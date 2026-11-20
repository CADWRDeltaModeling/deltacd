"""Pluggable meteorology providers and reference evapotranspiration (ETo) equations.

Implements the meteorological forcing and reference evapotranspiration equations
described in **Snyder et al. (2006), Appendix B (pp. 68-71)** and spatial
estimation rules (pp. 11-17):

- ``extraterrestrial_radiation`` (Duffie & Beckman, 1980; Eq. 1-4)
- ``hargreaves_samani`` (Hargreaves & Samani, 1982, 1985; Eq. 23)
- ``asce_penman_monteith`` (ASCE-EWRI, 2005 short-canopy daily PM; Eq. 5-22, 24-26)
- ``StationMeteorologyProvider`` (station ETo x subarea correction factor +
  Thiessen-weighted station precipitation)
- ``DatasetMeteorologyProvider`` (precomputed / gridded NetCDF or CSV ``et0`` &
  ``precip``)
"""

from __future__ import annotations

from pathlib import Path
from typing import Callable, Protocol

import numpy as np
import pandas as pd
import xarray as xr

from deltacd.detaw.config import DetawConfig
from deltacd.detaw.io import (
    _coerce_subarea_id,
    _extract_datetime_index,
    _find_column,
    read_station_precipitation,
    read_weather_station_csv,
)

# Physical constants (Snyder et al., 2006, Appendix B, pp. 69-71)
SOLAR_CONSTANT_MJ_M2_MIN: float = 0.082
STEFAN_BOLTZMANN_MJ_M2_D_K4: float = 4.90e-9
LATENT_HEAT_VAPORIZATION_MJ_KG: float = 2.45


def extraterrestrial_radiation(
    doy: xr.DataArray | np.ndarray,
    latitude_deg: float | xr.DataArray = 38.5,
) -> xr.DataArray:
    """Compute daily extraterrestrial solar radiation $R_a$ ($\\text{MJ m}^{-2}\\text{ d}^{-1}$).

    Follows Duffie & Beckman (1980) as presented in **Snyder et al. (2006),
    Appendix B, Step 1, Equations (1)-(4), p. 69**.

    Parameters
    ----------
    doy : xarray.DataArray or numpy.ndarray
        Calendar day of year (1..366).
    latitude_deg : float or xarray.DataArray, default 38.5
        Latitude in decimal degrees north (positive northern hemisphere).

    Returns
    -------
    xarray.DataArray
        Extraterrestrial radiation $R_a$ in $\\text{MJ m}^{-2}\\text{ d}^{-1}$.
    """
    doy_da = doy if isinstance(doy, xr.DataArray) else xr.DataArray(np.asarray(doy))
    phi = np.deg2rad(latitude_deg)

    # Eq. (1): Inverse relative Earth-sun distance factor d_r
    dr = 1.0 + 0.033 * np.cos((2.0 * np.pi / 365.0) * doy_da)
    # Eq. (2): Solar declination delta (radians)
    dec = 0.409 * np.sin((2.0 * np.pi / 365.0) * doy_da - 1.39)
    # Eq. (3): Sunrise hour angle omega_s (radians)
    arg = -np.tan(phi) * np.tan(dec)
    argsq = arg * arg
    # Match arccos(-tan(phi)*tan(dec)) safely across polar bounds
    ws = xr.where(
        argsq >= 1.0,
        np.pi / 2.0,
        (np.pi / 2.0) - np.arctan(arg / np.sqrt(np.maximum(1.0 - argsq, 1e-15))),
    )
    cosz = ws * np.sin(dec) * np.sin(phi) + np.cos(dec) * np.cos(phi) * np.sin(ws)
    # Eq. (4): Extraterrestrial radiation R_a
    ra = ((24.0 * 60.0) / np.pi) * SOLAR_CONSTANT_MJ_M2_MIN * dr * cosz
    ra.name = "ra"
    ra.attrs = {
        "units": "MJ m-2 day-1",
        "long_name": "extraterrestrial solar radiation",
    }
    return ra


def _prepare_temperatures(
    tmax: xr.DataArray,
    tmin: xr.DataArray,
    temp_qc_mode: str = "clamp",
) -> tuple[xr.DataArray, xr.DataArray]:
    """Apply quality-control rule when ``tmax < tmin``."""
    mode = temp_qc_mode.lower()
    if mode == "clamp":
        # Legacy DETAW baseline: set Tmax = Tmin where Tmax < Tmin
        tmax_clean = xr.where(tmax < tmin, tmin, tmax)
        tmin_clean = tmin
    elif mode == "swap":
        tmax_clean = np.maximum(tmax, tmin)
        tmin_clean = np.minimum(tmax, tmin)
    else:
        raise ValueError(f"Unsupported temp_qc_mode={temp_qc_mode!r}; expected 'clamp' or 'swap'")
    return tmax_clean, tmin_clean


def hargreaves_samani(
    ds_weather: xr.Dataset,
    latitude_deg: float | xr.DataArray = 38.5,
    temp_qc_mode: str = "clamp",
    **_: object,
) -> xr.DataArray:
    """Compute daily short-canopy reference ET ($ET_0$) via Hargreaves-Samani.

    Implements **Snyder et al. (2006), Appendix B, Step 4, Equation (23), p. 71**:

    $$ET_h = \\frac{1}{\\lambda} \\cdot 0.0023 \\cdot R_a \\cdot \\sqrt{T_x - T_n} \\cdot (T_m + 17.8)$$

    with $\\lambda = 2.45\\text{ MJ kg}^{-1}$.

    Parameters
    ----------
    ds_weather : xarray.Dataset
        Dataset with ``time`` coordinate and ``tmax``, ``tmin`` in °C.
    latitude_deg : float or xarray.DataArray, default 38.5
        Station/subarea latitude in decimal degrees.
    temp_qc_mode : {"clamp", "swap"}, default "clamp"
        Handling when ``tmax < tmin``.

    Returns
    -------
    xarray.DataArray
        Daily reference evapotranspiration $ET_0$ in $\\text{mm day}^{-1}$.
    """
    for req in ("tmax", "tmin"):
        if req not in ds_weather:
            raise ValueError(f"hargreaves_samani requires variable {req!r} in weather dataset")

    tmax, tmin = _prepare_temperatures(
        ds_weather["tmax"], ds_weather["tmin"], temp_qc_mode=temp_qc_mode
    )
    tdiff = (tmax - tmin).clip(min=0.0)
    tmean = 0.5 * (tmax + tmin)
    doy = ds_weather["time"].dt.dayofyear
    ra = extraterrestrial_radiation(doy, latitude_deg=latitude_deg)

    et0 = (0.0023 * ra * np.sqrt(tdiff) * (tmean + 17.8)) / LATENT_HEAT_VAPORIZATION_MJ_KG
    et0 = et0.clip(min=0.0)
    et0.name = "et0"
    et0.attrs = {
        "units": "mm day-1",
        "long_name": "reference evapotranspiration",
        "method": "hargreaves_samani",
    }
    return et0


def _saturation_vapor_pressure_kpa(temp_c: xr.DataArray) -> xr.DataArray:
    """Saturation vapor pressure $e_s(T)$ in kPa (Snyder et al., 2006, Eq. 8-9, p. 69)."""
    return 0.6108 * np.exp((17.27 * temp_c) / (temp_c + 237.3))


def asce_penman_monteith(
    ds_weather: xr.Dataset,
    latitude_deg: float | xr.DataArray = 38.5,
    elevation_m: float | xr.DataArray = 12.2,
    temp_qc_mode: str = "swap",
    **_: object,
) -> xr.DataArray:
    """Compute daily short-canopy (0.12 m) standardized Penman-Monteith $ET_0$.

    Implements **Snyder et al. (2006), Appendix B, Steps 1-3 & 5, Equations
    (5)-(22) and (24)-(26), pp. 69-71** (ASCE-EWRI, 2005).

    Required variables in ``ds_weather``:
    - ``tmax``, ``tmin`` (°C)
    - ``solar_rad`` ($R_s$, $\\text{MJ m}^{-2}\\text{ d}^{-1}$)
    - ``wind_speed`` ($u_2$ at 2 m height, $\\text{m s}^{-1}$)
    - Either ``dewpoint`` ($T_d$, °C; Eq. 11) or both ``rh_max`` and ``rh_min`` (%; Eq. 10)
    """
    for req in ("tmax", "tmin", "solar_rad", "wind_speed"):
        if req not in ds_weather:
            raise ValueError(f"asce_penman_monteith requires variable {req!r} in weather dataset")

    tmax, tmin = _prepare_temperatures(
        ds_weather["tmax"], ds_weather["tmin"], temp_qc_mode=temp_qc_mode
    )
    rs = ds_weather["solar_rad"].clip(min=0.0)
    u2 = ds_weather["wind_speed"].clip(min=0.0)
    doy = ds_weather["time"].dt.dayofyear
    ra = extraterrestrial_radiation(doy, latitude_deg=latitude_deg)

    # Step 2: Net radiation over grass (Eq. 5-14)
    rso = (0.75 + 2.0e-5 * elevation_m) * ra  # Eq. (5)
    rns = (1.0 - 0.23) * rs  # Eq. (6)
    rs_ratio = (rs / np.maximum(rso, 1e-6)).clip(min=0.3, max=1.0)
    f_cloud = 1.35 * rs_ratio - 0.35  # Eq. (7)

    es_tx = _saturation_vapor_pressure_kpa(tmax)  # Eq. (8)
    es_tn = _saturation_vapor_pressure_kpa(tmin)  # Eq. (9)

    if "dewpoint" in ds_weather:
        ea = _saturation_vapor_pressure_kpa(ds_weather["dewpoint"])  # Eq. (11)
    elif "rh_max" in ds_weather and "rh_min" in ds_weather:
        ea = 0.5 * (
            es_tn * (ds_weather["rh_max"].clip(0.0, 100.0) / 100.0)
            + es_tx * (ds_weather["rh_min"].clip(0.0, 100.0) / 100.0)
        )  # Eq. (10)
    else:
        raise ValueError(
            "asce_penman_monteith requires either 'dewpoint' or ('rh_max' and 'rh_min')"
        )

    eps_net = 0.34 - 0.14 * np.sqrt(ea.clip(min=0.0))  # Eq. (12)
    rnl = (
        -f_cloud
        * eps_net
        * STEFAN_BOLTZMANN_MJ_M2_D_K4
        * 0.5
        * ((tmax + 273.15) ** 4 + (tmin + 273.15) ** 4)
    )  # Eq. (13)
    rn = rns + rnl  # Eq. (14)

    # Step 3: Psychrometric & vapor-pressure slope variables (Eq. 15-22)
    pressure_kpa = 101.3 * ((293.0 - 0.0065 * elevation_m) / 293.0) ** 5.26  # Eq. (15)
    gamma = (0.00163 * pressure_kpa) / LATENT_HEAT_VAPORIZATION_MJ_KG  # Eq. (17)
    tmean = 0.5 * (tmax + tmin)  # Eq. (18)
    eo = _saturation_vapor_pressure_kpa(tmean)  # Eq. (19)
    delta = (4099.0 * eo) / ((tmean + 237.3) ** 2)  # Eq. (20)
    g_soil = 0.0  # Eq. (21)
    es = 0.5 * (es_tx + es_tn)  # Eq. (22)

    # Step 5: Radiation + Aerodynamic terms (Eq. 24-26)
    denom = delta + gamma * (1.0 + 0.34 * u2)
    inv_lambda = 1.0 / LATENT_HEAT_VAPORIZATION_MJ_KG  # ~0.408163
    r_o = (inv_lambda * delta * (rn - g_soil)) / denom  # Eq. (24)
    a_o = (gamma * (900.0 / (tmean + 273.0)) * u2 * (es - ea).clip(min=0.0)) / denom  # Eq. (25)

    et0 = (r_o + a_o).clip(min=0.0)  # Eq. (26)
    et0.name = "et0"
    et0.attrs = {
        "units": "mm day-1",
        "long_name": "reference evapotranspiration",
        "method": "asce_penman_monteith",
    }
    return et0


ETO_EQUATIONS: dict[str, Callable[..., xr.DataArray]] = {
    "hargreaves_samani": hargreaves_samani,
    "hs": hargreaves_samani,
    "asce_penman_monteith": asce_penman_monteith,
    "penman_monteith": asce_penman_monteith,
    "pm": asce_penman_monteith,
}


class MeteorologyProvider(Protocol):
    """Protocol for subarea meteorological forcing providers."""

    def load(
        self,
        start_water_year: int,
        end_water_year: int,
        ds_subarea: xr.Dataset,
    ) -> xr.Dataset:
        """Return ``xr.Dataset`` with ``et0(time, subarea)`` and ``precip(time, subarea)`` in mm/d."""
        ...


class StationMeteorologyProvider:
    """Station-based meteorology provider (default DETAW meteorology).

    Computes station $ET_0$ using a pluggable equation (Hargreaves-Samani or
    ASCE Penman-Monteith, or direct station ``et0`` column), scales by subarea
    ``eto_correction``, and computes subarea precipitation via ``xr.dot`` over
    normalized ``station`` coordinates with proper ``[0, 9990]`` clipping.

    Supports the optional 1-day spin-up on ``(start_water_year - 1)-09-30``
    automatically when present in the CSVs, without failing if CSVs start on
    ``(start_water_year - 1)-10-01``.
    """

    def __init__(self, config: DetawConfig):
        self.config = config

    def load(
        self,
        start_water_year: int,
        end_water_year: int,
        ds_subarea: xr.Dataset,
    ) -> xr.Dataset:
        met = self.config.meteorology
        if not met.input_pcp or not met.input_temperature:
            raise ValueError(
                "StationMeteorologyProvider requires 'input_pcp' and 'input_temperature' paths"
            )

        spinup_date = pd.Timestamp(f"{start_water_year - 1}-09-30")
        sim_start = pd.Timestamp(f"{start_water_year - 1}-10-01")
        sim_end = pd.Timestamp(f"{end_water_year}-09-30")

        da_station_pcp = read_station_precipitation(
            met.input_pcp,
            start_date=spinup_date,
            end_date=sim_end,
            station_aliases=met.station_aliases,
            column_map=self.config.column_map,
        )
        ds_weather = read_weather_station_csv(
            met.input_temperature,
            start_date=spinup_date,
            end_date=sim_end,
            column_map=self.config.column_map,
        )

        # Align shared time range (supports optional Sep 30 spin-up gracefully)
        has_spinup = spinup_date in pd.DatetimeIndex(
            da_station_pcp["time"].values
        ) and spinup_date in pd.DatetimeIndex(ds_weather["time"].values)
        effective_start = spinup_date if has_spinup else sim_start
        expected_index = pd.date_range(effective_start, sim_end, freq="D", name="time")

        if not expected_index.isin(pd.DatetimeIndex(da_station_pcp["time"].values)).all():
            missing = expected_index.difference(pd.DatetimeIndex(da_station_pcp["time"].values))
            raise ValueError(
                f"Precipitation CSV {met.input_pcp} is missing {len(missing)} dates "
                f"(first missing: {missing[0].date()})"
            )
        if not expected_index.isin(pd.DatetimeIndex(ds_weather["time"].values)).all():
            missing = expected_index.difference(pd.DatetimeIndex(ds_weather["time"].values))
            raise ValueError(
                f"Temperature/weather CSV {met.input_temperature} is missing {len(missing)} dates "
                f"(first missing: {missing[0].date()})"
            )

        da_station_pcp = da_station_pcp.sel(time=expected_index)
        ds_weather = ds_weather.sel(time=expected_index)

        # Validate station coordinate match before xr.dot
        weights = ds_subarea["precip_weights"]
        pcp_stations = set(da_station_pcp["station"].values.tolist())
        # Any station with non-zero weight must exist in the precipitation CSV
        nonzero_mask = (np.abs(weights) > 0.0).any(dim="subarea")
        required_stations = set(weights["station"].values[nonzero_mask.values].tolist())
        missing_stations = sorted(required_stations - pcp_stations)
        if missing_stations:
            raise ValueError(
                f"Stations with non-zero weights in et_correction are missing from "
                f"precipitation CSV: {missing_stations}. Available precip stations: {sorted(pcp_stations)}"
            )

        shared_stations = [s for s in weights["station"].values.tolist() if s in pcp_stations]
        da_subarea_pcp = xr.dot(
            da_station_pcp.sel(station=shared_stations),
            weights.sel(station=shared_stations),
            dim="station",
        ).transpose("time", "subarea")
        # Category A3 fix: actually assign the clipped precipitation array!
        if not self.config.options.emulate_legacy_bugs:
            da_subarea_pcp = da_subarea_pcp.clip(min=0.0, max=9990.0)
        da_subarea_pcp.name = "precip"
        da_subarea_pcp.attrs = {"units": "mm day-1", "long_name": "subarea daily precipitation"}

        # Compute reference ETo
        method = met.eto_method.lower()
        if method in ("direct", "column") and "et0" in ds_weather:
            da_station_et0 = ds_weather["et0"].clip(min=0.0)
        else:
            if method not in ETO_EQUATIONS:
                raise ValueError(
                    f"Unknown eto_method={met.eto_method!r}. Available: {sorted(ETO_EQUATIONS)}"
                )
            eto_fn = ETO_EQUATIONS[method]
            da_station_et0 = eto_fn(
                ds_weather,
                latitude_deg=met.latitude_deg,
                elevation_m=met.elevation_m,
                temp_qc_mode=met.temp_qc_mode,
            )

        da_subarea_et0 = (da_station_et0 * ds_subarea["eto_correction"]).transpose(
            "time", "subarea"
        )
        da_subarea_et0 = da_subarea_et0.clip(min=0.0)
        da_subarea_et0.name = "et0"
        da_subarea_et0.attrs = {
            "units": "mm day-1",
            "long_name": "subarea reference evapotranspiration",
            "method": method,
        }

        return xr.Dataset(
            data_vars={"et0": da_subarea_et0, "precip": da_subarea_pcp},
            coords={"time": expected_index.values, "subarea": ds_subarea["subarea"]},
            attrs={"has_spinup_day": int(has_spinup)},
        )


class DatasetMeteorologyProvider:
    """Precomputed or gridded subarea meteorology provider (NetCDF or CSV).

    Loads ``et0`` and ``precip`` directly by ``(time, subarea)`` coordinates,
    allowing users to drive DETAW with spatial CIMIS, GridMET, PRISM, or
    climate-change projections without station interpolation.
    """

    def __init__(self, config: DetawConfig, dataset: xr.Dataset | None = None):
        self.config = config
        self._dataset = dataset

    def load(
        self,
        start_water_year: int,
        end_water_year: int,
        ds_subarea: xr.Dataset,
    ) -> xr.Dataset:
        spinup_date = pd.Timestamp(f"{start_water_year - 1}-09-30")
        sim_start = pd.Timestamp(f"{start_water_year - 1}-10-01")
        sim_end = pd.Timestamp(f"{end_water_year}-09-30")

        if self._dataset is not None:
            ds_raw = self._dataset.copy()
        else:
            path_str = self.config.meteorology.dataset_path
            if not path_str:
                raise ValueError("DatasetMeteorologyProvider requires 'meteorology.dataset_path'")
            path = Path(path_str)
            if path.suffix.lower() in (".nc", ".nc4", ".netcdf"):
                with xr.open_dataset(path) as opened:
                    ds_raw = opened.load()
            else:
                df = pd.read_csv(path)
                dt_idx, time_cols = _extract_datetime_index(
                    df, path, column_map=self.config.column_map
                )
                area_col = _find_column(
                    df.columns,
                    ("subarea", "area_id", "subarea_id"),
                    column_map=self.config.column_map,
                    context=str(path),
                )
                et0_col = _find_column(
                    df.columns,
                    ("et0", "eto", "et_0"),
                    column_map=self.config.column_map,
                    context=str(path),
                )
                pcp_col = _find_column(
                    df.columns,
                    ("precip", "pcp", "ppt", "prcp"),
                    column_map=self.config.column_map,
                    context=str(path),
                )
                assert area_col is not None and et0_col is not None and pcp_col is not None
                work = df.copy()
                work["time"] = dt_idx
                work[area_col] = _coerce_subarea_id(work[area_col])
                p_et0 = work.pivot(index="time", columns=area_col, values=et0_col)
                p_pcp = work.pivot(index="time", columns=area_col, values=pcp_col)
                ds_raw = xr.Dataset(
                    data_vars={
                        "et0": (("time", "subarea"), p_et0.to_numpy(dtype=np.float64)),
                        "precip": (("time", "subarea"), p_pcp.to_numpy(dtype=np.float64)),
                    },
                    coords={
                        "time": p_et0.index.values,
                        "subarea": np.asarray(p_et0.columns.tolist()),
                    },
                )

        # Normalize variable names (e.g. ET0 -> et0)
        rename_map: dict[str, str] = {}
        for v in ds_raw.data_vars:
            if str(v).lower() in ("et0", "eto", "et_0") and v != "et0":
                rename_map[str(v)] = "et0"
            elif str(v).lower() in ("precip", "pcp", "ppt") and v != "precip":
                rename_map[str(v)] = "precip"
        if rename_map:
            ds_raw = ds_raw.rename(rename_map)

        for req in ("et0", "precip"):
            if req not in ds_raw:
                raise ValueError(f"Meteorology dataset is missing required variable {req!r}")

        target_subareas = ds_subarea["subarea"].values
        ds_aligned = ds_raw.sel(subarea=target_subareas)

        time_idx = pd.DatetimeIndex(ds_aligned["time"].values)
        has_spinup = spinup_date in time_idx
        effective_start = spinup_date if has_spinup else sim_start
        expected_index = pd.date_range(effective_start, sim_end, freq="D", name="time")
        ds_aligned = ds_aligned.sel(time=expected_index)

        et0 = ds_aligned["et0"].transpose("time", "subarea").clip(min=0.0)
        precip = ds_aligned["precip"].transpose("time", "subarea").clip(min=0.0, max=9990.0)
        return xr.Dataset(
            data_vars={"et0": et0, "precip": precip},
            coords={"time": expected_index.values, "subarea": ds_subarea["subarea"]},
            attrs={"has_spinup_day": int(has_spinup)},
        )


def create_meteorology_provider(config: DetawConfig) -> MeteorologyProvider:
    """Instantiate the configured ``MeteorologyProvider``."""
    provider_name = config.meteorology.provider.lower()
    if provider_name in ("station", "hargreaves_station", "legacy"):
        return StationMeteorologyProvider(config)
    if provider_name in ("dataset", "precomputed", "gridded", "netcdf"):
        return DatasetMeteorologyProvider(config)
    raise ValueError(
        f"Unknown meteorology provider {config.meteorology.provider!r}. "
        f"Supported: 'station', 'dataset'"
    )
