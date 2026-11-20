"""Crop coefficient ($K_c$), bare-soil evaporation ($OK_c$), and OpenET ($ET_c$) providers.

Implements:
1. Two-stage bare-soil evaporation model (`compute_bare_soil_okc`) from Ventura
   et al. (2006) & Stroosnijder (1987) as described in **Snyder et al. (2006),
   Appendix C, pp. 72–76**.
2. Modified Doorenbos & Pruitt (1977) seasonal crop coefficient curves (`IKc`,
   `Kc`) across growth stages A–B–C–D–E (**Snyder et al., 2006, pp. 18–26,
   76–83**) over real Gregorian calendar dates (including leap years).
3. ``KcETcProvider``: Standard internal ``ETc = Kc * ETo`` calculation.
4. ``ExternalETProvider``: Pluggable external evapotranspiration provider
   (e.g., **OpenET** daily or monthly NetCDF / CSV / Dataset for ``et_c`` or
   direct ``et_aw``).
"""

from __future__ import annotations

from pathlib import Path
from typing import Protocol

import numba
import numpy as np
import pandas as pd
import xarray as xr

from deltacd.detaw.config import DetawConfig
from deltacd.detaw.io import _coerce_subarea_id, _extract_datetime_index, _find_column
from deltacd.detaw.params import CropBehavior, resolve_crop_behavior


@numba.njit(cache=True)
def _bare_soil_okc_kernel(
    et0_2d: np.ndarray,
    pcp_2d: np.ndarray,
    beta: float,
    daily_increment: bool,
    two_eto_threshold: bool,
    truncate_2dec: bool,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Compute daily bare-soil evaporation coefficient ``OKc`` on ``(time, subarea)``.

    Parameters
    ----------
    et0_2d, pcp_2d : float64 array of shape (n_days, n_subareas)
        Daily reference ETo and precipitation in mm/d (including optional
        Sep 30 spin-up day at index 0).
    beta : float
        Stage-1 / Stage-2 hydraulic transition parameter ($\\beta = 2.6$,
        Snyder et al., 2006, p. 72).
    daily_increment : bool
        If False, uses legacy DETAW cumulative ratio ``CEs / CETo`` (with
        initial ``DCT = 1`` incremented to ``2`` on day 0). If True, uses
        Snyder et al. (pp. 73–74) daily increment ``(CEs[t] - CEs[t-1]) / ETo[t]``.
    two_eto_threshold : bool
        If False, significant wetting rainfall occurs when ``P > METo`` (legacy
        and Snyder et al. p. 73). If True, ``P > 2.0 * ETo`` (p. 76).
    truncate_2dec : bool
        Diagnostic legacy-bug emulation toggle reproducing C++ 2-decimal
        rounding ``int(x*100 + 0.5)/100`` and carrying previous ``ETo`` when
        ``ETo <= 0``.
    """
    n_days, n_subareas = et0_2d.shape
    okc = np.zeros((n_days, n_subareas), dtype=np.float64)
    et0_used = np.empty((n_days, n_subareas), dtype=np.float64)
    pcp_used = np.empty((n_days, n_subareas), dtype=np.float64)

    for s in range(n_subareas):
        prev_eto = 0.0
        c_eto = 0.0
        prev_ces = 0.0
        dct = 1 if not daily_increment else 0

        for t in range(n_days):
            e_val = et0_2d[t, s]
            p_val = pcp_2d[t, s]

            if truncate_2dec:
                if (e_val * 100.0 - int(e_val * 100.0)) > 0.5:
                    e_val = int(e_val * 100.0 + 1) / 100.0
                else:
                    e_val = int(e_val * 100.0) / 100.0

                if (p_val * 100.0 - int(p_val * 100.0)) > 0.5:
                    p_val = int(p_val * 100.0 + 1) / 100.0
                else:
                    p_val = int(p_val * 100.0) / 100.0

            pcp_used[t, s] = p_val

            dct += 1
            if e_val <= 0.0:
                e_val = prev_eto
            et0_used[t, s] = e_val

            c_eto += e_val
            m_eto = c_eto / dct

            is_wet = (p_val > 2.0 * e_val) if two_eto_threshold else (p_val > m_eto)
            if is_wet:
                c_eto = e_val
                dct = 1
                prev_ces = 0.0
                m_eto = c_eto / dct

            kx = 1.22 - 0.04 * m_eto
            c_ex = kx * c_eto
            if c_ex < 0.0:
                c_ex = 0.0

            sqrt_cex = np.sqrt(c_ex)
            if sqrt_cex < beta:
                c_es = c_ex
            else:
                c_es = beta * sqrt_cex

            if e_val != 0.0:
                if daily_increment:
                    es_day = c_es - prev_ces
                    if es_day < 0.0:
                        es_day = 0.0
                    okc[t, s] = es_day / e_val
                else:
                    okc[t, s] = c_es / c_eto

            prev_ces = c_es
            prev_eto = e_val

    return okc, et0_used, pcp_used


def compute_bare_soil_okc(
    ds_meteo: xr.Dataset,
    beta: float = 2.6,
    okc_mode: str = "legacy_cumulative",
    wet_threshold_mode: str = "mean_eto",
    truncate_2dec: bool = False,
) -> tuple[xr.DataArray, xr.DataArray, xr.DataArray]:
    """Compute daily bare-soil crop coefficient ``OKc(time, subarea)``.

    Returns ``(okc, et0_effective, precip_effective)`` on ``(time, subarea)``.
    """
    daily_inc = okc_mode.lower() == "daily_increment"
    if okc_mode.lower() not in ("legacy_cumulative", "daily_increment"):
        raise ValueError(
            f"Unsupported okc_mode={okc_mode!r}; expected 'legacy_cumulative' or 'daily_increment'"
        )
    two_eto = wet_threshold_mode.lower() == "two_eto"

    et0_in = np.ascontiguousarray(
        ds_meteo["et0"].transpose("time", "subarea").values, dtype=np.float64
    )
    pcp_in = np.ascontiguousarray(
        ds_meteo["precip"].transpose("time", "subarea").values, dtype=np.float64
    )

    okc_arr, et0_used, pcp_used = _bare_soil_okc_kernel(
        et0_in,
        pcp_in,
        float(beta),
        daily_inc,
        two_eto,
        bool(truncate_2dec),
    )

    coords = {"time": ds_meteo["time"], "subarea": ds_meteo["subarea"]}
    da_okc = xr.DataArray(
        okc_arr,
        dims=("time", "subarea"),
        coords=coords,
        name="okc",
        attrs={"units": "1", "long_name": "off-season bare-soil crop coefficient"},
    )
    da_et0 = xr.DataArray(
        et0_used,
        dims=("time", "subarea"),
        coords=coords,
        name="et0",
        attrs={"units": "mm day-1", "long_name": "effective reference evapotranspiration"},
    )
    da_pcp = xr.DataArray(
        pcp_used,
        dims=("time", "subarea"),
        coords=coords,
        name="precip",
        attrs={"units": "mm day-1", "long_name": "effective daily precipitation input"},
    )
    return da_okc, da_et0, da_pcp


def _is_critical_year_type(yt: str) -> bool:
    """Return True if a water-year classification is Critical ('C') or Dry ('D')."""
    return str(yt).strip().upper() in ("C", "D")


def _compute_piecewise_season_ikc(
    begin_doy: int,
    end_doy: int,
    ab_pct: float,
    ac_pct: float,
    ad_pct: float,
    kc_b: float,
    kc_c: float,
    kc_d: float,
    kc_e: float,
) -> tuple[np.ndarray, np.ndarray, int, int, int, int]:
    """Compute the in-season piecewise-linear ``IKc`` curve for one planting season.

    Follows **Snyder et al. (2006), Appendix C, pp. 76–78** and exact slope
    stepping from ``calc_ikc`` in DETAW:
    - Dates ``A..B``: constant ``kc_b``
    - Dates ``B+1..C``: linear from ``kc_b`` to ``kc_c``
    - Dates ``C+1..D``: linear from ``kc_c`` to ``kc_d``
    - Dates ``D+1..E``: linear from ``kc_d`` to ``kc_e``

    Returns ``(season_jj, season_ikc, B, C, D, end_adj)``.
    """
    end_adj = int(end_doy)
    if end_adj <= begin_doy:
        end_adj = 365 + end_adj
    len_date = end_adj - begin_doy

    b_day = int(0.01 * ab_pct * len_date + begin_doy)
    c_day = int(0.01 * ac_pct * len_date + begin_doy)
    d_day = int(0.01 * ad_pct * len_date + begin_doy)

    bc_slope = (kc_c - kc_b) / (c_day - b_day) if c_day > b_day else 0.0
    cd_slope = (kc_d - kc_c) / (d_day - c_day) if d_day > c_day else 0.0
    de_slope = (kc_e - kc_d) / (end_adj - d_day) if end_adj > d_day else 0.0

    jj = np.arange(begin_doy, end_adj + 1, dtype=np.int32)
    vals = np.empty(len(jj), dtype=np.float64)
    cur = kc_b
    for idx, day_num in enumerate(jj):
        if day_num <= d_day:
            if day_num <= c_day:
                if day_num <= b_day:
                    cur = kc_b
                else:
                    cur = cur + bc_slope
            else:
                cur = cur + cd_slope
        else:
            cur = cur + de_slope
        vals[idx] = cur

    return jj, vals, b_day, c_day, d_day, end_adj


def build_daily_ikc_and_kc(
    time_index: pd.DatetimeIndex,
    da_okc: xr.DataArray,
    da_et0: xr.DataArray,
    ds_crop_params: xr.Dataset,
    da_year_type: xr.DataArray,
    pad_year_type: str = "AN",
    season_mode: str = "legacy_wrapped",
    type1_initial_kc_mode: str = "legacy",
    beta: float = 2.6,
    emulate_legacy_leap_ikc_copy: bool = False,
) -> tuple[xr.DataArray, xr.DataArray]:
    """Build daily in-season ``IKc(time, crop)`` and combined ``Kc(time, subarea, crop)``.

    Parameters
    ----------
    time_index : pandas.DatetimeIndex
        Simulation daily timestamps (e.g., ``1921-10-01 .. 1932-09-30``).
    da_okc : xarray.DataArray
        Bare-soil evaporation coefficient ``OKc(time, subarea)`` aligned to ``time_index``.
    da_et0 : xarray.DataArray
        Daily reference ``ETo(time, subarea)`` aligned to ``time_index``.
    ds_crop_params : xarray.Dataset
        Tidy crop parameters & declarative behavior flags from ``read_crop_params``.
    da_year_type : xarray.DataArray
        Water-year classification ``year_type(water_year)`` from ``read_landuse``.
    pad_year_type : str, default "AN"
        Fallback water-year type for boundary calendar years before/after the
        land-use record (matches legacy ``"AN"`` padding).
    season_mode : {"legacy_wrapped", "continuous"}, default "legacy_wrapped"
        - ``"legacy_wrapped"``: Matches legacy DETAW calendar-year wrapping
          where wrapped seasons (`end_doy > 365`, e.g. fall-planted grain) wrap
          into Jan–May of the *same* calendar year.
        - ``"continuous"``: Continuous calendar timeline where a season planted
          in fall of calendar year ``Y`` continues seamlessly into Jan–May of
          calendar year ``Y + 1`` (**Category B1 option**).
    type1_initial_kc_mode : {"legacy", "irrigation_frequency"}, default "legacy"
        Optional SIMETAW Type-1 initial-growth ``KcB`` bump from irrigation
        frequency ``F`` (**Category B5 option**).
    emulate_legacy_leap_ikc_copy : bool, default False
        Diagnostic toggle reproducing legacy ``IKc[y, 366] = IKc[y, 365]`` on
        leap years.
    """
    n_days = len(time_index)
    n_crops = ds_crop_params.sizes["crop"]
    n_subareas = da_okc.sizes["subarea"]

    wy_lookup: dict[int, str] = {
        int(wy): str(yt) for wy, yt in zip(da_year_type["water_year"].values, da_year_type.values)
    }

    cal_years = sorted(time_index.year.unique().tolist())
    # Include cal_years[0] - 1 when running continuous cross-year seasons so
    # fall plantings from the preceding year populate Jan–May of the first year.
    start_cal_year = cal_years[0] - 1 if season_mode.lower() == "continuous" else cal_years[0]
    all_cal_years = list(range(start_cal_year, cal_years[-1] + 1))

    # Map each timestamp to integer position for O(1) slice assignment
    date_to_pos = {ts: idx for idx, ts in enumerate(time_index)}

    # Optional per-subarea Type-1 initial Kc requires (time, subarea, crop) IKc;
    # standard mode only depends on (time, crop), which is 168x faster and broadcasts cleanly!
    use_freq_bump = type1_initial_kc_mode.lower() == "irrigation_frequency"
    if use_freq_bump:
        ikc_full = np.zeros((n_days, n_subareas, n_crops), dtype=np.float64)
    else:
        ikc_2d = np.zeros((n_days, n_crops), dtype=np.float64)

    unshifted_flags = ds_crop_params["year_type_unshifted_in_legacy"].values
    open_water_flags = ds_crop_params["is_open_water"].values
    open_water_kcs = ds_crop_params["open_water_kc"].values

    for c_idx in range(n_crops):
        unshifted = bool(unshifted_flags[c_idx])

        for cy in all_cal_years:
            dpy = 366 if pd.Timestamp(year=cy, month=1, day=1).is_leap_year else 365

            if season_mode.lower() == "continuous":
                # Determine water year containing the planting date
                nc_bd = int(
                    ds_crop_params["begin_doy"].sel(regime="noncritical").isel(crop=c_idx).item()
                )
                plant_ts = pd.Timestamp(year=cy, month=1, day=1) + pd.Timedelta(days=nc_bd - 1)
                wy_for_season = plant_ts.year + 1 if plant_ts.month >= 10 else plant_ts.year
                yt = wy_lookup.get(wy_for_season, pad_year_type)
            else:
                # Legacy calendar-year water-year-type mapping:
                # year_i = 1 corresponds to model_start_year (e.g. 1921 for WY 1922).
                # yearType[1] is WY 1922 = cy + 1!
                # So unshifted (GR, DGR): yearType[y] -> WY (cy + 1);
                # shifted (all other crops): yearType[y - 1] -> WY cy.
                wy_for_season = (cy + 1) if unshifted else cy
                yt = wy_lookup.get(wy_for_season, pad_year_type)

            regime = "critical" if _is_critical_year_type(yt) else "noncritical"
            p_reg = ds_crop_params.sel(regime=regime).isel(crop=c_idx)

            crop_type = int(p_reg["crop_type"].item())
            bd = int(p_reg["begin_doy"].item())
            ed = int(p_reg["end_doy"].item())
            irrig_f = int(p_reg["irrig_freq"].item()) or 30
            kc1 = float(p_reg["kc1"].item())
            kc2 = float(p_reg["kc2"].item())
            kc3 = float(p_reg["kc3"].item())
            ab = float(p_reg["ab_pct"].item())
            ac = float(p_reg["ac_pct"].item())
            ad = float(p_reg["ad_pct"].item())

            # On non-leap years, full-year crops have ED=366 in legacy CSVs;
            # on leap years, PA non-critical has ED=365 in legacy CSVs.
            if not emulate_legacy_leap_ikc_copy and bd == 1 and ed in (365, 366):
                ed_eff = dpy
            else:
                ed_eff = ed

            if use_freq_bump and crop_type == 1:
                # Compute per-subarea initial growth KcB from irrigation frequency F
                end_tmp = ed_eff + 365 if ed_eff <= bd else ed_eff
                b_tmp = int(0.01 * ab * (end_tmp - bd) + bd)
                q_ts = pd.Timestamp(year=cy, month=1, day=1) + pd.Timedelta(days=bd - 1)
                r_ts = pd.Timestamp(year=cy, month=1, day=1) + pd.Timedelta(days=b_tmp - 1)
                et0_slice = da_et0.sel(time=slice(q_ts, r_ts))
                if et0_slice.sizes["time"] > 0:
                    ig_eto = et0_slice.mean(dim="time").values  # (subarea,)
                else:
                    ig_eto = np.zeros(n_subareas, dtype=np.float64)
                c_eto = irrig_f * ig_eto
                kx = 1.22 - 0.04 * ig_eto
                c_ex = np.maximum(kx * c_eto, 0.0)
                c_es = np.where(np.sqrt(c_ex) < beta, c_ex, beta * np.sqrt(c_ex))
                is_ikc = np.where(c_eto > 0.0, c_es / c_eto, 0.0)
                kc_b_per_sub = np.maximum(kc1, is_ikc)

                for s_idx in range(n_subareas):
                    jj, vals, _, _, _, _ = _compute_piecewise_season_ikc(
                        bd, ed_eff, ab, ac, ad, kc_b_per_sub[s_idx], kc2, kc2, kc3
                    )
                    for day_num, val in zip(jj, vals):
                        if season_mode.lower() == "continuous":
                            ts = pd.Timestamp(year=cy, month=1, day=1) + pd.Timedelta(
                                days=int(day_num) - 1
                            )
                        else:
                            doy_wrapped = int(day_num) - dpy if int(day_num) > dpy else int(day_num)
                            if doy_wrapped < 1 or doy_wrapped > dpy:
                                continue
                            ts = pd.Timestamp(year=cy, month=1, day=1) + pd.Timedelta(
                                days=doy_wrapped - 1
                            )
                        pos = date_to_pos.get(ts)
                        if pos is not None:
                            ikc_full[pos, s_idx, c_idx] = val
                continue

            jj, vals, _, _, _, _ = _compute_piecewise_season_ikc(
                bd, ed_eff, ab, ac, ad, kc1, kc2, kc2, kc3
            )

            if season_mode.lower() == "continuous":
                for day_num, val in zip(jj, vals):
                    ts = pd.Timestamp(year=cy, month=1, day=1) + pd.Timedelta(days=int(day_num) - 1)
                    pos = date_to_pos.get(ts)
                    if pos is not None:
                        ikc_2d[pos, c_idx] = val
            else:
                year_curve = np.zeros(dpy + 1, dtype=np.float64)
                for day_num, val in zip(jj, vals):
                    doy_wrapped = int(day_num) - dpy if int(day_num) > dpy else int(day_num)
                    if 1 <= doy_wrapped <= dpy:
                        year_curve[doy_wrapped] = val
                if emulate_legacy_leap_ikc_copy and dpy == 366:
                    year_curve[366] = year_curve[365]

                for doy_idx in range(1, dpy + 1):
                    ts = pd.Timestamp(year=cy, month=1, day=1) + pd.Timedelta(days=doy_idx - 1)
                    pos = date_to_pos.get(ts)
                    if pos is not None:
                        ikc_2d[pos, c_idx] = year_curve[doy_idx]

    okc_vals = da_okc.transpose("time", "subarea").values  # (time, subarea)
    if use_freq_bump:
        da_ikc = xr.DataArray(
            ikc_full,
            dims=("time", "subarea", "crop"),
            coords={
                "time": time_index,
                "subarea": da_okc["subarea"],
                "crop": ds_crop_params["crop"],
            },
            name="ikc",
            attrs={"units": "1", "long_name": "in-season crop coefficient"},
        )
        kc_3d = np.maximum(ikc_full, okc_vals[:, :, np.newaxis])
    else:
        da_ikc = xr.DataArray(
            ikc_2d,
            dims=("time", "crop"),
            coords={"time": time_index, "crop": ds_crop_params["crop"]},
            name="ikc",
            attrs={"units": "1", "long_name": "in-season crop coefficient"},
        )
        kc_3d = np.maximum(ikc_2d[:, np.newaxis, :], okc_vals[:, :, np.newaxis])

    for c_idx in range(n_crops):
        if open_water_flags[c_idx]:
            kc_3d[:, :, c_idx] = float(open_water_kcs[c_idx])

    da_kc = xr.DataArray(
        kc_3d,
        dims=("time", "subarea", "crop"),
        coords={
            "time": time_index,
            "subarea": da_okc["subarea"],
            "crop": ds_crop_params["crop"],
        },
        name="kc",
        attrs={"units": "1", "long_name": "combined daily crop/surface coefficient"},
    )
    return da_ikc, da_kc


class ETcProvider(Protocol):
    """Protocol for computing or supplying daily crop evapotranspiration ($ET_c$)."""

    def compute(
        self,
        time_index: pd.DatetimeIndex,
        ds_meteo_full: xr.Dataset,
        ds_crop_params: xr.Dataset,
        da_year_type: xr.DataArray,
    ) -> xr.Dataset:
        """Return ``xr.Dataset`` on ``time_index`` containing ``et_c_mm``, ``ikc``, ``kc``, ``okc``, ``et0``, ``precip``."""
        ...


class KcETcProvider:
    """Standard DETAW $ET_c = K_c \\times ET_0$ provider."""

    def __init__(self, config: DetawConfig):
        self.config = config

    def compute(
        self,
        time_index: pd.DatetimeIndex,
        ds_meteo_full: xr.Dataset,
        ds_crop_params: xr.Dataset,
        da_year_type: xr.DataArray,
    ) -> xr.Dataset:
        da_okc_full, da_et0_full, da_pcp_full = compute_bare_soil_okc(
            ds_meteo_full,
            beta=self.config.et.beta,
            okc_mode=self.config.et.okc_mode,
            wet_threshold_mode=self.config.et.wet_threshold_mode,
            truncate_2dec=self.config.options.emulate_legacy_bugs,
        )

        da_okc = da_okc_full.sel(time=time_index)
        da_et0 = da_et0_full.sel(time=time_index)
        da_pcp = da_pcp_full.sel(time=time_index)

        da_ikc, da_kc = build_daily_ikc_and_kc(
            time_index=time_index,
            da_okc=da_okc,
            da_et0=da_et0,
            ds_crop_params=ds_crop_params,
            da_year_type=da_year_type,
            pad_year_type=self.config.options.default_pad_year_type,
            season_mode=self.config.et.season_mode,
            type1_initial_kc_mode=self.config.et.type1_initial_kc_mode,
            beta=self.config.et.beta,
            emulate_legacy_leap_ikc_copy=self.config.options.emulate_legacy_bugs,
        )

        da_etc_mm = (da_kc * da_et0).transpose("time", "subarea", "crop")
        da_etc_mm.name = "et_c_mm"
        da_etc_mm.attrs = {
            "units": "mm day-1",
            "long_name": "daily crop evapotranspiration depth",
        }

        return xr.Dataset(
            data_vars={
                "et_c_mm": da_etc_mm,
                "ikc": da_ikc,
                "kc": da_kc,
                "okc": da_okc,
                "et0": da_et0,
                "precip": da_pcp,
            },
            coords={
                "time": time_index,
                "subarea": da_et0["subarea"],
                "crop": ds_crop_params["crop"],
                "crop_code": ds_crop_params["crop_code"],
            },
        )


class ExternalETProvider:
    """External / OpenET crop evapotranspiration ($ET_c$ or $ET_{aw}$) provider.

    Allows users to bypass internal $K_c \\times ET_0$ calculations and supply
    externally observed or remotely sensed ET (such as **OpenET** Ensemble,
    SSEBop, eeMETRIC, SIMS, PT-JPL, or DisALEXI) from a NetCDF file, CSV file,
    or in-memory ``xarray.Dataset``:
    - Supports ``(time, subarea, crop)`` crop-specific ET or ``(time, subarea)``
      subarea-composite ET (broadcast across crops).
    - Supports daily ET (``temporal_resolution="daily"``) or monthly ET
      (``temporal_resolution="monthly"``) disaggregated to daily steps either
      uniformly (``disaggregation="uniform"``) or scaled by daily $ET_0$
      pattern (``disaggregation="et0_scaled"``).
    - Supports ``mm/d``, ``mm/month`` (for monthly totals), or ``in/d`` units.
    - Also computes seasonal ``ikc`` mask from ``ds_crop_params`` soDETAW's
      root-zone water balance knows in-season vs. off-season irrigation windows
      (or accepts direct ``et_aw`` when ``external_variable == "et_aw"``).
    """

    def __init__(
        self,
        config: DetawConfig,
        dataset: xr.Dataset | None = None,
        extra_behaviors: tuple[CropBehavior, ...] = (),
    ):
        self.config = config
        self._dataset = dataset
        self._extra_behaviors = extra_behaviors

    def _load_external_dataarray(self) -> xr.DataArray:
        var_name = self.config.et.external_variable
        if self._dataset is not None:
            if var_name not in self._dataset:
                raise ValueError(
                    f"Variable {var_name!r} not found in external ET dataset: {list(self._dataset.data_vars)}"
                )
            return self._dataset[var_name].copy()

        ext_path_str = self.config.et.external_path
        if not ext_path_str:
            raise ValueError("ExternalETProvider requires 'et.path' (or 'et.external_path')")
        ext_path = Path(ext_path_str)

        if ext_path.suffix.lower() in (".nc", ".nc4", ".netcdf"):
            with xr.open_dataset(ext_path) as opened:
                if var_name not in opened:
                    raise ValueError(
                        f"Variable {var_name!r} not found in {ext_path}. Available: {list(opened.data_vars)}"
                    )
                return opened[var_name].load()

        # Order-independent CSV reader (supports long-table or wide crop columns)
        df = pd.read_csv(ext_path)
        dt_idx, time_cols = _extract_datetime_index(df, ext_path, column_map=self.config.column_map)
        area_col = _find_column(
            df.columns,
            ("subarea", "area_id", "subarea_id", "island"),
            column_map=self.config.column_map,
            context=str(ext_path),
        )
        crop_col = _find_column(
            df.columns,
            ("crop", "crop_code", "landuse"),
            column_map=self.config.column_map,
            required=False,
            context=str(ext_path),
        )
        val_col = _find_column(
            df.columns,
            (var_name, "et_c", "etc", "et", "openet", "et_aw", "etaw", "value"),
            column_map=self.config.column_map,
            required=False,
            context=str(ext_path),
        )
        assert area_col is not None

        work = df.copy()
        work["time"] = dt_idx
        work[area_col] = _coerce_subarea_id(work[area_col])

        if crop_col is not None and val_col is not None:
            # Tidy long table: (time, subarea, crop, value)
            work["_crop_label"] = [
                resolve_crop_behavior(str(c), extra_behaviors=self._extra_behaviors).label
                for c in work[crop_col]
            ]
            piv = (
                work.set_index(["time", area_col, "_crop_label"])[val_col]
                .astype(float)
                .to_xarray()
                .rename({area_col: "subarea", "_crop_label": "crop"})
            )
            return piv

        if val_col is not None and crop_col is None:
            # Subarea-composite table: (time, subarea, value)
            piv_df = work.pivot(index="time", columns=area_col, values=val_col)
            return xr.DataArray(
                piv_df.to_numpy(dtype=np.float64),
                dims=("time", "subarea"),
                coords={
                    "time": piv_df.index.values,
                    "subarea": np.asarray(piv_df.columns.tolist()),
                },
                name=var_name,
            )

        # Wide crop columns: (time, subarea, UR, PA, ..., TO, ...)
        used = {*time_cols, area_col, "time"}
        crop_cols = [c for c in df.columns if c not in used]
        if not crop_cols:
            raise ValueError(
                f"Could not identify ET value or crop columns in external CSV: {ext_path}"
            )
        labels = [
            resolve_crop_behavior(c, extra_behaviors=self._extra_behaviors).label for c in crop_cols
        ]
        subareas = sorted(work[area_col].unique().tolist(), key=lambda x: (isinstance(x, str), x))
        times = pd.DatetimeIndex(sorted(work["time"].unique()))
        arr = np.empty((len(times), len(subareas), len(crop_cols)), dtype=np.float64)
        for idx, col in enumerate(crop_cols):
            p = work.pivot(index="time", columns=area_col, values=col).reindex(
                index=times, columns=subareas
            )
            arr[:, :, idx] = p.to_numpy(dtype=np.float64)
        return xr.DataArray(
            arr,
            dims=("time", "subarea", "crop"),
            coords={
                "time": times.values,
                "subarea": np.asarray(subareas),
                "crop": np.asarray(labels, dtype=object),
            },
            name=var_name,
        )

    def _convert_units_to_mm(self, da: xr.DataArray, is_monthly_total: bool) -> xr.DataArray:
        u = self.config.et.external_units.strip().lower()
        if u in ("mm/d", "mm/day", "mm d-1", "mm day-1", "mm"):
            return da
        if u in ("mm/month", "mm/mon", "mm month-1"):
            if not is_monthly_total:
                raise ValueError("units='mm/month' requires temporal_resolution='monthly'")
            return da  # divided by days_in_month or scaled during disaggregation
        if u in ("in/d", "in/day", "inch/day", "inches/day", "in"):
            return da * 25.4
        if u in ("in/month", "inch/month", "inches/month"):
            if not is_monthly_total:
                raise ValueError("units='in/month' requires temporal_resolution='monthly'")
            return da * 25.4
        raise ValueError(
            f"Unsupported external ET units {self.config.et.external_units!r}; "
            f"expected 'mm/d', 'mm/month', 'in/d', or 'in/month'"
        )

    def compute(
        self,
        time_index: pd.DatetimeIndex,
        ds_meteo_full: xr.Dataset,
        ds_crop_params: xr.Dataset,
        da_year_type: xr.DataArray,
    ) -> xr.Dataset:
        da_okc_full, da_et0_full, da_pcp_full = compute_bare_soil_okc(
            ds_meteo_full,
            beta=self.config.et.beta,
            okc_mode=self.config.et.okc_mode,
            wet_threshold_mode=self.config.et.wet_threshold_mode,
            truncate_2dec=False,
        )
        da_okc = da_okc_full.sel(time=time_index)
        da_et0 = da_et0_full.sel(time=time_index)
        da_pcp = da_pcp_full.sel(time=time_index)

        da_ikc, _ = build_daily_ikc_and_kc(
            time_index=time_index,
            da_okc=da_okc,
            da_et0=da_et0,
            ds_crop_params=ds_crop_params,
            da_year_type=da_year_type,
            pad_year_type=self.config.options.default_pad_year_type,
            season_mode=self.config.et.season_mode,
            type1_initial_kc_mode="legacy",
            beta=self.config.et.beta,
        )

        da_ext = self._load_external_dataarray()
        # Normalize crop labels if 'crop' dim is present
        if "crop" in da_ext.dims:
            mapped_labels = [
                resolve_crop_behavior(str(c), extra_behaviors=self._extra_behaviors).label
                for c in da_ext["crop"].values
            ]
            da_ext = da_ext.assign_coords(crop=mapped_labels).sel(crop=ds_crop_params["crop"])
        else:
            # Broadcast subarea-composite ET across crops
            da_ext = da_ext.expand_dims(crop=ds_crop_params["crop"])

        da_ext = da_ext.sel(subarea=da_et0["subarea"])

        is_monthly = self.config.et.temporal_resolution == "monthly"
        u_lower = self.config.et.external_units.strip().lower()
        is_monthly_total = is_monthly and ("month" in u_lower or "mon" in u_lower)
        da_ext = self._convert_units_to_mm(da_ext, is_monthly_total=is_monthly_total)

        if is_monthly:
            # Map each timestamp in external data and daily time_index to (year, month) key
            ext_times = pd.DatetimeIndex(da_ext["time"].values)
            ext_ym = ext_times.year * 100 + ext_times.month
            da_ext_ym = da_ext.assign_coords(time=ext_ym).rename({"time": "ym"})

            day_ym = xr.DataArray(
                time_index.year * 100 + time_index.month,
                dims=("time",),
                coords={"time": time_index},
                name="ym",
            )
            days_in_mon = xr.DataArray(
                time_index.days_in_month.to_numpy(dtype=np.float64),
                dims=("time",),
                coords={"time": time_index},
            )

            missing_ym = sorted(set(np.unique(day_ym.values)) - set(da_ext_ym["ym"].values))
            if missing_ym:
                raise ValueError(
                    f"Monthly external ET dataset is missing year-months: {missing_ym[:6]}"
                )

            monthly_at_day = da_ext_ym.sel(ym=day_ym).drop_vars("ym")
            # Convert monthly mean rate (mm/d) to monthly total (mm/month) if needed
            monthly_total = monthly_at_day if is_monthly_total else (monthly_at_day * days_in_mon)

            if self.config.et.disaggregation == "et0_scaled":
                et0_ym_sum = da_et0.groupby(day_ym).sum("time").sel(ym=day_ym).drop_vars("ym")
                weight = xr.where(et0_ym_sum > 1e-9, da_et0 / et0_ym_sum, 1.0 / days_in_mon)
                da_daily_mm = monthly_total * weight
            elif self.config.et.disaggregation == "uniform":
                da_daily_mm = monthly_total / days_in_mon
            else:
                raise ValueError(
                    f"Unsupported monthly disaggregation={self.config.et.disaggregation!r}"
                )
        else:
            da_daily_mm = da_ext.sel(time=time_index)

        da_daily_mm = da_daily_mm.transpose("time", "subarea", "crop").clip(min=0.0)

        # Effective implied Kc = ETc / ETo for diagnostics
        da_kc_implied = xr.where(da_et0 > 1e-9, da_daily_mm / da_et0, 0.0).transpose(
            "time", "subarea", "crop"
        )

        is_direct_etaw = self.config.et.external_variable.lower() in ("et_aw", "etaw")
        out_vars: dict[str, xr.DataArray] = {
            "et_c_mm": da_daily_mm.rename("et_c_mm"),
            "ikc": da_ikc,
            "kc": da_kc_implied.rename("kc"),
            "okc": da_okc,
            "et0": da_et0,
            "precip": da_pcp,
        }
        if is_direct_etaw:
            out_vars["external_et_aw_mm"] = da_daily_mm.rename("external_et_aw_mm")

        return xr.Dataset(
            data_vars=out_vars,
            coords={
                "time": time_index,
                "subarea": da_et0["subarea"],
                "crop": ds_crop_params["crop"],
                "crop_code": ds_crop_params["crop_code"],
            },
        )


def create_etc_provider(
    config: DetawConfig,
    extra_behaviors: tuple[CropBehavior, ...] = (),
) -> ETcProvider:
    """Instantiate the configured ``ETcProvider``."""
    p = config.et.provider.lower()
    if p in ("kc_eto", "detaw", "internal", "legacy"):
        return KcETcProvider(config)
    if p in ("openet", "external", "dataset"):
        return ExternalETProvider(config, extra_behaviors=extra_behaviors)
    raise ValueError(
        f"Unknown ET provider {config.et.provider!r}. Supported: 'kc_eto', 'openet' / 'external'"
    )
