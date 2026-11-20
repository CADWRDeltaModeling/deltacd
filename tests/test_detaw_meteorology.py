"""Unit, equation, permutation, and legacy-equivalence tests for DETAW meteorology."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import xarray as xr

from deltacd.detaw.config import DetawConfig
from deltacd.detaw.io import read_subarea_params
from deltacd.detaw.meteorology import (
    DatasetMeteorologyProvider,
    StationMeteorologyProvider,
    asce_penman_monteith,
    hargreaves_samani,
)
from deltacd.detaw_legacy import (
    read_et_correction_factors,
    read_pcp,
    read_temperature,
    weatheroutput,
)

REPO_ROOT = Path(__file__).resolve().parent.parent
TEST_INPUTS = REPO_ROOT / "tests" / "testinputs"


@pytest.mark.parametrize(
    ("yaml_name", "start_str", "end_str"),
    [
        ("detaw_test_delta_dsm2.yaml", "1921-09-30", "1932-09-30"),
        ("detaw_test_suisun_schism.yaml", "1922-09-30", "1932-09-30"),
    ],
)
def test_station_meteorology_equivalence_with_legacy(yaml_name: str, start_str: str, end_str: str):
    """Verify StationMeteorologyProvider matches legacy weatheroutput() to ~1e-14 mm."""
    cfg = DetawConfig.from_yaml(TEST_INPUTS / yaml_name)
    ds_subarea = read_subarea_params(cfg.et_correction_path)
    provider = StationMeteorologyProvider(cfg)
    ds_met = provider.load(cfg.start_water_year, cfg.end_water_year, ds_subarea)

    ts_pcp = read_pcp(start_str, end_str, cfg.meteorology.input_pcp)
    ts_per, eto_corr, _ = read_et_correction_factors(cfg.et_correction_path)
    _, ts_days, tx, tn = read_temperature(start_str, end_str, cfg.meteorology.input_temperature)
    old_pcp, old_et0 = weatheroutput(ts_pcp, ts_per, ts_days, tx, tn, len(ts_days), eto_corr)

    np.testing.assert_allclose(ds_met["et0"].values, old_et0, rtol=1e-13, atol=1e-13)
    np.testing.assert_allclose(ds_met["precip"].values, old_pcp, rtol=1e-13, atol=1e-13)
    assert bool(ds_met.attrs["has_spinup_day"]) is True


def test_negative_precip_is_clipped_unlike_legacy_noop(tmp_path: Path):
    """Verify Category A3 fix: negative sentinel values in precip CSV are clipped to 0."""
    cfg = DetawConfig.from_yaml(TEST_INPUTS / "detaw_test_suisun_schism.yaml")
    df_pcp = pd.read_csv(cfg.meteorology.input_pcp)  # type: ignore[arg-type]
    df_pcp.loc[df_pcp.index[:5], "Brentwood"] = -999.0
    mod_pcp_path = tmp_path / "neg_pcp.csv"
    df_pcp.to_csv(mod_pcp_path, index=False)

    cfg.meteorology.input_pcp = str(mod_pcp_path)
    ds_subarea = read_subarea_params(cfg.et_correction_path)
    ds_met = StationMeteorologyProvider(cfg).load(
        cfg.start_water_year, cfg.end_water_year, ds_subarea
    )
    assert float(ds_met["precip"].min()) >= 0.0


def test_station_meteorology_works_without_sep30_spinup_day(tmp_path: Path):
    """Verify Category B9 improvement: weather starting on Oct 1 runs cleanly without crashing."""
    cfg = DetawConfig.from_yaml(TEST_INPUTS / "detaw_test_suisun_schism.yaml")
    cfg.end_water_year = 1924

    df_pcp = pd.read_csv(cfg.meteorology.input_pcp)  # type: ignore[arg-type]
    # Drop 1922-09-30 spin-up row so data starts on 1922-10-01
    df_pcp_no_spin = df_pcp[~((df_pcp["year"] == 1922) & (df_pcp["month"] == 9))].copy()
    pcp_path = tmp_path / "pcp_oct1.csv"
    df_pcp_no_spin.to_csv(pcp_path, index=False)

    df_tmp = pd.read_csv(cfg.meteorology.input_temperature)  # type: ignore[arg-type]
    df_tmp_no_spin = df_tmp[df_tmp["Date"] != "9/30/1922"].copy()
    tmp_csv_path = tmp_path / "temp_oct1.csv"
    df_tmp_no_spin.to_csv(tmp_csv_path, index=False)

    cfg.meteorology.input_pcp = str(pcp_path)
    cfg.meteorology.input_temperature = str(tmp_csv_path)

    ds_subarea = read_subarea_params(cfg.et_correction_path)
    ds_met = StationMeteorologyProvider(cfg).load(1923, 1924, ds_subarea)
    assert bool(ds_met.attrs["has_spinup_day"]) is False
    assert pd.Timestamp(ds_met["time"].values[0]) == pd.Timestamp("1922-10-01")
    assert pd.Timestamp(ds_met["time"].values[-1]) == pd.Timestamp("1924-09-30")


def test_temp_qc_swap_vs_clamp():
    """Verify temp_qc_mode='clamp' gives 0 ETo when Tmax < Tmin, while 'swap' gives positive ETo."""
    dates = pd.date_range("2020-07-01", periods=2, freq="D")
    ds = xr.Dataset(
        data_vars={
            "tmax": (("time",), np.array([15.0, 34.0])),
            "tmin": (("time",), np.array([30.0, 16.0])),  # Day 0 has Tmax < Tmin
        },
        coords={"time": dates},
    )
    et0_clamp = hargreaves_samani(ds, temp_qc_mode="clamp")
    et0_swap = hargreaves_samani(ds, temp_qc_mode="swap")
    assert float(et0_clamp.isel(time=0)) == 0.0
    assert float(et0_swap.isel(time=0)) > 4.0
    assert np.isclose(float(et0_clamp.isel(time=1)), float(et0_swap.isel(time=1)))


def test_asce_penman_monteith_equations():
    """Verify ASCE Standardized Penman-Monteith implementation (Snyder et al. Appendix B)."""
    dates = pd.date_range("2020-07-15", periods=3, freq="D")
    ds_dew = xr.Dataset(
        data_vars={
            "tmax": (("time",), np.array([33.0, 35.0, 31.0])),
            "tmin": (("time",), np.array([15.0, 16.0, 14.0])),
            "solar_rad": (("time",), np.array([28.0, 29.0, 26.5])),
            "wind_speed": (("time",), np.array([2.2, 2.5, 1.8])),
            "dewpoint": (("time",), np.array([12.0, 13.0, 11.5])),
        },
        coords={"time": dates},
    )
    et0_pm = asce_penman_monteith(ds_dew, latitude_deg=38.12, elevation_m=12.2)
    assert et0_pm.shape == (3,)
    # Typical mid-July Sacramento-San Joaquin Delta ETo is ~6.0-8.0 mm/d
    assert np.all((et0_pm.values > 5.0) & (et0_pm.values < 9.5))

    # Also test RHmax / RHmin vapor-pressure branch (Eq. 10)
    ds_rh = ds_dew.drop_vars("dewpoint").assign(
        rh_max=(("time",), np.array([82.0, 80.0, 84.0])),
        rh_min=(("time",), np.array([25.0, 22.0, 28.0])),
    )
    et0_pm_rh = asce_penman_monteith(ds_rh, latitude_deg=38.12, elevation_m=12.2)
    assert np.all((et0_pm_rh.values > 5.0) & (et0_pm_rh.values < 9.5))


def test_dataset_meteorology_provider_netcdf_and_csv(tmp_path: Path):
    """Verify DatasetMeteorologyProvider loads precomputed NetCDF and CSV with shuffled subareas."""
    cfg = DetawConfig.from_yaml(TEST_INPUTS / "detaw_test_suisun_schism.yaml")
    cfg.end_water_year = 1924
    ds_subarea = read_subarea_params(cfg.et_correction_path)

    ds_station = StationMeteorologyProvider(cfg).load(1923, 1924, ds_subarea)

    # Write NetCDF with reversed subarea order and reload via DatasetMeteorologyProvider
    nc_path = tmp_path / "precomp_meteo.nc"
    ds_station.isel(subarea=slice(None, None, -1)).to_netcdf(nc_path)

    cfg.meteorology.provider = "dataset"
    cfg.meteorology.dataset_path = str(nc_path)
    ds_from_nc = DatasetMeteorologyProvider(cfg).load(1923, 1924, ds_subarea)
    xr.testing.assert_allclose(ds_station, ds_from_nc)

    # Also test CSV round-trip with shuffled columns
    df_long = ds_station.to_dataframe().reset_index()[["subarea", "precip", "et0", "time"]]
    csv_path = tmp_path / "precomp_meteo.csv"
    df_long.to_csv(csv_path, index=False)
    cfg.meteorology.dataset_path = str(csv_path)
    ds_from_csv = DatasetMeteorologyProvider(cfg).load(1923, 1924, ds_subarea)
    xr.testing.assert_allclose(ds_station, ds_from_csv)
