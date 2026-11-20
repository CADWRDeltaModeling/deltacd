"""Unit and legacy-equivalence tests for DETAW OKc, IKc, Kc, and ETc."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

from deltacd.detaw.config import DetawConfig
from deltacd.detaw.et import KcETcProvider, build_daily_ikc_and_kc, compute_bare_soil_okc
from deltacd.detaw.io import read_crop_params, read_landuse, read_subarea_params
from deltacd.detaw.meteorology import StationMeteorologyProvider
from deltacd.detaw_legacy import calc_soil_evap

REPO_ROOT = Path(__file__).resolve().parent.parent
TEST_INPUTS = REPO_ROOT / "tests" / "testinputs"


def test_bare_soil_okc_exact_match_when_legacy_truncation_enabled():
    """With truncate_2dec=True, compute_bare_soil_okc matches legacy calc_soil_evap to machine precision."""
    cfg = DetawConfig.from_yaml(TEST_INPUTS / "detaw_test_suisun_schism.yaml")
    ds_sub = read_subarea_params(cfg.et_correction_path)
    ds_met = StationMeteorologyProvider(cfg).load(cfg.start_water_year, cfg.end_water_year, ds_sub)

    da_okc_trunc, _, _ = compute_bare_soil_okc(
        ds_met, beta=2.6, okc_mode="legacy_cumulative", truncate_2dec=True
    )
    da_okc_clean, _, _ = compute_bare_soil_okc(
        ds_met, beta=2.6, okc_mode="legacy_cumulative", truncate_2dec=False
    )

    # Compare subarea 0 against legacy calc_soil_evap
    times = pd.DatetimeIndex(ds_met["time"].values)
    idates = len(times)
    ts_year = times.year.to_numpy(dtype=np.int64)
    ts_days = times.dayofyear.to_numpy(dtype=np.int64)
    start1 = np.array([times[0].year, times[0].month, times[0].day, 23, 0], dtype=np.int32)
    n_years = times[-1].year - times[0].year + 1

    dpy_all = np.zeros(n_years + 2, dtype=np.int64)
    eto_daily = np.zeros((n_years + 2, 367), dtype=np.float64)
    pcp_daily = np.zeros((n_years + 2, 367), dtype=np.float64)
    okc_legacy = np.zeros((n_years + 2, 367), dtype=np.float64)

    calc_soil_evap(
        0,
        2.6,
        idates,
        ts_year,
        ts_days,
        start1,
        dpy_all,
        ds_met["et0"].values.T.copy(),
        ds_met["precip"].values.T.copy(),
        pcp_daily,
        eto_daily,
        okc_legacy,
    )

    # Extract continuous series from okc_legacy[iy, id]
    legacy_series = np.array(
        [okc_legacy[times[i].year - times[0].year + 1, times[i].dayofyear] for i in range(idates)]
    )
    np.testing.assert_allclose(da_okc_trunc.isel(subarea=0).values, legacy_series, atol=1e-14)
    # Without 2-decimal truncation, differences are tiny rounding differences (~O(1e-3))
    assert float(np.mean(np.abs(da_okc_clean.isel(subarea=0).values - legacy_series))) < 0.01


def test_bare_soil_okc_daily_increment_option():
    """Verify Snyder et al. (2006) daily_increment OKc decays faster than cumulative ratio after rain."""
    dates = pd.date_range("2020-06-01", periods=15, freq="D")
    pcp = np.zeros(15)
    pcp[0] = 25.0  # Significant rain on day 0
    ds_met = xr.Dataset(
        data_vars={
            "et0": (("time", "subarea"), np.full((15, 1), 6.0)),
            "precip": (("time", "subarea"), pcp[:, np.newaxis]),
        },
        coords={"time": dates, "subarea": [1]},
    )
    okc_cum, _, _ = compute_bare_soil_okc(ds_met, okc_mode="legacy_cumulative")
    okc_inc, _, _ = compute_bare_soil_okc(ds_met, okc_mode="daily_increment")

    # On day 0 (stage 1), both equal kx = 1.22 - 0.04 * 6.0 = 0.98
    assert np.isclose(float(okc_cum.isel(time=0, subarea=0)), 0.98)
    assert np.isclose(float(okc_inc.isel(time=0, subarea=0)), 0.98)
    # By day 14 (deep stage 2), daily Es/ETo is lower than cumulative CEs/CETo
    assert float(okc_inc.isel(time=14, subarea=0)) < float(okc_cum.isel(time=14, subarea=0))


def test_continuous_vs_legacy_wrapped_grain_seasons():
    """Verify season_mode='continuous' carries fall-planted grain Kc into Jan-May of next year."""
    cfg = DetawConfig.from_yaml(TEST_INPUTS / "detaw_test_delta_dsm2.yaml")
    ds_lu = read_landuse(cfg.landuse_path, water_years=range(1922, 1926))
    ds_cp = read_crop_params(
        cfg.critical_crop_path, cfg.noncritical_crop_path, target_crops=ds_lu["crop"].values
    )
    ds_sub = read_subarea_params(cfg.et_correction_path).isel(subarea=[0])
    ds_met = StationMeteorologyProvider(cfg).load(1922, 1925, ds_sub)
    time_idx = pd.date_range("1921-10-01", "1925-09-30", freq="D")
    da_okc, da_et0, _ = compute_bare_soil_okc(ds_met)

    # In Delta landuse: WY 1923 is BN (non-critical, Kc1=0.45 for GR), WY 1924 is C (critical, Kc1=0.33 for GR).
    # Calendar year 1923 fall planting (Oct 31, 1923) is in WY 1924 (C -> Kc1=0.33),
    # and its spring growth is Jan-May 1924.
    # Meanwhile WY 1923 planting (Oct 31, 1922) is in WY 1923 (BN -> Kc1=0.45),
    # so in 'continuous' mode Jan 1923 continues the BN planting (Kc1=0.45),
    # whereas in 'legacy_wrapped' mode Jan 1924 wraps from 1924's WY 1925 (D -> 0.33).
    ikc_wrap, _ = build_daily_ikc_and_kc(
        time_idx,
        da_okc.sel(time=time_idx),
        da_et0.sel(time=time_idx),
        ds_cp,
        ds_lu["year_type"],
        season_mode="legacy_wrapped",
    )
    ikc_cont, _ = build_daily_ikc_and_kc(
        time_idx,
        da_okc.sel(time=time_idx),
        da_et0.sel(time=time_idx),
        ds_cp,
        ds_lu["year_type"],
        season_mode="continuous",
    )

    # Check WY 1927 (W, non-critical, Kc1=0.45) vs WY 1926 (D, critical, Kc1=0.33):
    # Let's test on synthetic 2-year transition where WY 1923=C and WY 1924=W!
    yt_custom = xr.DataArray(
        ["AN", "C", "W", "W"],
        dims=("water_year",),
        coords={"water_year": [1922, 1923, 1924, 1925]},
    )
    ikc_w2, _ = build_daily_ikc_and_kc(
        time_idx,
        da_okc.sel(time=time_idx),
        da_et0.sel(time=time_idx),
        ds_cp,
        yt_custom,
        season_mode="legacy_wrapped",
    )
    ikc_c2, _ = build_daily_ikc_and_kc(
        time_idx,
        da_okc.sel(time=time_idx),
        da_et0.sel(time=time_idx),
        ds_cp,
        yt_custom,
        season_mode="continuous",
    )
    # Planting on Oct 31, 1922 (in WY 1923 = 'C' -> Kc1=0.33) extends to Nov 1922..Jan 1923:
    # On 1923-01-02 (still in initial stage A..B of Oct 31 planting, length 211 * 0.20 = 42 days -> Dec 12,
    # wait: on 1923-01-02 it's in rapid growth B..C, which starts from 0.33 vs 0.45!)
    val_cont_jan1923 = float(ikc_c2.sel(time="1923-01-02", crop="Irrig grain"))
    val_wrap_jan1923 = float(ikc_w2.sel(time="1923-01-02", crop="Irrig grain"))
    # In legacy_wrapped, cal year 1923 uses yearType[y] = WY 1924 ('W' -> higher Kc1=0.45)!
    # In continuous, Jan 1923 continues the Oct 1922 planting in WY 1923 ('C' -> lower Kc1=0.33)!
    assert val_cont_jan1923 < val_wrap_jan1923


def test_kc_etc_provider_shapes_and_open_water():
    """Verify KcETcProvider produces non-negative ETc and Kc=1.1 for Water surface."""
    cfg = DetawConfig.from_yaml(TEST_INPUTS / "detaw_test_suisun_schism.yaml")
    cfg.end_water_year = 1924
    ds_lu = read_landuse(cfg.landuse_path, water_years=[1923, 1924])
    ds_cp = read_crop_params(
        cfg.critical_crop_path, cfg.noncritical_crop_path, target_crops=ds_lu["crop"].values
    )
    ds_sub = read_subarea_params(cfg.et_correction_path)
    ds_met = StationMeteorologyProvider(cfg).load(1923, 1924, ds_sub)
    time_idx = pd.date_range("1922-10-01", "1924-09-30", freq="D")

    ds_et = KcETcProvider(cfg).compute(time_idx, ds_met, ds_cp, ds_lu["year_type"])
    assert ds_et["et_c_mm"].dims == ("time", "subarea", "crop")
    assert np.allclose(ds_et["kc"].sel(crop="Water surface").values, 1.1)
    assert float(ds_et["et_c_mm"].min()) >= 0.0
