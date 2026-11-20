"""Unit tests for ExternalETProvider (OpenET daily/monthly NetCDF & CSV integration)."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

from deltacd.detaw.config import DetawConfig
from deltacd.detaw.et import ExternalETProvider, KcETcProvider, create_etc_provider
from deltacd.detaw.io import read_crop_params, read_landuse, read_subarea_params
from deltacd.detaw.meteorology import StationMeteorologyProvider

REPO_ROOT = Path(__file__).resolve().parent.parent
TEST_INPUTS = REPO_ROOT / "tests" / "testinputs"


def _setup_two_year_suisun():
    cfg = DetawConfig.from_yaml(TEST_INPUTS / "detaw_test_suisun_schism.yaml")
    cfg.end_water_year = 1924
    ds_lu = read_landuse(cfg.landuse_path, water_years=[1923, 1924])
    ds_cp = read_crop_params(
        cfg.critical_crop_path, cfg.noncritical_crop_path, target_crops=ds_lu["crop"].values
    )
    ds_sub = read_subarea_params(cfg.et_correction_path).isel(subarea=[0, 1, 2])
    ds_met = StationMeteorologyProvider(cfg).load(1923, 1924, ds_sub)
    time_idx = pd.date_range("1922-10-01", "1924-09-30", freq="D")
    return cfg, ds_lu, ds_cp, ds_sub, ds_met, time_idx


def test_openet_daily_netcdf_exact_roundtrip(tmp_path: Path):
    """Supplying daily ETc via NetCDF (with permuted crops/subareas) reproduces ETc identically."""
    cfg, ds_lu, ds_cp, ds_sub, ds_met, time_idx = _setup_two_year_suisun()
    ds_internal = KcETcProvider(cfg).compute(time_idx, ds_met, ds_cp, ds_lu["year_type"])

    openet_nc = tmp_path / "openet_daily.nc"
    # Reverse crop and subarea coordinate order to prove order independence
    ds_to_save = xr.Dataset(
        {
            "et_c": ds_internal["et_c_mm"].isel(
                crop=slice(None, None, -1), subarea=slice(None, None, -1)
            )
        }
    )
    ds_to_save.to_netcdf(openet_nc)

    cfg.et.provider = "openet"
    cfg.et.external_path = str(openet_nc)
    cfg.et.external_variable = "et_c"
    cfg.et.external_units = "mm/d"
    cfg.et.temporal_resolution = "daily"

    provider = create_etc_provider(cfg)
    assert isinstance(provider, ExternalETProvider)
    ds_ext = provider.compute(time_idx, ds_met, ds_cp, ds_lu["year_type"])

    xr.testing.assert_allclose(ds_ext["et_c_mm"], ds_internal["et_c_mm"])


def test_openet_monthly_csv_disaggregation_conserves_monthly_volume(tmp_path: Path):
    """Monthly OpenET totals (mm/month) in shuffled long CSV conserve exact monthly totals when disaggregated."""
    cfg, ds_lu, ds_cp, ds_sub, ds_met, time_idx = _setup_two_year_suisun()

    months = pd.date_range("1922-10-01", "1924-09-01", freq="MS")
    records = []
    for m_ts in months:
        for sa in ds_sub["subarea"].values:
            # Subarea-composite OpenET monthly total (mm/month)
            seasonal_mm = 180.0 if m_ts.month in (6, 7, 8) else 45.0
            records.append(
                {
                    "ET_OpenET": seasonal_mm + float(sa),
                    "area_id": int(sa),
                    "Date": m_ts.strftime("%Y-%m-%d"),
                }
            )
    df_mon = pd.DataFrame(records).sample(frac=1.0, random_state=21)
    csv_path = tmp_path / "openet_monthly.csv"
    df_mon.to_csv(csv_path, index=False)

    cfg.et.provider = "openet"
    cfg.et.external_path = str(csv_path)
    cfg.et.external_variable = "ET_OpenET"
    cfg.et.external_units = "mm/month"
    cfg.et.temporal_resolution = "monthly"
    cfg.et.disaggregation = "et0_scaled"

    ds_ext = ExternalETProvider(cfg).compute(time_idx, ds_met, ds_cp, ds_lu["year_type"])

    # Sum daily et_c_mm back to monthly totals and check exact conservation
    july_1923_sum = float(
        ds_ext["et_c_mm"]
        .sel(time=slice("1923-07-01", "1923-07-31"), subarea=1, crop="Tomato")
        .sum()
    )
    assert np.isclose(july_1923_sum, 181.0, atol=1e-9)

    # Also verify daily pattern tracks daily ETo variations within July 1923
    july_et0 = ds_ext["et0"].sel(time=slice("1923-07-01", "1923-07-31"), subarea=1).values
    july_etc = (
        ds_ext["et_c_mm"]
        .sel(time=slice("1923-07-01", "1923-07-31"), subarea=1, crop="Tomato")
        .values
    )
    corr = np.corrcoef(july_et0, july_etc)[0, 1]
    assert corr > 0.999999
