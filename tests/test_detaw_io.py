"""Unit and permutation-invariance tests for DETAW config, params, and IO readers."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

from deltacd.detaw.config import DetawConfig
from deltacd.detaw.io import (
    read_crop_params,
    read_landuse,
    read_station_precipitation,
    read_subarea_params,
    read_weather_station_csv,
)

REPO_ROOT = Path(__file__).resolve().parent.parent
TEST_INPUTS = REPO_ROOT / "tests" / "testinputs"


def test_config_loads_legacy_yaml():
    cfg = DetawConfig.from_yaml(TEST_INPUTS / "detaw_test_delta_dsm2.yaml")
    assert cfg.start_water_year == 1922
    assert cfg.end_water_year == 1932
    assert cfg.model_start_year == 1921
    assert cfg.daily_output is True
    assert cfg.daily_output_unit == 1
    assert Path(cfg.landuse_path).exists()
    assert Path(cfg.meteorology.input_pcp).exists()  # type: ignore[arg-type]


def test_read_landuse_extracts_true_water_year_type_in_delta():
    """Verify Category A5 fix: Delta landuse CSV reads TYPE ('AN','C','D'), not DATE ('192')."""
    ds_lu = read_landuse(
        TEST_INPUTS / "inputs" / "landuse_dsm2.csv",
        water_years=range(1922, 1933),
    )
    assert ds_lu.sizes["water_year"] == 11
    assert ds_lu.sizes["subarea"] == 168
    assert ds_lu.sizes["crop"] == 15
    wy_types = dict(zip(ds_lu["water_year"].values.tolist(), ds_lu["year_type"].values.tolist()))
    assert wy_types[1922] == "AN"
    assert wy_types[1924] == "C"
    assert wy_types[1925] == "D"
    assert wy_types[1927] == "W"


def test_read_landuse_column_and_row_permutation_invariance(tmp_path: Path):
    """Shuffling columns (including area_id, DATE, TYPE, crops) and rows produces identical xarray values."""
    src_csv = TEST_INPUTS / "inputs" / "landuse_suisun_schism.csv"
    ds_orig = read_landuse(src_csv, water_years=range(1923, 1933))

    df = pd.read_csv(src_csv)
    rng = np.random.default_rng(42)
    shuffled_cols = rng.permutation(list(df.columns)).tolist()
    shuffled_df = df.loc[rng.permutation(df.index), shuffled_cols]

    perm_csv = tmp_path / "landuse_shuffled.csv"
    shuffled_df.to_csv(perm_csv, index=False)

    ds_perm = read_landuse(perm_csv, water_years=range(1923, 1933))
    # Align crop order to compare data identically
    ds_perm_aligned = ds_perm.sel(crop=ds_orig["crop"])
    xr.testing.assert_identical(ds_orig, ds_perm_aligned)


def test_read_crop_params_row_and_column_permutation_invariance(tmp_path: Path):
    """Legacy crop_info CSV reader resolves rows by label and crops by code/alias in any order."""
    crit_csv = TEST_INPUTS / "inputs" / "critical_suisun.csv"
    noncrit_csv = TEST_INPUTS / "inputs" / "noncritical_suisun.csv"

    ds_orig = read_crop_params(crit_csv, noncrit_csv)
    assert ds_orig.sizes["crop"] == 16
    assert "Duck pond" in ds_orig["crop"].values
    assert "All field" in ds_orig["crop"].values

    # Permute crop columns and parameter rows in the raw CSV
    raw = pd.read_csv(crit_csv, header=None, dtype=str)
    crop_col_indices = list(range(2, raw.shape[1]))
    rng = np.random.default_rng(7)
    perm_crop_cols = rng.permutation(crop_col_indices).tolist()
    raw_col_perm = raw.iloc[:, [0, 1] + perm_crop_cols].copy()

    # Shuffle parameter rows below the 'Code' row (rows 7..25)
    body_indices = list(range(7, 26))
    perm_body = rng.permutation(body_indices).tolist()
    row_order = list(range(0, 7)) + perm_body + list(range(26, len(raw)))
    raw_both_perm = raw_col_perm.iloc[row_order, :]

    perm_crit_csv = tmp_path / "crit_perm.csv"
    raw_both_perm.to_csv(perm_crit_csv, index=False, header=False)

    ds_perm = read_crop_params(
        perm_crit_csv,
        noncrit_csv,
        target_crops=ds_orig["crop"].values.tolist(),
    )
    xr.testing.assert_identical(ds_orig, ds_perm)


def test_read_tidy_crop_params_format(tmp_path: Path):
    """Verify read_crop_params also accepts a clean row-per-crop tidy CSV."""
    crit_csv = TEST_INPUTS / "inputs" / "crop_info_critical.csv"
    noncrit_csv = TEST_INPUTS / "inputs" / "crop_info_noncritical.csv"
    ds_legacy = read_crop_params(crit_csv, noncrit_csv)

    # Export critical regime to a tidy CSV with shuffled columns and read it back
    tidy_records = []
    for crop_idx, code in enumerate(ds_legacy["crop_code"].values):
        rec = {"crop_code": code}
        for var in ds_legacy.data_vars:
            if "regime" in ds_legacy[var].dims:
                rec[var] = ds_legacy[var].sel(regime="critical").isel(crop=crop_idx).item()
        tidy_records.append(rec)
    tidy_df = pd.DataFrame(tidy_records)
    tidy_path = tmp_path / "tidy_critical.csv"
    tidy_df.sample(frac=1.0, random_state=11)[tidy_df.columns[::-1]].to_csv(tidy_path, index=False)

    ds_from_tidy = read_crop_params(
        tidy_path,
        noncrit_csv,
        target_crops=ds_legacy["crop"].values.tolist(),
    )
    xr.testing.assert_allclose(ds_legacy, ds_from_tidy)


def test_station_matching_and_subarea_permutation_invariance(tmp_path: Path):
    """Station weights and precip match by normalized station aliases regardless of column/row order."""
    pcp_csv = TEST_INPUTS / "inputs" / "mm_pcp.csv"
    etc_csv = TEST_INPUTS / "inputs" / "et_correction_dsm2.csv"

    da_pcp = read_station_precipitation(pcp_csv, "1921-10-01", "1922-09-30")
    ds_sub = read_subarea_params(etc_csv)

    # Both should share identical canonical station sets despite 'RioVista' vs 'Rio Vista'
    # and 'Tracy' vs 'Tracy_Carbona'
    assert set(da_pcp["station"].values) == set(ds_sub["station"].values)

    # Shuffle station columns and subarea rows in both files
    df_pcp = pd.read_csv(pcp_csv)
    df_etc = pd.read_csv(etc_csv)
    rng = np.random.default_rng(99)

    pcp_perm_path = tmp_path / "pcp_perm.csv"
    df_pcp.loc[rng.permutation(df_pcp.index), rng.permutation(list(df_pcp.columns))].to_csv(
        pcp_perm_path, index=False
    )

    etc_perm_path = tmp_path / "etc_perm.csv"
    df_etc.loc[rng.permutation(df_etc.index), rng.permutation(list(df_etc.columns))].to_csv(
        etc_perm_path, index=False
    )

    da_pcp_perm = read_station_precipitation(pcp_perm_path, "1921-10-01", "1922-09-30")
    ds_sub_perm = read_subarea_params(etc_perm_path, subareas=ds_sub["subarea"].values)

    weighted_orig = xr.dot(da_pcp, ds_sub["precip_weights"], dim="station")
    weighted_perm = xr.dot(da_pcp_perm, ds_sub_perm["precip_weights"], dim="station")
    xr.testing.assert_allclose(weighted_orig, weighted_perm)


def test_read_subarea_params_handles_header_typo():
    """Verify dev/et_correction_calsim3.csv ('ETo Corection Factor') parses without error."""
    typo_csv = REPO_ROOT / "dev" / "et_correction_calsim3.csv"
    if typo_csv.exists():
        ds = read_subarea_params(typo_csv)
        assert "eto_correction" in ds
        assert ds.sizes["subarea"] > 0


def test_read_weather_station_column_order_invariance(tmp_path: Path):
    """Weather CSV reader parses Date and Tx/Tn columns in any column order."""
    src_csv = TEST_INPUTS / "inputs" / "LODI_PT.csv"
    ds_orig = read_weather_station_csv(src_csv, "1921-09-30", "1922-09-30")
    assert set(ds_orig.data_vars) >= {"tmax", "tmin", "precip"}
    assert ds_orig.sizes["time"] == 366  # includes Sep 30 spin-up day

    df = pd.read_csv(src_csv)
    perm_path = tmp_path / "lodi_perm.csv"
    df[df.columns[::-1]].iloc[::-1].to_csv(perm_path, index=False)

    ds_perm = read_weather_station_csv(perm_path, "1921-09-30", "1922-09-30")
    xr.testing.assert_identical(ds_orig, ds_perm)
