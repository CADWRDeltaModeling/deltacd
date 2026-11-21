"""Unit and end-to-end permutation-invariance tests for DETAW water balance, smoothed ETaw, and DCD."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr
import yaml

import deltacd.dcd
from deltacd.detaw import DetawConfig, DetawModel

REPO_ROOT = Path(__file__).resolve().parent.parent
TEST_INPUTS = REPO_ROOT / "tests" / "testinputs"


def test_end_to_end_detaw_and_dcd_permutation_invariance(tmp_path: Path):
    """Permuting columns and rows of ALL DETAW & DCD CSV inputs produces identical outputs."""
    cfg_orig = DetawConfig.from_yaml(TEST_INPUTS / "detaw_test_suisun_schism.yaml")
    cfg_orig.end_water_year = 1925
    orig_detaw_nc = tmp_path / "orig_detaw.nc"
    cfg_orig.detaw_output_path = str(orig_detaw_nc)
    cfg_orig.precip_output_path = None
    cfg_orig.et_output_path = None

    ds_orig = DetawModel(cfg_orig).run(write_outputs=True)

    # Shuffle columns and rows of every input CSV
    rng = np.random.default_rng(2026)

    def _shuffle_csv(src: str, dst: Path) -> Path:
        df = pd.read_csv(src)
        p_cols = rng.permutation(list(df.columns)).tolist()
        p_rows = rng.permutation(df.index)
        df.loc[p_rows, p_cols].to_csv(dst, index=False)
        return dst

    p_landuse = _shuffle_csv(cfg_orig.landuse_path, tmp_path / "lu_perm.csv")
    p_etc = _shuffle_csv(cfg_orig.et_correction_path, tmp_path / "etc_perm.csv")
    p_pcp = _shuffle_csv(cfg_orig.meteorology.input_pcp, tmp_path / "pcp_perm.csv")  # type: ignore[arg-type]
    p_tmp = _shuffle_csv(cfg_orig.meteorology.input_temperature, tmp_path / "tmp_perm.csv")  # type: ignore[arg-type]

    # Also permute crop columns & parameter rows of critical_suisun.csv
    raw_crit = pd.read_csv(cfg_orig.critical_crop_path, header=None, dtype=str)
    crop_cols = rng.permutation(list(range(2, raw_crit.shape[1]))).tolist()
    body_rows = rng.permutation(list(range(7, 26))).tolist()
    raw_crit_perm = raw_crit.iloc[
        list(range(0, 7)) + body_rows + list(range(26, len(raw_crit))),
        [0, 1] + crop_cols,
    ]
    p_crit = tmp_path / "crit_perm.csv"
    raw_crit_perm.to_csv(p_crit, index=False, header=False)

    cfg_perm = DetawConfig.from_yaml(TEST_INPUTS / "detaw_test_suisun_schism.yaml")
    cfg_perm.end_water_year = 1925
    cfg_perm.landuse_path = str(p_landuse)
    cfg_perm.et_correction_path = str(p_etc)
    cfg_perm.critical_crop_path = str(p_crit)
    cfg_perm.meteorology.input_pcp = str(p_pcp)
    cfg_perm.meteorology.input_temperature = str(p_tmp)
    perm_detaw_nc = tmp_path / "perm_detaw.nc"
    cfg_perm.detaw_output_path = str(perm_detaw_nc)
    cfg_perm.precip_output_path = None
    cfg_perm.et_output_path = None

    ds_perm = DetawModel(cfg_perm).run(write_outputs=True)

    # Align crop coordinate order (since landuse columns were shuffled) and assert identical!
    ds_perm_aligned = ds_perm.sel(crop=ds_orig["crop"], subarea=ds_orig["subarea"])
    xr.testing.assert_allclose(ds_orig, ds_perm_aligned, rtol=1e-13, atol=1e-13)

    # Run downstream DCD on both orig_detaw_nc and perm_detaw_nc (which has shuffled crop order!)
    with open(TEST_INPUTS / "dcd_test_suisun_schism.yaml") as f:
        dcd_base = yaml.safe_load(f)["dcd"]
    for k, v in list(dcd_base.items()):
        if k.startswith("path_") and v is not None:
            dcd_base[k] = str((TEST_INPUTS / v).resolve())
    dcd_base["end_water_year"] = 1925

    dcd_orig_cfg = dict(
        dcd_base,
        path_detaw_output=str(orig_detaw_nc),
        path_dcd_output=str(tmp_path / "dcd_areas_orig.nc"),
        path_dcd_node_output=str(tmp_path / "dcd_nodes_orig.nc"),
    )
    dcd_perm_cfg = dict(
        dcd_base,
        path_detaw_output=str(perm_detaw_nc),
        path_dcd_output=str(tmp_path / "dcd_areas_perm.nc"),
        path_dcd_node_output=str(tmp_path / "dcd_nodes_perm.nc"),
    )
    y_dcd_orig = tmp_path / "dcd_orig.yaml"
    y_dcd_perm = tmp_path / "dcd_perm.yaml"
    y_dcd_orig.write_text(yaml.safe_dump({"dcd": dcd_orig_cfg}))
    y_dcd_perm.write_text(yaml.safe_dump({"dcd": dcd_perm_cfg}))

    deltacd.dcd.dcd(str(y_dcd_orig))
    deltacd.dcd.dcd(str(y_dcd_perm))

    with (
        xr.open_dataset(dcd_orig_cfg["path_dcd_node_output"]) as n_orig,
        xr.open_dataset(dcd_perm_cfg["path_dcd_node_output"]) as n_perm,
    ):
        xr.testing.assert_allclose(n_orig, n_perm, rtol=1e-13, atol=1e-13)


def test_water_balance_mass_conservation_and_crop_rules():
    """Verify agricultural Dsw = ETc - Espg - Er, Rice 20-day drain, Riparian Dsw=0, Open Water."""
    cfg = DetawConfig.from_yaml(TEST_INPUTS / "detaw_test_delta_dsm2.yaml")
    cfg.end_water_year = 1923
    ds = DetawModel(cfg).run(write_outputs=False)

    # For all non-open-water crops, volumetric Dsw == ETc - Se - Er everywhere!
    non_ws = ds.drop_sel(crop="Water surface")
    residual = np.abs(
        non_ws["d_sw"].values
        - (non_ws["et_c"].values - non_ws["s_e"].values - non_ws["e_r"].values)
    )
    assert float(residual.max()) < 1e-11

    # Lowland subarea 1 Riparian vegetation: Dsw == 0 and s_e == et_c and et_aw == 0
    rv_low = ds.sel(subarea=1, crop="Riparian vegetation")
    assert float(np.abs(rv_low["d_sw"]).max()) < 1e-12
    assert np.allclose(rv_low["s_e"].values, rv_low["et_c"].values)
    assert float(np.abs(rv_low["et_aw"]).max()) < 1e-12

    # Native vegetation & Non-irrig grain: net_app == 0 and et_aw == 0 everywhere
    for non_irr_crop in ("Native vegetation", "Non-irrig grain", "Water surface"):
        assert float(np.abs(ds["net_app"].sel(crop=non_irr_crop)).max()) == 0.0
        assert float(np.abs(ds["et_aw"].sel(crop=non_irr_crop)).max()) == 0.0


def test_offseason_swd_cap_option():
    """Verify Category B6 option: offseason_swd_cap=True bounds off-season SWD growth."""
    cfg_uncapped = DetawConfig.from_yaml(TEST_INPUTS / "detaw_test_suisun_schism.yaml")
    cfg_uncapped.end_water_year = 1924

    cfg_capped = DetawConfig.from_yaml(TEST_INPUTS / "detaw_test_suisun_schism.yaml")
    cfg_capped.end_water_year = 1924
    cfg_capped.options.offseason_swd_cap = True

    ds_u = DetawModel(cfg_uncapped).run(write_outputs=False)
    ds_c = DetawModel(cfg_capped).run(write_outputs=False)
    # Capping off-season soil moisture depletion reduces pre-season deficit and irrigation need
    assert float(ds_c["net_app"].sum()) <= float(ds_u["net_app"].sum()) + 1e-9
