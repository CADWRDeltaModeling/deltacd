"""Regression tests for ETAW"""
from pathlib import Path
import os
import xarray as xr
import pytest
from pytest_regressions import ndarrays_regression  # noqa: F401
import deltacd.detaw


@pytest.fixture(scope="module")
def suisun_example():
    dirpath = Path(__file__).parent
    os.chdir(dirpath / "testinputs")
    fname_main_yaml = "detaw_test.yaml"
    deltacd.detaw.detaw(fname_main_yaml)


@pytest.fixture
def suisun_example_et0(suisun_example):
    dirpath = Path(__file__).parent
    filepath = dirpath / "testinputs/outputs/ET0.nc"
    return xr.open_dataset(filepath, decode_times=False)


@pytest.fixture
def suisun_example_precip(suisun_example):
    dirpath = Path(__file__).parent
    filepath = dirpath / "testinputs/outputs/precip_suisun.nc"
    return xr.open_dataset(filepath, decode_times=False)


@pytest.fixture
def suisun_example_etawoutput(suisun_example):
    dirpath = Path(__file__).parent
    filepath = dirpath / "testinputs/outputs/detawoutput_suisun_schism.nc"
    return xr.open_dataset(filepath, decode_times=False)


def test_regression_et0(suisun_example_et0, ndarrays_regression):  # noqa: F811
    ndarrays_regression.check(
        {
            "ET0": suisun_example_et0.ET0.values,
            "time": suisun_example_et0.time.values,
        }
    )


def test_regression_precip(suisun_example_precip, ndarrays_regression):  # noqa: F811
    ndarrays_regression.check(
        {
            "area": suisun_example_precip.area.values,
            "precip": suisun_example_precip.precip.values,
            "time": suisun_example_precip.time.values,
        }
    )


def test_regression_etaw(suisun_example_etawoutput, ndarrays_regression):  # noqa: F811
    ndarrays_regression.check(
        {
            "ET_c": suisun_example_etawoutput.et_c.values,
            "ET_aw": suisun_example_etawoutput.et_aw.values,
            "D_sw": suisun_example_etawoutput.d_sw.values,
            "E_r": suisun_example_etawoutput.e_r.values,
            "S_e": suisun_example_etawoutput.s_e.values,
            "precip": suisun_example_etawoutput.precip.values,
            "time": suisun_example_etawoutput.time.values,
            "area_id": suisun_example_etawoutput.area_id.values,
        }
    )
