"""Regression tests for DCD"""
from pathlib import Path
import os
import xarray as xr
import pytest
from pytest_regressions import ndarrays_regression  # noqa: F401
import deltacd.dcd


@pytest.fixture(scope="module")
def suisun_example():
    dirpath = Path(__file__).parent
    os.chdir(dirpath / "testinputs")
    fname_main_yaml = "dcd_test_suisun.yaml"
    deltacd.dcd.dcd(fname_main_yaml)


@pytest.fixture
def suisun_example_dcd_area_output(suisun_example):
    dirpath = Path(__file__).parent
    filepath = dirpath / "testinputs/outputs/dcd_areas_suisun_schism.nc"
    return xr.open_dataset(filepath, decode_times=False)


@pytest.fixture
def suisun_example_dcd_output(suisun_example):
    dirpath = Path(__file__).parent
    filepath = dirpath / "testinputs/outputs/dcd_suisun_schism.nc"
    return xr.open_dataset(filepath, decode_times=False)


def test_regression_dcd_area_output(
    suisun_example_dcd_area_output, ndarrays_regression  # noqa: F811
):
    ndarrays_regression.check(
        {
            "groundwater": suisun_example_dcd_area_output.groundwater.values,
            "runoff": suisun_example_dcd_area_output.runoff.values,
            "diversion": suisun_example_dcd_area_output.diversion.values,
            "drainage": suisun_example_dcd_area_output.drainage.values,
            "seepage": suisun_example_dcd_area_output.seepage.values,
        }
    )


def test_regression_dcd_output(
    suisun_example_dcd_output, ndarrays_regression  # noqa: F811
):
    ndarrays_regression.check(
        {
            "diversion": suisun_example_dcd_output.diversion.values,
            "seepage": suisun_example_dcd_output.seepage.values,
            "drainage": suisun_example_dcd_output.drainage.values,
        }
    )
