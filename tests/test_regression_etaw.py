"""Regression tests for ETAW"""

from pathlib import Path
import os
import xarray as xr
import pytest
from pytest_regressions import ndarrays_regression  # noqa: F401
import deltacd.detaw

EXAMPLES = ["suisun_schism", "delta_dsm2"]


@pytest.fixture(scope="module", params=EXAMPLES, ids=EXAMPLES)
def example_setup(request):
    """Parametrized fixture for running different test example configurations."""
    name = request.param
    dirpath = Path(__file__).parent
    output_dir = dirpath / "testinputs/outputs"

    # Change to test inputs directory and run DETAW
    os.chdir(dirpath / "testinputs")
    detaw_yaml_file = f"detaw_test_{name}.yaml"
    deltacd.detaw.detaw(detaw_yaml_file)

    # Yield the setup info for use by the tests
    yield name

    # Cleanup: Remove the output_dir
    if output_dir.exists():
        for file in output_dir.iterdir():
            if file.is_file():
                file.unlink()


@pytest.fixture
def example_et0(example_setup):
    """Parametrized fixture for ET0 data"""
    name = example_setup
    dirpath = Path(__file__).parent

    filepath = dirpath / f"testinputs/outputs/et0_{name}.nc"

    return name, xr.open_dataset(filepath, decode_times=False)


@pytest.fixture
def example_precip(example_setup):
    """Parametrized fixture for precipitation data"""
    name = example_setup
    dirpath = Path(__file__).parent
    filepath = dirpath / f"testinputs/outputs/precip_{name}.nc"
    return name, xr.open_dataset(filepath, decode_times=False)


@pytest.fixture
def example_etawoutput(example_setup):
    """Parametrized fixture for ETAW output data"""
    name = example_setup
    dirpath = Path(__file__).parent
    filepath = dirpath / f"testinputs/outputs/detawoutput_{name}.nc"
    return name, xr.open_dataset(filepath, decode_times=False)


def test_regression_et0(example_et0, ndarrays_regression):  # noqa: F811
    """Test ET0 regression for all examples"""
    name, data = example_et0
    ndarrays_regression.check(
        {
            "time": data.time.values,
            "ET0": data.ET0.values,
        },
        basename=f"test_regression_et0_{name}",
    )


def test_regression_precip(example_precip, ndarrays_regression):  # noqa: F811
    """Test precipitation regression for all examples"""
    name, data = example_precip
    ndarrays_regression.check(
        {
            "time": data.time.values,
            "subarea": data.subarea.values,
            "precip": data.precip.values,
        },
        basename=f"test_regression_precip_{name}",
    )


def test_regression_etaw(example_etawoutput, ndarrays_regression):  # noqa: F811
    """Test ETAW regression for all examples"""
    name, data = example_etawoutput
    ndarrays_regression.check(
        {
            "ET_c": data.et_c.values,
            "ET_aw": data.et_aw.values,
            "D_sw": data.d_sw.values,
            "E_r": data.e_r.values,
            "S_e": data.s_e.values,
            "precip": data.precip.values,
            "time": data.time.values,
            "subarea": data.subarea.values,
        },
        basename=f"test_regression_etaw_{name}",
    )
