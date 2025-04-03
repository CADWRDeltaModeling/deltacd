"""Regression tests for DCD"""

from pathlib import Path
import os
import xarray as xr
import pytest
from pytest_regressions import ndarrays_regression  # noqa: F401
import deltacd.detaw
import deltacd.dcd

EXAMPLES = ["suisun_schism", "delta_dsm2"]


@pytest.fixture(scope="module", params=EXAMPLES, ids=EXAMPLES)
def example_setup(request):
    """Run the DCD example for tests."""
    name = request.param
    dirpath = Path(__file__).parent
    output_dir = dirpath / "testinputs/outputs"

    # Change to test inputs directory and run DETAW and DCD
    os.chdir(dirpath / "testinputs")
    detaw_yaml_file = f"detaw_test_{name}.yaml"
    deltacd.detaw.detaw(detaw_yaml_file)
    dcd_yaml_file = f"dcd_test_{name}.yaml"
    deltacd.dcd.dcd(dcd_yaml_file)

    # Yield the setup info for use by the tests
    yield name

    # Cleanup: Remove the output_dir
    if output_dir.exists():
        for file in output_dir.iterdir():
            if file.is_file():
                file.unlink()


@pytest.fixture
def example_dcd_area_output(example_setup):
    """Get DCD area output for the example."""
    name = example_setup
    dirpath = Path(__file__).parent
    filepath = dirpath / f"testinputs/outputs/dcd_areas_{name}.nc"
    return name, xr.open_dataset(filepath, decode_times=False)


@pytest.fixture
def example_dcd_output(example_setup):
    """Get DCD output for the example."""
    name = example_setup
    dirpath = Path(__file__).parent
    filepath = dirpath / f"testinputs/outputs/dcd_{name}.nc"
    return name, xr.open_dataset(filepath, decode_times=False)


def test_regression_dcd_area_output(
    example_dcd_area_output, ndarrays_regression
):  # noqa: F811
    """Test the regression of DCD area output."""
    name, data = example_dcd_area_output
    ndarrays_regression.check(
        {
            "groundwater": data.groundwater.values,
            "runoff": data.runoff.values,
            "diversion": data.diversion.values,
            "drainage": data.drainage.values,
            "seepage": data.seepage.values,
        },
        basename=f"test_regression_dcd_areas_{name}.nc",
    )


def test_regression_dcd_output(example_dcd_output, ndarrays_regression):  # noqa: F811
    """Test the regression of DCD output."""
    name, data = example_dcd_output
    ndarrays_regression.check(
        {
            "diversion": data.diversion.values,
            "seepage": data.seepage.values,
            "drainage": data.drainage.values,
        },
        basename=f"test_regression_dcd_{name}.nc",
    )
