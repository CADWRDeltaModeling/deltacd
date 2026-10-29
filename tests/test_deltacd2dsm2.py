"""Unit tests for deltacd2dsm2 conversion tool."""

import numpy as np
import pandas as pd
import pytest
import xarray as xr

pyhecdss = pytest.importorskip("pyhecdss")
from deltacd.utils.deltacd2dsm2 import deltacd2dsm2


def test_deltacd2dsm2_time_shift(tmp_path):
    """Test that deltacd2dsm2 correctly shifts timestamps for PER-AVER DSS output."""
    nc_path = tmp_path / "test_deltacd.nc"
    dss_path = tmp_path / "test_output.dss"
    inp_path = tmp_path / "test_sources_sinks.inp"

    start_date = "2020-01-01"
    n_days = 5
    dates = pd.date_range(start_date, periods=n_days, freq="1D")
    nodes = ["1", "2"]

    div_data = np.array(
        [
            [10.0, 20.0],
            [11.0, 21.0],
            [12.0, 22.0],
            [13.0, 23.0],
            [14.0, 24.0],
        ]
    )
    drain_data = div_data * 0.5
    seep_data = div_data * 0.2

    ds = xr.Dataset(
        data_vars={
            "diversion": (("time", "node"), div_data),
            "drainage": (("time", "node"), drain_data),
            "seepage": (("time", "node"), seep_data),
        },
        coords={"time": dates, "node": nodes},
    )
    ds.to_netcdf(str(nc_path))

    deltacd2dsm2(nc_path, dss_path, inp_path)

    assert dss_path.exists()
    assert inp_path.exists()

    dcd = pyhecdss.DSSFile(str(dss_path))
    cat = dcd.read_catalog()

    # Check catalog records
    assert len(cat) == len(nodes) * 3  # 3 variables per node
    for d_part in cat["D"]:
        # The D-part in DSS catalog should begin in 2020, not shifted backward to 2019
        assert "2020" in d_part
        assert "2019" not in d_part

    # Read back each variable and check timestamps and values
    for node_idx, node in enumerate(nodes):
        cat_node = cat[cat.B == node]
        paths = dcd.get_pathnames(cat_node)
        for path in paths:
            rts = dcd.read_rts(path)
            data_series = rts.data.iloc[:, 0].dropna()

            # Ensure start and end period match original NetCDF dates
            expected_index = pd.period_range(start=start_date, periods=n_days, freq="D")
            pd.testing.assert_index_equal(data_series.index, expected_index)

            # Check values match
            c_part = path.split("/")[3]
            if c_part == "DIV-FLOW":
                np.testing.assert_allclose(data_series.values, div_data[:, node_idx])
            elif c_part == "DRAIN-FLOW":
                np.testing.assert_allclose(data_series.values, drain_data[:, node_idx])
            elif c_part == "SEEP-FLOW":
                np.testing.assert_allclose(data_series.values, seep_data[:, node_idx])

    dcd.close()

    # Verify .inp file content
    inp_content = inp_path.read_text()
    assert "SOURCE_FLOW" in inp_content
    assert "dicu_div_1" in inp_content
    assert "dicu_drain_2" in inp_content
