# DeltaCD, Delta channel depletion model
# version 0.9.0
#
# dcd.py
# This module calculates channel depletion time series based on outputs from
# Delta Evapotranspiration of Applied Water (DETAW).
# The code is based on the previous DCD v1.3.
#
# See the LICENSE file for the license of this software.

from pathlib import Path
import os
import argparse
import logging
import yaml
import json
import numpy as np
import pandas as pd
import xarray as xr
from deltacd.detaw import convert_to_absolute_path

# Set up a logger
logging.basicConfig(level=logging.INFO)


def read_water_year_types(path_file: str) -> pd.DataFrame:
    """Read water year type from the original text format

    The file has two columns. The first column is one-letter water year
    type, and the second for years.

    Parameters
    ----------
    path_file: str
        file name to read

    Returns
    -------
    pandas.DataFrame
        a DataFrame that has two columns, type and year, for water year
        types.
    """
    df = pd.read_csv(
        path_file, header=None, delim_whitespace=True, names=["type", "year"]
    )
    return df


def calculate_drained_seepage(df_subareas: pd.DataFrame) -> xr.DataArray:
    """Calculate drained seepage amount

    This function calculate drained seepage amount using subarea
    information. The result is drained seepage per area.

    Parameters
    ----------
    df_subareas: pandas.DataFrame
        subarea information

    Returns
    -------
    xarray.DataArray
        drained seepage per area

    """
    # FIXME Hard-wired areas in the lowland without seepage
    lowland_no_seepage = ["133", "134", "135", "136", "137", "140", "141", "142"]

    n_areas = len(df_subareas)
    areas = df_subareas["area_id"].values
    da_drnseep = xr.DataArray(
        data=np.zeros((n_areas)), dims=["area_id"], coords=dict(area_id=areas)
    )

    areas_selected = df_subareas.query("uplow == 'low' & docregion == 'Lower'")[
        "area_id"
    ]
    da_drnseep.loc[areas_selected.values] = (
        0.013 * df_subareas.loc[areas_selected.index, "acreage"]
    )

    areas_selected = df_subareas.query("uplow == 'low' & docregion == 'Midrange'")[
        "area_id"
    ]
    da_drnseep.loc[areas_selected.values] = (
        0.074 * df_subareas.loc[areas_selected.index, "acreage"]
    )

    areas_selected = df_subareas.query("uplow == 'low' & docregion == 'High'")[
        "area_id"
    ]
    da_drnseep.loc[areas_selected.values] = (
        0.095 * df_subareas.loc[areas_selected.index, "acreage"]
    )
    # FIXME Performing check so we only use for delta and not Suisun.
    # This can be removed when lowland_no_seegage is not hardwired.
    if n_areas > 15:
        da_drnseep.loc[lowland_no_seepage] = 0.0

    return da_drnseep


def calculate_groundwater_rates(df_gw_rates, dates):
    """Calculates daily groundwater rates based on yearly values

    Parameters
    ----------
    df_gw_rates: dataframe
        yearly groundwater rates for each area
    dates:
        dates used in the model run

    Returns
    -------
    xr.DataArray
        daily groundwater rates for model dates and nareas
    """
    # Assume water year
    df_gw_rates_daily = df_gw_rates.resample("D").ffill().shift(-92, freq="D")
    df_gw_rates_daily.index.rename("time", inplace=True)
    da_gwrates = xr.DataArray(
        data=df_gw_rates_daily.values,
        dims=["time", "area_id"],
        coords=dict(
            time=df_gw_rates_daily.index.to_timestamp(),
            area_id=df_gw_rates_daily.columns,
        ),
    )
    return da_gwrates.sel(time=dates)


def adjust_leach_water(da_lwa, da_lwd, da_ro):
    """Adjust (Reduce) leach water with runoff

    When there is runoff from precipitation, applied leach water can be
    reduced by the amount of the runoff. Even more, the leftover of the
    runoff can be applied as leach water for following days, if available.
    This logic is adopted from the previous Fortran version of DCD
    with some modification.

    Parameters
    ----------
    da_lwa: xarray.DataArray
        DataArray of daily applied leach water of all areas
    da_lwd: xarray.DataArray
        DataArray of daily drained leach water of all areas
    da_ro: xarray.DataArray
        DataArray of daily runoff of all areas

    Returns
    -------
    xarray.DataArray, xarray.DataArray
       DataArray of adjusted applied and drained leach water
    """
    da_lwa_adj = da_lwa.copy(deep=True)
    da_lwd_adj = da_lwd.copy(deep=True)
    dates = pd.date_range(da_lwa.time.values[0], da_lwa.time.values[-1])
    years = range(dates[0].year + 1, dates[-1].year + 1)
    areas = da_lwa.area_id.values
    n_areas = len(areas)
    for year in years:
        leach_saved = np.zeros((n_areas,))
        lwa = da_lwa.loc[f"{year - 1}-10-01":f"{year}-09-30", :].values
        lwd = da_lwd.loc[f"{year - 1}-10-01":f"{year}-09-30", :].values
        ro = da_ro.sel(time=slice(f"{year - 1}-10-01", f"{year}-09-30")).values.T
        lwa_adj = np.full_like(lwa, np.nan)
        lwd_adj = np.full_like(lwd, np.nan)
        for r_i in range(lwa.shape[0]):
            month = dates[r_i].month
            n_days = dates[r_i].days_in_month
            leach_saved, lwa_adj[r_i, :], lwd_adj[r_i, :] = np.vectorize(
                use_runoff_for_leach
            )(leach_saved, lwa[r_i, :], lwd[r_i, :], ro[r_i, :], month, n_days)
        da_lwa_adj.loc[f"{year - 1}-10-01":f"{year}-09-30", :] = lwa_adj
        da_lwd_adj.loc[f"{year - 1}-10-01":f"{year}-09-30", :] = lwd_adj

    return da_lwa_adj, da_lwd_adj


def use_runoff_for_leach(
    leach_saved_in: float, lwa: float, lwd: float, ro: float, month: int, n_days: int
):
    """Use runoff for leach water (daily)

    This function adjusts applied and drained leach water with runoff
    at each day at each area.
    """
    if lwa > 0.0:
        lwa_adj = lwa - ro - leach_saved_in
        if lwa_adj > 0.0:
            leach_saved = 0.0
        else:
            leach_saved = leach_saved_in + ro - lwa
            lwa_adj = 0.0
    else:
        leach_saved = leach_saved_in
        lwa_adj = 0.0
    if leach_saved > 0.0 and lwd > 0.0:
        # January
        if month == 1:
            lwd_adj = lwd - leach_saved * 0.56 / n_days
        elif month == 2:
            lwd_adj = lwd - leach_saved * 0.29 / n_days
        elif month == 3:
            lwd_adj = lwd - leach_saved * 0.14 / n_days
        else:
            lwd_adj = lwd - leach_saved * 0.01 / n_days
        if lwd_adj < 0.0:
            lwd_adj = 0.0
    else:
        lwd_adj = 0.0
    return leach_saved, lwa_adj, lwd_adj


def take_monthly_average(da: xr.DataArray) -> xr.DataArray:
    """Take monthly average of daily data

    Parameters
    ----------
    da: xarray.DataArray
        daily data

    Returns
    -------
    xarray.DataArray
        monthly average data
    """
    # FIXME Does not need to convert back and forth
    df = pd.DataFrame(data=da.values, index=da.date)
    df_mon = df.resample("ME").mean()
    da_mon = xr.DataArray(
        data=df_mon.values,
        dims=["month", "island"],
        coords=dict(month=df_mon.index, island=da.island),
        attrs=da.attrs,
    )
    da_mon.attrs["name"] = f"{da.attrs['name']}_mon"
    da_mon.attrs["part_e"] = "1MON"
    return da_mon


def calculate_model_depletion(ds_dcd, model_params, input_data):
    logging.info("Calculate monthly depletions...")
    ds_dcd_nodes = island_to_nodes(ds_dcd, input_data)
    ds_dcd_nodes = set_dcd_node_global_attributes(ds_dcd_nodes, model_params)
    return ds_dcd_nodes


def set_dcd_node_global_attributes(
    ds_dcd_nodes: xr.Dataset, model_params: dict
) -> xr.Dataset:
    """Set global attributes for the DCD node dataset

    Parameters
    ----------
    ds_dcd_nodes: xarray.Dataset
        dataset containing DCD node data
    """
    title = "DeltaCD outputs"
    leach_factor = model_params.get("leach_scale")
    water_surface_evaporation = model_params.get("is_adding_waterbody_evaporation")
    dcd_inputs = json.dumps(model_params)
    ds_dcd_nodes = ds_dcd_nodes.assign_attrs(
        title=title,
        leach_factor=leach_factor,
        water_surface_evaporation=str(water_surface_evaporation),
        dcd_inputs=dcd_inputs,
    )
    return ds_dcd_nodes


def read_distribution_ratios(path_file, df_split_table) -> xr.DataArray:
    """Read an island-to-node distribution file and create a distribution
    matrix.

    This function read a text file containing list of island-to-node
    distribution values and create a xarray.DataArray with a dense
    matrix size of (# of islands, # of distributed nodes).
    The coordinates of the DataArray are one-based indices of the
    islands and nodes.

    Parameters
    ----------
    path_file
        the input file name

    Returns
    -------
    xarray.DataArray
    """
    # Read the file, and drop the last column that is filled by
    # read_csv automatically.
    df_ratios = pd.read_csv(path_file, dtype={"area_id": str, "node": str})
    if df_split_table is not None:
        df_split_table_to_one = df_split_table.query(
            "subregion_id not in @df_ratios.node"
        )
        if not df_split_table_to_one.empty:
            df_ratios_w_one = pd.DataFrame(
                {
                    "area_id": df_split_table_to_one["area_id"],
                    "node": df_split_table_to_one["subregion_id"],
                    "factor": 1.0,
                }
            )
            df_ratios = pd.concat([df_ratios, df_ratios_w_one], ignore_index=True)

    # Get the list of unique islands
    # NOTE Need to make sure that this sorted indices work.
    area_ids = df_ratios["area_id"].unique()
    n_areas = len(area_ids)
    # Get the list of unique node numbers
    nodes = df_ratios["node"].unique()
    n_nodes = len(nodes)
    # Create an empty rate matrix
    rates = np.zeros((n_areas, n_nodes))
    # Loop through islands
    for ind, island in enumerate(area_ids):
        mask = df_ratios["area_id"] == island
        nodes_to_find = df_ratios[mask]["node"]
        node_args = np.vectorize(lambda x: np.argwhere(nodes == x))(nodes_to_find)
        rates[ind, node_args] = df_ratios[mask]["factor"]
    # Create a xarray.DataArray.
    da = xr.DataArray(
        data=rates, dims=["area_id", "node"], coords=dict(area_id=area_ids, node=nodes)
    )
    return da


def island_to_nodes(ds_dcd, input_data):
    """Distribute DCD results to model nodes

    Returns
    -------
    xarray.Dataset
    """
    logging.info("Distribute depletion to nodes...")
    # Read the updated DCD
    da_divrate = input_data.get("divrate")
    da_drnrate = input_data.get("drnrate")

    # Merge the two because their locations are not the same, and the
    # original approach has all three even though they do not need to be.
    ds_dist_rates = xr.merge(
        (da_divrate.to_dataset(name="divrate"), da_drnrate.to_dataset(name="drnrate")),
        fill_value=0.0,
    )

    # Diversion
    da_div_nodes = ds_dcd.diversion.dot(ds_dist_rates.divrate)
    da_div_nodes = da_div_nodes.assign_attrs(
        unit="cfs", description="diversion flow at nodes"
    )

    # Seepage
    da_spg_nodes = ds_dcd.seepage.dot(ds_dist_rates.divrate)
    da_spg_nodes = da_spg_nodes.assign_attrs(
        unit="cfs", description="seepage flow at nodes"
    )

    # Drainage
    da_drn_nodes = ds_dcd.drainage.dot(ds_dist_rates.drnrate) + ds_dcd.runoff.dot(
        ds_dist_rates.drnrate
    )
    da_drn_nodes = da_drn_nodes.assign_attrs(
        unit="cfs", description="drainage flow at nodes"
    )

    # Create a xarray.Dataset.
    # NOTE We may come up with better names. For now, use the old name for
    # convenience.
    # FIXME Add attributes.
    ds_dcd_nodes = da_div_nodes.to_dataset(name="diversion")
    ds_dcd_nodes["seepage"] = da_spg_nodes
    ds_dcd_nodes["drainage"] = da_drn_nodes

    return ds_dcd_nodes


def create_argparser() -> argparse.ArgumentParser:
    """Create an argument parser.

    Returns
    -------
    argparse.ArgumentParser
    """
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "input_yaml", type=str, help="An input YAML file to provide parameters"
    )
    return parser


def read_groundwater_rates(fpath: str) -> pd.DataFrame:
    """Read groundwater rate from a text file~

    Parameters
    ----------
    fpath: str
        file name to read

    Returns
    -------
    pandas.DataFrame
        groundwater rates
    """
    # The first column is the year.
    df_gw_rates = pd.read_csv(fpath, header=0, index_col=0, parse_dates=[0])
    df_gw_rates.index = df_gw_rates.index.to_period("Y")
    return df_gw_rates


def add_split_area_data(df_data, df_split_table, use_zeros=False):
    """Add data for split areas by copying data from the source areas

    Parameters
    ----------
    df_data: pandas.DataFrame
        data to split
    df_split_table: pandas.DataFrame
        split table
    use_zeros: bool, optional
        if True, set the values of the split areas to zero
        (default: False)

    Returns
    -------
    pandas.DataFrame
        augmented data with the split area information
    """
    df_data_updated = df_data.copy()
    # Need to convert area_id from int to str.
    # Select rows that have area_id_from in the base area_ids
    areas_to_add = df_split_table["area_id_from"].to_numpy()
    df_data_to_add = (
        df_data_updated.set_index("area_id").loc[areas_to_add].reset_index().copy()
    )
    if use_zeros:
        df_data_to_add.loc[:, df_data.columns[1:]] = 0.0
    df_data_to_add["area_id"] = df_split_table["area_id"]
    df_data_updated = pd.concat([df_data_updated, df_data_to_add], ignore_index=True)

    return df_data_updated


def calculate_depletion(model_params: dict, input_data: dict) -> xr.Dataset:
    """Calculate channel depletion per area

    Parameters
    ----------
    model_params: dict
        a dict of model parameters (from an input file)
    input_data: dict
        a dict of input data

    Returns
    -------
    xarray.DataSet
        a dataset containing daily time series per area, diversion,
        drainage, seepage, and runoff
    """
    logging.info("Reading input files ")

    # Get the calculation period
    start_water_year = model_params.get("start_water_year")
    end_water_year = model_params.get("end_water_year")
    start_date = pd.to_datetime(f"{start_water_year - 1}-10-01")
    end_date = pd.to_datetime(f"{end_water_year}-09-30")
    n_years = end_date.year - start_date.year
    # Create the dates in the date range.
    # This will be coordinates in many arrays later.
    dates = pd.date_range(start=start_date, end=end_date, freq="D")
    # Create months in the modeling period
    months = pd.period_range(start=start_date, end=end_date, freq="M")

    df_subareas = input_data.get("subarea_info")
    areas = df_subareas["area_id"].values

    # Read a DETAW NetCDF file
    path_detaw = model_params.get("path_detaw_output")
    ds_detaw = xr.open_dataset(path_detaw).sel(time=slice(start_date, end_date))
    # Change the data type of the area coordinates
    ds_detaw = ds_detaw.assign_coords(area_id=ds_detaw.area_id.astype(str))
    ds_detaw["et_aw"] = ds_detaw["et_aw"].clip(0.0)

    # --------------------------------------------------------
    # Preprocess DETAW outputs
    logging.info("Preprocessing DETAW output...")

    # Split the DETAW outputs
    df_split_table = input_data.get("split_table")
    # FIXME Crop names hardwired.
    crops = [
        "Urban",
        "Irrig pasture",
        "Alfalfa",
        "All field",
        "Sugar beets",
        "Irrig grain",
        "Rice",
        "Truck crops",
        "Tomato",
        "Orchard",
        "Vineyard",
        "Riparian vegetation",
        "Native vegetation",
        "Non-irrig grain",
        "Water surface",
        "Duck pond",
    ]
    # df_split_table_first = df_split_table.query("area_id_from in @df_subareas.area_id")
    # df_split_table_second = df_split_table.query("area_id_from in @df_split_table.area_id_to")

    ds_detaw_original = ds_detaw.copy(deep=True)
    if df_split_table is not None:
        df_landuse_split_ratios = input_data.get("landuse_split_ratios")
        col_ratios = df_landuse_split_ratios.columns[
            ~df_landuse_split_ratios.columns.isin(["area_id_from", "area_id_to"])
        ]
        n_crops = len(col_ratios)
        # Note that the code assumes the ordering of the crops in different datasets
        # This assumption is not safe.
        list_merged = []
        for varname, da in ds_detaw.data_vars.items():
            list_da = []
            for _, row in df_landuse_split_ratios.iterrows():
                area_id_from = row["area_id_from"]
                da_ratios = xr.DataArray(
                    data=row[col_ratios].astype(float),
                    dims=["crop"],
                    coords=dict(crop=crops[:n_crops]),
                )
                da_split = (
                    ds_detaw_original[varname].sel(area_id=area_id_from) * da_ratios
                )
                list_da.append(da_split)
            da_new = xr.concat(list_da, dim="area_id").assign_coords(
                area_id=df_split_table["area_id"].to_numpy()
            )
            da_merged = xr.concat([da, da_new], dim="area_id")
            list_merged.append(da_merged)
        ds_detaw_split = xr.merge(list_merged)
    else:
        ds_detaw_split = ds_detaw_original

    # Aggregating by areas
    # Need to omit the water surface components
    da_aw = ds_detaw_split.et_aw.drop_sel(crop="Water surface").clip(0.0).sum("crop")
    da_seepage = ds_detaw_split.s_e.clip(0.0).sum("crop")
    da_runoff = (
        (ds_detaw_split.precip - ds_detaw_split.e_r)
        .clip(0.0)
        .sum("crop")
        .rolling(time=5, min_periods=1)
        .mean()
    )
    cropname = "Water surface"
    da_depletion_waterbody = (
        ds_detaw_split.et_c.sel(crop=cropname)
        - ds_detaw_split.precip.sel(crop=cropname)
    ).clip(0.0)

    # Create an array of days in each month for later use
    days_in_month = dates.days_in_month
    da_days_in_month = xr.DataArray(
        data=days_in_month, dims=["time"], coords=dict(time=dates)
    )

    # Pre-processing the leach water data
    # FIXME Factor these out
    df_lwa = input_data.get("lwa")
    months_wateryear = [str((i + 9) % 12 + 1) for i in range(12)]
    df_lwa = pd.DataFrame(
        np.tile(
            df_lwa.set_index("area_id")[months_wateryear].to_numpy().transpose(),
            (n_years, 1),
        ),
        index=months,
    )
    df_lwa = df_lwa.resample("1D").ffill()
    da_lwa = xr.DataArray(
        data=df_lwa.values,
        dims=["time", "area_id"],
        coords=dict(time=dates, area_id=areas),
    )
    da_lwa /= da_days_in_month

    df_lwd = input_data.get("lwd")
    df_lwd = pd.DataFrame(
        np.tile(
            df_lwd.set_index("area_id")[months_wateryear].to_numpy().transpose(),
            (n_years, 1),
        ),
        index=months,
    )
    df_lwd = df_lwd.resample("1D").ffill()
    da_lwd = xr.DataArray(
        data=df_lwd.values,
        dims=["time", "area_id"],
        coords=dict(time=dates, area_id=areas),
    )
    da_lwd /= da_days_in_month

    # ---------------------------------------------------------------
    # Calculate channel depletion (DCD)

    # Calculate drained seepage (S_D) for each subarea
    # Simply multiplying by the area with coffecients to get monthly values
    # NOTE This can be improved by daily rates, not monthly rates.
    # NOTE This can be set per subareas as well.
    da_drnseep = calculate_drained_seepage(df_subareas)
    # Divide it by days in months
    da_drnseep = da_drnseep / da_days_in_month

    # Update seepage (S) by adding drained seepage (S_D)
    # S = S_E + S_D
    da_seepage += da_drnseep

    # Calculate daily groundwater rate per ar
    df_gw_rates = input_data.get("gw_rates")
    da_gwrates = calculate_groundwater_rates(df_gw_rates, dates)

    # Calculate groundwater component 1
    df_eta = input_data.get("eta")
    da_eta = xr.DataArray(
        data=df_eta["irrigation_efficiency"].values,
        dims=["area_id"],
        coords=dict(area_id=areas),
    )
    da_gw1 = da_gwrates * da_aw / da_eta

    # Calculate groundwater component 2
    # How much of seepage comes from groundwater
    # GW2 = GW_rate * S
    da_gw2 = da_gwrates * da_seepage

    # Calculate total groundwater
    # How much of water comes from groundwater
    da_gwf = da_gw1 + da_gw2

    # Apply leach water scale
    leach_scale = model_params.get("leach_scale")
    da_lwa *= leach_scale

    runoff_rate = model_params.get("runoff_rate")
    # Update the runoff (RO) with the coefficient
    da_runoff_after_percolation = da_runoff * runoff_rate

    # Adjust leach water with available runoff
    da_lwa, da_lwd = adjust_leach_water(da_lwa, da_lwd, da_runoff_after_percolation)

    # Calculate diversion without seepage
    # NOTE Routines to adjust runoff in leach water are not implemented.
    # NOTE FIXME 1977 adjustment of applied water is not implemented.
    da_diversion = da_aw / da_eta - da_gw1 + da_lwa
    if model_params.get("is_adding_waterbody_evaporation"):
        da_diversion += da_depletion_waterbody
        da_diversion = da_diversion.drop_vars("crop")

    # Calculate drainage without seepage
    da_drainage = da_aw * (1.0 - da_eta) / da_eta + da_lwd + da_drnseep

    # Calculate seepage, reducing seepage by the amount of the groundwater
    # S = S - GW2
    da_seepage -= da_gw2

    # -----------------------------------------
    # Write out results
    # Conversion factor.
    taf2cfs = 0.50417

    # Write out data
    # FIXME Add attributes
    da_gwf *= taf2cfs
    ds_dcd = da_gwf.to_dataset(name="groundwater")

    ds_dcd["groundwater_to_applied_water"] = da_gw1 * taf2cfs

    ds_dcd["applied_water"] = da_aw / da_eta * taf2cfs

    da_runoff_after_percolation *= taf2cfs
    ds_dcd["runoff"] = da_runoff_after_percolation

    # Calculate deep percolation
    factor_deep_percolation = model_params.get("deep_percolation_rate")
    ds_deep_percolation = da_runoff * factor_deep_percolation
    ds_dcd["deep_percolation"] = ds_deep_percolation

    da_diversion *= taf2cfs
    ds_dcd["diversion"] = da_diversion

    da_drainage *= taf2cfs
    ds_dcd["drainage"] = da_drainage

    da_seepage *= taf2cfs
    ds_dcd["seepage"] = da_seepage

    # Now we subtract the values from 'extension'
    ds_dcd_subtracted = ds_dcd.copy(deep=True)
    if df_split_table is not None:
        for _, row in df_split_table.iterrows():
            area_id_from = row["area_id_from"]
            area_id = row["area_id"]
            ds_dcd_subtracted.loc[dict(area_id=area_id_from)] -= ds_dcd.sel(
                area_id=area_id
            )

    return ds_dcd_subtracted


def convert_relpath_to_abspath(params: dict, dir_input_base: Path) -> dict:
    """Convert relative paths to absolute paths in a dict

    The function looks for items in the dict that start with "path_", check
    if it is a relative path, and then convert it to an absolute path.
    Note that this function modifies the input dict.

    Parameters
    ----------
    params: dict
        a dictionary of parameters
    dir_input_base: Path
        a base directory for relative paths

    Returns
    -------
    None
    """
    for key, value in params.items():
        if key.startswith("path_"):
            params[key] = convert_to_absolute_path(value, dir_input_base)


def read_dcd_input_data(params: dict) -> dict:
    """Read input data for DCD

    This functions read input data for DCD from files and split the information
    when necessary.

    Parameters
    ----------
    params: dict
        a dictionary of parameters

    Returns
    -------
    dict
        a dictionary of input data
    """
    input_data = {}
    # Read a split table
    path_split_table = params.get("path_split_table")
    df_split_table = None
    if path_split_table is not None:
        df_split_table = pd.read_csv(
            path_split_table,
            dtype={"area_id_from": str, "area_id_to": str, "subregion_id": str},
        )
        # Create and add new area ids for split-out areas
        df_split_table["area_id"] = (
            df_split_table["area_id_from"]
            + "_"
            + df_split_table["area_id_to"]
            + "_"
            + df_split_table["subregion_id"]
        )
        input_data["split_table"] = df_split_table

    # Read a subarea information file
    path_subareas = params.get("path_subarea_info")
    df_subareas = pd.read_csv(path_subareas, dtype={"area_id": str})
    if df_split_table is not None:
        df_subareas_updated = df_subareas.copy()
        path_subareas_split = params.get("path_subarea_info_split")
        df_subareas_split = pd.read_csv(
            path_subareas_split, dtype={"area_id_from": str, "area_id_to": str}
        )

        # Modify the split area information before concatenating
        df_subareas_split.drop(columns=["area_id_to"], inplace=True)
        df_subareas_split.rename(columns={"area_id_from": "area_id"}, inplace=True)
        # Select only the areas that are in the original subarea file
        df_subareas_split = df_subareas_split.query(
            "area_id in @df_subareas.area_id"
        ).copy()
        df_subareas_split["area_id"] = df_split_table.query(
            "area_id_from in @df_subareas.area_id"
        )["area_id"]

        df_subareas_updated = pd.concat(
            [df_subareas_updated, df_subareas_split], ignore_index=True
        )
        input_data["subarea_info"] = df_subareas_updated
    else:
        input_data["subarea_info"] = df_subareas

    if df_split_table is not None:
        path_landuse_split = params.get("path_landuse_split_ratios")
        df_landuse_split_ratios = pd.read_csv(
            path_landuse_split, dtype={"area_id_from": str, "area_id_to": str}
        )
        col_crops = df_landuse_split_ratios.columns[
            ~df_landuse_split_ratios.columns.isin(["area_id_from", "area_id_to"])
        ]
        df_landuse_split_ratios = df_landuse_split_ratios.astype(
            {c: float for c in col_crops}
        )
        df_landuse_split_ratios = df_landuse_split_ratios.query(
            "area_id_from in @df_subareas.area_id"
        ).copy()
        input_data["landuse_split_ratios"] = df_landuse_split_ratios

    # Read irrigation efficiency, \eta
    path_eta = params.get("path_irrigation_efficiency")
    df_eta = pd.read_csv(path_eta, dtype={"area_id": str})
    if df_split_table is not None:
        df_eta = add_split_area_data(df_eta, df_split_table)
    input_data["eta"] = df_eta

    # Read groundwater rates
    param_name = "path_groundwater_rates"
    path_groundwater_rates = params.get(param_name)
    if path_groundwater_rates is None:
        raise ValueError(f"The input does not include parameter, {param_name}.")
    if not os.path.exists(path_groundwater_rates):
        raise ValueError(f"File {path_groundwater_rates} not found.")
    df_gw_rates = read_groundwater_rates(path_groundwater_rates)
    if df_split_table is not None:
        df_split_table_to = df_split_table.query("area_id_from in @df_subareas.area_id")
        df_gw_rates_split = df_gw_rates.loc[:, df_split_table_to["area_id_from"]]
        df_gw_rates_split = df_gw_rates_split.set_axis(
            df_split_table_to["area_id"], axis=1
        )
        df_gw_rates = pd.concat([df_gw_rates, df_gw_rates_split], axis=1)
    input_data["gw_rates"] = df_gw_rates

    # Read monthly applied leach water (LW_A) and drained leach water (LW_D)
    path_lwa = params.get("path_leach_applied")
    df_lwa = pd.read_csv(path_lwa, dtype={"area_id": str})
    if df_split_table is not None:
        df_lwa = add_split_area_data(df_lwa, df_split_table, use_zeros=True)
    input_data["lwa"] = df_lwa

    path_lwd = params.get("path_leach_drained")
    df_lwd = pd.read_csv(path_lwd, dtype={"area_id": str})
    if df_split_table is not None:
        df_lwd = add_split_area_data(df_lwd, df_split_table, use_zeros=True)
    input_data["lwd"] = df_lwd

    # Read distribution rates
    divratefile = params.get("path_dsm2_diversion_rates")
    # Need to drop the last column with NaN
    da_divrate = read_distribution_ratios(divratefile, df_split_table)
    input_data["divrate"] = da_divrate

    drnratefile = params.get("path_dsm2_drainage_rates")
    da_drnrate = read_distribution_ratios(drnratefile, df_split_table)
    input_data["drnrate"] = da_drnrate

    return input_data


def dcd(fname_main_yaml: str) -> None:
    """Run DCD with an input yaml file

    Parameters
    ----------
    fname_main_yaml: str
        path to the main yaml input file
    """
    # Read the main yaml input file
    logging.info(f"Reading the main YAML input: {fname_main_yaml}")
    with open(fname_main_yaml, "r") as file_in:
        model_params = yaml.safe_load(file_in)
    # Get the direcotry of the main input file. This will be used to resolve
    # other relative path information.
    dir_input_base = Path(fname_main_yaml).resolve().parent
    model_params["dir_input_base"] = dir_input_base

    params = model_params.get("dcd")
    convert_relpath_to_abspath(params, dir_input_base)

    input_data = read_dcd_input_data(params)
    ds_dcd = calculate_depletion(params, input_data)

    logging.info("Write the depletion to a file...")
    path_dcd = params.get("path_dcd_output")
    head_tail = os.path.split(path_dcd)
    if not os.path.exists(head_tail[0]):
        os.mkdir(head_tail[0])
    ds_dcd.to_netcdf(path_dcd)

    ds_dcd_nodes = calculate_model_depletion(ds_dcd, params, input_data)

    logging.info("Write the depletion at nodes to a file...")
    path_dcd_nodes = params.get("path_dcd_node_output")
    head_tail = os.path.split(path_dcd_nodes)
    if not os.path.exists(head_tail[0]):
        os.mkdir(head_tail[0])
    ds_dcd_nodes.to_netcdf(path_dcd_nodes)

    logging.info("Finished. Exiting.")


def main() -> None:
    """DCD v2 main function"""
    parser = create_argparser()
    args = parser.parse_args()
    fname_main_yaml = args.input_yaml

    dcd(fname_main_yaml)


if __name__ == "__main__":
    main()
