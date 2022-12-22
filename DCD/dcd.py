# Delta Channel Depletion (DCD) model version 2.0.0
#
# DCD calculates channel depletion time series based on outputs from
# Delta Evapotranspiration of Applied Water (DETAW).
#
# <license>
#    Copyright (C) State of California, Department of Water Resources.
#
#    DCD v2.0.0 is free software:
#    you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    DCDv v2.0.0 is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with DCD v2.0.0.  If not, see <http://www.gnu.org/licenses>.
# </license>

from pathlib import Path
import os
import argparse
import logging
import yaml
import numpy as np
import pandas as pd
import xarray as xr

# Set up a logger
logging.basicConfig(level=logging.INFO)


def read_subarea_info_file(path_file: str) -> pd.DataFrame:
    """ Read the subarea info text file and return it in a DataFrame

        The file has fixed format and not have headers.

        Parameters
        ----------
        path_file
            file path to read

        Returns
        -------
        pandas.DataFrame
    """
    colnames = ["area", "name", "uplow", "acreage", "docregion"]
    colspecs = ((0, 3), (5, 56), (57, 58), (60, 75), (81, 90))

    df = pd.read_fwf(path_file,
                     colspecs=colspecs,
                     skiprows=6,
                     names=colnames)
    return df


def read_water_year_types(path_file: str) -> pd.DataFrame:
    """ Read water year type from the original text format

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
    df = pd.read_csv(path_file,
                     header=None,
                     delim_whitespace=True,
                     names=["type", "year"])
    return df


def calculate_drained_seepage(df_subareas: pd.DataFrame) -> xr.DataArray:
    """ Calculate drained seepage amount

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
    # Hard-wired areas in the lowland without seepage
    lowland_no_seepage = [133, 134, 135, 136, 137, 140, 141, 142]

    n_areas = len(df_subareas)
    areas = df_subareas["area"].values
    da_drnseep = xr.DataArray(data=np.zeros((n_areas)), dims=["area"],
                              coords=dict(area=areas))

    areas_selected = df_subareas.query(
        "uplow == 1 & docregion == 'Lower'")["area"]
    da_drnseep.loc[areas_selected.values] = 0.013 * \
        df_subareas.loc[areas_selected.index, "acreage"]

    areas_selected = df_subareas.query(
        "uplow == 1 & docregion == 'Midrange'")["area"]
    da_drnseep.loc[areas_selected.values] = 0.074 * \
        df_subareas.loc[areas_selected.index, "acreage"]

    areas_selected = df_subareas.query(
        "uplow == 1 & docregion == 'High'")["area"]
    da_drnseep.loc[areas_selected.values] = 0.095 * \
        df_subareas.loc[areas_selected.index, "acreage"]

    da_drnseep.loc[lowland_no_seepage] = 0.
    return da_drnseep

# XXX This commented function is no longer needed and can be deleted as the gwrates now have to be preprocessed
# def calculate_groundwater_rates(df_subareas, df_gw_rates, dates):
#     """
#     """
#     # Create an empty groundwater array
#     n_dates = len(dates)
#     areas = df_subareas["area"].values
#     n_areas = len(areas)
#     da_gwrates = xr.DataArray(data=np.full((n_dates, n_areas), np.nan),
#                               dims=["time", "area"],
#                               coords=dict(time=dates, area=areas))

#     areas_lower = df_subareas.query("uplow == 1 & docregion == 'Lower'")[
#         'area'].values
#     areas_mid = df_subareas.query("uplow == 1 & docregion == 'Midrange'")[
#         'area'].values
#     areas_high = df_subareas.query("uplow == 1 & docregion == 'High'")[
#         'area'].values

#     for year in range(dates[0].year + 1, dates[-1].year + 1):
#         da_gwrates.loc[f"{year - 1}-10-01":f"{year}-09-30",
#                        :] = df_gw_rates.query("year == @year")["gw_rate"].values
#     da_gwrates.loc[:, areas_lower] = 0.35
#     da_gwrates.loc[:, areas_mid] = 0.30
#     da_gwrates.loc[:, areas_high] = 0.25
#     return da_gwrates


def adjust_leach_water(da_lwa, da_lwd, da_ro):
    """ Adjust (Reduce) leach water with runoff

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
    areas = da_lwa.area.values
    n_areas = len(areas)
    for year in years:
        leach_saved = np.zeros((n_areas, ))
        lwa = da_lwa.loc[f"{year - 1}-10-01":f"{year}-09-30", :].values
        lwd = da_lwd.loc[f"{year - 1}-10-01":f"{year}-09-30", :].values
        ro = da_ro.sel(time=slice(
            f"{year - 1}-10-01", f"{year}-09-30")).values.T
        lwa_adj = np.full_like(lwa, np.nan)
        lwd_adj = np.full_like(lwd, np.nan)
        for r_i in range(lwa.shape[0]):
            month = dates[r_i].month
            n_days = dates[r_i].days_in_month
            leach_saved, lwa_adj[r_i, :], lwd_adj[r_i, :] = \
                np.vectorize(use_runoff_for_leach)(leach_saved,
                                                   lwa[r_i, :],
                                                   lwd[r_i, :],
                                                   ro[r_i, :],
                                                   month,
                                                   n_days)
        da_lwa_adj.loc[f"{year - 1}-10-01":f"{year}-09-30", :] = lwa_adj
        da_lwd_adj.loc[f"{year - 1}-10-01":f"{year}-09-30", :] = lwd_adj

    return da_lwa_adj, da_lwd_adj


def use_runoff_for_leach(leach_saved, lwa, lwd, ro, month, n_days):
    """ Use runoff for leach water (daily)

        This function adjusts applied and drained leach water with runoff
        at each day at each area.
    """
    if lwa > 0.:
        lwa_adj = lwa - ro - leach_saved
        if lwa_adj > 0.:
            leach_saved = 0.
        else:
            leach_saved = -lwa_adj
            lwa_adj = 0.
    else:
        lwa_adj = lwa
    if leach_saved > 0. and lwd > 0.:
        if month == 1:
            lwd_adj = lwd - leach_saved * 0.56 / n_days
        elif month == 2:
            lwd_adj = lwd - leach_saved * 0.29 / n_days
        elif month == 3:
            lwd_adj = lwd - leach_saved * 0.14 / n_days
        else:
            lwd_adj = lwd - leach_saved * 0.01 / n_days
        if lwd_adj < 0.:
            lwd_adj = 0.
    else:
        lwd_adj = 0.
    return leach_saved, lwa_adj, lwd_adj


def take_monthly_average(da):
    df = pd.DataFrame(data=da.values, index=da.date)
    df_mon = df.resample('M').mean()
    da_mon = xr.DataArray(data=df_mon.values,
                          dims=["month", "island"],
                          coords=dict(month=df_mon.index,
                                      island=da.island),
                          attrs=da.attrs)
    da_mon.attrs["name"] = f"{da.attrs['name']}_mon"
    da_mon.attrs["part_e"] = "1MON"
    return da_mon


def calculate_model_depletion(ds_dcd, model_params):
    logging.info("Calculate monthly depletions...")
    # Assuming 25% of the surface runoff goes to the deep percolation
    # See Annual Report 2017, Chapter 3.
    dp_factor = model_params.get("deep_percolation_rate")
    da_percolation = ds_dcd.runoff * dp_factor
    da_percolation.attrs = ds_dcd.runoff.attrs
    da_percolation.attrs["name"] = "DP_island"
    da_percolation.attrs["part_c"] = "DP-FLOW"
    # ds_all[da_percolation.attrs["name"]] = da_percolation
    # Save the data for later uses
    # suffix = f"_{extension}" if extension != "" else ""
    # path_out = f"DCD_islands{suffix}.nc"
    # ds_all.to_netcdf(path_out)

    # Take monthly averages of the NODCU outputs
    # names_to_process = ["drn_wo_ro_island", "div_wo_spg_island",
    #                     "GW_per_island", "RO_island",
    #                     "DP_island", "spg_island"]
    # for name in names_to_process:
    #     da = ds_dcd[name]
    #     da_mon = take_monthly_average(da)
    #     name_mon = f"{name}_mon"
    #     ds_dcd[name_mon] = da_mon

    # FIXME DSM2 only
    extension = ""
    ds_dcd_nodes = island_to_nodes(ds_dcd, model_params, extension)
    return ds_dcd_nodes


def read_distribution_rates(path_file):
    """ Read an island-to-node distribution file and create a distribution
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
    df = pd.read_csv(path_file, delim_whitespace=True, skiprows=2).iloc[:, :3]
    df["I"] = df["I"]
    df["N"] = df["N"]
    # Convert percentage to decimal numbers.
    colname_rate = df.columns[2]
    df[colname_rate] = df[colname_rate] * 0.01
    # Get the list of unique islands
    islands = np.sort(df["I"].unique())
    n_islands = len(islands)
    # Get the list of unique node numbers
    nodes = np.sort(df["N"].unique())
    n_nodes = len(nodes)
    # Create an empty rate matrix
    rates = np.zeros((n_islands, n_nodes))
    # Loop through islands
    for island in islands:
        mask = (df["I"] == island)
        nodes_to_find = df[mask]["N"]
        node_args = np.vectorize(
            lambda x: np.argwhere(nodes == x))(nodes_to_find)
        rates[island - 1, node_args] = df[mask][colname_rate]
    # Create a xarray.DataArray.
    da = xr.DataArray(data=rates, dims=["area", "node"],
                      coords=dict(area=islands,
                                  node=nodes))
    return da


def island_to_nodes(ds_dcd, model_params, extension):
    """ Distribute DCD results to model nodes

        FIXME WIP. Needs to be polished up.

        Returns
        -------
        xarray.Dataset
    """
    # FIXME Some dead code can be removed.
    logging.info("Distribute depletion to nodes...")
    # Read the updated DCD
    # FIXME Do not use hard-wired file names, but OK for now?
    # ext = "" if extension == "" else f"_{ext}.nc"
    # path_dcd_adj = f"DCD_islands_adj.nc"
    # ds_dcd_adj = xr.open_dataset(path_dcd_adj)

    # Read distribution rates
    divratefile = model_params.get("path_dsm2_diversion_rates")
    # Need to drop the last column with NaN
    da_divrate = read_distribution_rates(divratefile)
    drnratefile = model_params.get("path_dsm2_drainage_rates")
    da_drnrate = read_distribution_rates(drnratefile)
    # Merge the two because their locations are not the same, and the
    # original approach has all three even though they do not need to be.
    ds_dist_rates = xr.merge((da_divrate.to_dataset(name='divrate'),
                              da_drnrate.to_dataset(name='drnrate')),
                             fill_value=0.)

    # Diversion
    da_div_nodes = ds_dcd.diversion.dot(ds_dist_rates.divrate)

    # Seepage
    da_spg_nodes = ds_dcd.seepage.dot(ds_dist_rates.divrate)

    # Drainage
    da_drn_nodes = ds_dcd.drainage.dot(ds_dist_rates.drnrate) + \
        ds_dcd.runoff.dot(ds_dist_rates.drnrate)

    # BBID
    # FIXME I think that this can be folded with the code above, but
    #       I am not sure the best way to deal with the BBID island index.
    # FIXME a hard-wired path information. Not good.
    # FIXME Temporaryly commenting out the BBID section
    # path_ext = "../../path_ext.csv"
    # df_path_ext = pd.read_csv(path_ext)
    # extensions = df_path_ext['ext'].unique()
    # dss = []
    # for ext in extensions:
    #     # Read NODCU outputs that are saved for each extension previously
    #     # FIXME Better not to use intermediate files
    #     # FIXME A similar job is done in `split_bbid_new` already.
    #     path_nodcu_ext_nc = f"DCD_islands_{ext}.nc"
    #     ds_dcd_ext = xr.open_dataset(path_nodcu_ext_nc)
    #     df_islands = df_path_ext.query("ext == @ext & part_a == 'BBID'")
    #     if not df_islands.empty:
    #         ds = ds_dcd_ext.sel(island=df_islands["part_b"].values)
    #         dss.append(ds)
    # ds_bbid_sum = xr.concat(dss, dim='island').sum(dim='island')
    # # Add back run off to the drainage
    # # FIXME this is not done completely.
    # ds_bbid_sum['DRAIN-FLOW'] = ds_bbid_sum['drn_wo_ro_island'] + \
    #                             ds_bbid_sum['RO_island']

    # path_dcd_base = f"DCD_islands.nc"
    # ds_dcd_base = xr.open_dataset(path_dcd_base)
    # bbid_islands = df_path_ext[df_path_ext["part_a"] == 'BBID']["part_b"].values
    # ds_dcd_base.sel(island=bbid_islands)

    # Create a xarray.Dataset.
    # NOTE We may come up with better names. For now, use the old name for
    # convenience.
    # FIXME Add attributes.
    ds_dcd_nodes = da_div_nodes.to_dataset(name="DIV-FLOW")
    ds_dcd_nodes["SEEP-FLOW"] = da_spg_nodes
    ds_dcd_nodes["DRAIN-FLOW"] = da_drn_nodes

    return ds_dcd_nodes


def create_argparser() -> argparse.ArgumentParser:
    """ Create an argument parser.

        Returns
        -------
        argepase.ArgumentParser
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("input_yaml", type=str,
                        help="An input YAML file to provide parameters")
    return parser


def read_groundwater_rates(fpath, dates):
    """ Read groundwater rate from a text file

        Parameters
        ----------
        fpath: str
            file name to read
        dates:
            Dates used in the run

        Returns
        -------
        xr.DataArray
            groundwater rates
    """
    start_date = dates[0]
    end_date = dates[-1]
    df_gw_rates = pd.read_csv(fpath,header=0,index_col=0,parse_dates=True)
    df_gw_rates = df_gw_rates.query("time >= @start_date & time <= @end_date")
    n_dates = len(dates)
    areas = np.linspace(1,df_gw_rates.shape[1],df_gw_rates.shape[1],dtype=int)
    n_areas = len(areas)
    da_gwrates = xr.DataArray(data=np.full((n_dates, n_areas), df_gw_rates.values),
                              dims=["time", "area"],
                              coords=dict(time=dates, area=areas))
    return da_gwrates


def calculate_depletion(model_params: dict) -> xr.Dataset:
    """ Calculate channel depletion per area

        Parameters
        ----------
        model_params: dict
           a dict of model parameters (from an input file)

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
    start_date = pd.to_datetime(f'{start_water_year - 1}-10-01')
    end_date = pd.to_datetime(f'{end_water_year}-09-30')
    n_years = end_date.year - start_date.year
    # Create the dates in the date range.
    # This will be coordinates in many arrays later.
    dates = pd.date_range(start=start_date, end=end_date, freq='D')
    # Create months in the modeling period
    months = pd.period_range(start=start_date, end=end_date, freq='M')

    # Read a water year types
    # path_wy_types = model_params.get("path_wateryear_types")
    # df_wy_types = read_water_year_types(path_wy_types)

    # Read a subarea information file
    path_subareas = model_params.get("path_subarea_info")
    df_subareas = read_subarea_info_file(path_subareas)
    # Create an area array.
    # This will be an coordinates in many arrays.
    areas = df_subareas["area"].values
    # n_areas = len(df_subareas)

    # Read irrigation efficiency, \eta
    path_eta = model_params.get("path_irrigation_efficiency")
    df_eta = pd.read_csv(path_eta, header=None)
    da_eta = xr.DataArray(data=df_eta[0].values, dims=[
                          "area"], coords=dict(area=areas))

    # Read groundwater rates
    param_name = "path_groundwater_rates"
    path_groundwater_rates = model_params.get(param_name)
    if path_groundwater_rates is None:
        raise ValueError(f"The input does not include parameter, {param_name}.")
    if not os.path.exists(path_groundwater_rates):
        raise ValueError(f"File {path_groundwater_rates} not found.")
    da_gwrates = read_groundwater_rates(path_groundwater_rates,
                                         dates)

    # Read monthly applied leach water (LW_A) and drained leach water (LW_D)
    path_lwam = model_params.get("path_leach_applied")
    df_lwam = pd.read_csv(path_lwam, delim_whitespace=True,
                          header=None).loc[:, 1:].transpose()
    path_lwdm = model_params.get("path_leach_drained")
    df_lwdm = pd.read_csv(path_lwdm, delim_whitespace=True,
                          header=None).loc[:, 1:].transpose()

    # Read a DETAW NetCDF file
    path_detaw = model_params.get("path_detaw_output")
    ds_detaw = xr.open_dataset(path_detaw).sel(time=slice(start_date, end_date))

    #--------------------------------------------------------
    # Preprocess DETAW outputs
    logging.info("Preprocessing DETAW output...")
    # Aggregating by areas
    da_aw = ds_detaw.et_aw.sum('crop').clip(0.)
    da_seepage = ds_detaw.s_e.sum('crop')
    da_runoff = (ds_detaw.precip - ds_detaw.e_r).sum(
        'crop').rolling(time=5, min_periods=1).mean()
    # da_depletion = ds_detaw.et_c.sum('crop')
    cropname = 'Water Surface'
    da_waterbody = (ds_detaw.et_c.sel(
        crop=cropname) - ds_detaw.precip.sel(crop=cropname)).clip(0.)


    # n_dates = len(dates)
    # Create an array of days in each month for later use
    days_in_month = dates.days_in_month
    da_days_in_month = xr.DataArray(data=days_in_month,
                                    dims=["time"],
                                    coords=dict(time=dates))


    # Months in the ordering in a water year, which is from October
    # to September, 10, 11, 12, 1, ..., 9
    # wy_months = np.roll(np.arange(1, 13), 3)

    # Pre-processing the leach water data
    # FIXME Factor these out
    df_lwam = pd.DataFrame(np.tile(df_lwam.values, (n_years, 1)), index=months)
    df_lwam.rename(columns={i: i+1 for i in df_lwam.columns}, inplace=True)
    df_lwa = df_lwam.resample('1D').ffill()
    da_lwa = xr.DataArray(data=df_lwa.values, dims=["time", "area"],
                          coords=dict(time=dates, area=areas))
    da_lwa /= da_days_in_month

    n_years = end_date.year - start_date.year
    df_lwdm = pd.DataFrame(np.tile(df_lwdm.values, (n_years, 1)), index=months)
    df_lwdm.rename(columns={i: i+1 for i in df_lwdm.columns}, inplace=True)
    df_lwd = df_lwdm.resample('1D').ffill()
    da_lwd = xr.DataArray(data=df_lwd.values, dims=["time", "area"],
                          coords=dict(time=dates, area=areas))
    da_lwd /= da_days_in_month

    # ---------------------------------------------------------------
    # Calculate channel depletion (DCD)

    # Calculate drained seepage (S_D) for each subarea
    da_drnseep = calculate_drained_seepage(df_subareas)
    # Divide it by days in months
    da_drnseep = da_drnseep / da_days_in_month

    # Update seepage (S) by adding drained seepage (S_D)
    # Calculate daily seepage values by dividing by days in months
    da_seepage += da_drnseep

    # Calculate daily groundwater rate per area
    # da_gwrates = calculate_groundwater_rates(df_subareas, df_gw_rates, dates)

    # Calculate groundwater component 1
    da_gw1 = da_gwrates * da_aw / da_eta

    # Calculate groundwater component 2
    da_gw2 = da_gwrates * da_seepage

    # Calculate total groundwater
    da_gwf = da_gw1 + da_gw2

    leach_scale = model_params.get("leach_scale")
    da_lwa *= leach_scale

    runoff_rate = model_params.get("runoff_rate")
    # Update the runoff (RO) with the coefficient
    da_runoff *= runoff_rate

    # Adjust leach water
    da_lwa, da_lwd = adjust_leach_water(da_lwa, da_lwd, da_runoff)

    # Calculate diversion without seepage
    # NOTE Routines to adjust runoff in leach water are not implemented.
    # NOTE FIXME 1977 adjustment of applied water is not implemented.
    da_diversion = da_aw / da_eta - da_gw1 + da_lwa
    if model_params.get("is_adding_waterbody_evaporation"):
        da_diversion += da_waterbody
        da_diversion = da_diversion.drop_vars('crop')

    # Calculate drainage without seepage
    da_drainage = da_aw * (1. - da_eta) / da_eta + da_lwd + \
        da_drnseep

    # Calculate seepage
    da_seepage -= da_gw2

    # -----------------------------------------
    # Write out results
    # Conversion factor.
    taf2cfs = 0.50417

    # Write out data
    # FIXME Add attributes
    da_gwf *= taf2cfs
    ds_dcd = da_gwf.to_dataset(name="groundwater")

    da_runoff *= taf2cfs
    ds_dcd["runoff"] = da_runoff

    da_diversion *= taf2cfs
    ds_dcd["diversion"] = da_diversion

    da_drainage *= taf2cfs
    ds_dcd["drainage"] = da_drainage

    da_seepage *= taf2cfs
    ds_dcd["seepage"] = da_seepage

    return ds_dcd


def main() -> None:
    """ DCD v2 main function
    """
    parser = create_argparser()
    args = parser.parse_args()
    fname_main_yaml = args.input_yaml

    # Read the main yaml input file
    logging.info(f"Reading the main YAML input: {fname_main_yaml}")
    with open(fname_main_yaml, 'r') as file_in:
        model_params = yaml.safe_load(file_in)

    params = model_params.get("dcd")
    ds_dcd = calculate_depletion(params)

    logging.info("Write the depletion to a file...")
    path_dcd = params.get("path_dcd_output")
    ds_dcd.to_netcdf(path_dcd)

    ds_dcd_nodes = calculate_model_depletion(ds_dcd, params)

    logging.info("Write the depletion at nodes to a file...")
    path_dcd_nodes = params.get("path_dcd_node_output")
    ds_dcd_nodes.to_netcdf(path_dcd_nodes)

    logging.info("Finished. Exiting.")


if __name__ == "__main__":
    main()
