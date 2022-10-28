#!/usr/bin/env python
# coding: utf-8
# A script to calculate Delta Channel Depletion from the DETAW outputs.
#
# The script is a Python implementation of the former NODCU Fortran codes in
# DETAW-DCD.
#
# Note that this version does not include parts that are not being used in
# DSM2, SCHISM, and CalSim3, such as concentration.
# Also note that some details are not converted for simplicity. For example,
# routines to adjust the amount of leach water by the runoff is not implemented
# yet.
# FIXME Input file names are hardwired and assumed in the current directory.

import logging
import numpy as np
import pandas as pd
import xarray as xr

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


def sniff_n_years_from_detaw_output(path_file: str) -> int:
    """ Read a DETAW text output file and find out how many years of
        data it has.

        How many years of data in a DETAW text output file is necessary
        to read the whole file. This function sniffs it out.

        Parameters
        ----------
        path_file: str
            file name to read

        Returns
        -------
        int
            the number of years in the data

    """
    with open(path_file) as fh:
        # Skip the first two
        fh.readline()
        fh.readline()
        # Test up to 200 lines
        for line_no in range(200):
            line = fh.readline()
            n_items = len(line.strip().split())
            if n_items < 355:
                break
        return line_no


def read_detaw_output_to_dataframe(path_file: str) -> pd.DataFrame:
    """ Read a DETAW output file in text such as DICU5.* into
        pandas.DataFrame.

        The outputs from the current DETAW are in text.
        Each island (subarea) section has two-line header and following data.
        Each data line contains series of daily values for a water year,
        for example,
        1A 1922 35.3 36.0 34.0 ...
        If a year is a leap year, the number of data is 366 plus two (for
        the area name and the year.)
        This pattern repeats for areas.

        This function removes letter parts from the first column.
        The last number of non-leap years will be NaN.

        Parameters
        ----------
        path_file
            file path to read

        Returns
        -------
        pandas.DataFrame
            DETAW output data
    """
    # Find out how many years of data the file has.
    n_years = sniff_n_years_from_detaw_output(path_file)
    # Read data into a DataFrame first
    colnames = ["area", "year"] + list(range(1, 367))
    df = pd.read_csv(path_file,
                     header=None,
                     skiprows=lambda x: x % (n_years + 2) < 2,
                     names=colnames,
                     delim_whitespace=True)
    return df


def read_detaw_output_to_dataarray(path_file):
    """ Read a DETAW output file in text such as DICU5.* into
        xarray.DataAarray.

        This function reads a DICU5.* file into a DataFrame first. The
        DataFrame is converted into a xarray.DataArray.
        The DataArray contains a two-dimensional array, dates in rows, and
        areas in columns. The column names are integer area numbers.

        NOTE: This may be slow due to rearranging data structure. This
        will not be necessary once we move away from the text outputs from
        DETAW.

        Parameters
        ----------
        path_file
            file path to read

        Returns
        -------
        pandas.DataFrame
            DETAW output data
    """
    # Read the data into a DataFrame first.
    df = read_detaw_output_to_dataframe(path_file)

    # Remove a trailing character from area names and convert them into
    # integers.
    # This assumes that the names of area are always like '1A'.
    areas = np.array([int(x[:-1]) for x in df["area"].unique()])
    n_areas = len(areas)

    # Make a one long time series from the DataFrame
    df_ts = df.loc[:, 1:].transpose().melt().drop(
        'variable', axis=1).dropna()
    years = df["year"].unique()
    # Years are water years, so one needs to be subtracted to get the calendar
    # year.
    start_date = pd.to_datetime(f"{years[0] - 1}-10-01")
    end_date = pd.to_datetime(f"{years[-1]}-09-30")
    # Add areas and dates columns and make them indices
    df_ts["area"] = pd.Categorical(
        areas.repeat(len(pd.date_range(start=start_date, end=end_date))),
        categories=areas)
    dates = pd.date_range(start=start_date, end=end_date, freq='D')
    df_ts["date"] = np.tile(dates, n_areas)
    da = xr.DataArray(data=df_ts.pivot(columns="area",
                                       index="date",
                                       values="value").values,
                      dims=["date", "area"],
                      coords=dict(date=dates, area=areas))
    return da


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


def calculate_groundwater_rates(df_subareas, dates):
    # Read groundwater ratios
    # FIXME hard-wired file name
    # FIXME The groundwater rate data from the file may need to be clipped.
    path_gwrates = "GW_RATES.TXT"
    df_gw_rates = pd.read_csv(path_gwrates, delim_whitespace=True, header=None)
    df_gw_rates.rename(columns={0: "year", 1: "gw_rate"}, inplace=True)

    # Create an empty groundwater array
    n_dates = len(dates)
    areas = df_subareas["area"].values
    n_areas = len(areas)
    da_gwrates = xr.DataArray(data=np.full((n_dates, n_areas), np.nan),
                              dims=["date", "area"],
                              coords=dict(date=dates, area=areas))

    areas_lower = df_subareas.query("uplow == 1 & docregion == 'Lower'")[
        'area'].values
    areas_mid = df_subareas.query("uplow == 1 & docregion == 'Midrange'")[
        'area'].values
    areas_high = df_subareas.query("uplow == 1 & docregion == 'High'")[
        'area'].values

    for year in range(dates[0].year + 1, dates[-1].year + 1):
        da_gwrates.loc[f"{year - 1}-10-01":f"{year}-09-30",
                       :] = df_gw_rates.query("year == @year")["gw_rate"].values
    da_gwrates.loc[:, areas_lower] = 0.35
    da_gwrates.loc[:, areas_mid] = 0.30
    da_gwrates.loc[:, areas_high] = 0.25
    return da_gwrates


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
    dates = pd.date_range(da_lwa.date.values[0], da_lwa.date.values[-1])
    years = range(dates[0].year + 1, dates[-1].year + 1)
    areas = da_lwa.area.values
    n_areas = len(areas)
    for year in years:
        leach_saved = np.zeros((n_areas, ))
        lwa = da_lwa.loc[f"{year - 1}-10-01":f"{year}-09-30", :].values
        lwd = da_lwd.loc[f"{year - 1}-10-01":f"{year}-09-30", :].values
        ro = da_ro.loc[f"{year - 1}-10-01":f"{year}-09-30", :].values
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


def main():
    """ A main function to run DCD

        FIXME this will be renamed and refactored later.
    """
    # Read a water year types
    # FIXME hard-wired file name
    path_wy_types = "WYTYPES"
    df_wy_types = read_water_year_types(path_wy_types)

    # Read a subarea information file
    # FIXME hard-wired file name
    path_subareas = "subarea-info.TXT"
    df_subareas = read_subarea_info_file(path_subareas)

    # Create an area array.
    # This will be an coordinates in many arrays.
    areas = df_subareas["area"].values
    n_areas = len(df_subareas)

    # Read irrigation efficiency, \eta
    # FIXME hard-wired file name
    path_eta = "IRREFF-3MWQIregions"
    df_eta = pd.read_csv(path_eta, header=None)
    da_eta = xr.DataArray(data=df_eta[0].values, dims=[
                          "area"], coords=dict(area=areas))

    # Read Applied Water (I_A) from DETAW outputs
    # FIXME hard-wired file name
    path_aw = "DICU5.17"
    da_aw = read_detaw_output_to_dataarray(path_aw)

    # Read (effective) seepage (S_E) from DETAW outputs
    # FIXME hard-wired file name
    path_seepage = "DICU5.14"
    da_seepage = read_detaw_output_to_dataarray(path_seepage)

    # Read runoff (RO) from DETAW outputs
    path_runoff = "DICU5.12"
    da_runoff = read_detaw_output_to_dataarray(path_runoff)

    # Read depletion from waterbody
    path_wncu = "DICU5.30"
    da_wncu = read_detaw_output_to_dataarray(path_wncu)

    # FIXME hard-wired dates.
    start_date = pd.to_datetime("2015-10-01")
    end_date = pd.to_datetime("2017-09-30")

    # Create the dates in the date range.
    # This will be coordinates in many arrays later.
    dates = pd.date_range(start=start_date, end=end_date, freq='D')
    n_dates = len(dates)
    # Create an array of days in each month for later use
    days_in_month = dates.days_in_month
    da_days_in_month = xr.DataArray(data=days_in_month,
                                    dims=["date"],
                                    coords=dict(date=dates))

    # Create months in the modeling period
    months = pd.period_range(start=start_date, end=end_date, freq='M')

    # Months in the ordering in a water year, which is from October
    # to September, 10, 11, 12, 1, ..., 9
    wy_months = np.roll(np.arange(1, 13), 3)

    # Read monthly applied leach water (LW_A) and drained leach water (LW_D)
    # FIXME hard-wired file name
    # FIXME Factor out the block
    path_lwam = "LEACHAPL.DAT"
    df_lwam = pd.read_csv(path_lwam, delim_whitespace=True,
                          header=None).loc[:, 1:].transpose()
    # Repeat the leach water data for convenience.
    n_years = end_date.year - start_date.year
    df_lwam = pd.DataFrame(np.tile(df_lwam.values, (n_years, 1)), index=months)
    df_lwam.rename(columns={i: i+1 for i in df_lwam.columns}, inplace=True)
    df_lwa = df_lwam.resample('1D').ffill()
    da_lwa = xr.DataArray(data=df_lwa.values, dims=["date", "area"],
                          coords=dict(date=dates, area=areas))
    da_lwa /= da_days_in_month

    path_lwdm = "LEACHDRN.DAT"
    df_lwdm = pd.read_csv(path_lwdm, delim_whitespace=True,
                          header=None).loc[:, 1:].transpose()
    # Repeat the leach water data for convenience.
    n_years = end_date.year - start_date.year
    df_lwdm = pd.DataFrame(np.tile(df_lwdm.values, (n_years, 1)), index=months)
    df_lwdm.rename(columns={i: i+1 for i in df_lwdm.columns}, inplace=True)
    df_lwd = df_lwdm.resample('1D').ffill()
    da_lwd = xr.DataArray(data=df_lwd.values, dims=["date", "area"],
                          coords=dict(date=dates, area=areas))
    da_lwd /= da_days_in_month

    # ---------------------------------------------------------------
    # Calculate channel depletion

    # Calculate drained seepage (S_D) for each subarea
    da_drnseep = calculate_drained_seepage(df_subareas)
    # Divide it by days in months
    da_drnseep = da_drnseep / da_days_in_month

    # Update seepage (S) by adding drained seepage (S_D)
    # Calculate daily seepage values by dividing by days in months
    da_seepage += da_drnseep

    # Calculate daily groundwater rate per area
    da_gwrates = calculate_groundwater_rates(df_subareas, dates)

    # Calculate groundwater component 1
    da_gw1 = da_gwrates * da_aw / da_eta

    # Calculate groundwater component 2
    da_gw2 = da_gwrates * da_seepage

    # Calculate total groundwater
    da_gwf = da_gw1 + da_gw2

    # FIXME hard-wired value
    leach_scale = 5.0
    da_lwa *= leach_scale

    # FIXME hard-wired value
    runoff_rate = 0.75
    # Update the runoff (RO) with the coefficient
    da_runoff *= runoff_rate

    # Adjust leach water
    da_lwa, da_lwd = adjust_leach_water(da_lwa, da_lwd, da_runoff)

    # Calculate diversion without seepage
    # NOTE Routines to adjust runoff in leach water are not implemented.
    # NOTE FIXME 1977 adjustment of applied water is not implemented.
    da_diversion = da_aw / da_eta - da_gw1 + da_wncu + da_lwa

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

    path_dcd = "dcd.nc"
    ds_dcd.to_netcdf(path_dcd)


if __name__ == "__main__":
    main()
