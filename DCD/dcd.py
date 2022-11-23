# Delta Channel Depletion Model (DCDv1.2)
# <license>
#    Copyright (C) State of California, Department of Water Resources.
#    This file is part of DCDv1.2.

#    DCDv1.1 is free software:
#    you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    DCDv1.1 is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with DCDv1.2.  If not, see <http://www.gnu.org/licenses>.
# </license>
#
# This version 1.2 was developed for supporting CALSIM3,DSM2 and SCHISM studies.
#
# Enter the following line in the command window to run the model:
#    Python DCD1.2.py
#

import dcd_postprocess
from pathlib import Path
import argparse
import logging
import yaml

import pandas as pd
import pyhecdss
import os
import sys
import shutil
import platform
import numpy as np
import xarray as xr
# FIXME Remove this path kludge later
sys.path.append("../DETAW")

# Set up a logger
logging.basicConfig(level=logging.INFO)


def callDETAW(supmodel, leachoption):
    owd = os.getcwd()
    dir_dst = "../DETAW/"
    os.chdir(dir_dst)
    months = ["Jan", "Feb", "Mar", "Apr", "May", "Jun",
              "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]

    if supmodel == 1 or supmodel == 2:
        inputfile = "./Input/historical_study/LODI_PT.csv"
    elif supmodel == 3:
        inputfile = "./Input/planning_study/LODI_PT.csv"

    f0 = open(inputfile, 'r')
    templ = ""
    endyear = 0
    for line in f0:
        if line:
            templ = line
            if line.split(",")[5].strip() == "0" and line.split(",")[6].strip() == "0":
                break
    endyear = templ.split(",")[1]
    endmonth = int(templ.split(",")[2])
    outputfile = "DCD_"+months[endmonth-1] + \
        endyear+"_Lch"+str(leachoption)+".dss"
    SKIP_DETAW = True  # FIXME make this an option
    if not SKIP_DETAW:
        # FIXME Using a hardwired file name for now...
        path_detaw_yaml = "../DETAW/Input/detaw.yaml"
        detaw.detaw(path_detaw_yaml)
    os.chdir(owd)
    return(endyear, outputfile)


def get_kernel_exe():
    csys = platform.system()
    if csys == 'Windows':
        return '.\DCD_kernel.exe'
    elif csys == 'Linux':
        return './DCD_kernel'
    else:
        raise 'Unsupported platform: %s' % csys


# module level global for temp working space for DCD_kernel
DCD_OUTPUT_TEMP_DIR = './DCD_outputs'


def ensure_output_dirs(owd):
    out_dir = os.path.join(owd, 'Output')
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)
    subdirs = [os.path.join(out_dir, dir)
               for dir in ['DSM2', 'SCHISM', 'CALSIM']]
    for sdir in subdirs:
        if not os.path.exists(sdir):
            os.mkdir(sdir)


def change_dir_and_copy_files():
    """ Create a temporary working directory and move into it. """
    dir_dst = "./NODCU/DCD_Cal/"
    os.chdir(dir_dst)
    if not os.path.exists(DCD_OUTPUT_TEMP_DIR):
        os.mkdir(DCD_OUTPUT_TEMP_DIR)
    shutil.copy(get_kernel_exe(), DCD_OUTPUT_TEMP_DIR)
    shutil.copy('WYTYPES', DCD_OUTPUT_TEMP_DIR)
    os.chdir(DCD_OUTPUT_TEMP_DIR)


def set_env_vars_for_nodcu(supmodel, leachoption, endyear, outputfile,
                           extension):
    """ Set environment variables for the Fortran code, NODCU

        Parameters
        ----------
        supmodel: str
            Target model for processing.
        leachoption: str
            Leach scale factor option
        endyear: str
            End year of the processing
        outputfile: str
            Output file from DETAW through this process
        extenion: str
            Extension name, empty string for the base case

        Returns
        -------
        None
    """
    ext = f'_{extension}' if extension else ''
    os.environ['DICU5_14'] = f'../../../../DETAW/Output/DICU5{ext}.14'
    os.environ['DICU5_17'] = f'../../../../DETAW/Output/DICU5{ext}.17'
    os.environ['DICU5_12'] = f'../../../../DETAW/Output/DICU5{ext}.12'
    os.environ['DICU5_27'] = f'../../../../DETAW/Output/DICU5{ext}.27'
    os.environ['DICU5_30'] = f'../../../../DETAW/Output/DICU5{ext}.30'

    if supmodel == 3:
        os.environ['GW_RATES_TXT'] = '../../../NODCU/GW_RATES_CALSIM3.TXT'
    else:
        # update data in the file each year ----no adjustment-GW_RATES.TXT
        os.environ['GW_RATES_TXT'] = '../../../NODCU/GW_RATES.TXT'
    # set for DETAW-CD
    os.environ['GW_LOWLANDS_TXT'] = '../../../NODCU/GW_LOWLANDS.TXT'
    # os.environ['DIVFCTR_RMA']='../../../NODCU/DIVFCTR.2020'
    # os.environ['DRNFCTR_RMA']='../../../NODCU/DRNFCTR.2020'
    os.environ['LEACHAPL_DAT'] = f'../../../NODCU/LEACHAPL{ext}.DAT'
    os.environ['LEACHDRN_DAT'] = f'../../../NODCU/LEACHDRN{ext}.DAT'
    os.environ['IDRNTDS_DAT'] = '../../../NODCU/IDRNTDS.DAT'
    os.environ['DRNCL_123'] = '../../../NODCU/DRNCL.123'
    os.environ['GEOM_NODES'] = '../../../NODCU/GEOM-NODES-1.5'

    os.environ['IRREFF_DAT'] = '../../../NODCU/IRREFF-3MWQIregions'
    os.environ['subarea_info'] = f'../../../NODCU/subarea-info{ext}.TXT'

    # Runtime variables
    # The years assumed are incorrect, so 'N'
    os.environ['years_ok'] = 'N'
    # The correct beginning year to run is
    os.environ['begwy'] = '2016'
    # The correct last year to run is
    os.environ['endwy'] = endyear
    # Type of drainage concentration data (1 for TDS, 2 for chloride)
    os.environ['datatype'] = '1'
    # Do you want to creat an ascii file?
    os.environ['ascii'] = 'Y'
    # The dss file to save output
    os.environ['dssfile'] = outputfile
    # The leach scale factor
    os.environ['leachscale'] = str(leachoption)


def run_nodcu():
    """ Run NODCU (the Fortran code) with the environment variables
        set in the previous steps.
    """
    status = os.system(get_kernel_exe())
    files_to_convert = ['roisl.txt', 'gwbyisl.txt', 'drn_wo_ro_isl.txt',
                        'div_wo_spg_isl.txt', 'spgisl.txt']
    list_datasets = []
    for pathfile in files_to_convert:
        da = read_nodcu_out(pathfile)
        ds = da.to_dataset(name=da.attrs['name'])
        list_datasets.append(ds)
    ds_all = xr.merge(list_datasets)
    return ds_all


def read_nodcu_out(filepath):
    """ Read NODCU output and return in a pandas.DataFrame.

        The function reads a NODCU output file. It is assumed that the
        ordering of the data blocks (time series) are in an ascending order
        of the island indices.

        Parameters
        ----------
        filepath
            file name to read

        Returns
        -------
        xarray.DataArray
    """
    # Read the file to figure out the number of years in the file
    df = pd.read_csv(filepath, skiprows=5, header=None, nrows=200*365)
    n_count = df[df[0] == 'END'].index[0]
    # Get some info from the header
    with open(filepath, 'r') as fh:
        # First line for the file name
        l_name = fh.readline().strip()
        # Second line for DSS parts
        l_parts = fh.readline()
        parts = [x.split()[0] for x in l_parts.split("=")[1:]]
        l_unit = fh.readline().strip()
        l_period_type = fh.readline().strip()
        l_start_time = fh.readline().strip()

    # Read the file again while skipping non-data parts
    df = pd.read_csv(filepath,
                     skiprows=lambda x: x % (n_count + 5) < 5,
                     dtype='f4',
                     header=None)
    # Create a 2-d array
    n_islands = df.shape[0] // n_count
    data = np.hstack(
        [df.loc[i * n_count:(i + 1) * n_count - 1, :].values
         for i in range(n_islands)])

    # Put dates column
    start_date = pd.to_datetime(l_start_time.split()[0])
    end_data = start_date + pd.to_timedelta(n_count - 1, unit='D')
    dates = pd.date_range(start=start_date, end=end_data)

    # Create a DataArray
    da = xr.DataArray(data,
                      dims=["date", "island"],
                      coords=dict(date=dates,
                                  island=np.arange(n_islands) + 1))

    # Add metadata
    da.attrs["name"] = l_name.split(".")[0]
    da.attrs["part_a"] = parts[0]
    da.attrs["part_c"] = parts[2]
    da.attrs["part_e"] = parts[4]
    da.attrs["part_f"] = parts[5]
    da.attrs["unit"] = l_unit
    da.attrs["period_type"] = l_period_type
    return da


def callDCD(supmodel, leachoption, endyear, outputfile, extension):
    """ Run NODCU (DCD_kernel)

        Parameters
        ----------
        supmodel
            Target model for processing.
        leachoption
            Leach scale factor option
        endyear
            End year of the processing
        outputfile
            Output file from DETAW through this process
        extension
            Extension option name

        Returns
        -------
        None
            The function creates a few processed DSS files
    """
    owd = os.getcwd()
    change_dir_and_copy_files()

    set_env_vars_for_nodcu(supmodel, leachoption, endyear, outputfile,
                           extension)

    ds_all = run_nodcu()

    # First, calculate deep percolation
    da_ro_island = ds_all.RO_island
    # Assuming 25% of the surface runoff goes to the deep percolation
    # See Annual Report 2017, Chapter 3.
    dp_factor = 0.25
    da_dp_island = da_ro_island * dp_factor
    da_dp_island.attrs = da_ro_island.attrs
    da_dp_island.attrs["name"] = "DP_island"
    da_dp_island.attrs["part_c"] = "DP-FLOW"
    ds_all[da_dp_island.attrs["name"]] = da_dp_island
    # Save the data for later uses
    suffix = f"_{extension}" if extension != "" else ""
    path_out = f"DCD_islands{suffix}.nc"
    ds_all.to_netcdf(path_out)

    # Take monthly averages of the NODCU outputs
    names_to_process = ["drn_wo_ro_island", "div_wo_spg_island",
                        "GW_per_island", "RO_island",
                        "DP_island", "spg_island"]
    for name in names_to_process:
        da = ds_all[name]
        da_mon = take_monthly_average(da)
        name_mon = f"{name}_mon"
        ds_all[name_mon] = da_mon

    # Write monthly data
    name_out = "DCD_island_month"
    if extension == "":
        ensure_output_dirs(owd)
        if supmodel == 1:  # DSM2
            model = "DSM2"
        elif supmodel == 2:
            model = "SCHISM"
        elif supmodel == 3:
            model = "CALSISM3"
        else:
            raise ValueError()
        path_output = Path(owd) / f"Output/{model}/{name_out}.nc"
        names_to_write = [f"{x}_mon" for x in names_to_process]
        ds_all[names_to_write].to_netcdf(path_output)
    else:
        if extension == "ex3":
            if supmodel == 3:
                path_output = Path(owd) / f"Output/CALSIM3/{name_out}.nc"
                # FIXME The final output writing is better not to be here.
                names_to_write = [f"ext_{x}_mon" for x in names_to_process]
                ds_all[names_to_write].to_netcdf(path_output)
            dcd_postprocess.split_bbid_new()
    os.chdir(owd)


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


def get_pathname(dssfh, path):
    '''
    reads catalog to get full pathname if it exists
    The assumption is that the path provided does not have a time window
    returns the first matching pathname with exactly the A,B,C,E and F parts
    '''
    pathparts = path.split('/')
    dfcat = dssfh.read_catalog()
    # @dfpath=dfcat[(dfcat.A==pathparts[1]) & (dfcat.B==pathparts[2])
    # @                & (dfcat.C==pathparts[3]) & (dfcat.E==pathparts[5])
    # @                & (dfcat.F==pathparts[6])]
    dfpath = dfcat[(dfcat.B == pathparts[2])]
    pathname = dssfh.get_pathnames(dfpath)[0]
    return pathname


def daytomonth(inputfile, outputfile):
    d = pyhecdss.DSSFile(inputfile)
    #outputfile = os.getcwd()+"/"+inputfile.split(".")[0]+"_mon.dss"
    do = pyhecdss.DSSFile(outputfile, create_new=True)
    plist = d.get_pathnames()
    for p in plist:
        df, u, p = d.read_rts(p)
        do.write_rts(df.columns[0].replace('1DAY', '1MON'),
                     df.resample('M').mean(), u, 'PER-AVER')
    d.close()
    do.close()


def changepaths(inDSSfile, pathfile, outDSSfile, EPART):
    f0 = open(pathfile, 'r')
    islands = []
    for line in f0:
        if line:
            islands.append(line)

    dssifh = pyhecdss.DSSFile(inDSSfile)
    dssofh = pyhecdss.DSSFile(outDSSfile, create_new=True)
    for i in range(len(islands)):
        templ = islands[i]
        pathin = "//"+templ.split(",")[0].strip()+"/////"
        cpart = get_pathname(dssifh, pathin).split("/")[3]
        tdss, cunits, ctype = dssifh.read_rts(get_pathname(dssifh, pathin))
        pathout = "/"+templ.split(",")[1].strip()+"/"+templ.split(
            ",")[0].strip()+"/"+cpart+"//"+EPART+"/"+templ.split(",")[2].strip()+"/"
        dssofh.write_rts(pathout, tdss.shift(freq='D'), cunits, ctype)
    dssifh.close()
    dssofh.close()


def callDCD_ext(supmodel, leachoption, endyear, outputfile, extension):
    owd = os.getcwd()
    change_dir_and_copy_files()

    set_env_vars_for_nodcu(supmodel, leachoption,
                           endyear, outputfile, extension)

    run_nodcu()

    dcd_postprocess(outputfile, "base")

    filestocopy_day = ["spg_island.dss", "RO_island.dss", "drn_wo_ro_island.dss",
                       "div_wo_spg_island.dss", "DP_island.dss", "GW_per_island.dss"]
    os.chdir(owd)


def createoutputs(outputfile, modeloption):
    owd = os.getcwd()

    dir_dst = "./NODCU/DCD_Cal/"
    os.chdir(dir_dst)
    if not os.path.exists(DCD_OUTPUT_TEMP_DIR):
        os.mkdir(DCD_OUTPUT_TEMP_DIR)
    shutil.copy(get_kernel_exe(), DCD_OUTPUT_TEMP_DIR)
    shutil.copy('WYTYPES', DCD_OUTPUT_TEMP_DIR)
    os.chdir(DCD_OUTPUT_TEMP_DIR)

    extension = f"out_{modeloption}"
    dcd_postprocess.dcd_postprocess(outputfile, extension)
    # shutil.rmtree(DCD_OUTPUT_TEMP_DIR)
    os.chdir(owd)


def main_dcd_old():
    """ (Deprecated) Old DCD main function
    """
    pyhecdss.set_message_level(0)
    owd = os.getcwd()
    modelparafile = "./NODCU/DCD_parameters.inp"
    fmp = open(modelparafile, "r")
    modeloption = 0
    outputfile = ""
    for line in fmp:
        if line:
            if not("#" in line):
                if "Model to streamline" in line:
                    modeloption = int(line.split("=")[1])
                if "Leach scale factor" in line:
                    leachoption = int(line.split("=")[1])
    (endyear, outputfile) = callDETAW(modeloption, leachoption)
    callDCD(modeloption, leachoption, endyear, outputfile, "")
    callDCD(modeloption, leachoption, endyear, outputfile, "ex1")
    callDCD(modeloption, leachoption, endyear, outputfile, "ex2")
    callDCD(modeloption, leachoption, endyear, outputfile, "ex3")
    createoutputs(outputfile, modeloption)

    # shutil.rmtree("./NODCU/DCD_Cal/DCD_outputs", ignore_errors=True)
    os.chdir(owd)


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


def calculate_groundwater_rates(df_subareas, df_gw_rates, dates):
    """
    """
    # Create an empty groundwater array
    n_dates = len(dates)
    areas = df_subareas["area"].values
    n_areas = len(areas)
    da_gwrates = xr.DataArray(data=np.full((n_dates, n_areas), np.nan),
                              dims=["time", "area"],
                              coords=dict(time=dates, area=areas))

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


def calculate_depletion(model_params: dict) -> None:
    """ Calculate channel depletion

        Parameters
        ----------
        model_params: dict
           a dict of model parameters (from an input file)

        Returns
        -------
        None
    """
    logging.info("Reading input files ")
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
    path_groundwater_rates = model_params.get("path_groundwater_rates")
    df_gw_rates = pd.read_csv(path_groundwater_rates,
                              delim_whitespace=True, header=None)
    df_gw_rates.rename(columns={0: "year", 1: "gw_rate"}, inplace=True)

    # Read monthly applied leach water (LW_A) and drained leach water (LW_D)
    path_lwam = model_params.get("path_leach_applied")
    df_lwam = pd.read_csv(path_lwam, delim_whitespace=True,
                          header=None).loc[:, 1:].transpose()
    path_lwdm = model_params.get("path_leach_drained")
    df_lwdm = pd.read_csv(path_lwdm, delim_whitespace=True,
                          header=None).loc[:, 1:].transpose()

    # Read a DETAW NetCDF file
    path_detaw = model_params.get("path_detaw_output")
    ds_detaw = xr.open_dataset(path_detaw)

    start_water_year = model_params.get("start_water_year")
    end_water_year = model_params.get("end_water_year")
    start_date = pd.to_datetime(f'{start_water_year - 1}-10-01')
    end_date = pd.to_datetime(f'{end_water_year}-09-30')
    n_years = end_date.year - start_date.year

    # Preproces DETAW outputs
    logging.info("Preprocessing DETAW output...")
    # Aggregating by areas
    da_aw = ds_detaw.et_aw.sum('crop')
    da_seepage = ds_detaw.s_e.sum('crop')
    da_runoff = (ds_detaw.precip - ds_detaw.e_r).sum(
        'crop').rolling(time=5, min_periods=1).mean()
    # da_depletion = ds_detaw.et_c.sum('crop')
    cropname = 'Water Surface'
    da_waterbody = (ds_detaw.et_c.sel(
        crop=cropname) - ds_detaw.precip.sel(crop=cropname)).clip(0.)

    # Create the dates in the date range.
    # This will be coordinates in many arrays later.
    dates = pd.date_range(start=start_date, end=end_date, freq='D')
    # n_dates = len(dates)
    # Create an array of days in each month for later use
    days_in_month = dates.days_in_month
    da_days_in_month = xr.DataArray(data=days_in_month,
                                    dims=["time"],
                                    coords=dict(time=dates))

    # Create months in the modeling period
    months = pd.period_range(start=start_date, end=end_date, freq='M')

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
    # Calculate daily seepagnnn values by dividing by days in months
    da_seepage += da_drnseep

    # Calculate daily groundwater rate per area
    da_gwrates = calculate_groundwater_rates(df_subareas, df_gw_rates, dates)

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
