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

from pathlib import Path
import pandas as pd
import pyhecdss
import os
import sys
import shutil
import platform
import numpy as np
import xarray as xr
import dcd_postprocess


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
    SKIP_DETAW = False  # FIXME make this an option
    if not SKIP_DETAW:
        status = os.system('python DETAW.py')
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


def main():
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

    shutil.rmtree("./NODCU/DCD_Cal/DCD_outputs", ignore_errors=True)
    os.chdir(owd)


if __name__ == "__main__":
    main()
