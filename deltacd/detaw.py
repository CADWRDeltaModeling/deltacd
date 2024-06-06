# Delta Evapotranspiration of Applied Water(DETAWv2.0)
# <license>
#    Copyright (C) State of California, Department of Water Resources.
#    This file is part of Delta Evapotranspiration of Applied Water
#    (DETAWv2.0).

#    DETAWv2.0 is free software:
#    you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    DETAWv2.0 is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with DETAWv2.0.  If not, see <http://www.gnu.org/licenses>.
# </license>
#
# Enter the following line in the command window to run the model:
#    Python detaw.py
#

import sys
from pathlib import Path

if not sys.warnoptions:
    import warnings
    warnings.simplefilter("ignore")

import argparse
import logging
import pandas as pd
import xarray as xr
from numpy import add, pi, array, zeros

import os
import string
import math
import numpy
from os import listdir
from math import cos, sin, tan, atan, sqrt, pi, pow

import timeit
import numba
import yaml


DEBUG_TIMING = True
DEBUG_OUTPUT = False
NO_OUTPUT = True

# Set up a logger
# FIXME The log level can be managed by the user input.
logging.basicConfig(level=logging.INFO)

def read_and_clean_crop_info(source):
    '''
    reads critical.csv or non-critical.csv for crop information for critical and non-critical types
    return the cleaned up dataframe
    '''
    #equivalent in pandas
    df = pd.read_csv(source, skiprows=6)
    df = df.rename(columns={'Crop': 'Data Type'})
    df = df.drop(columns=[df.columns[1]])
    df = df.drop(range(20, 27))
    return df


def fill_zeros_with_last(arr):
    '''
    fill in zeros with previous value.
    https://stackoverflow.com/questions/30488961/fill-zero-values-of-1d-numpy-array-with-last-non-zero-values
    '''
    prev = numpy.arange(len(arr))
    prev[arr <= 0] = 0
    prev = numpy.maximum.accumulate(prev)
    return arr[prev]


def write_to_netcdf(detawoutput, model_start_year, fn_detaw_output):
    ''' Write DETAW output array into a NetCDF

        This function writes out DETAW outputs in a numpy array into a NetCDF
        file after converting the array into an xarray.Dataset. The converted
        dataset is returned.

        Parameters
        ----------
        detawoutput: numpy.ndarray
            Output array from `historicalETAW`. The dimensions of the array are
            (etvars, areas, crops, times).
        model_start_year: int
            start year in water year
        fn_detaw_output: str
            filename to write the output

        Returns
        -------
        xarray.Dataset
            DETAW output converted into xarray.Dataset.
    '''
    etvars = ["et_c", "s_e", "precip", "et_aw", "d_sw", "e_r"]
    dims = ['area_id', 'crop', 'time']
    if detawoutput.shape[2]-1 == 15: # The Delta has 15 landuse crop categories
        coords = {'area_id': numpy.arange(detawoutput.shape[1]-1, dtype='i4')+1,
              'crop': ["Urban", "Irrig pasture", "Alfalfa", "All field",
                       "Sugar beets", "Irrig grain", "Rice", "Truck crops",
                       "Tomato", "Orchard", "Vineyard", "Riparian vegetation",
                       "Native vegetation", "Non-irrig grain", "Water surface"],
              'time': pd.date_range(str(model_start_year) + '-10-01',
                                    periods=detawoutput.shape[-1])}  # last dimension is time
    elif detawoutput.shape[2]-1 == 16: # Suisun Marsh has one more landuse crop category than the Delta
        coords = {'area_id': numpy.arange(detawoutput.shape[1]-1, dtype='i4')+1,
              'crop': ["Urban", "Irrig pasture", "Alfalfa", "All field",
                       "Sugar beets", "Irrig grain", "Rice", "Truck crops",
                       "Tomato", "Orchard", "Vineyard", "Riparian vegetation",
                       "Native vegetation", "Non-irrig grain", "Water surface","Duck Pond"],
              'time': pd.date_range(str(model_start_year) + '-10-01',
                                    periods=detawoutput.shape[-1])}  # last dimension is time

    ds = xr.Dataset({etvar: xr.DataArray(detawoutput[i, :-1, :-1, :],
                                         dims=dims,
                                         coords=coords, attrs={'units': 'Acre-feet'})
                     for i, etvar in enumerate(etvars)})

    head_tail = os.path.split(fn_detaw_output)
    if not os.path.exists(head_tail[0]):
        os.mkdir(head_tail[0])
    ds.to_netcdf(fn_detaw_output)
    return ds


def weatheroutput_to_netcdf(pcp, ET0, model_start_year, fn_precip_output, fn_et_output):
    '''
    write the precip and ET0 for all areas to netcdf4 format
    '''
    dpcp = xr.DataArray(pcp,
                        dims=['time', 'area'],
                        coords={'time': pd.date_range(
                            str(model_start_year) + '-10-01', periods=pcp.shape[0], freq='D'), 'area': numpy.arange(pcp.shape[-1], dtype='i4')+1, },
                        attrs={'units': 'mm'},
                        name='precip')
    det0 = xr.DataArray(ET0,
                        dims=['time', 'area'],
                        coords={'time': pd.date_range(
                            str(model_start_year) + '-10-01', periods=pcp.shape[0], freq='D'), 'area': numpy.arange(pcp.shape[-1], dtype='i4')+1, },
                        attrs={'units': 'mm'},
                        name='ET0')

    head_tail = os.path.split(fn_precip_output)
    if not os.path.exists(head_tail[0]):
        os.mkdir(head_tail[0])
    head_tail = os.path.split(fn_et_output)
    if not os.path.exists(head_tail[0]):
        os.mkdir(head_tail[0])
    dpcp.to_netcdf(fn_precip_output)
    det0.to_netcdf(fn_et_output)


def write_to_dss(dssfh, arr, path, startdatetime, cunits, ctype):
    '''
    write to the pyhecdss.DSSFile for an array with starttime and assuming
    daily data with the pathname path, cunits and ctype
    '''
    if NO_OUTPUT:
        return
    fstr = '1D'
    epart = path.split('/')[5]
    if epart == '1DAY':
        fstr = '1D'
    elif epart == '1MONTH':
        fstr = '1M'
    else:
        raise RuntimeError('Not recognized frequency in path: %s' % path)
    # df=pd.DataFrame(arr,index=pd.date_range(startdatetime,periods=len(arr),freq=fstr))
    #write_dataframe(dssfh, df, path, cunits, ctype)
    sp = pd.to_datetime(startdatetime)
    darr = numpy.array(arr, dtype='d')
    pyhecdss.pyheclib.hec_zsrtsxd(dssfh.ifltab, path,
                                  sp.strftime("%d%b%Y").upper(), sp.round(
                                      freq='T').strftime("%H%M"),
                                  darr, cunits[:8], ctype[:8])


def write_dataframe(dssfh, df, path, cunits, ctype='INST-VAL'):
    ''' write data frame to DSS file handle '''
    dssfh.write_rts(path, df, cunits, ctype)


def write_weather_dss(dftmax, dftmin, dfptotal, dfet0, outputfile):
    ctype = "INST-VAL"  # "PER-AVER"
    pyhecdss.set_message_level(0)
    pyhecdss.set_program_name('DETAW')
    dssfh = pyhecdss.DSSFile(outputfile, create_new=True)
    path = "/detaw/LODI_Tmax/Temp//1DAY/detaw/"
    write_dataframe(dssfh, dftmax, path, 'oC', ctype)
    path = "/detaw/LODI_Tmin/Temp//1DAY/detaw/"
    write_dataframe(dssfh, dftmin, path, 'oC', ctype)
    path = "/detaw/LODI_Tmax/Temp//1MONTH/detaw/"
    write_dataframe(dssfh, dftmax.resample('M').mean(), path, 'oC', ctype)
    path = "/detaw/LODI_Tmin/Temp//1MONTH/detaw/"
    write_dataframe(dssfh, dftmin.resample('M').mean(), path, 'oC', ctype)
    for j in range(0, ilands):
        precip_area = dfptotal.iloc[:, j].to_frame()
        et0_area = dfet0.iloc[:, j].to_frame()
        path = "/detaw/island_"+str(j+1)+"/precipitation//1DAY/detaw/"
        write_dataframe(dssfh, precip_area, path, 'mm', ctype)
        path = "/detaw/island_"+str(j+1)+"/ET0//1DAY/detaw/"
        write_dataframe(dssfh, et0_area, path, 'mm', ctype)
        path = "/detaw/island_"+str(j+1)+"/precipitation//1MONTH/detaw/"
        write_dataframe(dssfh, precip_area.resample(
            'M').mean(), path, 'mm', ctype)
        path = "/detaw/island_"+str(j+1)+"/ET0//1MONTH/detaw/"
        write_dataframe(dssfh, et0_area.resample(
            'M').mean(), path, "mm", ctype)
    dssfh.close()


def weatheroutput(ts_pcp, ts_per, ts_mon, ts_days, Tmax, Tmin, ilands, idates, isites, ETo_corrector, filepath, start1):
    """
        calculate the precipitation and reference evapotranspiration for each island

    input: daily precipitation at seven sites
           percentages of seven sites for 168 islands
    output: Nothing,
            Generate daily file and monthly weather file.
    """
    monthname = ["JAN", "FEB", "MAR", "APR", "MAY",
                 "JUN", "JUL", "AUG", "SEP", "OCT", "NOV", "DEC"]
    # precipitation is the product of precip stations times the distribution percentages
    ts_ptotal = numpy.dot(numpy.transpose(ts_pcp), ts_per)
    # limit to > 0 and < 9990 ?
    numpy.clip(ts_ptotal, 0, 9990)
    # set Tmax to Tmin where Tmax < Tmin (why should this be?)
    # should we not just exchange the values of min and max ?
    bad_vals = numpy.where(Tmax < Tmin)
    Tmax[bad_vals] = Tmin[bad_vals]
    TDiffTemp = Tmax-Tmin
    Tm = 0.5*(Tmax+Tmin)
    # Global constants for calcs
    Lat = 38.5  # Mean latitude
    PHI = math.pi*Lat/180  # phi angle with equator line ?
    GSC = 0.082  # Global Solar constant
    Lam = 2.45  # Latent heat of Vaporization (units?)
    # Calculate ETo with Hargres-Samani equation
    df = 1 + 0.033*numpy.cos(2*math.pi/365*ts_days)  # Earth-sun distance
    dec1 = 0.409*numpy.sin(2*math.pi/365*ts_days-1.39)  # Delination of sun
    arg = -numpy.tan(PHI)*numpy.tan(dec1)  # argument for tangant in radians
    argsq = arg*arg
    # WS = sunrise hour angle in radians
    WS = pi/2.-numpy.arctan(arg/numpy.sqrt(1-argsq))
    WS[numpy.where(argsq >= 1)] = pi/2
    cosz = WS*numpy.sin(dec1)*math.sin(PHI) + \
        (numpy.cos(dec1)*math.cos(PHI)*numpy.sin(WS))
    Ra = (24*60/math.pi)*GSC*df*cosz  # Extrater Radiation
    ET0 = (0.0023*Ra*numpy.sqrt(TDiffTemp)*(Tm+17.8))/Lam
    numpy.clip(ET0, 0, None, out=ET0)
    ET0Daily = ETo_corrector*ET0.reshape(idates, 1)
    # start data & time
    # The start time for monthly interval must be the first day of each month.
    startdate = str(start1[2])+monthname[start1[1]-1]+str(start1[0])
    starttime = str(start1[3])+"00"
    # create data frames with time index
    dtindex = pd.date_range(startdate+'T'+starttime,
                            periods=len(Tmax), freq='D')
    dftmax = pd.DataFrame(Tmax, index=dtindex)
    dftmin = pd.DataFrame(Tmin, index=dtindex)
    dfptotal = pd.DataFrame(ts_ptotal, index=dtindex)
    dfet0 = pd.DataFrame(ET0Daily, index=dtindex)
    # --debug output only output to dss
    if DEBUG_OUTPUT:
        outputfile = os.path.join(filepath, 'Output', 'weather.dss')
        write_weather_dss(dftmax, dftmin, dfptotal, dfet0, outputfile)
    return(ts_ptotal, ET0Daily)


@numba.jit(nopython=True, cache=True)
def calc_soil_evap(k, Beta1, idates, idays, ts_year, ts_days, start1, dpyAll, ET0, pcp, PcpDaily, EToDaily, OKc):
    PETo = 0.0
    CETo = 0.0
    METo = 0.0
    CEx = 0.0
    CEs = 0.0
    DCT = 1
    dpyAll[0] = 365
    for kk in range(0, idates):
        iy = ts_year[kk] - start1[0] + 1
        id = ts_days[kk]

        if (ET0[k, kk]*100-int(ET0[k, kk]*100)) > 0.5:
            EToDaily[iy, id] = int(ET0[k, kk]*100+1)/100.
        else:
            EToDaily[iy, id] = int(ET0[k, kk]*100)/100.

        if (pcp[k, kk]*100-int(pcp[k, kk]*100)) > 0.5:
            PcpDaily[iy, id] = int(pcp[k, kk]*100+1)/100.
        else:
            PcpDaily[iy, id] = int(pcp[k, kk]*100)/100.

        dpy = 365
        if ts_year[kk] % 4 == 0:
            dpy = 366
        dpyAll[iy] = dpy

        DCT = DCT + 1
        if EToDaily[iy, id] <= 0:
            EToDaily[iy, id] = PETo
        # Cumulative ETo
        CETo = CETo+EToDaily[iy, id]
        METo = CETo/DCT
        if PcpDaily[iy, id] > METo:
            CETo = EToDaily[iy, id]
            DCT = 1
            PEs = 0
        METo = CETo/DCT
        kx = 1.22 - 0.04*METo
        CEx = kx*CETo
        if CEx < 0:
            CEx = 0
        if pow(CEx, 0.5) < Beta1:
            CEs = CEx
        else:
            CEs = Beta1*pow(CEx, 0.5)
        if EToDaily[iy, id] != 0:
            OKc[iy, id] = CEs/CETo
        PETo = EToDaily[iy, id]
        PEs = CEs


@numba.jit(nopython=True, cache=True)
def main_calc_loop(iyears, j, yearType, CBeginDate, CEndDate, Ckc1, Ckc2, Ckc3, CAB, CAC, CAD, NCBeginDate, NCEndDate,
                   NCkc1, NCkc2, NCkc3, NCAB, NCAC, NCAD, kkc1, kkc2, kkc3, EndDate1, BeginDate1, AB1, AC1, AD1,
                   EToMonthly, ETcMonthly, PcpMonthly, ERnMonthly, SpgMonthly, EspgMonthly, MonNetApp, MonDsw, MonDswPos,
                   dpyAll, NumDaysPerMon, erd, Region, k, SWD, OKc, IKc, Kc, EToDaily, PcpDaily, ETcDaily,
                   NA1, BIYear, NA2, LIYear, NA3, osSWDx, SWD0, Dsw, NetApp, NII, NI, osCETc, HAcre,
                   osMCETc, isCETc, isMCETc, osCERn, osCSpg, osMCERn, osMCSpg, isCERn, isPCP, isETaw, isCSpg, isMCERn, isMCSpg,
                   aw1, ADep1, Espg, WSCESpg, WSEspgMonthly, ytemp, DOYtemp, IrrigYear, HAcre_temp, HAcretemp, OKctemp, IKctemp,
                   CCKc, CCKctemp, ETotemp, Kctemp, ETctemp, Pcptemp, Ertemp, Spgtemp, ESpgtemp, Dsw0temp, SWDtemp, SWDxtemp,
                   FC0temp, PWPtemp, SWC0temp, YTDtemp, NAtemp, DOYLIrrig, DOYGrainLIrrig, CPcptemp, CErtemp, CESpgtemp, CETctemp, CDswtemp,
                   HAcreDaily, yDaily, DOY, OKcDaily, IKcDaily, CCKcDaily, EToDaily2, KcDaily, ETcDaily2, PcpDaily2, ErDaily, SpgDaily,
                   ESpgDaily, DswDaily, SWDDaily, SWDxDaily, FCDaily, PWPDaily, SWCDaily, YTDDaily, NADaily, CPcpDaily, CErDaily, CESpgDaily, CETcDaily, CDswDaily,
                   ETAWMonDay, ETAWDaily, MonETAW, idayoutput, dailyunit, dataday, date_index, data2day, data6day, data7day, data8day,
                   data9day, data10day, data11day, data12day, data13day, data14day, data15day, data16day, data17day, data18day, data26day,
                   data27day, data28day, data29day, data30day, data31day, CETAWDaily, data32day, imonthoutput, datamon, data2mon, data3mon,
                   data4mon, data5mon, data6mon, data7mon, data8mon, data9mon, data12mon, data13mon, data14mon, data15mon, data16mon, data17mon,
                   data18mon, data21mon, data10mon, data19mon, data11mon, data20mon,model_start_year):

    nmonth = 0  # pointer for datamon* arrays
    for y in range(1, iyears+1):
        if j == 6 or j == 14:
            yearTypeCal = yearType[y]
        else:
            yearTypeCal = yearType[y-1]

        if yearTypeCal.upper() == "C" or yearTypeCal.upper() == "D":
            BeginDate1 = CBeginDate[j]
            EndDate1 = CEndDate[j]
            kkc1 = Ckc1[j]
            kkc2 = Ckc2[j]
            kkc3 = Ckc3[j]
            AB1 = CAB[j]
            AC1 = CAC[j]
            AD1 = CAD[j]
            # end of critical years
        else:
            BeginDate1 = NCBeginDate[j]
            EndDate1 = NCEndDate[j]
            kkc1 = NCkc1[j]
            kkc2 = NCkc2[j]
            kkc3 = NCkc3[j]
            AB1 = NCAB[j]
            AC1 = NCAC[j]
            AD1 = NCAD[j]
            # end of noncritical years
        if j == 7:
            RiceIrrigEndDate = EndDate1-20
        # Initial Kc values for dates B,C,D,E
        KcB = kkc1
        KcC = kkc2
        KcD = kkc2
        KcE = kkc3

        # ..........................................................
        # .. Calculate day of year corresponding to A,B,C, and E ...
        # .. Note that B, C, D, and E can be bigger than 365   ....
        # .. Note that BB,CC,DD, and EE are always <=365        ...
        # ..........................................................
        if EndDate1 <= BeginDate1:
            EndDate1 = 365+EndDate1
        lenDate = EndDate1-BeginDate1
        B = int(0.01*AB1*lenDate+BeginDate1)
        C = int(0.01*AC1*lenDate+BeginDate1)
        D = int(0.01*AD1*lenDate+BeginDate1)

        BB = B
        if B > 365:
            BB = B-365
        CC = C
        if C > 365:
            CC = C-365
        DD = D
        if D > 365:
            DD = D-365
        EE = EndDate1
        if EndDate1 > 365:
            EE = EndDate1-365

        for Mon in range(0, 12):
            EToMonthly[Mon] = 0
            ETcMonthly[Mon] = 0
            PcpMonthly[Mon] = 0
            ERnMonthly[Mon] = 0
            SpgMonthly[Mon] = 0
            EspgMonthly[Mon] = 0
            MonNetApp[Mon] = 0
            MonDsw[Mon] = 0
            MonDswPos[Mon] = 0
        dpy = dpyAll[y]
        FinalIrrig = 0.0
        Mon = 0
        for ii in range(1, dpy+1):
            # initialize cumulative variable for the first day
            if (y == iyears) and ((y % 4 != 0 and ii > 273) or (y % 4 == 0 and ii > 274)):
                break
            if (y == 1 and ii == 274) or (y == 1 and ii == 1):
                SWD = 0
                PSWD = 0
                NetApp = 0
                CPcp = 0
                CERn = 0
                CESpg = 0
                CDsw = 0
                DCETo = 0
                DCETc = 0
            if Mon < 12:
                # Spg=0.15/NumDaysPerMon[Mon+1]*erd
                if j == 13 or j == 12:  # for native and riparian vegetation
                    Spg = 0.15/NumDaysPerMon[Mon+1]*erd  # 0.025 0.05 0.15
                elif j == 16: # Duck pond
                    Spg = 0.06 / NumDaysPerMon[Mon + 1] * erd
                else:
                    Spg = 0.025/NumDaysPerMon[Mon+1]*erd  # 0.025 0.15
            if Region[k] == 1:
                Spg = 0.0
            PSWD = SWD
            OKc1 = OKc[y, ii]
            IKc1 = IKc[y, ii]
            Kc11 = Kc[y, ii]
            ETo = EToDaily[y, ii]
            PCP = PcpDaily[y, ii]
            # This will reduce the value for Spg when the seepage is greater than the Dsw

            if (SWD+ETcDaily[y, ii]-Spg) >= 0:
                Espg = Spg
            else:
                Espg = SWD+ETcDaily[y, ii]
                if Region[k] == 1:
                    Espg = 0
            # adjust for the rice
            if j == 7 and IKc1 != 0:
                Espg = 0
            # water surface and RV
            # water surface is now cropNum 15
            # if((j==12 or j==15)) Espg=ETcDaily[y,i]
            if j == 12:
                Espg = ETcDaily[y, ii]
                if Region[k] == 1:
                    Espg = 0
            if j == 15:
                Espg = 0

            # if y == 33  and ii>250 and j==1:
            # print y,ii,SWD,ETcDaily[y,ii],Espg,PCP

            if (SWD+ETcDaily[y, ii]-Espg-PCP) >= 0:
                ERn = PCP
            else:
                ERn = SWD + ETcDaily[y, ii]-Espg

            # rice
            if j == 7 and IKc1 != 0:
                ERn = 0
            # if((j==12 or j==15)) ERn=0
            if j == 12:
                ERn = 0
            if j == 15:
                ERn = PCP
            # This calculates the Dsw adjusted for Espg and ERn
            # only for water surface dsw=+espg

            # if (y == 34 or y ==35) and ii<120 and j==1:
            # print y,ii,ETcDaily[y,ii], ERn, Espg, ETcDaily[y,ii]-ERn+Espg

            if j != 15:
                Dsw = ETcDaily[y, ii]-ERn-Espg
            else:  # for water surface
                Dsw = ETcDaily[y, ii]-ERn+Espg
            # Dswp=Dsw

            if ii >= BeginDate1:
                SWDx = NA1[y]

            if ii >= BIYear[y]:
                SWDx = NA2[y]
            if ii >= LIYear[y]:
                SWDx = NA3[y]
            if ii == EE:
                SWDx = NA3[y]
            if EE < BeginDate1 and ii < BIYear[y]:
                SWDx = NA2[y]
            # if j== 1 and k == 0 and (y==1 or y==2):
            # print NA1[y],NA2[y],NA3[y],BeginDate1,BIYear[y],LIYear[y],EE

            # rice-water surface-Riparian
            if j == 7 or j == 12 or j == 15 or j == 16:
                SWDx = osSWDx
            # off season
            if IKc1 == 0:
                SWDx = SWD0
            if IKc1 != 0:
                if (SWD+Dsw) > SWDx:
                    NetApp = SWD+Dsw
                else:
                    NetApp = 0
                # rice - every day we have irrigation except the last 20 days
                if j == 7:
                    NetApp = SWD+Dsw
                    if ii > RiceIrrigEndDate:
                        NetApp = 0
                # ***********add for Native Vegetation, not in the DETAW-UCD***************
                if j == 13:
                    NetApp = 0
                # *************************************************************
            # end of in-season
            if IKc1 == 0:
                NetApp = 0
            SWD = SWD + Dsw - NetApp
            if SWD < 0:
                SWD = 0

            if dpyAll[y] == 366:
                if ii > NII[Mon]:
                    Mon = Mon + 1
            else:
                if ii > NI[Mon]:
                    Mon = Mon + 1

            if IKc1 == 0:
                if SWD0 < osSWDx:
                    SWD0 = osSWDx
                Diff = 0
                if (SWD+Dsw) >= SWD0:
                    Diff = SWD0-SWD
                if (SWD+Dsw) < SWD0:
                    if (y % 4 != 0 and ii > 273) or (y % 4 == 0 and ii > 274):
                        osCETc[y] = osCETc[y]+ETcDaily[y, ii] * \
                            HAcre[k, y, j]*0.0081071
                    else:
                        osCETc[y-1] = osCETc[y-1]+ETcDaily[y, ii] * \
                            HAcre[k, y-1, j]*0.0081071
                    osMCETc[Mon] = osMCETc[Mon]+ETcDaily[y, ii]
                else:
                    if((y % 4 != 0 and ii > 273) or (y % 4 == 0 and ii > 274)):
                        osCETc[y] = osCETc[y]+Diff*HAcre[k, y, j]*0.0081071
                    else:
                        osCETc[y-1] = osCETc[y-1]+Diff * \
                            HAcre[k, y-1, j]*0.0081071
                    osMCETc[Mon] = osMCETc[Mon]+Diff
            else:
                if (y % 4 != 0 and ii > 273) or (y % 4 == 0 and ii > 274):
                    isCETc[y] = isCETc[y]+ETcDaily[y, ii] * \
                        HAcre[k, y, j]*0.0081071
                else:
                    isCETc[y-1] = isCETc[y-1]+ETcDaily[y, ii] * \
                        HAcre[k, y-1, j]*0.0081071
                isMCETc[Mon] = isMCETc[Mon]+ETcDaily[y, ii]

            # ISCERn and OsCERn are in and off-season cum effect rainfall
            if IKc1 == 0:
                if (y % 4 != 0 and ii > 273) or (y % 4 == 0 and ii > 274):
                    osCERn[y] = osCERn[y]+ERn*HAcre[k, y, j]*0.0081071
                    osCSpg[y] = osCSpg[y]+Spg*HAcre[k, y, j]*0.0081071
                else:
                    osCERn[y-1] = osCERn[y-1]+ERn*HAcre[k, y-1, j]*0.0081071
                    osCSpg[y-1] = osCSpg[y-1]+Spg*HAcre[k, y-1, j]*0.0081071
                osMCERn[Mon] = osMCERn[Mon]+ERn
                osMCSpg[Mon] = osMCSpg[Mon]+Spg
            else:
                if (y % 4 != 0 and ii > 273) or (y % 4 == 0 and ii > 274):
                    isCERn[y] = isCERn[y]+ERn*HAcre[k, y, j]*0.0081071
                    isPCP[y] = isPCP[y]+PCP*HAcre[k, y, j]*0.0081071
                    isETaw[y] = isETaw[y]+NetApp*HAcre[k, y, j]*0.0081071

                    isCSpg[y] = isCSpg[y]+Spg*HAcre[k, y, j]*0.0081071
                else:
                    isCERn[y-1] = isCERn[y-1]+ERn*HAcre[k, y-1, j]*0.0081071
                    isPCP[y-1] = isPCP[y-1]+PCP*HAcre[k, y-1, j]*0.0081071
                    isETaw[y-1] = isETaw[y-1]+NetApp*HAcre[k, y-1, j]*0.0081071
                    isCSpg[y-1] = isCSpg[y-1]+Spg*HAcre[k, y-1, j]*0.0081071
                isMCERn[Mon] = isMCERn[Mon]+ERn
                isMCSpg[Mon] = isMCSpg[Mon]+Spg
                # end of Pcp

            # the following 3 lines: no use now written in original code
            # if IKc1 != 0 and flag == "s":  ##off season, the following 3 lines are modififed
            # SWD0=SWD
            ##flag = " "

            MaxSWD = erd*aw1*ADep1/100
            FC = MaxSWD*4
            PWP = MaxSWD*2
            SWC = FC-SWD
            # YTDD is for YTD in column
            YTDD = FC-SWDx
            # calculate yeartypeCal for oct to oct year
            if (y % 4 != 0 and ii >= 274) or (y % 4 == 0 and ii >= 275):
                yearTypeCal = yearType[y]
            else:
                yearTypeCal = yearType[y-1]
            # for first year only print ftom oct 1(day274)
            CPcp = CPcp+PCP
            CERn = CERn+ERn
            CESpg = CESpg+Espg
            # for water surface we use espg* hacre each crop
            if j != 15:
                if(y % 4 != 0 and ii > 273) or (y % 4 == 0 and ii > 274):
                    WSCESpg[y, ii] = WSCESpg[y, ii] + \
                        CESpg*HAcre[k, y, j]*0.0081071
                    WSEspgMonthly[y, Mon] = WSEspgMonthly[y,
                                                          Mon]+Espg*HAcre[k, y, j]*0.0081071
                else:
                    WSCESpg[y, ii] = WSCESpg[y, ii] + \
                        CESpg*HAcre[k, y-1, j]*0.0081071
                    WSEspgMonthly[y, Mon] = WSEspgMonthly[y,
                                                          Mon]+Espg*HAcre[k, y-1, j]*0.0081071
                # before oct 1
            # end of if crop is not water surface
            DCETo = DCETo + EToDaily[y, ii]
            DCETc = DCETc+ETcDaily[y, ii]
            CDsw = CDsw+Dsw

            EToMonthly[Mon] = EToMonthly[Mon]+EToDaily[y, ii]
            ETcMonthly[Mon] = ETcMonthly[Mon] + ETcDaily[y, ii]
            ERnMonthly[Mon] = ERnMonthly[Mon]+ERn
            SpgMonthly[Mon] = SpgMonthly[Mon]+Spg
            EspgMonthly[Mon] = EspgMonthly[Mon]+Espg

            PcpMonthly[Mon] = PcpMonthly[Mon]+PcpDaily[y, ii]
            MonNetApp[Mon] = MonNetApp[Mon]+NetApp
            MonDsw[Mon] = ETcMonthly[Mon] - (ERnMonthly[Mon]+EspgMonthly[Mon])
            # if(MonACETAW[Mon]<0) MonACETAW[Mon]=0
            MonDswPos[Mon] = MonDsw[Mon]
            if MonDsw[Mon] < 0:
                MonDswPos[Mon] = 0
            if(y != 1 and y != iyears) or (y == 1 and ii > 273) or (y == iyears and ii < 274 and y % 4 != 0) or (y == iyears and ii < 275 and y % 4 == 0):
                if j != 15:
                    if j == 14:
                        NetApp = 0
                    if (y % 4 != 0 and ii > 273) or (y % 4 == 0 and ii > 274):
                        HAcre_temp = HAcre[k, y, j]*2.471
                    else:
                        HAcre_temp = HAcre[k, y-1, j]*2.471
                    # crop is not water surface
                else:
                    Spg = 0
                    Espg = 0
                    NetApp = 0
                    if (y % 4 != 0 and ii > 273) or (y % 4 == 0 and ii > 274):
                        HAcre_temp = HAcre[k, y, j]*2.471
                    else:
                        HAcre_temp = HAcre[k, y-1, j]*2.471
                # crop is water surface
                if y == 1 and ii == 274:
                    ik = 1

                ytemp[ik] = y + model_start_year - 1
                DOYtemp[ik] = ii
                if DOYtemp[ik] == 1:
                    IrrigYear = IrrigYear + 1

                HAcretemp[ik] = HAcre_temp
                OKctemp[ik] = OKc[y, ii]
                IKctemp[ik] = IKc[y, ii]
                CCKctemp[ik] = CCKc[ii]
                ETotemp[ik] = EToDaily[y, ii]
                Kctemp[ik] = Kc[y, ii]
                ETctemp[ik] = ETcDaily[y, ii]
                Pcptemp[ik] = PcpDaily[y, ii]

                Ertemp[ik] = ERn
                Spgtemp[ik] = Spg
                ESpgtemp[ik] = Espg
                Dsw0temp[ik] = Dsw
                SWDtemp[ik] = SWD
                SWDxtemp[ik] = SWDx
                FC0temp[ik] = FC
                PWPtemp[ik] = PWP
                SWC0temp[ik] = SWC
                YTDtemp[ik] = YTDD
                NAtemp[ik] = NetApp

                if NAtemp[ik] > 0.000001:
                    DOYLIrrig[IrrigYear] = DOYtemp[ik]
                    if j == 6:
                        if DOYtemp[ik] < 152:
                            DOYGrainLIrrig[IrrigYear] = DOYtemp[ik]
                CPcptemp[ik] = CPcp
                CErtemp[ik] = CERn
                CESpgtemp[ik] = CESpg
                CETctemp[ik] = DCETc
                CDswtemp[ik] = CDsw
                ik = ik + 1
                if ik > dpy:
                    ik = 1

                if y > 1 and ik == 1:
                    ic = slice(1, dpy+1)
                    HAcreDaily[ic] = HAcretemp[ic]
                    yDaily[ic] = ytemp[ic]
                    DOY[ic] = DOYtemp[ic]
                    OKcDaily[ic] = OKctemp[ic]
                    IKcDaily[ic] = IKctemp[ic]
                    CCKcDaily[ic] = CCKctemp[ic]
                    EToDaily2[ic] = ETotemp[ic]
                    KcDaily[ic] = Kctemp[ic]
                    ETcDaily2[ic] = ETctemp[ic]
                    PcpDaily2[ic] = Pcptemp[ic]
                    ErDaily[ic] = Ertemp[ic]
                    SpgDaily[ic] = Spgtemp[ic]
                    ESpgDaily[ic] = ESpgtemp[ic]
                    DswDaily[ic] = Dsw0temp[ic]
                    SWDDaily[ic] = SWDtemp[ic]
                    SWDxDaily[ic] = SWDxtemp[ic]
                    FCDaily[ic] = FC0temp[ic]
                    PWPDaily[ic] = PWPtemp[ic]
                    SWCDaily[ic] = SWC0temp[ic]
                    YTDDaily[ic] = YTDtemp[ic]
                    NADaily[ic] = NAtemp[ic]
                    CPcpDaily[ic] = CPcptemp[ic]
                    CErDaily[ic] = CErtemp[ic]
                    CESpgDaily[ic] = CESpgtemp[ic]
                    CETcDaily[ic] = CETctemp[ic]
                    CDswDaily[ic] = CDswtemp[ic]

                ########Combine OLDHSA and HSA together#################################################

                # if y > 1 and ik==1:
                if (y > 1 and ii == dpy) or (y == iyears and ik == 1):
                    for iq in range(1, dpy+1):  # add 04/29/09
                        # initialize cumulative for start of water year

                        # @if iq ==1:
                        # @    SumDelSWC = 0
                        if yDaily[iq] == model_start_year and iq == 1:
                            FCtemp = int(FCDaily[iq]*10)/10.0
                            ##PSWC = FCDaily[iq]
                            PSWC = FCtemp
                            SumDelSWC = 0
                        SWCtemp = int(SWCDaily[iq]*1000)/1000.0
                        ##DelSWC = PSWC - SWCDaily[iq]
                        DelSWC = PSWC - SWCtemp
                        SumDelSWC = SumDelSWC+DelSWC
                        PSWC = SWCDaily[iq]
                        # do the calculation on the sept 30 of each year  and print
                        yy = yDaily[iq]

                        dpy = 365
                        if yy % 4 == 0:
                            dpy = 366

                        if (yy % 4 != 0 and iq == 365) or (yy % 4 == 0 and iq == 366):
                            MonthlySWC = SumDelSWC/12.0
                            # convert etwawdaily to mon day
                            etMon = 10
                            etDay = 1
                            etdd = calc_etdd(
                                dpy, etMon, yy, NII, NI, DswDaily, MonthlySWC, ETAWMonDay, etDay, NumDaysPerMon)
                            # end of for etdd

                            # calculate aveage ETaw
                            SumNegETAW = 0.0
                            CountNegETAW = 0
                            calc_etawavg(ETAWMonDay, SumNegETAW,
                                         CountNegETAW, NumDaysPerMon)
                            # convert etawdmon day to etaw daily
                            # ETAWDaily[1:dpy+1]=ETAWMonDay[(10:10+12)%12,(273+leap_year(yearCal):)]
                            etMon = 10
                            etDay = 1
                            calc_etaw_month(
                                dpy, etMon, yy, NII, etdd, NI, ETAWMonDay, etDay, ETAWDaily)
                            ##ii = 0
                            CETAWDaily = 0
                            Mon = 10
                            CETAWDaily = calc_etaw_daily(
                                dpy, Mon, yy, DOY, NII, NI, j, DOYLIrrig, yDaily, ETAWDaily, DOYGrainLIrrig, CETAWDaily, MonETAW,model_start_year)
                            # split loop
                            temp_vector = HAcreDaily[1:dpy+1]*0.00328084
                            if idayoutput == 1:
                                if dailyunit == 1:
                                    dataday[date_index:date_index +
                                            dpy] = ETcDaily2[1:dpy+1]*temp_vector
                                    dataday[date_index:date_index +
                                            dpy] = (ETcDaily2[1:dpy+1]*temp_vector)
                                    data2day[date_index:date_index +
                                             dpy] = (PcpDaily2[1:dpy+1]*temp_vector)
                                    data6day[date_index:date_index +
                                             dpy] = (EToDaily2[1:dpy+1]*temp_vector)
                                    data7day[date_index:date_index +
                                             dpy] = (KcDaily[1:dpy+1]*temp_vector)
                                    data8day[date_index:date_index +
                                             dpy] = (SWCDaily[1:dpy+1]*temp_vector)
                                    data9day[date_index:date_index +
                                             dpy] = (SpgDaily[1:dpy+1]*temp_vector)
                                    data10day[date_index:date_index +
                                              dpy] = (ESpgDaily[1:dpy+1]*temp_vector)
                                    data11day[date_index:date_index +
                                              dpy] = (DswDaily[1:dpy+1]*temp_vector)
                                    data12day[date_index:date_index +
                                              dpy] = (ETAWDaily[1:dpy+1]*temp_vector)
                                    data13day[date_index:date_index +
                                              dpy] = (ErDaily[1:dpy+1]*temp_vector)
                                    data14day[date_index:date_index +
                                              dpy] = (SWDxDaily[1:dpy+1]*temp_vector)
                                    data15day[date_index:date_index +
                                              dpy] = (FCDaily[1:dpy+1]*temp_vector)
                                    data16day[date_index:date_index +
                                              dpy] = (PWPDaily[1:dpy+1]*temp_vector)
                                    data17day[date_index:date_index +
                                              dpy] = (SWDDaily[1:dpy+1]*temp_vector)
                                    data18day[date_index:date_index +
                                              dpy] = (YTDDaily[1:dpy+1]*temp_vector)
                                else:
                                    dataday[date_index:date_index +
                                            dpy] = (ETcDaily2[1:dpy+1])
                                    data2day[date_index:date_index +
                                             dpy] = (PcpDaily2[1:dpy+1])
                                    data6day[date_index:date_index +
                                             dpy] = (EToDaily2[1:dpy+1])
                                    data7day[date_index:date_index +
                                             dpy] = (KcDaily[1:dpy+1])
                                    data8day[date_index:date_index +
                                             dpy] = (SWCDaily[1:dpy+1])
                                    data9day[date_index:date_index +
                                             dpy] = (SpgDaily[1:dpy+1])
                                    data10day[date_index:date_index +
                                              dpy] = (ESpgDaily[1:dpy+1])
                                    data11day[date_index:date_index +
                                              dpy] = (DswDaily[1:dpy+1])
                                    data12day[date_index:date_index +
                                              dpy] = (ETAWDaily[1:dpy+1])
                                    data13day[date_index:date_index +
                                              dpy] = (ErDaily[1:dpy+1])
                                    data14day[date_index:date_index +
                                              dpy] = (SWDxDaily[1:dpy+1])
                                    data15day[date_index:date_index +
                                              dpy] = (FCDaily[1:dpy+1])
                                    data16day[date_index:date_index +
                                              dpy] = (PWPDaily[1:dpy+1])
                                    data17day[date_index:date_index +
                                              dpy] = (SWDDaily[1:dpy+1])
                                    data18day[date_index:date_index +
                                              dpy] = (YTDDaily[1:dpy+1])
                                # data19day[date_index:date_index+dpy]=(NADaily[1:dpy+1])
                                # data20day[date_index:date_index+dpy]=(CPcpDaily[1:dpy+1])
                                # data21day[date_index:date_index+dpy]=(CErDaily[1:dpy+1])
                                # data22day[date_index:date_index+dpy]=(CESpgDaily[1:dpy+1])
                                # data23day[date_index:date_index+dpy]=(CETcDaily[1:dpy+1])
                                # data24day[date_index:date_index+dpy]=(CDswDaily[1:dpy+1])
                                # data25day[date_index:date_index+dpy]=(CETAWDaily)
                                data26day[date_index:date_index +
                                          dpy] = (NADaily[1:dpy+1]*temp_vector)
                                data27day[date_index:date_index +
                                          dpy] = (CPcpDaily[1:dpy+1]*temp_vector)
                                data28day[date_index:date_index +
                                          dpy] = (CErDaily[1:dpy+1]*temp_vector)
                                data29day[date_index:date_index +
                                          dpy] = (CESpgDaily[1:dpy+1]*temp_vector)
                                data30day[date_index:date_index +
                                          dpy] = (CETcDaily[1:dpy+1]*temp_vector)
                                data31day[date_index:date_index +
                                          dpy] = (CDswDaily[1:dpy+1]*temp_vector)
                                data32day[date_index:date_index +
                                          dpy] = (CETAWDaily*temp_vector)
                                date_index = date_index+dpy
                                SumDelSWC = 0
                if (y % 4 != 0 and ii == 273) or (y % 4 == 0 and ii == 274):
                    NetApp = 0
                    CPcp = 0
                    CERn = 0
                    CESpg = 0
                    DCETc = 0
                    CDsw = 0

                if ii == dpy or (y == iyears and ii == 273 and y % 4 != 0) or (y == iyears and ii == 274 and y % 4 == 0):
                    for Mon in range(0, 12):
                        if(y != 1 and y != iyears) or (y == 1 and Mon >= 9) or (y == iyears and Mon < 9):

                            if j != 15:
                                # NetApp = 0 for non-irrig Grain
                                if j == 14:
                                    MonNetApp[Mon] = 0
                                if Mon >= 9:
                                    temp_scalar = HAcre[k, y, j]*0.0081071
                                else:
                                    temp_scalar = HAcre[k, y-1, j]*0.0081071
                            else:
                                MonNetApp[Mon] = 0
                                SpgMonthly[Mon] = 0
                                EspgMonthly[Mon] = 0
                                if Mon >= 9:
                                    temp_scalar = HAcre[k, y, j]*0.0081071
                                else:
                                    temp_scalar = HAcre[k, y-1, j]*0.0081071
                            if imonthoutput == 1:
                                nmonth = nmonth+1
                                datamon[nmonth] = (ETcMonthly[Mon])
                                data2mon[nmonth] = (EToMonthly[Mon])
                                data3mon[nmonth] = (MonNetApp[Mon])
                                data4mon[nmonth] = (PcpMonthly[Mon])
                                data5mon[nmonth] = (ERnMonthly[Mon])
                                data6mon[nmonth] = (SpgMonthly[Mon])
                                data7mon[nmonth] = (EspgMonthly[Mon])
                                data8mon[nmonth] = (MonDsw[Mon])
                                data9mon[nmonth] = (MonDswPos[Mon])
                                ##data10mon[nmonth] = (MonETAW[y,Mon])
                                ##data11mon[nmonth] = (MonETAWPos)
                                data12mon[nmonth] = (
                                    MonNetApp[Mon]*temp_scalar)
                                data13mon[nmonth] = (
                                    PcpMonthly[Mon]*temp_scalar)
                                data14mon[nmonth] = (
                                    ERnMonthly[Mon]*temp_scalar)
                                if j != 15:
                                    data15mon[nmonth] = (
                                        EspgMonthly[Mon]*temp_scalar)
                                else:
                                    data15mon[nmonth] = (WSEspgMonthly[y, Mon])
                                data16mon[nmonth] = (
                                    ETcMonthly[Mon]*temp_scalar)
                                data17mon[nmonth] = (MonDsw[Mon]*temp_scalar)
                                data18mon[nmonth] = (
                                    MonDswPos[Mon]*temp_scalar)
                                ##data19mon[nmonth] = (MonETAW[y,Mon]*temp_scalar)
                                ##data20mon[nmonth] = (MonETAWPos*temp_scalar)
                                data21mon[nmonth] = (
                                    EToMonthly[Mon]*temp_scalar)

                            if y > 1:
                                if Mon == 0:
                                    for imtemp in range(10, 13):
                                        data10mon[nmonth] = (
                                            MonETAW[y-1, imtemp])
                                        data19mon[nmonth] = (
                                            MonETAW[y-1, imtemp]*temp_scalar)
                                        # 5/1/09 revised
                                        MonETAWPos = MonETAW[y-1, imtemp]
                                        if MonETAWPos < 0:
                                            MonETAWPos = 0
                                        if imonthoutput == 1:
                                            data11mon[nmonth] = (MonETAWPos)
                                            data20mon[nmonth] = (
                                                MonETAWPos*temp_scalar)
                                if Mon < 9:
                                    # 5/1/09 revised
                                    MonETAWPos = MonETAW[y, Mon+1]
                                    if MonETAWPos < 0:
                                        MonETAWPos = 0
                                    if imonthoutput == 1:
                                        if j == 15 and k == 0 and MonETAW[y, Mon+1] > 0.0:
                                            print(y, Mon+1, " MonETAW=",
                                                  MonETAW[y, Mon+1])
                                        data10mon[nmonth] = (MonETAW[y, Mon+1])
                                        data19mon[nmonth] = (
                                            MonETAW[y, Mon+1]*temp_scalar)
                                        data11mon[nmonth] = (MonETAWPos)
                                        data20mon[nmonth] = (
                                            MonETAWPos*temp_scalar)


@numba.jit(nopython=True, cache=True)
def calc_kc_vals(iyears, yearType, CCropType, j, CBeginDate, CEndDate, Cf, Ckc1, Ckc2, Ckc3, CAB, CAC, CAD, CSDx, CRDxU, CRDxL, CawL, CawU, CADep, Region, k, NCCropType, NCBeginDate, NCEndDate, NCf, NCkc1, NCkc2, NCkc3, NCAB, NCAC, NCAD, NCSDx, NCRDxU, NCRDxL, NCawL, NCawU, NCADep, osIkc, OKc, EToDaily, IGETo, osFkc, isIkc, Beta1):
    for i in range(0, iyears):
        SWD = 0.0
        PSWD = 0.0
        if yearType[i].upper() == "C" or yearType[i].upper() == "D":
            CropType1 = CCropType[j]
            BeginDate1 = CBeginDate[j]
            EndDate1 = CEndDate[j]
            f1 = Cf[j]
            kkc1 = Ckc1[j]
            kkc2 = Ckc2[j]
            kkc3 = Ckc3[j]
            AB1 = CAB[j]
            AC1 = CAC[j]
            AD1 = CAD[j]
            SDx1 = CSDx[j]
            CRDxU1 = CRDxU[j]
            CRDxL1 = CRDxL[j]
            CawL1 = CawL[j]
            CawU1 = CawU[j]
            ADep1 = CADep[j]
            # Set RDX1 to RDxU or RDxL
            if Region[k] == 0:
                RDx1 = CRDxL1
                aw1 = CawL1
            else:
                RDx1 = CRDxU1
                aw1 = CawU1
        else:
            CropType1 = NCCropType[j]
            BeginDate1 = NCBeginDate[j]
            EndDate1 = NCEndDate[j]
            f1 = NCf[j]
            kkc1 = NCkc1[j]
            kkc2 = NCkc2[j]
            kkc3 = NCkc3[j]
            AB1 = NCAB[j]
            AC1 = NCAC[j]
            AD1 = NCAD[j]
            SDx1 = NCSDx[j]
            NCRDxU1 = NCRDxU[j]
            NCRDxL1 = NCRDxL[j]
            NCawL1 = NCawL[j]
            NCawU1 = NCawU[j]
            ADep1 = NCADep[j]
            # Set RDX1 to RDxU or RDxL
            if Region[k] == 0:
                RDx1 = NCRDxL1
                aw1 = NCawL1
            else:
                RDx1 = NCRDxU1
                aw1 = NCawU1
        # set the 20 days before end date for rice
        if j == 7:
            RiceIrrigEndDate = EndDate1-20
        # set the default irrigation frequency to 30
        if RDx1 < SDx1:
            erd = RDx1
        else:
            erd = SDx1
        if f1 == 0:
            f1 = 30
        PAW1 = erd*aw1
        # YTD=PAW*(ADep/100)
        YTD = PAW1*(ADep1/100)
        osSWDx = aw1*300*ADep1/100  # off season max soil water depletion
        # Initilal Kc values for dates B,C,D,E
        KcB = kkc1
        KcC = kkc2
        KcD = kkc2
        KcE = kkc3
        # ..........................................................
        # .. Calculate day of year corresponding to A,B,C, and E ...
        # .. Note that B, C, D, and E can be bigger than 365   ....
        # .. Note that BB,CC,DD, and EE are always <=365        ...
        # ..........................................................
        if EndDate1 <= BeginDate1:
            EndDate1 = 365+EndDate1
        lenDate = EndDate1-BeginDate1
        B = int(0.01*AB1*lenDate+BeginDate1)
        C = int(0.01*AC1*lenDate+BeginDate1)
        D = int(0.01*AD1*lenDate+BeginDate1)

        BB = B
        if B > 365:
            BB = B-365
        CC = C
        if C > 365:
            CC = C-365
        DD = D
        if D > 365:
            DD = D-365
        EE = EndDate1
        if EndDate1 > 365:
            EE = EndDate1-365
        # ...............................................................
        # The following loop calculates Kc values
        # ...............................................................
        # for y in range(1, iyears+1):
        # Jan29-2007
        # starting date of grain and Non-irrig Greain is october
        y = i+1
        if j == 6 or j == 14:
            yearTypeCal = yearType[y]
        else:
            yearTypeCal = yearType[y-1]
        # critical year
        if yearTypeCal.upper() == "C" or yearTypeCal.upper() == "D":
            BeginDate1 = CBeginDate[j]
            EndDate1 = CEndDate[j]
            kkc1 = Ckc1[j]
            kkc2 = Ckc2[j]
            kkc3 = Ckc3[j]
            AB1 = CAB[j]
            AC1 = CAC[j]
            AD1 = CAD[j]
        else:
            BeginDate1 = NCBeginDate[j]
            EndDate1 = NCEndDate[j]
            kkc1 = NCkc1[j]
            kkc2 = NCkc2[j]
            kkc3 = NCkc3[j]
            AB1 = NCAB[j]
            AC1 = NCAC[j]
            AD1 = NCAD[j]

        if j == 7:
            RiceIrrigEndDate = EndDate1-20
        # Initial Kc values for dates B,C,D,E
        KcB = kkc1
        KcC = kkc2
        KcD = kkc2
        KcE = kkc3
        # ..........................................................
        # .. Calculate day of year corresponding to A,B,C, and E ...
        # .. Note that B, C, D, and E can be bigger than 365   ....
        # .. Note that BB,CC,DD, and EE are always <=365        ...
        if EndDate1 <= BeginDate1:
            EndDate1 = 365+EndDate1
        lenDate = EndDate1-BeginDate1
        B = int(0.01*AB1*lenDate+BeginDate1)
        C = int(0.01*AC1*lenDate+BeginDate1)
        D = int(0.01*AD1*lenDate+BeginDate1)

        BB = B
        if B > 365:
            BB = B-365
        CC = C
        if C > 365:
            CC = C-365
        DD = D
        if D > 365:
            DD = D-365
        EE = EndDate1
        if EndDate1 > 365:
            EE = EndDate1-365

        # end of jan 29-2007
        osIkc[y] = 0
        ctn = 0
        IGKc = 0
        IGETo1 = 0
        Q = BeginDate1
        R = B
        # because in this program I use int number for croptypw I change the
        # condition from cropType>2 to cropType>=2
        if CropType1 > 2:
            R = B+10
        ctn = R+1-Q
        IGKc = numpy.sum(OKc[y, Q:R+1])
        IGETo1 = numpy.sum(EToDaily[y, Q:R+1])
        osIkc[y] = IGKc/ctn  # initial growth Kc from off-season
        IGETo[y] = IGETo1/ctn  # initial growth mean Eto rate
        # Identify osFkc for Kc on date E
        osFkc[y] = 0
        ctn = 0
        EKc = 0
        R = EndDate1
        if EndDate1 > 365:
            R = EndDate1-365
        if R == 1 and EndDate1 >= 365:
            R = 365
        Q = R-10
        ctn = 11
        EKc = numpy.sum(OKc[y, R-10:R+1])
        osFkc[y] = EKc/ctn  # final Kc on date E from off-season
        # .. Identify initial growth Kc from irrig freq(F)
        isIkc[y] = 0
        CETo = f1*IGETo[y]

        # RichEdit1->Lines->Add("CETo="+FloatToStr(CETo))
        # kx=1.05-0.03*IGETo[y]
        kx = 1.22-0.04*IGETo[y]
        CEx = kx*CETo
        if pow(CEx, 0.5) < Beta1:
            CEs = CEx
        else:
            CEs = Beta1*pow(CEx, 0.5)

        if CETo != 0:
            isIkc[y] = CEs/CETo

        if CropType1 > 1:
            isIkc[y] = 0
    # end of for  y<iyears  ********************************
    return kkc1, kkc2, kkc3, EndDate1, BeginDate1, AB1, AC1, AD1, CropType1, y, f1, YTD, osSWDx, erd, aw1, ADep1


@numba.jit(nopython=True, cache=True)
def calc_kc_daily(iyears, j, yearType, CBeginDate, CEndDate, Ckc1, Ckc2, Ckc3, CAB, CAC, CAD, NCBeginDate, NCEndDate, NCkc1, NCkc2, NCkc3, NCAB, NCAC, NCAD, kkc1, kkc2, kkc3, EndDate1, BeginDate1, AB1, AC1, AD1, CropType1, LowIkc, isIkc, KcByr, y, KcCyr, KcDyr, KcEyr):
    for y in range(1, iyears+1):
        # Jan29-2007

        # starting date of grain and Non-irrig Greain is october
        if j == 6 or j == 14:
            yearTypeCal = yearType[y]
        else:
            yearTypeCal = yearType[y-1]

        # critical years
        if yearTypeCal.upper() == "C" or yearTypeCal.upper() == "D":
            BeginDate1 = CBeginDate[j]
            EndDate1 = CEndDate[j]
            kkc1 = Ckc1[j]
            kkc2 = Ckc2[j]
            kkc3 = Ckc3[j]
            AB1 = CAB[j]
            AC1 = CAC[j]
            AD1 = CAD[j]
            # end of critical years
        else:  # non-critical yeras
            BeginDate1 = NCBeginDate[j]
            EndDate1 = NCEndDate[j]
            kkc1 = NCkc1[j]
            kkc2 = NCkc2[j]
            kkc3 = NCkc3[j]
            AB1 = NCAB[j]
            AC1 = NCAC[j]
            AD1 = NCAD[j]
            # we have to set RDX1 to RDxU or RDxL
            # end of non-critical yeras

        # set the 20 days before end date for rice
        if j == 7:
            RiceIrrigEndDate = EndDate1-20

        # Initilal Kc values for dates B,C,D,E
        KcB = kkc1
        KcC = kkc2
        KcD = kkc2
        KcE = kkc3

        # ..........................................................
        # .. Calculate day of year corresponding to A,B,C, and E ...
        # .. Note that B, C, D, and E can be bigger than 365   ....
        # .. Note that BB,CC,DD, and EE are always <=365        ...
        # ..........................................................
        if EndDate1 <= BeginDate1:
            EndDate1 = 365+EndDate1
        lenDate = EndDate1-BeginDate1
        B = int(0.01*AB1*lenDate+BeginDate1)
        C = int(0.01*AC1*lenDate+BeginDate1)
        D = int(0.01*AD1*lenDate+BeginDate1)

        BB = B
        if B > 365:
            BB = B-365
        CC = C
        if C > 365:
            CC = C-365
        DD = D
        if D > 365:
            DD = D-365
        EE = EndDate1
        if EndDate1 > 365:
            EE = EndDate1-365

        if CropType1 > 1:
            if CropType1 > 2:
                if CropType1 <= 3:
                    KcB = LowIkc
                    KcE = kkc3
            # end of if CropType>2
            else:  # CropType>2
                KcB = kkc1
                KcC = kkc2
                KcD = kkc2
                KcE = kkc3
                # RichEdit1->Lines->Add("SaraKcBKcc")
            # end of CropType >2
        # end of if croptype>1
        else:  # CropType>1
            KcB = LowIkc
            KcC = kkc2
            KcD = kkc2
            KcE = kkc3

            if isIkc[y] > KcB:
                KcB = isIkc[y]

        # end of else CropType >1
        KcByr[y] = KcB
        KcCyr[y] = KcC
        KcDyr[y] = KcD
        KcEyr[y] = KcE
    # end of for y<YCt
    return kkc1, kkc2, kkc3, EndDate1, BeginDate1, AB1, AC1, AD1, y


@numba.jit(nopython=True, cache=True)
def calc_ikc_daily(iyears, j, yearType, CBeginDate, CEndDate, Ckc1, Ckc2, Ckc3, CAB, CAC, CAD, NCBeginDate, NCEndDate, NCkc1, NCkc2, NCkc3, NCAB, NCAC, NCAD, kkc1, kkc2, kkc3, EndDate1, BeginDate1, AB1, AC1, AD1, IKc, y, idays, dpyAll, BeginDateYear, IKc1, OKc, CropType1, CCKc, CC1, Kc, EToDaily, ETcDaily, f1, BIYear, NA1, YTD, LIYear, NA2, NA3):
    for y in range(1, iyears+1):

        # Jan29-2007

        # starting date of grain and Non-irrig Greain is october
        if j == 6 or j == 14:
            yearTypeCal = yearType[y]
        else:
            yearTypeCal = yearType[y-1]
        # critical years
        if yearTypeCal.upper() == "C" or yearTypeCal.upper() == "D":
            BeginDate1 = CBeginDate[j]
            EndDate1 = CEndDate[j]
            kkc1 = Ckc1[j]
            kkc2 = Ckc2[j]
            kkc3 = Ckc3[j]
            AB1 = CAB[j]
            AC1 = CAC[j]
            AD1 = CAD[j]
            # end of critical years
        else:  # non-critical yeras
            BeginDate1 = NCBeginDate[j]
            EndDate1 = NCEndDate[j]
            kkc1 = NCkc1[j]
            kkc2 = NCkc2[j]
            kkc3 = NCkc3[j]
            AB1 = NCAB[j]
            AC1 = NCAC[j]
            AD1 = NCAD[j]
            # we have to set RDX1 to RDxU or RDxL
            # end of non-critical yeras

        # set the 20 days before end date for rice
        if j == 7:
            RiceIrrigEndDate = EndDate1-20

        # Initilal Kc values for dates B,C,D,E
        KcB = kkc1
        KcC = kkc2
        KcD = kkc2
        KcE = kkc3

        # ..........................................................
        # .. Calculate day of year corresponding to A,B,C, and E ...
        # .. Note that B, C, D, and E can be bigger than 365   ....
        # .. Note that BB,CC,DD, and EE are always <=365        ...
        # ..........................................................
        if EndDate1 <= BeginDate1:
            EndDate1 = 365+EndDate1
        lenDate = EndDate1-BeginDate1
        B = int(0.01*AB1*lenDate+BeginDate1)
        C = int(0.01*AC1*lenDate+BeginDate1)
        D = int(0.01*AD1*lenDate+BeginDate1)

        BB = B
        if B > 365:
            BB = B-365
        CC = C
        if C > 365:
            CC = C-365
        DD = D
        if D > 365:
            DD = D-365
        EE = EndDate1
        if EndDate1 > 365:
            EE = EndDate1-365
        # end of jan 29-2007
        IKc[y, 0:idays+1] = 0
        dpy = dpyAll[y]
        # KcE=KcEyr[y]
        KcB = kkc1
        KcC = kkc2
        KcD = kkc2
        KcE = kkc3
        BCslope = (KcC-KcB)/(C-B)
        CDslope = (KcD-KcC)/(D-C)
        DEslope = (KcE-KcD)/(EndDate1-D)
        BeginDateYear[y] = BeginDate1  # project doesn't have this line

        #? loop eliminate ?#
        ii = 0
        IKc1, ii, jj = calc_ikc(BeginDate1, EndDate1, dpy, D, C,
                                B, KcB, BCslope, IKc1, CDslope, DEslope, IKc, y, ii)

        if dpyAll[y] == 366:
            IKc[y, 366] = IKc[y, 365]

        # This section calculates daily Kc and print the results....
        # apply startig and ending dates of Stress factor

        # IKcs = KcD*ks   ##Ikcs max kc with stress coeficcient based on date D
            # RichEdit1->Lines->Add(dpyAll[y])
        #? loop eliminate ?#
        calc_etc_daily(dpyAll, y, IKc, OKc, CropType1, CCKc,
                       CC1, IKc1, j, Kc, ii, EToDaily, ETcDaily)

        # end of for i<=dpyAll

        # ....... Determine application number and amounts
        CETc = 0
        #? loop eliminate ?#
        CETc = calc_cetc(BeginDate1, EndDate1, dpyAll, y,
                         EToDaily, jj, Kc, ETcDaily, CETc)

        ET1 = 0
        NumI1 = 0
        # if (PIrr.upperCase()=="Y") NumI1=1
        # modified on jan 8
        # NumI=floor((B-A)/F)
        NumI = int((B-BeginDate1)/f1)
        NumI1 = NumI1+NumI
        BI = BeginDate1+f1*NumI

        if NumI1 == 0:
            BI = B
        BIYear[y] = BI
        #? eliminate loop ?#
        ET1, jj = calc_ET1(BeginDate1, BI, dpyAll, y, ET1, ETcDaily, jj)

        if NumI1 != 0:
            NA1[y] = ET1/NumI1
        ET2 = 0
        #? eliminate loop ?#
        calc_ET2(BI, EndDate1, dpyAll, y, ET2, ETcDaily, jj,
                 ET1, CETc, YTD, LIYear, NA2, NumI1, NA1, NA3)
        ##ii = EndDate1 + 1
        # to get out of loop i<=E
        # end of for i<=E
    # end of y<YCI


@numba.jit(nopython=True, cache=True)
def calc_ET1(BeginDate1, BI, dpyAll, y, ET1, ETcDaily, jj):
    #? eliminate loop ?#
    for ii in range(BeginDate1, BI+1):
        jj = ii
        if ii > dpyAll[y-1]:
            jj = ii-dpyAll[y-1]
        # ET1=ET1+EToDaily[y,j]*Kc[y,j]
        ET1 = ET1+ETcDaily[y, jj]
    return ET1, jj


@numba.jit(nopython=True, cache=True)
def calc_ikc(BeginDate1, EndDate1, dpy, D, C, B, KcB, BCslope, IKc1, CDslope, DEslope, IKc, y, ii):
    for jj in range(BeginDate1, EndDate1+1):
        ii = jj
        if jj > dpy:
            ii = jj-dpy
        if jj <= D:
            if jj <= C:
                if jj <= B:
                    IKc1 = KcB
                else:
                    IKc1 = IKc1+BCslope
            else:
                IKc1 = IKc1 + CDslope
        else:
            IKc1 = IKc1 + DEslope

        IKc[y, ii] = IKc1
    return IKc1, ii, jj


@numba.jit(nopython=True, cache=True)
def calc_ET2(BI, EndDate1, dpyAll, y, ET2, ETcDaily, jj, ET1, CETc, YTD, LIYear, NA2, NumI1, NA1, NA3):
    #? eliminate loop ?#
    for ii in range(BI+1, EndDate1+1):
        jj = ii
        if ii > dpyAll[y-1]:
            jj = ii-dpyAll[y-1]

        ET2 = ET2 + ETcDaily[y, jj]
        if (ET2+ET1) > (CETc-YTD):
            LI = ii
            LIYear[y] = LI
            NI2 = int(ET2/YTD)+1
            NA2[y] = ET2/NI2
            if NumI1 == 0:
                NA1[y] = NA2[y]
            NA3[y] = YTD
            break
            ##ii = EndDate1 + 1
            # to get out of loop i<=E
    # end of for i<=E


@numba.jit(nopython=True, cache=True)
def calc_cetc(BeginDate1, EndDate1, dpyAll, y, EToDaily, jj, Kc, ETcDaily, CETc):
    #? loop eliminate ?#
    for ii in range(BeginDate1, EndDate1+1):
        jj = ii
        if ii > dpyAll[y-1]:
            jj = ii-dpyAll[y]
        # modified for test on oct19-2005
        # CETc=CETc+ETodaily[y,j]*Kc[y,j]
        ETcDaily[y, jj] = EToDaily[y, jj]*Kc[y, jj]
        CETc = CETc+ETcDaily[y, jj]
    return CETc


@numba.jit(nopython=True, cache=True)
def calc_etc_daily(dpyAll, y, IKc, OKc, CropType1, CCKc, CC1, IKc1, j, Kc, ii, EToDaily, ETcDaily):
    for ii in range(1, dpyAll[y]+1):

        IKc1 = IKc[y, ii]
        OKc1 = OKc[y, ii]

        Kc11 = IKc1
        if OKc1 < Kc11:
            Kc11 = OKc1

        if CropType1 >= 3:
            ##Kc11 = Kc11 + CCKc[ii]
            if IKc[y, ii] != 0:
                if Kc11 > 1.15:
                    Kc11 = 1.15
            else:
                if Kc11 > 1.05:
                    Kc11 = 1.05

            if CCKc[ii] == CC1 and Kc11 < 0.9:
                Kc11 = 0.9  # with Cover crop min Kc=0.9
            # apply stress factor and start and end date to crop number 3

            # if CropType1 == 3:
                # if ii>KsStartDate and ii<ksEndDate:
                ##    Kc11 = Kc11*ks
        # End of Kc11 adjustment for type 3 and 4
        Kc11 = IKc1
        if OKc1 > Kc11:
            Kc11 = OKc1

        # May29-2007
        if j == 15:
            Kc11 = 1.1
        Kc[y, ii] = Kc11
        ETcDaily[y, ii] = EToDaily[y, ii]*Kc[y, ii]

    # end of for i<=dpyAll


@numba.jit(nopython=True, cache=True)
def calc_etdd(dpy, etMon, yy, NII, NI, DswDaily, MonthlySWC, ETAWMonDay, etDay, NumDay):
    for etdd in range(1, dpy+1):
        if etMon < 10:
            yearCal = yy
        else:
            yearCal = yy-1
        if yearCal % 4 == 0:
            etDOY = etdd+274
            if etDOY > 366:
                etDOY = etDOY-366
        if yearCal % 4 != 0:
            etDOY = etdd+273
            if etDOY > 365:
                etDOY = etDOY-365

        if yearCal % 4 == 0:
            if etDOY > NII[etMon-1]:
                etMon = etMon + 1
                etDay = 1
        else:
            if etDOY > NI[etMon-1]:
                etMon = etMon + 1
                etDay = 1

        # ETAW=DswDaily[dd]-MonthlySWC/NumDay[Mon]
        Dswtemp = int(DswDaily[etdd]*1000)/1000.0
        ETAWMonDay[etMon, etDay] = Dswtemp-MonthlySWC/NumDay[etMon]
        # @ETAWMonDay[etMon,etDay]=Dswtemp-abs(MonthlySWC/NumDay[etMon])

        etDay = etDay+1
        if (yearCal % 4 != 0 and etDOY == 365) or (yearCal % 4 == 0 and etDOY == 366):
            etMon = 1
            etDay = 1
    # end of for etdd
    return etdd


@numba.jit(nopython=True, cache=True)
def calc_etawavg(ETAWMonDay, SumNegETAW, CountNegETAW, NumDay):
    for etMon in range(1, 13):
        flagETAW = "true"
        while flagETAW == "true":
            for etDay in range(1, NumDay[etMon]+1):
                if ETAWMonDay[etMon, etDay] < 0:
                    SumNegETAW = SumNegETAW+ETAWMonDay[etMon, etDay]
                    CountNegETAW = CountNegETAW+1
                    ETAWMonDay[etMon, etDay] = 0
            if abs(SumNegETAW) > 0:
                flagETAW = "true"
            else:
                flagETAW = "false"
            if NumDay[etMon] > CountNegETAW:
                avgNeg = SumNegETAW/(NumDay[etMon]-CountNegETAW)
            SumNegETAW = 0.0
            CountNegETAW = 0
            for etDay in range(1, NumDay[etMon]+1):
                if ETAWMonDay[etMon, etDay] > 0:
                    ETAWMonDay[etMon, etDay] = ETAWMonDay[etMon, etDay]+avgNeg


@numba.jit(nopython=True, cache=True)
def calc_etaw_month(dpy, etMon, yy, NII, etdd, NI, ETAWMonDay, etDay, ETAWDaily):
    for etdd in range(1, dpy+1):
        if etMon < 10:
            yearCal = yy
        else:
            yearCal = yy - 1
        if yearCal % 4 == 0:
            etDOY = etdd + 274
            if etDOY > 366:
                etDOY = etDOY - 366
            if etDOY > NII[etMon-1]:
                etMon = etMon + 1
                etDay = 1
        else:
            etDOY = etdd + 273
            if etDOY > 365:
                etDOY = etDOY - 365
            if etDOY > NI[etMon-1]:
                etMon = etMon + 1
                etDay = 1
        ETAWDaily[etdd] = ETAWMonDay[etMon, etDay]
        etDay = etDay + 1
        if (yearCal % 4 != 0 and etDOY == 365) or (yearCal % 4 == 0 and etDOY == 366):
            etMon = 1
            etDay = 1


@numba.jit(nopython=True, cache=True)
def calc_etaw_daily(dpy, Mon, yy, DOY, NII, NI, j, DOYLIrrig, yDaily, ETAWDaily, DOYGrainLIrrig, CETAWDaily, MonETAW,model_start_year):
    for dd in range(1, dpy+1):
        if Mon < 10:
            yearCal = yy
        else:
            yearCal = yy-1
        if yearCal % 4 == 0:
            if DOY[dd] > NII[Mon-1]:
                Mon = Mon+1
                ##MonETAW[yearCal-1920,Mon]= 0.0
        else:
            if DOY[dd] > NI[Mon-1]:
                Mon = Mon + 1
                ##MonETAW[yearCal-1920,Mon]= 0.0
        if (yearCal % 4 != 0 and DOY[dd] == 365) or (yearCal % 4 == 0 and DOY[dd] == 366):
            Mon = 1
        if(j != 6):
            if DOY[dd] >= DOYLIrrig[yDaily[dd]-model_start_year]:
                ETAWDaily[dd] = 0
        else:
            if DOYLIrrig[yDaily[dd]-model_start_year] < 275:
                if DOY[dd] >= DOYGrainLIrrig[yDaily[dd]-model_start_year] and DOY[dd] < 275:
                    ETAWDaily[dd] = 0
            else:
                if DOY[dd] >= DOYLIrrig[yDaily[dd]-model_start_year]:
                    ETAWDaily[dd] = 0
        CETAWDaily = CETAWDaily + ETAWDaily[dd]
        # y is for sep 30 each year, so for previous oct-Dec we subtract y by 1
        if Mon < 10:
            MonETAW[yy-model_start_year-1, Mon] = MonETAW[yy-model_start_year-1, Mon] + ETAWDaily[dd]
        else:
            MonETAW[yy-model_start_year, Mon] = MonETAW[yy-model_start_year, Mon] + ETAWDaily[dd]
    return CETAWDaily
# _______________________________________________________________________________
# _______________________________________________________________________________


def historicalETAW(ts_per, ETo_corrector, Region, pcp, ET0, tmax, tmin, ilands, idates, isites, ts_year, ts_mon,
                   ts_days, start1, filepath, NI, NII, NumDaysPerMon, iyears, idayoutput, imonthoutput,
                   iyearoutput, itotaloutput, dailyunit, forDSM2_daily, streamlinemodel,model_start_year,yearType,HAcre,icroptype,crdf,ncrdf):
    InpHSACrop = "  "
    SACropDaily = "  "
    HSACropDailyMean = "  "
    HSACropMonMean = "  "
    Date = "  "
    cpartt = "  "
    monthname = ["JAN", "FEB", "MAR", "APR", "MAY",
                 "JUN", "JUL", "AUG", "SEP", "OCT", "NOV", "DEC"]
    # icroptype is now read in from read_landuse can delete this commented block
    # #? why not read these from input file ?#
    # CropName = ["Urban", "Irrig pasture", "Alfalfa", "All Field", "Sugar beets",
    #             "Irrig Grain", "Rice", "Truck Crops", "Tomato", "Orchard",
    #             "Vineyard", "Riparian Vegetation", "Native Vegetation",
    #             "Non-irrig Grain", "Water Surface"]
    # icroptype = len(CropName)
    idays = 366
    imonths = 12
    # pyhecdss.set_message_level(0)
    # pyhecdss.set_program_name('DETAW')
    ctype = "INST-VAL"  # "PER-AVER"
    iplan = int(0)
    # FIXME start2 is different from start1. Tried to make it the same as start1,
    # but outputs were not the same. My guess is due to numpy indexing issues.
    start2 = [model_start_year,10,1,23,0]
    startdate = str(start2[2])+monthname[start2[1]-1]+str(start2[0])
    starttime = str(start2[3])+"00"
    dpyAll = zeros((iyears+2), int)
    Beta1 = 2.6
    IKc = zeros((iyears+2, idays+1), float)
    Kc = zeros((iyears+2, idays+1), float)

    OKc = zeros((iyears+2, idays+1), float)
    EToDaily = zeros((iyears+2, idays+1), float)
    PcpDaily = zeros((iyears+2, idays+1), float)
    _OKc = zeros((iyears+2, idays+1), float)
    _EToDaily = zeros((iyears+2, idays+1), float)
    _PcpDaily = zeros((iyears+2, idays+1), float)

    ETcDaily = zeros((iyears+2, idays+1), float)
    Date1Daily = "  "*(iyears+2)*(idays+1)
    WSCESpg = zeros((iyears+2, idays+1), float)

    # HAcre = zeros((ilands+1, iyears+2, icroptype+1), float)
    Kc1 = zeros((icroptype+1), float)
    Kc2 = zeros((icroptype+1), float)
    Kc3 = zeros((icroptype+1), float)
    AB = zeros((icroptype+1), float)
    AC = zeros((icroptype+1), float)
    AD = zeros((icroptype+1), float)
    SDx = zeros((icroptype+1), float)
    RDxU = zeros((icroptype+1), float)
    RDxL = zeros((icroptype+1), float)
    awL = zeros((icroptype+1), float)
    awU = zeros((icroptype+1), float)
    ADep = zeros((icroptype+1), float)
    # crop info for critical years
    CBeginDate = zeros((icroptype+1), int)
    CEndDate = zeros((icroptype+1), int)
    CCropType = zeros((icroptype+1), int)
    Cf = zeros((icroptype+1), int)
    Ckc1 = zeros((icroptype+1), float)
    Ckc2 = zeros((icroptype+1), float)
    Ckc3 = zeros((icroptype+1), float)
    CAB = zeros((icroptype+1), float)
    CAC = zeros((icroptype+1), float)
    CAD = zeros((icroptype+1), float)
    CSDx = zeros((icroptype+1), float)
    CRDxU = zeros((icroptype+1), float)
    CRDxL = zeros((icroptype+1), float)
    CawL = zeros((icroptype+1), float)
    CawU = zeros((icroptype+1), float)
    CADep = zeros((icroptype+1), float)
    # crop info for non-critical years
    NCBeginDate = zeros((icroptype+1), int)
    NCEndDate = zeros((icroptype+1), int)
    NCCropType = zeros((icroptype+1), int)
    NCf = zeros((icroptype+1), int)
    NCkc1 = zeros((icroptype+1), float)
    NCkc2 = zeros((icroptype+1), float)
    NCkc3 = zeros((icroptype+1), float)
    NCAB = zeros((icroptype+1), float)
    NCAC = zeros((icroptype+1), float)
    NCAD = zeros((icroptype+1), float)
    NCSDx = zeros((icroptype+1), float)
    NCRDxU = zeros((icroptype+1), float)
    NCRDxL = zeros((icroptype+1), float)
    NCawL = zeros((icroptype+1), float)
    NCawU = zeros((icroptype+1), float)
    NCADep = zeros((icroptype+1), float)

    BeginDateYear = zeros((iyears+2), int)

    osIkc = zeros((iyears+2), float)
    isIkc = zeros((iyears+2), float)
    IGETo = zeros((iyears+2), float)
    osFkc = zeros((iyears+2), float)

    KcByr = zeros((iyears+2), float)
    KcCyr = zeros((iyears+2), float)
    KcDyr = zeros((iyears+2), float)
    KcEyr = zeros((iyears+2), float)

    CCKc = zeros((idays+1), float)
    LIYear = zeros((iyears+2), int)
    BIYear = zeros((iyears+2), int)
    NA1 = zeros((iyears+2), float)
    NA2 = zeros((iyears+2), float)
    NA3 = zeros((iyears+2), float)
    isCETc = zeros((iyears+2), float)
    osCETc = zeros((iyears+2), float)
    isCERn = zeros((iyears+2), float)
    osCERn = zeros((iyears+2), float)
    isPCP = zeros((iyears+2), float)
    isETaw = zeros((iyears+2), float)

    isMCETc = zeros((imonths+1), float)
    osMCETc = zeros((imonths+1), float)
    isMCERn = zeros((imonths+1), float)
    osMCERn = zeros((imonths+1), float)
    isCSpg = zeros((iyears+2), float)
    osCSpg = zeros((iyears+2), float)
    isMCSpg = zeros((imonths+1), float)
    osMCSpg = zeros((imonths+1), float)
    EToMonthly = zeros((imonths+1), float)
    PcpMonthly = zeros((imonths+1), float)
    ETcMonthly = zeros((imonths+1), float)
    ERnMonthly = zeros((imonths+1), float)
    SpgMonthly = zeros((imonths+1), float)
    EspgMonthly = zeros((imonths+1), float)
    WSEspgMonthly = zeros((iyears+2, imonths+1), float)
    MonDsw = zeros((imonths+1), float)
    MonDswPos = zeros((imonths+1), float)
    MonNetApp = zeros((imonths+1), float)
    # yearType = numpy.empty(iyears+2, dtype='<U3')  # "  "*(iyears+1)

    # for HSA****.csv (not for OLDHSA****.csv)
    SACropDaily = "  "
    HSACropMonMean = "  "
    OldHSACropMonMean = "  "
    InpOldHSACropMonMean = "  "
    fpOldHSACrpMonMean = "  "
    # 4/20/09 yearTypeDaily = [] ##"  "*(idays+1)
    DateDaily = "  "*(idays+1)
    ##yearType = "  "
    MonDay = "  "
    flagETAW = "  "

    yDaily = zeros((idays+1), int)
    DOY = zeros((idays+1), int)
    EToDaily2 = zeros((idays+1), float)
    PcpDaily2 = zeros((idays+1), float)
    ETcDaily2 = zeros((idays+1), float)
    HAcreDaily = zeros((idays+1), float)
    OKcDaily = zeros((idays+1), float)
    IKcDaily = zeros((idays+1), float)
    CCKcDaily = zeros((idays+1), float)
    KcDaily = zeros((idays+1), float)
    ErDaily = zeros((idays+1), float)
    SpgDaily = zeros((idays+1), float)
    ESpgDaily = zeros((idays+1), float)
    DswDaily = zeros((idays+1), float)
    SWDDaily = zeros((idays+1), float)
    SWDxDaily = zeros((idays+1), float)
    FCDaily = zeros((idays+1), float)
    PWPDaily = zeros((idays+1), float)
    SWCDaily = zeros((idays+1), float)
    YTDDaily = zeros((idays+1), float)
    NADaily = zeros((idays+1), float)
    CPcpDaily = zeros((idays+1), float)
    CErDaily = zeros((idays+1), float)
    CESpgDaily = zeros((idays+1), float)
    CETcDaily = zeros((idays+1), float)
    CDswDaily = zeros((idays+1), float)
    ETAWDaily = zeros((idays+1), float)
    DETAWOUTPUT = zeros((6, ilands+1, icroptype+1, idates-1), float)

 # for irrigation and hydrology year convertion (((((
    ytemp = zeros((idays+1), int)
    DOYtemp = zeros((idays+1), int)
    ETotemp = zeros((idays+1), float)
    Pcptemp = zeros((idays+1), float)
    ETctemp = zeros((idays+1), float)
    HAcretemp = zeros((idays+1), float)
    OKctemp = zeros((idays+1), float)
    IKctemp = zeros((idays+1), float)
    CCKctemp = zeros((idays+1), float)
    Kctemp = zeros((idays+1), float)
    Ertemp = zeros((idays+1), float)
    Spgtemp = zeros((idays+1), float)
    ESpgtemp = zeros((idays+1), float)
    Dsw0temp = zeros((idays+1), float)
    SWDtemp = zeros((idays+1), float)
    SWDxtemp = zeros((idays+1), float)
    FC0temp = zeros((idays+1), float)
    PWPtemp = zeros((idays+1), float)
    SWC0temp = zeros((idays+1), float)
    YTDtemp = zeros((idays+1), float)
    NAtemp = zeros((idays+1), float)
    CPcptemp = zeros((idays+1), float)
    CErtemp = zeros((idays+1), float)
    CESpgtemp = zeros((idays+1), float)
    CETctemp = zeros((idays+1), float)
    CDswtemp = zeros((idays+1), float)
    HAcretemp = zeros((idays+1), float)
    # for irrigation and hydrology year convertion))))))
    EToMonthly = zeros((imonths+1), float)
    PcpMonthly = zeros((imonths+1), float)
    ETcMonthly = zeros((imonths+1), float)
    ERnMonthly = zeros((imonths+1), float)
    SpgMonthly = zeros((imonths+1), float)
    EspgMonthly = zeros((imonths+1), float)
    MonETAW = zeros((iyears+2, imonths+1), float)
    ETAWMonDay = zeros((imonths+1, 32), float)
    DOYLIrrig = zeros((iyears+2), float)
    DOYGrainLIrrig = zeros((iyears+2), float)
    DOYLIrrig[0] = 0
    DOYGrainLIrrig[0] = 0
    # for HSA****.csv  (not OLDHSA***.csv)
    # Add more variables for DETAW python program to efficient work
    ##cropareas = zeros((ilands,icroptype,iyears),float)
    # end:Add more variables for DETAW python program to efficient work
    PEs = 0.0
    PETo = 0.0
    DCT = 0
    ##j = 0
    CETo = 0.0
    METo = 0.0
    CEx = 0.0
    CEs = 0.0
    Es = 0.0
    EKc = 0.0
    Beta1 = 2.6
    LowIkc = 0.0
    LowFkc = 0.0

    IKc1 = 0.0
    OKc1 = 0.0
    Kc11 = 0.0
    ##IKcs = 0.0
    CC1 = 0.35
    ET1 = 0.0
    ET2 = 0.0
    CETc = 0.0
    NumI = 0
    NumI1 = 0
    NI2 = 0
    BI = 0
    LI = 0
    Dsw = 0.0
    Dswp = 0.0
    C = 0
    D = 0
    B = 0
    EndDate1 = 0
    BB = 0
    CC = 0
    DD = 0
    EE = 0
    HAcre_temp = 0.0
    Date1 = "  "
    yearTypeCal = "  "
    ##yearType[0] = "AN"
    NetApp = 0.0
    CCKc[:idays+1] = 0
    osIkc[:iyears+1] = 0
    isIkc[:iyears+1] = 0
    IGETo[:iyears+1] = 0
    osFkc[:iyears+1] = 0
    KcByr[:iyears+1] = 0
    KcCyr[:iyears+1] = 0
    KcDyr[:iyears+1] = 0
    KcEyr[:iyears+1] = 0
    dpyAll[:iyears+1] = 0
    # yrs[:iyears+1]=0
    NA1[:iyears+1] = 0
    NA2[:iyears+1] = 0
    NA3[:iyears+1] = 0
    EToDaily[:iyears+1, :idays+1] = 0
    # Date1Daily[:iyears+1,:idays+1]=0
    PcpDaily[:iyears+1, :idays+1] = 0
    OKc[:iyears+1, :idays+1] = 0
    IKc[:iyears+1, :idays+1] = 0
    Kc[:iyears+1, :idays+1] = 0

    # ++++++++++++++++++++++++++++++++++++++++++++++++++++
    ##
    # Determine ETo Rain base bare soil evaporation
    ##
    # ++++++++++++++++++++++++++++++++++++++++++++++++++++
    # step1: read percentage file: percentage, ET0_corrector, and region

    # step2: read \weatheroutput\ files, data:DayOfYear, TMax, TMin, PCP, ETo
    # The weather data come from the main program

    # step3: read crop information for critical years

    # if streamlinemodel == "CALSIM3":
    #     source = os.path.join(
    #         filepath, 'Input', 'planning_study', 'critical.csv')
    # else:
    #     source = os.path.join(
    #         filepath, 'Input', 'historical_study', 'critical.csv')
    ts_type = "rts"
    # crdf = read_and_clean_crop_info(source)
    CCropType[1:icroptype+1] = crdf.iloc[2, 1:].values
    CBeginDate[1:icroptype+1] = crdf.iloc[5, 1:].values
    CEndDate[1:icroptype+1] = crdf.iloc[6, 1:].values
    Cf[1:icroptype+1] = crdf.iloc[7, 1:].values
    Ckc1[1:icroptype+1] = crdf.iloc[8, 1:].values
    Ckc2[1:icroptype+1] = crdf.iloc[9, 1:].values
    Ckc3[1:icroptype+1] = crdf.iloc[10, 1:].values
    CAB[1:icroptype+1] = crdf.iloc[11, 1:].values
    CAC[1:icroptype+1] = crdf.iloc[12, 1:].values
    CAD[1:icroptype+1] = crdf.iloc[13, 1:].values
    CSDx[1:icroptype+1] = crdf.iloc[14, 1:].values
    CRDxL[1:icroptype+1] = crdf.iloc[15, 1:].values
    CRDxU[1:icroptype+1] = crdf.iloc[16, 1:].values
    CawL[1:icroptype+1] = crdf.iloc[17, 1:].values
    CawU[1:icroptype+1] = crdf.iloc[18, 1:].values
    CADep[1:icroptype+1] = crdf.iloc[19, 1:].values

    # step4: read crop information for non-critical years
    # if streamlinemodel == "CALSIM3":
    #     source = os.path.join(
    #         filepath, 'Input', 'planning_study', 'noncritical.csv')
    # else:
    #     source = os.path.join(
    #         filepath, 'Input', 'historical_study', 'noncritical.csv')
    # ncrdf = read_and_clean_crop_info(source)
    NCCropType[1:icroptype+1] = ncrdf.iloc[2, 1:].values
    NCBeginDate[1:icroptype+1] = ncrdf.iloc[5, 1:].values
    NCEndDate[1:icroptype+1] = ncrdf.iloc[6, 1:].values
    NCf[1:icroptype+1] = ncrdf.iloc[7, 1:].values
    NCkc1[1:icroptype+1] = ncrdf.iloc[8, 1:].values
    NCkc2[1:icroptype+1] = ncrdf.iloc[9, 1:].values
    NCkc3[1:icroptype+1] = ncrdf.iloc[10, 1:].values
    NCAB[1:icroptype+1] = ncrdf.iloc[11, 1:].values
    NCAC[1:icroptype+1] = ncrdf.iloc[12, 1:].values
    NCAD[1:icroptype+1] = ncrdf.iloc[13, 1:].values
    NCSDx[1:icroptype+1] = ncrdf.iloc[14, 1:].values
    NCRDxL[1:icroptype+1] = ncrdf.iloc[15, 1:].values
    NCRDxU[1:icroptype+1] = ncrdf.iloc[16, 1:].values
    NCawL[1:icroptype+1] = ncrdf.iloc[17, 1:].values
    NCawU[1:icroptype+1] = ncrdf.iloc[18, 1:].values
    NCADep[1:icroptype+1] = ncrdf.iloc[19, 1:].values

    # step5: read land use from .\Landuse folder !!!!!!!Not checked 3/13/2009
    # get year type of each year

    # if streamlinemodel == "CALSIM3":
    #     source = os.path.join(
    #         filepath, 'Input', 'planning_study', 'Landuse', 'SA0001.csv')
    # else:
    #     source = os.path.join(
    #         filepath, 'Input', 'historical_study', 'Landuse', 'SA0001.csv')
    # f0 = open(source)
    # iline0 = 1
    # ncount = 0
    # #? XXX why pad the first year with AN type?
    # yearType[ncount] = "AN"
    # for line in f0:
    #     if line:
    #         if line[0] == "1" or line[0] == "2":
    #             ncount = ncount+1
    #             yearType[ncount] = line.split(",")[1].strip()
    # ncount = ncount+1
    # #? XXX why pad the last year with AN type?
    # yearType[ncount] = "AN"  # for last year
    # # get Hectares of each crop type, year and island
    # if streamlinemodel == "CALSIM3":
    #     hist_path = os.path.join(
    #         filepath, 'Input', 'planning_study', 'Landuse')  # ---08/02/2010
    # else:
    #     hist_path = os.path.join(
    #         filepath, 'Input', 'historical_study', 'Landuse')
    # files = listdir(hist_path)
    # ts_type = "rts"
    # for file in files:
    #     if ".csv" in file:
    #         sheetname = file.replace(".csv", "")
    #         ilandno = int(sheetname.split("A0")[1])  # Landuse ---08/02/2010
    #         source = os.path.join(hist_path, file)
    #         df = pd.read_csv(source)
    #         if ilandno > ilands:
    #             break
    #         HAcre[ilandno-1, 1:iyears, 1:icroptype+1] = df.iloc[1:iyears,
    #                                                             2:icroptype+2].astype('float').values
    # End of data input
    # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ##
    # Calculating ETAW for each crop in SA
    ##
    # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    for k in range(0, ilands):
        print("island =", k+1)
        calc_soil_evap(k, Beta1, idates, idays, ts_year, ts_days,
                       start1, dpyAll, ET0, pcp, PcpDaily, EToDaily, OKc)
        #import pytest
        #PcpDaily == pytest.approx(_PcpDaily,rel=1e-2)
        #EToDaily == pytest.approx(_EToDaily,rel=1e-2)
        #OKc == pytest.approx(_OKc,rel=1e-2)
        WSCESpg[:, :] = 0
        WSEspgMonthly[:, :] = 0
        for j in range(1, icroptype+1):
            # for HSA*** (not for OLDHSA***)   +++++++++++++++++++++
            # initalize Monthly ETaw
            #print("island =",k+1, " croptype =",j)
            MonETAW[:, :] = 0
            # initialize cumulative ETAW
            CETAWDaily = 0
            # initializa DoyLast irrigation
            DOYLIrrig[:] = 0
            DOYGrainLIrrig[:] = 0
            # for HSA*** (not for OLDHSA***)  +++++++++++++++++++++
            datadays = zeros((32, idates-1), float)
            dataday = zeros(idates-1, float)
            data2day = zeros(idates-1, float)
            data3day = zeros(idates-1, float)
            data4day = zeros(idates-1, float)
            data5day = zeros(idates-1, float)
            data6day = zeros(idates-1, float)
            data7day = zeros(idates-1, float)
            data8day = zeros(idates-1, float)
            data9day = zeros(idates-1, float)
            data10day = zeros(idates-1, float)
            data11day = zeros(idates-1, float)
            data12day = zeros(idates-1, float)
            data13day = zeros(idates-1, float)
            data14day = zeros(idates-1, float)
            data15day = zeros(idates-1, float)
            data16day = zeros(idates-1, float)
            data17day = zeros(idates-1, float)
            data18day = zeros(idates-1, float)
            data19day = zeros(idates-1, float)
            data20day = zeros(idates-1, float)
            data21day = zeros(idates-1, float)
            data22day = zeros(idates-1, float)
            data23day = zeros(idates-1, float)
            data24day = zeros(idates-1, float)
            data25day = zeros(idates-1, float)
            data26day = zeros(idates-1, float)
            data27day = zeros(idates-1, float)
            data28day = zeros(idates-1, float)
            data29day = zeros(idates-1, float)
            data30day = zeros(idates-1, float)
            data31day = zeros(idates-1, float)
            data32day = zeros(idates-1, float)
            data1day_14crops = zeros(idates-1, float)
            data10day_14crops = zeros(idates-1, float)
            data11day_14crops = zeros(idates-1, float)
            data12day_14crops = zeros(idates-1, float)
            data13day_14crops = zeros(idates-1, float)
            data1day_water = zeros(idates-1, float)
            data10day_water = zeros(idates-1, float)
            data11day_water = zeros(idates-1, float)
            data12day_water = zeros(idates-1, float)
            data13day_water = zeros(idates-1, float)

            nmonths = 12*(iyears-1)
            datamon = zeros(nmonths, float)
            data2mon = zeros(nmonths, float)
            data3mon = zeros(nmonths, float)
            data4mon = zeros(nmonths, float)
            data5mon = zeros(nmonths, float)
            data6mon = zeros(nmonths, float)
            data7mon = zeros(nmonths, float)
            data8mon = zeros(nmonths, float)
            data9mon = zeros(nmonths, float)
            data10mon = zeros(nmonths, float)
            data11mon = zeros(nmonths, float)
            data12mon = zeros(nmonths, float)
            data13mon = zeros(nmonths, float)
            data14mon = zeros(nmonths, float)
            data15mon = zeros(nmonths, float)
            data16mon = zeros(nmonths, float)
            data17mon = zeros(nmonths, float)
            data18mon = zeros(nmonths, float)
            data19mon = zeros(nmonths, float)
            data20mon = zeros(nmonths, float)
            data21mon = zeros(nmonths, float)

            datayr = zeros(iyears-1, float)
            data2yr = zeros(iyears-1, float)
            data3yr = zeros(iyears-1, float)
            data4yr = zeros(iyears-1, float)
            data5yr = zeros(iyears-1, float)
            data6yr = zeros(iyears-1, float)
            data7yr = zeros(iyears-1, float)
            data8yr = zeros(iyears-1, float)
            kkc1, kkc2, kkc3, EndDate1, BeginDate1, AB1, AC1, AD1, CropType1, y, f1, YTD, osSWDx, erd, aw1, ADep1 = calc_kc_vals(iyears, yearType, CCropType, j, CBeginDate, CEndDate, Cf, Ckc1, Ckc2, Ckc3, CAB, CAC, CAD, CSDx, CRDxU, CRDxL, CawL, CawU, CADep, Region, k, NCCropType, NCBeginDate, NCEndDate, NCf, NCkc1, NCkc2, NCkc3, NCAB, NCAC, NCAD, NCSDx, NCRDxU, NCRDxL, NCawL, NCawU, NCADep, osIkc, OKc, EToDaily, IGETo, osFkc, isIkc, Beta1)
            # end of for  y<iyears  ********************************
            # calculates values of osIkc, isIkc, osFkc and IGETo
            # Identify initial growth Kc(KcB) and final Kc(KcE)
            LowIkc = min(2, numpy.amin(osIkc[1:iyears+1]))
            LowFkc = min(2, numpy.amin(osFkc[1:iyears+1]))
            # ......... Print Kc value selection process
            # ......   Note that Kc's are not printed in final version
            # ........  Print rows were changed to remarks. however, the
            # loop is retained to assign Kc's to subscripts.
            # ...........................................................

            # Loop through years to calculate daily Kc values
            kkc1, kkc2, kkc3, EndDate1, BeginDate1, AB1, AC1, AD1, y = calc_kc_daily(iyears, j, yearType, CBeginDate, CEndDate, Ckc1, Ckc2, Ckc3, CAB, CAC, CAD, NCBeginDate,
                                                                                     NCEndDate, NCkc1, NCkc2, NCkc3, NCAB, NCAC, NCAD, kkc1, kkc2, kkc3, EndDate1, BeginDate1, AB1, AC1, AD1, CropType1, LowIkc, isIkc, KcByr, y, KcCyr, KcDyr, KcEyr)
            # end of for y<YCt

            # *****************************************************
            ##
            # Write Daily results in OldHSA****C***.csv (OutSACropDaily)
            # Write yearly results(Daily mean) in HSA****C***.eaw.csv (OutHSACropDailyMean)
            # Write Monthly results in OldHSA****C***.mtv.csv (OutHSACropMonMean)
            ##
            # *****************************************************

            # This loop determines daily IKc for each year and subscripts
            # the results by year and day(ends after 5020)
            #import pdb; pdb.set_trace()
            calc_ikc_daily(iyears, j, yearType, CBeginDate, CEndDate, Ckc1, Ckc2, Ckc3, CAB, CAC, CAD, NCBeginDate, NCEndDate, NCkc1, NCkc2, NCkc3, NCAB, NCAC, NCAD, kkc1, kkc2, kkc3,
                           EndDate1, BeginDate1, AB1, AC1, AD1, IKc, y, idays, dpyAll, BeginDateYear, IKc1, OKc, CropType1, CCKc, CC1, Kc, EToDaily, ETcDaily, f1, BIYear, NA1, YTD, LIYear, NA2, NA3)
            ##ii = EndDate1 + 1
            # to get out of loop i<=E
            # end of for i<=E
            # end of y<YCI

            #? eliminate loop ?#
            isCERn[1:iyears+1] = 0
            isPCP[1:iyears+1] = 0
            isETaw[1:iyears+1] = 0
            osCERn[1:iyears+1] = 0
            isCSpg[1:iyears+1] = 0
            osCSpg[1:iyears+1] = 0
            isCETc[1:iyears+1] = 0
            osCETc[1:iyears+1] = 0
            #? eliminate loop ?#
            Mon = slice(0, 12)
            isMCERn[Mon] = 0
            osMCERn[Mon] = 0
            isMCSpg[Mon] = 0
            osMCSpg[Mon] = 0
            isMCETc[Mon] = 0
            osMCETc[Mon] = 0
            EToMonthly[Mon] = 0
            ETcMonthly[Mon] = 0
            PcpMonthly[Mon] = 0
            ERnMonthly[Mon] = 0
            SpgMonthly[Mon] = 0
            EspgMonthly[Mon] = 0
            MonDsw[Mon] = 0
            MonDswPos[Mon] = 0
            MonNetApp[Mon] = 0

            # Loop to calculate Kc's,Etc, &SWD for sceduling
            PSW = 0.0
            SWD = 0.0
            SWD0 = osSWDx
            Espg = 0

            # for HSA**** (not for OLDHSA****)
            IrrigYear = 0
            #? The most costly loop below ~ 2.8s?#
            date_index = 0
            main_calc_loop(iyears, j, yearType, CBeginDate, CEndDate, Ckc1, Ckc2, Ckc3, CAB, CAC, CAD, NCBeginDate, NCEndDate, NCkc1, NCkc2, NCkc3, NCAB, NCAC, NCAD,
                           kkc1, kkc2, kkc3, EndDate1, BeginDate1, AB1, AC1, AD1, EToMonthly, ETcMonthly, PcpMonthly, ERnMonthly, SpgMonthly, EspgMonthly, MonNetApp, MonDsw, MonDswPos,
                           dpyAll, NumDaysPerMon, erd, Region, k, SWD, OKc, IKc, Kc, EToDaily, PcpDaily, ETcDaily, NA1, BIYear, NA2, LIYear, NA3, osSWDx, SWD0, Dsw, NetApp, NII, NI, osCETc,
                           HAcre, osMCETc, isCETc, isMCETc, osCERn, osCSpg, osMCERn, osMCSpg, isCERn, isPCP, isETaw, isCSpg, isMCERn, isMCSpg, aw1, ADep1, Espg, WSCESpg, WSEspgMonthly,
                           ytemp, DOYtemp, IrrigYear, HAcre_temp, HAcretemp, OKctemp, IKctemp, CCKc, CCKctemp, ETotemp, Kctemp, ETctemp, Pcptemp, Ertemp, Spgtemp, ESpgtemp, Dsw0temp,
                           SWDtemp, SWDxtemp, FC0temp, PWPtemp, SWC0temp, YTDtemp, NAtemp, DOYLIrrig, DOYGrainLIrrig, CPcptemp, CErtemp, CESpgtemp, CETctemp, CDswtemp, HAcreDaily,
                           yDaily, DOY, OKcDaily, IKcDaily, CCKcDaily, EToDaily2, KcDaily, ETcDaily2, PcpDaily2, ErDaily, SpgDaily, ESpgDaily, DswDaily, SWDDaily, SWDxDaily, FCDaily,
                           PWPDaily, SWCDaily, YTDDaily, NADaily, CPcpDaily, CErDaily, CESpgDaily, CETcDaily, CDswDaily, ETAWMonDay, ETAWDaily, MonETAW, idayoutput, dailyunit,
                           dataday, date_index, data2day, data6day, data7day, data8day, data9day, data10day, data11day, data12day, data13day, data14day, data15day, data16day,
                           data17day, data18day, data26day, data27day, data28day, data29day, data30day, data31day, CETAWDaily, data32day, imonthoutput,
                           datamon, data2mon, data3mon, data4mon, data5mon, data6mon, data7mon, data8mon, data9mon, data12mon, data13mon, data14mon, data15mon,
                           data16mon, data17mon, data18mon, data21mon, data10mon, data19mon, data11mon, data20mon,model_start_year)
            SACETC = 0
            SACERN = 0
            SACSpg = 0
            SISCETC = 0
            SISCERN = 0
            SISCSpg = 0
            SOSCETC = 0
            SOSCERN = 0
            SOSCSpg = 0
            for y in range(1, iyears+1):
                SISCETC = SISCETC+isCETc[y]
                SISCERN = SISCERN+isCERn[y]
                SISCSpg = SISCSpg+isCSpg[y]
                SOSCETC = SOSCETC+osCETc[y]
                SOSCERN = SOSCERN+osCERn[y]
                SOSCSpg = SOSCSpg+osCSpg[y]

                # **save: y+1921,yearType[y],HAcre[k,y,j]*2.471,isPCP[y],isCETc[y],isCERn[y],isCSpg[y],
                # **isETaw[y],osCETc[y],osCERn[y],osCSpg[y]
                if iyearoutput == 1:
                    datayr[y-1] = (isCETc[y])
                    data2yr[y-1] = (isPCP[y])
                    data3yr[y-1] = (isCERn[y])
                    data4yr[y-1] = (isCSpg[y])
                    data5yr[y-1] = (isETaw[y])
                    data6yr[y-1] = (osCETc[y])
                    data7yr[y-1] = (osCERn[y])
                    data8yr[y-1] = (osCSpg[y])
            # calculate &print mean over years for CETc,CERn, ETAW
            MACETC = SACETC/iyears
            MACERN = SACERN/iyears
            MACSpg = SACSpg/iyears
            MISCETC = SISCETC/iyears
            MISCERN = SISCERN/iyears
            MISCSpg = SISCSpg/iyears
            MOSCETC = SOSCETC/iyears
            MOSCERN = SOSCERN/iyears
            MOSCSpg = SOSCSpg/iyears

            ###-------------OUTPUT CODE BELOW HERE----------------###
            #start = start1+days(2)
            start = str(start1[2]+1)+monthname[start1[1]-1] + \
                str(start1[0])+" "+"2300"
            dt = "1DAY"  # days(1)
            unit = ""
            if idayoutput == 1:
                DETAWOUTPUT[0, k, j-1, :] = dataday[:]
                DETAWOUTPUT[1, k, j-1, :] = data10day[:]
                DETAWOUTPUT[2, k, j-1, :] = data2day[:]
                DETAWOUTPUT[3, k, j-1, :] = data12day[:]
                DETAWOUTPUT[4, k, j-1, :] = data11day[:]
                DETAWOUTPUT[5, k, j-1, :] = data13day[:]
                if forDSM2_daily == 1:
                    ddatalist = dataday, data10day, data2day, data12day, data11day, data13day
                    dcpartlist = "ETc", "ESpg", "PCP", "ETAW", "Dsw", "ER"
                    dunitlist = "A-ft", "A-ft", "A-ft", "A-ft", "A-ft", "A-ft"
                else:
                    ddatalist = dataday, data2day, data6day, data7day, data8day, data9day,  \
                        data10day, data11day, data12day, data13day,    \
                        data14day, data15day, data16day, data17day, data18day,  \
                        data26day, data27day, data28day, data29day, data30day, data31day, data32day
                    dcpartlist = "ETc", "PCP", "ETo", "Kc", "SWC", "Spg",   \
                        "ESpg", "Dsw", "ETAW", "ER",    \
                        "SWDxDaily", "FCDaily", "PWPDaily", "SWDDaily", "YTDDaily", \
                        "NA", "CPcp", "CEr", "CEspg", "CETc", "CDsw", "CETaw"
                    if dailyunit == 1:
                        dunitlist = "A-ft", "A-ft", "A-ft", "A-ft", "A-ft", "A-ft",   \
                            "A-ft", "A-ft", "A-ft", "A-ft",    \
                            "A-ft", "A-ft", "A-ft", "A-ft", "A-ft",    \
                            "A-ft", "A-ft", "A-ft", "A-ft", "A-ft", "A-ft", "A-ft"
                    else:
                        dunitlist = "mm", "mm", "mm", "mm", "mm", "mm",   \
                            "mm", "mm", "mm", "mm",    \
                            "mm", "mm", "mm", "mm", "mm",    \
                            "A-ft", "A-ft", "A-ft", "A-ft", "A-ft", "A-ft", "A-ft"
                if itotaloutput == 1:
                    if j == 1 and k == 0:
                        data1day_14crops = dataday
                        data10day_14crops = data10day
                        data11day_14crops = data11day
                        data12day_14crops = data12day
                        data13day_14crops = data13day
                    elif j < 15:
                        data1day_14crops = data1day_14crops + dataday
                        data10day_14crops = data10day_14crops + data10day
                        data11day_14crops = data11day_14crops + data11day
                        data12day_14crops = data12day_14crops + data12day
                        data13day_14crops = data13day_14crops + data13day
                    if j == 15:
                        if k == 0:
                            data1day_water = dataday
                            data10day_water = data10day
                            data11day_water = data11day
                            data12day_water = data12day
                            data13day_water = data13day
                        else:
                            data1day_water = data1day_water + dataday
                            data10day_water = data10day_water + data10day
                            data11day_water = data11day_water + data11day
                            data12day_water = data12day_water + data12day
                            data13day_water = data13day_water + data13day
            ##start = start1
            dt2 = "1MONTH"  # months(1)
            if imonthoutput == 1:
                #props = {"unit":"mm"}
                mdatalist = datamon, data2mon, data3mon, data4mon, data5mon, data6mon, \
                    data7mon, data8mon, data9mon, data10mon, data11mon, \
                    data12mon, data13mon, data14mon, data15mon, data16mon,  \
                    data17mon, data18mon, data19mon, data20mon, data21mon
                mcpartlist = "ETc", "ETo", "NetApp", "Pcp", "ERn", "Spg",  \
                             "Espg", "Dsw", "Dsw+ve", "ETAW", "ETAW+ve",  \
                             "NA_A_ft", "Pcp_A_ft", "Er_A_ft", "Espg_A_ft", "ETc_A_ft",  \
                             "Dsw_A_ft", "Dsw+ve_A_ft", "ETaw_A_ft", "ETaw+ve_A_ft", "ETo_A_ft"
                munitlist = "mm", "mm", "mm", "mm", "mm", "mm",  \
                    "mm", "mm", "mm", "mm", "mm",  \
                    "A-ft", "A-ft", "A-ft", "A-ft", "A-ft",  \
                    "A-ft", "A-ft", "A-ft", "A-ft", "A-ft"

                if itotaloutput != 1:
                    if j == 2:
                        data12mon_veg = data12mon
                        data13mon_veg = data13mon
                        data14mon_veg = data14mon
                        data15mon_veg = data15mon
                        data16mon_veg = data16mon
                        data17mon_veg = data17mon
                        data18mon_veg = data18mon
                        data19mon_veg = data19mon
                        data20mon_veg = data20mon
                        data21mon_veg = data21mon
                    if j > 2 and j < 12:
                        data12mon_veg = add(data12mon_veg, data12mon)
                        data13mon_veg = add(data13mon_veg, data13mon)
                        data14mon_veg = add(data14mon_veg, data14mon)
                        data15mon_veg = add(data15mon_veg, data15mon)
                        data16mon_veg = add(data16mon_veg, data16mon)
                        data17mon_veg = add(data17mon_veg, data17mon)
                        data18mon_veg = add(data18mon_veg, data18mon)
                        data19mon_veg = add(data19mon_veg, data19mon)
                        data20mon_veg = add(data20mon_veg, data20mon)
                        data21mon_veg = add(data21mon_veg, data21mon)
                if itotaloutput == 1:
                    if j == 1 and k == 0:
                        data14mon_total = data14mon
                        data15mon_total = data15mon
                        data16mon_total = data16mon
                        data20mon_total = data20mon
                    elif j < 15:
                        data14mon_total = add(data14mon_total, data14mon)
                        data15mon_total = add(data15mon_total, data15mon)
                        data16mon_total = add(data16mon_total, data16mon)
                        data20mon_total = add(data20mon_total, data20mon)
                    if j == 15:
                        if k == 0:
                            data14mon_water = data14mon
                            data15mon_water = data15mon
                            data16mon_water = data16mon
                            data20mon_water = data20mon
                        else:
                            data14mon_water = add(data14mon_water, data14mon)
                            data15mon_water = add(data15mon_water, data15mon)
                            data16mon_water = add(data16mon_water, data16mon)
                            data20mon_water = add(data20mon_water, data20mon)
            ##start = start1
            dt3 = "1YEAR"  # years(1)
            if iyearoutput == 1:
                ydatalist = datayr, data2yr, data3yr, data4yr, data5yr, data6yr, data7yr, data8yr
                ycpartlist = "CETc", "CPcp", "CEr", "CSpg", "CETaw", "OCETc", "OCEr", "OCSpg"

            Apart = "DWRBDO-DETAW"
            if (k+1) < 10:
                Bpart = "HSA_000" + str((k+1))
                Bveg = "HSA_000" + str((k+1)) + "_veg"
            elif (k+1) < 100:
                Bpart = "HSA_00" + str((k+1))
                Bveg = "HSA_00" + str((k+1)) + "_veg"
            else:
                Bpart = "HSA_0" + str((k+1))
                Bveg = "HSA_0" + str((k+1)) + "_veg"

            Cpart = cpartt
            Epart = "1DAY"
            Fpart = "Crop_"+str(j)

            ktemp = int(k/10)

            if idayoutput == 1 and not NO_OUTPUT:
                destination = os.path.join(
                    filepath, 'Output', 'DETAW_day_'+str(ktemp)+'.dss')
                dssfh = pyhecdss.DSSFile(destination, create_new=True)
                for ilist in range(0, len(ddatalist)):
                    path = "/"+Apart+"/"+Bpart+"/" + \
                        dcpartlist[ilist]+"//"+Epart+"/"+Fpart+"/"
                    listlen = len(ddatalist[ilist])
                    write_to_dss(
                        dssfh, ddatalist[ilist], path, startdate + " "+starttime, dunitlist[ilist], ctype)
                if itotaloutput == 1:
                    if j == 15 and k == (ilands-1):
                        print("daily output")
                        path = "/"+Apart+"/"+Bpart+"/ETc//"+Epart+"/"+Fpart+"/"
                        listlen = len(data1day_14crops)
                        if dailyunit == 1:
                            dunits = "A-FT"
                        else:
                            dunits = "mm"
                        write_to_dss(dssfh, data1day_14crops, path,
                                     startdate + " "+starttime, dunits, ctype)
                        path = "/"+Apart+"/"+Bpart+"/ESpg//"+Epart+"/"+Fpart+"/"
                        write_to_dss(dssfh, data10day_14crops, path,
                                     startdate + " "+starttime, dunits, ctype)
                        path = "/"+Apart+"/"+Bpart+"/Dsw//"+Epart+"/"+Fpart+"/"
                        write_to_dss(dssfh, data11day_14crops, path,
                                     startdate + " "+starttime, dunits, ctype)
                        path = "/"+Apart+"/"+Bpart+"/ETaw//"+Epart+"/"+Fpart+"/"
                        write_to_dss(dssfh, data12day_14crops, path,
                                     startdate + " "+starttime, dunits, ctype)
                        path = "/"+Apart+"/"+Bpart+"/Er//"+Epart+"/"+Fpart+"/"
                        write_to_dss(dssfh, data13day_14crops, path,
                                     startdate + " "+starttime, dunits, ctype)
                        path = "/"+Apart+"/"+Bpart+"/ETc//"+Epart+"/"+Fpart+"/"
                        write_to_dss(dssfh, data1day_water, path,
                                     startdate + " "+starttime, dunits, ctype)
                        path = "/"+Apart+"/"+Bpart+"/ESpg//"+Epart+"/"+Fpart+"/"
                        write_to_dss(dssfh, data1day_water, path,
                                     startdate + " "+starttime, dunits, ctype)
                        path = "/"+Apart+"/"+Bpart+"/Dsw//"+Epart+"/"+Fpart+"/"
                        write_to_dss(dssfh, data11day_water, path,
                                     startdate + " "+starttime, dunits, ctype)
                        path = "/"+Apart+"/"+Bpart+"/ETaw//"+Epart+"/"+Fpart+"/"
                        write_to_dss(dssfh, data13day_14crops, path,
                                     startdate + " "+starttime, dunits, ctype)
                        path = "/"+Apart+"/"+Bpart+"/Er//"+Epart+"/"+Fpart+"/"
                        write_to_dss(dssfh, data13day_water, path,
                                     startdate + " "+starttime, dunits, ctype)
                dssfh.close()

            Epart = "1MONTH"
            if imonthoutput == 1 and not NO_OUTPUT:
                destination = os.path.join(
                    filepath, 'Output', 'DETAW_month.dss')
                dssfh = pyhecdss.DSSFile(destination, create_new=True)

                for ilist in range(0, 21):
                    path = "/"+Apart+"/"+Bpart+"/" + \
                        mcpartlist[ilist]+"//"+Epart+"/"+Fpart+"/"
                    listlen = len(mdatalist[ilist])
                    templist = list(mdatalist[ilist])
                    tempunit = munitlist[ilist]
                    write_to_dss(dssfh, templist, path, startdate +
                                 " "+starttime, tempunit, ctype)

                if itotaloutput == 1:
                    if j == 15 and k == (ilands-1):
                        print("monthly output")
                        listlen = len(data14mon_total)
                        path = "/"+Apart+"/total_no_w/Er//"+Epart+"/"+Fpart+"/"
                        write_to_dss(dssfh, data14mon_total, path,
                                     startdate + " "+starttime, "A-ft", ctype)
                        path = "/"+Apart+"/total_no_w/Espg//"+Epart+"/"+Fpart+"/"
                        write_to_dss(dssfh, data15mon_total, path,
                                     startdate + " "+starttime, "A-ft", ctype)
                        path = "/"+Apart+"/total_no_w/ETc//"+Epart+"/"+Fpart+"/"
                        write_to_dss(dssfh, data16mon_total, path,
                                     startdate + " "+starttime, "A-ft", ctype)
                        path = "/"+Apart+"/total_no_w/ETaw+ve//"+Epart+"/"+Fpart+"/"
                        write_to_dss(dssfh, data20mon_total, path,
                                     startdate + " "+starttime, "A-ft", ctype)
                        path = "/"+Apart+"/water/Er//"+Epart+"/"+Fpart+"/"
                        write_to_dss(dssfh, data14mon_water, path,
                                     startdate + " "+starttime, "A-ft", ctype)
                        path = "/"+Apart+"/water/Espg//"+Epart+"/"+Fpart+"/"
                        write_to_dss(dssfh, data15mon_water, path,
                                     startdate + " "+starttime, "A-ft", ctype)
                        path = "/"+Apart+"/water/ETc//"+Epart+"/"+Fpart+"/"
                        write_to_dss(dssfh, data16mon_water, path,
                                     startdate + " "+starttime, "A-ft", ctype)
                        path = "/"+Apart+"/water/ETaw+ve//"+Epart+"/"+Fpart+"/"
                        write_to_dss(dssfh, data20mon_water, path,
                                     startdate + " "+starttime, "A-ft", ctype)
                else:
                    if j == 11:
                        listlen = len(data12mon_veg)
                        path = "/"+Apart+"/veg/NA//"+Epart+"/"+Fpart+"/"
                        write_to_dss(dssfh, data12mon_veg, path,
                                     startdate + " "+starttime, "A-ft", ctype)
                        path = "/"+Apart+"/veg/Pcp//"+Epart+"/"+Fpart+"/"
                        write_to_dss(dssfh, data13mon_veg, path,
                                     startdate + " "+starttime, "A-ft", ctype)
                        path = "/"+Apart+"/veg/Er//"+Epart+"/"+Fpart+"/"
                        write_to_dss(dssfh, data14mon_veg, path,
                                     startdate + " "+starttime, "A-ft", ctype)
                        path = "/"+Apart+"/veg/Espg//"+Epart+"/"+Fpart+"/"
                        write_to_dss(dssfh, data15mon_veg, path,
                                     startdate + " "+starttime, "A-ft", ctype)
                        path = "/"+Apart+"/veg/ETc//"+Epart+"/"+Fpart+"/"
                        write_to_dss(dssfh, data16mon_veg, path,
                                     startdate + " "+starttime, "A-ft", ctype)
                        path = "/"+Apart+"/veg/Dsw//"+Epart+"/"+Fpart+"/"
                        write_to_dss(dssfh, data17mon_veg, path,
                                     startdate + " "+starttime, "A-ft", ctype)
                        path = "/"+Apart+"/veg/Dsw+ve//"+Epart+"/"+Fpart+"/"
                        write_to_dss(dssfh, data18mon_veg, path,
                                     startdate + " "+starttime, "A-ft", ctype)
                        path = "/"+Apart+"/veg/ETaw//"+Epart+"/"+Fpart+"/"
                        write_to_dss(dssfh, data19mon_veg, path,
                                     startdate + " "+starttime, "A-ft", ctype)
                        path = "/"+Apart+"/veg/ETaw+ve//"+Epart+"/"+Fpart+"/"
                        write_to_dss(dssfh, data20mon_veg, path,
                                     startdate + " "+starttime, "A-ft", ctype)
                        path = "/"+Apart+"/veg/ETo//"+Epart+"/"+Fpart+"/"
                        write_to_dss(dssfh, data21mon_veg, path,
                                     startdate + " "+starttime, "A-ft", ctype)
                dssfh.close()

            Epart = "1YEAR"
            if iyearoutput == 1:
                destination = os.path.join(
                    filepath, 'Output', 'DETAW_year.dss')
                dssfh = pyhecdss.DSSFile(destination, create_new=True)
                for ilist in range(0, len(ydatalist)):
                    path = "/"+Apart+"/"+Bpart+"/" + \
                        ycpartlist[ilist]+"//"+Epart+"/"+Fpart+"/"
                    listlen = len(ydatalist[ilist])
                    write_to_dss(
                        dssfh, ydatalist[ilist], path, startdate + " "+starttime, yunitlist[ilist], ctype)
                dssfh.close()
    return(DETAWOUTPUT)


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


def read_pcp(start_date_str, end_date_str, fn):
    """ Read pcp input data

        Parameters
        ----------
            streamlinemodel: str
            start_date_str: str
            end_date_str: str
            pcplocs: str
            fn: str
        Returns
        -------
            ts_pcp: timeseries of station precipitation
    """
    # FIXME Avoid to use a current directory for jobs
    # filepath = os.getcwd()

    # if streamlinemodel == "CALSIM3":
    #     source = os.path.join(filepath, 'Input', 'planning_study', fn)
    # else:
    #     source = os.path.join(filepath, 'Input', 'historical_study', fn)

    pcp_df = pd.read_csv(fn,index_col=False)
    # add date column to dataframe
    pcp_df['file_dates'] = pd.to_datetime(pcp_df[['year', 'month', 'day']])
    mask = pcp_df['file_dates'].between(start_date_str,end_date_str)
    # get the column_names from the dataframe
    column_names = pcp_df.columns.values.tolist()
    skip_cols =['year','month','day','DOY','file_dates']
    pcplocs = list((filter(lambda val: val not in skip_cols, column_names)))
    # subset the precip values and transpose before converting to array
    ts_pcp = pcp_df.loc[mask,pcplocs].T.to_numpy()
    return(ts_pcp)


def read_temperature(streamlinemodel, start_date_str,end_date_str,fn):
    """ Read temperature input data

        Parameters
        ----------
            streamlinemodel: str
            start_date_str: str
            end_date_str: str
            fn: str
        Returns
        -------
            ts_year: array of years
            ts_mon: array of mon
            ts_days: array of days
            ts_LODI_tx: array tmax
            ts_LODI_tn: array tmin
    """
    # FIXME Avoid to use a current directory for jobs
    # filepath = os.getcwd()

    # if streamlinemodel == "CALSIM3":
    #     source = os.path.join(filepath, 'Input', 'planning_study', fn)
    # else:
    #     source = os.path.join(filepath, 'Input', 'historical_study', fn)

    temp_df = pd.read_csv(fn,parse_dates=[0],index_col=0,header=0)

    ts_year = temp_df[start_date_str:end_date_str]['Year'].T.to_numpy()
    ts_mon = temp_df[start_date_str:end_date_str]['Month'].T.to_numpy()
    ts_days = temp_df[start_date_str:end_date_str]['DOY'].T.to_numpy()
    ts_LODI_tx = temp_df[start_date_str:end_date_str]['Tx(oC)'].T.to_numpy()
    ts_LODI_tn = temp_df[start_date_str:end_date_str]['Tn(oC)'].T.to_numpy()

    return(ts_year,ts_mon,ts_days,ts_LODI_tx,ts_LODI_tn)


def read_landuse(fn_landuse, iyears, water_years, n_areas):
    """ Read landuse acreages

        Parameters
        ----------
            fn_landuse: str
            iyears: int
            water_years: array
            n_areas: int
                number of areas in the model

        Returns
        -------
            yearType: array
            Landuse Area: array
            icroptype: int
    """
    # # FIXME Avoid to use a current directory for jobs
    # filepath = os.getcwd()

    lu_comb_df = pd.read_csv(fn_landuse,  header=[0])
    mask = lu_comb_df['DATE'].isin(water_years)
    # clipped based on water years
    sub_df = lu_comb_df.loc[mask, :]

    column_names = list(lu_comb_df.columns)
    skip_cols =['DATE','TYPE','area_id']
    lucols = list((filter(lambda val: val not in skip_cols, column_names)))
    icroptype = len(lucols)

    # FIXME year_type is set to iyears but should be water_years like
    year_type = numpy.empty(iyears+1, dtype='<U3')
    # FIXME Why do we need n_areas from the argument? Can it be inferred from the landuse?
    landuse_area_hectare = zeros((n_areas, iyears + 2, icroptype + 1), float)

    lu_hectare_lf = sub_df.loc[:,lucols].values
    landuse_area_hectare[:, 1:iyears, 1:] = lu_hectare_lf.reshape((n_areas,len(water_years),icroptype))

    # FIXME the padding  with AN on either end is carry forward from legacy code
    year_type[0]="AN"
    year_type[1:len(water_years)+1] = sub_df.iloc[0:len(water_years),1].values
    year_type[-1] = "AN"
    year_type = year_type.astype(dtype='<U3')

    return (year_type, landuse_area_hectare, icroptype)


def read_et_correction_factors(fn):
    """ Read et correction factors

        Parameters
        ----------
            streamlinemodel: str
            fn: str
        Returns
        -------
            ts_per: array
            ETo_corrector: array
            Region: array
    """
    # FIXME Avoid to use a current directory for jobs
    # filepath = os.getcwd()

    # if streamlinemodel == "CALSIM3":
    #     source = os.path.join(filepath, 'Input', 'planning_study', fn)
    # else:
    #     source = os.path.join(filepath, 'Input', 'historical_study', fn)

    # read data from the csv file
    data = pd.read_csv(fn, index_col=False)
    skip_cols = ['area_id', 'SubArea', 'extension', 'ETo Correction Factor', 'REGION', 'REGION.1']
    ts_per = data[data.columns[~data.columns.isin(skip_cols)]].T.to_numpy()
    ETo_corrector = data.loc[:, 'ETo Correction Factor'].T.to_numpy()
    Region = data.loc[:, 'REGION'].T.to_numpy()
    return (ts_per, ETo_corrector, Region)


def convert_to_absolute_path(filename, dir_input_base):
    """Check if a file name is in an absolute path and convert it into one
    by adding the base directory of the main input

    Parameters
    ----------
    filename: str
        the file name to be checked
    dir_input_base: str
        The base directory name to add

    Returns
    -------
    pathlib.Path
        An absolute path
    """
    path = Path(filename)
    if not path.is_absolute():
        path = dir_input_base / Path(path)
    return str(path)


def detaw(fname_main_yaml: str) -> None:
    """ Run DETAW and DCD with a YAML file

        Parameters
        ----------
        fname_main_yaml: str
            the main input file name
    """
    logging.info(f"Reading the main YAML input: {fname_main_yaml}")
    with open(fname_main_yaml, 'r') as file_in:
        model_params = yaml.safe_load(file_in)
    dir_input_base = Path(fname_main_yaml).resolve().parent

    # FIXME For now, passing the yaml information to the current model
    #       parameters. They can be passed as a dict directly.
    detaw_params = model_params["detaw"]
    streamlinemodel = detaw_params["target_model"]
    idayoutput = detaw_params["daily_output"]
    imonthoutput = detaw_params["monthly_output"]
    iyearoutput = detaw_params["yearly_output"]
    itotaloutput = detaw_params["delta_output"]
    dailyunit = detaw_params["daily_output_unit"]
    forDSM2_daily = detaw_params["for_dsm2_only"]
    start_water_year = detaw_params["start_water_year"]
    end_water_year = detaw_params['end_water_year']
    fn_input_pcp = convert_to_absolute_path(detaw_params['input_pcp'], dir_input_base)
    fn_input_temperature = convert_to_absolute_path(detaw_params['input_temperature'], dir_input_base)
    fn_landuse = convert_to_absolute_path(detaw_params['landuse'], dir_input_base)
    fn_et_correction = convert_to_absolute_path(detaw_params['et_correction'], dir_input_base)
    fn_detaw_output = convert_to_absolute_path(detaw_params['detaw_output'], dir_input_base)
    fn_precip_output = convert_to_absolute_path(detaw_params['precip_output'], dir_input_base)
    fn_et_output = convert_to_absolute_path(detaw_params['et_output'], dir_input_base)
    fn_critical = convert_to_absolute_path(detaw_params['critical'], dir_input_base)
    fn_noncritical = convert_to_absolute_path(detaw_params['noncritical'], dir_input_base)

    # FIXME Avoid to use a current directory for jobs
    filepath = os.getcwd()

    model_start_year = int(start_water_year)-1
    # FIXME the start water year date of start_water_year-09-30 is a carry forward from the old code
    # setting this to year-10-01 and removing the extra input data results in numpy errors. I tried fixing
    # this by changing the dimensions of the declared variables but the results are different and hence I
    # rolled back to the old version.
    start_date_str = str(model_start_year) + '-09-30'
    end_date_str = str(end_water_year) + '-09-30'
    # convert string to datetime
    start_water_year_dt = pd.to_datetime(start_date_str)
    end_water_year_dt = pd.to_datetime(end_date_str)


    water_years = numpy.arange(start_water_year,end_water_year+1)
    # get endyear for the model run
    endyear = end_water_year_dt.year
    # get length in days for the model run
    idates = len(pd.date_range(start_water_year_dt, end_water_year_dt, freq='D'))
    # print("endyear =",endyear)
    print("idates =", idates)

    start1 = numpy.array([start_water_year_dt.year, start_water_year_dt.month, start_water_year_dt.day, 23, 0], dtype='i4')
    iyears = endyear-start1[0]+1
    print("iyears =", iyears)

    # Setting the value of ilands from the landuse file
    tmp_landuse_df = pd.read_csv(fn_landuse,header=[0])
    ilands = len(tmp_landuse_df['area_id'].unique())

    # Reading the input pcp file to get the number of pcp stations (isites). Reading the pcp values is performed in read_pcp
    tmp_df = pd.read_csv(fn_input_pcp,header=[0])
    # subtract columns by 4 to ingnore year, month, day, doy to get the number of pcp stations
    isites = tmp_df.shape[1]-4

    NumDay = numpy.array([0, 31, 28, 31, 30, 31, 30, 31,
                         31, 30, 31, 30, 31], dtype='i4')
    NI = numpy.array([31, 59, 90, 120, 151, 181, 212,
                     243, 273, 304, 334, 365], dtype='i4')
    NII = numpy.array([31, 60, 91, 121, 152, 182, 213,
                      244, 274, 305, 335, 366], dtype='i4')

    # XXX Need to fix hardwired stations
    #? all the above are hardwired. why?#
    ts_pcp = zeros((isites, idates), float)
    ts_per = zeros((isites, ilands), float)
    ETo_corrector = zeros((ilands), float)
    Region = zeros((ilands), int)
    ts_year = zeros((idates), int)
    ts_mon = zeros((idates), int)
    ts_days = zeros((idates), int)
    ts_LODI_tx = zeros((idates), float)
    ts_LODI_tn = zeros((idates), float)

    ts_pcp = read_pcp(start_date_str, end_date_str, fn_input_pcp)

    [ts_per, ETo_corrector, Region] = read_et_correction_factors(fn_et_correction)

    [ts_year, ts_mon, ts_days, ts_LODI_tx, ts_LODI_tn] = read_temperature(streamlinemodel, start_date_str, end_date_str, fn_input_temperature)

    [yearType, HAcre, icroptype] = read_landuse(fn_landuse, iyears, water_years, ilands)

    crdf = read_and_clean_crop_info(fn_critical)

    ncrdf = read_and_clean_crop_info(fn_noncritical)

    if DEBUG_TIMING:
        st = timeit.default_timer()
    (pcp, ET0) = weatheroutput(ts_pcp, ts_per, ts_mon, ts_days, ts_LODI_tx,
                               ts_LODI_tn, ilands, idates, isites, ETo_corrector, filepath, start1)
    weatheroutput_to_netcdf(pcp, ET0, model_start_year, fn_precip_output, fn_et_output)
    if DEBUG_TIMING:
        print('weather output took', timeit.default_timer()-st, ' seconds')
    pcp = pcp.T
    ET0 = ET0.T
    if DEBUG_TIMING:
        st = timeit.default_timer()
    # output dimensioned by (var, island,landuse,time)
    (DETAWOUTPUT) = historicalETAW(ts_per, ETo_corrector, Region, pcp, ET0, ts_LODI_tx, ts_LODI_tn,
                                   ilands, idates, isites, ts_year, ts_mon, ts_days, start1, filepath, NI, NII, NumDay, iyears,
                                   idayoutput, imonthoutput, iyearoutput, itotaloutput, dailyunit, forDSM2_daily, streamlinemodel,model_start_year,
                                   yearType,HAcre,icroptype,crdf,ncrdf)

    if DEBUG_TIMING:
        print('historical etaw calculations took ',
              timeit.default_timer()-st, ' seconds')
    if DEBUG_TIMING:
        st = timeit.default_timer()
    dx = write_to_netcdf(DETAWOUTPUT,model_start_year,fn_detaw_output)
    if DEBUG_TIMING:
        print('detaw output to netcdf4 took',
              timeit.default_timer()-st, ' seconds')

    print("done")


def main() -> None:
    """ Run ETAW and post-process ETAW results for DCD

        This function accepts one command line argument for model parameters,
        e.g.

        $ python detaw.py detaw.yaml

        Returns
        -------
        None
    """
    # Create an argument parser to get command line input
    parser = create_argparser()
    args = parser.parse_args()
    # Read the main yaml input file
    fname_main_yaml = args.input_yaml
    detaw(fname_main_yaml)


# _______________________________________________________________________________
if __name__ == "__main__":
    main()
