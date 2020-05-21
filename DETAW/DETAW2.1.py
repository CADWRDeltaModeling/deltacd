# Delta Evapotranspiration of Applied Water(DETAWv2.0)
#<license>
#    Copyright (C) State of California, Department of Water Resources.
#    This file is part of Delta Evapotranspiration of Applied Water
#    (DETAWv2.1).

#    DETAWv2.1 is free software: 
#    you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    DETAWv2.1 is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with DETAWv2.1.  If not, see <http://www.gnu.org/licenses>.
#</license>
#
# This version 2.1 was developed only for supporting CALSIM3 studies.
#
# Enter the following line in the command window to run the model:
#    Python detaw.py 
# 

import pyhecdss
import pandas as pd
from numpy import add,pi,array,zeros

import os, sys, string, math, numpy
from os import listdir
from math import cos,sin,tan,atan,sqrt,pi,pow

import for_DCD
from for_DCD import timeseries_combine, forNODCU

def write_to_dss(dssfh, arr, path, startdatetime, cunits, ctype):
    '''
    write to the pyhecdss.DSSFile for an array with starttime and assuming
    daily data with the pathname path, cunits and ctype
    '''
    fstr='1D'
    epart=path.split('/')[5]
    if epart == '1DAY':
        fstr='1D'
    elif epart == '1MONTH':
        fstr='1M'
    else:
        raise RuntimeError('Not recognized frequency in path: %s'%path)
    df=pd.DataFrame(arr,index=pd.date_range(startdatetime,periods=len(arr),freq=fstr))
    dssfh.write_rts(path,df,cunits,'INST-VAL')


#def list_add_list(list1,list2):
#    #list1 and list2 have the same length
#    listtemp = []
#    for i in range(0, len(list1)):
#        listtemp.append(list1[i]+list2[i])
#    return(listtemp)
    

def weatheroutput(ts_pcp,ts_per,ts_mon,ts_days,Tmax,Tmin,ilands,idates,isites,ETo_corrector,filepath,start1):
    """
        calculate the precipitation and reference evapotranpiration for each island

    input: daily precipitation at seven sites
           precentages of seven sites for 168 islands
    output: Nothing,
            Generate daily file and monthly weather file.
    """
    monthname = ["JAN","FEB","MAR","APR","MAY","JUN","JUL","AUG","SEP","OCT","NOV","DEC"] 
    
    istat = int(1)
    ifltab = zeros(600,"i")
    outputfile = filepath+"\\Output\\weather.dss"
    pyhecdss.set_message_level(0)
    pyhecdss.set_program_name('DETAW')
    dssfh=pyhecdss.DSSFile(outputfile)
    ctype = "INST-VAL" #"PER-AVER"
    iplan = int(0)
    
    ts_ptotal = zeros((ilands,idates),float)
    ET0Daily = zeros((ilands,idates),float)
    #start1 = [1921,9,30,23,0]  #datetime(1921,9,1,0,0)
    startdate = str(start1[2])+monthname[start1[1]-1]+str(start1[0])
    starttime = str(start1[3])+"00"
    ##The start time for monthly interval must be the first day of each month.
    
    for j in range(0,ilands):
        for k in range(0,idates):
                           
            for i in range(0,isites):
                ts_ptotal[j,k] = ts_ptotal[j,k] + ts_pcp[i,k]*ts_per[i,j]
                if i == (isites-1):
                    ## if ts_ptotal[j,k] < 0.0:
                    if ts_ptotal[j,k] < 0.0 or ts_ptotal[j,k] > 9990:  ##revised 3/19/2009
                        ts_ptotal[j,k] = 0.0
                        
            if (Tmax[k] < Tmin[k]) and (j==0):
                Tmax[k] = Tmin[k]

            ##Calculate ETo with Hargres-Samani equation
            ##PI = 3.1415926
            Lat = 38.5
            PHI = pi*Lat/180
            GSC = 0.082  ##Global Solar constant
            Lam = 2.45   ##Latent heat of Vaporization
            TDiffTemp = Tmax[k] - Tmin[k]
            Tm = (Tmax[k]+Tmin[k])/2

            df = 1 + 0.033*cos(2*pi/365*ts_days[k])  ##Earth-sun distance
            dec1 = 0.409*sin(2*pi/365*ts_days[k]-1.39)   ##Delination of sun
            arg = -tan(PHI)*tan(dec1)  ##argument for tangant in radians
            if arg*arg >= 1:           ## WS = sunrise hour angle in radians
                WS = pi/2.
            else:
                WS = pi/2.-atan(arg/sqrt(1-arg*arg))
            cosz = WS*sin(dec1)*sin(PHI)+(cos(dec1)*cos(PHI)*sin(WS))
            Ra = (24*60/pi)*GSC*df*cosz    ##Extrater Radiation
            ET0 = (0.0023*Ra*sqrt(TDiffTemp)*(Tm+17.8))/Lam
            if ET0 < 0.0:
                ET0 = 0.0
            ET0Daily[j,k] = ET0*ETo_corrector[j]
        
        path = "/detaw/island_"+str(j+1)+"/precipitation//1DAY/detaw/"
        write_to_dss(dssfh, ts_ptotal[j], path, startdate+" "+starttime, 'mm', ctype)
        #props={TIMESTAMP:PERIOD_START,AGGREGATION:MEAN,"UNIT":"mm"}    
        #tss_P = rts(ts_ptotal[j],start,dt,props)
        #path = "/detaw/island_"+str(j+1)+"/precipitation//1DAY/detaw/"
        #dss_store_ts(tss_P,outputfile,path)
        
        #tss_ET0 = rts(ET0Daily[j],start,days(1),props)
        path = "/detaw/island_"+str(j+1)+"/ET0//1DAY/detaw/"
        write_to_dss(dssfh, ET0Daily[j], path, startdate+" "+starttime, 'mm', ctype)
        #dss_store_ts(tss_ET0,outputfile,path)
        
        if j == 0:
            #props={TIMESTAMP:PERIOD_START,AGGREGATION:MEAN,"UNIT":"oC"}
            #tss_tx = rts(Tmax,start,days(1),props)
            path = "/detaw/LODI_Tmax/Temp//1DAY/detaw/"
            write_to_dss(dssfh, Tmax, path, startdate+" "+starttime, 'oC', ctype)
            #dss_store_ts(tss_tx,outputfile,path)
                       
            #tss_tn = rts(Tmin,start,days(1),props)
            path = "/detaw/LODI_Tmin/Temp//1DAY/detaw/"
            write_to_dss(dssfh, Tmin, path, startdate+" "+starttime, "oC", ctype)
            #dss_store_ts(tss_tn,outputfile,path)


        ##Generate monthly weather output: Tmax, Tmin, Pcp, ET0
        tempmonth = start1[1]
        mon_pcp = []
        mon_tx = []
        mon_tn = []
        mon_ET0 = []
        imon = 0
        sum_pcp = 0.0
        sum_tx = 0.0
        sum_tn = 0.0
        sum_ET0 = 0.0
        itemp = 0
        for k in range(0,idates):
            
            if ts_mon[k] == tempmonth:
                if k<(idates-1):
                    sum_pcp = sum_pcp+ts_ptotal[j,k]
                    sum_tx = sum_tx+Tmax[k]
                    sum_tn = sum_tn + Tmin[k]
                    sum_ET0 = sum_ET0 + ET0Daily[j,k]
                    itemp = itemp + 1
                else:
                    ##mon_pcp[imon].append(sum_pcp/itemp)
                    mon_pcp.append(sum_pcp)
                    mon_tx.append(sum_tx/itemp)
                    mon_tn.append(sum_tn/itemp)
                    mon_ET0.append(sum_ET0/itemp)
            else:
                mon_pcp.append(sum_pcp)
                mon_tx.append(sum_tx/itemp)
                mon_tn.append(sum_tn/itemp)
                mon_ET0.append(sum_ET0/itemp)
                imon = imon + 1
                if k < (idates-1):
                    itemp = 1
                    sum_pcp = ts_ptotal[j,k]
                    sum_tx = Tmax[k]
                    sum_tn = Tmin[k]
                    sum_ET0 = ET0Daily[j,k]
                    tempmonth = ts_mon[k]
                else:
                    mon_pcp.append(ts_ptotal[j,k])
                    mon_tx.append(Tmax[k])
                    mon_tn.append(Tmin[k])
                    mon_ET0.append(ET0Daily[j,k])
                    
        #props={TIMESTAMP:PERIOD_START,AGGREGATION:MEAN,"UNIT":"mm"}    
        #tss_mon_P = rts(mon_pcp,start1,dt_mon,props)
        path = "/detaw/island_"+str(j+1)+"/precipitation//1MONTH/detaw/"
        write_to_dss(dssfh, mon_pcp, path, startdate+" "+starttime, 'MM', ctype)
        #--delete--[ifltab,istat] = hecdss.zsrts(ifltab,path,startdate,starttime,len(mon_pcp),mon_pcp,"MM",ctype,iplan)
        #dss_store_ts(tss_mon_P,outputfile,path)
        
        #tss_mon_ET0 = rts(mon_ET0,start1,dt_mon,props)
        path = "/detaw/island_"+str(j+1)+"/ET0//1MONTH/detaw/"
        write_to_dss(dssfh, mon_ET0, path, startdate+" "+starttime, "MM", ctype)
        #dss_store_ts(tss_mon_ET0,outputfile,path)
        
        if j == 0:
            #props={TIMESTAMP:PERIOD_START,AGGREGATION:MEAN,"UNIT":"oC"}
            #tss_mon_tx = rts(mon_tx,start1,dt_mon,props)
            path = "/detaw/LODI_Tmax/Temp//1MONTH/detaw/"
            write_to_dss(dssfh, mon_tx, path, startdate+" "+starttime, "oC", ctype)
            #dss_store_ts(tss_mon_tx,outputfile,path)
            
            #tss_mon_tn = rts(mon_tn,start1,dt_mon,props)
            path = "/detaw/LODI_Tmin/Temp//1MONTH/detaw/"
            write_to_dss(dssfh, mon_tn, path, startdate+" "+starttime, "oC", ctype)
            #dss_store_ts(tss_mon_tn,outputfile,path)
    
    dssfh.close()        
    return(ts_ptotal, ET0Daily)
        
                          
##_______________________________________________________________________________
##_______________________________________________________________________________

def historicalETAW(ts_per,ETo_corrector,Region,pcp,ET0,tmax,tmin,ilands,idates,isites,ts_year,ts_mon,\
                   ts_days,start1,filepath,NI,NII,NumDaysPerMon,iyears,idayoutput,imonthoutput,\
                   iyearoutput,itotaloutput,dailyunit,forDSM2_daily,streamlinemodel):
    InpHSACrop = "  "
    SACropDaily = "  "
    HSACropDailyMean = "  "
    HSACropMonMean = "  "
    Date = "  "
    cpartt = "  "
    monthname = ["JAN","FEB","MAR","APR","MAY","JUN","JUL","AUG","SEP","OCT","NOV","DEC"] 
    CropName = ["Urban","Irrig pasture","Alfalfa","All Field","Sugar beets",  \
                "Irrig Grain","Rice", "Truck Crops", "Tomato","Orchard",  \
                "Vineyard", "Riparian Vegetation","Native Vegetation", \
                "Non-irrig Grain","Water Surface"]
    icroptype = 15
    ##iyears = ts_year[len(ts_year)-1]-start1.year
    ##iyears = 2008-start1.year+1
    idays = 366
    imonths = 12
    ##idays = int(idates/iyears)
    
    
    istatd = int(1)
    istatm = int(1)
    istaty = int(1)
    ifltabd = zeros(600,"i")
    ifltabm = zeros(600,"i")
    ifltaby = zeros(600,"i")
    pyhecdss.set_message_level(0)
    pyhecdss.set_program_name('DETAW')
    ctype = "INST-VAL" # "PER-AVER"
    iplan = int(0)
    start2 = [1921,10,1,23,0]
    #startdate = str(start1[2])+monthname[start1[1]-1]+str(start1[0])
    startdate = str(start2[2])+monthname[start2[1]-1]+str(start2[0])
    #starttime = str(start1[3])+"00"
    starttime = str(start2[3])+"00"
    
    
    dpyAll = zeros((iyears+2),int)
    IKc = zeros((iyears+2,idays+1),float)
    Kc = zeros((iyears+2,idays+1),float)
    OKc = zeros((iyears+2,idays+1),float)

    EToDaily = zeros((iyears+2,idays+1),float)
    PcpDaily = zeros((iyears+2,idays+1),float)
    ETcDaily = zeros((iyears+2,idays+1),float)
    Date1Daily = "  "*(iyears+2)*(idays+1)
    
    
    ##DavisPcp = zeros((iyears,idays),float)
    WSCESpg = zeros((iyears+2,idays+1),float)
    ##NumDaysPerMon = [0,31,28,31,30,31,30,31,31,30,31,30,31]
    
    HAcre = zeros((ilands+1,iyears+2,icroptype+1),float)
    Kc1 = zeros((icroptype+1),float)
    Kc2 = zeros((icroptype+1),float)
    Kc3 = zeros((icroptype+1),float)
    AB = zeros((icroptype+1),float)
    AC = zeros((icroptype+1),float)
    AD = zeros((icroptype+1),float)
    SDx = zeros((icroptype+1),float)
    RDxU = zeros((icroptype+1),float)
    RDxL = zeros((icroptype+1),float)
    awL = zeros((icroptype+1),float)
    awU = zeros((icroptype+1),float)
    ADep = zeros((icroptype+1),float)
    ## crop info for critical years
    CBeginDate = zeros((icroptype+1),int)
    CEndDate = zeros((icroptype+1),int)
    CCropType = zeros((icroptype+1),int)
    Cf = zeros((icroptype+1),int)
    Ckc1 = zeros((icroptype+1),float)
    Ckc2 = zeros((icroptype+1),float)
    Ckc3 = zeros((icroptype+1),float)
    CAB = zeros((icroptype+1),float)
    CAC = zeros((icroptype+1),float)
    CAD = zeros((icroptype+1),float)
    CSDx = zeros((icroptype+1),float)
    CRDxU = zeros((icroptype+1),float)
    CRDxL = zeros((icroptype+1),float)
    CawL = zeros((icroptype+1),float)
    CawU = zeros((icroptype+1),float)
    CADep = zeros((icroptype+1),float)
    ## crop info for non-critical years
    NCBeginDate = zeros((icroptype+1),int)
    NCEndDate = zeros((icroptype+1),int)
    NCCropType = zeros((icroptype+1),int)
    NCf = zeros((icroptype+1),int)
    NCkc1 = zeros((icroptype+1),float)
    NCkc2 = zeros((icroptype+1),float)
    NCkc3 = zeros((icroptype+1),float)
    NCAB = zeros((icroptype+1),float)
    NCAC = zeros((icroptype+1),float)
    NCAD = zeros((icroptype+1),float)
    NCSDx = zeros((icroptype+1),float)
    NCRDxU = zeros((icroptype+1),float)
    NCRDxL = zeros((icroptype+1),float)
    NCawL = zeros((icroptype+1),float)
    NCawU = zeros((icroptype+1),float)
    NCADep = zeros((icroptype+1),float)

    BeginDateYear = zeros((iyears+2),int)    

    osIkc = zeros((iyears+2),float)
    isIkc = zeros((iyears+2),float)
    IGETo = zeros((iyears+2),float)
    osFkc = zeros((iyears+2),float)

    KcByr = zeros((iyears+2),float)
    KcCyr = zeros((iyears+2),float)
    KcDyr = zeros((iyears+2),float)
    KcEyr = zeros((iyears+2),float)
    
    CCKc = zeros((idays+1),float)
    LIYear = zeros((iyears+2),int)
    BIYear = zeros((iyears+2),int)
    NA1 = zeros((iyears+2),float)
    NA2 = zeros((iyears+2),float)
    NA3 = zeros((iyears+2),float)
    isCETc = zeros((iyears+2),float)
    osCETc = zeros((iyears+2),float)
    isCERn = zeros((iyears+2),float)
    osCERn = zeros((iyears+2),float)
    isPCP = zeros((iyears+2),float)
    isETaw = zeros((iyears+2),float)
    
    isMCETc = zeros((imonths+1),float)
    osMCETc = zeros((imonths+1),float)
    isMCERn = zeros((imonths+1),float)
    osMCERn = zeros((imonths+1),float)
    isCSpg= zeros((iyears+2),float)
    osCSpg= zeros((iyears+2),float)
    isMCSpg= zeros((imonths+1),float)
    osMCSpg= zeros((imonths+1),float)
    EToMonthly= zeros((imonths+1),float)
    PcpMonthly= zeros((imonths+1),float)
    ETcMonthly= zeros((imonths+1),float)
    ERnMonthly= zeros((imonths+1),float)
    SpgMonthly= zeros((imonths+1),float)
    EspgMonthly= zeros((imonths+1),float)
    WSEspgMonthly= zeros((iyears+2,imonths+1),float)
    MonDsw = zeros((imonths+1),float)
    MonDswPos = zeros((imonths+1),float)
    MonNetApp = zeros((imonths+1),float)
    yearType = []  ##"  "*(iyears+1)
    
    
    ## for HSA****.csv (not for OLDHSA****.csv)
    SACropDaily = "  "
    HSACropMonMean = "  "
    OldHSACropMonMean = "  "
    InpOldHSACropMonMean = "  "
    fpOldHSACrpMonMean = "  "
    ##4/20/09 yearTypeDaily = [] ##"  "*(idays+1)
    DateDaily = "  "*(idays+1)
    ##yearType = "  "
    MonDay = "  "
    flagETAW = "  "
    
    yDaily = zeros((idays+1),int)
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
    DETAWOUTPUT = zeros((6,ilands+1,icroptype+1,idates-1),float)
    
 ## for irrigation and hydrology year convertion (((((
    ytemp = zeros((idays+1),int)
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
    HAcretemp = zeros((idays+1),float)
## for irrigation and hydrology year convertion))))))



    EToMonthly = zeros((imonths+1), float)
    PcpMonthly = zeros((imonths+1), float)
    ETcMonthly = zeros((imonths+1), float)
    ERnMonthly = zeros((imonths+1), float)
    SpgMonthly = zeros((imonths+1), float)
    EspgMonthly = zeros((imonths+1), float)
    MonETAW = zeros((iyears+2,imonths+1), float)
    ETAWMonDay = zeros((imonths+1,32), float)


    
    
    DOYLIrrig = zeros((iyears+2), float)
    DOYGrainLIrrig = zeros((iyears+2), float)
    DOYLIrrig[0] = 0
    DOYGrainLIrrig[0] = 0
    ##for HSA****.csv  (not OLDHSA***.csv)

    ## Add more variables for DETAW python program to efficient work    
    ##cropareas = zeros((ilands,icroptype,iyears),float)
    ## end:Add more variables for DETAW python program to efficient work
    
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
    
    for i in range(0,idays+1):
        CCKc[i] = 0.0

    for i in range(0,iyears+1):
        osIkc[i]=0
        isIkc[i]=0
        IGETo[i]=0
        osFkc[i]=0
        KcByr[i]=0
        KcCyr[i]=0
        KcDyr[i]=0
        KcEyr[i]=0
        dpyAll[i]=0
        ## yrs[i]=0
        NA1[i]=0
        NA2[i]=0
        NA3[i]=0
        for j in range(0,idays+1):
            EToDaily[i,j]=0
            ##Date1Daily[i,j]=0
            PcpDaily[i,j]=0
            OKc[i,j]=0
            IKc[i,j]=0
            Kc[i,j]=0
            
    ##++++++++++++++++++++++++++++++++++++++++++++++++++++
    ##
    ## Determine ETo Rain base bare soil evaporation
    ##
    ##++++++++++++++++++++++++++++++++++++++++++++++++++++
    ##step1: read percentage file: percentage, ET0_corrector, and region
            
    ##step2: read \weatheroutput\ files, data:DayOfYear, TMax, TMin, PCP, ETo            
    ## The weather data come from the main program

    ##step3: read crop information for critical years
    if streamlinemodel == "CALSIM3":
        source = filepath+"\\Input\\planning_study\\critical.csv"
    else:
        source = filepath+"\\Input\\historical_study\\critical.csv"
    ts_type = "rts"
    ##start1 = datetime(1921,9,29,23,0)
    ##intl = time_interval(days=1)
    
    fint = open(source,'r')
    icon = 0
    for oneline in fint:
        if icon == 9:
            for j in range(1,icroptype+1):
                CCropType[j] = int(oneline.split(",")[j+1])
        #selector = "critical$C10:Q10"
        #tss_temp = excel_retrieve_ts(source,selector,ts_type,start=start1,interval=intl)
        #for j in range(1,icroptype+1):
        #    CCropType[j] = int(tss_temp[j-1].data[0])
               
        if icon == 12:
            for j in range(1,icroptype+1):
                CBeginDate[j] = int(oneline.split(",")[j+1])
        if icon == 13:
            for j in range(1,icroptype+1):
                CEndDate[j] = int(oneline.split(",")[j+1])
        if icon == 14:
            for j in range(1,icroptype+1):
                Cf[j] = int(oneline.split(",")[j+1])
        if icon == 15:
            for j in range(1,icroptype+1):
                Ckc1[j] = float(oneline.split(",")[j+1])
        if icon == 16:
            for j in range(1,icroptype+1):
                Ckc2[j] = float(oneline.split(",")[j+1])
        if icon == 17:
            for j in range(1,icroptype+1):
                Ckc3[j] = float(oneline.split(",")[j+1])
        if icon == 18:
            for j in range(1,icroptype+1):
                CAB[j] = float(oneline.split(",")[j+1])
        if icon == 19:
            for j in range(1,icroptype+1):
                CAC[j] = float(oneline.split(",")[j+1])
        if icon == 20:
            for j in range(1,icroptype+1):
                CAD[j] = float(oneline.split(",")[j+1])
        if icon == 21:
            for j in range(1,icroptype+1):
                CSDx[j] = float(oneline.split(",")[j+1])
        if icon == 22:
            for j in range(1,icroptype+1):
                CRDxL[j] = float(oneline.split(",")[j+1])
        if icon == 23:
            for j in range(1,icroptype+1):
                CRDxU[j] = float(oneline.split(",")[j+1])
        if icon == 24:
            for j in range(1,icroptype+1):
                CawL[j] = float(oneline.split(",")[j+1])
        if icon == 25:
            for j in range(1,icroptype+1):
                CawU[j] = float(oneline.split(",")[j+1])
        if icon == 26:
            for j in range(1,icroptype+1):
                CADep[j] = float(oneline.split(",")[j+1])
        icon += 1    
    fint.close()
    
    ##step4: read crop information for non-critical years
    if streamlinemodel == "CALSIM3":
        source = filepath+"\\Input\\planning_study\\noncritical.csv"
    else:
        source = filepath+"\\Input\\historical_study\\noncritical.csv"
    fint = open(source,'r')
    icon = 0
    for oneline in fint:
        if icon == 9:
            for j in range(1,icroptype+1):
                NCCropType[j] = int(oneline.split(",")[j+1])
        
        if icon == 12:
            for j in range(1,icroptype+1):
                NCBeginDate[j] = int(oneline.split(",")[j+1])
        if icon == 13:
            for j in range(1,icroptype+1):
                NCEndDate[j] = int(oneline.split(",")[j+1])
        if icon == 14:
            for j in range(1,icroptype+1):
                NCf[j] = int(oneline.split(",")[j+1])
        if icon == 15:
            for j in range(1,icroptype+1):
                NCkc1[j] = float(oneline.split(",")[j+1])
        if icon == 16:
            for j in range(1,icroptype+1):
                NCkc2[j] = float(oneline.split(",")[j+1])
        if icon == 17:
            for j in range(1,icroptype+1):
                NCkc3[j] = float(oneline.split(",")[j+1])
        if icon == 18:
            for j in range(1,icroptype+1):
                NCAB[j] = float(oneline.split(",")[j+1])
        if icon == 19:
            for j in range(1,icroptype+1):
                NCAC[j] = float(oneline.split(",")[j+1])
        if icon == 20:
            for j in range(1,icroptype+1):
                NCAD[j] = float(oneline.split(",")[j+1])
        if icon == 21:
            for j in range(1,icroptype+1):
                NCSDx[j] = float(oneline.split(",")[j+1])
        if icon == 22:
            for j in range(1,icroptype+1):
                NCRDxL[j] = float(oneline.split(",")[j+1])
        if icon == 23:
            for j in range(1,icroptype+1):
                NCRDxU[j] = float(oneline.split(",")[j+1])
        if icon == 24:
            for j in range(1,icroptype+1):
                NCawL[j] = float(oneline.split(",")[j+1])
        if icon == 25:
            for j in range(1,icroptype+1):
                NCawU[j] = float(oneline.split(",")[j+1])
        if icon == 26:
            for j in range(1,icroptype+1):
                NCADep[j] = float(oneline.split(",")[j+1])
        icon += 1    
    fint.close()
   
    ##step5: read land use from .\Landuse folder !!!!!!!Not checked 3/13/2009
             ##get year type of each year
    
    if streamlinemodel == "CALSIM3":
        source = filepath +"\\Input\\planning_study\\Landuse\\SA0001.csv"    
    else:
        source = filepath +"\\Input\\historical_study\\Landuse\\SA0001.csv" 
    f0 = open(source)
    iline0 = 1
    yearType.append("AN")
    for line in f0:
        if line:
            if line[0] == "1" or line[0] == "2":
                yearType.append(line.split(",")[1].strip())                
    yearType.append("AN")  ##for last year
                
            
    ## get Hectares of each crop type, year and island        
    
    if streamlinemodel == "CALSIM3":
        hist_path = filepath + "\\Input\\planning_study\\Landuse\\"      ##---08/02/2010
    else:
        hist_path = filepath + "\\Input\\historical_study\\Landuse\\"
    files = listdir(hist_path)
    ts_type = "rts"
    ##intl = time_interval(years=1)
    ##start0 = datetime(1922,1,1,1,0)
    for file in files:
        if ".csv" in file:
            sheetname = file.replace(".csv","")
            ilandno = int(sheetname.split("A0")[1])  ##Landuse ---08/02/2010
            source = hist_path + file
            fint = open(source,'r')
            icon = 0
            for oneline in fint:
                #print icon, oneline
                if icon > 1:
                    for j in range(1,icroptype+1):
                        HAcre[ilandno-1,icon-1,j] = float(oneline.split(",")[j+1])
                        #print "HAcre =", HAcre[ilandno-1,icon-1,j]
                icon += 1
            fint.close()
    ##End of data input
    #--for itt in range(0,icon):
    #--    print HAcre[167,itt,3]
                    
    ##+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ##
    ##Calculating ETAW for each crop in SA
    ##
    ##+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    
    for k in range(0,ilands):
        print("island =",k+1)
        PETo = 0.0
        DCT = 1
        
        dpyAll[0] = 365
        for kk in range(0,idates):
            iy = ts_year[kk] - start1[0] + 1
            id = ts_days[kk]
            
            if (ET0[k,kk]*100-int(ET0[k,kk]*100))>0.5:
                EToDaily[iy,id] = int(ET0[k,kk]*100+1)/100.
            else:
                EToDaily[iy,id] = int(ET0[k,kk]*100)/100.
            
            if (pcp[k,kk]*100-int(pcp[k,kk]*100))>0.5:
                PcpDaily[iy,id] = int(pcp[k,kk]*100+1)/100.
            else:
                PcpDaily[iy,id] = int(pcp[k,kk]*100)/100.
            
            dpy=365
            if ts_year[kk]%4==0:
                dpy=366
            dpyAll[iy]= dpy
                        
            DCT = DCT + 1
            if EToDaily[iy,id] <= 0:
                EToDaily[iy,id] = PETo
            CETo = CETo+EToDaily[iy,id]
            METo = CETo/DCT
            if PcpDaily[iy,id] > METo:
                CETo = EToDaily[iy,id] 
                DCT = 1
                PEs = 0
            METo = CETo/DCT
            kx = 1.22 - 0.04*METo
            CEx = kx*CETo
            if CEx < 0:
                CEx = 0
            if pow(CEx,0.5) < Beta1:
                CEs = CEx
            else:
                CEs = Beta1*pow(CEx,0.5)
            Es = CEs - PEs
            
            if EToDaily[iy,id] != 0:
                OKc[iy,id] = CEs/CETo
            
            PETo = EToDaily[iy,id]
            PEs = CEs

        for y in range(0, iyears+1):
            for dd in range(0, idays+1):
                WSCESpg[y,dd] = 0
            for Mon in range(0, imonths):
                WSEspgMonthly[y,Mon] = 0
            
              
        for j in range (1,icroptype+1):
            ##for HSA*** (not for OLDHSA***)   +++++++++++++++++++++  
            ##initalize Monthly ETaw
            #print("island =",k+1, " croptype =",j)
            
            for y in range(0,iyears+2):
                for Mon in range(0, 13):
                    MonETAW[y,Mon]=0
            ##initialize cumulative ETAW
            CETAWDaily=0
            ##initializa DoyLast irrigation
            for y in range(0,iyears+1):
                DOYLIrrig[y]=0
                DOYGrainLIrrig[y]=0
            ##for HSA*** (not for OLDHSA***)  +++++++++++++++++++++
            dataday = []            
            data2day = []
            data3day = []
            data4day = []
            data5day = []
            data6day = []
            data7day = []
            data8day = []
            data9day = []
            data10day = []
            data11day = []
            data12day = []
            data13day = []
            data14day = []
            data15day = []
            data16day = []
            data17day = []
            data18day = []
            data19day = []
            data20day = []
            data21day = []
            data22day = []
            data23day = []
            data24day = []
            data25day = []
            data26day = []
            data27day = []
            data28day = []
            data29day = []
            data30day = []
            data31day = []
            data32day = []
            data1day_14crops = []
            data10day_14crops = []
            data11day_14crops = []
            data12day_14crops = []
            data13day_14crops = []
            data1day_water = []
            data10day_water = []
            data11day_water = []
            data12day_water = []
            data13day_water = []
            
            datamon = []
            data2mon = []
            data3mon = []
            data4mon = []
            data5mon = []
            data6mon = []
            data7mon = []
            data8mon = []
            data9mon = []
            data10mon = []
            data11mon = []
            data12mon = []
            data13mon = []
            data14mon = []
            data15mon = []
            data16mon = []
            data17mon = []
            data18mon = []
            data19mon = []
            data20mon = []
            data21mon = []
            
            
            datayr = []
            data2yr = []
            data3yr = []
            data4yr = []
            data5yr = []
            data6yr = []
            data7yr = []
            data8yr = []

            for i in range(0,iyears):
                SWD = 0.0
                PSWD = 0.0               
                if yearType[i].upper()=="C" or yearType[i].upper()=="D":
                    CropType1=CCropType[j]
                    BeginDate1=CBeginDate[j]
                    EndDate1=CEndDate[j]
                    f1=Cf[j]
                    kkc1=Ckc1[j]
                    kkc2=Ckc2[j]
                    kkc3=Ckc3[j]
                    AB1=CAB[j]
                    AC1=CAC[j]
                    AD1=CAD[j]
                    SDx1=CSDx[j]
                    CRDxU1=CRDxU[j]
                    CRDxL1=CRDxL[j]
                    CawL1=CawL[j]
                    CawU1=CawU[j]
                    ADep1=CADep[j]
                    ##Set RDX1 to RDxU or RDxL
                    if Region[k] == 0:
                        RDx1=CRDxL1
                        aw1=CawL1
                    else:
                        RDx1=CRDxU1
                        aw1=CawU1
                else:
                    CropType1=NCCropType[j]
                    BeginDate1=NCBeginDate[j]
                    EndDate1=NCEndDate[j]
                    f1=NCf[j]
                    kkc1=NCkc1[j]
                    kkc2=NCkc2[j]
                    kkc3=NCkc3[j]
                    AB1=NCAB[j]
                    AC1=NCAC[j]
                    AD1=NCAD[j]
                    SDx1=NCSDx[j]
                    NCRDxU1=NCRDxU[j]
                    NCRDxL1=NCRDxL[j]
                    NCawL1=NCawL[j]
                    NCawU1=NCawU[j]
                    ADep1=NCADep[j]
                    ##Set RDX1 to RDxU or RDxL
                    if Region[k] == 0:
                        RDx1=NCRDxL1
                        aw1=NCawL1
                    else:
                        RDx1=NCRDxU1
                        aw1=NCawU1
                ##set the 20 days before end date for rice
                if j==7:
                    RiceIrrigEndDate=EndDate1-20
                ##set the default irrigation frequency to 30
                if RDx1<SDx1:
                    erd=RDx1
                else:
                    erd=SDx1
                if f1==0:
                    f1=30
                PAW1=erd*aw1
                ##YTD=PAW*(ADep/100)
                YTD=PAW1*(ADep1/100)

                osSWDx=aw1*300*ADep1/100      ##off season max soil water depletion

                ##Initilal Kc values for dates B,C,D,E
                KcB=kkc1
                KcC=kkc2
                KcD=kkc2
                KcE=kkc3

                ##..........................................................
                ##.. Calculate day of year corresponding to A,B,C, and E ...
                ## .. Note that B, C, D, and E can be bigger than 365   ....
                ## .. Note that BB,CC,DD, and EE are always <=365        ...
                ##..........................................................
                if EndDate1<=BeginDate1:
                    EndDate1=365+EndDate1
                lenDate=EndDate1-BeginDate1
                B=int(0.01*AB1*lenDate+BeginDate1)
                C=int(0.01*AC1*lenDate+BeginDate1)
                D=int(0.01*AD1*lenDate+BeginDate1)

                BB=B
                if B>365:
                    BB=B-365
                CC=C
                if C>365:
                    CC=C-365
                DD=D
                if D>365:
                    DD=D-365
                EE=EndDate1
                if EndDate1>365:
                    EE=EndDate1-365
                ##...............................................................
                ## The following loop calculates Kc values
                ##...............................................................
                ## for y in range(1, iyears+1):
                ## Jan29-2007
                ## starting date of grain and Non-irrig Greain is october
                y = i+1
                if j==6 or j==14:
                    yearTypeCal=yearType[y]
                else:
                    yearTypeCal=yearType[y-1]
                ##critical year
                if yearTypeCal.upper()=="C" or yearTypeCal.upper()=="D":
                    BeginDate1=CBeginDate[j]
                    EndDate1=CEndDate[j]
                    kkc1=Ckc1[j]
                    kkc2=Ckc2[j]
                    kkc3=Ckc3[j]
                    AB1=CAB[j]
                    AC1=CAC[j]
                    AD1=CAD[j]
                else:
                    BeginDate1=NCBeginDate[j]
                    EndDate1=NCEndDate[j]
                    kkc1=NCkc1[j]
                    kkc2=NCkc2[j]
                    kkc3=NCkc3[j]
                    AB1=NCAB[j]
                    AC1=NCAC[j]
                    AD1=NCAD[j]

                if j==7:
                    RiceIrrigEndDate=EndDate1-20
                ##Initial Kc values for dates B,C,D,E
                KcB=kkc1
                KcC=kkc2
                KcD=kkc2
                KcE=kkc3
                ##..........................................................
                ##.. Calculate day of year corresponding to A,B,C, and E ...
                ## .. Note that B, C, D, and E can be bigger than 365   ....
                ## .. Note that BB,CC,DD, and EE are always <=365        ...
                if EndDate1<=BeginDate1:
                    EndDate1=365+EndDate1
                lenDate=EndDate1-BeginDate1
                B=int(0.01*AB1*lenDate+BeginDate1)
                C=int(0.01*AC1*lenDate+BeginDate1)
                D=int(0.01*AD1*lenDate+BeginDate1)

                BB=B
                if B>365:
                    BB=B-365
                CC=C
                if C>365:
                    CC=C-365
                DD=D
                if D>365:
                    DD=D-365
                EE=EndDate1
                if EndDate1>365:
                    EE=EndDate1-365

                ##end of jan 29-2007
                osIkc[y]=0
                ctn=0
                IGKc=0
                IGETo1=0
                Q=BeginDate1
                R=B
                ##because in this program I use int number for croptypw I change the
                ##condition from cropType>2 to cropType>=2
                if CropType1 > 2:
                    R=B+10
                for ii in range(Q, R+1):
                    ctn=ctn+1
                    IGKc=IGKc+OKc[y,ii]
                    IGETo1=IGETo1+EToDaily[y,ii]

                             
                osIkc[y]=IGKc/ctn        ##initial growth Kc from off-season
                IGETo[y]=IGETo1/ctn      ##initial growth mean Eto rate
                ## Identify osFkc for Kc on date E
                osFkc[y]=0
                ctn=0
                EKc=0
                R=EndDate1
                if EndDate1>365:
                    R=EndDate1-365
                if R==1 and EndDate1>=365:
                    R=365
                Q=R-10

                for ii in range(Q,R+1):
                    ctn=ctn+1
                    EKc=EKc+OKc[y,ii]

              
                osFkc[y]=EKc/ctn  ## final Kc on date E from off-season

                ##.. Identify initial growth Kc from irrig freq(F)
                isIkc[y]=0
                CETo=f1*IGETo[y]

                ##RichEdit1->Lines->Add("CETo="+FloatToStr(CETo))
                ##kx=1.05-0.03*IGETo[y]
                kx=1.22-0.04*IGETo[y]
                CEx=kx*CETo
                if pow(CEx,0.5) <Beta1:
                    CEs=CEx 
                else:
                    CEs=Beta1*pow(CEx,0.5)

                if CETo != 0:
                    isIkc[y]=CEs/CETo

                if CropType1 >1:
                    isIkc[y]=0

            ## end of for  y<iyears  ********************************

            ## Identify initial growth Kc(KcB) and final Kc(KcE)
            LowIkc=2
            LowFkc=2
            for y in range(1,iyears+1):
                if LowIkc>osIkc[y]:
                    LowIkc=osIkc[y]
                
                if LowFkc>osFkc[y]:
                    LowFkc=osFkc[y]
                    
            ##......... Print Kc value selection process
            ## ......   Note that Kc's are not printed in final version
            ##........  Print rows were changed to remarks. however, the
            ## loop is retained to assign Kc's to subscripts.
            ##...........................................................

            ## Loop through years to calculate daily Kc values
            for y in range(1,iyears+1):
                ## Jan29-2007

                ## starting date of grain and Non-irrig Greain is october
                if j==6 or j==14:
                    yearTypeCal=yearType[y]                    
                else:
                    yearTypeCal=yearType[y-1]
                    
                ##critical years
                if yearTypeCal.upper()=="C" or yearTypeCal.upper()=="D": 
                    BeginDate1=CBeginDate[j]
                    EndDate1=CEndDate[j]
                    kkc1=Ckc1[j]
                    kkc2=Ckc2[j]
                    kkc3=Ckc3[j]
                    AB1=CAB[j]
                    AC1=CAC[j]
                    AD1=CAD[j]
                    ##end of critical years
                else:   ##non-critical yeras
                    BeginDate1=NCBeginDate[j]
                    EndDate1=NCEndDate[j]
                    kkc1=NCkc1[j]
                    kkc2=NCkc2[j]
                    kkc3=NCkc3[j]
                    AB1=NCAB[j]
                    AC1=NCAC[j]
                    AD1=NCAD[j]
                    ##we have to set RDX1 to RDxU or RDxL
                    ##end of non-critical yeras

                ##set the 20 days before end date for rice
                if j==7:
                    RiceIrrigEndDate=EndDate1-20

                ## Initilal Kc values for dates B,C,D,E
                KcB=kkc1
                KcC=kkc2
                KcD=kkc2
                KcE=kkc3

                ##..........................................................
                ##.. Calculate day of year corresponding to A,B,C, and E ...
                ## .. Note that B, C, D, and E can be bigger than 365   ....
                ## .. Note that BB,CC,DD, and EE are always <=365        ...
                ##..........................................................
                if EndDate1<=BeginDate1:
                    EndDate1=365+EndDate1
                lenDate=EndDate1-BeginDate1
                B=int(0.01*AB1*lenDate+BeginDate1)
                C=int(0.01*AC1*lenDate+BeginDate1)
                D=int(0.01*AD1*lenDate+BeginDate1)

                BB=B
                if B>365:
                    BB=B-365
                CC=C
                if C>365:
                    CC=C-365
                DD=D
                if D>365:
                    DD=D-365
                EE=EndDate1
                if EndDate1>365:
                    EE=EndDate1-365
                         
                if CropType1 > 1:
                    if CropType1 > 2:
                        if CropType1 <= 3:
                            KcB=LowIkc
                            KcE=kkc3
                    ## end of if CropType>2
                    else:   ## CropType>2
                        KcB=kkc1
                        KcC=kkc2
                        KcD=kkc2
                        KcE=kkc3
                        ##RichEdit1->Lines->Add("SaraKcBKcc")
                    ## end of CropType >2
                ## end of if croptype>1
                else:   ## CropType>1
                    KcB=LowIkc
                    KcC=kkc2
                    KcD=kkc2
                    KcE=kkc3

                    if isIkc[y] >KcB:
                        KcB=isIkc[y]

                ## end of else CropType >1
                KcByr[y]=KcB
                KcCyr[y]=KcC
                KcDyr[y]=KcD
                KcEyr[y]=KcE
            ## end of for y<YCt



            ##*****************************************************
            ##
            ##Write Daily results in OldHSA****C***.csv (OutSACropDaily)
            ##Write yearly results(Daily mean) in HSA****C***.eaw.csv (OutHSACropDailyMean)
            ##Write Monthly results in OldHSA****C***.mtv.csv (OutHSACropMonMean)
            ##
            ##*****************************************************
            
            ## This loop determines daily IKc for each year and subscripts
            ## the results by year and day(ends after 5020)
            for y in range(1,iyears+1):
                
                ## Jan29-2007

                ## starting date of grain and Non-irrig Greain is october
                if j==6 or j==14:
                    yearTypeCal=yearType[y]
                else:
                    yearTypeCal=yearType[y-1]
                ##critical years
                if yearTypeCal.upper()=="C" or yearTypeCal.upper()=="D":
                    BeginDate1=CBeginDate[j]
                    EndDate1=CEndDate[j]
                    kkc1=Ckc1[j]
                    kkc2=Ckc2[j]
                    kkc3=Ckc3[j]
                    AB1=CAB[j]
                    AC1=CAC[j]
                    AD1=CAD[j]
                    ##end of critical years
                else:   ##non-critical yeras
                    BeginDate1=NCBeginDate[j]
                    EndDate1=NCEndDate[j]
                    kkc1=NCkc1[j]
                    kkc2=NCkc2[j]
                    kkc3=NCkc3[j]
                    AB1=NCAB[j]
                    AC1=NCAC[j]
                    AD1=NCAD[j]
                    ##we have to set RDX1 to RDxU or RDxL
                    ##end of non-critical yeras
                    
                ##set the 20 days before end date for rice
                if j==7:
                    RiceIrrigEndDate=EndDate1-20

                ## Initilal Kc values for dates B,C,D,E
                KcB=kkc1
                KcC=kkc2
                KcD=kkc2
                KcE=kkc3

                ##..........................................................
                ##.. Calculate day of year corresponding to A,B,C, and E ...
                ## .. Note that B, C, D, and E can be bigger than 365   ....
                ## .. Note that BB,CC,DD, and EE are always <=365        ...
                ##..........................................................
                if EndDate1<=BeginDate1:
                    EndDate1=365+EndDate1
                lenDate=EndDate1-BeginDate1
                B=int(0.01*AB1*lenDate+BeginDate1)
                C=int(0.01*AC1*lenDate+BeginDate1)
                D=int(0.01*AD1*lenDate+BeginDate1)

                BB=B
                if B>365:
                    BB=B-365
                CC=C
                if C>365:
                    CC=C-365
                DD=D
                if D>365:
                    DD=D-365
                EE=EndDate1
                if EndDate1>365:
                    EE=EndDate1-365
                ##end of jan 29-2007
          
                for ii in range(0, idays+1):
                    IKc[y,ii] = 0

                dpy = dpyAll[y]
                ##KcE=KcEyr[y]
                KcB=kkc1
                KcC=kkc2
                KcD=kkc2
                KcE=kkc3
                BCslope=(KcC-KcB)/(C-B)
                CDslope=(KcD-KcC)/(D-C)
                DEslope=(KcE-KcD)/(EndDate1-D)
                BeginDateYear[y]=BeginDate1  ##project doesn't have this line
                
               
                for jj in range(BeginDate1,EndDate1+1):
                    ii = jj 
                    if jj>dpy:
                        ii=jj-dpy
                    if jj<=D:
                        if jj<=C:
                            if jj<=B:
                                IKc1 = KcB
                            else:
                                IKc1 = IKc1+BCslope
                        else:
                            IKc1 = IKc1 + CDslope
                    else:
                        IKc1 = IKc1 + DEslope

                    IKc[y,ii] = IKc1
                    
                if dpyAll[y]==366:
                    IKc[y,366]=IKc[y,365] 
                
                ## This section calculates daily Kc and print the results....
                ##apply startig and ending dates of Stress factor    

                ##IKcs = KcD*ks   ##Ikcs max kc with stress coeficcient based on date D
                                ## RichEdit1->Lines->Add(dpyAll[y])               
                
                for ii in range(1,dpyAll[y]+1):
                    
                    IKc1 = IKc[y,ii]
                    OKc1 = OKc[y,ii]
                                       
                    Kc11 = IKc1
                    if OKc1 < Kc11:
                        Kc11 = OKc1

                    if CropType1 >= 3:
                        ##Kc11 = Kc11 + CCKc[ii]
                        if IKc[y,ii] != 0:
                            if Kc11>1.15:
                                Kc11 = 1.15
                        else:
                            if Kc11>1.05:
                                Kc11 = 1.05

                        if CCKc[ii]==CC1 and Kc11<0.9:
                            Kc11 = 0.9  ## with Cover crop min Kc=0.9
                        ##apply stress factor and start and end date to crop number 3

                        ##if CropType1 == 3:
                            ##if ii>KsStartDate and ii<ksEndDate:
                            ##    Kc11 = Kc11*ks
                    ##End of Kc11 adjustment for type 3 and 4
                    Kc11 = IKc1
                    if OKc1>Kc11:
                        Kc11=OKc1

                    ##May29-2007
                    if j==15:
                        Kc11=1.1                                                 
                    Kc[y,ii]=Kc11 
                    ETcDaily[y,ii]=EToDaily[y,ii]*Kc[y,ii]
                   
                ##end of for i<=dpyAll
                    
                ##....... Determine application number and amounts
                CETc = 0
                for ii in range(BeginDate1, EndDate1+1):
                    jj=ii 
                    if ii>dpyAll[y-1]:
                        jj=ii-dpyAll[y] 
                    ##modified for test on oct19-2005
                    ##CETc=CETc+ETodaily[y,j]*Kc[y,j] 
                    ETcDaily[y,jj]= EToDaily[y,jj]*Kc[y,jj] 
                    CETc=CETc+ETcDaily[y,jj] 
                    
                ET1 = 0
                NumI1 = 0
                ##if (PIrr.upperCase()=="Y") NumI1=1 
                ## modified on jan 8
                ##NumI=floor((B-A)/F) 
                NumI=int((B-BeginDate1)/f1) 
                NumI1=NumI1+NumI 
                BI=BeginDate1+f1*NumI 

                if NumI1==0:
                    BI=B
                BIYear[y]=BI
                for ii in range(BeginDate1,BI+1):
                    jj=ii 
                    if ii>dpyAll[y-1]:
                        jj=ii-dpyAll[y-1]
                    ##ET1=ET1+EToDaily[y,j]*Kc[y,j] 
                    ET1=ET1+ETcDaily[y,jj] 

                if NumI1!=0:
                    NA1[y]=ET1/NumI1 
                ET2=0 
                for ii in range(BI+1,EndDate1+1):
                    jj = ii
                    if ii>dpyAll[y-1]:
                        jj = ii-dpyAll[y-1]
                    
                    ET2 = ET2 + ETcDaily[y,jj]
                    
                        
                    if (ET2+ET1)>(CETc-YTD):
                        LI = ii
                        LIYear[y] = LI
                        NI2 = int(ET2/YTD)+1
                        NA2[y] = ET2/NI2
                        if NumI1 == 0:
                            NA1[y] = NA2[y]
                        NA3[y] = YTD
                        break
                        ##ii = EndDate1 + 1
                      
                        
                        ##to get out of loop i<=E
                ##end of for i<=E
            ##end of y<YCI
                
                
            for y in range(1,iyears+1):
                isCERn[y]=0 
                isPCP[y]=0 
                isETaw[y]=0 
                osCERn[y]=0 
                isCSpg[y]=0 
                osCSpg[y]=0 
                isCETc[y]=0 
                osCETc[y]=0 

            for Mon in range(0,12):
                isMCERn[Mon]=0 
                osMCERn[Mon]=0 
                isMCSpg[Mon]=0 
                osMCSpg[Mon]=0 
                isMCETc[Mon]=0 
                osMCETc[Mon]=0 
                EToMonthly[Mon]=0 
                ETcMonthly[Mon]=0 
                PcpMonthly[Mon]=0 
                ERnMonthly[Mon]=0 
                SpgMonthly[Mon]=0 
                EspgMonthly[Mon]=0 
                MonDsw[Mon]=0 
                MonDswPos[Mon]=0 
                MonNetApp[Mon]=0 
                
            ## Loop to calculate Kc's,Etc, &SWD for sceduling                     
            PSW = 0.0
            SWD = 0.0
            SWD0= osSWDx
            Espg=0
            
            ## for HSA**** (not for OLDHSA****)
            IrrigYear=0
                        
            for y in range(1,iyears+1):
                
                if j==6 or j==14:
                    yearTypeCal=yearType[y]
                else:
                    yearTypeCal=yearType[y-1]

                if yearTypeCal.upper()=="C" or yearTypeCal.upper()=="D":
                    BeginDate1=CBeginDate[j] 
                    EndDate1=CEndDate[j] 
                    kkc1=Ckc1[j] 
                    kkc2=Ckc2[j] 
                    kkc3=Ckc3[j] 
                    AB1=CAB[j] 
                    AC1=CAC[j] 
                    AD1=CAD[j] 
                    ##end of critical years
                else:
                    BeginDate1=NCBeginDate[j] 
                    EndDate1=NCEndDate[j] 
                    kkc1=NCkc1[j] 
                    kkc2=NCkc2[j] 
                    kkc3=NCkc3[j] 
                    AB1=NCAB[j] 
                    AC1=NCAC[j] 
                    AD1=NCAD[j] 
                    ##end of noncritical years
                if j==7:
                    RiceIrrigEndDate=EndDate1-20
                ##Initial Kc values for dates B,C,D,E
                KcB = kkc1
                KcC = kkc2 
                KcD = kkc2 
                KcE = kkc3

                ##..........................................................
                ##.. Calculate day of year corresponding to A,B,C, and E ...
                ## .. Note that B, C, D, and E can be bigger than 365   ....
                ## .. Note that BB,CC,DD, and EE are always <=365        ...
                ##..........................................................
                if EndDate1<=BeginDate1:
                    EndDate1=365+EndDate1 
                lenDate=EndDate1-BeginDate1 
                B=int(0.01*AB1*lenDate+BeginDate1) 
                C=int(0.01*AC1*lenDate+BeginDate1)
                D=int(0.01*AD1*lenDate+BeginDate1) 

                BB=B 
                if B>365:
                    BB=B-365 
                CC=C 
                if C>365:
                    CC=C-365 
                DD=D 
                if D>365:
                    DD=D-365 
                EE=EndDate1 
                if EndDate1>365:
                    EE=EndDate1-365
                    

                for Mon in range(0,12):
                    EToMonthly[Mon]=0 
                    ETcMonthly[Mon]=0 
                    PcpMonthly[Mon]=0 
                    ERnMonthly[Mon]=0 
                    SpgMonthly[Mon]=0 
                    EspgMonthly[Mon]=0 
                    MonNetApp[Mon]=0 
                    MonDsw[Mon]=0 
                    MonDswPos[Mon]=0 
                dpy = dpyAll[y]
                FinalIrrig = 0.0 
                Mon = 0
                
                ##CPcp = 0.0
                ##CERn = 0.0
                ##CESpg = 0.0
                ##CDsw = 0.0
                ##4/20/09 yearTypeDaily.append(" first")
                ##4/20/09 yearTypeDaily[0] = "First"
                
                for ii in range(1,dpy+1):
                    
                    ##initialize cumulative variable for the first day
                    if (y==iyears) and ((y%4!=0 and ii>273) or (y%4==0 and ii>274)):
                        break
                    if (y==1 and ii==274) or (y==1 and ii == 1):
                    ##if y==1 and ii = 274:
                        SWD=0 
                        PSWD=0 
                        NetApp=0 
                        CPcp=0 
                        CERn=0 
                        CESpg=0 
                        CDsw=0 
                        DCETo=0 
                        DCETc=0 
                    if Mon<12:
                        ##Spg=0.15/NumDaysPerMon[Mon+1]*erd
                        if j == 13 or j==12:   ## for native and riparian vegetation
                            Spg=0.15/NumDaysPerMon[Mon+1]*erd   ##0.025 0.05 0.15
                        else:
                            Spg=0.025/NumDaysPerMon[Mon+1]*erd   ##0.025 0.15
                    if Region[k] == 1:
                        Spg = 0.0
                    PSWD=SWD                           
                    OKc1=OKc[y,ii] 
                    IKc1=IKc[y,ii] 
                    Kc11=Kc[y,ii] 
                    ETo=EToDaily[y,ii]
                    PCP=PcpDaily[y,ii] 
                    ##This will reduce the value for Spg when the seepage is greater than the Dsw
                    
                    if (SWD+ETcDaily[y,ii]-Spg)>=0:
                        Espg=Spg 
                    else:
                        Espg=SWD+ETcDaily[y,ii]
                        if Region[k] == 1:     
                            Espg = 0
                    ## adjust for the rice
                    if j==7 and IKc1!=0:
                        Espg=0 
                    ##water surface and RV
                    ##water surface is now cropNum 15
                    ## if((j==12 or j==15)) Espg=ETcDaily[y,i] 
                    if j==12:
                        Espg=ETcDaily[y,ii]
                        if Region[k] == 1:     
                            Espg = 0
                    if j==15:
                        Espg=0
                        
                    ##if y == 33  and ii>250 and j==1:
                    ##    print y,ii,SWD,ETcDaily[y,ii],Espg,PCP
                     
                    if (SWD+ETcDaily[y,ii]-Espg-PCP)>=0:
                        ERn=PCP 
                    else:
                        ERn=SWD+ ETcDaily[y,ii]-Espg 

                    ##rice
                    if j==7 and IKc1!=0:
                        ERn=0 
                    ##if((j==12 or j==15)) ERn=0 
                    if j==12:
                        ERn=0 
                    if j==15:
                        ERn=PCP 
                    ## This calculates the Dsw adjusted for Espg and ERn
                    ##only for water surface dsw=+espg
                        
                    ##if (y == 34 or y ==35) and ii<120 and j==1:
                    ##    print y,ii,ETcDaily[y,ii], ERn, Espg, ETcDaily[y,ii]-ERn+Espg
                    
                        
                    if j!=15:
                        Dsw=ETcDaily[y,ii]-ERn-Espg 
                    else:          ##for water surface
                        Dsw=ETcDaily[y,ii]-ERn+Espg 
                    ## Dswp=Dsw 
                    
                    if ii>=BeginDate1:
                        SWDx=NA1[y] 
                        
                    if ii>=BIYear[y]:
                        SWDx=NA2[y] 
                    if ii>=LIYear[y]:
                        SWDx=NA3[y] 
                    if ii==EE:
                        SWDx=NA3[y] 
                    if EE<BeginDate1 and ii<BIYear[y]:
                        SWDx=NA2[y] 
                    ##if j== 1 and k == 0 and (y==1 or y==2):
                    ##    print NA1[y],NA2[y],NA3[y],BeginDate1,BIYear[y],LIYear[y],EE
                        
                    ##rice-water surface-Riparian
                    if j==7 or j==12 or j==15:
                        SWDx=osSWDx 
                    ##off season
                    if IKc1==0:
                        SWDx=SWD0
                        
                        
                    
                    if IKc1!=0:
                        if (SWD+Dsw)>SWDx:
                            NetApp=SWD+Dsw 
                        else:
                            NetApp=0 
                        ##rice - every day we have irrigation except the last 20 days
                        if j==7:
                            NetApp=SWD+Dsw 
                            if ii>RiceIrrigEndDate:
                                NetApp=0
                        ##***********add for Native Vegetation, not in the DETAW-UCD***************
                        if j==13:
                            NetApp = 0
                        ##*************************************************************
                    ## end of in-season
                    if IKc1 == 0:
                        NetApp = 0
                    SWD = SWD + Dsw - NetApp
                    if SWD<0:
                        SWD = 0 

                    if dpyAll[y] == 366:
                        if ii > NII[Mon]:
                            Mon = Mon + 1
                    else:
                        if ii>NI[Mon]:
                            Mon = Mon + 1
                                
                    if IKc1 == 0:
                        if SWD0<osSWDx:
                            SWD0=osSWDx
                        Diff = 0
                        if (SWD+Dsw) >=SWD0:
                            Diff=SWD0-SWD
                        if (SWD+Dsw) <SWD0:
                            if (y%4!=0 and ii>273) or (y%4==0 and ii>274):
                                osCETc[y]=osCETc[y]+ETcDaily[y,ii]*HAcre[k,y,j]*0.0081071
                            else:
                                osCETc[y-1]=osCETc[y-1]+ETcDaily[y,ii]*HAcre[k,y-1,j]*0.0081071
                            osMCETc[Mon]=osMCETc[Mon]+ETcDaily[y,ii] 
                        else:
                            if((y%4!=0 and ii>273)or (y%4==0 and ii>274)):
                                osCETc[y]=osCETc[y]+Diff*HAcre[k,y,j]*0.0081071 
                            else:
                                osCETc[y-1]=osCETc[y-1]+Diff*HAcre[k,y-1,j]*0.0081071 

                            osMCETc[Mon]=osMCETc[Mon]+Diff
                    else:
                        if (y%4!=0 and ii>273) or (y%4==0 and ii>274):
                            isCETc[y]=isCETc[y]+ETcDaily[y,ii]*HAcre[k,y,j]*0.0081071 
                        else:
                            isCETc[y-1]=isCETc[y-1]+ETcDaily[y,ii]*HAcre[k,y-1,j]*0.0081071 

                        isMCETc[Mon]=isMCETc[Mon]+ETcDaily[y,ii] 
                            
                    ##ISCERn and OsCERn are in and off-season cum effect rainfall
                    if IKc1 == 0:
                        if (y%4!=0 and ii>273) or (y%4==0 and ii>274):
                            osCERn[y]=osCERn[y]+ERn*HAcre[k,y,j]*0.0081071 
                            osCSpg[y]=osCSpg[y]+Spg*HAcre[k,y,j]*0.0081071
                        else:
                            osCERn[y-1]=osCERn[y-1]+ERn*HAcre[k,y-1,j]*0.0081071 
                            osCSpg[y-1]=osCSpg[y-1]+Spg*HAcre[k,y-1,j]*0.0081071 
                        osMCERn[Mon]=osMCERn[Mon]+ERn 
                        osMCSpg[Mon]=osMCSpg[Mon]+Spg
                    else:
                        if (y%4!=0 and ii>273) or (y%4==0 and ii>274):
                            isCERn[y]=isCERn[y]+ERn*HAcre[k,y,j]*0.0081071 
                            isPCP[y]=isPCP[y]+PCP*HAcre[k,y,j]*0.0081071 
                            isETaw[y]=isETaw[y]+NetApp*HAcre[k,y,j]*0.0081071 

                            isCSpg[y]=isCSpg[y]+Spg*HAcre[k,y,j]*0.0081071 
                        else:
                            isCERn[y-1]=isCERn[y-1]+ERn*HAcre[k,y-1,j]*0.0081071 
                            isPCP[y-1]=isPCP[y-1]+PCP*HAcre[k,y-1,j]*0.0081071 
                            isETaw[y-1]=isETaw[y-1]+NetApp*HAcre[k,y-1,j]*0.0081071 

                            isCSpg[y-1]=isCSpg[y-1]+Spg*HAcre[k,y-1,j]*0.0081071 
      
                        isMCERn[Mon]=isMCERn[Mon]+ERn  
                        isMCSpg[Mon]=isMCSpg[Mon]+Spg
                        ##end of Pcp
                        
                    ##the following 3 lines: no use now written in original code 
                    ##if IKc1 != 0 and flag == "s":  ##off season, the following 3 lines are modififed                         
                    ##    SWD0=SWD
                    ##flag = " "
                        
                    MaxSWD=erd*aw1*ADep1/100 
                    FC=MaxSWD*4 
                    PWP=MaxSWD*2 
                    SWC=FC-SWD 
                    ##YTDD is for YTD in column
                    YTDD=FC-SWDx    
                    ##calculate yeartypeCal for oct to oct year
                    if (y%4!=0 and ii>=274) or (y%4==0 and ii>=275):
                        yearTypeCal= yearType[y]
                    else:
                        yearTypeCal= yearType[y-1]
                    ##for first year only print ftom oct 1(day274)
                    CPcp=CPcp+PCP 
                    CERn=CERn+ERn 
                    CESpg=CESpg+Espg
                    ##for water surface we use espg* hacre each crop
                    if j!=15:
                        if(y%4!=0 and ii>273) or (y%4==0 and ii>274):
                            WSCESpg[y,ii]=WSCESpg[y,ii]+CESpg*HAcre[k,y,j]*0.0081071 
                            WSEspgMonthly[y,Mon]=WSEspgMonthly[y,Mon]+Espg*HAcre[k,y,j]*0.0081071 
                        else:
                            WSCESpg[y,ii]=WSCESpg[y,ii]+CESpg*HAcre[k,y-1,j]*0.0081071 
                            WSEspgMonthly[y,Mon]=WSEspgMonthly[y,Mon]+Espg*HAcre[k,y-1,j]*0.0081071 
                        ##before oct 1
                    ##end of if crop is not water surface
                    DCETo=DCETo+ EToDaily[y,ii] 
                    DCETc=DCETc+ETcDaily[y,ii] 
                    CDsw=CDsw+Dsw 
                
                    EToMonthly[Mon]=EToMonthly[Mon]+EToDaily[y,ii] 
                    ETcMonthly[Mon]=ETcMonthly[Mon]+ ETcDaily[y,ii] 
                    ERnMonthly[Mon]=ERnMonthly[Mon]+ERn 
                    SpgMonthly[Mon]=SpgMonthly[Mon]+Spg 
                    EspgMonthly[Mon]=EspgMonthly[Mon]+Espg 
                    
                    PcpMonthly[Mon]=PcpMonthly[Mon]+PcpDaily[y,ii] 
                    MonNetApp[Mon]=MonNetApp[Mon]+NetApp 
                    MonDsw[Mon]= ETcMonthly[Mon]- (ERnMonthly[Mon]+EspgMonthly[Mon]) 
                    ##if(MonACETAW[Mon]<0) MonACETAW[Mon]=0 
                    MonDswPos[Mon]=MonDsw[Mon] 
                    if MonDsw[Mon]<0:
                        MonDswPos[Mon]=0
                    
                    if(y!=1 and y!=iyears) or (y==1 and ii>273) or (y==iyears and ii<274 and y%4!=0) or (y==iyears and ii<275 and y%4==0):
                        
                        if j!=15:
                            if j == 14:
                                NetApp = 0
                            if (y%4!=0 and ii>273) or (y%4==0 and ii>274):
                                HAcre_temp = HAcre[k,y,j]*2.471
                              
                            else:
                                HAcre_temp = HAcre[k,y-1,j]*2.471
                            ##crop is not water surface
                        else:
                            Spg = 0
                            Espg = 0
                            NetApp = 0
                            if (y%4!=0 and ii>273) or (y%4==0 and ii>274):
                                HAcre_temp = HAcre[k,y,j]*2.471
                            else:
                                HAcre_temp = HAcre[k,y-1,j]*2.471
                        ##crop is water surface
                        if y==1 and ii==274:
                            ik = 1                        

                        ytemp[ik] = y+1920
                        DOYtemp[ik] = ii
                        if DOYtemp[ik] == 1:
                            IrrigYear = IrrigYear +1
                            
                        HAcretemp[ik] = HAcre_temp               
                        OKctemp[ik] = OKc[y,ii]
                        IKctemp[ik] = IKc[y,ii]
                        CCKctemp[ik] = CCKc[ii]
                        ETotemp[ik] = EToDaily[y,ii]
                        Kctemp[ik] = Kc[y,ii]
                        ETctemp[ik] = ETcDaily[y,ii]
                        Pcptemp[ik] = PcpDaily[y,ii]
                        
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
                            if j==6:
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

                        if y > 1 and ik==1:
                            for ic in range(1,dpy+1):
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
   
                        ##if y > 1 and ik==1:
                        if (y > 1 and ii == dpy) or (y == iyears and ik==1):
                            for iq in range(1,dpy+1):  ##add 04/29/09    
                                ##initialize cumulative for start of water year

                                ##@if iq ==1:
                                ##@    SumDelSWC = 0           
                                if yDaily[iq] == 1921 and iq ==1:
                                    FCtemp = int(FCDaily[iq]*10)/10.0
                                    ##PSWC = FCDaily[iq]
                                    PSWC = FCtemp
                                    SumDelSWC = 0
                                SWCtemp = int(SWCDaily[iq]*1000)/1000.0
                                ##DelSWC = PSWC - SWCDaily[iq]
                                DelSWC = PSWC - SWCtemp
                                SumDelSWC=SumDelSWC+DelSWC
                                PSWC=SWCDaily[iq]
                                ##do the calculation on the sept 30 of each year  and print
                                yy=yDaily[iq]
                                
                                dpy=365
                                if yy%4 == 0:
                                    dpy = 366
                    
                                if (yy%4!=0 and iq==365) or (yy%4==0 and iq==366):
                                    MonthlySWC = SumDelSWC/12.0
                                    ##convert etwawdaily to mon day
                                    etMon=10
                                    etDay=1
                                    for etdd in range(1, dpy+1):
                                        if etMon<10:
                                            yearCal=yy
                                        else:
                                            yearCal=yy-1
                                        if yearCal%4==0:
                                            etDOY=etdd+274
                                            if etDOY>366:
                                                etDOY=etDOY-366
                                        if yearCal%4!=0:
                                            etDOY=etdd+273
                                            if etDOY>365:
                                                etDOY=etDOY-365
                                    
                                        if yearCal%4==0:
                                            if etDOY>NII[etMon-1]:
                                                etMon = etMon + 1
                                                etDay=1
                                        else:
                                            if etDOY>NI[etMon-1]:
                                                etMon = etMon + 1
                                                etDay=1
                        
                                        ##ETAW=DswDaily[dd]-MonthlySWC/NumDay[Mon]
                                        Dswtemp = int(DswDaily[etdd]*1000)/1000.0
                                        ETAWMonDay[etMon,etDay]=Dswtemp-MonthlySWC/NumDay[etMon]                                        
                                        ##@ETAWMonDay[etMon,etDay]=Dswtemp-abs(MonthlySWC/NumDay[etMon])
                                        
                                        etDay=etDay+1
                                        if (yearCal%4!=0 and etDOY==365) or (yearCal%4==0 and etDOY==366):
                                            etMon=1
                                            etDay=1       
                                    ##end of for etdd

                                    ##calculate aveage ETaw
                                    SumNegETAW = 0.0
                                    CountNegETAW = 0
                                    for etMon in range(1,13):
                                        flagETAW = "true"
                                        while flagETAW == "true":
                                            for etDay in range(1,NumDay[etMon]+1):
                                                if ETAWMonDay[etMon,etDay] < 0:
                                                    SumNegETAW=SumNegETAW+ETAWMonDay[etMon,etDay]
                                                    CountNegETAW=CountNegETAW+1
                                                    ETAWMonDay[etMon,etDay]=0
                                            if abs(SumNegETAW) > 0:
                                                flagETAW = "true"
                                            else:
                                                flagETAW = "false"
                                            if NumDay[etMon] > CountNegETAW:
                                                avgNeg = SumNegETAW/(NumDay[etMon]-CountNegETAW)
                                            SumNegETAW = 0.0
                                            CountNegETAW = 0
                                            for etDay in range(1,NumDay[etMon]+1):
                                                if ETAWMonDay[etMon,etDay] > 0:
                                                    ##@ETAWMonDay[etMon,etDay] = ETAWMonDay[etMon,etDay]
                                                    ETAWMonDay[etMon,etDay] = ETAWMonDay[etMon,etDay]+avgNeg
                                                    
                                                    
                                    ##convert etawdmon day to etaw daily
                                    etMon = 10
                                    etDay = 1
                                    for etdd in range(1,dpy+1):
                                        if etMon < 10:
                                            yearCal = yy
                                        else:
                                            yearCal = yy -1
                                        if yearCal%4==0:
                                            etDOY = etdd + 274
                                            if etDOY > 366:
                                                etDOY = etDOY - 366
                                        if yearCal%4!=0:
                                            etDOY = etdd + 273
                                            if etDOY > 365:
                                                etDOY = etDOY - 365
                                        if yearCal%4==0:
                                            if etDOY>NII[etMon-1]:
                                                etMon = etMon + 1
                                                etDay = 1
                                        else:
                                            if etDOY > NI[etMon-1]:
                                                etMon = etMon + 1
                                                etDay = 1
                                        ETAWDaily[etdd] = ETAWMonDay[etMon,etDay]
                                        
                            
                                        etDay = etDay + 1
                                        if (yearCal%4!=0 and etDOY==365) or (yearCal%4==0 and etDOY==366):
                                            etMon = 1
                                            etDay = 1
                                    ##ii = 0

                                                                                    
                                    CETAWDaily= 0
                                    Mon = 10
                                    for dd in range(1,dpy+1):
                            
                                        if Mon<10:
                                            yearCal = yy
                                        else:
                                            yearCal = yy-1
                                        if yearCal%4==0:
                                            if DOY[dd]>NII[Mon-1]:
                                                Mon = Mon+1
                                                ##MonETAW[yearCal-1920,Mon]= 0.0
                                        else:
                                            if DOY[dd]>NI[Mon-1]:
                                                Mon = Mon + 1
                                                ##MonETAW[yearCal-1920,Mon]= 0.0
                                        if (yearCal%4!=0 and DOY[dd]==365) or (yearCal%4==0 and DOY[dd]==366):
                                            Mon = 1

                                                                                       
                                        if(j!=6):
                                            if DOY[dd] >= DOYLIrrig[yDaily[dd]-1921]:
                                                ETAWDaily[dd] = 0
                                        else:
                                            if DOYLIrrig[yDaily[dd]-1921] < 275:
                                                if DOY[dd]>=DOYGrainLIrrig[yDaily[dd]-1921] and DOY[dd]<275:
                                                    ETAWDaily[dd] = 0
                                            else:
                                                if DOY[dd]>=DOYLIrrig[yDaily[dd]-1921]:
                                                    ETAWDaily[dd] = 0
                                        
                                        CETAWDaily = CETAWDaily + ETAWDaily[dd]
                                        ## y is for sep 30 each year, so for previous oct-Dec we subtract y by 1
                                        if Mon<10:                                        
                                            MonETAW[yy-1920,Mon] = MonETAW[yy-1920,Mon] + ETAWDaily[dd]                                          
                                       
                                        else:                                            
                                            MonETAW[yy-1921,Mon] = MonETAW[yy-1921,Mon]+ ETAWDaily[dd]
                                          
                                        temp = HAcreDaily[dd]*0.00328084

                                        if idayoutput==1:
                                            if dailyunit == 1:
                                                dataday.append(ETcDaily2[dd]*temp)
                                                data2day.append(PcpDaily2[dd]*temp)
                                                data6day.append(EToDaily2[dd]*temp)                           
                                                data7day.append(KcDaily[dd]*temp)                                            
                                                data8day.append(SWCDaily[dd]*temp)
                                                data9day.append(SpgDaily[dd]*temp)
                                                data10day.append(ESpgDaily[dd]*temp)
                                                data11day.append(DswDaily[dd]*temp)
                                                data12day.append(ETAWDaily[dd]*temp)
                                                data13day.append(ErDaily[dd]*temp)
                                                data14day.append(SWDxDaily[dd]*temp)
                                                data15day.append(FCDaily[dd]*temp)
                                                data16day.append(PWPDaily[dd]*temp)
                                                data17day.append(SWDDaily[dd]*temp)
                                                data18day.append(YTDDaily[dd]*temp)
                                            else: 
                                                dataday.append(ETcDaily2[dd])
                                                data2day.append(PcpDaily2[dd])
                                                data6day.append(EToDaily2[dd])                           
                                                data7day.append(KcDaily[dd])                                            
                                                data8day.append(SWCDaily[dd])
                                                data9day.append(SpgDaily[dd])
                                                data10day.append(ESpgDaily[dd])
                                                data11day.append(DswDaily[dd])
                                                data12day.append(ETAWDaily[dd])
                                                data13day.append(ErDaily[dd])
                                                data14day.append(SWDxDaily[dd])
                                                data15day.append(FCDaily[dd])
                                                data16day.append(PWPDaily[dd])
                                                data17day.append(SWDDaily[dd])
                                                data18day.append(YTDDaily[dd])
                                            ##data19day.append(NADaily[dd])
                                            ##data20day.append(CPcpDaily[dd])
                                            ##data21day.append(CErDaily[dd])
                                            ##data22day.append(CESpgDaily[dd])
                                            ##data23day.append(CETcDaily[dd])
                                            ##data24day.append(CDswDaily[dd])
                                            ##data25day.append(CETAWDaily)
                                            data26day.append(NADaily[dd]*temp)
                                            data27day.append(CPcpDaily[dd]*temp)
                                            data28day.append(CErDaily[dd]*temp)
                                            data29day.append(CESpgDaily[dd]*temp)
                                            data30day.append(CETcDaily[dd]*temp)
                                            data31day.append(CDswDaily[dd]*temp)
                                            data32day.append(CETAWDaily*temp)
                            
                                        SumDelSWC = 0
                    ##print "   dpy=", dpy, "  y=",y, "  ii=",ii
                        if (y%4!=0 and ii==273) or (y%4==0 and ii==274):
                            NetApp=0 
                            CPcp=0 
                            CERn=0 
                            CESpg=0 
                            DCETc=0 
                            CDsw=0 
                             
                        
                        if ii == dpy or (y==iyears and ii==273 and y%4!=0) or (y==iyears and ii==274 and y%4==0):
                            for Mon in range(0,12):
                                if(y!=1 and y!=iyears) or (y==1 and Mon>=9) or (y==iyears and Mon<9):                                                                    
                                    
##____________________________OldHSA [
                                    if j!=15:
                                        ##NetApp = 0 for non-irrig Grain
                                        if j==14:
                                            MonNetApp[Mon] = 0
                                        if Mon >= 9:
##____________________________OldHSA ]
                                            temp = HAcre[k,y,j]*0.0081071
                                        else:                             
                                            temp = HAcre[k,y-1,j]*0.0081071
                                    else:
                                        MonNetApp[Mon] = 0
                                        SpgMonthly[Mon]=0 
                                        EspgMonthly[Mon]=0
                                        if Mon>=9:
                                            temp = HAcre[k,y,j]*0.0081071
                                        else:
                                            temp = HAcre[k,y-1,j]*0.0081071
                                    if imonthoutput == 1:        
                                        datamon.append(ETcMonthly[Mon])
                                        data2mon.append(EToMonthly[Mon])
                                        data3mon.append(MonNetApp[Mon])
                                        data4mon.append(PcpMonthly[Mon])
                                        data5mon.append(ERnMonthly[Mon])
                                        data6mon.append(SpgMonthly[Mon])
                                        data7mon.append(EspgMonthly[Mon])
                                        data8mon.append(MonDsw[Mon])
                                        data9mon.append(MonDswPos[Mon])
                                        ##data10mon.append(MonETAW[y,Mon])
                                        ##data11mon.append(MonETAWPos)
                                        data12mon.append(MonNetApp[Mon]*temp)
                                        data13mon.append(PcpMonthly[Mon]*temp)
                                        data14mon.append(ERnMonthly[Mon]*temp)
                                        if j!=15:
                                            data15mon.append(EspgMonthly[Mon]*temp)
                                        else:
                                            data15mon.append(WSEspgMonthly[y,Mon])
                                        data16mon.append(ETcMonthly[Mon]*temp)
                                        data17mon.append(MonDsw[Mon]*temp)
                                        data18mon.append(MonDswPos[Mon]*temp)
                                        ##data19mon.append(MonETAW[y,Mon]*temp)
                                        ##data20mon.append(MonETAWPos*temp)
                                        data21mon.append(EToMonthly[Mon]*temp)
                                        
                                    if y > 1:
                                        if Mon == 0:
                                            for imtemp in range(10,13):
                                                data10mon.append(MonETAW[y-1,imtemp])
                                                data19mon.append(MonETAW[y-1,imtemp]*temp)
                                                MonETAWPos=MonETAW[y-1,imtemp]    ##5/1/09 revised
                                                if MonETAWPos< 0:
                                                    MonETAWPos=0
                                                if imonthoutput == 1:
                                                    data11mon.append(MonETAWPos)
                                                    data20mon.append(MonETAWPos*temp)
                                        if Mon < 9:        
                                            MonETAWPos=MonETAW[y,Mon+1]    ##5/1/09 revised
                                            if MonETAWPos< 0:
                                                MonETAWPos=0
                                            if imonthoutput == 1:
                                                if j == 15 and k==0 and MonETAW[y,Mon+1]>0.0:
                                                    print(y,Mon+1, " MonETAW=",MonETAW[y,Mon+1])
                                                data10mon.append(MonETAW[y,Mon+1])
                                                data19mon.append(MonETAW[y,Mon+1]*temp)
                                                data11mon.append(MonETAWPos)                                  
                                                data20mon.append(MonETAWPos*temp)                                                
                        
                        
########Combine OLDHSA and HSA together#################################################                        
                            
            
            
            SACETC=0 
            SACERN=0 
            SACSpg=0 
            SISCETC=0 
            SISCERN=0 
            SISCSpg=0 
            SOSCETC=0 
            SOSCERN=0 
            SOSCSpg=0
            for y in range(1,iyears+1):
                SISCETC=SISCETC+isCETc[y] 
                SISCERN=SISCERN+isCERn[y] 
                SISCSpg=SISCSpg+isCSpg[y] 
                SOSCETC=SOSCETC+osCETc[y] 
                SOSCERN=SOSCERN+osCERn[y] 
                SOSCSpg=SOSCSpg+osCSpg[y]

                ##**save: y+1921,yearType[y],HAcre[k,y,j]*2.471,isPCP[y],isCETc[y],isCERn[y],isCSpg[y],
                ##**isETaw[y],osCETc[y],osCERn[y],osCSpg[y]
                if iyearoutput == 1:
                    datayr.append(isCETc[y])
                    data2yr.append(isPCP[y])
                    data3yr.append(isCERn[y])
                    data4yr.append(isCSpg[y])
                    data5yr.append(isETaw[y])
                    data6yr.append(osCETc[y])
                    data7yr.append(osCERn[y])
                    data8yr.append(osCSpg[y])
            ## calculate &print mean over years for CETc,CERn, ETAW
            MACETC=SACETC/iyears 
            MACERN=SACERN/iyears 
            MACSpg=SACSpg/iyears 
            MISCETC=SISCETC/iyears 
            MISCERN=SISCERN/iyears 
            MISCSpg=SISCSpg/iyears 
            MOSCETC=SOSCETC/iyears 
            MOSCERN=SOSCERN/iyears 
            MOSCSpg=SOSCSpg/iyears
            
            
            #start = start1+days(2)
            start = str(start1[2]+1)+monthname[start1[1]-1]+str(start1[0])+" "+"2300"
            dt = "1DAY"  #days(1)
            unit = ""
            if idayoutput==1:
                #cpartt = "ETc"
                #props = {"unit":"mm"}
                DETAWOUTPUT[0,k,j-1,:] = dataday[:] 
                DETAWOUTPUT[1,k,j-1,:] = data10day[:]
                DETAWOUTPUT[2,k,j-1,:] = data2day[:]
                DETAWOUTPUT[3,k,j-1,:] = data12day[:]
                DETAWOUTPUT[4,k,j-1,:] = data11day[:]
                DETAWOUTPUT[5,k,j-1,:] = data13day[:]
                if forDSM2_daily == 1:
                    #ddatalist = dataday,data2day,data10day,data11day,data12day,data13day
                    ddatalist = dataday,data10day,data2day,data12day,data11day,data13day
                    #dcpartlist = "ETc","PCP","ESpg","Dsw","ETAW","ER"
                    dcpartlist = "ETc","ESpg","PCP","ETAW","Dsw","ER"
                    dunitlist = "A-ft","A-ft","A-ft","A-ft","A-ft","A-ft"   
                else:
                    ddatalist = dataday,data2day,data6day,data7day,data8day,data9day,  \
                           data10day,data11day,data12day,data13day,    \
                           data14day,data15day,data16day,data17day,data18day,  \
                           data26day,data27day,data28day,data29day,data30day,data31day,data32day
                    dcpartlist = "ETc","PCP","ETo","Kc","SWC","Spg",   \
                            "ESpg","Dsw","ETAW","ER",    \
                            "SWDxDaily","FCDaily","PWPDaily","SWDDaily","YTDDaily", \
                            "NA","CPcp","CEr","CEspg","CETc","CDsw","CETaw"
                    if dailyunit == 1:
                        dunitlist = "A-ft","A-ft","A-ft","A-ft","A-ft","A-ft",   \
                            "A-ft","A-ft","A-ft","A-ft",    \
                            "A-ft","A-ft","A-ft","A-ft","A-ft",    \
                            "A-ft","A-ft","A-ft","A-ft","A-ft","A-ft","A-ft"
                    else:
                        dunitlist = "mm","mm","mm","mm","mm","mm",   \
                            "mm","mm","mm","mm",    \
                            "mm","mm","mm","mm","mm",    \
                            "A-ft","A-ft","A-ft","A-ft","A-ft","A-ft","A-ft"
                if itotaloutput == 1:
                    if j==1 and k==0:
                        data1day_14crops = dataday
                        data10day_14crops = data10day
                        data11day_14crops = data11day
                        data12day_14crops = data12day
                        data13day_14crops = data13day
                    elif j<15:
                        data1day_14crops = add(data1day_14crops, dataday)
                        data10day_14crops = add(data10day_14crops, data10day)
                        data11day_14crops = add(data11day_14crops, data11day)
                        data12day_14crops = add(data12day_14crops, data12day)
                        data13day_14crops = add(data13day_14crops, data13day)
                    if j==15:
                        if k==0:
                            data1day_water = dataday
                            data10day_water = data10day
                            data11day_water = data11day
                            data12day_water = data12day
                            data13day_water = data13day
                        else:
                            data1day_water = add(data1day_water, dataday)
                            data10day_water = add(data10day_water, data10day)
                            data11day_water = add(data11day_water, data11day)
                            data12day_water = add(data12day_water, data12day)
                            data13day_water = add(data13day_water, data13day)
            
                        
            ##start = start1            
            dt2 = "1MONTH"  #months(1)            
            if imonthoutput == 1:
                #props = {"unit":"mm"}                
                mdatalist= datamon,data2mon,data3mon,data4mon,data5mon,data6mon, \
                           data7mon,data8mon,data9mon,data10mon,data11mon, \
                           data12mon,data13mon,data14mon,data15mon,data16mon,  \
                           data17mon,data18mon,data19mon,data20mon,data21mon
                mcpartlist = "ETc","ETo","NetApp","Pcp","ERn","Spg",  \
                             "Espg","Dsw","Dsw+ve","ETAW","ETAW+ve",  \
                             "NA_A_ft","Pcp_A_ft","Er_A_ft","Espg_A_ft","ETc_A_ft",  \
                             "Dsw_A_ft","Dsw+ve_A_ft","ETaw_A_ft","ETaw+ve_A_ft","ETo_A_ft"
                munitlist = "mm","mm","mm","mm","mm","mm",  \
                             "mm","mm","mm","mm","mm",  \
                             "A-ft","A-ft","A-ft","A-ft","A-ft",  \
                             "A-ft","A-ft","A-ft","A-ft","A-ft"
                
                if itotaloutput != 1:
                    if j==2:
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
                    if j>2 and j<12:
                        data12mon_veg = add(data12mon_veg,data12mon)
                        data13mon_veg = add(data13mon_veg,data13mon)
                        data14mon_veg = add(data14mon_veg,data14mon)
                        data15mon_veg = add(data15mon_veg,data15mon)
                        data16mon_veg = add(data16mon_veg,data16mon)
                        data17mon_veg = add(data17mon_veg,data17mon)
                        data18mon_veg = add(data18mon_veg,data18mon)
                        data19mon_veg = add(data19mon_veg,data19mon)
                        data20mon_veg = add(data20mon_veg,data20mon)
                        data21mon_veg = add(data21mon_veg,data21mon)
                if itotaloutput == 1:
                    if j==1 and k==0:
                        data14mon_total = data14mon
                        data15mon_total = data15mon
                        data16mon_total = data16mon
                        data20mon_total = data20mon
                    elif j<15:
                        data14mon_total = add(data14mon_total,data14mon)
                        data15mon_total = add(data15mon_total,data15mon)
                        data16mon_total = add(data16mon_total,data16mon)
                        data20mon_total = add(data20mon_total,data20mon)
                    if j==15:
                        if k==0:
                            data14mon_water = data14mon
                            data15mon_water = data15mon
                            data16mon_water = data16mon
                            data20mon_water = data20mon
                        else:
                            data14mon_water = add(data14mon_water,data14mon)
                            data15mon_water = add(data15mon_water,data15mon)
                            data16mon_water = add(data16mon_water,data16mon)
                            data20mon_water = add(data20mon_water,data20mon)
            
            ##start = start1
            dt3 = "1YEAR"  #years(1)
            if iyearoutput == 1:
                ydatalist = datayr, data2yr, data3yr, data4yr,data5yr,data6yr,data7yr,data8yr
                ycpartlist="CETc","CPcp","CEr","CSpg","CETaw","OCETc","OCEr","OCSpg"
                
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
            
            if idayoutput == 1:
                destination = filepath+"\\Output\\DETAW_day_"+str(ktemp)+".dss"
                dssfh=pyhecdss.DSSFile(destination)
                for ilist in range(0,len(ddatalist)):
                    path = "/"+Apart+"/"+Bpart+"/"+dcpartlist[ilist]+"//"+Epart+"/"+Fpart+"/"
                    listlen = len(ddatalist[ilist])
                    write_to_dss(dssfh, ddatalist[ilist], path,startdate +" "+starttime, dunitlist[ilist], ctype)
                    #temptext += write_1ts_to_txt(start,ddatalist[ilist],Apart,Bpart,dcpartlist[ilist], \
                    #           Epart,Fpart,dunitlist[ilist],dt)
                
                if itotaloutput == 1:
                    if j==15 and k == (ilands-1):
                        print("daily output")
                        path = "/"+Apart+"/"+Bpart+"/ETc//"+Epart+"/"+Fpart+"/"
                        listlen = len(data1day_14crops)
                        if dailyunit == 1:
                            dunits = "A-FT"
                        else:
                            dunits = "mm"
                        write_to_dss(dssfh, data1day_14crops, path,startdate +" "+starttime, dunits, ctype)
                        #temptext += write_1ts_to_txt(start,data1day_14crops,Apart,Bpart,"ETc", \
                        #       Epart,Fpart,"mm",dt)
                        path = "/"+Apart+"/"+Bpart+"/ESpg//"+Epart+"/"+Fpart+"/"
                        write_to_dss(dssfh, data10day_14crops, path,startdate +" "+starttime, dunits, ctype)
                        #temptext += write_1ts_to_txt(start,data10day_14crops,Apart,Bpart,"ESpg", \
                        #       Epart,Fpart,"mm",dt)
                        path = "/"+Apart+"/"+Bpart+"/Dsw//"+Epart+"/"+Fpart+"/"
                        write_to_dss(dssfh, data11day_14crops, path,startdate +" "+starttime, dunits, ctype)
                        #temptext += write_1ts_to_txt(start,data11day_14crops,Apart,Bpart,"Dsw", \
                        #       Epart,Fpart,"mm",dt)
                        path = "/"+Apart+"/"+Bpart+"/ETaw//"+Epart+"/"+Fpart+"/"
                        write_to_dss(dssfh, data12day_14crops, path,startdate +" "+starttime, dunits, ctype)
                        #temptext += write_1ts_to_txt(start,data12day_14crops,Apart,Bpart,"ETaw", \
                        #       Epart,Fpart,"mm",dt)
                        path = "/"+Apart+"/"+Bpart+"/Er//"+Epart+"/"+Fpart+"/"
                        write_to_dss(dssfh, data13day_14crops, path,startdate +" "+starttime, dunits, ctype)
                        #temptext += write_1ts_to_txt(start,data13day_14crops,Apart,Bpart,"Er", \
                        #       Epart,Fpart,"mm",dt)
                        path = "/"+Apart+"/"+Bpart+"/ETc//"+Epart+"/"+Fpart+"/"
                        write_to_dss(dssfh, data1day_water, path,startdate +" "+starttime, dunits, ctype)
                        #temptext += write_1ts_to_txt(start,data1day_water,Apart,Bpart,"ETc", \
                        #       Epart,Fpart,"mm",dt)
                        path = "/"+Apart+"/"+Bpart+"/ESpg//"+Epart+"/"+Fpart+"/"
                        write_to_dss(dssfh, data1day_water, path,startdate +" "+starttime, dunits, ctype)
                        #temptext += write_1ts_to_txt(start,data10day_water,Apart,Bpart,"ESpg", \
                        #       Epart,Fpart,"mm",dt)
                        path = "/"+Apart+"/"+Bpart+"/Dsw//"+Epart+"/"+Fpart+"/"
                        write_to_dss(dssfh, data11day_water, path,startdate +" "+starttime, dunits, ctype)
                        #temptext += write_1ts_to_txt(start,data11day_water,Apart,Bpart,"Dsw", \
                        #       Epart,Fpart,"mm",dt)
                        path = "/"+Apart+"/"+Bpart+"/ETaw//"+Epart+"/"+Fpart+"/"
                        write_to_dss(dssfh, data13day_14crops, path,startdate +" "+starttime, dunits, ctype)
                        #temptext += write_1ts_to_txt(start,data12day_water,Apart,Bpart,"ETaw", \
                        #       Epart,Fpart,"mm",dt)
                        path = "/"+Apart+"/"+Bpart+"/Er//"+Epart+"/"+Fpart+"/"
                        write_to_dss(dssfh, data13day_water, path,startdate +" "+starttime, dunits, ctype)
                        #temptext += write_1ts_to_txt(start,data13day_water,Apart,Bpart,"Er", \
                        #       Epart,Fpart,"mm",dt)       
                dssfh.close()
                
            Epart = "1MONTH"
            if imonthoutput == 1:
                destination = filepath+"\\Output\\DETAW_month.dss"
                dssfh=pyhecdss.DSSFile(destination)
                
                for ilist in range(0,21):
                    path = "/"+Apart+"/"+Bpart+"/"+mcpartlist[ilist]+"//"+Epart+"/"+Fpart+"/"
                    listlen = len(mdatalist[ilist])
                    templist = list(mdatalist[ilist])
                    tempunit = munitlist[ilist]
                    write_to_dss(dssfh, templist, path,startdate +" "+starttime, tempunit, ctype)
                    #temptext += write_1ts_to_txt(start,mdatalist[ilist],Apart,Bpart,mcpartlist[ilist], \
                    #           Epart,Fpart,munitlist[ilist],dt2)
                                            
                if itotaloutput == 1:
                    if j==15 and k == (ilands-1):
                        print("monthly output")
                        listlen = len(data14mon_total)
                        path = "/"+Apart+"/total_no_w/Er//"+Epart+"/"+Fpart+"/"
                        write_to_dss(dssfh, data14mon_total, path,startdate +" "+starttime, "A-ft", ctype)
                        #temptext += write_1ts_to_txt(start,data14mon_total,Apart,"total_no_w","Er", \
                        #       Epart,Fpart,"A-ft",dt2)
                        path = "/"+Apart+"/total_no_w/Espg//"+Epart+"/"+Fpart+"/"
                        write_to_dss(dssfh, data15mon_total, path,startdate +" "+starttime, "A-ft", ctype)
                        #temptext += write_1ts_to_txt(start,data15mon_total,Apart,"total_no_w","Espg", \
                        #       Epart,Fpart,"A-ft",dt2)
                        path = "/"+Apart+"/total_no_w/ETc//"+Epart+"/"+Fpart+"/"
                        write_to_dss(dssfh, data16mon_total, path,startdate +" "+starttime, "A-ft", ctype)
                        #temptext += write_1ts_to_txt(start,data16mon_total,Apart,"total_no_w","ETc", \
                        #       Epart,Fpart,"A-ft",dt2)
                        path = "/"+Apart+"/total_no_w/ETaw+ve//"+Epart+"/"+Fpart+"/"
                        write_to_dss(dssfh, data20mon_total, path,startdate +" "+starttime, "A-ft", ctype)
                        #temptext += write_1ts_to_txt(start,data20mon_total,Apart,"total_no_w","ETaw+ve", \
                        #       Epart,Fpart,"A-ft",dt2)
            
                        path = "/"+Apart+"/water/Er//"+Epart+"/"+Fpart+"/"
                        write_to_dss(dssfh, data14mon_water, path,startdate +" "+starttime, "A-ft", ctype)
                        #temptext += write_1ts_to_txt(start,data14mon_water,Apart,"water","Er", \
                        #       Epart,Fpart,"A-ft",dt2)
                        path = "/"+Apart+"/water/Espg//"+Epart+"/"+Fpart+"/"
                        write_to_dss(dssfh, data15mon_water, path,startdate +" "+starttime, "A-ft", ctype)
                        #temptext += write_1ts_to_txt(start,data15mon_water,Apart,"water","Espg", \
                        #       Epart,Fpart,"A-ft",dt2)
                        path = "/"+Apart+"/water/ETc//"+Epart+"/"+Fpart+"/"
                        write_to_dss(dssfh, data16mon_water, path,startdate +" "+starttime, "A-ft", ctype)
                        #temptext += write_1ts_to_txt(start,data16mon_water,Apart,"water","ETc", \
                        #      Epart,Fpart,"A-ft",dt2)
                        path = "/"+Apart+"/water/ETaw+ve//"+Epart+"/"+Fpart+"/"
                        write_to_dss(dssfh, data20mon_water, path,startdate +" "+starttime, "A-ft", ctype)
                        #temptext += write_1ts_to_txt(start,data20mon_water,Apart,"water","ETaw+ve", \
                        #       Epart,Fpart,"A-ft",dt2)
                
                else:
                    if j == 11:
                        listlen = len(data12mon_veg)
                        path = "/"+Apart+"/veg/NA//"+Epart+"/"+Fpart+"/"
                        write_to_dss(dssfh, data12mon_veg, path,startdate +" "+starttime, "A-ft", ctype)
                        #temptext += write_1ts_to_txt(start,data12mon_veg,Apart,"veg","NA", \
                        #       Epart,Fpart,"A-ft",dt2)
                        path = "/"+Apart+"/veg/Pcp//"+Epart+"/"+Fpart+"/"
                        write_to_dss(dssfh, data13mon_veg, path,startdate +" "+starttime, "A-ft", ctype)
                        #temptext += write_1ts_to_txt(start,data13mon_veg,Apart,"veg","Pcp", \
                        #       Epart,Fpart,"A-ft",dt2)
                        path = "/"+Apart+"/veg/Er//"+Epart+"/"+Fpart+"/"
                        write_to_dss(dssfh, data14mon_veg, path,startdate +" "+starttime, "A-ft", ctype)
                        #temptext += write_1ts_to_txt(start,data14mon_veg,Apart,"veg","Er", \
                        #       Epart,Fpart,"A-ft",dt2)
                        path = "/"+Apart+"/veg/Espg//"+Epart+"/"+Fpart+"/"
                        write_to_dss(dssfh, data15mon_veg, path,startdate +" "+starttime, "A-ft", ctype)
                        #temptext += write_1ts_to_txt(start,data15mon_veg,Apart,"veg","Espg", \
                        #       Epart,Fpart,"A-ft",dt2)
                        path = "/"+Apart+"/veg/ETc//"+Epart+"/"+Fpart+"/"
                        write_to_dss(dssfh, data16mon_veg, path,startdate +" "+starttime, "A-ft", ctype)
                        #temptext += write_1ts_to_txt(start,data16mon_veg,Apart,"veg","ETc", \
                        #       Epart,Fpart,"A-ft",dt2)
                        path = "/"+Apart+"/veg/Dsw//"+Epart+"/"+Fpart+"/"
                        write_to_dss(dssfh, data17mon_veg, path,startdate +" "+starttime, "A-ft", ctype)
                        #temptext += write_1ts_to_txt(start,data17mon_veg,Apart,"veg","Dsw", \
                        #       Epart,Fpart,"A-ft",dt2)
                        path = "/"+Apart+"/veg/Dsw+ve//"+Epart+"/"+Fpart+"/"
                        write_to_dss(dssfh, data18mon_veg, path,startdate +" "+starttime, "A-ft", ctype)
                        #temptext += write_1ts_to_txt(start,data18mon_veg,Apart,"veg","Dsw+ve", \
                        #       Epart,Fpart,"A-ft",dt2)
                        path = "/"+Apart+"/veg/ETaw//"+Epart+"/"+Fpart+"/"
                        write_to_dss(dssfh, data19mon_veg, path,startdate +" "+starttime, "A-ft", ctype)
                        #temptext += write_1ts_to_txt(start,data19mon_veg,Apart,"veg","ETaw", \
                        #       Epart,Fpart,"A-ft",dt2)
                        path = "/"+Apart+"/veg/ETaw+ve//"+Epart+"/"+Fpart+"/"
                        write_to_dss(dssfh, data20mon_veg, path,startdate +" "+starttime, "A-ft", ctype)
                        #temptext += write_1ts_to_txt(start,data20mon_veg,Apart,"veg","ETaw+ve", \
                        #       Epart,Fpart,"A-ft",dt2)
                        path = "/"+Apart+"/veg/ETo//"+Epart+"/"+Fpart+"/"
                        write_to_dss(dssfh, data21mon_veg, path,startdate +" "+starttime, "A-ft", ctype)
                        #temptext += write_1ts_to_txt(start,data21mon_veg,Apart,"veg","ETo", \
                        #       Epart,Fpart,"A-ft",dt2)
                dssfh.close()
                
                
            Epart = "1YEAR"
            if iyearoutput == 1:
                destination = filepath+"\\Output\\DETAW_year.dss"
                dssfh=pyhecdss.DSSFile(destination)
                for ilist in range(0,len(ydatalist)):
                    path = "/"+Apart+"/"+Bpart+"/"+ycpartlist[ilist]+"//"+Epart+"/"+Fpart+"/"
                    listlen = len(ydatalist[ilist])
                    write_to_dss(dssfh, ydatalist[ilist], path,startdate +" "+starttime, yunitlist[ilist], ctype)
                    #temptext += write_1ts_to_txt(start,ydatalist[ilist],Apart,Bpart,ycpartlist[ilist], \
                    #           Epart,Fpart,yunitlist[ilist],dt3)
                dssfh.close()
    return(DETAWOUTPUT)
##_______________________________________________________________________________

if __name__ == "__main__":

    files = listdir(".")
    filepath = os.getcwd()
    
    inputfile = "./Input/DETAW_para.inp"
    f0 = open(inputfile, 'r')
    itemp = 0
    endyear = 0
    for line in f0:
        if line:
            if not("#" in line):
                itemp = itemp + 1
                if "Model to streamline" in line:
                    modelno = int(line.split("=")[1])
                    if modelno == 1:
                        streamlinemodel = "DSM2"
                    elif modelno == 2:
                        streamlinemodel = "SCHISM"
                    elif modelno == 3:
                        streamlinemodel = "CALSIM3"
                if "Daily output" in line:
                    idayoutput = int(line.split("=")[1])
                if "Monthly output" in line:
                    imonthoutput = int(line.split("=")[1])
                if "Yearly output" in line:
                    iyearoutput = int(line.split("=")[1])
                if  "Delta output" in line:
                    itotaloutput = int(line.split("=")[1])
                if  "Daily output unit" in line:
                    dailyunit = int(line.split("=")[1])
                if "forDSM2_daily" in line:
                    forDSM2_daily = int(line.split("=")[1])
                if "End year" in line:
                    endyear = int(line.split("=")[1])
                if "Days" in line:
                    idates = int(line.split("=")[1])
    
    file = ["  "]*3
    
    file[0] = "mm_Pcp.csv"
    file[1] = "Percentage.csv"
    file[2] = "LODI_PT.csv"

    fileout = "weather.dss"
    

    ts_type = "rts"
    start1 = [1921,9,30,23,0]
    iyears = endyear-start1[0]+1
    #iyears = 2015-start1[0]+1 #2011-start1[0]+1  ##2008  1924
    
    ilands = 168
    #idates = 34219 #32873  ##31870  1096
    isites = 7
    NumDay=[0,31,28,31,30,31,30,31,31,30,31,30,31]
    NumDayL=[0,31,29,31,30,31,30,31,31,30,31,30,31]
    NI=[31,59,90,120,151,181,212,243,273,304,334,365]
    NII=[31,60,91,121,152,182,213,244,274,305,335,366]
    pcplocs = ["Brentwood","Davis","Galt","Lodi","RioVista","Stockton","Tracy"]
    perclocs = ["Davis","Stockton","Lodi","Tracy Carbona","Rio Vista","Brentwood","Galt"]
    ts_pcp = zeros((isites,idates),float)
    ts_per = zeros((isites,ilands),float)
    ETo_corrector = zeros((ilands),float)
    Region = zeros((ilands),int)
    ts_year = zeros((idates),int)
    ts_mon = zeros((idates),int)
    ts_days = zeros((idates),int)
    ts_LODI_tx = zeros((idates),float)
    ts_LODI_tn = zeros((idates),float)
    for ifile in range(0,3):
        print(file[ifile])
        if streamlinemodel == "CALSIM3":
            source = filepath+"\\Input\\planning_study\\"+file[ifile]
        else:
            source = filepath+"\\Input\\historical_study\\"+file[ifile]
        ff = open(source,"r")
        if ifile == 0:
            icon = 0
            for line in ff:
                if icon > 0:
                    for j in range(0,isites):
                        ts_pcp[j,icon-1] = float(line.split(",")[4+j])
                icon += 1            
        if ifile == 1:
            icon = 0
            for line in ff:
                if icon > 0:
                    for j in range(0,isites):
                        for k in range(0,isites):
                            if perclocs[j][0:2] in pcplocs[k]:
                                ts_per[k,icon-1] = float(line.split(",")[j+1])
                    ETo_corrector[icon-1] = float(line.split(",")[8])
                    Region[icon-1] = int(line.split(",")[9])
                icon += 1            
        if ifile == 2:
            icon = 0
            for line in ff:
                if icon > 0:
                    ts_year[icon-1] = int(line.split(",")[1])
                    ts_mon[icon-1] = int(line.split(",")[2])
                    ts_days[icon-1] = int(line.split(",")[3])
                    ts_LODI_tx[icon-1] = float(line.split(",")[5])
                    ts_LODI_tn[icon-1] = float(line.split(",")[6])
                icon += 1
        ff.close()
    
    (pcp,ET0) = weatheroutput(ts_pcp,ts_per,ts_mon,ts_days,ts_LODI_tx,ts_LODI_tn,ilands,idates,isites,ETo_corrector,filepath,start1)
    
    (DETAWOUTPUT) = historicalETAW(ts_per,ETo_corrector,Region,pcp,ET0,ts_LODI_tx,ts_LODI_tn, \
            ilands,idates,isites,ts_year,ts_mon,ts_days,start1,filepath,NI,NII,NumDay,iyears,  \
            idayoutput,imonthoutput,iyearoutput,itotaloutput,dailyunit,forDSM2_daily,streamlinemodel)
    
    (DETAWISL168) = timeseries_combine(DETAWOUTPUT,ilands,ilands,15,idates-1,"")
    forNODCU(DETAWISL168,streamlinemodel,endyear,ilands,"")
    if streamlinemodel == "CALSIM3":
        print("in the double-counting process", idates)
        tempfile = filepath+"\\Input\\planning_study\\"+"CS3_DCD_rate1.txt"
        (DETAWISL168) = timeseries_combine(DETAWOUTPUT,ilands,ilands,15,idates-1,tempfile)
        forNODCU(DETAWISL168,streamlinemodel,endyear,ilands,"_ex1")
        
        tempfile = filepath+"\\Input\\planning_study\\"+"CS3_DCD_rate2.txt"
        (DETAWISL168) = timeseries_combine(DETAWOUTPUT,ilands,ilands,15,idates-1,tempfile)
        forNODCU(DETAWISL168,streamlinemodel,endyear,ilands,"_ex2")
        
        tempfile = filepath+"\\Input\\planning_study\\"+"CS3_DCD_rate3.txt"
        (DETAWISL168) = timeseries_combine(DETAWOUTPUT,ilands,ilands,15,idates-1,tempfile)
        forNODCU(DETAWISL168,streamlinemodel,endyear,ilands,"_ex3")
        
    print("done")