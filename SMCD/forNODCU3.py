# the runoff (PCP) has been distributed evenly into the following 5 days.  _7/22/2013
# Espg has been double counted. It has to be reduced. 8/7/2012
# extra rainfall is spreaded as an assumed hydrograph, 1 day 0.32, 2 day 0.38, 3 day, 0.2, 4 day 0.1


# Import VTools time and numpy array creation function
from vtools.data.vtime import *
from vtools.data.constants import *
from vtools.data.timeseries import *   
from datetime import datetime         
from numpy import arange,sin,pi,array, zeros

# Import VTools dss utility function
from vtools.datastore.dss.api import *
# Import VTools Excel utility function
from vtools.datastore.excel.api import *
import os, sys, string
from os import listdir
import calendar
 
## giving a existing datasource.
##source = open("SA0011.xls", "r")

def daysofayr(n):
    if calendar.isleap(n):  #leap year
        days = 366
    else:     #not leap year
        days = 365
    return days

def daysinyr_currentday(n,m,d):
    leapdays = [31,29,31,30,31,30,31,31,30,31,30,31]
    noleapdays = [31,28,31,30,31,30,31,31,30,31,30,31]
    days = 0
    for i in range(0,12):
        if m == (i+1):
            if i == 0:
                days = d
            else:
                for j in range(0,i):
                    if calendar.isleap(n):
                        days += leapdays[j]
                    else:
                        days += noleapdays[j]
                days += d
    return days

def main(argv1,argv2):
    Cpart = ["  "]*5
    Cpart[0] = "ETAW"
    Cpart[1] = "ESPG"
    Cpart[2] = "ETC"
    Cpart[3] = "PCP"
    Cpart[4] = "DSW"
    
    beginyear = 1922
    endyear = int(argv1)    
    tyr = endyear-beginyear+1
    num_islands = int(argv2)
    #tempts = zeros((tyr,366),float)
    
    inputfile = "detaw_168to142.dss"
    f1 = open("DICU5.17","w")
    f2 = open("DICU5.14","w")
    f3 = open("DICU5.27","w")
    f4 = open("DICU5.12","w")
    f5 = open("DICU5.30","w")
    
    daysofyear = zeros((366,4,tyr+1),int)
    leapdays = [31,29,31,30,31,30,31,31,30,31,30,31]
    noleapdays = [31,28,31,30,31,30,31,31,30,31,30,31]
    kyr = 0
    for iyr in range(beginyear-1,endyear+1):
        if iyr == (beginyear-1):
            for m in range(10,13):
                for d in range(1,noleapdays[m-1]+1):
                    if m==10 and d==1:
                        idays = 0
                    daysofyear[idays,0,kyr]=iyr
                    daysofyear[idays,1,kyr]=m
                    daysofyear[idays,2,kyr]=d
                    daysofyear[idays,3,kyr]=daysinyr_currentday(iyr,m,d)
                    idays += 1        
        elif iyr == endyear:
            for m in range(1,10):
                if calendar.isleap(iyr):
                    for d in range(1,leapdays[m-1]+1):
                        if m==10 and d==1:
                            idays = 0
                            kyr += 1
                        daysofyear[idays,0,kyr]=iyr
                        daysofyear[idays,1,kyr]=m
                        daysofyear[idays,2,kyr]=d
                        daysofyear[idays,3,kyr]=daysinyr_currentday(iyr,m,d)
                        idays += 1
                else:
                    for d in range(1,noleapdays[m-1]+1):
                        if m==10 and d==1:
                            if iyr == endyear:
                                break
                            idays = 0
                            kyr += 1    
                        #@print idays, kyr,iyr,m,d,daysinyr_currentday(iyr,m,d)
                        daysofyear[idays,0,kyr]=iyr
                        daysofyear[idays,1,kyr]=m
                        daysofyear[idays,2,kyr]=d
                        daysofyear[idays,3,kyr]=daysinyr_currentday(iyr,m,d)
                        idays += 1
        else:
            for m in range(1,13):
                if calendar.isleap(iyr):
                    for d in range(1,leapdays[m-1]+1):
                        if m==10 and d==1:
                            idays = 0
                            kyr += 1
                        daysofyear[idays,0,kyr]=iyr
                        daysofyear[idays,1,kyr]=m
                        daysofyear[idays,2,kyr]=d
                        daysofyear[idays,3,kyr]=daysinyr_currentday(iyr,m,d)
                        idays += 1
                else:
                    for d in range(1,noleapdays[m-1]+1):
                        if m==10 and d==1:
                            if iyr == endyear:
                                break
                            idays = 0
                            kyr += 1    
                        #print idays, kyr,iyr,m,d,daysinyr_currentday(iyr,m,d)
                        daysofyear[idays,0,kyr]=iyr
                        daysofyear[idays,1,kyr]=m
                        daysofyear[idays,2,kyr]=d
                        daysofyear[idays,3,kyr]=daysinyr_currentday(iyr,m,d)
                        idays += 1
                        
    for ifile in range(0,5): #:(3,4)
        print "which file =", ifile
        for iland in range(1,num_islands+1):
            print  "island ",iland
            if ifile == 0:
                strt = str("%3i" % iland)+"AT1 4 HISTORIC DEPLETION OF APPLIED WATER BY IRR. AND URBAN, AREA   "
                strt = strt + str(iland)+ "\n"
                f1.writelines(strt)
                strt = str("%3i" % iland)+"AT2            FROM DETAW STUDY \n"
                f1.writelines(strt)
            elif ifile == 1:
                strt = str("%3i" % iland)+"AT1 4  TOTAL HISTORIC CU OF SEEPAGE  "
                strt = strt + str(iland)+ "\n"
                f2.writelines(strt)
                strt = str("%3i" % iland)+"AT2            FROM DETAW STUDY \n"
                f2.writelines(strt)
            elif ifile == 2:
                strt = str("%3i" % iland)+"AT1 4   HISTORIC DEPLETION, AREA   "
                strt = strt + str(iland)+ "\n"
                f3.writelines(strt)
                strt = str("%3i" % iland)+"AT2            FROM DETAW STUDY \n"
                f3.writelines(strt)
            elif ifile == 3:
                strt = str("%3i" % iland)+"AT1 4 TOTAL BASIN RUNOFF, AREA   "
                strt = strt + str(iland)+ "\n"
                f4.writelines(strt)
                strt = str("%3i" % iland)+"AT2            FROM DETAW STUDY \n"
                f4.writelines(strt)
            else:
                strt = str("%3i" % iland)+"AT1 4 TOTAL NET DEPLETION of WATERBODY, AREA   "
                strt = strt + str(iland)+ "\n"
                f5.writelines(strt)
                strt = str("%3i" % iland)+"AT2            FROM DETAW STUDY \n"
                f5.writelines(strt)    
                
            if iland < 10:
                Bpart = "HSA000" + str(iland)
            elif iland < 100:
                Bpart = "HSA00" + str(iland)
            else:
                Bpart = "HSA0" + str(iland)
            path = "/DETAW/"+Bpart+"/"+Cpart[ifile]+"//1DAY/DWR/"
            tss = dss_retrieve_ts(inputfile,path)
            ndays = 0
            tempts = zeros((tyr,366),float)
            for iyr in range(0,tyr):
                for iday in range(0,366):
                    if (iday == 365 and daysofyear[iday,0,iyr] <> 0) or (iday<365):
                        if ifile == 3:
                            if iyr == 0 and iday < 5:
                                for kday in range(0, iday+1):
                                    temp_c = 0.2
                                    tempts[iyr][iday] += tss.data[iday-kday]*temp_c
                            else:
                                for kday in range(0,5):
                                    temp_c = 0.2
                                    tempts[iyr][iday] += tss.data[ndays-kday]*temp_c                            
                        else:
                            #print iyr, iday, ndays, ifile
                            #print tss.data[ndays]
                            tempts[iyr][iday] = tss.data[ndays]
                        # Water doesn't has daily ESPG. Monthly data double counts ESPG because of water body.
                        #if ifile == 1:   
                        #    tempts[iyr][iday] = tempts[iyr][iday]/2.0
                        ndays+=1
            
            for iyr in range(0,tyr):
                strt = str("%3i" % iland)+"A   "+ str(beginyear+iyr)
                for iday in range(0,366):
                    if (iday == 365 and daysofyear[iday,0,iyr] <> 0) or (iday<365):
                        strt = strt + str("%8.1f" % tempts[iyr][iday])
                strt = strt + "\n"
                if ifile == 0:
                    f1.writelines(strt)
                elif ifile == 1:
                    f2.writelines(strt)
                elif ifile == 2:
                    f3.writelines(strt)
                elif ifile == 3:
                    f4.writelines(strt)
                else:
                    f5.writelines(strt)
    f1.close()
    f2.close()
    f3.close()
    f4.close()
    f5.close()
    
if __name__ == "__main__":
    main(sys.argv[1],sys.argv[2])