#Edit the effective precipitation for islands > 143, found the bug 1/16/2015

# Import VTools time and numpy array creation function
from vtools.data.vtime import *
from vtools.data.constants import *
from vtools.data.timeseries import *   
from datetime import datetime         
from numpy import arange,sin,pi,array, zeros
from vtools.functions.shift import *

# Import VTools dss utility function
from vtools.datastore.dss.api import *
# Import VTools Excel utility function
from vtools.datastore.excel.api import *
import os, sys, string
from os import listdir
import time

## giving a existing datasource.
##source = open("SA0011.xls", "r")

if __name__ == "__main__":

    Cpart = ["  "]*6
    Cpart[0] = "ETC"
    Cpart[1] = "ESPG"
    Cpart[2] = "PCP"
    Cpart[3] = "DSW"
    Cpart[4] = "ETAW"
    Cpart[5] = "ER"  #"ER_A_FT"
  
    inputfile =".\hist_isl_crop.dss"
    outputfile = ".\detaw_168to142.dss"
    
    islid = zeros((26,2),int)
    
    
    for iterm in range(0,5): #5):
    ##for iterm in range(2,2):
        ##destination="R:\DETAW 8-8-07\Run_Jan3\HistoricalSAOutput\hist_mon_isl.dss"
        ##destination="Z:\lliang\detaw\croptoisland\hist_mon.dss"

        for iland in range(0,15):
            print iland
            #time.sleep(2)
            if (iland+1) < 143:
                if (iland+1)<10:
                    Bpart = "HSA000" + str(iland+1)
                elif (iland+1) < 100:
                    Bpart = "HSA00" + str(iland+1)
                else:
                    Bpart = "HSA0" + str(iland+1)
                    
                path = "/DETAW/"+Bpart+"/"+Cpart[iterm]+"//1DAY/DWR/"
                if iterm == 2:
                    path2 = "/DETAW/"+Bpart+"/"+Cpart[5]+"//1DAY/DWR/"
                    tss = dss_retrieve_ts(inputfile,path) - dss_retrieve_ts(inputfile,path2)
                    for ig in range(0,len(tss)):
                        if tss.data[ig] < 0:
                            tss.data[ig] = 0.0
                    #dss_store_ts(tss, outputfile2, path)
                    tss.props[AGGREGATION] = "MEAN"
                    tss=shift(tss, days(1))  
                else:
                    tss = dss_retrieve_ts(inputfile,path)
                    
            else:
                Bpart = "HSA0" + str(iland+1)
                path = "/DETAW/"+Bpart+"/"+Cpart[iterm]+"//1DAY/DWR/"
                if iterm == 2:
                    path2 =  "/DETAW/"+Bpart+"/"+Cpart[5]+"//1DAY/DWR/"
                    tss1 = dss_retrieve_ts(inputfile,path) - dss_retrieve_ts(inputfile,path2)
                    for ig in range(0,len(tss1)):
                        if tss1.data[ig] < 0:
                            tss1.data[ig] = 0.0
                    #tss1.props[AGGREGATION] = "MEAN"
                    #tss1=shift(tss1, days(1)) 
                else:
                    tss1 = dss_retrieve_ts(inputfile,path)
                
                for id in range(0, isl):
                    if islid[id,0] == (iland+1):
                        if islid[id,1]<10:
                            Bpart = "HSA000" + str(islid[id,1])
                        elif islid[id,1]< 100:
                            Bpart = "HSA00" + str(islid[id,1])
                        else:
                            Bpart = "HSA0" + str(islid[id,1])
                path = "/DETAW/"+Bpart+"/"+Cpart[iterm]+"//1DAY/DWR/"
                
                tss = dss_retrieve_ts(outputfile,path)
                tss = tss + tss1
                tss.props[AGGREGATION] = "MEAN"
                tss=shift(tss, days(1))
                
                ##path = "/DETAW/"+Bpart+"/"+Cpart[iterm]+"//1DAY/MTV/"
            
            dss_store_ts(tss, outputfile, path)