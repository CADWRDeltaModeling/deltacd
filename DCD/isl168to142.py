#Edit the effective precipitation for islands > 143, found the bug 1/16/2015

import pandas as pd
import pyhecdss
from datetime import datetime         
from numpy import arange,sin,pi,array, zeros
import os, sys, string
from os import listdir
import time
from dss2dss_combine import get_pathname
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
    pyhecdss.set_message_level(0)
    pyhecdss.set_program_name('DCD')

    inputfile =".\hist_isl_crop.dss"
    dssifh=pyhecdss.DSSFile(inputfile)
    outputfile = ".\detaw_168to142.dss"
    dssofh=pyhecdss.DSSFile(outputfile)
    islandfile = ".\island_id.txt"
    f1 = open(islandfile)
    islid = zeros((26,2),int)
    

    isl = 0
    for line in f1:
        if line:
            ll = line.split()
            islid[isl,0] = int(ll[0])
            islid[isl,1] = int(ll[1])
            isl = isl + 1
    f1.close()
    
    for iterm in range(0,5): #5):
    ##for iterm in range(2,2):
        ##destination="R:\DETAW 8-8-07\Run_Jan3\HistoricalSAOutput\hist_mon_isl.dss"
        ##destination="Z:\lliang\detaw\croptoisland\hist_mon.dss"

        for iland in range(0,168):
            print(iland)
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
                    tss,cunits,ctype = dssifh.read_rts(get_pathname(dssifh,path))
                    tssx,cux,ctx = dssifh.read_rts(get_pathname(dssifh,path2))
                    tss.iloc[:,0]=tss.iloc[:,0]-tssx.iloc[:,0]
                    tss.clip(lower=0.0,inplace=True)
                    #dss_store_ts(tss, outputfile2, path)
                else:
                    tss,cunits,ctype = dssifh.read_rts(get_pathname(dssifh,path))   
            else:
                Bpart = "HSA0" + str(iland+1)
                path = "/DETAW/"+Bpart+"/"+Cpart[iterm]+"//1DAY/DWR/"
                if iterm == 2:
                    path2 =  "/DETAW/"+Bpart+"/"+Cpart[5]+"//1DAY/DWR/"
                    tss1, cunits, ctype = dssifh.read_rts(get_pathname(dssifh,path)) 
                    tss1x, cx, ct = dssifh.read_rts(get_pathname(dssifh,path2))
                    tss1.iloc[:,0] = tss1.iloc[:,0] - tss1x.iloc[:,0]
                    tss1.clip(lower=0.0, inplace=True)
                else:
                    tss1, cunits, ctype = dssifh.read_rts(get_pathname(dssifh,path))
                
                for id in range(0, isl):
                    if islid[id,0] == (iland+1):
                        if islid[id,1]<10:
                            Bpart = "HSA000" + str(islid[id,1])
                        elif islid[id,1]< 100:
                            Bpart = "HSA00" + str(islid[id,1])
                        else:
                            Bpart = "HSA0" + str(islid[id,1])
                path = "/DETAW/"+Bpart+"/"+Cpart[iterm]+"//1DAY/DWR/"
                
                tss, cunits, ctype = dssofh.read_rts(get_pathname(dssofh,path))
                tss.iloc[:,0] = tss.iloc[:,0] + tss1.iloc[:,0]
            
            dssofh.write_rts(path, tss, cunits, ctype)