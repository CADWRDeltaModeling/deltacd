# After the ETAW for water has been calculated, it has to be removed in this file. 12/22/22015

# Import VTools time and numpy array creation function
from vtools.data.vtime import *
from vtools.data.constants import *
from vtools.data.timeseries import *   
from datetime import datetime         
from numpy import arange,sin,pi
from vtools.functions.shift import *

# Import VTools dss utility function
from vtools.datastore.dss.api import *
# Import VTools Excel utility function
from vtools.datastore.excel.api import *
import os, sys, string
from os import listdir


## giving a existing datasource.
##source = open("SA0011.xls", "r")



def main(argv):
    Cpart = ["  "]*6
    Cpart[0] = "ETC"  #"ETC_A_FT"
    Cpart[1] = "ESPG"  #"ESPG_A_FT"
    Cpart[2] = "PCP"  #"PCP_A_FT"
    Cpart[3] = "ETAW"   #"ETAW+VE_A_FT"
    Cpart[4] = "DSW"   #"DSW+VE_A_FT"
    Cpart[5] = "ER"  #"ER_A_FT"
    
    
    outputfile = ".\hist_isl_crop.dss"
    fileno = int(sys.argv[1])
    print fileno
    owd = os.getcwd()
    
    for ifile in range(fileno-1,fileno):
        inputfile = "./DETAW_day_"+str(ifile)+".dss"        
        for iterm in range(0,6):  
            isl_beg = 10*ifile
            isl_end = 10*(ifile+1)
            if ifile == 16:
                isl_end = 168
            for iland in range(isl_beg,isl_end):     #168):
                if (iland+1)<10:
                    Bpart = "HSA_000" + str(iland+1)
                elif (iland+1) < 100:
                    Bpart = "HSA_00" + str(iland+1)
                else:
                    Bpart = "HSA_0" + str(iland+1)
                
                if iterm <> 4:
                    if iterm == 3:
                        tempcrop = 14
                    else:
                        tempcrop = 15
                    for icrop in range(0,tempcrop):
                        Fpart = "CROP_"+str(icrop+1)
                        path = "/DWRBDO-DETAW/"+Bpart+"/"+Cpart[iterm]+"//1DAY/"+Fpart+"/"
                        tss = dss_retrieve_ts(inputfile,path)
                        print iland, " crop=", icrop
                        
                        if icrop == 0:
                            ttss = tss
                        else:
                            ttss = ttss + tss
                        del tss
                else:
                    Fpart = "CROP_15"
                    path1 = "/DWRBDO-DETAW/"+Bpart+"/"+Cpart[0]+"//1DAY/"+Fpart+"/"
                    path2 = "/DWRBDO-DETAW/"+Bpart+"/"+Cpart[2]+"//1DAY/"+Fpart+"/"
                    ttss = dss_retrieve_ts(inputfile,path1) - dss_retrieve_ts(inputfile,path2)                   
                
                #print Bpart
                if iterm == 3 or iterm==4:
                    for i in range(0, len(ttss)):
                        if ttss.data[i]<0.0:
                            ttss.data[i] = 0.0
                
                if (iland+1)<10:
                    Bpart = "HSA000" + str(iland+1)
                elif (iland+1) < 100:
                    Bpart = "HSA00" + str(iland+1)
                else:
                    Bpart = "HSA0" + str(iland+1)
                path = "/DETAW/"+Bpart+"/"+Cpart[iterm]+"//1DAY/DWR/"
                #if iterm == 5:
                #    path = "/DETAW/"+Bpart+"/W_"+Cpart[iterm]+"//1DAY/DWR/"
                ttss.props[AGGREGATION] = "MEAN"
                ttss=shift(ttss, days(1))
                dss_store_ts(ttss, outputfile, path)
                del ttss             
                
                
if __name__ == "__main__":
    main(sys.argv[1:])