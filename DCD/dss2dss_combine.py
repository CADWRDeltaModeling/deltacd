# After the ETAW for water has been calculated, it has to be removed in this file. 12/22/22015

import pandas as pd
import pyhecdss
from datetime import datetime         
from numpy import arange,sin,pi
import os, sys, string
from os import listdir


## giving a existing datasource.
##source = open("SA0011.xls", "r")
def get_pathname(dssfh, path):
    '''
    reads catalog to get full pathname if it exists
    The assumption is that the path provided does not have a time window
    returns the first matching pathname with exactly the A,B,C,E and F parts
    '''
    pathparts=path.split('/')
    dfcat=dssfh.read_catalog()
    dfpath=dfcat[(dfcat.A==pathparts[1]) & (dfcat.B==pathparts[2]) 
                    & (dfcat.C==pathparts[3]) & (dfcat.E==pathparts[5])
                    & (dfcat.F==pathparts[6])]
    pathname=dssfh.get_pathnames(dfpath)[0]
    return pathname



def main(argv):
    Cpart = ["  "]*6
    Cpart[0] = "ETC"  #"ETC_A_FT"
    Cpart[1] = "ESPG"  #"ESPG_A_FT"
    Cpart[2] = "PCP"  #"PCP_A_FT"
    Cpart[3] = "ETAW"   #"ETAW+VE_A_FT"
    Cpart[4] = "DSW"   #"DSW+VE_A_FT"
    Cpart[5] = "ER"  #"ER_A_FT"
    
    pyhecdss.set_message_level(0)
    pyhecdss.set_program_name('DCD')
    outputfile = ".\hist_isl_crop.dss"
    dssofh=pyhecdss.DSSFile(outputfile)
    fileno = int(sys.argv[1])
    print(fileno)
    owd = os.getcwd()
    
    for ifile in range(fileno-1,fileno):
        inputfile = "./DETAW_day_"+str(ifile)+".dss"
        dssfh=pyhecdss.DSSFile(inputfile)      
        dfcat=dssfh.read_catalog()
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
                
                if iterm != 4:
                    if iterm == 3:
                        tempcrop = 14
                    else:
                        tempcrop = 15
                    for icrop in range(0,tempcrop):
                        Fpart = "CROP_"+str(icrop+1)
                        path = "/DWRBDO-DETAW/"+Bpart+"/"+Cpart[iterm]+"//1DAY/"+Fpart+"/"
                        tss,cunits,ctype = dssfh.read_rts(get_pathname(dssfh,path))
                        print(iland, " crop=", icrop)
                        if icrop == 0:
                            ttss = tss
                        else:
                            ttss.iloc[:,0] = ttss.iloc[:,0] + tss.iloc[:,0]
                        del tss
                else:
                    Fpart = "CROP_15"
                    path1 = "/DWRBDO-DETAW/"+Bpart+"/"+Cpart[0]+"//1DAY/"+Fpart+"/"
                    path2 = "/DWRBDO-DETAW/"+Bpart+"/"+Cpart[2]+"//1DAY/"+Fpart+"/"
                    ttss,cunits,ctype = dssfh.read_rts(get_pathname(dssfh,path1))
                    ttss2,cunits2,ctype2 = dssfh.read_rts(get_pathname(dssfh,path2)) 
                    ttss.iloc[:,0]=ttss.iloc[:,0]-ttss2.iloc[:,0]
                #print Bpart
                if iterm == 3 or iterm==4:
                    ttss.clip(lower=0.0,inplace=True)
                if (iland+1)<10:
                    Bpart = "HSA000" + str(iland+1)
                elif (iland+1) < 100:
                    Bpart = "HSA00" + str(iland+1)
                else:
                    Bpart = "HSA0" + str(iland+1)
                path = "/DETAW/"+Bpart+"/"+Cpart[iterm]+"//1DAY/DWR/"
                #if iterm == 5:
                #    path = "/DETAW/"+Bpart+"/W_"+Cpart[iterm]+"//1DAY/DWR/"
                dssofh.write_rts(path,ttss,cunits,ctype)
                del ttss             
                
                
if __name__ == "__main__":
    main(sys.argv[1:])