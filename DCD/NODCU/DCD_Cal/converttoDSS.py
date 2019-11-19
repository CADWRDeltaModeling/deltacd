import pyhecdss
import pandas as pd
from numpy import arange,pi,array,zeros

import os, sys, string, math, numpy
from os import listdir
from math import cos,sin,tan,atan,sqrt,pi,pow


def write_to_dss(dssfh, arr, path, startdatetime, cunits, ctype):
    '''
    write to the pyhecdss.DSSFile for an array with starttime and assuming
    daily data with the pathname path, cunits and ctype
    '''
    if "2400" in startdatetime:
        startdatetime=startdatetime.replace("2400","2300")
        
    fstr='1D'
    epart=path.split('/')[5]
    if epart == '1DAY':
        fstr='1D'
    elif epart == '1MONTH':
        fstr='1M'
    else:
        raise RuntimeError('Not recognized frequency in path: %s'%path)
    df=pd.DataFrame(arr,index=pd.date_range(startdatetime,periods=len(arr),freq=fstr))
    dssfh.write_rts(path,df.shift(freq='H'),cunits,'PER-AVER')
    
if __name__ == "__main__":    
    pyhecdss.set_message_level(0)
    
    f1 = open(sys.argv[1])
    
    ileng = 0
    tempseries = []
    oneseries = 0
    for line in f1:
        if ileng == 0 and ".dss" in line:
            outputfile = line.strip()
            dssfh=pyhecdss.DSSFile(outputfile)
        else:
            if "A=" in line and oneseries == 0:
                pathnames = line.split("=")
                Apart = pathnames[1].split("  ")[0].strip()
                Bpart = pathnames[2][:5].strip()
                Cpart = pathnames[3].split("  ")[0].strip()
                Epart = pathnames[5].split("  ")[0].strip()
                Fpart = pathnames[6].strip()
                path = "/"+Apart+"/"+Bpart+"/"+Cpart+"//"+Epart+"/"+Fpart+"/"
                oneseries = 1 
            elif oneseries == 1:
                cunit = line.strip()
                oneseries += 1
            elif oneseries == 2:
                ctype = line.strip()
                oneseries += 1
            elif oneseries == 3:
                starttime = line.strip()
                oneseries += 1
            elif oneseries > 3:
                if "END" in line:
                    write_to_dss(dssfh, tempseries, path,starttime, cunit, ctype)
                    oneseries = 0
                    tempseries = []
                elif "FINISH" in line:
                    print("Done")         
                else:
                    tempseries.append(float(line.strip()))
                    oneseries += 1
    f1.close()
    
    