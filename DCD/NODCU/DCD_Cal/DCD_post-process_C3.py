
# Import VTools time and numpy array creation function
from datetime import datetime,timedelta  
import pyhecdss
import pandas as pd
import numpy as np
import pyhecdss
import os, sys, string
import calendar

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
    
def PP_DSS():
    inputfile =".\RO_island.dss"
    outputfile = ".\DP_island.dss"
    
    pyhecdss.set_message_level(0)
    dssfh=pyhecdss.DSSFile(inputfile)
    dssofh=pyhecdss.DSSFile(outputfile)
    plist=dssfh.get_pathnames()
    for p in plist:
        df,u,p=dssfh.read_rts(p)
        df.values[:]=df.values[:]*0.25
        dssofh.write_rts(df.columns[0].replace("RO-FLOW","DP-FLOW"),df.shift(freq='D'),u,p)  #'PER-AVER')
    dssfh.close()
    dssofh.close()
    print("Deep percolation")    
    
def write_to_dss(dssfh, arr, path, startdatetime, cunits, ctype):
    '''
    write to the pyhecdss.DSSFile for an array with starttime and assuming
    daily data with the pathname path, cunits and ctype
    '''
    fstr='1D'
    epart=path.split('/')[5]
    if epart == '1DAY':
        fstr='1D'
    elif epart == '1MON':
        fstr='1M'
    else:
        raise RuntimeError('Not recognized frequency in path: %s'%path)
    #df=pd.DataFrame(arr,index=pd.date_range(startdatetime,periods=len(arr),freq=fstr))
    df=pd.DataFrame(arr,index=pd.period_range(startdatetime,periods=len(arr),freq=fstr))
    dssfh.write_rts(path,df,cunits,ctype)
    
def daytomonth(inputfile):
    d=pyhecdss.DSSFile(inputfile)
    outputfile = os.getcwd()+"\\"+inputfile.split(".")[0]+"_mon.dss"
    do=pyhecdss.DSSFile(outputfile)
    plist=d.get_pathnames()
    for p in plist:
        df,u,p=d.read_rts(p)
        do.write_rts(df.columns[0].replace('1DAY','1MON'),df.resample('M').mean(),u,'PER-AVER')
    d.close()
    do.close()    

def DCD_to_CALSIM(inputfile):
    inputfile = inputfile.split(".")[0]+"_mon.dss"
    dssifh=pyhecdss.DSSFile(inputfile)
    outputfile = inputfile.split(".")[0]+"_C3.dss"
    dssofh=pyhecdss.DSSFile(outputfile)
    DCD_C3_nodes = os.path.join("..","DCD_CALSIM3_nodes.csv")
    C3_nodes = ["OMR","SJR_EAST","SJR_WEST","SAC_WEST","MOK","SAC_SOUTH","SAC_NORTH","50_PA2"]
    C3_paths = ["IRR","SEEP","DRN"]
    DCD_paths = ["DIV-FLOW","SEEP-FLOW","DRAIN-FLOW"]
    
    f0 = open(DCD_C3_nodes)
    DCD_nodes = []
    ili = 0
    for line in f0:
        ili += 1
        if line:
            if ili > 1:
                DCD_nodes.append(line)
    f0.close()
    
    for ipath in range(0,len(C3_paths)):
        for c3j in range(0,len(C3_nodes)):
            inode = 0
            for i in range(0,len(DCD_nodes)):
                if C3_nodes[c3j] == DCD_nodes[i].split(",")[1].strip():
                    inode += 1
                    if DCD_nodes[i].split(",")[0].strip()=="BBID":
                        path = "/DICU-HIST+RSVR/"+DCD_nodes[i].split(",")[0].strip()+"/"+DCD_paths[ipath]+"//1MON/DWR-BDO/"
                    else:
                        path = "/DICU-HIST+NODE/"+DCD_nodes[i].split(",")[0].strip()+"/"+DCD_paths[ipath]+"//1MON/DWR-BDO/"
                    if inode == 1:
                        tdss,cunits,ctype = dssifh.read_rts(get_pathname(dssifh,path))
                    else:
                        ttss2,cunits,ctype = dssifh.read_rts(get_pathname(dssifh,path)) 
                        tdss.iloc[:,0]=tdss.iloc[:,0]+ttss2.iloc[:,0]
            path = "/CALSIM/"+C3_paths[ipath]+"_"+C3_nodes[c3j]+"/"+C3_paths[ipath]+"//1MON/L2015A/"
            dssofh.write_rts(path,tdss,cunits,ctype)

if __name__ == "__main__":
    pyhecdss.set_message_level(0)
    PP_DSS()
    daytomonth(sys.argv[1])  #DCD_Sep2016.dss
    daytomonth(sys.argv[2])  #"DP_island.dss"
    daytomonth(sys.argv[3])   #"GW_per_island.dss"
    DCD_to_CALSIM(sys.argv[1])
