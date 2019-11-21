# Delta Channel Depletion Model (DCDv1.0)
#<license>
#    Copyright (C) State of California, Department of Water Resources.
#    This file is part of DCDv1.0.

#    DCDv1.0 is free software: 
#    you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    DCDv1.0 is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with DCDv1.0.  If not, see <http://www.gnu.org/licenses>.
#</license>
#
# Enter the following line in the command window to run the model:
#    Python DCD.py 
# 

import pandas as pd
import pyhecdss
import os,sys
import shutil


def main():
    owd = os.getcwd()
    modelparafile = os.path.join("NODCU","DCD_parameters.inp")
    fmp = open(modelparafile,"r")
    for line in fmp:
        if line:
            if not("#" in line):
                supmodel = int(line.split("=")[1]) 
                
    dir_dst = os.path.join("..","DETAW")
    os.chdir(dir_dst)
    months = ["Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"]
    
    if supmodel == 1:
        inputfile = os.path.join("Input","historical_study","LODI_PT.csv")
        f0 = open(inputfile, 'r')
        templ = ""
        endyear = 0
        for line in f0:
            if line:
                templ = line
                if line.split(",")[5].strip()=="0" and line.split(",")[6].strip()=="0":
                    break
        endyear = templ.split(",")[1]  
        endmonth = int(templ.split(",")[2])
        outputfile = "DCD_"+months[endmonth-1]+endyear+".dss"        
    elif supmodel == 3:
        inputfile = os.path.join("Input","planning_study","LODI_PT.csv")
        f0 = open(inputfile, 'r')
        templ = ""
        endyear = 0
        for line in f0:
            if line:
                templ = line
                if line.split(",")[5].strip()=="0" and line.split(",")[6].strip()=="0":
                    break
        endyear = templ.split(",")[1]  
        endmonth = int(templ.split(",")[2])
        outputfile = "DCD_"+months[endmonth-1]+endyear+".dss"        
    elif supmodel == 2:
        inputfile = os.path.join("Input","historical_study","LODI_PT.csv")
        f0 = open(inputfile, 'r')
        templ = ""
        endyear = 0
        for line in f0:
            if line:
                templ = line
                if line.split(",")[5].strip()=="0" and line.split(",")[6].strip()=="0":
                    break
        endyear = templ.split(",")[1]  
        endmonth = int(templ.split(",")[2])
        outputfile = "DCD_noWS_"+months[endmonth-1]+endyear+".dss"
    status=os.system('python DETAW.py')
    print("output file =", outputfile)
        
    os.chdir(owd)
    exe_file='./DETAW_CD' if os.name=='posix' else 'DETAW_CD.exe'
    dir_dst = os.path.join("NODCU","DCD_Cal")
    os.chdir(dir_dst)
    if not os.path.exists("DCD_outputs"): os.mkdir("DCD_outputs")
    shutil.copy(exe_file,"DCD_outputs")
    shutil.copy('WYTYPES',"DCD_outputs")
    os.chdir("DCD_outputs")
    
    detaw_output_dir=os.path.join('..','..','..','..','DETAW','Output')
    os.environ['DICU5.14']=os.path.join(detaw_output_dir,'DICU5.14')
    os.environ['DICU5.17']=os.path.join(detaw_output_dir,'DICU5.17')
    os.environ['DICU5.12']=os.path.join(detaw_output_dir,'DICU5.12')
    os.environ['DICU5.27']=os.path.join(detaw_output_dir,'DICU5.27')
    os.environ['DICU5.30']=os.path.join(detaw_output_dir,'DICU5.30')
    
    nodcu_dir=os.path.join('..','..','..','NODCU')
    os.environ['GW_RATES.TXT'] = os.path.join(nodcu_dir,'GW_RATES.txt')   #update data in the file each year ----no adjustment-GW_RATES.TXT
    os.environ['GW_LOWLANDS.TXT']=os.path.join(nodcu_dir,'GW_LOWLANDS.txt') #set for DETAW-CD
    os.environ['DIVFCTR.RMA']=os.path.join(nodcu_dir,'DIVFCTR.DSM.2-92')
    os.environ['DRNFCTR.RMA']=os.path.join(nodcu_dir,'DRNFCTR.DSM.2-92')
    os.environ['LEACHAPL.DAT']=os.path.join(nodcu_dir,'LEACHAPL.DAT')
    os.environ['LEACHDRN.DAT']=os.path.join(nodcu_dir,'LEACHDRN.DAT')
    os.environ['IDRNTDS.DAT']=os.path.join(nodcu_dir,'IDRNTDS.DAT')
    os.environ['DRNCL.123']=os.path.join(nodcu_dir,'DRNCL.123')
    os.environ['GEOM-NODES']=os.path.join(nodcu_dir,'GEOM-NODES-1.5')
    os.environ['IRREFF.DAT']=os.path.join(nodcu_dir,'IRREFF-3MWQIregions')
    os.environ['subarea-info']=os.path.join(nodcu_dir,'subarea-info')
    
    
    # Runtime variables
    # The years assumed are incorrect, so 'N'
    os.environ['years_ok']='N'
    # The correct beginning year to run is
    os.environ['begwy']='1922'
    # The correct last year to run is
    os.environ['endwy']=endyear                        
    # Type of drainage concentration data (1 for TDS, 2 for chloride)
    os.environ['datatype']='1'
    # Do you want to creat an ascii file?
    os.environ['ascii']='Y'
    # The dss file to save output
    os.environ['dssfile']=outputfile        
    
    status=os.system(exe_file)
    
    status=os.system('python ../converttoDSS.py junk1_1.txt')
    status=os.system('python ../converttoDSS.py junk1_2.txt')
    status=os.system('python ../converttoDSS.py junk2_1.txt')
    status=os.system('python ../converttoDSS.py junk2_2.txt')
    status=os.system('python ../converttoDSS.py junk3_1.txt')
    status=os.system('python ../converttoDSS.py junk3_2.txt')
    if supmodel == 1:
        shutil.copy(outputfile,os.path.join(owd,"Output","DSM2"))
    elif supmodel == 2:
        shutil.copy(outputfile,os.path.join(owd,"Output","SCHISM"))
    elif supmodel == 3:
        status=os.system('python ../converttoDSS.py roisl.txt')
        status=os.system('python ../converttoDSS.py gwbyisl.txt')
        tempstr = "python "+os.path.join("..","DCD_post-process_C3.py")+" "+outputfile+" DP_island.dss GW_per_island.dss"
        status = os.system(tempstr)
        tempfile = outputfile.split(".")[0].strip()+"_mon_C3.dss"
        shutil.copy(tempfile,os.path.join(owd,"Output","CALSIM3"))
        tempfile = outputfile.split(".")[0].strip()+"_mon.dss"
        shutil.copy(tempfile,os.path.join(owd,"Output","CALSIM3"))
        tempfile = "DP_island_mon.dss"
        shutil.copy(tempfile,os.path.join(owd,"Output","CALSIM3"))
        tempfile = "GW_per_island_mon.dss"
        shutil.copy(tempfile,os.path.join(owd,"Output","CALSIM3"))
    os.chdir("..")
    shutil.rmtree("DCD_outputs")

    print("finish")

if __name__ == "__main__":
    main()
