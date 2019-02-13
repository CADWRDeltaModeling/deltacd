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

import os,sys
import shutil


def main():
    owd = os.getcwd()
    dir_dst = "../DETAWv2.0/"
    os.chdir(dir_dst)
    
    inputfile = "LODI_PT.csv"
    months = ["Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"]
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
    outputfile = outputfile = "DCD_"+months[endmonth-1]+endyear+".dss"
    print "output file =", outputfile
    
    status=os.system('python DETAW.py')
    for ifile in range(0,17):
        filename = "DETAW_day_"+str(ifile)+".dss"
        shutil.copy(filename,owd)
    os.chdir(owd)
    for ifile in range(0,17):
        tempstr = 'python dss2dss_combine.py '+str(ifile+1)
        status=os.system(tempstr)
    for ifile in range(0,17):
        filename = "DETAW_day_"+str(ifile)+".dss"
        os.remove(filename)
        filename = "DETAW_day_"+str(ifile)+".dsc"
        os.remove(filename)
        filename = "DETAW_day_"+str(ifile)+".dsd"
        os.remove(filename)
        filename = "DETAW_day_"+str(ifile)+".dsk"
        os.remove(filename)
    status=os.system('python isl168to142.py')
    status=os.system('python forNODCU3.py')
    
    dir_dst =".\\DCD_inputs\\"
    shutil.copy("DICU5.12",dir_dst)
    shutil.copy("DICU5.14",dir_dst)
    shutil.copy("DICU5.17",dir_dst)
    shutil.copy("DICU5.27",dir_dst)
    shutil.copy("DICU5.30",dir_dst)
    
    
    dir_dst = ".\\NODCU\\DCD_Cal\\DCD_outputs\\"
    os.chdir(dir_dst)
    
    os.environ['DICU5.14']='../../../DCD_inputs/DICU5.14' 
    os.environ['DICU5.17']='../../../DCD_inputs/DICU5.17' 
    os.environ['DICU5.12']='../../../DCD_inputs/DICU5.12' 
    os.environ['DICU5.27']='../../../DCD_inputs/DICU5.27' 
    os.environ['DICU5.30']='../../../DCD_inputs/DICU5.30' 
    
    os.environ['GW_RATES.TXT'] = '../../../NODCU/GW_RATES.TXT'   #update data in the file each year ----no adjustment-GW_RATES.TXT
    os.environ['GW_LOWLANDS.TXT']='../../../NODCU/GW_LOWLANDS.TXT' #set for DETAW-CD
    os.environ['DIVFCTR.RMA']='../../../NODCU/DIVFCTR.DSM.2-92'
    os.environ['DRNFCTR.RMA']='../../../NODCU/DRNFCTR.DSM.2-92'
    os.environ['LEACHAPL.DAT']='../../../NODCU/LEACHAPL.DAT'
    os.environ['LEACHDRN.DAT']='../../../NODCU/LEACHDRN.DAT'
    os.environ['IDRNTDS.DAT']='../../../NODCU/IDRNTDS.DAT'
    os.environ['DRNCL.123']='../../../NODCU/DRNCL.123'
    os.environ['GEOM-NODES']='../../../NODCU/GEOM-NODES-1.5'
    
    os.environ['IRREFF.DAT']='../../../NODCU/IRREFF-3MWQIregions'
    os.environ['subarea-info']='../../../NODCU/subarea-info'
    
    
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
    
    #date
    status=os.system('DETAW_CD.exe')
    
    status=os.system('dssts < junk1_1.txt')
    status=os.system('dssts < junk1_2.txt')
    status=os.system('dssts < junk2_1.txt')
    status=os.system('dssts < junk2_2.txt')
    status=os.system('dssts < junk3_1.txt')
    status=os.system('dssts < junk3_2.txt')
    
    shutil.copy(outputfile,owd)
    print "finish"

if __name__ == "__main__":
    main()