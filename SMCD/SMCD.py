#********************************************************************************
#            Suisun Marsh Channel Depletion Model (SMCD)
#********************************************************************************
#<license>
#    Copyright (C) State of California, Department of Water Resources.
#    This file is part of SMCD.
#
#    SMCD is free software: 
#    you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    SMCD is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with SMCD.  If not, see <http://www.gnu.org/licenses>.
#</license>
#
#*********************************************************************************
# Enter the following line in the command window to run the model:
#    Python SMCD.py 
#********************************************************************************* 
# 
import os,sys
import shutil


def main():
    owd = os.getcwd()
    dir_dst = "../ETAW/"
    os.chdir(dir_dst)
    
    list = os.listdir(dir_dst+"Historical/") 
    Num_islands = len(list)
    
    inputfile = "Suisun_PT.csv"
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
    outputfile = "SMCD_pond_"+months[endmonth-1]+endyear+"_1.dss"
    print "output file =", outputfile, " endyear =",endyear
    
    status=os.system('python SMETAW.py')
    for ifile in range(0,2):
        filename = "ETAW_day_"+str(ifile)+".dss"
        shutil.copy(filename,owd)
    os.chdir(owd)
    for ifile in range(0,2):
       tempstr = 'python dss2dss_combine.py '+str(ifile+1)
       status=os.system(tempstr)
    for ifile in range(0,2):
        filename = "ETAW_day_"+str(ifile)+".dss"
        os.remove(filename)
        filename = "ETAW_day_"+str(ifile)+".dsc"
        os.remove(filename)
        filename = "ETAW_day_"+str(ifile)+".dsd"
        os.remove(filename)
        filename = "ETAW_day_"+str(ifile)+".dsk"
        os.remove(filename)
    status=os.system('python isl168to142.py')
    
    tempstr = 'python forNODCU3.py '+endyear+" "+str(Num_islands)    
    status=os.system(tempstr)
    
    dir_dst =".\\CD_inputs\\"
    shutil.copy("DICU5.12",dir_dst)
    shutil.copy("DICU5.14",dir_dst)
    shutil.copy("DICU5.17",dir_dst)
    shutil.copy("DICU5.27",dir_dst)
    shutil.copy("DICU5.30",dir_dst)
    
    
    dir_dst = ".\\NODCU\\CD_Cal\\DCD_outputs\\"
    os.chdir(dir_dst)
    
    os.environ['DICU5.14']='../../../CD_inputs/DICU5.14' 
    os.environ['DICU5.17']='../../../CD_inputs/DICU5.17' 
    os.environ['DICU5.12']='../../../CD_inputs/DICU5.12' 
    os.environ['DICU5.27']='../../../CD_inputs/DICU5.27' 
    os.environ['DICU5.30']='../../../CD_inputs/DICU5.30' 
    
    os.environ['GW_RATES.TXT'] = '../../../NODCU/GW_RATES.txt'   
    os.environ['DIVFCTR.RMA']='../../../NODCU/DIVFCTR.SS.2019'
    os.environ['DRNFCTR.RMA']='../../../NODCU/DRNFCTR.SS.2019'
    os.environ['LEACHAPL.DAT']='../../../NODCU/LEACHAPL_SS.DAT'
    os.environ['LEACHDRN.DAT']='../../../NODCU/LEACHDRN_SS.DAT'
    os.environ['IDRNTDS.DAT']='../../../NODCU/IDRNTDS_SS.DAT'
    os.environ['DRNCL.123']='../../../NODCU/DRNCL.SS.123'
    os.environ['GEOM-NODES']='../../../NODCU/GEOM-NODES-SS.txt'
    
    os.environ['IRREFF.DAT']='../../../NODCU/IRREFF-SS.txt'
    os.environ['subarea-info']='../../../NODCU/subarea-info-SS.txt'
    
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
    # The leach scale factor   
    os.environ['leachscale']=str(1)
    
    status = os.system('DCD_kernel.exe')
    
    status=os.system('dssts < junk1_1.txt')
    status=os.system('dssts < junk2_1.txt')
    status=os.system('dssts < junk3_1.txt')
    
    shutil.copy(outputfile,owd)
    print "finish" 

if __name__ == "__main__":
    main()