# After the ETAW for water has been calculated, it has to be removed in this file. 12/22/22015

import numpy as np
from numpy import zeros,sum,clip
import os, sys, string

def timeseries_combine(DETAWOUTPUT,Oldisl, Newisl,Numcrop,idates,ratefile):
    DETAWISL168 = zeros((6,Oldisl,idates),float)
    DETAWISLnew = zeros((5,Newisl,idates),float)
    DETAWISLold = np.copy(DETAWOUTPUT)
    ratelines = []
    if ratefile != "":
        f0 = open(ratefile, 'r')
        for line in f0:
            if line:
                if len(line) > 2:
                    ratelines.append(line)
        rates = zeros((len(ratelines),15),float)
        ex_isl = zeros((len(ratelines)),int)
        for i in range(len(ratelines)):
            temp = int(ratelines[i].split(",")[0].strip())
            if temp > 0: 
                ex_isl[i]=temp
                for j in range(Numcrop):
                    rates[i,j] = float(ratelines[i].split(",")[3+j].strip())
        for isl in range(0,Oldisl):
            isign = 0
            for ii in range(len(ex_isl)):
                if (isl+1) == ex_isl[ii]:
                    isign = 1
                    for j in range(Numcrop):
                        DETAWISLold[:,isl,j,:] = DETAWISLold[:,isl,j,:]*rates[ii,j]
                    break
            if isign == 0:
                DETAWISLold[:,isl,:,:] = 0.0        
    #sum up the island ET of all 15 crops for each island
    for iterm in range(0,6):  #len(DETAWOUTPUT)):  
        for iland in range(0,Oldisl):     #168):
            if iterm != 4:
                if iterm == 3:
                    tempcrop = Numcrop-1
                else:
                    tempcrop = Numcrop
                for idd in range(0,idates):
                    DETAWISL168[iterm,iland,idd] = sum(DETAWISLold[iterm,iland,0:tempcrop,idd])
            else:
                DETAWISL168[4,iland,:] = DETAWISLold[0,iland,Numcrop-1,:]-DETAWISLold[2,iland,Numcrop-1,:]
                DETAWISL168[4,iland,:] = clip(DETAWISL168[4,iland,:],0.0,None)
                
    #DETAWISL168 = clip(DETAWISL168,0.0,None)
    #convert 168 islands ET to 142 islands ET
    for iterm in range(0,5): #5):
        for iland in range(0,Newisl):
            if iterm == 2:
                DETAWISLnew[iterm,iland,:] = DETAWISL168[iterm,iland,:]-DETAWISL168[5,iland,:]   
                DETAWISLnew[iterm,iland,:] = clip(DETAWISLnew[iterm,iland,:],0.0,None)                   
            else:
                DETAWISLnew[iterm,iland,:] = clip(DETAWISL168[iterm,iland,:],0.0,None)
    return DETAWISLnew
    
def forNODCU(DETAWISL168,inputversion,endyear,ilands,outfilenames):     
    #prepare the text input files for DCD-NODCU
    #Cpart[0] = "ETAW"
    #Cpart[1] = "ESPG"
    #Cpart[2] = "ETC"
    #Cpart[3] = "PCP"
    #Cpart[4] = "DSW"
    
    beginyear = 1922
    tyr = endyear-beginyear+1
    
    #dssfh=pyhecdss.DSSFile(inputfile)
    tempname = "./Output/DICU5"+outfilenames
    f1 = open(tempname+".27","w")
    f2 = open(tempname+".14","w")
    f3 = open(tempname+".12","w")
    f4 = open(tempname+".17","w")
    f5 = open(tempname+".30","w")
    
    if inputversion.strip() == "CALSIM3":
        daysfile = "./Input/planning_study/calender.txt"    #update the txt file too!!!!
    else:
        daysfile = "./Input/historical_study/calender.txt"
    f0 = open(daysfile)
    daysofyear = zeros((366,4,tyr),int)
    isl = 0
    idays = 0
    iyr = 0
    for line in f0:
        if line:
            if isl>0:
                ll = line.split()
                if int(ll[1])==10 and int(ll[2])==1:
                    idays = 0
                    iyr = int(ll[0])+1-beginyear
                daysofyear[idays,0,iyr] = int(ll[0])
                daysofyear[idays,1,iyr] = int(ll[1])
                daysofyear[idays,2,iyr] = int(ll[2])
                daysofyear[idays,3,iyr] = int(ll[3])
            isl = isl + 1
            idays += 1
    f0.close()
    
    
    for ifile in range(0,5): 
        print(ifile)
        for iland in range(1,ilands+1):
            if ifile == 3:
                strt = str("%3i" % iland)+"AT1 4 HISTORIC DEPLETION OF APPLIED WATER BY IRR. AND URBAN, AREA   "
                strt = strt + str(iland)+ "\n"
                f4.writelines(strt)
                strt = str("%3i" % iland)+"AT2            FROM DETAW STUDY \n"
                f4.writelines(strt)                
            elif ifile == 1:
                strt = str("%3i" % iland)+"AT1 4  TOTAL HISTORIC CU OF SEEPAGE  "
                strt = strt + str(iland)+ "\n"
                f2.writelines(strt)
                strt = str("%3i" % iland)+"AT2            FROM DETAW STUDY \n"
                f2.writelines(strt)
            elif ifile == 0:
                strt = str("%3i" % iland)+"AT1 4   HISTORIC DEPLETION, AREA   "
                strt = strt + str(iland)+ "\n"
                f1.writelines(strt)
                strt = str("%3i" % iland)+"AT2            FROM DETAW STUDY \n"
                f1.writelines(strt)
            elif ifile == 2:
                strt = str("%3i" % iland)+"AT1 4 TOTAL BASIN RUNOFF, AREA   "
                strt = strt + str(iland)+ "\n"
                f3.writelines(strt)
                strt = str("%3i" % iland)+"AT2            FROM DETAW STUDY \n"
                f3.writelines(strt)
            elif ifile == 4:
                strt = str("%3i" % iland)+"AT1 4 TOTAL NET DEPLETION of WATERBODY, AREA   "
                strt = strt + str(iland)+ "\n"
                f5.writelines(strt)
                strt = str("%3i" % iland)+"AT2            FROM DETAW STUDY \n"
                f5.writelines(strt) 
                
            ndays = 0
            tempts = zeros((tyr,366),float)
            for iyr in range(0,tyr):
                for iday in range(0,366):
                    if (iday == 365 and daysofyear[iday,0,iyr] != 0) or (iday<365):
                        if ifile == 2:
                            if iyr == 0 and iday < 5:
                                for kday in range(0, iday+1):
                                    temp_c = 0.2
                                    tempts[iyr,iday] += DETAWISL168[ifile,iland-1,iday-kday]*temp_c
                            else:
                                for kday in range(0,5):
                                    temp_c = 0.2
                                    tempts[iyr,iday] += DETAWISL168[ifile,iland-1,ndays-kday]*temp_c
                        else:
                            tempts[iyr,iday] = DETAWISL168[ifile,iland-1,ndays]
                        ndays+=1  
            for iyr in range(0,tyr):
                strt = str("%3i" % iland)+"A   "+ str(beginyear+iyr)
                for iday in range(0,366):
                    if (iday == 365 and daysofyear[iday,0,iyr] != 0) or (iday<365):
                        if ifile==4 and inputversion.strip() == "SCHISM":
                            tempts[iyr][iday] = 0.0  # Set surface water evaporation as zero for SCHISM  11/28/2017
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
