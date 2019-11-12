# After the ETAW for water has been calculated, it has to be removed in this file. 12/22/22015

from numpy import zeros,sum,clip
import os, sys, string

def timeseries_combine(DETAWOUTPUT,Oldisl, Newisl,Numcrop,idates,inputversion):
    #DETAWOUTPUT [0-ETC,1-ESPG,2-PCP,3-ETAW,4-DSW,5-ER] unit: A-FT 
    if inputversion.strip() == "CALSIM3":
        islandfile = ".\Input\planning_study\island_id.txt"
    else:
        islandfile = ".\Input\historical_study\island_id.txt"
    DETAWISL168 = zeros((6,Oldisl,idates),float)
    DETAWISL142 = zeros((5,Newisl,idates),float)
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
    
    #sum up the island ET of all 15 crops for each island
    for iterm in range(0,6):  #len(DETAWOUTPUT)):  
        for iland in range(0,Oldisl):     #168):
            if iterm != 4:
                if iterm == 3:
                    tempcrop = Numcrop-1
                else:
                    tempcrop = Numcrop
                for idd in range(0,idates):
                    DETAWISL168[iterm,iland,idd] = sum(DETAWOUTPUT[iterm,iland,:,idd])
            else:
                DETAWISL168[4,iland,:] = DETAWOUTPUT[0,iland,Numcrop-1,:]-DETAWOUTPUT[2,iland,Numcrop-1,:]
                DETAWISL168[4,iland,:] = clip(DETAWISL168[4,iland,:],0.0,None)
    #DETAWISL168 = clip(DETAWISL168,0.0,None)
    #convert 168 islands ET to 142 islands ET
    for iterm in range(0,5): #5):
        for iland in range(0,Oldisl):
            if (iland+1) < Newisl+1:
                if iterm == 2:
                    DETAWISL142[iterm,iland,:]=DETAWISL168[iterm,iland,:]-DETAWISL168[5,iland,:]   
                    DETAWISL142[iterm,iland,:] = clip(DETAWISL142[iterm,iland,:],0.0,None)                   
                else:
                    DETAWISL142[iterm,iland,:] = DETAWISL168[iterm,iland,:]
            else:
                for id in range(0, isl):
                    if islid[id,0] == (iland+1):
                        tempisl = islid[id,1]-1
                        if iterm == 2:
                            DETAWISL142[iterm,tempisl,:]+=(DETAWISL168[iterm,iland,:]-DETAWISL168[5,iland,:])
                            DETAWISL142[iterm,tempisl,:]=clip(DETAWISL142[iterm,tempisl,:],0.0,None)
                        else:
                            DETAWISL142[iterm,tempisl,:]+=DETAWISL168[iterm,iland,:]
    #DETAWISL142=clip(DETAWISL142,0.0,None)
    return DETAWISL142
    
def forNODCU(DETAWISL142,inputversion,endyear):     
    #prepare the text input files for DCD-NODCU
    #Cpart[0] = "ETAW"
    #Cpart[1] = "ESPG"
    #Cpart[2] = "ETC"
    #Cpart[3] = "PCP"
    #Cpart[4] = "DSW"
    
    beginyear = 1922
    tyr = endyear-beginyear+1
    
    #dssfh=pyhecdss.DSSFile(inputfile)
    f1 = open(".\Output\DICU5.27","w")
    f2 = open(".\Output\DICU5.14","w")
    f3 = open(".\Output\DICU5.12","w")
    f4 = open(".\Output\DICU5.17","w")
    f5 = open(".\Output\DICU5.30","w")
    
    if inputversion.strip() == "CALSIM3":
        daysfile = ".\Input\planning_study\calender.txt"    #update the txt file too!!!!
    else:
        daysfile = ".\Input\historical_study\calender.txt"
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
        for iland in range(1,143):
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
                                    tempts[iyr,iday] += DETAWISL142[ifile,iland-1,iday-kday]*temp_c
                            else:
                                for kday in range(0,5):
                                    temp_c = 0.2
                                    tempts[iyr,iday] += DETAWISL142[ifile,iland-1,ndays-kday]*temp_c
                        else:
                            tempts[iyr,iday] = DETAWISL142[ifile,iland-1,ndays]
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
