
# Import VTools time and numpy array creation function
from datetime import datetime, timedelta
import pyhecdss
import pandas as pd
import numpy as np
import os
import sys
import string
import calendar
import shutil


def get_rts(file, path):
    gen = pyhecdss.get_ts(file, path)
    return list(gen)


def get_pathname(dssfh, path, partno):
    '''
    reads catalog to get full pathname if it exists
    The assumption is that the path provided does not have a time window
    returns the first matching pathname with exactly the A,B,C,E and F parts
    '''
    pathparts = path.split('/')
    dfcat = dssfh.read_catalog()
    if partno == 5:
        dfpath = dfcat[(dfcat.A == pathparts[1]) & (dfcat.B == pathparts[2])
                       & (dfcat.C == pathparts[3]) & (dfcat.E == pathparts[5])
                       & (dfcat.F == pathparts[6])]
    elif partno == 1:
        dfpath = dfcat[(dfcat.B == pathparts[2])]
    elif partno == 2:
        dfpath = dfcat[(dfcat.A == pathparts[1]) & (dfcat.B == pathparts[2])]
    elif partno == 4:
        dfpath = dfcat[(dfcat.A == pathparts[1]) & (dfcat.B == pathparts[2])
                       & (dfcat.E == pathparts[5]) & (dfcat.F == pathparts[6])]
    pathname = dssfh.get_pathnames(dfpath)[0]
    return pathname


def PP_DSS():
    inputfile = "./RO_island.dss"
    outputfile = "./DP_island.dss"

    pyhecdss.set_message_level(0)
    dssfh = pyhecdss.DSSFile(inputfile)
    dssofh = pyhecdss.DSSFile(outputfile, create_new=True)
    plist = dssfh.get_pathnames()
    for p in plist:
        df, u, p = dssfh.read_rts(p)
        df.values[:] = df.values[:]*0.25
        dssofh.write_rts(df.columns[0].replace(
            "RO-FLOW", "DP-FLOW"), df, u, p)  # 'PER-AVER')
    dssfh.close()
    dssofh.close()


def daytomonth(inputfile):
    d = pyhecdss.DSSFile(inputfile)
    outputfile = os.getcwd()+"/"+inputfile.split(".")[0]+"_mon.dss"
    do = pyhecdss.DSSFile(outputfile, create_new=True)
    plist = d.get_pathnames()
    for p in plist:
        df, u, p = d.read_rts(p)
        do.write_rts(df.columns[0].replace('1DAY', '1MON'),
                     df.resample('M').mean(), u, 'PER-AVER')
    d.close()
    do.close()


def DCD_to_CALSIM_ISLAND(divfile, spgfile, drnfile, rofile, inputfile):
    inputfile = inputfile.split(".")[0]+"_mon.dss"
    outputfile = inputfile.split(".")[0]+"_C3.dss"
    dssofh = pyhecdss.DSSFile(outputfile, create_new=True)

    DCD_C3_islands = "../DCD_CALSIM3_islands_N.csv"
    C3_nodes = ["OMR", "SJR_EAST", "SJR_WEST", "SAC_WEST",
                "MOK", "SAC_SOUTH", "SAC_NORTH", "50_PA2"]
    C3_paths = ["IRR", "SEEP", "DRN"]
    DSM2N_paths = ["DIV-FLOW", "SEEP-FLOW", "DRAIN-FLOW"]
    DCD_paths = ["DIV-WO-SPG-FLOW", "SPG-FLOW", "DRN-WO-RO-FLOW", "RO-FLOW"]

    f0 = open(DCD_C3_islands)
    DCD_islands = []
    ili = 0
    for line in f0:
        ili += 1
        if line:
            if ili > 1:
                DCD_islands.append(line)
    f0.close()
    for ipath in range(0, len(C3_paths)):
        if ipath == 0:
            dssifh = pyhecdss.DSSFile(divfile.split(".")[0]+"_mon.dss")
        elif ipath == 1:
            dssifh = pyhecdss.DSSFile(spgfile.split(".")[0]+"_mon.dss")
        elif ipath == 2:
            dssifh = pyhecdss.DSSFile(drnfile.split(".")[0]+"_mon.dss")
            dssifh2 = pyhecdss.DSSFile(rofile.split(".")[0]+"_mon.dss")
        for c3j in range(0, len(C3_nodes)-1):
            iisland = 0
            for i in range(0, len(DCD_islands)):
                if C3_nodes[c3j] == DCD_islands[i].split(",")[1].strip():
                    iisland += 1
                    tempIsl = int(DCD_islands[i].split(",")[0].strip())

                    if ipath == 0 or ipath == 1:
                        path = "/DICU-ISLAND/" + \
                            DCD_islands[i].split(",")[0].strip(
                            )+"/"+DCD_paths[ipath]+"//1MON/DWR-BDO/"
                        if iisland == 1:
                            tdss, cunits, ctype = dssifh.read_rts(
                                get_pathname(dssifh, path, 5))
                        else:
                            ttss2, cunits, ctype = dssifh.read_rts(
                                get_pathname(dssifh, path, 5))
                            tdss.iloc[:, 0] += ttss2.iloc[:, 0]
                    elif ipath == 2:
                        path = "/DICU-ISLAND/" + \
                            DCD_islands[i].split(",")[0].strip(
                            )+"/"+DCD_paths[ipath]+"//1MON/DWR-BDO/"
                        path2 = "/DICU-ISLAND/" + \
                            DCD_islands[i].split(",")[0].strip(
                            )+"/"+DCD_paths[ipath+1]+"//1MON/DWR-BDO/"
                        if iisland == 1:
                            tdss, cunits, ctype = dssifh.read_rts(
                                get_pathname(dssifh, path, 5))
                            tdss_ro, cunits, ctype = dssifh2.read_rts(
                                get_pathname(dssifh2, path2, 5))
                            tdss.iloc[:, 0] += tdss_ro.iloc[:, 0]
                        else:
                            ttss2, cunits, ctype = dssifh.read_rts(
                                get_pathname(dssifh, path, 5))
                            ttss_ro2, cunits, ctype = dssifh2.read_rts(
                                get_pathname(dssifh2, path2, 5))
                            tdss.iloc[:, 0] += ttss2.iloc[:, 0] + \
                                ttss_ro2.iloc[:, 0]
            path = "/CALSIM/" + \
                C3_paths[ipath]+"_"+C3_nodes[c3j] + \
                "/"+C3_paths[ipath]+"//1MON/L2015A/"
            dssofh.write_rts(path, tdss, cunits, ctype)
        dssifh = pyhecdss.DSSFile(inputfile)
        pathin = "/DICU-HIST+RSVR/BBID/"+DSM2N_paths[ipath]+"//1MON/DWR-BDO/"
        tdssb, cunits, ctype = dssifh.read_rts(get_pathname(dssifh, pathin, 5))
        pathout = "/CALSIM/" + \
            C3_paths[ipath]+"_" + \
            C3_nodes[len(C3_nodes)-1]+"/"+C3_paths[ipath]+"//1MON/L2015A/"
        dssofh.write_rts(pathout, tdssb, cunits, ctype)
        dssifh.close()
        if ipath == 2:
            dssifh2.close()
    dssofh.close()


def split_BBID(divfile, spgfile, drnfile, rofile, outputfile, option):
    Tisland = 168
    DCD_paths = ["DIV-WO-SPG-FLOW", "SPG-FLOW", "DRN-WO-RO-FLOW", "RO-FLOW"]

    inputfiles = [divfile, spgfile, drnfile, rofile]
    # Reduce BBID amounts from the island outputs and add BBID into the island outputs
    BBIDisl = [33, 34, 41, 103, 128, 130]
    for ifile in range(len(inputfiles)):
        extfile = "ext_" + inputfiles[ifile][2::]
        orgfile = inputfiles[ifile]
        dssout = pyhecdss.DSSFile(orgfile, create_new=True)
        for i in range(len(BBIDisl)):
            path1 = "/DICU-ISLAND/" + \
                str(BBIDisl[i])+"/"+DCD_paths[ifile]+"//1DAY/DWR-BDO/"
            path2 = "/BBID/"+str(BBIDisl[i])+"/"+DCD_paths[ifile]+"//1DAY//"
            tdss1 = get_rts(orgfile, path1)
            tdss2 = get_rts(extfile, path2)
            pathout = "/BBID/"+str(BBIDisl[i]) + \
                "/"+DCD_paths[ifile]+"//1DAY/DWR-BDO/"
            dssout.write_rts(pathout, tdss2[0][0], tdss2[0][1], tdss2[0][2])
            pathout = path1.replace(str(BBIDisl[i]), str(BBIDisl[i])+"_w_BBID")
            dssout.write_rts(pathout, tdss1[0][0], tdss1[0][1], tdss1[0][2])
            tdss1[0][0].iloc[:, 0] = tdss1[0][0].iloc[:, 0] - \
                tdss2[0][0].iloc[:, 0]
            dssout.write_rts(path1, tdss1[0][0], tdss1[0][1], tdss1[0][2])

        dssifh2 = pyhecdss.DSSFile(extfile, create_new=True)
        dfcat = dssifh2.read_catalog()
        dfpath = dfcat[(dfcat.A != "BBID")]
        pathnames = dssifh2.get_pathnames(dfpath)
        for i in range(len(pathnames)):
            path1 = "/DICU-ISLAND/" + \
                pathnames[i].split("/")[2]+"/" + \
                DCD_paths[ifile]+"//1DAY/DWR-BDO/"
            tdss1, cunits, ctype = dssout.read_rts(
                get_pathname(dssout, path1, 5))
            tdss2, cunits, ctype = dssifh2.read_rts(pathnames[i])
            if option == 2:
                tdss1.iloc[:, 0] = tdss1.iloc[:, 0]-tdss2.iloc[:, 0]
            dssout.write_rts(path1, tdss1, cunits, ctype)
        dssifh2.close()
        dssout.close()


def islandtoDSM2node(divfile, spgfile, drnfile, rofile, outputfile):
    Tisland = 168
    DCD_paths = ["DIV-WO-SPG-FLOW", "SPG-FLOW", "DRN-WO-RO-FLOW", "RO-FLOW"]
    divratefile = '../../../NODCU/DIVFCTR_CS3_NorOMR.2020'
    drnratefile = '../../../NODCU/DRNFCTR_CS3_NorOMR.2020'

    inputfiles = [divfile, spgfile, drnfile, rofile]
    # Reduce BBID amounts from the island outputs and add BBID into the island outputs
    BBIDisl = [33, 34, 41, 103, 128, 130]

    # Allocate island values to DSM2 nodes

    divalloc = 687
    drnalloc = 427
    divisl = [0]*divalloc
    divnode = [0]*divalloc
    divrate = [0.0]*divalloc
    drnisl = [0]*drnalloc
    drnnode = [0]*drnalloc
    drnrate = [0.0]*drnalloc
    f1 = open(divratefile)
    ili = 0
    maxnode = 0
    for line in f1:
        ili += 1
        if line:
            if ili > 4:
                if line.strip() != "":
                    if int(line[0:5]) > 0 and int(line[5:11]) > 0:
                        divisl[ili-5] = int(line[0:5])
                        divnode[ili-5] = int(line[5:11])
                        divrate[ili-5] = float(line[12:len(line)])*0.01
    f1.close()
    f2 = open(drnratefile)
    ili = 0
    for line in f2:
        ili += 1
        if line:
            if ili > 4:
                if line.strip() != "":
                    if int(line[0:5]) > 0 and int(line[5:11]) > 0:
                        drnisl[ili-5] = int(line[0:5])
                        drnnode[ili-5] = int(line[5:11])
                        drnrate[ili-5] = float(line[12:len(line)])*0.01
    f2.close()
    nodes = []
    for i in range(len(divnode)):
        if divnode[i] not in nodes:
            nodes.append(divnode[i])
    for i in range(len(drnnode)):
        if drnnode[i] not in nodes:
            nodes.append(drnnode[i])
    Sortednodes = np.sort(nodes)
    dssout = pyhecdss.DSSFile(outputfile, create_new=True)
    for ifile in range(len(inputfiles)-1):
        orgfile = inputfiles[ifile]
        dssinputf = pyhecdss.DSSFile(orgfile)
        if ifile == 0 or ifile == 1:
            for i in range(len(Sortednodes)):
                nonode = 0
                if Sortednodes[i] > 0:
                    for j in range(len(divnode)):
                        if Sortednodes[i] == divnode[j]:
                            nonode += 1
                            pathisl = "/DICU-ISLAND/"+str(divisl[j])+"/////"
                            if nonode == 1:
                                tdss1, cunits, ctype = dssinputf.read_rts(
                                    get_pathname(dssinputf, pathisl, 2))
                                tdss1.iloc[:, 0] = tdss1.iloc[:, 0]*divrate[j]
                            else:
                                tdss2, cunits, ctype = dssinputf.read_rts(
                                    get_pathname(dssinputf, pathisl, 2))
                                tdss1.iloc[:, 0] = tdss1.iloc[:, 0] + \
                                    tdss2.iloc[:, 0]*divrate[j]
                    if nonode == 0:
                        pathisl = "/DICU-ISLAND/1/////"
                        tdss1, cunits, ctype = dssinputf.read_rts(
                            get_pathname(dssinputf, pathisl, 2))
                        tdss1.iloc[:, 0] = tdss1.iloc[:, 0]*0.0
                    if ifile == 0:
                        pathout = "/DICU-HIST+NODE/" + \
                            str(Sortednodes[i])+"/DIV-FLOW//1DAY/DWR-BDO/"
                    elif ifile == 1:
                        pathout = "/DICU-HIST+NODE/" + \
                            str(Sortednodes[i])+"/SEEP-FLOW//1DAY/DWR-BDO/"
                    dssout.write_rts(pathout, tdss1, cunits, ctype)
        elif ifile == 2:
            orgfile2 = inputfiles[ifile+1]
            dssinputf2 = pyhecdss.DSSFile(orgfile2)
            for i in range(len(Sortednodes)):
                nonode = 0
                if Sortednodes[i] > 0:
                    for j in range(len(drnnode)):
                        if Sortednodes[i] == drnnode[j]:
                            nonode += 1
                            pathisl = "/DICU-ISLAND/"+str(drnisl[j])+"/////"
                            if nonode == 1:
                                tdss1, cunits, ctype = dssinputf.read_rts(
                                    get_pathname(dssinputf, pathisl, 2))
                                tdssro, cunits, ctype = dssinputf2.read_rts(
                                    get_pathname(dssinputf2, pathisl, 2))
                                tdss1.iloc[:, 0] = tdss1.iloc[:, 0] * \
                                    drnrate[j]+tdssro.iloc[:, 0]*drnrate[j]
                            else:
                                tdss2, cunits, ctype = dssinputf.read_rts(
                                    get_pathname(dssinputf, pathisl, 2))
                                tdssro, cunits, ctype = dssinputf2.read_rts(
                                    get_pathname(dssinputf2, pathisl, 2))
                                tdss1.iloc[:, 0] = tdss1.iloc[:, 0]+tdss2.iloc[:,
                                                                               0]*drnrate[j]+tdssro.iloc[:, 0]*drnrate[j]
                    if nonode == 0:
                        pathisl = "/DICU-ISLAND/1/////"
                        tdss1, cunits, ctype = dssinputf.read_rts(
                            get_pathname(dssinputf, pathisl, 2))
                        tdss1.iloc[:, 0] = tdss1.iloc[:, 0]*0.0
                    pathout = "/DICU-HIST+NODE/" + \
                        str(Sortednodes[i])+"/DRAIN-FLOW//1DAY/DWR-BDO/"
                    dssout.write_rts(pathout, tdss1, cunits, ctype)
            dssinputf2.close()
        for i in range(len(BBIDisl)):
            pathname = "/BBID/" + \
                str(BBIDisl[i])+"/"+DCD_paths[ifile]+"//1DAY/DWR-BDO/"
            if i == 0:
                print(pathname, inputfiles[ifile])
                tdssb = get_rts(inputfiles[ifile], pathname)
                #print("1st Tdssb =",tdssb)
            else:
                tdsst = get_rts(inputfiles[ifile], pathname)
                tdssb[0][0].iloc[:, 0] += tdsst[0][0].iloc[:, 0]
            if ifile == 2:
                pathname = "/BBID/" + \
                    str(BBIDisl[i])+"/"+DCD_paths[ifile+1]+"//1DAY/DWR-BDO/"
                tdssro = get_rts(inputfiles[ifile+1], pathname)
                tdssb[0][0].iloc[:, 0] += tdssro[0][0].iloc[:, 0]

        if ifile == 0:
            pathout = "/DICU-HIST+RSVR/BBID/DIV-FLOW//1DAY/DWR-BDO/"
        elif ifile == 1:
            pathout = "/DICU-HIST+RSVR/BBID/SEEP-FLOW//1DAY/DWR-BDO/"
        elif ifile == 2:
            pathout = "/DICU-HIST+RSVR/BBID/DRAIN-FLOW//1DAY/DWR-BDO/"
        # FIXME: why shift by 1 day? to match older results
        dssout.write_rts(
            pathout, tdssb[0][0].shift(-1), tdssb[0][1], tdssb[0][2])
        dssinputf.close()
    dssout.close()


def dcd_postprocess_base():
    """ Process the base extension """
    PP_DSS()
    daytomonth("DP_island.dss")
    daytomonth("GW_per_island.dss")


def dcd_postprocess_ex3(outputfile):
    """ Postprocess with the ex3 extension option """
    split_BBID("D_div_wo_spg_island.dss", "D_spg_island.dss", "D_drn_wo_ro_island.dss",
               "D_RO_island.dss", outputfile, 2)  # DSM2 node daily output
    split_BBID("C_div_wo_spg_island.dss", "C_spg_island.dss",
               "C_drn_wo_ro_island.dss", "C_RO_island.dss", "delta_"+outputfile, 1)


def dcd_postprocess_out1(outputfile):
    """ Postprocess with the out1 extension option, DSM2 """
    islandtoDSM2node("D_div_wo_spg_island.dss", "D_spg_island.dss", "D_drn_wo_ro_island.dss",
                     "D_RO_island.dss", outputfile)  # DSM2 node daily output
    islandtoDSM2node("C_div_wo_spg_island.dss", "C_spg_island.dss",
                     "C_drn_wo_ro_island.dss", "C_RO_island.dss", "delta_"+outputfile)
    shutil.copy("delta_"+outputfile, "../../../Output/DSM2/")


def dcd_postprocess_out2(outputfile):
    """ Postprocess with the out2 extension option, SCHISM """
    islandtoDSM2node("D_div_wo_spg_island.dss", "D_spg_island.dss", "D_drn_wo_ro_island.dss",
                     "D_RO_island.dss", outputfile)  # DSM2 node daily output
    islandtoDSM2node("C_div_wo_spg_island.dss", "C_spg_island.dss",
                     "C_drn_wo_ro_island.dss", "C_RO_island.dss", "delta_"+outputfile)
    shutil.copy("delta_"+outputfile, "../../../Output/SCHISM/")


def dcd_postprocess_out3(outputfile, extension):
    """ Postprocess with the out3 extension option, CalSim3 """
    islandtoDSM2node("D_div_wo_spg_island.dss", "D_spg_island.dss", "D_drn_wo_ro_island.dss",
                     "D_RO_island.dss", outputfile)  # DSM2 node daily output
    daytomonth(outputfile)  # DSM2 node monthly output
    daytomonth("D_drn_wo_ro_island.dss")
    daytomonth("D_RO_island.dss")
    daytomonth("D_div_wo_spg_island.dss")
    daytomonth("D_spg_island.dss")
    DCD_to_CALSIM_ISLAND("D_div_wo_spg_island.dss", "D_spg_island.dss", "D_drn_wo_ro_island.dss",
                         "D_RO_island.dss", outputfile)  # combine island outputs to CS3 Delta grids

    islandtoDSM2node("C_div_wo_spg_island.dss", "C_spg_island.dss",
                     "C_drn_wo_ro_island.dss", "C_RO_island.dss", "delta_"+outputfile)
    daytomonth("delta_"+outputfile)
    daytomonth("C_drn_wo_ro_island.dss")
    daytomonth("C_RO_island.dss")
    daytomonth("C_div_wo_spg_island.dss")
    daytomonth("C_spg_island.dss")
    DCD_to_CALSIM_ISLAND("C_div_wo_spg_island.dss", "C_spg_island.dss", "C_drn_wo_ro_island.dss",
                         "C_RO_island.dss", "delta_" + outputfile)  # combine island outputs to CS3 Delta grids

    finalfiles = [outputfile.split(
        ".")[0]+"_mon.dss", outputfile.split(".")[0]+"_mon_C3.dss"]
    for i in range(len(finalfiles)):
        shutil.copy(finalfiles[i], "../../../Output/CALSIM3/")


def dcd_postprocess(outputfile, extension):
    if extension == "base":
        dcd_postprocess_base()
    if extension == "ex3":
        dcd_postprocess_ex3(outputfile)
    if extension == "out_1":
        dcd_postprocess_out1(outputfile)
    if extension == "out_2":
        dcd_postprocess_out2(outputfile)
    if extension == "out_3":
        dcd_postprocess_out3(outputfile, extension)


if __name__ == "__main__":
    pyhecdss.set_message_level(0)
    extension = sys.argv[2].strip()
    outputfile = sys.argv[1].strip()
    dcd_postprocess(outputfile, extension)
