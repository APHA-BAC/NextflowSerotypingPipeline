#!/usr/bin/env python3

# Originally by:    Javier Nunez, CSU, APHA
#                   Dolapo Ajayi, Bacterial Genomic Epidemiology, APHA
# Edits by:         Jaromir Guzinski, Bacterial Genomic Epidemiology, APHA
#                   Aaron Fishman, CSU, APHA
# Reworked by:      Joshua Potter, Core Bioinformatics, APHA, October 2021

# This script extracts results from respective pipeline tools and presents them in a summary table

import csv
import os
import sys
import re
import numpy
import fnmatch
import glob
import zipfile
import argparse
import pandas as pd

####################################################################################################

tableHeader = [
"Consensus",
"#Reads_raw",
"#Reads_filtered",
"GC%",
"KmerID",
"Contam_Flag",
"MOST",
"MOST_Light",
"MOST_ST",
"MLST",
"MLST_meanCov",
"SeqSero",
"SeqSero_comment",
"N50",
"sistr_Serogroup",
"sistr_Serovar",
"serovar_antigen",
"serovar_cgmlst",
"vaccine",
"mono",
"sseJ",
"ReadLenRange",
"#Contigs",
"#Contigs>25Kbp",
"#Contigs>50Kbp",
"AssemblySize",
"AssemblyGC",
"L50"
"EBG"
]

# tableHeader = ["StrainID", "Consensus", "#ReadsR1", "GC%R1", "R1Kmerid", "ContaminationFlag", "MOST", "Most_light", "st", "MLST", "MLST mean cov", "SeqSero", "SS comment", "N50", "serogroup", "serovar", "serovar_antigen", "serovar_cgmlst", "vaccine", "mono", "sseJ", "ReadLenRange", "#Contigs", "#Contigs>25Kbp", "#Contigs>50Kbp", "AssemblySize", "AssemblyGC", "L50"]

####################################################################################################

# def writeCSV(fname, matrix):
#     with open(fname, "w") as fileOut:
#         writer = csv.writer(fileOut)
#         writer.writerows(matrix)

def readTable(fname):
    infile = open(fname, "r")
    a = infile.readlines()
    infile.close()
    infile = open(fname, "r")
    if any(["\t" in x for x in a]):
        data = csv.reader(infile, delimiter="\t")
    else:
        data = csv.reader(infile, delimiter=",")
    dataOut = [row for row in data]
    infile.close()
    return(dataOut)

def get_subdirs(direc):
    subDirs = [os.path.abspath(x) for x in [os.path.join(direc, subpath) for subpath in os.listdir(direc)] if os.path.isdir(x)]
    return subDirs


####################################################################################################
# Get original readcount
def raw_count(sampleDir, sampleID):
    readcountDir = os.path.join(sampleDir, "readcount")
    readcountFile = os.path.join(readcountDir, "{}_readcount.txt".format(sampleID))
    if not os.path.exists(readcountFile) or os.path.getsize(readcountFile) == 0:
        return None
    with open(readcountFile) as inFile:
        rawCount = int(inFile.readline().strip())
    return rawCount


####################################################################################################
# Get QC metrics from FastQC

def fastqc_summary(sampleDir, sampleID):
    fastqcDir = os.path.join(sampleDir, "FASTQC_Reports")
    fastqcZip = [x for x in glob.glob(os.path.join(fastqcDir, "{}*_R1_fastqc.zip".format(sampleID)))]
    if len(fastqcZip) == 0 or os.path.getsize(fastqcZip[0]) == 0:
        return None, None, None, None
    with zipfile.ZipFile(fastqcZip[0], 'r') as archive:
        fastqcOutFile = [x for x in archive.namelist() if 'fastqc_data.txt' in x]
        if len(fastqcOutFile) == 0:
            return None, None, None, None
        with archive.open(fastqcOutFile[0]) as fastqcOut:
            for line in fastqcOut:
                line = line.decode('utf-8').strip()
                if "Total Sequences" in line:
                    readCount = line.split()[-1]
                if "Sequences flagged as poor quality" in line:
                    poorReads = line.split()[-1]
                if "Sequence length" in line:
                    readLen = line.split()[-1]
                if "%GC" in line:
                    gcR1 = line.split()[-1]
    return readCount, poorReads, readLen, gcR1


####################################################################################################
# kmerid analysis - can now account for variable tsv file names from kmerid

def kmerid_summary(sampleDir, sampleID):
    kmeridDir = os.path.join(sampleDir, "Kmerid")
    kmeridOut = [x for x in glob.glob(os.path.join(kmeridDir, "{}*_R1.tsv".format(sampleID)))]
    if len(kmeridOut) == 0 or os.path.getsize(kmeridOut[0]) == 0:
        return None, None
    kmeridOut = kmeridOut[0]
    kmeridResults = readTable(kmeridOut)
    kmeridResults = kmeridResults[2:]
    pathogens = [x[1] for x in kmeridResults]
    UniPathogens = list(set(pathogens))
    lista = []
    newkmeridResults = []
    for pathogen in UniPathogens:
        pathogenTable = [[round(float(x[0]), 2), x[1]] for x in kmeridResults if x[1] == pathogen]
        pathogenTable.sort(key=lambda x: -x[0])
        newkmeridResults = newkmeridResults+[pathogenTable[0]]
        newkmeridResults.sort(key=lambda x: -x[0])
    # Check if KmerID result file has any results by checking if array is filled
    if newkmeridResults:
        if newkmeridResults[1][0] > 20:
            kmeridFlag = "On"
        else:
            kmeridFlag = "Off"
        lista = [str(x[0])+"-"+x[1] for x in newkmeridResults]
        R1Kmerid = "--".join(lista)
    else:
        R1Kmerid = "NoResults"
        kmeridFlag = "NoResults"
    return R1Kmerid, kmeridFlag

####################################################################################################
# MOST analysis

def most_summary(sampleDir):
    mostDir = os.path.join(sampleDir, "MOST")
    mostFile = [x for x in glob.glob(os.path.join(mostDir, "*MLST_result.csv"))]
    lightFile = [x for x in glob.glob(os.path.join(mostDir, "*.results.xml"))]
    if len(lightFile) == 0 or os.path.getsize(lightFile[0]) == 0:
        light = None
    else:
        lightF = open(lightFile[0], "r")
        lightText = lightF.readlines()
        lightF.close()
        if any("AMBER" in x for x in lightText):
            light = "AMBER"
        if any("GREEN" in x for x in lightText):
            light = "GREEN"
        if any("RED" in x for x in lightText):
            light = "RED"
    if len(mostFile) == 0 or os.path.getsize(mostFile[0]) == 0:
        mostType = None
        st = None
        meanMLSTCov = None
        mlst = None
    else:
        mostFileName = mostFile[0]
        mostResults = readTable(mostFileName)
        mostType = [x for x in mostResults if "Predicted Serotype" in x][0][1].replace(
            "(", "").replace(")", "").replace("'", "").replace(",", "")  # .split(" ")
        st = [x for x in mostResults if "st value:" in x][0][1]
        AROC = [x for x in mostResults if "AROC" in x][0][1]
        DNAN = [x for x in mostResults if "DNAN" in x][0][1]
        HEMD = [x for x in mostResults if "HEMD" in x][0][1]
        HISD = [x for x in mostResults if "HISD" in x][0][1]
        PURE = [x for x in mostResults if "PURE" in x][0][1]
        SUCA = [x for x in mostResults if "SUCA" in x][0][1]
        THRA = [x for x in mostResults if "THRA" in x][0][1]
        AROCCov = [x for x in mostResults if "AROC" in x][0][5]
        DNANCov = [x for x in mostResults if "DNAN" in x][0][5]
        HEMDCov = [x for x in mostResults if "HEMD" in x][0][5]
        HISDCov = [x for x in mostResults if "HISD" in x][0][5]
        PURECov = [x for x in mostResults if "PURE" in x][0][5]
        SUCACov = [x for x in mostResults if "SUCA" in x][0][5]
        THRACov = [x for x in mostResults if "THRA" in x][0][5]
        # print([float(x) for x in [AROCCov, DNANCov, HEMDCov, HISDCov, PURECov, SUCACov, THRACov]])
        meanMLSTCov = round(numpy.mean(list(
            map(float, [AROCCov, DNANCov, HEMDCov, HISDCov, PURECov, SUCACov, THRACov]))), 2)
        mlst = ",".join([AROC, DNAN, HEMD, HISD, PURE, SUCA, THRA])
        if mostType != "no ST-serotype":
            mostType = mostType.split(" ")[0]
        else:
            mostType = "no type"
    return light, mostType, st, meanMLSTCov, mlst

####################################################################################################
# Using the most sequence type data, we can indentify the eBurst group for each sample

def ebgs(sampleDir):
    mostDir = os.path.join(sampleDir, "MOST")
    mostFile = [x for x in glob.glob(os.path.join(mostDir, "*MLST_result.csv"))]

    if len(mostFile) == 0 or os.path.getsize(mostFile[0]) == 0:
        ebg = None

    else:
        mostFileName = mostFile[0]
        mostResults = readTable(mostFileName)

        ebgFile = os.path.expanduser("~/summary/ebgs.csv")
        ebgData = readTable(ebgFile)

        st = [x for x in mostResults if "st value:" in x][0][1]
        st_str = str(st)
        st_str = st_str[1:]

        ebg = None
        for item in ebgData:
            if st_str in item:
                ebg = item[0]
    if ebg:
        return ebg
    else:
        return "No ebg"

####################################################################################################
# SeqSero2 analysis

def seqsero_summary(sampleDir):
    seqseroDir = os.path.join(sampleDir, "SeqSero2")
    seqseroFile = [x for x in glob.glob(os.path.join(seqseroDir, "SeqSero_result.txt"))]
    if len(seqseroFile) == 0 or os.path.getsize(seqseroFile[0]) == 0:
        seqseroType = None
        seqseroComment = None
    else:
        seqseroFileName = seqseroFile[0]
        seqseroFile = open(seqseroFileName, "r")
        seqseroResults = seqseroFile.readlines()
        seqseroFile.close()
        seqseroType = [x for x in seqseroResults if "Predicted serotype" in x][0].split("\t")[
            1].replace("\n", "")
        if "N/A" in seqseroType or "See comments below*" in seqseroType:
            seqseroType = "no type"
        seqseroComment = [x for x in seqseroResults if "Note:" in x][0].split("\t")[1].replace("\n", "")
    return seqseroType, seqseroComment

####################################################################################################
# quast analysis
    # 5/6/18 - Additional check - of whether N50 is empty or not

def quast_summary(sampleDir):
    quastDir = os.path.join(sampleDir, "quast")
    quastFile = [x for x in glob.glob(os.path.join(quastDir, "transposed_report.tsv"))]
    if len(quastFile) == 0 or os.path.getsize(quastFile[0]) == 0:
        contigs25k = contigs50k = contigs = assemLen = assemGC = N50 = L50 = None
    else:
        quastOut = readTable(quastFile[0])
        if len(quastOut) > 1:
            quastOut[0] = [x.strip() for x in quastOut[0]]
            contigs25k = quastOut[1][quastOut[0].index("# contigs (>= 25000 bp)")]
            contigs50k = quastOut[1][quastOut[0].index("# contigs (>= 50000 bp)")]
            contigs = quastOut[1][quastOut[0].index("# contigs")]
            assemLen = quastOut[1][quastOut[0].index("Total length")]
            assemGC = quastOut[1][quastOut[0].index("GC (%)")]
            N50 = quastOut[1][quastOut[0].index("N50")]
            L50 = quastOut[1][quastOut[0].index("L50")]
    return contigs25k, contigs50k, contigs, assemLen, assemGC, N50, L50


####################################################################################################
# sistr analysis

def sistr_summary(sampleDir):
    sistrDir = os.path.join(sampleDir, "sistr")
    sistrFile = [x for x in glob.glob(os.path.join(sistrDir, "sistr_prediction.csv"))]
    if len(sistrFile) == 0 or os.path.getsize(sistrFile[0]) == 0:
        serogroup = None
        serovar = None
        serovar_antigen = None
        serovar_cgmlst = None
    else:
        tab = readTable(sistrFile[0])
        serogroupind = tab[0].index("serogroup")
        serogroup = tab[1][serogroupind]

        serovarind = tab[0].index("serovar")
        serovar = tab[1][serovarind]

        serovar_antigenind = tab[0].index("serovar_antigen")
        serovar_antigen = tab[1][serovar_antigenind]

        serovar_cgmlstind = tab[0].index("serovar_cgmlst")
        serovar_cgmlst = tab[1][serovar_cgmlstind]
    return serogroup, serovar, serovar_antigen, serovar_cgmlst

####################################################################################################
# consensus calculation
# Additions - Reconiliation of seqsero and sistr results - 07/06/18
#           - Reconciliation of no type results - 08/06/18

def calc_consensus(seqseroType, mostType, serovar):
    # Set default consensus:
    clear_seqseroType = seqseroType
    clear_mostType = mostType
    clear_serovar = serovar
    # If the string that forms the MOST result is found in the SeqSero result and the string monophasic is not in
    # the SeqSero result then make the SeqSero and MOST result the same to show consensus, because MOST and SeqSero
    # express the same thing in different ways.
    if mostType.lower() in seqseroType.lower() and "monophasic" not in seqseroType:
        clear_seqseroType = mostType
    # To reconcile Seqsero's unassertive return of results with Sistr's assertive result e.g. Seqsero returning Senftenberg or Dessau* and Sistr returning Senftenberg.
    if serovar.lower() in seqseroType.lower():
        clear_serovar = clear_seqseroType = serovar
    # To account for the varying ways in which each of the programs express no type this line has been added.
    if "no " in mostType.lower():
        clear_mostType = "No Type"
    if seqseroType == "- -:-:-":
        clear_seqseroType = "No Type"
    if serovar == "-:-:-":
        clear_serovar = "No Type"
    # To account for the different ways in which SeqSero and Sistr express Monophasic Typhimurium and the fact that MOST
    # cannot report Monophasics but nevertheless reports as best as it can this line has been added to conclude that
    # the sample is Monophasic if the following conditions are met:
    if "typhimurium" == mostType.lower() and "I 4,[5],12:i:-" == seqseroType and "I 1,4,[5],12:i:-" == serovar:
        clear_mostType = clear_seqseroType = clear_serovar = "Monophasic Typhimurium"
    consensusTable = [clear_mostType, clear_seqseroType, clear_serovar]
    uniqueTypes = list(set(consensusTable))
    consensus = []
    for sero in uniqueTypes:
        count = consensusTable.count(sero)
        consensus.append([count, sero])
    consensus.sort(key=lambda x: -x[0])
    lista = [str(x[0])+"-"+x[1] for x in consensus]
    consensus = "--".join(lista)
    if "Enteritidis" in mostType and "See comments" in seqseroType:
        consensus = "Enteritidis"
    return consensus

####################################################################################################
# vaccine extraction

def vaccine_diff(sampleDir, serovar):
    if serovar.lower() in ["enteritidis", "typhimurium"]:
        srst2Dir = os.path.join(sampleDir, "srst2")
        srst2File = [x for x in glob.glob(os.path.join(srst2Dir, "__genes__*results.txt"))]
        if len(srst2File) == 0 or os.path.getsize(srst2File[0]) == 0:
            vaccine = "srst2 result file not found"
        else:
            tab = readTable(srst2File[0])[1][1:]
            tab = [x for x in tab if "WT" not in x and "not" not in x and "Nobilis" not in x and "sseJ" not in x]
            if len(tab) > 0:
                # print(tab)
                vaccine = "_".join(tab)
                if vaccine.count('Salmoporc') == 2:
                    vaccine = '2-SalmoporcSTM'
                elif vaccine.count('Salmoporc') == 1:
                    vaccine = '1-SalmoporcSTM--1-Wild_Type'
                elif vaccine.count('AviproE') == 2:
                    vaccine = '2-AviproE'
                elif vaccine.count('AviproE') == 1:
                    vaccine.count('1-AviproE--1-Wild_Type')
                elif vaccine.count('AviproT') == 2:
                    vaccine = '2-AviproT'
                elif vaccine.count('AviproT') == 1:
                    vaccine.count('1-AviproT--1-Wild_Type')
                elif vaccine.count('Gallivac') == 2:
                    vaccine = '2-Salmovac440'           # Previously known as Gallivac
                elif vaccine.count('Gallivac') == 1:
                    vaccine = '1-Salmovac440--1-Wild_Type'
            else:
                vaccine = "2-Wild_Type"
            if len(tab) == 0:
                vaccine = "2-Wild_Type"
    elif serovar.lower() in ["gallinarum", "pullorum"]:
        srst2Dir = os.path.join(sampleDir, "srst2")
        srst2File = [x for x in glob.glob(os.path.join(srst2Dir, "__genes__*results.txt"))]
        if len(srst2File) == 0 or os.path.getsize(srst2File[0]) == 0:
            vaccine = "srst2 result file not found"
        else:
            tab = readTable(srst2File[0])[1][1:]
            tab = [x for x in tab if "WT" not in x and "not" not in x and "sseJ" not in x] #and "like" not in x
            if len(tab) > 0:
                # print(tab)
                vaccine = "_".join(tab)
                if vaccine.count('Pullorum') == 2:
                    vaccine = '2-Pullorum'
                elif vaccine.count('Pullorum') == 1:
                    vaccine = '1-Pullorum--1-Wild_Type'
                elif vaccine.count('Gallinarum') == 2 and vaccine.count('Nobilis') == 0:
                    vaccine = '2-Gallinarum'
                elif vaccine.count('Gallinarum') == 2 and vaccine.count('Nobilis') == 1:
                    vaccine = '2-Gallinarum'
                elif vaccine.count('Gallinarum') == 2 and vaccine.count('Nobilis_A') == 1 and vaccine.count('Nobilis_B') == 1 and vaccine.count('Nobilis_C') == 1:
                    vaccine = 'Nobilis'
                elif vaccine.count('Gallinarum') == 2 and vaccine.count('Nob_C_like') == 1:
                    vaccine = 'Nobilis-like'
                elif vaccine.count('Gallinarum') == 2 and vaccine.count('Nobilis_A') == 1 and vaccine.count('Nobilis_B') == 1:
                    vaccine = 'NobilisA_NobillisB'
                elif vaccine.count('Gallinarum') == 2 and vaccine.count('Nobilis_A') == 1 and vaccine.count('Nobilis_C') == 1:
                    vaccine = 'NobilisA_NobillisC'
                elif vaccine.count('Gallinarum') == 2 and vaccine.count('Nobilis_B') == 1 and vaccine.count('Nobilis_C') == 1:
                    vaccine = 'NobilisB_NobilisC'
            else:
                vaccine = "2-Wild_Type"
            if len(tab) == 0:
                vaccine = "2-Wild_Type"
    else:
        vaccine = "NA"
    return vaccine


####################################################################################################
# 13,23i differentiation

def thirteen23i_diff(sampleDir, mostType):
    if mostType.lower() in ["idikan","kedougou"]:
        srst2Dir = os.path.join(sampleDir, "srst2")
        srst2File = [x for x in glob.glob(os.path.join(srst2Dir, "__genes__*results.txt"))]
        if len(srst2File) == 0 or os.path.getsize(srst2File[0]) == 0:
            mono = "srst2 results file not found"
        else:
            tab = readTable(srst2File[0])[1][1:]
            tab = [x for x in tab if "WT" not in x and "Not" not in x and "Nobilis_C" not in x and "sseJ" not in x]
            if len(tab) > 0:
                # print(tab)
                mono = "_".join(tab)
                if mono.count('monoIdikanA') == 1 and mono.count('monoIdikanB') == 1:
                    mono = 'MonophasicIdikan'
                elif mono.count('monoKedougouA') == 1 and mono.count('monoKedougouB') == 1:
                    mono = 'MonophasicKedougou'
                elif mono.count('monoIdikanA') == 0 and mono.count('monoIdikanB') == 0:
                    mono = 'NA'
                elif mono.count('monoKedougouA') == 0 and mono.count('monoKedougouB') == 0:
                    mono = 'NA'
            else:
                mono = "NA"
    else:
        mono = "NA"
    return mono

####################################################################################################
# ParatyphiB/Java differentiation

def paratyphiB_java_diff(sampleDir, mostType):
    if mostType.lower() in ["java","paratyphi"]:
        srst2Dir = os.path.join(sampleDir, "srst2")
        srst2File = [x for x in glob.glob(os.path.join(srst2Dir, "__genes__*results.txt"))]
        if len(srst2File) == 0 or os.path.getsize(srst2File[0]) == 0:
            sseJ = "srst2 results file not found"
        else:
            tab = readTable(srst2File[0])[1][1:]
            tab = [x for x in tab if "WT" not in x and "Not" not in x and "Nob" not in x and "Nobilis" not in x]
            if len(tab) > 0:
                # print(tab)
                sseJ = "_".join(tab)
                if sseJ.count('sseJ') == 1:
                    sseJ = 'Java'
                elif sseJ.count('sseJ') == 0:
                    sseJ = 'ParatyphiB'
            else:
                sseJ = "NA"
    else:
        sseJ = "NA"
    return sseJ

####################################################################################################

def instantiate_summary(resultsDir, runID):
    sampleDirs = get_subdirs(resultsDir)
    sampleNames = [os.path.basename(x) for x in sampleDirs]
    df = pd.DataFrame(columns = tableHeader, index = sampleNames)
    df.index.rename('Isolate_ID', inplace=True)
    df = df.fillna('no_result')
    for sampleDir in sampleDirs:
        sampleID = os.path.basename(sampleDir)
        rawCount = raw_count(sampleDir, sampleID)
        if rawCount:
            df.loc[sampleID, "#Reads_raw"] = rawCount
            df.loc[sampleID, "SeqSero_comment"] = ""
    summaryFileName = os.path.join(resultsDir, runID + "_SummaryTable.csv")
    print(df)
    df.to_csv(summaryFileName)

def fill_summary(resultsDir, runID):
    summaryFileName = os.path.join(resultsDir, runID + "_SummaryTable.csv")
    df = pd.read_csv(summaryFileName, index_col=0)
    sampleDirs = get_subdirs(resultsDir)
    for sampleDir in sampleDirs:
        sampleID = os.path.basename(sampleDir)
        readCount, poorReads, readLen, gcR1 = fastqc_summary(sampleDir, sampleID)
        if readCount:
            df.loc[sampleID, "#Reads_filtered"] = readCount
        if readLen:
            df.loc[sampleID, "ReadLenRange"] = readLen
        if gcR1:
            df.loc[sampleID, "GC%"] = gcR1
        R1Kmerid, kmeridFlag = kmerid_summary(sampleDir, sampleID)
        if R1Kmerid:
            df.loc[sampleID, "KmerID"] = R1Kmerid
        if kmeridFlag:
            df.loc[sampleID, "Contam_Flag"] = kmeridFlag
        light, mostType, st, meanMLSTCov, mlst = most_summary(sampleDir)
        if light:
            df.loc[sampleID, "MOST_Light"] = light
        if mostType:
            df.loc[sampleID, "MOST"] = mostType
        if st:
            df.loc[sampleID, "MOST_ST"] = st
        if meanMLSTCov:
            df.loc[sampleID, "MLST_meanCov"] = meanMLSTCov
        if mlst:
            df.loc[sampleID, "MLST"] = mlst
        seqseroType, seqseroComment = seqsero_summary(sampleDir)
        if seqseroType:
            df.loc[sampleID, "SeqSero"] = seqseroType
        if seqseroComment:
            df.loc[sampleID, "SeqSero_comment"] = seqseroComment
        else:
            df.loc[sampleID, "SeqSero_comment"] = ""
        contigs25k, contigs50k, contigs, assemLen, assemGC, N50, L50 = quast_summary(sampleDir)
        if contigs25k:
            df.loc[sampleID, "#Contigs>25Kbp"] = contigs25k
        if contigs50k:
            df.loc[sampleID, "#Contigs>50Kbp"] = contigs50k
        if contigs:
            df.loc[sampleID, "#Contigs"] = contigs
        if assemLen:
            df.loc[sampleID, "AssemblySize"] = assemLen
        if assemGC:
            df.loc[sampleID, "AssemblyGC"] = assemGC
        if N50:
            df.loc[sampleID, "N50"] = N50
        if L50:
            df.loc[sampleID, "L50"] = L50
        serogroup, serovar, serovar_antigen, serovar_cgmlst = sistr_summary(sampleDir)
        if serogroup:
            df.loc[sampleID, "sistr_Serogroup"] = serogroup
        if serovar:
            df.loc[sampleID, "sistr_Serovar"] = serovar
        if serovar_antigen:
            df.loc[sampleID, "serovar_antigen"] = serovar_antigen
        if serovar_cgmlst:
            df.loc[sampleID, "serovar_cgmlst"] = serovar_cgmlst
        if not seqseroType:
            seqseroType = "no_result"
        if not mostType:
            mostType = "no_result"
        if not serovar:
            serovar = "no_result"
        consensus = calc_consensus(seqseroType, mostType, serovar)
        if consensus:
            df.loc[sampleID, "Consensus"] = consensus
        vaccine = vaccine_diff(sampleDir, serovar)
        if vaccine:
            df.loc[sampleID, "vaccine"] = vaccine
        mono = thirteen23i_diff(sampleDir, mostType)
        if mono:
            df.loc[sampleID, "mono"] = mono
        sseJ = paratyphiB_java_diff(sampleDir, mostType)
        if sseJ:
            df.loc[sampleID, "sseJ"] = sseJ
        ebg = ebgs(sampleDir)
        if ebg:

            df.loc[sampleID, "EBG"] = ebg

    print(df)
    df.to_csv(summaryFileName)


# Command line args
def main():
    parser = argparse.ArgumentParser(description="Create summary table from Salmonella pipeline run.")
    parser.add_argument("runID", help="The run ID of the specific pipeline output to be analysed, e.g. 1234 would be the run ID for ~/WGS_Results/1234")
    parser.add_argument("--instantiate", default=False, action='store_true')
    args = parser.parse_args()
    return args

if __name__ == '__main__':
    args = main()
    runID = args.runID
    print(runID)
    resultsDir = os.path.join(os.path.expanduser("~/WGS_Results"), runID)
    print(resultsDir)
    if args.instantiate:
        print("Instantiating summary table...")
        instantiate_summary(resultsDir, runID)
    else:
        fill_summary(resultsDir, runID)

quit()

# del sys.argv[0]

# if not len(sys.argv) == 1:
#     sys.exit("Error, no run name specified")


# runID = sys.argv[0]
# print(runID)
# resultsDir = os.path.join(os.path.expanduser("~/WGS_Results"), runID)
# print(resultsDir)
# outTable = [tableHeader]
# sampleDirs = get_subdirs(resultsDir)

# for sampleDir in sampleDirs:
#     sampleID = os.path.basename(sampleDir)
#     print(sampleID)
#     readCount, poorReads, readLen, gcR1 = fastqc_summary(sampleDir)
#     R1Kmerid, kmeridFlag = kmerid_summary(sampleDir)
#     light, mostType, st, meanMLSTCov, mlst = most_summary(sampleDir)
#     seqseroType, seqseroComment = seqsero_summary(sampleDir)
#     contigs25k, contigs50k, contigs, assemLen, assemGC, N50, L50 = quast_summary(sampleDir)
#     serogroup, serovar, serovar_antigen, serovar_cgmlst = sistr_summary(sampleDir)
#     consensus = calc_consensus(seqseroType, mostType, serovar)
#     vaccine = vaccine_diff(sampleDir, serovar)
#     mono = thirteen23i_diff(sampleDir, mostType)
#     sseJ = paratyphiB_java_diff(sampleDir, mostType)
#     outTable.append([sampleID, consensus, readCount, gcR1, R1Kmerid, kmeridFlag, mostType, light, st, mlst, meanMLSTCov,
#         seqseroType, seqseroComment, N50, serogroup, serovar, serovar_antigen, serovar_cgmlst,
#         vaccine, mono, sseJ, readLen, contigs, contigs25k, contigs50k, assemLen, assemGC, L50])

# writeCSV(os.path.join(resultsDir, runID+"_SummaryTable.csv"), outTable)

# def filterDirs(folder):
#     return [d for d in (os.path.join(folder, d1) for d1 in os.listdir(folder)) if os.path.isdir(d)]

# def nreads(fName):
#     os.system("unpigz -stdout "+fName+" > "+fName[:-3])
#     os.system("wc -l "+fName[:-3]+" > "+fName[:-3]+".wcl.txt")
#     f = open(fName[:-3]+".wcl.txt", "rb")
#     val = f.readlines()[0].split(" ")[0]
#     print(val)
#     val = int(val)/4
#     f.close()
#     os.system("rm "+fName[:-3]+".wcl.txt")
#     os.system("rm "+fName[:-3])
#     return val

# def GCcontent(fName):
#     os.system("unpigz --stdout "+fName+" > "+fName[:-3])
#     os.system("wc -l "+fName[:-3]+" > "+fName[:-3]+".wcl.txt")
#     f = open(fName[:-3]+".wcl.txt", "r")
#     val = int(f.readlines()[0].split(" ")[0])/4
#     print(val)
#     f.close()
#     fileIn = open(fName[:-3], 'r')
#     As = 0
#     Ts = 0
#     Cs = 0
#     Gs = 0
#     allb = 0
#     line = fileIn.readline()
#     cont2 = 0
#     while line and cont2 < 200000:
#         if line[0] == "@":
#             cont2 = cont2+1
#             seq = fileIn.readline()
#             As = As+seq.count("A")+seq.count("a")
#             Ts = Ts+seq.count("T")+seq.count("t")
#             Cs = Cs+seq.count("C")+seq.count("c")
#             Gs = Gs+seq.count("G")+seq.count("g")
#             allb = allb+len(seq[:-1])
#             fileIn.readline()  # +
#             fileIn.readline()  # Qualities
#             line = fileIn.readline()
#     fileIn.close()
#     os.system("rm "+fName[:-3])
#     os.system("rm "+fName[:-3]+".wcl.txt")
#     if allb == 0:
#         pass
#     else:
#         return val, round(100*float(Cs+Gs)/allb, 2)

# def find(pattern, path):
#     result = []
#     for root, dirs, files in os.walk(path):
#         for name in files:
#             if fnmatch.fnmatch(name, pattern):
#                 result.append(os.path.join(root, name))
#     return result
