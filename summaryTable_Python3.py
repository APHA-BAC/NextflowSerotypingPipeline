#!/usr/bin/env python

# Javier Nunez, CSU, APHA
# This script extraction the results from respective pipeline programs and presents them in a summary table

import csv
import os
import os.path
import sys
import re
import numpy
import fnmatch
import glob


def writeCSV(fname, matrix):
    with open(fname, "w") as fileOut:
        writer = csv.writer(fileOut)
        writer.writerows(matrix)


def readTable(fname, ch):
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


def find(pattern, path):
    result = []
    for root, dirs, files in os.walk(path):
        for name in files:
            if fnmatch.fnmatch(name, pattern):
                result.append(os.path.join(root, name))
    return result


def filterDirs(folder):
    return [d for d in (os.path.join(folder, d1) for d1 in os.listdir(folder)) if os.path.isdir(d)]


def nreads(fName):
    os.system("unpigz -stdout "+fName+" > "+fName[:-3])
    os.system("wc -l "+fName[:-3]+" > "+fName[:-3]+".wcl.txt")
    f = open(fName[:-3]+".wcl.txt", "rb")
    val = f.readlines()[0].split(" ")[0]
    print(val)
    val = int(val)/4
    f.close()
    os.system("rm "+fName[:-3]+".wcl.txt")
    os.system("rm "+fName[:-3])
    return val


def GCcontent(fName):
    os.system("unpigz --stdout "+fName+" > "+fName[:-3])
    os.system("wc -l "+fName[:-3]+" > "+fName[:-3]+".wcl.txt")
    f = open(fName[:-3]+".wcl.txt", "r")
    val = int(f.readlines()[0].split(" ")[0])/4
    print(val)
    f.close()
    fileIn = open(fName[:-3], 'r')
    As = 0
    Ts = 0
    Cs = 0
    Gs = 0
    allb = 0
    line = fileIn.readline()
    cont2 = 0
    while line and cont2 < 200000:
        if line[0] == "@":
            cont2 = cont2+1
            seq = fileIn.readline()
            As = As+seq.count("A")+seq.count("a")
            Ts = Ts+seq.count("T")+seq.count("t")
            Cs = Cs+seq.count("C")+seq.count("c")
            Gs = Gs+seq.count("G")+seq.count("g")
            allb = allb+len(seq[:-1])
            fileIn.readline()  # +
            fileIn.readline()  # Qualities
            line = fileIn.readline()
    fileIn.close()
    os.system("rm "+fName[:-3])
    os.system("rm "+fName[:-3]+".wcl.txt")
    if allb == 0:
        pass
    else:
        return val, round(100*float(Cs+Gs)/allb, 2)


##########################
##########################
##########################


#### HACK: passing directory in as input ###

input_dir = sys.argv[1]

###

pathoData = glob.glob(input_dir + "*")
pathoData = str(pathoData[0])
print(pathoData)

pathoResults = pathoData.replace("WGS_Data", "WGS_Results")
runName = pathoData.split(os.sep)[-1]
dirs = filterDirs(pathoResults)


#pathoData = glob.glob("/home/*/WGS_Data/*")
#pathoData = str(pathoData[0])
#pathoResults = pathoData.replace("WGS_Data", "WGS_Results")
#runName = pathoData.split(os.sep)[-1]
#dirs = filterDirs(pathoResults)
#print (input_dir)
#print(pathoData)
#print(pathoResults)
#print(dirs)


tabo = [["StrainID", "Consensus", "#ReadsR1", "GC%R1", "R1Kmerid", "ContaminationFlag", "MOST", "Most_light", "st", "MLST",
         "MLST mean cov", "SeqSero", "SS comment", "N50", "serogroup", "serovar", "serovar_antigen", "serovar_cgmlst", "vaccine", "mono", "sseJ"]]
for diro in dirs:
    strainName = diro.split(os.sep)[-1]
    print("Summarizing sample "+strainName)
    fastqgz = sorted(find(strainName+"*.fastq.gz", pathoData))
    if len(fastqgz) > 0:
        readsR1, GCR1 = GCcontent(fastqgz[0])
    else:
        readsR1 = "NaN"
        GCR1 = "NaN"

    ### kmerid analysis - can now account for variable tsv file names from kmerid  ###
    kmeridR1File = [f for f in os.listdir(os.path.join(
        diro, "Kmerid")) if strainName in f and "_R1.tsv" in f]
    if len(kmeridR1File) > 0:
        kmeridFileName = kmeridR1File[0]
        print(kmeridFileName)
        kmeridResults = readTable(os.path.join(
            diro, "Kmerid", kmeridFileName), "\t")[2:]
        pathogens = [x[1] for x in kmeridResults]
        UniPathogens = list(set(pathogens))
        lista = []
        print(UniPathogens)
        newkmeridResults = []
        for pathogen in UniPathogens:
            pathogenTable = [[round(float(x[0]), 2), x[1]]
                             for x in kmeridResults if x[1] == pathogen]
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
    else:
        R1Kmerid = "NoResults"
        kmeridFlag = "NoResults"

    # MOST analysis
    mostFile = find("*MLST_result.csv", diro)
    lightFile = find("*.results.xml", diro)
    light = "NA"
    if len(lightFile) > 0:
        lightF = open(lightFile[0], "r")
        lightText = lightF.readlines()
        lightF.close()
        if any("AMBER" in x for x in lightText):
            light = "AMBER"
        if any("GREEN" in x for x in lightText):
            light = "GREEN"
        if any("RED" in x for x in lightText):
            light = "RED"

    if len(mostFile) > 0:
        mostFileName = mostFile[0]
        mostResults = readTable(mostFileName, ",")
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
        print(list(map(float, [AROCCov, DNANCov, HEMDCov,
                               HISDCov, PURECov, SUCACov, THRACov])))
        meanMLSTCov = round(numpy.mean(list(
            map(float, [AROCCov, DNANCov, HEMDCov, HISDCov, PURECov, SUCACov, THRACov]))), 2)
        mlst = ",".join([AROC, DNAN, HEMD, HISD, PURE, SUCA, THRA])
        if mostType != "no ST-serotype":
            mostType = mostType.split(" ")[0]
        else:
            mostType = "no type"
    else:
        mostType = "NoResults"
        st = "NoResults"
        meanMLSTCov = "NoResults"
        mlst = "NoResults"
    print(mostType)

    # Seqsero analysis
    seqseroFile = find("*Sero_result.txt", diro)
    print(seqseroFile)
    if len(seqseroFile) > 0:
        seqseroFileName = seqseroFile[0]
        seqseroFile = open(seqseroFileName, "r")
        seqseroResults = seqseroFile.readlines()
        seqseroFile.close()
        seqseroType = [x for x in seqseroResults if "Predicted serotype" in x][0].split("\t")[
            1].replace("\n", "")
        if "N/A" in seqseroType or "See comments below*" in seqseroType:
            seqseroType = "no type"
    else:
        seqseroType = "NoResults"
    print(seqseroType)

    comment = [x for x in seqseroResults if "Note:" in x][0].split("\t")[
        1].replace("\n", "")

    # quast analysis
    # 5/6/18 - Additional check - of whether N50 is empty or not
    quatFile = find("*transposed_report.tsv", diro)
    print(quatFile)
    if len(quatFile) > 0:
        N50 = readTable(quatFile[0], "\t")
        if len(N50) > 0:
            N50ind = N50[0].index("N50")
            N50 = N50[1][N50ind]
    else:
        N50ind = "NoQuast"
        N50 = "NoQuast"

    # sistr analysis
    sistrFile = find("*sistr_prediction.csv", diro)
    if len(sistrFile) > 0:
        tab = readTable(sistrFile[0], ",")
        serogroupind = tab[0].index("serogroup")
        serogroup = tab[1][serogroupind]

        serovarind = tab[0].index("serovar")
        serovar = tab[1][serovarind]

        serovar_antigenind = tab[0].index("serovar_antigen")
        serovar_antigen = tab[1][serovar_antigenind]

        serovar_cgmlstind = tab[0].index("serovar_cgmlst")
        serovar_cgmlst = tab[1][serovar_cgmlstind]
    else:
        serogroup = "NoResults"
        serovar = "NoResults"
        serovar_antigen = "NoResults"
        serovar_cgmlst = "NoResults"

    # consensus calculation
    # Additions - Reconiliation of seqsero and sistr results - 07/06/18
    #           - Reconciliation of no type results - 08/06/18

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

    # new rules to added on 21/10/2020

    # vaccine extraction
    #if mostType.lower() in ["enteritidis","typhimurium"]:
    if serovar.lower() in ["enteritidis", "typhimurium"]:
        srst2File = find("*__genes__*", diro)
        if len(srst2File) > 0:
            tab = readTable(srst2File[0], "/t")[1][1:]
            tab = [x for x in tab if "WT" not in x and "not" not in x and "Nobilis" not in x and "sseJ" not in x]
            if len(tab) > 0:
                print(tab)
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

        else:
            vaccine = "srst2 file result not found"


    #elif mostType.lower() in ["gallinarum","pullorum"]:
    elif serovar.lower() in ["gallinarum", "pullorum"]:
        srst2File = find("*__genes__*", diro)
        if len(srst2File) > 0:
            tab = readTable(srst2File[0], "/t")[1][1:]
            tab = [x for x in tab if "WT" not in x and "not" not in x and "sseJ" not in x] #and "like" not in x
            if len(tab) > 0:
                print(tab)
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
            vaccine = "srst2 file result not found"
    else:
        vaccine = "NA"


    if "Enteritidis" in mostType and "See comments" in seqseroType:
        consensus = "Enteritidis"


    # 13,23i differentiation
    if mostType.lower() in ["idikan","kedougou"]:
    #if serovar.lower() in ["idikan","kedougou"]:
        srst2File = find("*__genes__*", diro)
        if len(srst2File) > 0:
            tab = readTable(srst2File[0], "/t")[1][1:]
            tab = [x for x in tab if "WT" not in x and "Not" not in x and "Nobilis_C" not in x and "sseJ" not in x]
            if len(tab) > 0:
                print(tab)
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
    else:
        mono = "NA"


    # ParatyphiB/Java differentiation
    if mostType.lower() in ["java","paratyphi"]:
    #if serovar.lower() in ["java","paratyphi"]:
        srst2File = find("*__genes__*", diro)
        if len(srst2File) > 0:
            tab = readTable(srst2File[0], "/t")[1][1:]
            tab = [x for x in tab if "WT" not in x and "Not" not in x and "Nob" not in x and "Nobilis" not in x]
            if len(tab) > 0:
                print(tab)
                sseJ = "_".join(tab)

                if sseJ.count('sseJ') == 1:
                    sseJ = 'Java'
                elif sseJ.count('sseJ') == 0:
                    sseJ = 'ParatyphiB'

            else:
                sseJ = "NA"
        else:
            sseJ = "NA"
    else:
        sseJ = "NA"

    tabo = tabo+[[strainName, consensus, readsR1, GCR1, R1Kmerid, kmeridFlag, mostType, light, st, mlst,
                  meanMLSTCov, seqseroType, comment, N50, serogroup, serovar, serovar_antigen, serovar_cgmlst, vaccine, mono, sseJ]]
    writeCSV(os.path.join(pathoResults, runName+"_SummaryTable.csv"), tabo)
