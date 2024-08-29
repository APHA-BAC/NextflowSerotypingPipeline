#!/usr/bin/python3
# lims_rules.py by Josh
# Edits by: Arslan Hussaini

import sys
import csv
import re
import os
import argparse
import pandas as pd

consensusRegex = re.compile("^([123]-(.+?))(--1-(.+?)){0,1}(--1-(.+)){0,1}$")
subgenusRegex = re.compile("^[IVab]+ (.+)$")
anySpaceRegex = re.compile(" +")

def safe_int(someValue):
    try:
        someValue = int(someValue)
    except:
        pass
    return someValue

def safe_float(someValue):
    try:
        someValue = float(someValue)
    except:
        pass
    return someValue

def build_sero_dict():
    seroDict = {}
    with open(os.path.expanduser('~/summary/serogroup_lookup_table.tsv'), 'r') as seroLookup:
        reader = csv.reader(seroLookup, delimiter='\t')
        for idx, row in enumerate(reader):
            if idx == 0:
                continue
            serotype, seroDesc, serogroup = row[0], row[1], row[6]
            serotype = serotype.strip()
            seroDesc = seroDesc.strip()
            serotype = anySpaceRegex.sub("_", serotype)
            seroDesc = anySpaceRegex.sub("_", seroDesc)
            serotype = serotype.upper()

            if not seroDesc.upper() == serotype:
                seroDict[seroDesc.upper()] = [seroDesc, serogroup, "undetermined"]
            if "JAVA" in serotype:
                serotype = serotype.strip()
                seroDesc = seroDesc.strip()
                serotype = anySpaceRegex.sub("_", serotype)
                seroDesc = anySpaceRegex.sub("_", seroDesc)
                serotype = serotype.upper()
                seroDesc = seroDesc.replace("_", " ")
                seroDict[seroDesc.upper()] = [seroDesc, "B", "undetermined"]

            seroDict[serotype] = [seroDesc, serogroup, "undetermined"]

    with open(os.path.expanduser('~/summary/subgenus_lookup_table.tsv'), 'r') as subgLookup:
        reader = csv.reader(subgLookup, delimiter='\t')
        for idx, row in enumerate(reader):
            if idx == 0:
                continue
            serotype, subgenus = row[0:2]
            serotype = serotype.strip()
            serotype = anySpaceRegex.sub("_", serotype)
            serotype = serotype.upper()

            if serotype == "PARATYPHI_B_JAVA" or serotype == "PARATYPHI_B_VAR._JAVA" or serotype == "PARATYPHI_B_VARIANT_JAVA":

                serotype = "PARATYPHI B VAR. JAVA"

            if serotype in seroDict:
                seroDict[serotype][2] = subgenus
            else:
                seroDict[serotype] = [serotype, "undetermined", subgenus]
    return seroDict

def parse_consensus(consensus):
    if consensus == "no_result":
        serotypes = ["no_result"]
    else:
        if not consensusRegex.match(consensus):
            print(consensus)
            # sys.exit("Error, consensus not in recognisable format!")
            return False
        serotypes = [x for x in [consensusRegex.match(consensus).groups()[i] for i in (1, 3, 5)] if x is not None]
    return serotypes

def parse_seros(serotypes):
    serogroups = []
    subgenera = []
    limsSerotypes = []
    for serotype in serotypes:

        if serotype == "Monophasic Typhimurium":
            lookupSero = "Typhimurium"
        elif serotype in ["Java", "Paratyphi B var. L(+) tartrate+", "Paratyphi", "1-I 1,4,[5],12:b:", "1,4,[5],12:B:-"] or serotype in "1-Paratyphi--1-I 1,4,[5],12:b:---1-I 4:b:-":
            lookupSero = "PARATYPHI B VAR. JAVA"
        elif serotype in ["4:D:-", "I 1,4,[5],12:d:-"]:
            lookupSero = "4,12:D:-"
        elif subgenusRegex.match(serotype):
            lookupSero = subgenusRegex.match(serotype).group(1)
        else:
            lookupSero = serotype
        
        lookupSero = lookupSero.strip()
        lookupSero = lookupSero.replace("_", " ")
        lookupSero = lookupSero.upper()
        # print(lookupSero)
        if "JAVA" in lookupSero:
            lookupSero = "PARATYPHI B VAR. JAVA"
        if lookupSero not in seroDict:
            # print("Warning, serotype not in lookup dict:", lookupSero)
            serogroup = "undetermined"
            subgenus = "undetermined"
            limsSerotype = serotype
        else:
            limsSerotype = seroDict[lookupSero][0]
            serogroup = seroDict[lookupSero][1]
            subgenus = seroDict[lookupSero][2]
        serogroups.append(serogroup)
        subgenera.append(subgenus)
        if serotype == "Monophasic Typhimurium":
            limsSerotype = serotype

        limsSerotypes.append(limsSerotype)

    if len(serogroups) == 3:
        limsSerogroup = [x for x in serogroups if serogroups.count(x) > 1]
        if len(limsSerogroup) > 0:

            limsSerogroup = limsSerogroup[0]
        else:

            limsSerogroup = "undetermined"
    else:

        limsSerogroup = serogroups[0]
    
    if len(subgenera) == 3:
        limsSubgenus = [x for x in subgenera if subgenera.count(x) > 1]
        if len(limsSubgenus) > 0:
            limsSubgenus = limsSubgenus[0]

        else:
            limsSubgenus = "undetermined"
    else:
        limsSubgenus = subgenera[0]

    if serotype in ["JAVA", "Java", "Paratyphi B Variant Java"]:
        limsSubgenus = "I"
        limsSerogroup = "B"
    if "I 1,4,[5],12:d:-" in serotypes:
        limsSubgenus = "I"
        limsSerogroup = "B"
    if "I 1,4,[5],12:b:-" in serotypes:
        limsSubgenus = "I"
        limsSerogroup = "B"
    if "Paratyphi B var. Java" in serotypes:
        limsSubgenus = "I"
        limsSerogroup = "B"
    if "IIIb 61:k:1,5,(7)" in serotypes:
        limsSubgenus = "IIIb"
        limsSerogroup = "61"
    if "I G:i:-" in serotypes:
        limsSubgenus = "I"
        limsSerogroup = "G"

    return limsSerotypes, limsSerogroup, limsSubgenus


def parse_row(row):
    sampleID = row['Isolate_ID']
    consensus = row['Consensus']
    rawCount = row['#Reads_raw']
    readCount = row['#Reads_filtered']
    gc = row['GC%']
    kmerid = row['KmerID']
    contamFlag = row['Contam_Flag']
    most = row['MOST']
    mostLight = row['MOST_Light']
    st = row['MOST_ST']
    eBG = row["eBG"]
    mlst = row['MLST']
    mlstMeanCov = row['MLST_meanCov']
    seqsero = row['SeqSero']
    seqseroComment = row['SeqSero_comment']
    n50 = row['N50']
    serogr = row['sistr_Serogroup']
    serovar = row['sistr_Serovar']
    seroAnt = row['serovar_antigen']
    seroCGMLST = row['serovar_cgmlst']
    vaccine = row['vaccine']
    mono = row['mono']
    sseJ = row['sseJ']
    readRange = row['ReadLenRange']
    numContigs = row['#Contigs']
    numContigs25Kb = row['#Contigs>25Kbp']
    numContigs50Kb = row['#Contigs>50Kbp']
    assemblySize = row['AssemblySize']
    assemblyGC = row['AssemblyGC']
    L50 = row['L50']
    return sampleID, consensus, rawCount, readCount, gc, kmerid, contamFlag, most, mostLight, st, eBG, mlst, mlstMeanCov, seqsero, seqseroComment, n50, serogr, serovar, seroAnt, seroCGMLST, vaccine, mono, sseJ, readRange, numContigs, numContigs25Kb, numContigs50Kb, assemblySize, assemblyGC, L50

def apply_rules(limsSerotypes, limsSerogroup, limsSubgenus, row):
    sampleID, consensus, rawCount, readCount, gc, kmerid, contamFlag, most, mostLight, st, eBG, mlst, mlstMeanCov, seqsero, seqseroComment, n50, serogr, serovar, seroAnt, seroCGMLST, vaccine, mono, sseJ, readRange, numContigs, numContigs25Kb, numContigs50Kb, assemblySize, assemblyGC, L50 = parse_row(row)
    numReads = safe_int(rawCount)
    assemblySize = safe_int(assemblySize)
    copy_status = "Yes"
    try:
        n50.repalce(",","")
    except:
        pass
    n50 = safe_int(n50)
    numContigs = safe_int(numContigs)
    mlstMeanCov = safe_float(mlstMeanCov)
    kmerids = kmerid.split("--")
    salmPercent = [x.split("-")[0] for x in kmerids if "salmonella" in x]

    if len(salmPercent) == 1:
        salmPercent = safe_float(salmPercent[0])
    else:
        salmPercent = 0
    limsSerotypes = [subgenusRegex.match(x).group(1) if subgenusRegex.match(x) else x for x in limsSerotypes]
    if len(limsSerotypes) == 1:
        limsSerotype = limsSerotypes[0]
    elif len(limsSerotypes) == 2:
        limsSerotype = limsSerotypes[0]+" (2/3consensus)"
    else:
        limsSerotype = "no_consensus"
    limsReason = ""
    limsVariant = ""
    limsVaccine = ""
    limsStatus = ""

    # Paratyphi B rules
    if consensus == "2-Paratyphi--1-Paratyphi B var. Java":
        limsSerotype = "Paratyphi B Variant Java"
        limsVariant = "Variant Java"
        limsStatus = "Pass"
        limsReason = ""
    elif ("1-Paratyphi B var. Java" in consensus and "1-Java" in consensus and "1-Paratyphi B var. L(+) tartrate+" in consensus):
        limsSerotype = "Paratyphi B Variant Java"
        limsVariant = "Variant Java"
        limsStatus = "Pass"
    elif ("1-I 4:b:-" in consensus and "1-I 1,4,[5],12:b:-" in consensus and "1-Paratyphi" in consensus) and sseJ == "Java" and limsSubgenus == 'I' and salmPercent > 75:
        LIMS_SerotypeID = "Paratyphi B Variant Java"
        limsSerotype = "Paratyphi B Variant Java"
        limsVariant = "Monophasic Java"
        limsStatus = "Pass"
        limsSerogroup = "B"
    elif len([x for x in limsSerotypes if x in ("4:b:-", "1,4,[5],12:b:-", "Paratyphi", "Paratyphi B var. Java")]) == len(limsSerotypes) and sseJ == 'Java':
        limsSerotype = "Paratyphi B Variant Java"
        limsVariant = "Monophasic Java"
        limsStatus = "Pass"
    elif len([x for x in limsSerotypes if x in ("4:b:-", "1,4,[5],12:b:-", "Paratyphi", "Paratyphi B")]) == len(limsSerotypes) and sseJ == 'Paratyphi':
        limsSerotype = "Paratyphi B"
        limsStatus = "Pass"
    
    # Monophasic Typhi rule
    elif ("1-I 1,4,[5],12:i:-" in consensus and "2-Typhimurium" in consensus) or ("1-Typhimurium" in consensus and "1-I 1,4,[5],12:i:" in consensus and "1-I 4,[5],12:i:-" in consensus):
        limsSerotype = "Monophasic Typhimurium"
        limsStatus = "Pass"
    elif "Detected a deletion that causes O5- variant of Typhimurium" in seqseroComment:
        limsVariant = "O5 negative"
        
    # Idekan rule checked
    elif "3-Idikan" == consensus and mono == "MonophasicIdikan":
        limsVariant = "Monophasic Idikan"
        limsStatus = "Pass"
        limsReason = ""

    elif "3-Kedougou" == consensus and mono == "MonophasicKedougou" and salmPercent > 75:
        limsVariant = "Monophasic Kedougou"
        limsStatus = "Pass"
        limsReason = ""
    
    # Kedogou rules
    # elif ("1-I G:i:" in consensus or "1-I 13:i" in consensus or "1-Kedogou" in consensus) and mono == "MonophasicIdikan":
    #     limsSerotype = "Kedougou"
    #     limsVariant = "Monophasic Kedougou"
    #     limsStatus = "Pass"
    #     limsReason = ""
    elif ("1-I G:i:" in consensus and "1-I 13:i:" in consensus and "1-Kedougou" in consensus):
        limsSerotype = "Kedougou"
        limsStatus = "Pass"
        limsReason = ""
    
    
    #ARIZONAE IIIA 18:Z4:Z32
    elif consensus == "2-IIIa 44:z4,z23:---1-Arizonae":
        limsSerotype = "44:z4,z23:-"
        limsStatus == "CheckRequired"
    
    # RULE DIARIZONAE
    elif ("1-Arizonae" in consensus and "1-IIIb 61:k:1,5,(7)" in consensus and "1-IIIb O:61:k:1,5,7" in consensus) and salmPercent > 38:
        print("***********")
        print(consensus)
        print("Is it this one")
        print("********")
        limsSerotype = "61:k:1,5,7"
        limsStatus = "Pass"
    elif ("1-diarizonae" in consensus and "1-IIIb 61:k:1,5,(7)" in consensus and "1-IIIb O:61:k:1,5,7" in consensus) and salmPercent > 38:
        print("***********")
        print(consensus)
        print("or this one")
        print("********")
        limsSerotype = "61:k:1,5,7"
        limsStatus = "Pass"
    elif ("1-No Type" in consensus and "1-Arizonae" in consensus and "1-O:61:k:1,5,7" in consensus) and salmPercent > 38:
        print("***********")
        print(consensus)
        print("maybe this one")
        print("********")
        limsSerotype = "61:k:1,5,7"
        limsStatus = "Pass"

    # Clarify what this should be
    elif ("1-I 4:d:" in consensus and "1-Unnamed" in consensus and "1-I 1,4,[5],12:d:" in consensus):
        limsSerotype = "4,12:d:-"
        limsStatus = "Pass"
        limsReason = ""
    
    # Houtenae RULE
    elif limsSubgenus == "IV" and salmPercent > 38 and consensus == "2-IV 44:z4,z23:---1-Unnamed":
        limsSerotype = "Houtenae IV 44:z4,z23:-"
        limsStatus = "Pass"
    # RULE 17 HOUTENAE N/A
    elif consensus == "2-IV 44:z4,z23:---1-Unnamed":
        limsReason = "Check Serovar"
        limsSerotype = "Houtenae IV 44:z4,z23:-"
        limsStatus = "CheckRequired"
    
    # Uncommon enteritidis
    elif ("1-No Type--" in consensus and "1-Gallinarum or Enteritidis" in consensus and "1-Berta|Pensacola|Sangalkam" in consensus) or "2-No Type--1-Pensacola|Sangalkam" == consensus:
        limsSerotype = "Enteritidis"
        limsStatus = "CheckRequired"
        # limsReason = ""
    
    # RULE 22 4,12:D:-
    elif ("1-No Type" in consensus and "1-Unnamed" in consensus and "1-I 4,[5],12:d:-" in consensus):
        limsSerotype = "I 4,[5],12:d:-"
        limsStatus = "Pass"

    # RULE 24 BOVISMORBIFICANS
    elif "2-Bovismorbificans" in consensus and ("Bovis-Morbificans" in consensus or "Bovis-morbificans" in consensus):
        limsSerotype = "Bovismorbificans"
        limsStatus = "Pass"
        limsReason = ""
        copy_status = "Yes"
        consensus = "3-Bovismorbificans"

    # Goldcoast rule
    elif "2-Goldcoast" in consensus and "Gold-Coast" in consensus:
        limsSerotype = "Goldcoast"
        limsStatus = "Pass"
        limsReason= ""
        consensus = "3-Goldcoast"

    elif consensus == "2-IIIb 61:i:z53--1-No Type":
        limsStatus = "Pass"
        limsSerotype = "61:i:z53 "

    # 2-Pullorum vaccine RULE
    if vaccine == "2-Pullorum":
        limsSerotype = "Pullorum"
        limsVariant = "Variant Pullorum"
        limsVaccine = ""

    # Vaccine rules
    if vaccine not in ("NA", "srst2 result file not found"):
        if "2-AviproE" in vaccine or "8321_AviproE*_8322_AviproE*" in vaccine:
            limsVaccine = "Salmonella Enteritidis AVIPRO VAC E Vaccine Strain"
            limsStatus = "Pass"
            limsReason = ""
        elif "2-Salmovac440" in vaccine:
            limsVaccine = "Salmonella Enteritidis Salmovac 440 Vaccine Strain"
            if consensus == "2-Enteritidis--1-No Type" or consensus == "3-Enteritidis":
                limsStatus = "Pass"
                limsReason = ""
            else:
                limsStatus = "Inconclusive"
        elif "2-AviproT" in vaccine or vaccine == "527E_AviproT*" or "526_21b_AviproT*_527E_AviproT" in vaccine:
            limsVaccine = "Salmonella Typhimurium AVIPRO VAC T Vaccine Strain"
            limsStatus = "Pass"
            limsReason = ""
        elif "2-SalmoporcSTM" in vaccine:
            limsVaccine = "Salmonella Typhimurium Salmoporc STM Vaccine Strain"
            limsStatus = "Pass"
            limsReason = ""
        elif "2-Nobilis" in vaccine:
            limsVaccine = "Salmonella Gallinarum Noblis SG 9R Vaccine Strain"
            limsStatus = "Pass"
            limsReason = ""
        elif "2-Wild_Type" in vaccine:
            limsVaccine = "NOT VACCINE STRAIN"
        elif "Wild_Type" in vaccine:
            limsVaccine = "NOT VACCINE STRAIN"
        else:
            limsVaccine = vaccine
    elif "srst2 result file not found" in vaccine:
        limsVaccine = "Not typed by srst2"

    limsSerotype = limsSerotype.replace("(2/3consensus)", "")

    # if "3-No Type" in consensus:
    #     limsStatus = "Inconclusive"
    #     limsReason = "Contaminated: noIDedSerotypes"
    # if "2-No Type" in consensus and "e,h:-" in consensus:
    #     limsStatus = "Inconclusive"
    #     limsReason = "Contaminated: noIDedSerotypes"
    
    if limsSerotype == "no_consensus":
        limsStatus = "CheckRequired"
        limsReason = "Check Serovar"
    # Subgenus/KmerID % check
    if limsSubgenus == "I" and salmPercent > 75 and "3-" in consensus and "--" not in consensus:
        limsStatus = limsStatus.replace(" ", "")
        if len(limsStatus) <= 0:
            limsStatus = "Pass"
    elif limsSubgenus in ("II", 'IIIa', 'IIIb', 'IV', 'V') and salmPercent > 38 and "3-" in consensus and "--" not in consensus:
        if len(limsStatus) <= 0:
            limsStatus = "Pass"
    

    exceptions_list = ["Paratyphi B Variant Java","Paratyphi B","Monophasic Typhimurium",
    "Kedougou","4,12:d:-","S. enterica 4,12:d:-","Houtenae IV 44:z4,z23:-",
    "I 4,[5],12:d:-","Bovismorbificans", "61:k:1,5,7"]

    #  QUALITY CHECKS
    if numReads == "no_result" or isinstance(numReads, int) and numReads < 500000:
        limsReason = "InsufficientData: readCount<500K"
        copy_status = "No"
        limsStatus = "Inconclusive"
    elif (limsSubgenus == "I" or limsSubgenus == "undetermined") and salmPercent < 75:
        limsReason = "Contaminated: EntericaKmerID<75%"
        limsStatus = "Inconclusive"
        copy_status = "No"
    elif limsSubgenus in ("II", "IIIa", "IIIb", "IV", "V") and salmPercent < 38:
        limsReason = "Contaminated: non-EntericaKmerID<38%"
        limsStatus = "Inconclusive"
        copy_status = "No"
    elif isinstance(assemblySize, int) and assemblySize > 5800000:
        limsReason = "Contaminated: assembly>5.8Mbp"
        limsStatus = "Inconclusive"
        copy_status = "No"
    elif mostLight == "RED":
        limsReason = "PoorQuality: MOSTlightRED"
        limsStatus = "Inconclusive"
        copy_status = "No"
    elif isinstance(mlstMeanCov, float) and mlstMeanCov < 30:
        limsReason = "PoorQuality: MLSTcov<30x"
        limsStatus = "Inconclusive"
        copy_status = "No"
    elif st == "Failed(incomplete locus coverage)":
        limsReason = "PoorQuality: incomplSTcov(MOST)"
        limsStatus = "Inconclusive"
        copy_status = "No"
    elif isinstance(numContigs, int) and numContigs > 600:
        limsReason = "PoorQuality: contigCount>600"
        limsStatus = "Inconclusive"
        copy_status = "No"
    elif isinstance(n50, int) and n50 < 20000:
        limsReason = "PoorQuality: N50<20Kbp"
        limsStatus = "Inconclusive"
        copy_status = "No"
    elif isinstance(assemblySize, int) and assemblySize < 4000000:
        limsReason = "PoorQuality: assembly<4Mbp"
        limsStatus = "Inconclusive"
        copy_status = "No"
    elif "Co-existence of multiple serotypes detected" in seqseroComment and "3-" not in consensus:
        limsStatus = "Inconclusive"
        limsReason = "Contaminated: multipleSerotypesDetected(SeqSero2)"
        copy_status = "No"
    elif "No serotype antigens were detected. This is an atypical result that should be further investigated" in seqseroComment:
        limsStatus = "Inconclusive"
        limsReason = "PoorQuality: No serotype antigens detected by SeqSero"
        copy_status = "No"
    elif len([x for x in limsSerotypes if x in ('No Type', 'No Results')]) == len(limsSerotypes):
        limsReason = "Contaminated: noIDedSerotypes"
        limsStatus = "Inconclusive"
        copy_status = "No"
    if limsSubgenus not in ("II", "IIIa", "IIIb", "IV", "V") and salmPercent < 75 and limsReason != "InsufficientData: readCount<500K":
        limsStatus = "Inconclusive"
        limsReason = "Contaminated: EntericaKmerID<75%"
        

    if "3-" in consensus and ((salmPercent > 75 and (limsSubgenus == "I" or limsSubgenus == "undetermined")) or (salmPercent > 38 and limsSubgenus in ("II", 'IIIa', 'IIIb', 'IV', 'V'))):
        if "No Type" in consensus or "no_result" in consensus:
            pass 
        else:
            limsStatus = "Pass"
            limsReason = ""
    
    if limsSerotype.strip() in exceptions_list and "3-" not in consensus:
        limsStatus = "Pass"
        limsReason = ""

    elif "2-" in consensus and limsStatus == "":
        limsStatus = "CheckRequired"
        limsReason = "Check Serovar"
        copy_status = "No"

    return limsStatus, limsReason, limsSerotype, limsVariant, limsVaccine, copy_status

def parse_table(summaryTable):
    firstColNames = ['Isolate_ID', 'Consensus', 'LIMS_Status', 'LIMS_Reason', 'LIMS_SerotypeID', 'LIMS_Subgenus', 
    'LIMS_Serogroup', 'LIMS_Variant', 'LIMS_Vaccine', 'Copy_Status']
    outTable = []
    outFileName = os.path.basename(summaryTable).replace(".csv", "_plusLIMS.csv")
    df = pd.read_csv(summaryTable, keep_default_na=False)
    for idx, row in df.iterrows():
        
        if idx == 0:
            otherColNames = [x for x in list(df.columns) if x not in ('Isolate_ID', 'Consensus')]
            outTable.append(firstColNames + otherColNames)
        sampleID = row['Isolate_ID']
        consensus = row['Consensus']
        serotypes = parse_consensus(consensus)
        if not serotypes:
            print("Problem on this row: {}".format(idx))
            continue
        limsSerotypes, limsSerogroup, limsSubgenus = parse_seros(serotypes)

        if "Iiia" in consensus or "IIIa" in consensus:
            limsSubgenus = "IIIa"
        if "IIIb" in consensus:
            limsSubgenus = "IIIb"
        if "2-IV" in consensus:
            limsSubgenus = "IV"
        
        limsStatus, limsReason, limsSerotype, limsVariant, limsVaccine, copy_status = apply_rules(limsSerotypes, limsSerogroup, limsSubgenus, row)
        

        outRow = [sampleID, consensus, limsStatus, limsReason, limsSerotype, limsSubgenus, limsSerogroup, limsVariant, limsVaccine, copy_status] + list(row[otherColNames])
        outTable.append([str(x) for x in outRow])
    with open(outFileName, 'w') as outFile:
        writer = csv.writer(outFile, delimiter=',')
        writer.writerows(outTable)

if __name__ == '__main__':
    # Parse
    parser = argparse.ArgumentParser(description="Apply LIMS rule to pipeline summary table")
    parser.add_argument("summary_table", help="Pipeline summary CSV to be parsed")
    args = parser.parse_args()
    # Run
    seroDict = build_sero_dict()
    parse_table(args.summary_table)