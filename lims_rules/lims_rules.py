#!/usr/bin/python3
# lims_rules.py by Josh

import sys
import csv
import re
import os
import argparse

headerRow = ["Isolate_ID",
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
]

outTableHeader = ["Isolate_ID",
"Consensus", 
"LIMS_Status", "LIMS_Reason", "LIMS_SerotypeID", "LIMS_Subgenus", "LIMS_Serogroup", "LIMS_Variant", "LIMS_Vaccine",
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
]

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
                print("\t".join([serotype, seroDesc]))
                seroDict[seroDesc.upper()] = [seroDesc, serogroup, "undetermined"]
            seroDict[serotype] = [seroDesc, serogroup, "undetermined"]
    print("\n\n")
    with open(os.path.expanduser('~/summary/subgenus_lookup_table.tsv'), 'r') as subgLookup:
        reader = csv.reader(subgLookup, delimiter='\t')
        for idx, row in enumerate(reader):
            if idx == 0:
                continue
            serotype, subgenus = row[0:2]
            serotype = serotype.strip()
            serotype = anySpaceRegex.sub("_", serotype)
            serotype = serotype.upper()
            if serotype == "PARATYPHI_B_JAVA":
                serotype = "PARATYPHI_B_VAR._JAVA"
            if serotype in seroDict:
                seroDict[serotype][2] = subgenus
            else:
                seroDict[serotype] = [serotype, "undetermined", subgenus]
    return seroDict

def parse_consensus(consensus):
    if not consensusRegex.match(consensus):
        print(consensus)
        sys.exit("Error, consensus not in recognisable format!")
    serotypes = [x for x in [consensusRegex.match(consensus).groups()[i] for i in (1, 3, 5)] if x is not None]
    return serotypes

def parse_seros(serotypes):
    serogroups = []
    subgenera = []
    limsSerotypes = []
    for serotype in serotypes:
        if serotype == "Monophasic Typhimurium":
            lookupSero = "Typhimurium"
        elif serotype in ["Java", "Paratyphi B var. L(+) tartrate+"]:
            lookupSero = "Paratyphi B var. Java"
        elif subgenusRegex.match(serotype):
            lookupSero = subgenusRegex.match(serotype).group(1)
        else:
            lookupSero = serotype
        lookupSero = lookupSero.strip()
        lookupSero = anySpaceRegex.sub("_", lookupSero)
        lookupSero = lookupSero.upper()
        if lookupSero not in seroDict:
            print("Warning, serotype not in lookup dict:", lookupSero)
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
    return limsSerotypes, limsSerogroup, limsSubgenus

def apply_rules(limsSerotypes, limsSerogroup, limsSubgenus, row):
    sampleID, consensus, rawCount, readCount, gc, kmerid, contamFlag, most, mostLight, st, mlst, mlstMeanCov, seqsero, seqseroComment, n50, serogr, serovar, seroAnt, seroCGMLST, vaccine, mono, sseJ, readRange, numContigs, numContigs25Kb, numContigs50Kb, assemblySize, assemblyGC, L50 = row
    numReads = safe_int(rawCount)
    assemblySize = safe_int(assemblySize)
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
    limsStatus = "CheckRequired"
    # RULE 1 READCOUNT
    if numReads < 500000 or numReads == "no_result":
        limsReason = "InsufficientData: readCount<50K"
        print("Low read count:", numReads)
        limsStatus = "Inconclusive"
    # RULE 2 LOW KMERID ENTERICA
    elif limsSubgenus == "I" and salmPercent < 75:
        limsReason = "Contaminated: EntericaKmerID<75%"
        print("Low KmerID score for Enterica (< 75%):", salmPercent)
        limsStatus = "Inconclusive"
    # RULE 3 LOW KMERID NON-ENTERICA
    elif limsSubgenus in ("II", "IIIa", "IIIb", "IV", "V") and salmPercent < 38:
        limsReason = "Contaminated: non-EntericaKmerID<38%"
        print("Low KmerID score for non-Enterica (< 38%):", salmPercent)
        limsStatus = "Inconclusive"
    # RULE 4 RED LIGHT
    elif mostLight == "RED":
        limsReason = "Contaminated: MOSTlightRED"
        print("Most light:", mostLight)
        limsStatus = "Inconclusive"
    # RULE 5 SEQSERO2 COMMENT
    elif "Co-existence of multiple serotypes detected" in seqseroComment:
        if assemblySize < 5800000 and assemblySize > 4000000 and n50 > 20000 and numContigs < 600 and mlstMeanCov > 20 and len(limsSerotypes) == 1:
            limsStatus = "CheckRequired"
            limsReason = "multipleSerotypesDetected(SeqSero2)"
        else:
            limsStatus = "Inconclusive"
            limsReason = "Contaminated: multiSerotypes(SeqSero2)"
        print("SeqSero2 comment:", seqseroComment)
    # RULE 6 ASSEMBLY SIZE UPPER
    elif assemblySize > 5800000:
        limsReason = "Contaminated: assembly>5.8Mbp"
        print("Assembly too large:", assemblySize)
        limsStatus = "Inconclusive"
    # RULE 7 N50
    elif n50 < 20000:
        limsReason = "PoorAssembly: N50<20Kbp"
        print("N50 too small:", n50)
        limsStatus = "Inconclusive"
    # RULE 8 CONTIG COUNT
    elif numContigs > 600:
        limsReason = "PoorAssembly: contigCount>600"
        print("Too many contigs:", numContigs)
        limsStatus = "Inconclusive"
    # RULE 9 LOW MLST COVERAGE TYPH
    elif limsSerotypes == ["Typhimurium"] and mlstMeanCov < 20:
        limsReason = "PoorAssembly: MLSTcov<20xTyphimurium"
        print("MLST coverage:", mlstMeanCov)
        limsStatus = "Inconclusive"
    # RULE 10 LOW MLST COVERAGE
    elif limsSerotypes != ["Typhimurium"] and mlstMeanCov < 30:
        limsReason = "PoorAssembly: MLSTcov<30x"
        print("MLST coverage:", mlstMeanCov)
        limsStatus = "Inconclusive"
    # RULE 11 ASSEMBLY SIZE LOWER
    elif assemblySize < 4000000:
        limsReason = "PoorAssembly: assembly<4Mbp"
        print("Assembly too small:", assemblySize)
        limsStatus = "Inconclusive"
    # RULE 12 NO RESULTS
    elif len([x for x in limsSerotypes if x in ('No Type', 'No Results')]) == len(limsSerotypes):
        limsReason = "Contaminated: noIDedSerotypes"
        print("No identified serotypes:", limsSerotypes)
        limsStatus = "Inconclusive"
    # RULE 13 LOW MLST COVERAGE
    elif st == "Failed(incomplete locus coverage)":
        limsReason = "PoorAssembly: incomplSTcov(MOST)"
        print("Incomplete ST locus coverage")
        limsStatus = "Inconclusive"
    # RULE 14 ARIZONAE IIIA 18:Z4:Z32
    elif consensus in ("2-IIIa--1-IIIa 18:z4,z23:-", "2-IIIa 18:z4,z23:---1-Arizonae"):
        limsReason = "check Serovar"
        limsSerotype = "Arizonae IIIa 18:z4,z32:-"
        limsStatus = "CheckRequired"
    # RULE 15 ARIZONAE IIIA 44:Z4:Z23
    elif consensus == "2-IIIa 44:z4,z23:---1-Arizonae":
        limsReason = "check Serovar"
        limsSerotype = "Arizonae IIIa 44:z4,z23:-"
        limsStatus = "CheckRequired"
    # RULE 16 HOUTENAE
    elif consensus == "2-IV 44:z4,z23:---1-Unnamed":
        limsReason = "check Serovar"
        limsSerotype = "Houtenae IV 44:z4,z23:-"
        limsStatus = "CheckRequired"
    # RULE 17 PARATYPHI B JAVA
    elif len([x for x in limsSerotypes if x in ("4:b:-", "1,4,[5],12:b:-", "Paratyphi", "Paratyphi B var. Java")]) == len(limsSerotypes) and sseJ == 'Java':
        limsSerotype = "Paratyphi B var. Java"
        limsVariant = "Paratyphi B var. Java"
        limsStatus = "Pass"
    # RULE 18 PARATYPHI B
    elif len([x for x in limsSerotypes if x in ("4:b:-", "1,4,[5],12:b:-", "Paratyphi", "Paratyphi B")]) == len(limsSerotypes) and sseJ == 'Paratyphi':
        limsSerotype = "Paratyphi B"
        limsStatus = "Pass"
    # RULE 19 MONO IDIKAN
    elif mono == "MonophasicIdikan":
        limsVariant = "Monophasic Idikan"
        limsStatus = "Pass"
    # RULE 20 MONO KEDOUGOU
    elif mono == "MonophasicKedougou":
        limsVariant = "Monophasic Kedougou"
        limsStatus = "Pass"
    # RULE 21 4,12:D:-
    elif consensus == "1-No Type--1-Unnamed--1-I 4,[5],12:d:-":
        limsSerotype = "I 4,[5],12:d:-"
        limsStatus = "Pass"
    # RULE 22 DIARIZONAE
    elif consensus in ("1-IIIb 61:k:1,5,(7)--1-IIIb O:61:k:1,5,7--1-Arizonae", "1-No Type--1-Arizonae--1-O61:k:1,5,7"):
        limsSerotype = "Diarizonae IIIb O:61:k:1,5,7"
        limsStatus = "Pass"
    # RULE 23 BOVISMORBIFICANS
    elif len([x for x in limsSerotypes if x in ("Bovismorbificans", "Bovis-Morbificans")]) == len(limsSerotypes):
        limsSerotype = "Bovismorbificans"
        limsStatus = "Pass"
    # RULE 24 ALL GOOD
    elif len(limsSerotypes) == 1:
        limsStatus = "Pass"
    # RULE 25 EVERYTHING ELSE
    else:
        limsReason = "unknown"
        limsStatus = "CheckRequired"
    if vaccine not in ("NA", "srst2 result file not found"):
        if "2-AviproE" in vaccine:
            limsVaccine = "AVIPRO VAC E Vaccine strain"
        elif "2-Salmovac440" in vaccine:
            limsVaccine = "Salmovac 440 Vaccine strain"
        elif "2-AviproT" in vaccine:
            limsVaccine = "AVIPRO VAC T Vaccine strain"
        elif "2-SalmoporcSTM" in vaccine:
            limsVaccine = "Salmoporc STM Vaccine strain"
        elif "2-Nobilis" in vaccine:
            limsVaccine = "Nobilis SG 9R Vaccine strain"
        elif "2-Wild_Type" in vaccine:
            limsVaccine = "Field strain"
        elif vaccine == "527E_AviproT*":
            limsVaccine = "AVIPRO VAC T Vaccine strain"
        elif "526_21b_AviproT*_527E_AviproT" in vaccine:
            limsVaccine = "AVIPRO VAC T Vaccine strain"
        elif "8321_AviproE*_8322_AviproE*" in vaccine:
            limsVaccine = "AVIPRO VAC E Vaccine strain"
        else:
            limsVaccine = vaccine
    # print("LIMS serotype:", limsSerotype)
    # print("LIMS status:", limsStatus.upper())
    # print("LIMS reason:", limsReason, "\n")
    return limsStatus, limsReason, limsSerotype, limsVariant, limsVaccine
    # numReads, assemblySize, n50, numContigs, mostLight, kmerid, st, mlstMeanCov, contamFlag, vaccine, mono, sseJ


def parse_table(summaryTable):
    outTable = []
    outTable.append(outTableHeader)
    outFileName = os.path.basename(summaryTable).replace(".csv", "_plusLIMS.csv")
    with open(summaryTable, 'r') as inFile:
        reader = csv.reader(inFile, delimiter=',')
        for idx, row in enumerate(reader):
            if idx == 0:
                if row != headerRow:
                    print(headerRow)
                    print(row)
                    sys.exit("Error, table header not in expected format")
                continue
            sampleID, consensus = row[0:2]
            serotypes = parse_consensus(consensus)
            # print("Consensus serotypes:", serotypes)
            limsSerotypes, limsSerogroup, limsSubgenus = parse_seros(serotypes)
            limsStatus, limsReason, limsSerotype, limsVariant, limsVaccine = apply_rules(limsSerotypes, limsSerogroup, limsSubgenus, row)
            outTable.append(str(x) for x in [sampleID, consensus, limsStatus, limsReason, limsSerotype, limsSubgenus, limsSerogroup, limsVariant, limsVaccine] + 
                row[2:])
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
