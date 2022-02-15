#!/usr/bin/python3
# lims_rules.py by Josh

import sys
import csv
import re
import os
import argparse

headerRow = ['StrainID', 'Consensus', '#ReadsR1', 'GC%R1', 'R1Kmerid', 'ContaminationFlag', 'MOST', 
'Most_light', 'st', 'MLST', 'MLST mean cov', 'SeqSero', 'SS comment', 'N50', 'serogroup', 'serovar', 
'serovar_antigen', 'serovar_cgmlst', 'vaccine', 'mono', 'sseJ', 
"ReadLenRange", "#Contigs", "#Contigs>25Kbp", "#Contigs>50Kbp", "AssemblySize", "AssemblyGC", "L50"]

consensusRegex = re.compile("^([123]-(.+?))(--1-(.+?)){0,1}(--1-(.+)){0,1}$")
subgenusRegex = re.compile("^[IVab]+ (.+)$")

def build_sero_dict():
    seroDict = {}
    with open('serogroup_lookup_table.tsv', 'r') as seroLookup:
        reader = csv.reader(seroLookup, delimiter='\t')
        for idx, row in enumerate(reader):
            if idx == 0:
                continue
            serotype, serogroup = row[1], row[6]
            seroDict[serotype.upper()] = [serogroup]
    with open('subgenus_lookup_table.tsv', 'r') as subgLookup:
        reader = csv.reader(subgLookup, delimiter='\t')
        for idx, row in enumerate(reader):
            if idx == 0:
                continue
            serotype, subgenus = row[0:2]
            if serotype in seroDict:
                seroDict[serotype].append(subgenus)
            else:
                seroDict[serotype] = ["undetermined",subgenus]
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
    for serotype in serotypes:
        if serotype == "Monophasic Typhimurium":
            lookupSero = "Typhimurium"
        elif subgenusRegex.match(serotype):
            lookupSero = subgenusRegex.match(serotype).group(1)
        else:
            lookupSero = serotype
        if lookupSero.upper() not in seroDict:
            print("Warning, serotype not in lookup dict:", lookupSero)
            serogroup = "undetermined"
            subgenus = "undetermined"
        else:
            serogroup = seroDict[lookupSero.upper()][0]
            if len(seroDict[lookupSero.upper()]) == 2:
                subgenus = seroDict[lookupSero.upper()][1]
            else:
                subgenus = "undetermined"
        serogroups.append(serogroup)
        subgenera.append(subgenus)
    if len(serogroups) == 3:
        limsSerogroup = [x for x in serogroups if serogroups.count(x) > 1]
        if len(limsSerogroup) == 1:
            limsSerogroup = limsSerogroup[0]
        else:
            limsSerogroup = "undetermined"
    else:
        limsSerogroup = serogroups[0]
    if len(subgenera) == 3:
        limsSubgenus = [x for x in subgenera if subgenera.count(x) > 1]
        if len(limsSubgenus) == 1:
            limsSubgenus = limsSubgenus[0]
        else:
            limsSubgenus = "undetermined"
    else:
        limsSubgenus = subgenera[0]
    return limsSerogroup, limsSubgenus

def apply_rules(serotypes, limsSerogroup, limsSubgenus, row):
    sampleID, consensus, numReads, gc, kmerid, contamFlag, most, mostLight, st, mlst, mlstMeanCov, seqsero, seqseroComment, n50, serogr, serovar, seroAnt, seroCGMLST, vaccine, mono, sseJ, readRange, numContigs, numContigs25Kb, numContigs50Kb, assemblySize, assemblyGC, L50 = row
    numReads = int(numReads)
    assemblySize = int(assemblySize)
    n50 = int(n50)
    numContigs = int(numContigs)
    mlstMeanCov = float(mlstMeanCov)
    kmerids = kmerid.split("--")
    salmPercent = [x.split("-")[0] for x in kmerids if "salmonella" in x]
    if len(salmPercent) == 1:
        salmPercent = float(salmPercent[0])
    else:
        salmPercent = 0
    serotypes = [subgenusRegex.match(x).group(1) if subgenusRegex.match(x) else x for x in serotypes]
    if len(serotypes) == 1:
        limsSerotype = serotypes[0]
    elif len(serotypes) == 2:
        limsSerotype = serotypes[0]+" (2/3consen)"
    else:
        limsSerotype = "no_consensus"
    limsVariant = ""
    limsVaccine = ""
    limsStatus = "CheckRequired"
    # RULE 1 READCOUNT
    if numReads < 500000:
        print("Low read count:", numReads)
        limsStatus = "InsufficientData/PoorAssembly"
    # RULE 2 ASSEMBLY SIZE LOWER
    elif assemblySize < 4000000:
        print("Assembly too small:", assemblySize)
        limsStatus = "InsufficientData/PoorAssembly"
    # RULE 3 ASSEMBLY SIZE UPPER
    elif assemblySize > 5800000:
        print("Assembly too large:", assemblySize)
        limsStatus = "InsufficientData/PoorAssembly"
    # RULE 4 N50
    elif n50 < 20000:
        print("N50 too small:", n50)
        limsStatus = "InsufficientData/PoorAssembly"
    # RULE 5 CONTIG COUNT
    elif numContigs > 600:
        print("Too many contigs:", numContigs)
        limsStatus = "InsufficientData/PoorAssembly"
    # RULE 6 NO RESULTS
    elif len([x for x in serotypes if x in ('No Type', 'No Results')]) == len(serotypes):
        print("No identified serotypes:", serotypes)
        limsStatus = "Contamination"
    # RULE 7 RED LIGHT
    elif mostLight == "RED":
        print("Most light:", mostLight)
        limsStatus = "Contamination"
    # RULE 8 LOW MLST COVERAGE
    elif st == "Failed(incomplete locus coverage)":
        print("Incomplete ST locus coverage")
        limsStatus = "Contamination"
    # RULE 9 LOW KMERID ENTERICA
    elif limsSubgenus == "I" and salmPercent < 75:
        print("Low KmerID score for Enterica (< 75%):", salmPercent)
        limsStatus = "Contamination"
    # RULE 10 LOW KMERID NON-ENTERICA
    elif limsSubgenus in ("II", "IIIa", "IIIb", "IV", "V") and salmPercent < 38:
        print("Low KmerID score for non-Enterica (< 38%):", salmPercent)
        limsStatus = "Contamination"
    # RULE 11 ARIZONAE IIIA 18:Z4:Z32
    elif consensus in ("2-IIIa--1-IIIa 18:z4,z23:-", "2-IIIa 18:z4,z23:---1-Arizonae"):
        limsSerotype = "Arizonae IIIa 18:z4,z32:-"
        limsStatus = "CheckRequired"
    # RULE 12 ARIZONAE IIIA 44:Z4:Z23
    elif consensus == "2-IIIa 44:z4,z23:---1-Arizonae":
        limsSerotype = "Arizonae IIIa 44:z4,z23:-"
        limsStatus = "CheckRequired"
    # RULE 13 HOUTENAE
    elif consensus == "2-IV 44:z4,z23:---1-Unnamed":
        limsSerotype = "Houtenae IV 44:z4,z23:-"
        limsStatus = "CheckRequired"
    # RULE 14 LOW MLST COVERAGE TYPH
    elif serotypes == ["Typhimurium"] and mlstMeanCov < 20:
        limsStatus = "CheckRequired"
    # RULE 15 LOW MLST COVERAGE
    elif serotypes in [["Enteritidis"], ["Infantis"]] and mlstMeanCov < 30:
        limsStatus = "CheckRequired"
    # RULE 16 PARATYPHI B JAVA
    elif len([x for x in serotypes if x in ("4:b:-", "1,4,[5],12:b:-", "Paratyphi", "Paratyphi B var. Java")]) == len(serotypes) and sseJ == 'Java':
        limsSerotype = "Paratyphi B var. Java"
        limsStatus = "Pass"
    # RULE 17 PARATYPHI B
    elif len([x for x in serotypes if x in ("4:b:-", "1,4,[5],12:b:-", "Paratyphi", "Paratyphi B")]) == len(serotypes) and sseJ == 'Paratyphi':
        limsVariant = "Paratyphi B"
        limsStatus = "Pass"
    # RULE 18 MONO IDIKAN
    elif mono == "MonophasicIdikan":
        limsVariant = "Monophasic Idikan"
        limsStatus = "Pass"
    # RULE 19 MONO KEDOUGOU
    elif mono == "MonophasicKedougou":
        limsVariant = "Monophasic Kedougou"
        limsStatus = "Pass"
    # RULE 20 4,12:D:-
    elif consensus == "1-No Type--1-Unnamed--1-I 4,[5],12:d:-":
        limsSerotype = "I 4,[5],12:d:-"
        limsStatus = "Pass"
    # RULE 21 DIARIZONAE
    elif consensus in ("1-IIIb 61:k:1,5,(7)--1-IIIb O:61:k:1,5,7--1-Arizonae", "1-No Type--1-Arizonae--1-O61:k:1,5,7"):
        limsSerotype = "Diarizonae IIIb O:61:k:1,5,7"
        limsStatus = "Pass"
    # RULE 22 ALL GOOD
    elif len(serotypes) == 1:
        limsStatus = "Pass"
    # RULE 23 EVERYTHING ELSE
    else:
        limsStatus = "CheckRequired"
    if vaccine not in ("NA", "srst2 result file not found"):
        limsVaccine = vaccine
    print("LIMS serotype:", limsSerotype)
    print("LIMS status:", limsStatus.upper() + "\n")
    return limsStatus, limsSerotype, limsVariant, limsVaccine, numReads, assemblySize, n50, numContigs, mostLight, kmerid, st, mlstMeanCov, contamFlag, vaccine, mono, sseJ


def parse_table(summaryTable):
    outTable = []
    outTable.append(["StrainID", "Consensus", "LIMS_Status", "LIMS_SerotypeID", "LIMS_Subgenus", "LIMS_Serogroup", "LIMS_Variant", "LIMS_Vaccine", 
    "#ReadsR1", "AssemblySize", "N50", "#Contigs", "Most_light", "R1Kmerid", "st", "MLST_mean_cov", "ContaminationFlag"])
    outFileName = os.path.basename(summaryTable).replace("_SummaryTable.csv", "_LIMS.csv")
    with open(summaryTable, 'r') as inFile:
        reader = csv.reader(inFile, delimiter=',')
        for idx, row in enumerate(reader):
            if idx == 0:
                if row != headerRow:
                    print(row)
                    sys.exit("Error, table header not in expected format")
                continue
            sampleID, consensus = row[0:2]
            serotypes = parse_consensus(consensus)
            print("Consensus serotypes:", serotypes)
            limsSerogroup, limsSubgenus = parse_seros(serotypes)
            limsStatus, limsSerotype, limsVariant, limsVaccine, numReads, assemblySize, n50, numContigs, mostLight, kmerid, st, mlstMeanCov, contamFlag, vaccine, mono, sseJ = apply_rules(serotypes, limsSerogroup, limsSubgenus, row)
            outTable.append(str(x) for x in [sampleID, consensus, limsStatus, limsSerotype, limsSubgenus, limsSerogroup, limsVariant, limsVaccine, numReads, assemblySize, n50, numContigs, mostLight, kmerid, st, mlstMeanCov, contamFlag, vaccine, mono, sseJ])
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
