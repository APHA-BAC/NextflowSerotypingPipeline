#!/bin/bash
#
#================================================================
# inclusivity.bash
#================================================================
#
#% DESCRIPTION
#%    Tests the nextflow pipeline on a minimal dataset

# Import
source test_isolates/utils/aliases.bash

# Args
TESTCASE=$1

serovar=$(print_csv_value './test_isolates/data/inclusivity_cases.csv' $TESTCASE serovar)
accession=$(print_csv_value './test_isolates/data/inclusivity_cases.csv' $TESTCASE accession)

# Fetch SRA Data
prefetch $accession -O ./
fasterq-dump ./$accession
rm ./$accession/*.sra
rm -r ./$accession

gzip ${accession}_1.fastq ${accession}_2.fastq
mv ${accession}_1.fastq.gz /home/WGS_Data/TestIsolates/${accession}_R1.fastq.gz
mv ${accession}_2.fastq.gz /home/WGS_Data/TestIsolates/${accession}_R2.fastq.gz

# Run nextflow
nextflowtest

#Check results
SUMMARY_TABLE=$(print_path_to_summary_table)
assert_first_csv_row $SUMMARY_TABLE "Consensus" "$serovar"
