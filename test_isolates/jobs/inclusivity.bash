#!/bin/bash
#
#================================================================
# inclusivity.bash
#================================================================
#
#% DESCRIPTION
#%    Tests the nextflow pipeline on a minimal dataset

set -e

# Import
source test_isolates/utils/aliases.bash

# Args
TESTCASE=$1

serovar=$(print_csv_value './test_isolates/data/inclusivity_cases.csv' $TESTCASE serovar)
accession=$(print_csv_value './test_isolates/data/inclusivity_cases.csv' $TESTCASE accession)

# Fetch SRA Data
# TODO: call prefetch from $PATH rather than /usr/local/bin/
#       it is being explicitly caleed in this way to disambiguate from the conda install
#       which is prioritised in the path.
/usr/local/bin/prefetch $accession -O ./
/usr/local/bin/fasterq-dump ./$accession
rm ./$accession/*.sra
rm -r ./$accession

gzip ${accession}_1.fastq ${accession}_2.fastq
mv ${accession}_1.fastq.gz $HOME/WGS_Data/TestIsolates/${accession}_R1.fastq.gz
mv ${accession}_2.fastq.gz $HOME/WGS_Data/TestIsolates/${accession}_R2.fastq.gz

# Run nextflow
nextflowtest

#Check results
SUMMARY_TABLE=$(print_path_to_summary_table)
assert_first_csv_row $SUMMARY_TABLE "Consensus" "$serovar"
