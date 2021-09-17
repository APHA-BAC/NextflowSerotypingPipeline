#!/bin/bash
#
#================================================================
# inclusivity.bash
#================================================================
#
#% DESCRIPTION
#%    Tests the nextflow pipeline on a minimal dataset

# Import
source /NextflowSerotypingPipeline/test_isolates/utils/aliases.bash

# Args
TESTCASE=$1

serovar=$(print_csv_value './NextflowSerotypingPipeline/test_isolates/data/inclusivity_cases.csv' $TESTCASE serovar)
accession=$(print_csv_value './NextflowSerotypingPipeline/test_isolates/data/inclusivity_cases.csv' $TESTCASE accession)

# Fetch SRA Data
prefetch $accession -O ./
fasterq-dump ./$accession
rm ./$accession/*.sra
rm -r ./$accession

gzip ${accession}_1.fastq ${accession}_2.fastq
mv ${accession}_1.fastq.gz /WGS_Data/${accession}_R1.fastq.gz
mv ${accession}_2.fastq.gz /WGS_Data/${accession}_R2.fastq.gz

# Run nextflow
nextflowtest

# Check results
WGS_CLUSTER_CSV=$(print_todays_wgs_cluster)
assert_first_csv_row $WGS_CLUSTER_CSV "Outcome" "Pass"
assert_first_csv_row $WGS_CLUSTER_CSV "serovar" "$serovar"
