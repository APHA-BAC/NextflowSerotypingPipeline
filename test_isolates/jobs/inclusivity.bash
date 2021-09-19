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
mv ${accession}_1.fastq.gz /WGS_Data/${accession}_R1.fastq.gz
mv ${accession}_2.fastq.gz /WGS_Data/${accession}_R2.fastq.gz

# Run nextflow
nextflowtest


#check dir contents
for file in /NextflowSerotypingPipeline/WGS_Data/*; do
  echo "${file##*/}"
done

for file in /NextflowSerotypingPipeline/WGS_Results/*; do
  echo "${file##*/}"
done

for entry in "/NextflowSerotypingPipeline/WGS_Results"/*
do
  echo "$entry"
done


for file in /WGS_Data/*; do
  echo "${file##*/}"
done

for file in /WGS_Results/*; do
  echo "${file##*/}"
done

for entry in "/WGS_Results"/*
do
  echo "$entry"
done


for file in /WGS_Data/test_isolates/*; do
  echo "${file##*/}"
done

for file in /WGS_Results/test_isolates/*; do
  echo "${file##*/}"
done

for entry in "/WGS_Results/test_isolates"/*
do
  echo "$entry"
done


for file in /home/WGS_Data/*; do
  echo "${file##*/}"
done

for file in /home/WGS_Results/*; do
  echo "${file##*/}"
done

for entry in "/home/WGS_Results"/*
do
  echo "$entry"
done


for file in /home/$USER/WGS_Data/*; do
  echo "${file##*/}"
done

for file in /home/$USER/WGS_Results/*; do
  echo "${file##*/}"
done

for entry in "/home/$USER/WGS_Results"/*
do
  echo "$entry"
done


for file in /home/$USER/WGS_Data/test_isolates*; do
  echo "${file##*/}"
done

for file in /home/$USER/WGS_Results/test_isolates*; do
  echo "${file##*/}"
done

for entry in "/home/$USER/WGS_Results/test_isolates"/*
do
  echo "$entry"
done


for entry in "/home/WGS_Results/test_isolates"/*
do
  echo "$entry"
done


for entry in "/home/WGS_Results/test_isolates/ERR2230776"/*
do
  echo "$entry"
done
#check dir contents


# Check results
WGS_CLUSTER_CSV=$(print_todays_wgs_cluster)
assert_first_csv_row $WGS_CLUSTER_CSV "Outcome" "Pass"
assert_first_csv_row $WGS_CLUSTER_CSV "serovar" "$serovar"
