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


for entry in "/home/WGS_Data/TestIsolates"/*
do
  echo "$entry"
done

for entry in "/home/WGS_Data/TestIsolates/ERR2230776"/*
do
  echo "$entry"
done

for entry in "/home/WGS_Data/TestIsolates/ERR2231037"/*
do
  echo "$entry"
done

for entry in "/home/WGS_Data/TestIsolates/ERR2235662"/*
do
  echo "$entry"
done

for entry in "/home/WGS_Data/TestIsolates/ERR2208776"/*
do
  echo "$entry"
done

for entry in "/WGS_Data"/*
do
  echo "$entry"
done

for entry in "/home/WGS_Data"/*
do
  echo "$entry"
done

for entry in "/home/WGS_Data/test_isolates"/*
do
  echo "$entry"
done

for entry in "/home/WGS_Data/test_isolates/ERR2230776"/*
do
  echo "$entry"
done

for entry in "/home/WGS_Data/test_isolates/ERR2231037"/*
do
  echo "$entry"
done

for entry in "/home/WGS_Data/test_isolates/ERR2235662"/*
do
  echo "$entry"
done

for entry in "/home/WGS_Data/test_isolates/ERR2208776"/*
do
  echo "$entry"
done



for entry in "/home/WGS_Results/TestIsolates"/*
do
  echo "$entry"
done

for entry in "/home/WGS_Results/TestIsolates/ERR2230776"/*
do
  echo "$entry"
done

for entry in "/home/WGS_Results/TestIsolates/ERR2231037"/*
do
  echo "$entry"
done

for entry in "/home/WGS_Results/TestIsolates/ERR2235662"/*
do
  echo "$entry"
done

for entry in "/home/WGS_Results/TestIsolates/ERR2208776"/*
do
  echo "$entry"
done

for entry in "/home/WGS_Results"/*
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

for entry in "/home/WGS_Results/test_isolates/ERR2231037"/*
do
  echo "$entry"
done

for entry in "/home/WGS_Results/test_isolates/ERR2235662"/*
do
  echo "$entry"
done

for entry in "/home/WGS_Results/test_isolates/ERR2208776"/*
do
  echo "$entry"
done

#Check results
WGS_CLUSTER_CSV=$(print_todays_wgs_cluster)
assert_first_csv_row $WGS_CLUSTER_CSV "Consensus" "$serovar"
assert_first_csv_row $WGS_CLUSTER_CSV "serovar" "$serovar"
