#!/bin/bash
#
#================================================================
# nextflowtest.bash
#================================================================
#
#% DESCRIPTION
#%    Integration tests for nextflow should be called this way. 

$HOME/nextflow/nextflow $HOME/NextflowSerotypingPipeline/SCE3_pipeline_update.nf -with-report /artifacts/report.html
# \
# --outdir "$HOME/WGS_Results/TestIsolates/" \
# --reads "$HOME/WGS_Data/TestIsolates/*_{R1,R2}.fastq.gz" \
