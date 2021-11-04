#!/bin/bash
#
#================================================================
# nextflowtest.bash
#================================================================
#
#% DESCRIPTION
#%    Integration tests for nextflow should be called this way. 

$HOME/nextflow/nextflow run $HOME/nextflow/SCE3_pipeline_update.nf \
--outdir "$HOME/WGS_Results/TestIsolates/" \
--reads "$HOME/WGS_Data/TestIsolates/*_{R1,R2}.fastq.gz" \
-with-report /artifacts/report.html
