#!/bin/bash
#
#================================================================
# nextflowtest.bash
#================================================================
#
#% DESCRIPTION
#%    Integration tests for nextflow should be called this way. 

/NextflowSerotypingPipeline/nextflow run SCE3_pipeline_update.nf \
--outdir "/WGS_Results/" \
--reads "/WGS_Data/*_{R1,R2}.fastq.gz" \
-with-report /artifacts/report.html
