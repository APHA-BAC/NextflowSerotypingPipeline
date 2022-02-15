#!/bin/bash
#
#================================================================
# nextflowtest.bash
#================================================================
#
#% DESCRIPTION
#%    Integration tests for nextflow should be called this way. 

$HOME/nextflow/nextflow run ./SCE3_pipeline_update.nf -with-report /artifacts/report.html
