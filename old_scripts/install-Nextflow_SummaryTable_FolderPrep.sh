#!/bin/bash

set -e

mkdir /home/$USER/nextflow
cd /home/$USER/nextflow
curl -s https://get.nextflow.io | bash


cp ~/NextflowSerotypingPipeline/SCE3_pipeline_update.nf ./

mkdir /home/$USER/summary
cd /home/$USER/summary

cp ~/NextflowSerotypingPipeline/summaryTable_Python3.py ./

mkdir /home/$USER/WGS_Data
mkdir /home/$USER/WGS_Results

