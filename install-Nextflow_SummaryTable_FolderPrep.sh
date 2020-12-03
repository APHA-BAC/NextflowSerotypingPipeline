#!/bin/bash

mkdir /home/$USER/nextflow
cd /home/$USER/nextflow
curl -s https://get.nextflow.io | bash
wget --no-check-certificate --content-disposition https://raw.githubusercontent.com/JaromirGuzinski/NextflowSerotypingPipeline/master/SCE3_pipeline_update.nf

mkdir /home/$USER/summary
cd /home/$USER/summary
wget --no-check-certificate --content-disposition https://raw.githubusercontent.com/JaromirGuzinski/NextflowSerotypingPipeline/master/summaryTable_Python3.py

mkdir /home/$USER/WGS_Data
mkdir /home/$USER/WGS_Results




