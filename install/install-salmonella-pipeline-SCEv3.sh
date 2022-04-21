#!/bin/bash

#fastqc

cp ~/NextflowSerotypingPipeline/install/install-fastqc.sh ./
bash install-fastqc.sh

#shovill

cp ~/NextflowSerotypingPipeline/install/install-shovill.sh ./
bash install-shovill.sh

#quast

cp ~/NextflowSerotypingPipeline/install/install-quast.sh ./
bash install-quast.sh

#most

cp ~/NextflowSerotypingPipeline/install/install-most.sh ./
bash install-most.sh

#kmerid

cp ~/NextflowSerotypingPipeline/install/install-kmerid.sh ./
bash install-kmerid.sh

#seqsero2

cp ~/NextflowSerotypingPipeline/install/install=seqsero2.sh ./
bash install-seqsero2.sh

#sistr

cp ~/NextflowSerotypingPipeline/install/install-sistr.sh ./
bash install-sistr.sh

#srst2

cp ~/NextflowSerotypingPipeline/install/install-srst2.sh ./
bash install-srst2.sh

#Nextflow_SummaryTable_FolderPrep

cp ~/NextflowSerotypingPipeline/install/install-Nextflow_SummaryTable_FolderPrep.sh ./
bash install-Nextflow_SummaryTable_FolderPrep.sh

sudo reboot
