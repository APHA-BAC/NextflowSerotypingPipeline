#!/bin/bash

#fastqc
wget --no-check-certificate --content-disposition https://raw.githubusercontent.com/JaromirGuzinski/NextflowSerotypingPipeline/master/install-fastqc.sh
bash install-fastqc.sh

#shovill 
wget --no-check-certificate --content-disposition https://raw.githubusercontent.com/JaromirGuzinski/NextflowSerotypingPipeline/master/install-shovill.sh
bash install-shovill.sh

#quast
wget --no-check-certificate --content-disposition https://raw.githubusercontent.com/JaromirGuzinski/NextflowSerotypingPipeline/master/install-quast.sh
bash install-quast.sh

#most
wget --no-check-certificate --content-disposition https://raw.githubusercontent.com/JaromirGuzinski/NextflowSerotypingPipeline/master/install-most.sh
bash install-most.sh

#kmerid
wget --no-check-certificate --content-disposition https://raw.githubusercontent.com/JaromirGuzinski/NextflowSerotypingPipeline/master/install-kmerid.sh
bash install-kmerid.sh

#seqsero2
wget --no-check-certificate --content-disposition https://raw.githubusercontent.com/JaromirGuzinski/NextflowSerotypingPipeline/master/install-seqsero2.sh
bash install-seqsero2.sh

#sistr
wget --no-check-certificate --content-disposition https://raw.githubusercontent.com/JaromirGuzinski/NextflowSerotypingPipeline/master/install-sistr.sh
bash install-sistr.sh

#srst2
wget --no-check-certificate --content-disposition https://raw.githubusercontent.com/JaromirGuzinski/NextflowSerotypingPipeline/master/install-srst2.sh
bash install-srst2.sh

#Nextflow_SummaryTable_FolderPrep
wget --no-check-certificate --content-disposition https://raw.githubusercontent.com/JaromirGuzinski/NextflowSerotypingPipeline/master/install-Nextflow_SummaryTable_FolderPrep.sh
bash install-Nextflow_SummaryTable_FolderPrep.sh

sudo reboot
