#!/bin/bash

wget https://repo.anaconda.com/miniconda/Miniconda3-py37_4.8.3-Linux-x86_64.sh
sudo chmod +x Miniconda3-py37_4.8.3-Linux-x86_64.sh 
./Miniconda3-py37_4.8.3-Linux-x86_64.sh -b
rm -rf Miniconda3-py37_4.8.3-Linux-x86_64.sh

conda install -c conda-forge -c bioconda -c defaults shovill=0.9.0


