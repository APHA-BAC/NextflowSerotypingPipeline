#!/bin/bash

wget --no-check-certificate --content-disposition https://raw.githubusercontent.com/JaromirGuzinski/NextflowSerotypingPipeline/master/pipeline_component_software/most/most.zip
unzip most.zip
sudo cp -r most /opt
sudo chmod +x /opt/most/MOST-master/MOST.py
rm -rf most
rm most.zip

pip3 install lxml==4.5.2
pip3 install biopython==1.73
conda install -c bioconda emboss
conda install -c kantorlab blastn
sudo apt-get install libncurses5






