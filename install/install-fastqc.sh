#!/bin/bash

set -e

######################################################################################
# FastQC
cd $HOME

sudo DEBIAN_FRONTEND=noninteractive apt-get -y install unzip
sudo DEBIAN_FRONTEND=noninteractive apt-get -y install openjdk-11-jdk
sudo DEBIAN_FRONTEND=noninteractive apt-get -y install default-jre

wget --no-check-certificate --content-disposition https://raw.githubusercontent.com/APHA-BAC/NextflowSerotypingPipeline/master/pipeline_component_software/FastQC/fastqc_v0.11.9.zip
unzip fastqc_v0.11.9.zip
sudo cp -r FastQC /opt
sudo chmod +x /opt/FastQC/fastqc
sudo ln -s /opt/FastQC/fastqc /usr/local/bin
rm -rf FastQC
rm fastqc_v0.11.9.zip

# wget --no-check-certificate --content-disposition https://github.com/s-andrews/FastQC/archive/refs/heads/master.zip
