#!/bin/bash
set -e

######################################################################################
# fastp
cd $HOME

wget --no-check-certificate --content-disposition http://opengene.org/fastp/fastp.0.23.1
mv fastp.0.23.1 fastp
sudo cp fastp /opt
sudo chmod +x /opt/fastp
sudo ln -s /opt/fastp /usr/local/bin
rm fastp
