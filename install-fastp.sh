#!/bin/bash
set -e

######################################################################################
# fastp
cd $HOME

wget --no-check-certificate --content-disposition http://opengene.org/fastp/fastp
sudo cp fastp /opt
sudo chmod +x /opt/fastp
sudo ln -s /opt/fastp /usr/local/bin
rm fastp
