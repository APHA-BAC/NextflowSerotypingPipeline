#!/bin/bash

wget --no-check-certificate --content-disposition https://raw.githubusercontent.com/JaromirGuzinski/NextflowSerotypingPipeline/master/pipeline_component_software/srst2/srst2.zip
unzip srst2.zip
sudo cp -r srst2 /opt
sudo chmod +x /opt/srst2/scripts/srst2.py
sudo ln -s /opt/srst2/scripts/srst2.py /usr/local/bin
rm -rf srst2
rm srst2.zip

pip3 install scipy
