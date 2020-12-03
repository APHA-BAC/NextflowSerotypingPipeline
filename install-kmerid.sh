#!/bin/bash

wget --no-check-certificate --content-disposition https://raw.githubusercontent.com/JaromirGuzinski/NextflowSerotypingPipeline/master/pipeline_component_software/KmerID/kmerid.zip
sudo unzip kmerid.zip -d /opt/ 
rm kmerid.zip

aws s3 cp s3://s3-staging-area/kmerid_ref/kmerid_ref.zip kmerid_ref.zip
sudo unzip kmerid_ref.zip
sudo rm -rf kmerid_ref.zip
sudo mv ref /opt/kmerid/ref

cd /opt/kmerid
sudo make all
sudo chmod +x /opt/kmerid/setup_refs.py
sudo chmod +x /opt/kmerid/kmerid_python3.py
