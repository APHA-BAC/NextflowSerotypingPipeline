#!/bin/bash
sudo ln -s /usr/bin/python3 /usr/bin/python
sudo apt-get install python3-pip
sudo python -m pip install -U pip
sudo python -m pip install -U matplotlib

wget --no-check-certificate --content-disposition https://raw.githubusercontent.com/APHA-BAC/NextflowSerotypingPipeline/master/pipeline_component_software/quast/quast-5.1.0rc1.tar.gz
tar -xvf quast-5.1.0rc1.tar.gz

sudo cp -r quast-5.1.0rc1 /opt
cd /opt/quast-5.1.0rc1
sudo ./install_full.sh
cd
sudo ln -s /opt/quast-5.1.0rc1/quast.py /usr/local/bin
rm -fr quast-5.1.0rc1
rm quast-5.1.0rc1.tar.gz

