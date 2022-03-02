#!/bin/bash

####################################################################################
# The following section is for Docker ONLY
# Skip to FastQC section if installing pipeline directly into an SCE3 VM

apt-get -y update
DEBIAN_FRONTEND=noninteractive apt-get install -y -q sudo
DEBIAN_FRONTEND=noninteractive apt-get install -y -q apt-utils
DEBIAN_FRONTEND=noninteractive apt-get install -y -q make
DEBIAN_FRONTEND=noninteractive apt-get install -y -q gcc

DEBIAN_FRONTEND=noninteractive apt-get install -y -q git
DEBIAN_FRONTEND=noninteractive apt-get install -y -q wget
DEBIAN_FRONTEND=noninteractive apt-get install -y -q curl

DEBIAN_FRONTEND=noninteractive apt-get install -y -q libz-dev
DEBIAN_FRONTEND=noninteractive apt-get install -y -q libbz2-dev
DEBIAN_FRONTEND=noninteractive apt-get install -y -q libncurses5-dev
DEBIAN_FRONTEND=noninteractive apt-get install -y -q libncursesw5-dev

DEBIAN_FRONTEND=noninteractive apt-get install -y -q libghc-bzlib-prof
DEBIAN_FRONTEND=noninteractive apt-get install -y -q zlib1g-dev
DEBIAN_FRONTEND=noninteractive apt-get install -y -q libcurl4-openssl-dev

DEBIAN_FRONTEND=noninteractive apt-get install -y -q libfile-copy-recursive-perl
DEBIAN_FRONTEND=noninteractive apt-get install -y -q libio-socket-ssl-perl
DEBIAN_FRONTEND=noninteractive apt-get install -y -q libio-tee-perl
DEBIAN_FRONTEND=noninteractive apt-get install -y -q libunicode-string-perl

DEBIAN_FRONTEND=noninteractive apt-get install -y -q nano
DEBIAN_FRONTEND=noninteractive apt-get install -y -q python3

# The following original line is no longernecessary. Conda gives us pip later
# and we should try to avoid root pip installs!
# DEBIAN_FRONTEND=noninteractive apt-get install -y -q python3-pip

######################################################################################
# FastQC
cd $HOME

sudo apt-get -y install unzip
sudo apt-get -y install openjdk-11-jdk
sudo apt-get -y install default-jre

wget --no-check-certificate --content-disposition https://raw.githubusercontent.com/APHA-BAC/NextflowSerotypingPipeline/master/pipeline_component_software/FastQC/fastqc_v0.11.9.zip
unzip fastqc_v0.11.9.zip
sudo cp -r FastQC /opt
sudo chmod +x /opt/FastQC/fastqc
sudo ln -s /opt/FastQC/fastqc /usr/local/bin
rm -rf FastQC
rm fastqc_v0.11.9.zip

# wget --no-check-certificate --content-disposition https://github.com/s-andrews/FastQC/archive/refs/heads/master.zip

#########################################################################################
# Miniconda
cd $HOME

wget https://repo.anaconda.com/miniconda/Miniconda3-py37_4.8.3-Linux-x86_64.sh
sudo chmod +x Miniconda3-py37_4.8.3-Linux-x86_64.sh
sudo ./Miniconda3-py37_4.8.3-Linux-x86_64.sh -b -p /opt/conda
# -b flag is for silent batch installation with no direct PATH modification
rm Miniconda3-py37_4.8.3-Linux-x86_64.sh

export PATH=/opt/conda/bin:$PATH
/opt/conda/bin/conda init bash
source ~/.bashrc

# Conda environment
# Now this is where things change. I create a specific Conda environment for the pipeline
# Any other projects on the SCE3 machine will thus be unaffected by conda and pip installs within this environment
conda create -n salmpipe python=3.7
conda activate salmpipe

# Matplotlib
pip install matplotlib

###############################################################################################

# The following original lines are not needed; we already have Python3 and pip from Conda:
# sudo apt-get -y install python3-pip
# sudo python -m pip install -U pip
# sudo python -m pip install -U matplotlib

# This also avoids using pip as root with sudo, which is strongly discouraged as it may break dependencies in other projects
# Using the pip provided by Conda ensures only our current Conda environment is updated
# We also use "pip install" (Conda) lower down for MOST, so axing "sudo python -m pip" makes sense

###########################################################################################
# Shovill
conda install -c conda-forge -c bioconda -c defaults shovill=0.9.0
shovill --check

###########################################################################################
# Quast
cd $HOME
sudo ln -s /usr/bin/python3 /usr/bin/python

# I have changed this to follow the style of other installations into /opt
wget --no-check-certificate --content-disposition https://raw.githubusercontent.com/APHA-BAC/NextflowSerotypingPipeline/master/pipeline_component_software/quast/quast-5.1.0rc1.tar.gz
tar -xzf quast-5.1.0rc1.tar.gz
sudo cp -r quast-5.1.0rc1 /opt
sudo /opt/quast-5.1.0rc1/install.sh
sudo ln -s /opt/quast-5.1.0rc1/quast.py /usr/local/bin
rm -rf quast-5.1.0rc1
sudo rm -rf quast_test_output
rm quast-5.1.0rc1.tar.gz

# wget --no-check-certificate --content-disposition https://github.com/ablab/quast/archive/refs/heads/master.zip

############################################################################################
# MOST
cd $HOME

sudo apt-get -y install libncurses5
conda install -c bioconda emboss
conda install -c kantorlab blastn
pip install lxml==4.5.2
pip install biopython==1.73

wget --no-check-certificate --content-disposition https://raw.githubusercontent.com/APHA-BAC/NextflowSerotypingPipeline/master/pipeline_component_software/most/most.zip
unzip most.zip
sudo cp -r most /opt
sudo chmod +x /opt/most/MOST-master/MOST.py
rm -rf most
rm most.zip

#############################################################################
# KmerID
cd $HOME

wget --no-check-certificate --content-disposition https://raw.githubusercontent.com/APHA-BAC/NextflowSerotypingPipeline/master/pipeline_component_software/KmerID/kmerid.zip
unzip kmerid.zip
sudo cp -r kmerid /opt
sudo make -C /opt/kmerid all # Changed this so no need to change directory
sudo chmod +x /opt/kmerid/setup_refs.py
sudo chmod +x /opt/kmerid/kmerid_python3.py
rm -rf kmerid
rm kmerid.zip

sudo apt-get -y install nodejs
sudo apt-get -y install npm
sudo npm install -g github-files-fetcher

# fetcher --url="https://github.com/APHA-BAC/NextflowSerotypingPipeline/tree/circleci-project-setup/kmerid_ref/Acinetobacter" --out=/opt/kmerid/ref
sudo fetcher --url="https://github.com/APHA-BAC/NextflowSerotypingPipeline/tree/circleci-project-setup/kmerid_ref/Citrobacter" --out=/opt/kmerid/ref
sudo fetcher --url="https://github.com/APHA-BAC/NextflowSerotypingPipeline/tree/circleci-project-setup/kmerid_ref/Salmonella" --out=/opt/kmerid/ref

###########################################################################################
# SeqSero2
conda install -c bioconda seqsero2=1.1.1

############################################################################################
# Sistr
pip install numpy pandas==1.1.5
# pip install wheel # Requirements already satisfied

conda config --add channels conda-forge
conda config --add channels r
conda config --add channels bioconda
conda install sistr_cmd

#############################################################################################
# Srst2
cd $HOME

wget --no-check-certificate --content-disposition https://raw.githubusercontent.com/APHA-BAC/NextflowSerotypingPipeline/master/pipeline_component_software/srst2/srst2.zip
unzip srst2.zip
sudo cp -r srst2 /opt
sudo chmod +x /opt/srst2/scripts/srst2.py
sudo ln -s /opt/srst2/scripts/srst2.py /usr/local/bin
rm -rf srst2
rm srst2.zip

#################################################################################################
# Nextflow
mkdir $HOME/nextflow
cd $HOME/nextflow
curl -L0 https://github.com/nextflow-io/nextflow/releases/download/v21.04.1/nextflow-21.04.1-all | bash
# This ^ gives a warning, "curl: (23) Failed writing body (4096 != 16384)", but seems to work?

#################################################################################################
# Summary table script
mkdir $HOME/summary # Changed to $HOME for SC3 - We don't want summary installed in /home/summary
cd $HOME/summary
wget --no-check-certificate --content-disposition https://raw.githubusercontent.com/APHA-BAC/NextflowSerotypingPipeline/master/summaryTable_Python3.py

#################################################################################################
# Remaining setup
mkdir $HOME/WGS_Data
mkdir $HOME/WGS_Results

# pip install lxml # Requirement already satisfied
# pip install scipy # Requirement already satisfied

###################################################################################################
# Retrieve test isolates
cd $HOME

fetcher --url="https://github.com/APHA-BAC/NextflowSerotypingPipeline/tree/circleci-project-setup/test_isolates" --out=$HOME/WGS_Data

cd $HOME/nextflow
# Updated this next line, which wasn't working:
fetcher --url="https://github.com/APHA-BAC/NextflowSerotypingPipeline/blob/master/SCE3_pipeline_update.nf" --out=$HOME/nextflow

