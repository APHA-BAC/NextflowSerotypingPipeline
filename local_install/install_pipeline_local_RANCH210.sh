#!/usr/bin/bash
# LOCAL INSTALL SCRIPT FOR RANCH-210

####################################################################################
# The following section is primarily for blank VMs (or Docker) only
# Skip to next section if installing pipeline locally on Ranch-210
cd $HOME

sudo apt-get update
# sudo apt-get install -y -q sudo
# sudo apt-get install -y -q apt-utils
# sudo apt-get install -y -q make
# sudo apt-get install -y -q gcc

# sudo apt-get install -y -q git
# sudo apt-get install -y -q wget
# sudo apt-get install -y -q curl

# sudo apt-get install -y -q libz-dev
# sudo apt-get install -y -q libbz2-dev
# sudo apt-get install -y -q libncurses5-dev
# sudo apt-get install -y -q libncursesw5-dev

# sudo apt-get install -y -q libghc-bzlib-prof
# sudo apt-get install -y -q zlib1g-dev
# sudo apt-get install -y -q libcurl4-openssl-dev

# sudo apt-get install -y -q libfile-copy-recursive-perl
# sudo apt-get install -y -q libio-socket-ssl-perl
# sudo apt-get install -y -q libio-tee-perl
# sudo apt-get install -y -q libunicode-string-perl

# sudo apt-get install -y -q nano

#########################################################################################
# Miniconda
cd $HOME

wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
chmod +x Miniconda3-latest-Linux-x86_64.sh
./Miniconda3-latest-Linux-x86_64.sh
rm Miniconda3-latest-Linux-x86_64.sh

source ~/.bashrc

# Conda environment
# Now this is where things change. I create a specific Conda environment for the pipeline
# Any other projects on the SCE3 machine will thus be unaffected by conda and pip installs within this environment
conda create -n salmpipe python=3.7
conda activate salmpipe

# Matplotlib
pip install matplotlib

# SPECIFIC VERSION, OLD INSTALL, DO NOT RUN
# wget https://repo.anaconda.com/miniconda/Miniconda3-py37_4.8.3-Linux-x86_64.sh

######################################################################################
# fastp
# cd $HOME

# wget --no-check-certificate --content-disposition http://opengene.org/fastp/fastp.0.23.1
# mv fastp.0.23.1 fastp
# sudo cp fastp /opt
# sudo chmod +x /opt/fastp
# sudo ln -s /opt/fastp /usr/local/bin
# rm fastp


######################################################################################
# FastQC
# cd $HOME

# # sudo apt-get -y install unzip
# # sudo apt-get -y install openjdk-11-jdk
# # sudo apt-get -y install default-jre

# wget --no-check-certificate https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.9.zip
# unzip fastqc_v0.11.9.zip
# sudo cp -r FastQC /opt
# sudo chmod +x /opt/FastQC/fastqc
# sudo ln -s /opt/FastQC/fastqc /usr/local/bin
# rm -r FastQC
# rm fastqc_v0.11.9.zip


######################################################################################
# seqtk
cd $HOME

conda install -c bioconda seqtk


###########################################################################################
# Shovill
cd $HOME

conda install -c conda-forge -c bioconda -c defaults shovill=0.9.0
shovill --check


###########################################################################################
# Quast
# cd $HOME

# wget --no-check-certificate --content-disposition https://github.com/ablab/quast/releases/download/quast_5.1.0rc1/quast-5.1.0rc1.tar.gz
# tar -xzvf quast-5.1.0rc1.tar.gz

# sudo cp -r quast-5.1.0rc1 /opt
# sudo /opt/quast-5.1.0rc1/install.sh

# sudo ln -s /opt/quast-5.1.0rc1/quast.py /usr/local/bin
# rm -r quast-5.1.0rc1
# sudo rm -rf quast_test_output
# rm quast-5.1.0rc1.tar.gz


###########################################################################################
# MOST
cd $HOME

# sudo apt-get -y install libncurses5

conda install -c bioconda emboss
conda install -c kantorlab blastn
pip install lxml==4.5.2
pip install biopython==1.73

# wget --no-check-certificate --content-disposition https://raw.githubusercontent.com/APHA-BAC/MOST-PHE/main/pipeline_component_software/most/most.zip
# unzip most.zip
# sudo cp -r most /opt
# sudo chmod +x /opt/most/MOST-master/MOST.py

# rm -r most
# rm most.zip


#############################################################################
# KmerID
# cd $HOME

# wget --no-check-certificate --content-disposition https://raw.githubusercontent.com/JaromirGuzinski/NextflowSerotypingPipeline/master/pipeline_component_software/KmerID/kmerid.zip
# unzip kmerid.zip
# sudo cp -r kmerid /opt

# aws s3 cp s3://s3-staging-area/kmerid_ref/kmerid_ref.zip kmerid_ref.zip
# jar xvf kmerid_ref.zip # Unzip can fail for unknown reasons; this is a workaround
# sudo mv ref /opt/kmerid

# sudo make -C /opt/kmerid all # No need to change directory when the -C option is given
# sudo chmod +x /opt/kmerid/setup_refs.py
# sudo chmod +x /opt/kmerid/kmerid_python3.py

# rm -r kmerid
# rm kmerid.zip
# rm kmerid_ref.zip

# sudo apt-get -y install nodejs
# sudo apt-get -y install npm
# sudo npm install -g github-files-fetcher


###########################################################################################
# SeqSero2
conda install -c bioconda seqsero2=1.1.1


############################################################################################
# Sistr
pip install numpy

conda config --add channels conda-forge
conda config --add channels r
conda config --add channels bioconda
conda install sistr_cmd

pip install pandas==1.0.5

#############################################################################################
# Srst2
# cd $HOME

# wget --no-check-certificate --content-disposition https://raw.githubusercontent.com/APHA-BAC/SRST2/master/pipeline_component_software/srst2/srst2.zip
# unzip srst2.zip
# sudo cp -r srst2 /opt
# sudo chmod +x /opt/srst2/scripts/srst2.py
# sudo ln -s /opt/srst2/scripts/srst2.py /usr/local/bin
# rm -rf srst2
# rm srst2.zip


#################################################################################################
# Nextflow
mkdir $HOME/nextflow
cd $HOME/nextflow
curl -L0 https://github.com/nextflow-io/nextflow/releases/download/v21.04.1/nextflow-21.04.1-all | bash
# This ^ gives a warning, "curl: (23) Failed writing body (4096 != 16384)", but seems to work?


#################################################################################################
# Remaining setup
mkdir $HOME/WGS_Data
mkdir $HOME/WGS_Results

#################################################################################################
# Summary table script
mkdir $HOME/summary
cd $HOME/summary

wget --no-check-certificate --content-disposition https://raw.githubusercontent.com/APHA-BAC/NextflowSerotypingPipeline/master/summaryTable_reworked.py

fetcher --url="https://github.com/APHA-BAC/NextflowSerotypingPipeline/blob/master/local_install/SCE3_pipeline_update_LOCAL.nf" --out=$HOME/nextflow
