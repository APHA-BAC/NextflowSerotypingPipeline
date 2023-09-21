# MOST
set -e
cd $HOME

sudo apt-get -y install libncurses5
conda install -c bioconda emboss
conda install -c kantorlab blastn
pip install lxml==4.5.2
pip install biopython==1.73

# TODO: download directly from authour
wget --no-check-certificate --content-disposition https://raw.githubusercontent.com/APHA-BAC/MOST-PHE/most_consensus_change/pipeline_component_software/most/most.zip
unzip most.zip
sudo cp -r most /opt
sudo chmod +x /opt/most/MOST-master/MOST.py
rm -rf most
rm most.zip
