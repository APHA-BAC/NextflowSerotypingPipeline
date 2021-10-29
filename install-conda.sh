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

# Link to path
sudo ln -s /usr/bin/python3 /usr/bin/python

# Matplotlib
pip install matplotlib
