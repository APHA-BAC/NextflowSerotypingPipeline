set -e

#################################################################################################
# Remaining setup
mkdir $HOME/WGS_Data
mkdir $HOME/WGS_Results

#################################################################################################
# Summary table script
mkdir $HOME/summary # Changed to $HOME for SC3 - We don't want summary installed in /home/summary
cd $HOME/summary
wget --no-check-certificate --content-disposition https://raw.githubusercontent.com/APHA-BAC/NextflowSerotypingPipeline/master/summaryTable_reworked.py
