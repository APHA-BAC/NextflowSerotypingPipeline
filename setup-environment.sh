set -e

#################################################################################################
# Remaining setup
mkdir $HOME/WGS_Data
mkdir $HOME/WGS_Results

###################################################################################################
# Retrieve test isolates
cd $HOME

fetcher --url="https://github.com/APHA-BAC/NextflowSerotypingPipeline/tree/circleci-project-setup/test_isolates" --out=$HOME/WGS_Data

cd $HOME/nextflow
# Updated this next line, which wasn't working:
fetcher --url="https://github.com/APHA-BAC/NextflowSerotypingPipeline/blob/master/SCE3_pipeline_update.nf" --out=$HOME/nextflow

