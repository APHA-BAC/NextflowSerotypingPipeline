#################################################################################################
set -e

# Nextflow
mkdir $HOME/nextflow
cd $HOME/nextflow
curl -L0 https://github.com/nextflow-io/nextflow/releases/download/v21.04.1/nextflow-21.04.1-all | bash
# This ^ gives a warning, "curl: (23) Failed writing body (4096 != 16384)", but seems to work?

ln -s $HOME/nextflow /usr/local/bin/nextflow

fetcher --url="https://github.com/APHA-BAC/NextflowSerotypingPipeline/blob/seqtk_subsampling/SCE3_pipeline_update.nf" --out=$HOME/nextflow
# fetcher --url="https://github.com/APHA-BAC/NextflowSerotypingPipeline/blob/master/SCE3_pipeline_update.nf" --out=$HOME/nextflow
