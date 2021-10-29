#################################################################################################
set -e

# Nextflow
mkdir $HOME/nextflow
cd $HOME/nextflow
curl -L0 https://github.com/nextflow-io/nextflow/releases/download/v21.04.1/nextflow-21.04.1-all | bash
# This ^ gives a warning, "curl: (23) Failed writing body (4096 != 16384)", but seems to work?

ln -s $PWD/nextflow /usr/local/bin/nextflow