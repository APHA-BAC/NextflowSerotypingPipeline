###########################################################################################
# Shovill
set -e

conda install -c conda-forge -c bioconda -c defaults shovill
shovill --check
