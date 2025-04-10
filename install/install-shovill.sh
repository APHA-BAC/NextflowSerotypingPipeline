###########################################################################################
# Shovill
set -e

# conda install -c conda-forge -c bioconda -c defaults shovill
conda install -c "bioconda/label/main" shovill
shovill --check
