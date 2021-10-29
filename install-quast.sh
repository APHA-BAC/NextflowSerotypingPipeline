###########################################################################################
# Quast
cd $HOME
set -e

# TODO: don't download from master branch, download directly from the authour

# I have changed this to follow the style of other installations into /opt
wget --no-check-certificate --content-disposition https://raw.githubusercontent.com/APHA-BAC/NextflowSerotypingPipeline/master/pipeline_component_software/quast/quast-5.1.0rc1.tar.gz
tar -xzf quast-5.1.0rc1.tar.gz

sudo cp -r quast-5.1.0rc1 /opt
sudo /opt/quast-5.1.0rc1/install.sh

sudo ln -s /opt/quast-5.1.0rc1/quast.py /usr/local/bin
rm -rf quast-5.1.0rc1
sudo rm -rf quast_test_output
rm quast-5.1.0rc1.tar.gz
