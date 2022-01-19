###########################################################################################
# Quast
set -e

cd $HOME

# I have changed this to follow the style of other installations into /opt
wget --no-check-certificate --content-disposition https://github.com/ablab/quast/releases/download/quast-5.1.0rc1/quast-5.1.0rc1.tar.gz
tar -xzf quast-5.1.0rc1.tar.gz

sudo cp -r quast-5.1.0rc1 /opt
sudo /opt/quast-5.1.0rc1/install.sh

sudo ln -s /opt/quast-5.1.0rc1/quast.py /usr/local/bin
rm -rf quast-5.1.0rc1
sudo rm -rf quast_test_output
rm quast-5.1.0rc1.tar.gz
