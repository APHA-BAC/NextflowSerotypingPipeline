# KmerID
set -e
cd $HOME

# TODO: download directly from authour
wget --no-check-certificate --content-disposition https://raw.githubusercontent.com/APHA-BAC/NextflowSerotypingPipeline/master/pipeline_component_software/KmerID/kmerid.zip
unzip kmerid.zip
sudo cp -r kmerid /opt
sudo make -C /opt/kmerid all # Changed this so no need to change directory
sudo chmod +x /opt/kmerid/setup_refs.py
sudo chmod +x /opt/kmerid/kmerid_python3.py
rm -rf kmerid
rm kmerid.zip

sudo apt-get -y install nodejs
sudo apt-get -y install npm
sudo npm install -g github-files-fetcher

# TODO: download from somewhere that's not this repo's branch 
# fetcher --url="https://github.com/APHA-BAC/NextflowSerotypingPipeline/tree/circleci-project-setup/kmerid_ref/Acinetobacter" --out=/opt/kmerid/ref
sudo fetcher --url="https://github.com/APHA-BAC/NextflowSerotypingPipeline/tree/circleci-project-setup/kmerid_ref/Citrobacter" --out=/opt/kmerid/ref
sudo fetcher --url="https://github.com/APHA-BAC/NextflowSerotypingPipeline/tree/circleci-project-setup/kmerid_ref/Salmonella" --out=/opt/kmerid/ref
