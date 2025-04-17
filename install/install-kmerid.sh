# KmerID
set -e
cd $HOME

wget --no-check-certificate --content-disposition https://raw.githubusercontent.com/APHA-BAC/KMERID-PHE/master/pipeline_component_software/KmerID/kmerid.zip
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

# sudo fetcher --url="https://github.com/APHA-BAC/KMERID-PHE/tree/master/kmerid_ref/Citrobacter" --out=/opt/kmerid/ref
# sudo fetcher --url="https://github.com/APHA-BAC/KMERID-PHE/tree/master/kmerid_ref/Salmonella" --out=/opt/kmerid/ref


curl "https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip" -o "awscliv2.zip"
unzip awscliv2.zip
sudo ./aws/install

# aws s3 cp --recursive s3://s3-ranch-046/KmerID_Ref_Genomes/ /opt/kmerid/
# aws s3 cp --recursive s3://s3-ranch-046/KmerID_Ref_Genomes/ref/Salmonella/ /opt/kmerid/ref/Salmonella
