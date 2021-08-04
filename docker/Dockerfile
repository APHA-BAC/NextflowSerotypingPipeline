FROM ubuntu:20.04

################## METADATA ##########################

LABEL base.image=ubuntu:20.04
LABEL software="Salmonella Serotyping Pipeline Image"
LABEL about.summary="Salmonella Serotyping Pipeline"
LABEL about.documentation="https://github.com/APHA-BAC/NextflowSerotypingPipeline"
LABEL about.tags="Genomics, WGS"

# Sudo
RUN apt-get -y update && DEBIAN_FRONTEND=noninteractive apt-get install -y sudo
RUN sudo apt-get -y install apt-utils
RUN sudo apt-get -y install wget
RUN sudo apt-get -y install curl
RUN sudo apt-get -y install make
RUN sudo apt-get -y install openjdk-11-jdk
RUN sudo apt-get -y install libz-dev
RUN sudo apt-get -y install libbz2-dev
RUN sudo apt-get -y install libncurses5-dev
RUN sudo apt-get -y install libncursesw5-dev
RUN sudo apt-get -y install libghc-bzlib-prof
RUN sudo apt-get -y install gcc
RUN sudo apt-get -y install zlib1g-dev
RUN sudo apt-get -y install libcurl4-openssl-dev
RUN sudo apt-get -y install python3
RUN sudo apt-get -y install default-jre
RUN sudo apt-get -y install default-jdk

# Pipeline
WORKDIR /NextflowSerotypingPipeline/
COPY ./ ./

# FastQC
RUN sh install-fastqc.sh

# Shovill
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-py37_4.8.3-Linux-x86_64.sh  && \
     sudo chmod +x Miniconda3-py37_4.8.3-Linux-x86_64.sh && \
     ./Miniconda3-py37_4.8.3-Linux-x86_64.sh -b -p /opt/conda && \     
     rm Miniconda3-py37_4.8.3-Linux-x86_64.sh

RUN /opt/conda/bin/conda install -c conda-forge -c bioconda -c defaults shovill=0.9.0
RUN /opt/conda/bin/shovill --check

# Quast
RUN sudo ln -s /usr/bin/python3 /usr/bin/python
RUN sudo apt-get -y install python3-pip
RUN sudo python -m pip install -U pip
RUN sudo python -m pip install -U matplotlib

RUN wget --no-check-certificate --content-disposition https://raw.githubusercontent.com/APHA-BAC/NextflowSerotypingPipeline/master/pipeline_component_software/quast/quast-5.1.0rc1.tar.gz
RUN tar -xvf quast-5.1.0rc1.tar.gz
RUN rm quast-5.1.0rc1.tar.gz

RUN ./quast-5.1.0rc1/install.sh
RUN sudo ln -s /quast-5.1.0rc1/quast.py /usr/local/bin

# MOST 
RUN wget --no-check-certificate --content-disposition https://raw.githubusercontent.com/APHA-BAC/NextflowSerotypingPipeline/master/pipeline_component_software/most/most.zip
RUN unzip most.zip
RUN sudo cp -r most /opt
RUN sudo chmod +x /opt/most/MOST-master/MOST.py
RUN rm -rf most
RUN rm most.zip

RUN pip install lxml==4.5.2
RUN pip install biopython==1.73
RUN /opt/conda/bin/conda install -c bioconda emboss
RUN /opt/conda/bin/conda install -c kantorlab blastn
RUN sudo apt-get -y install libncurses5

# Kmerid
RUN wget --no-check-certificate --content-disposition https://raw.githubusercontent.com/APHA-BAC/NextflowSerotypingPipeline/master/pipeline_component_software/KmerID/kmerid.zip
RUN sudo unzip kmerid.zip -d /opt/ 
RUN rm kmerid.zip

RUN cd /opt/kmerid && \
     sudo make all
RUN sudo chmod +x /opt/kmerid/setup_refs.py
RUN sudo chmod +x /opt/kmerid/kmerid_python3.py

RUN sudo apt-get -y install nodejs
RUN sudo apt-get -y install npm
RUN npm install -g github-files-fetcher
RUN fetcher --url="https://github.com/APHA-BAC/NextflowSerotypingPipeline/tree/circleci-project-setup/kmerid_ref/Acinetobacter" --out=/opt/kmerid/ref
RUN fetcher --url="https://github.com/APHA-BAC/NextflowSerotypingPipeline/tree/circleci-project-setup/kmerid_ref/Salmonella" --out=/opt/kmerid/ref

# SeqSero2
RUN /opt/conda/bin/conda install -c bioconda seqsero2=1.1.1

# Sistr
RUN pip install --upgrade pip
RUN pip install wheel
RUN pip install numpy pandas==1.1.5

RUN /opt/conda/bin/conda config --add channels defaults
RUN /opt/conda/bin/conda config --add channels conda-forge
RUN /opt/conda/bin/conda config --add channels r
RUN /opt/conda/bin/conda config --add channels bioconda
RUN /opt/conda/bin/conda install sistr_cmd

# Srst2
RUN wget --no-check-certificate --content-disposition https://raw.githubusercontent.com/APHA-BAC/NextflowSerotypingPipeline/master/pipeline_component_software/srst2/srst2.zip
RUN unzip srst2.zip
RUN sudo cp -r srst2 /opt
RUN sudo chmod +x /opt/srst2/scripts/srst2.py
RUN sudo ln -s /opt/srst2/scripts/srst2.py /usr/local/bin
RUN rm -rf srst2
RUN rm srst2.zip

RUN pip install scipy

# Nextflow
RUN mkdir /home/$USER/nextflow
RUN cd /home/$USER/nextflow
RUN curl -s https://get.nextflow.io | bash
RUN fetcher --url="https://github.com/APHA-BAC/NextflowSerotypingPipeline/master/SCE3_pipeline_update.nf" --out=/home/$USER/nextflow

RUN mkdir /home/summary
RUN cd /home/summary
RUN wget --no-check-certificate --content-disposition https://raw.githubusercontent.com/APHA-BAC/NextflowSerotypingPipeline/master/summaryTable_Python3.py

RUN mkdir /home/WGS_Data
RUN mkdir /home/WGS_Results

# Pipeline run
RUN fetcher --url="https://github.com/APHA-BAC/NextflowSerotypingPipeline/tree/circleci-project-setup/test_isolates" --out=/home/WGS_Data
RUN cd /home/$USER/nextflow
RUN ./nextflow run SCE3_pipeline_update.nf
