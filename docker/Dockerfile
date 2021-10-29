FROM ubuntu:20.04

################## METADATA ###########################

LABEL base.image=ubuntu:20.04
LABEL software="Salmonella Serotyping Pipeline Image"
LABEL about.summary="Salmonella Serotyping Pipeline"
LABEL about.documentation="https://github.com/APHA-BAC/NextflowSerotypingPipeline"
LABEL about.tags="Genomics, WGS"


################## DEPENDENCIES ######################

# Copy repository
WORKDIR /root/pipeline/
COPY ./install/ ./install/

# RETURN
RUN bash install/install.bash

COPY ./install-fastqc.sh ./install/install-fastqc.sh
RUN bash install/install-fastqc.sh

ENV PATH=/opt/conda/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin
COPY ./install-conda.sh ./install/install-conda.sh
RUN bash install/install-conda.sh

COPY ./install-shovill.sh ./install/install-shovill.sh
RUN bash install/install-shovill.sh

COPY ./install-quast.sh ./install/install-quast.sh
RUN bash install/install-quast.sh

COPY ./install-most.sh ./install/install-most.sh
RUN bash install/install-most.sh

COPY ./install-kmerid.sh ./install/install-kmerid.sh
RUN bash install/install-kmerid.sh

#TODO set -e

# # Sudo
# RUN apt-get -y update 
# RUN DEBIAN_FRONTEND=noninteractive apt-get install -y -q sudo
# RUN DEBIAN_FRONTEND=noninteractive apt-get install -y -q apt-utils
# RUN DEBIAN_FRONTEND=noninteractive apt-get install -y -q wget
# RUN DEBIAN_FRONTEND=noninteractive apt-get install -y -q git
# RUN DEBIAN_FRONTEND=noninteractive apt-get install -y -q curl
# RUN DEBIAN_FRONTEND=noninteractive apt-get install -y -q make
# RUN DEBIAN_FRONTEND=noninteractive apt-get install -y -q openjdk-11-jdk
# RUN DEBIAN_FRONTEND=noninteractive apt-get install -y -q libz-dev
# RUN DEBIAN_FRONTEND=noninteractive apt-get install -y -q libbz2-dev
# RUN DEBIAN_FRONTEND=noninteractive apt-get install -y -q libncurses5-dev
# RUN DEBIAN_FRONTEND=noninteractive apt-get install -y -q libncursesw5-dev
# RUN DEBIAN_FRONTEND=noninteractive apt-get install -y -q libghc-bzlib-prof
# RUN DEBIAN_FRONTEND=noninteractive apt-get install -y -q gcc
# RUN DEBIAN_FRONTEND=noninteractive apt-get install -y -q zlib1g-dev
# RUN DEBIAN_FRONTEND=noninteractive apt-get install -y -q libcurl4-openssl-dev
# RUN DEBIAN_FRONTEND=noninteractive apt-get install -y -q python3
# RUN DEBIAN_FRONTEND=noninteractive apt-get install -y -q default-jre
# RUN DEBIAN_FRONTEND=noninteractive apt-get install -y -q libfile-copy-recursive-perl
# RUN DEBIAN_FRONTEND=noninteractive apt-get install -y -q libio-socket-ssl-perl
# RUN DEBIAN_FRONTEND=noninteractive apt-get install -y -q libio-tee-perl
# RUN DEBIAN_FRONTEND=noninteractive apt-get install -y -q libunicode-string-perl
# RUN DEBIAN_FRONTEND=noninteractive apt-get install -y -q nano
# RUN DEBIAN_FRONTEND=noninteractive apt-get install -y -q gzip
# RUN DEBIAN_FRONTEND=noninteractive apt-get install -y -q rename

# # Pipeline
# WORKDIR /NextflowSerotypingPipeline/
# COPY ./ ./

# # FastQC
# RUN sh install-fastqc.sh

# # Shovill
# RUN wget https://repo.anaconda.com/miniconda/Miniconda3-py37_4.8.3-Linux-x86_64.sh  && \
#      sudo chmod +x Miniconda3-py37_4.8.3-Linux-x86_64.sh && \
#      ./Miniconda3-py37_4.8.3-Linux-x86_64.sh -b -p /opt/conda && \     
#      rm Miniconda3-py37_4.8.3-Linux-x86_64.sh

# # Add conda binaries to PATH
# ENV PATH=/opt/conda/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin

# RUN /opt/conda/bin/conda init bash && \
#       . ~/.bashrc && \
#       /opt/conda/bin/conda install -c conda-forge -c bioconda -c defaults shovill=0.9.0 && \
#       /opt/conda/bin/shovill --check

# # Quast
# RUN sudo ln -s /usr/bin/python3 /usr/bin/python
# RUN sudo apt-get -y install python3-pip
# RUN sudo python -m pip install -U pip
# RUN sudo python -m pip install -U matplotlib

# RUN wget --no-check-certificate --content-disposition https://raw.githubusercontent.com/APHA-BAC/NextflowSerotypingPipeline/master/pipeline_component_software/quast/quast-5.1.0rc1.tar.gz
# RUN tar -xvf quast-5.1.0rc1.tar.gz
# RUN rm quast-5.1.0rc1.tar.gz

# RUN ./quast-5.1.0rc1/install.sh
# RUN sudo ln -s /quast-5.1.0rc1/quast.py /usr/local/bin

# # MOST 
# RUN wget --no-check-certificate --content-disposition https://raw.githubusercontent.com/APHA-BAC/NextflowSerotypingPipeline/master/pipeline_component_software/most/most.zip
# RUN unzip most.zip
# RUN sudo cp -r most /opt
# RUN sudo chmod +x /opt/most/MOST-master/MOST.py
# RUN rm -rf most
# RUN rm most.zip

# RUN pip install lxml==4.5.2
# RUN pip install biopython==1.73
# RUN /opt/conda/bin/conda install -c bioconda emboss
# RUN /opt/conda/bin/conda install -c kantorlab blastn
# RUN sudo apt-get -y install libncurses5

# # Kmerid
# RUN wget --no-check-certificate --content-disposition https://raw.githubusercontent.com/APHA-BAC/NextflowSerotypingPipeline/master/pipeline_component_software/KmerID/kmerid.zip
# RUN sudo unzip kmerid.zip -d /opt/ 
# RUN rm kmerid.zip

# RUN cd /opt/kmerid && \
#      sudo make all
# RUN sudo chmod +x /opt/kmerid/setup_refs.py
# RUN sudo chmod +x /opt/kmerid/kmerid_python3.py

# RUN sudo apt-get -y install nodejs
# RUN sudo apt-get -y install npm
# RUN npm install -g github-files-fetcher
# RUN fetcher --url="https://github.com/APHA-BAC/NextflowSerotypingPipeline/tree/master/kmerid_ref/Citrobacter" --out=/opt/kmerid/ref
# RUN fetcher --url="https://github.com/APHA-BAC/NextflowSerotypingPipeline/tree/master/kmerid_ref/Salmonella" --out=/opt/kmerid/ref

# # SeqSero2
# RUN /opt/conda/bin/conda install -c bioconda seqsero2=1.1.1

# # Sistr
# RUN pip install --upgrade pip
# RUN pip install wheel
# RUN pip install numpy pandas==1.1.5

# RUN /opt/conda/bin/conda config --add channels defaults
# RUN /opt/conda/bin/conda config --add channels conda-forge
# RUN /opt/conda/bin/conda config --add channels r
# RUN /opt/conda/bin/conda config --add channels bioconda
# RUN /opt/conda/bin/conda install sistr_cmd

# # Srst2
# RUN wget --no-check-certificate --content-disposition https://raw.githubusercontent.com/APHA-BAC/NextflowSerotypingPipeline/master/pipeline_component_software/srst2/srst2.zip
# RUN unzip srst2.zip
# RUN sudo cp -r srst2 /opt
# RUN sudo chmod +x /opt/srst2/scripts/srst2.py
# RUN sudo ln -s /opt/srst2/scripts/srst2.py /usr/local/bin
# RUN rm -rf srst2
# RUN rm srst2.zip

# RUN pip install scipy

# # Nextflow
# RUN mkdir $HOME/nextflow
# RUN cd $HOME/nextflow
# RUN curl -L0 https://github.com/nextflow-io/nextflow/releases/download/v21.04.1/nextflow-21.04.1-all | bash

# RUN mkdir $HOME/summary
# RUN cd $HOME/summary
# RUN wget --no-check-certificate --content-disposition https://raw.githubusercontent.com/APHA-BAC/NextflowSerotypingPipeline/master/summaryTable_Python3.py

# RUN mkdir $HOME/WGS_Data
# RUN mkdir $HOME/WGS_Results

# RUN pip install lxml

# # Pipeline run
# #RUN fetcher --url="https://github.com/APHA-BAC/NextflowSerotypingPipeline/tree/circleci-project-setup/test_isolates" --out=$HOME/WGS_Data
# RUN cd $HOME/nextflow

# RUN pip install lxml
