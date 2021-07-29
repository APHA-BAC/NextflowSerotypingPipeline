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
RUN sudo apt-get -y install python3

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



