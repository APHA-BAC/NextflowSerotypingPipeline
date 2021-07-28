FROM ubuntu:20.04

################## METADATA ##########################

LABEL base.image=ubuntu:20.04
LABEL software="Salmonella Serotyping Pipeline Image"
LABEL about.summary="Salmonella Serotyping Pipeline"
LABEL about.documentation="https://github.com/APHA-BAC/NextflowSerotypingPipeline"
LABEL about.tags="Genomics, WGS"

# Sudo
RUN apt-get -y update && DEBIAN_FRONTEND=noninteractive apt-get install -y sudo

WORKDIR /NextflowSerotypingPipeline/
COPY ./ ./
RUN sudo apt-get -y install apt-utils
RUN sudo apt-get -y install wget
RUN sh install-fastqc.sh
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-py37_4.8.3-Linux-x86_64.sh
RUN sudo chmod +x Miniconda3-py37_4.8.3-Linux-x86_64.sh
RUN sh Miniconda3-py37_4.8.3-Linux-x86_64.sh -b
RUN conda create --name shovill
#RUN conda install -c conda-forge bazel
