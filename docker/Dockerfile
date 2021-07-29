FROM ubuntu:20.04

################## METADATA ##########################

LABEL base.image=ubuntu:20.04
LABEL software="Salmonella Serotyping Pipeline Image"
LABEL about.summary="Salmonella Serotyping Pipeline"
LABEL about.documentation="https://github.com/APHA-BAC/NextflowSerotypingPipeline"
LABEL about.tags="Genomics, WGS"

#Sudo
RUN apt-get -y update && DEBIAN_FRONTEND=noninteractive apt-get install -y sudo

RUN sudo apt-get -y install apt-utils
RUN sudo apt-get -y install wget


# Pipeline
WORKDIR /NextflowSerotypingPipeline/
COPY ./ ./
RUN chmod +x ./bin/*

RUN sh install-fastqc.sh

RUN wget https://repo.anaconda.com/miniconda/Miniconda3-py37_4.8.3-Linux-x86_64.sh \
&& sudo chmod +x Miniconda3-py37_4.8.3-Linux-x86_64.sh \
&& ./Miniconda3-py37_4.8.3-Linux-x86_64.sh -b -p /miniconda \
&& rm Miniconda3-py37_4.8.3-Linux-x86_64.sh

RUN conda --version

#RUN wget https://repo.anaconda.com/miniconda/Miniconda3-py37_4.8.3-Linux-x86_64.sh

#RUN bash Miniconda3-py37_4.8.3-Linux-x86_64.sh -b && \
#    echo "export PATH="/root/anaconda3/bin:$PATH"" >> ~/.bashrc && \
#    /bin/bash -c "source ~/.bashrc"
#ENV PATH /root/anaconda3/bin:$PATH

#RUN sudo apt-get -y update

#RUN conda --version



#RUN PATH=/root/anaconda/bin:$PATH

#RUN conda create --name shovill

#RUN conda init bash

#ENV PATH=/root/anaconda/bin:$PATH

#RUN sudo /usr/sbin/shutdown -r now
#RUN conda create --name shovill
#RUN conda install -c conda-forge bazel
