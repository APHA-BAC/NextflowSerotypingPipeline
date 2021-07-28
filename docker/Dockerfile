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

#RUN INSTALL_PATH=~/anaconda \
#    && wget https://repo.anaconda.com/miniconda/Miniconda3-py37_4.8.3-Linux-x86_64.sh \
#    && bash Miniconda3-py37_4.8.3-Linux-x86_64.sh -fbp $INSTALL_PATH \
#    && PATH=$INSTALL_PATH/bin:$PATH \
#ENV PATH="/root/anaconda/bin:${PATH}"

RUN wget https://repo.anaconda.com/miniconda/Miniconda3-py37_4.8.3-Linux-x86_64.sh

RUN bash Miniconda3-py37_4.8.3-Linux-x86_64.sh -b && \
    echo "export PATH="/root/anaconda3/bin:$PATH"" >> ~/.bashrc && \
    /bin/bash -c "source ~/.bashrc"
ENV PATH /root/anaconda3/bin:$PATH
RUN conda update --all



#RUN PATH=/root/anaconda/bin:$PATH

#RUN conda create --name shovill

#RUN conda init bash

#ENV PATH=/root/anaconda/bin:$PATH

#RUN sudo /usr/sbin/shutdown -r now
#RUN conda create --name shovill
#RUN conda install -c conda-forge bazel
