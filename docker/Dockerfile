FROM ubuntu:20.04

################## METADATA ##########################

LABEL base.image=ubuntu:20.04
LABEL software="Salmonella Serotyping Pipeline"
LABEL about.summary="Bioinformatics Pipeline for post-processing of Bovine TB fastq reads"
LABEL about.documentation="https://github.com/APHA-CSU/BovTB-nf"
LABEL about.tags="Genomics, WGS"

# Sudo
RUN apt-get -y update && DEBIAN_FRONTEND=noninteractive apt-get install -y sudo

WORKDIR /NextflowSerotypingPipeline/
COPY ./ ./
RUN sudo bash install.bash
