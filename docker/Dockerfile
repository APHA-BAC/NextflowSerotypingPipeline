FROM ubuntu:20.04

################## METADATA ##########################

LABEL base.image=ubuntu:20.04
LABEL software="Bovine-TB Pipeline Image"
LABEL about.summary="Bioinformatics Pipeline for post-processing of Bovine TB fastq reads"
LABEL about.documentation="https://github.com/APHA-CSU/BovTB-nf"
LABEL about.tags="Genomics, WGS"

RUN apt-get update
RUN apt-get -y install git