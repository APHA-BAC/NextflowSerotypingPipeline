FROM ubuntu:20.04

################## METADATA ###########################

LABEL base.image=ubuntu:20.04
LABEL software="Salmonella Serotyping Pipeline Image"
LABEL about.summary="Salmonella Serotyping Pipeline"
LABEL about.documentation="https://github.com/APHA-BAC/NextflowSerotypingPipeline"
LABEL about.tags="Genomics, WGS"


################## DEPENDENCIES ######################

# Copy repository
WORKDIR /root/NextflowSerotypingPipeline/
COPY ./install/ ./install/

# RETURN
RUN bash install/install.bash

COPY ./install-fastp.sh ./install/install-fastp.sh
RUN bash install/install-fastp.sh

COPY ./install-seqtk.sh ./install/install-seqtk.sh
RUN bash install/install-seqtk.sh

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

COPY ./install-seqsero2.sh ./install/install-seqsero2.sh
RUN bash install/install-seqsero2.sh

COPY ./install-sistr.sh ./install/install-sistr.sh
RUN bash install/install-sistr.sh

COPY ./install-srst2.sh ./install/install-srst2.sh
RUN bash install/install-srst2.sh

COPY ./install-nextflow.sh ./install/install-nextflow.sh
RUN bash install/install-nextflow.sh

COPY ./install-sra-toolkit.sh ./install/install-sra-toolkit.sh
RUN bash install/install-sra-toolkit.sh

COPY ./setup-environment.sh ./install/setup-environment.sh
RUN bash install/setup-environment.sh

COPY ./ ./
