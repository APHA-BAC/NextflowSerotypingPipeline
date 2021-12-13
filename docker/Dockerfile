FROM ubuntu:20.04

################## METADATA ###########################

LABEL base.image=ubuntu:20.04
LABEL software="Salmonella Serotyping Pipeline Image"
LABEL about.summary="Salmonella Serotyping Pipeline"
LABEL about.documentation="https://github.com/APHA-BAC/NextflowSerotypingPipeline"
LABEL about.tags="Genomics, WGS"


################## ARGS #############################

ARG SALMO_PATH="/NextflowSerotypingPipeline/"


################## DEPENDENCIES ######################
# Copy repository
WORKDIR $SALMO_PATH
COPY ./ ./
#RUN chmod +x ./bin/*

# RETURN
RUN bash install/install.bash



#COPY ./ ./
