set -e

# Apt-get
sudo apt-get -y update && DEBIAN_FRONTEND=noninteractive apt-get install -y \
    openjdk-11-jdk \
    wget \
    make \
    git \
    curl \
    libz-dev \
    libbz2-dev \
    libncurses5-dev \
    libncursesw5-dev \
    libghc-bzlib-prof \
    gcc \
    unzip \
    zlib1g-dev

sudo apt-get install -y libcurl4-openssl-dev
sudo apt-get install -y python3 \
    python3-numpy \
    python3-pip \
    vim \
    nano \
    bc

# python 
# TODO: Version pin pandas to 1.15
pip3 install biopython pandas gitpython
sudo ln -s /usr/bin/python3 /usr/bin/python

# Install the dependencies
# sh install-fastqc.sh