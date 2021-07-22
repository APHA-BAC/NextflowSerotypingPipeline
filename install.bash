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

# Install the dependencies
# sh install-fastqc.sh