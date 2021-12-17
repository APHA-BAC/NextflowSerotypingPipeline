set -e

# Apt Packages
apt-get -y update
DEBIAN_FRONTEND=noninteractive apt-get install -y -q \
    sudo \
    apt-utils \
    make \
    gcc \
    bc \
    git \
    wget \
    curl \
    libz-dev \
    libbz2-dev \
    libncurses5-dev \
    libncursesw5-dev \
    libghc-bzlib-prof \
    zlib1g-dev \
    libcurl4-openssl-dev \
    libfile-copy-recursive-perl \
    libio-socket-ssl-perl \
    libio-tee-perl \
    libunicode-string-perl \
    nano \
    python3

# Bioinformatics Tools
bash install/install-conda.sh
bash install/install-fastp.sh
bash install/install-fastqc.sh
bash install/install-seqtk.sh
bash install/install-shovill.sh
bash install/install-quast.sh
bash install/install-most.sh
bash install/install-kmerid.sh
bash install/install-seqsero2.sh
bash install/install-sistr.sh
bash install/install-srst2.sh
bash install/install-nextflow.sh
bash install/install-sra-toolkit.sh
bash install/setup-environment.sh
