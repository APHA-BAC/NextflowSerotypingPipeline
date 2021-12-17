set -e

# TODO: DRY up the apt-get so it's in one call

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


bash -e install/install-conda.sh
bash -e install/install-fastp.sh
bash -e install/install-fastqc.sh
bash -e install/install-seqtk.sh
bash -e install/install-shovill.sh
bash -e install/install-quast.sh
bash -e install/install-most.sh
bash -e install/install-kmerid.sh
bash -e install/install-seqsero2.sh
bash -e install/install-sistr.sh
bash -e install/install-srst2.sh
bash -e install/install-nextflow.sh
bash -e install/install-sra-toolkit.sh
bash -e install/setup-environment.sh



