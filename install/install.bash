set -e

# TODO: DRY up the apt-get so it's in one call

apt-get -y update
DEBIAN_FRONTEND=noninteractive apt-get install -y -q sudo
DEBIAN_FRONTEND=noninteractive apt-get install -y -q apt-utils
DEBIAN_FRONTEND=noninteractive apt-get install -y -q make
DEBIAN_FRONTEND=noninteractive apt-get install -y -q gcc
DEBIAN_FRONTEND=noninteractive apt-get install -y -q bc

DEBIAN_FRONTEND=noninteractive apt-get install -y -q git
DEBIAN_FRONTEND=noninteractive apt-get install -y -q wget
DEBIAN_FRONTEND=noninteractive apt-get install -y -q curl

DEBIAN_FRONTEND=noninteractive apt-get install -y -q libz-dev
DEBIAN_FRONTEND=noninteractive apt-get install -y -q libbz2-dev
DEBIAN_FRONTEND=noninteractive apt-get install -y -q libncurses5-dev
DEBIAN_FRONTEND=noninteractive apt-get install -y -q libncursesw5-dev

DEBIAN_FRONTEND=noninteractive apt-get install -y -q libghc-bzlib-prof
DEBIAN_FRONTEND=noninteractive apt-get install -y -q zlib1g-dev
DEBIAN_FRONTEND=noninteractive apt-get install -y -q libcurl4-openssl-dev

DEBIAN_FRONTEND=noninteractive apt-get install -y -q libfile-copy-recursive-perl
DEBIAN_FRONTEND=noninteractive apt-get install -y -q libio-socket-ssl-perl
DEBIAN_FRONTEND=noninteractive apt-get install -y -q libio-tee-perl
DEBIAN_FRONTEND=noninteractive apt-get install -y -q libunicode-string-perl

DEBIAN_FRONTEND=noninteractive apt-get install -y -q nano
DEBIAN_FRONTEND=noninteractive apt-get install -y -q python3

cp ./install-fastp.sh ./install
bash -e ./install/install-fastp.sh

cp ./install-fastqc.sh ./install/install-fastqc.sh
bash -e install/install-fastqc.sh

#export 
PATH=$PATH:/opt/conda/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin
cp ./install-conda.sh ./install/install-conda.sh
bash -e install/install-conda.sh

cp ./install-seqtk.sh ./install/install-seqtk.sh
bash -e install/install-seqtk.sh

cp ./install-shovill.sh ./install/install-shovill.sh
bash -e install/install-shovill.sh

cp ./install-quast.sh ./install/install-quast.sh
bash -e install/install-quast.sh

cp ./install-most.sh ./install/install-most.sh
bash -e install/install-most.sh

cp ./install-kmerid.sh ./install/install-kmerid.sh
bash -e install/install-kmerid.sh

cp ./install-seqsero2.sh ./install/install-seqsero2.sh
bash -e install/install-seqsero2.sh

cp ./install-sistr.sh ./install/install-sistr.sh
bash -e install/install-sistr.sh

cp ./install-srst2.sh ./install/install-srst2.sh
bash -e install/install-srst2.sh

cp ./install-nextflow.sh ./install/install-nextflow.sh
bash -e install/install-nextflow.sh

cp ./install-sra-toolkit.sh ./install/install-sra-toolkit.sh
bash -e install/install-sra-toolkit.sh

cp ./setup-environment.sh ./install/setup-environment.sh
bash -e install/setup-environment.sh



