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

#cp ./install-fastp.sh ./install/install-fastp.sh
bash -e ./install-fastp.sh

cp ./install-fastqc.sh ./install/install-fastqc.sh
bash -e install/install-fastqc.sh
