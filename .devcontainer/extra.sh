#!/bin/bash

set -uxe
set -o pipefail

VER_SAMTOOLS=1.21
VER_BCFTOOLS=1.21

CPU=`grep -c ^processor /proc/cpuinfo`
if [ $? -eq 0 ]; then
  if [ "$CPU" -gt "6" ]; then
    CPU=6
  fi
else
  CPU=1
fi

# install htslib & samtools if not already present
set +e
hash samtools >& /dev/null || (
    set -e
    sudo apt-get -y update
    sudo apt-get -y install autoconf automake make gcc perl zlib1g-dev libbz2-dev liblzma-dev libcurl4-gnutls-dev libssl-dev libncurses5-dev libdeflate-dev
    cd /tmp
    curl -sSL https://github.com/samtools/samtools/releases/download/${VER_SAMTOOLS}/samtools-${VER_SAMTOOLS}.tar.bz2 | tar -jx
    cd samtools-*
    ./configure
    make -j${CPU}
    sudo make install-htslib
    sudo make install
    cd /tmp
    rm -rf samtools-*
)

# install bcftools if not already present
set +e
hash bcftools >& /dev/null || (
    set -e
    cd /tmp
    curl -sSL https://github.com/samtools/bcftools/releases/download/${VER_BCFTOOLS}/bcftools-${VER_BCFTOOLS}.tar.bz2 | tar -jx
    cd bcftools-*
    ./configure
    make -j${CPU}
    sudo make install
    cd /tmp
    rm -rf bcftools-*
)

pip install --upgrade pip

## global install of pre-commit
pip install pre-commit
pre-commit install --install-hooks

(cd /tmp && sudo rm -rf gitflow && git clone https://github.com/datasift/gitflow && cd gitflow && sudo ./install.sh && sudo git hf upgrade)
sudo rm -rf /tmp/gitflow


# set +e
# # this fails if unstaged changes, however they should only exist when resuming
# git hf init -f
# hash nextflow >& /dev/null || (
#     set -e
#     cd /tmp
#     curl -s https://get.nextflow.io | bash
#     mkdir -p $HOME/.local/bin
#     mv nextflow $HOME/.local/bin/
# )

# pip install nf-core

exit 0
