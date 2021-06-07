#!/bin/bash

HTSLIB_VERSION=${1}
SAMTOOLS_VERSION=${HTSLIB_VERSION}
BCFTOOLS_VERSION=${HTSLIB_VERSION}

curl -LO https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2
tar xvf samtools-${SAMTOOLS_VERSION}.tar.bz2
cd samtools-${SAMTOOLS_VERSION}

# build htslib (including tabix and bgzip)
cd htslib-${HTSLIB_VERSION}
./configure --without-curses
make -j4
strip tabix
strip bgzip

# build samtools
cd ..
./configure --without-curses
make -j4
strip samtools

cd /
curl -LO https://github.com/samtools/bcftools/releases/download/${BCFTOOLS_VERSION}/bcftools-${BCFTOOLS_VERSION}.tar.bz2
tar xvf bcftools-${BCFTOOLS_VERSION}.tar.bz2
cd bcftools-${BCFTOOLS_VERSION}
./configure --without-curses
make -j4
strip bcftools
