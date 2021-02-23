#!/bin/bash

HTSLIB_VERSION=1.11

# get tool & YQ htslib wrapper
mkdir -p /build
cd /build
git clone --recursive https://github.com/yiq/gtenricher.git
cd gtenricher/contrib 

# get htslib & build
curl -LO https://github.com/samtools/htslib/releases/download/${HTSLIB_VERSION}/htslib-${HTSLIB_VERSION}.tar.bz2
tar -xf htslib-${HTSLIB_VERSION}.tar.bz2 
mv ./htslib-${HTSLIB_VERSION} ./htslib
cd htslib
./configure --prefix=/usr/local
make
make install

# build
cd ../..
./configure --with-htslib=/build/gtenricher/contrib/htslib
make
make install

# put tool in path & slim down
strip /usr/local/bin/gtenricher
# cp /usr/local/bin/gtenricher /usr/bin/

