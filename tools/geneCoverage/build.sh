#!/bin/bash

BAMTOOLS_VERSION=2.5.1

curl -LO https://github.com/pezmaster31/bamtools/archive/v${BAMTOOLS_VERSION}.tar.gz
tar xf v${BAMTOOLS_VERSION}.tar.gz

cd bamtools-${BAMTOOLS_VERSION}
mkdir build
cd build
cmake -DEnableNodeJS=true ..
make -j4
make install

cd /
git clone https://github.com/AlistairNWard/gene_coverage.git
cd gene_coverage/src
g++ -I/bamtools-${BAMTOOLS_VERSION}/src -o /genecoverage *.cpp -lbamtools -lz -lm
strip /genecoverage
