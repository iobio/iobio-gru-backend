#!/bin/bash

BAMTOOLS_VERSION=2.5.1

curl -LO https://github.com/pezmaster31/bamtools/archive/v${BAMTOOLS_VERSION}.tar.gz
tar xf v${BAMTOOLS_VERSION}.tar.gz
pushd bamtools-${BAMTOOLS_VERSION}
mkdir build
pushd build
cmake -DEnableNodeJS=true ..
make -j4
make install
popd
popd

git clone --depth 1 https://github.com/yiq/bamstatsAlive

pushd bamstatsAlive
export BAMTOOLS=../bamtools-${BAMTOOLS_VERSION}
make
popd
