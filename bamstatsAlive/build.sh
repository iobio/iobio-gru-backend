#!/bin/bash

BAMTOOLS_VERSION=2.5.1

apt-get update && apt-get -y install wget curl git cmake build-essential zlib1g-dev

wget -q https://github.com/pezmaster31/bamtools/archive/v${BAMTOOLS_VERSION}.tar.gz
tar xzf v${BAMTOOLS_VERSION}.tar.gz
pushd bamtools-${BAMTOOLS_VERSION}
mkdir build
pushd build
cmake -DEnableNodeJS=true ..
make -j4
make install
popd
popd

git clone https://github.com/yiq/bamstatsAlive

pushd bamstatsAlive
# TODO: I'm sure there's a way to override just CFLAGS=-static without
# resorting to a git patch, but I couldn't get it to work in 5 minutes so
# gave up.
git apply ../build/static_link.patch
export BAMTOOLS=../bamtools-${BAMTOOLS_VERSION}
make
popd

cp bamstatsAlive/bamstatsAlive build/
