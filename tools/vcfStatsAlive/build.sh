#!/bin/bash

mkdir -p /build
cd /build
git clone --recursive https://github.com/yiq/vcfstatsalive.git
cd vcfstatsalive/lib/htslib
make -j4
mv libhts.a ../
cd ../../
HTSLIB_HOME=. make -j4
strip vcfstatsalive
cp vcfstatsalive /vcfStatsAlive
