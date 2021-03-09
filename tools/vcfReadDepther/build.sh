#!/bin/bash

curl -LO https://raw.githubusercontent.com/tonydisera/vcfReadDepther/master/sampler.cpp
g++ -o vcfReadDepther sampler.cpp
strip vcfReadDepther

mkdir -p build/lib/vcfReadDepther/
mv vcfReadDepther build/lib/vcfReadDepther/

cp vcfReadDeptherHelper.sh build/vcfReadDepther
