#!/bin/bash

curl -LO https://raw.githubusercontent.com/tonydisera/vcfReadDepther/master/sampler.cpp
g++ -o vcfReadDepther sampler.cpp
strip vcfReadDepther

mkdir -p lib/vcfReadDepther/
mv vcfReadDepther lib/vcfReadDepther/

cp vcfReadDeptherHelper.sh vcfReadDepther
