#!/bin/bash

curl -LO https://raw.githubusercontent.com/tonydisera/vcfReadDepther/master/sampler.cpp
g++ -o vcfReadDepther sampler.cpp
strip vcfReadDepther
