#!/bin/bash

apt-get update && apt-get -y install wget clang zlib1g-dev

wget -q https://raw.githubusercontent.com/chmille4/bamReadDepther/master/bamReadDepther.cpp

clang++ -o bamReadDepther -static -O3 bamReadDepther.cpp

mv bamReadDepther build/bamReadDepther
