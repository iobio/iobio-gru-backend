#!/bin/bash

curl -LO https://raw.githubusercontent.com/chmille4/bamReadDepther/master/bamReadDepther.cpp

clang++ -o bamReadDepther -static -O3 bamReadDepther.cpp
