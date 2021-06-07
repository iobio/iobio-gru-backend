#!/bin/bash

curl -LO https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh

bash Miniconda3-latest-Linux-x86_64.sh -b -p $PWD/miniconda3

rm Miniconda3-latest-Linux-x86_64.sh

./miniconda3/bin/pip3 install pyvcf
