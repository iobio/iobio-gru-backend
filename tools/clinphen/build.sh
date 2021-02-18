#!/bin/bash

curl -LO https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh

bash Miniconda3-latest-Linux-x86_64.sh -b -p $PWD/miniconda3

./miniconda3/bin/pip3 install clinphen

./miniconda3/bin/python3 -m nltk.downloader -d nltk_data wordnet
