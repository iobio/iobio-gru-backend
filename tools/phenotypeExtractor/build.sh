#!/bin/bash

# TODO: Figure out how to activate from the Dockerfile
source miniconda3/bin/activate

pip3 install nltk==3.4.5

python3 -m nltk.downloader -d nltk_data wordnet

git clone https://github.com/adityaekawade/Phenotype-extractor.git

cd Phenotype-extractor

git checkout 198bf0a55816ae172aa27619f9e1e02c1c8396e8

npm install
