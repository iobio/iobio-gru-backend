#!/bin/bash

# TODO: Figure out how to activate from the Dockerfile
source miniconda3/bin/activate

pip3 install nltk==3.4.5

python3 -m nltk.downloader -d nltk_data wordnet

git clone https://github.com/adityaekawade/Phenotype-extractor.git

cd Phenotype-extractor

git checkout 9cb4c0c7de850f7c839d1e1f229a02e9be04e3d3

npm install
