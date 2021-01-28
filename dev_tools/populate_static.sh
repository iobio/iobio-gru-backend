#!/bin/bash

mkdir -p static/clinvar/GRCh37
wget -P static/clinvar/GRCh37 https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar.vcf.gz
wget -P static/clinvar/GRCh37 https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar.vcf.gz.tbi
mkdir -p static/clinvar/GRCh38
wget -P static/clinvar/GRCh38 https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz
wget -P static/clinvar/GRCh38 https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz.tbi
