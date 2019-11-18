#!/bin/bash

mkdir -p static/GRCh37/2018
wget -P static/GRCh37/2018 ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/archive_2.0/2018/clinvar_20181202.vcf.gz
wget -P static/GRCh37/2018 ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/archive_2.0/2018/clinvar_20181202.vcf.gz.tbi
mkdir -p static/GRCh37/2019
wget -P static/GRCh37/2019 ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/archive_2.0/2019/clinvar_20191021.vcf.gz
wget -P static/GRCh37/2019 ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/archive_2.0/2019/clinvar_20191021.vcf.gz.tbi

mkdir -p static/GRCh38/2018
wget -P static/GRCh38/2018 ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/archive_2.0/2018/clinvar_20181202.vcf.gz
wget -P static/GRCh38/2018 ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/archive_2.0/2018/clinvar_20181202.vcf.gz.tbi
mkdir -p static/GRCh38/2019
wget -P static/GRCh38/2019 ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/archive_2.0/2019/clinvar_20191021.vcf.gz
wget -P static/GRCh38/2019 ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/archive_2.0/2019/clinvar_20191021.vcf.gz.tbi
