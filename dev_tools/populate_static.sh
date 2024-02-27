#!/bin/bash

set -xeuo pipefail

mkdir -p static/clinvar/GRCh37
curl --fail -o static/clinvar/GRCh37/clinvar.vcf.gz https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/weekly/clinvar_20230930.vcf.gz
curl --fail -o static/clinvar/GRCh37/clinvar.vcf.gz.tbi https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/weekly/clinvar_20230930.vcf.gz.tbi
mkdir -p static/clinvar/GRCh38
curl --fail -o static/clinvar/GRCh38/clinvar.vcf.gz https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/weekly/clinvar_20230930.vcf.gz
curl --fail -o static/clinvar/GRCh38/clinvar.vcf.gz.tbi https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/weekly/clinvar_20230930.vcf.gz.tbi

mkdir -p static/ensembl/GRCh37
curl --fail -o static/ensembl/GRCh37/geneSet37.gff3.gz https://ftp.ensembl.org/pub/grch37/current/gff3/homo_sapiens/Homo_sapiens.GRCh37.87.chr.gff3.gz
mkdir -p static/ensembl/GRCh38
curl --fail https://ftp.ensembl.org/pub/release-109/gff3/homo_sapiens/Homo_sapiens.GRCh38.109.chr.gff3.gz | gzip -d - | perl -pe 's/^([0-9]+|[X]|[Y]|[M])/chr$1/' - | gzip - >> static/ensembl/GRCh38/geneSet38.gff3.gz 
