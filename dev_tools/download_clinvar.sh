#!/bin/bash

set -xeuo pipefail

clinvar_version=20241103

mkdir -p static/clinvar/GRCh37
curl --fail -o static/clinvar/GRCh37/clinvar.vcf.gz https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/weekly/clinvar_${clinvar_version}.vcf.gz
curl --fail -o static/clinvar/GRCh37/clinvar.vcf.gz.tbi https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/weekly/clinvar_${clinvar_version}.vcf.gz.tbi
mkdir -p static/clinvar/GRCh38
curl --fail -o static/clinvar/GRCh38/clinvar.vcf.gz https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/weekly/clinvar_${clinvar_version}.vcf.gz
curl --fail -o static/clinvar/GRCh38/clinvar.vcf.gz.tbi https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/weekly/clinvar_${clinvar_version}.vcf.gz.tbi
