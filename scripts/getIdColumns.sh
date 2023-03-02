#!/bin/bash
set -euo pipefail

vcfUrl=$1
region=$2

bcftools view -r $region $vcfUrl
