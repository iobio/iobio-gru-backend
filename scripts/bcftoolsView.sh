#!/bin/bash
# SJG updated Aug2023
# General purpose record retrieval script
# Can add further functionality as needed in future
set -euo pipefail

vcfUrl=$1
region=$2
lineNumber=$3

regionArg=$region
if [ ! -z "$region" ]; then
    regionArg="-r $region"
fi

headStage=cat
if [ ! -z "$lineNumber" ]; then
    function headFunc {
        head -n $lineNumber
    }
    headStage=headFunc
fi

bcftools view $regionArg $vcfUrl | \
    headStage
