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
headerArg=$""
if [ ! -z "$lineNumber" ]; then
    function headFunc {
        head -n $lineNumber
    }
    headStage=headFunc
    headerArg="-H" #Don't include header
fi

bcftools view $headerArg $regionArg $vcfUrl | \
    headStage
