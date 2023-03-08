#!/bin/bash
set -euo pipefail

vcfUrl=$1
tbiUrl=$2

tabixVcfArg=$vcfUrl
if [ -n "${tbiUrl}" ]; then
    tabixVcfArg="$vcfUrl##idx##$tbiUrl"
fi

tabix -l $tabixVcfArg
