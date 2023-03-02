#!/bin/bash
set -euo pipefail

vcfUrl=$1
indexUrl=$2


args=$vcfUrl
if [ -n "${indexUrl}" ]; then
    args="$vcfUrl##idx##$indexUrl"
fi

tabix -H $args
