#!/bin/bash
set -euo pipefail

vcfUrl=$1
indexUrl=$2

tempDir=$(mktemp -d)
cd $tempDir

args=$vcfUrl
if [ -n "${indexUrl}" ]; then
    args="$vcfUrl##idx##$indexUrl"
fi

tabix -H $args

rm -rf $tempDir
