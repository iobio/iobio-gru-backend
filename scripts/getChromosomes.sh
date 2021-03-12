#!/bin/bash

vcfUrl=$1
tbiUrl=$2

tempDir=$(mktemp -d)
cd $tempDir

tabixVcfArg=$vcfUrl
if [ -n "${tbiUrl}" ]; then
    tabixVcfArg="$vcfUrl##idx##$tbiUrl"
fi

tabix -l $tabixVcfArg

rm -rf $tempDir
