#!/bin/bash
set -euo pipefail

url=$1
indexUrl=$2
regionsStr=$3

tabixBedArg=$url
if [ -n "${indexUrl}" ]; then
    tabixBedArg="$url##idx##$indexUrl"
fi

runDir=$PWD
tempDir=$(mktemp -d)
cd $tempDir


tabix $tabixBedArg $regionsStr

cd $runDir
rm -rf $tempDir
