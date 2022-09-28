#!/bin/bash
set -euo pipefail

url=$1
indexUrl=$2
samtoolsRegion=$3
bamstatsRegions=$4
dataDir=$5

tempDir=$(mktemp -d)
cd $tempDir

dataOpts=$url
if [ -n "${indexUrl}" ]; then
    dataOpts="-X $url $indexUrl"
fi

export REF_CACHE=$dataDir/md5_reference_cache/%2s/%2s/%s

samtools-1.11 view -b $dataOpts $samtoolsRegion | bamstatsAlive -u 500 -k 1 -r $bamstatsRegions

rm -rf $tempDir
