#!/bin/bash
set -euo pipefail

url=$1
indexUrl=$2
samtoolsRegion=$3
dataDir=$4

dataOpts=$url
if [ -n "${indexUrl}" ]; then
    dataOpts="-X $url $indexUrl"
fi

export REF_CACHE=$dataDir/md5_reference_cache/%2s/%2s/%s

samtools depth -a $dataOpts -r $samtoolsRegion
