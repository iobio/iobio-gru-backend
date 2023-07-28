#!/bin/bash
set -euo pipefail

vcfUrl=$1
tbiUrl=$2
refName=$3
regions=$4
contigStr=$5
refFastaFile=$6


# TODO: it doesn't seem to work unless I have a file extension on the end...
printf "$contigStr" > contigs.txt

#echo vcfUrl: $vcfUrl tbiUrl: $tbiUrl refName: $refName regions: $regions contigStr: $contigStr refFastaFile: $refFastaFile

tabixVcfArg=$vcfUrl
if [ -n "${tbiUrl}" ]; then
    tabixVcfArg="$vcfUrl##idx##$tbiUrl"
fi

tabix -h $tabixVcfArg $regions | \
    bcftools annotate -h contigs.txt - | \
    bcftools norm -m - -w 10000 -f $refFastaFile -
