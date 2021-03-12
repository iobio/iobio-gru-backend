#!/bin/bash

url=$1
indexUrl=$2
regionsStr=$3
contigStr=$4
sampleNamesStr=$5

tabixVcfArg=$url
if [ -n "${indexUrl}" ]; then
    tabixVcfArg="$url##idx##$indexUrl"
fi

tempDir=$(mktemp -d)
cd $tempDir

echo -e "$contigStr" > contigs.txt

if [ "$sampleNamesStr" ]; then
        echo -e "$sampleNamesStr" > sample_names.txt
        tabix -h $tabixVcfArg $regionsStr | \
                bcftools annotate -h contigs.txt - | \
                vt subset -s sample_names.txt - | \
                vcfStatsAlive -u 1000 -Q 1000
else
        tabix -h $tabixVcfArg $regionsStr | \
                bcftools annotate -h contigs.txt - | \
                vcfStatsAlive -u 1000 -Q 1000
fi

cd $runDir
rm -rf $tempDir
