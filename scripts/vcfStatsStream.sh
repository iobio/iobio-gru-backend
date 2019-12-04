#!/bin/bash

url=$1
indexUrl=$2
regionsStr=$3
contigStr=$4

tempDir=$(mktemp -d)
cd $tempDir

echo -e "$contigStr" > contigs.txt

tabix_od -h $url $regionsStr $indexUrl | \
        bcftools annotate -h contigs.txt - | \
        vcfStatsAlive -u 1000 -Q 1000

cd $runDir
rm -rf $tempDir
