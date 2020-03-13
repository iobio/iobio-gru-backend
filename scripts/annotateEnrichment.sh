#!/bin/bash

vcfUrl=$1
tbiUrl=$2
region=$3
contigStr=$4
refFastaFile=$5
genomeBuildName=$6
filterArgs=$7
expSamples=$8
controlNamesFile=$9

echo -e "$contigStr" > contigs.txt

tabix_od -h $vcfUrl $region $tbiUrl | \
    bcftools annotate -h contigs.txt - | \
    vt subset -s $controlNamesFile - | \
    $filterCmd | \
    gtenricher expSamples

rm -rf $tempDir
cd $runDir
