#!/bin/bash
#SJG updated Sept2022
set -euo pipefail

vcfUrl=$1
selectedSamples=$2
regions=$3
somaticFilterPhrase=$4
refFastaFile=$5
gffFile=$6

runDir=$PWD
tempDir=$(mktemp -d)
cd $tempDir

# Add region filter stage if we have regions
regionFilterStage=cat
if [ ! -z "$regions" ]; then
    function regionFilterFunc {
        bcftools filter -t $regions -
    }
    regionFilterStage=regionFilterFunc
fi

#Add somatic filter stage if we have a filter
somFilterStage=cat
if [ ! -z "$somaticFilterPhrase" ]; then
    function somFilterFunc {
        bcftools filter -i $somaticFilterPhrase -
    }
    somFilterStage=somFilterFunc
fi

#Do work
bcftools view -s $selectedSamples $vcfUrl | \
    $regionFilterStage | \
    bcftools norm -m - -w 10000 -f $refFastaFile - | \
    $somFilterStage | \
    bcftools csq -f $refFastaFile -g $gffFile - -Ov -l

rm -rf $tempDir
cd $runDir

