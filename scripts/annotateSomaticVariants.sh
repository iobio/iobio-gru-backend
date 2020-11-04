#!/bin/bash
#SJG updated Jun2020

vcfUrl=$1
selectedSamples=$2
regions=$3
somaticFilterPhrase=$4
genomeBuildName=$5
vepCacheDir=$6
fastaPath=$7

runDir=$PWD
tempDir=$(mktemp -d)
cd $tempDir

#Add somatic filter stage if we have a filter
somFilterStage=cat
if [ ! -z "$somaticFilterPhrase" ]; then
    function somFilterFunc {
        bcftools filter -i $somaticFilterPhrase -
    }
    somFilterStage=somFilterFunc
fi

#Compose args
vepArgs="--assembly $genomeBuildName --format vcf --allele_number --dir_cache $vepCacheDir"
echo -e "$regions" >> regions.txt

#Do work
bcftools view -s $selectedSamples $vcfUrl | \
    bcftools filter -t $regions - | \
    bcftools norm -m - -w 10000 -f $fastaPath - | \
    $somFilterStage | \
    vep $vepArgs

rm -rf $tempDir
cd $runDir

