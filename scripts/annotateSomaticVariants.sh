#!/bin/bash
#SJG updated Jun2020

vcfUrl=$1
selectedSamples=$2
regions=$3
somaticFilterPhrase=$4
genomeBuildName=$5
vepCacheDir=$6
refFastaFile=$7

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

#Compose args
vepBaseArgs="-i STDIN --format vcf --cache --dir_cache $vepCacheDir --offline --vcf -o STDOUT --no_stats --no_escape --sift b --polyphen b --regulatory --fork 4 --merged --fasta $refFastaFile"
vepArgs="$vepBaseArgs --assembly $genomeBuildName --allele_number"
echo -e "$regions" >> regions.txt

#Do work
bcftools view -s $selectedSamples $vcfUrl | \
    $regionFilterStage | \
    bcftools norm -m - -w 10000 -f $refFastaFile - | \
    $somFilterStage | \
    vep $vepArgs

rm -rf $tempDir
cd $runDir

