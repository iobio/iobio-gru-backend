#!/bin/bash

vcfUrl=$1
tbiUrl=$2
region=$3
contigStr=$4
refFastaFile=$5
genomeBuildName=$6
gnomadMergeAnnots=$7
pathoFilterPhrase=${8}

runDir=$PWD
tempDir=$(mktemp -d)
cd $tempDir

echo -e "$contigStr" > contigs.txt

tabixVcfArg=$vcfUrl
if [ -n "${tbiUrl}" ]; then
    tabixVcfArg="$vcfUrl##idx##$tbiUrl"
fi

#optional gnomad stage
gnomadAnnotStage=cat
if [ "$gnomadMergeAnnots" ]; then
    
    if [ "$genomeBuildName" == "GRCh38" ]; then
        toml="/data/gnomad/vcfanno_gnomad_3.1_grch38.toml"
    fi
    	toml="/data/gnomad/vcfanno_gnomad_2.1_grch37.toml"

    function gnomadAnnotFunc {
        # Add the gnomAD INFO fields to the input vcf
        vcfanno $toml /dev/stdin
    }

    gnomadAnnotStage=gnomadAnnotFunc
fi

tabix -h $tabixVcfArg $region | \
    bcftools annotate -h contigs.txt - | \
    vt normalize -n -r $refFastaFile - | \
    bcftools filter -i $pathoFilterPhrase - | \
    $gnomadAnnotStage

#echo $tempDir
rm -rf $tempDir
cd $runDir
