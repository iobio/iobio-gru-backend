#!/bin/bash

vcfUrl=$1
tbiUrl=$2
region=$3
contigStr=$4
vcfSampleNamesStr=$5
refFastaFile=$6
genomeBuildName=$7
vepCacheDir=$8
vepREVELFile=$9
vepAF=${10}
vepPluginDir=${11}
hgvsNotation=${12}
getRsId=${13}
gnomadMergeAnnots=${14}
decompose=${15}

# default optional stages to no-op
subsetStage=cat
gnomadAnnotStage=cat
decomposeStage=cat

runDir=$PWD
tempDir=$(mktemp -d)
cd $tempDir

if [ "$vcfSampleNamesStr" ]; then
    echo -e "$vcfSampleNamesStr" > samples.txt

    function subsetFunc {
        vt subset -s samples.txt -
    }

    subsetStage=subsetFunc
fi


if [ "$decompose" == "true" ]; then
    function decomposeFunc {
        vt decompose -s  -
    }

    decomposeStage=decomposeFunc
fi

if [ "$gnomadMergeAnnots" ]; then
    
    if [ "$genomeBuildName" == "GRCh38" ]; then
        toml="/data/gnomad/vcfanno_gnomad_3.1_grch38.toml"
    else
    	toml="/data/gnomad/vcfanno_gnomad_2.1_grch37.toml"
    fi

    function gnomadAnnotFunc {
        # Add the gnomAD INFO fields to the input vcf
        vcfanno $toml /dev/stdin
    }

    gnomadAnnotStage=gnomadAnnotFunc
fi


echo -e "$contigStr" > contigs.txt

vepBaseArgs="-i STDIN --format vcf --cache --dir_cache $vepCacheDir --offline --vcf -o STDOUT --no_stats --no_escape --sift b --polyphen b --regulatory --fork 4 --merged --fasta $refFastaFile"

vepArgs="$vepBaseArgs --assembly $genomeBuildName --allele_number"

if [ "$vepREVELFile" ]; then
    vepArgs="$vepArgs --dir_plugins $vepPluginDir --plugin REVEL,$vepREVELFile"
fi

if [ "$vepAF" == "true" ]; then
    vepArgs="$vepArgs --af --af_gnomad --af_esp --af_1kg --max_af"
fi

if [ "$hgvsNotation" == "true" ]; then
    vepArgs="$vepArgs --hgvs"
fi

if [ "$getRsId" == "true" ]; then
    vepArgs="$vepArgs --check_existing"
fi

tabixVcfArg=$vcfUrl
if [ -n "${tbiUrl}" ]; then
    tabixVcfArg="$vcfUrl##idx##$tbiUrl"
fi

tabix -h $tabixVcfArg $region | \
    bcftools annotate -h contigs.txt - | \
    $subsetStage | \
    $decomposeStage | \
    vt normalize -n -r $refFastaFile - | \
    vep $vepArgs | \
    $gnomadAnnotStage

#echo $tempDir
rm -rf $tempDir
cd $runDir
