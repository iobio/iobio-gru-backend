#!/bin/bash
set -euo pipefail

vcfUrl=$1
tbiUrl=$2
region=$3
contigStr=$4
vcfSampleNamesStr=$5
refFastaFile=$6
genomeBuildName=$7
vepCacheDir=$8
vepREVELFile=$9
vepPluginDir=${10}
hgvsNotation=${11}
getRsId=${12}
decompose=${13}

# default optional stages to no-op
subsetStage=cat
decomposeStage=cat

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

if [ "$genomeBuildName" == "GRCh38" ]; then
    toml="/gru_data/vcfanno/GRCh38/annotate_variants.toml"
else
    toml="/gru_data/vcfanno/GRCh37/annotate_variants.toml"
fi

custom_lua="/gru_data/vcfanno/custom.lua"

vcfanno="vcfanno --lua $custom_lua $toml /dev/stdin"

echo -e "$contigStr" > contigs.txt

vepBaseArgs="-i STDIN --format vcf --cache --dir_cache $vepCacheDir --offline --vcf -o STDOUT --no_stats --no_escape --sift b --polyphen b --regulatory --fork 4 --merged --fasta $refFastaFile"

vepArgs="$vepBaseArgs --assembly $genomeBuildName --allele_number"

if [ "$hgvsNotation" == "true" ]; then
    vepArgs="$vepArgs --hgvs"
fi

if [ "$getRsId" == "true" ]; then
    vepArgs="$vepArgs --check_existing"
fi

if [ "$vepREVELFile" ]; then
    vepArgs="$vepArgs --dir_plugins $vepPluginDir --plugin REVEL,$vepREVELFile"
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
    $vcfanno
