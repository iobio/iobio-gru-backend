#!/bin/bash

vcfUrl=$1
tbiUrl=$2
region=$3
contigStr=$4
vcfSampleNamesStr=$5
refFastaFile=$6
genomeBuildName=$7
vepREVELFile=$8
vepAF=$9
isRefSeq=${10}
hgvsNotation=${11}
getRsId=${12}
#globalGetRsId=${13}

if [ "$tbiUrl" == "null" ]; then
    tbiUrl=
fi

# default to no-op
subsetStage=cat

if [ "$vcfSampleNamesStr" ]; then
    samplesFile=/tmp/iobio_annotateVariants_samplesFile_$(cat /proc/sys/kernel/random/uuid)
    echo -e "$vcfSampleNamesStr" > $samplesFile

    subsetStage="vt subset -s $samplesFile -"
fi

contigFile=/tmp/iobio_annotateVariants_contigFile_$(cat /proc/sys/kernel/random/uuid)
echo -e "$contigStr" > $contigFile


vepArgs="--assembly $genomeBuildName --format vcf --allele_number"

if [ "$vepREVELFile" ]; then
    vepArgs="$vepArgs --plugin REVEL,$vepREVELFile"
fi

if [ "$vepAF" ]; then
    vepArgs="$vepArgs --af --af_gnomad --af_esp --af_1kg --max_af"
fi

if [ "$isRefSeq" ]; then
    vepArgs="$vepArgs --refseq"
fi

if [ "$hgvsNotation" ]; then
    vepArgs="$vepArgs --hgvs"
fi

if [ "$hgvsNotation" ]; then
    vepArgs="$vepArgs --check_existing"
fi

# TODO: compare globalGetRsId to what gene is doing. It's returning a
# function. Not sure what the semantics should be.
#if [ "$hgvsNotation" ] || [ "$globalGetRsId" ] || [ "$isRefSeq" ]; then
if [ "$hgvsNotation" ] || [ "$isRefSeq" ]; then
    vepArgs="$vepArgs --fasta $refFastaFile"
fi

echo vepArgs: $vepArgs > vepArgs.txt

tabix_od -h $vcfUrl $region $tbiUrl | \
    bcftools annotate -h $contigFile - | \
    $subsetStage | \
    vt normalize -n -r $refFastaFile - | \
    vep $vepArgs 

rm $samplesFile
rm $contigFile
