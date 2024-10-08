#!/bin/bash

# TODO: add 'set -u' once we fix unbound variables
#set -euo pipefail
set -eo pipefail

alignmentUrls=$1
alignmentIndices=$2
region=$3
refFastaFile=$4
useSuggestedVariants=$5
clinvarUrl=$6
genomeBuildName=$7
vepREVELFile=$8
samplesFileStr=${9}
extraArgs=${10}
vepCacheDir=${11}
vepPluginDir=${12}
decompose=${13}
contigStr=${14}
dataDir=${15}

contigFile="contig.txt"
echo -e "$contigStr" > $contigFile

samplesFile="samples.txt"
echo -e "$samplesFileStr" > $samplesFile

export REF_CACHE=$dataDir/md5_reference_cache/%2s/%2s/%s

freebayesArgs="-s $samplesFile"

# split alignments by ','
IFS=','
read -ra urls <<< "$alignmentUrls"
read -ra indices <<< "$alignmentIndices"
for i in "${!urls[@]}"; do
    url=${urls[$i]}
    indexUrl=${indices[$i]}
    alignmentFile=$(mktemp --tmpdir=./)

    if [ -n "${indexUrl}" ]; then
        samtools view -b -X $url $indexUrl $region > $alignmentFile &
    else
        samtools view -b $url $region > $alignmentFile &
    fi

    freebayesArgs="$freebayesArgs -b $alignmentFile"
done
IFS=' '

tabixCommand=''
if [ "$useSuggestedVariants" == "true" ]; then
    suggFile="sugg.vcf"
    tabix -h $clinvarUrl $region | vt view -f "INFO.CLNSIG=~'5|4'" - > $suggFile
    freebayesArgs="$freebayesArgs -@ $suggFile"
fi

decomposeStage=''
if [ "$decompose" == "true" ]; then
    function decomposeFunc {
        vt decompose -s  -
    }
    decomposeStage=decomposeFunc
fi


if [ "$genomeBuildName" == "GRCh38" ]; then
    toml="/gru_data/annotations/GRCh38/vcfanno_annotate_variants.toml"
else
    toml="/gru_data/annotations/GRCh37/vcfanno_annotate_variants.toml"
fi

custom_lua="/gru_data/annotations/vcfanno_custom.lua"

vcfanno="vcfanno --lua $custom_lua $toml /dev/stdin"


freebayesArgs="$freebayesArgs $extraArgs"

vepBaseArgs="-i STDIN --format vcf --cache --dir_cache $vepCacheDir --offline --vcf -o STDOUT --no_stats --no_escape --sift b --polyphen b --regulatory --fork 4 --merged --fasta $refFastaFile"

vepArgs="$vepBaseArgs --assembly $genomeBuildName --allele_number --hgvs --check_existing"

if [ "$vepREVELFile" ]; then
    vepArgs="$vepArgs --dir_plugins $vepPluginDir --plugin REVEL,$vepREVELFile"
fi

wait

# TODO: vt filter -d used to be "Variants called by iobio" when invoked by
# minion, but for some reason I can't get it to accept more than a single word
# directly in the bash script
freebayes -f $refFastaFile $freebayesArgs | \
    $decomposeStage | \
    bcftools norm -m - -w 10000 -f $refFastaFile - \ |
    vt filter -f 'QUAL>1' -t 'PASS' -d 'iobio' - | \
    bcftools annotate -h $contigFile | \
    vep $vepArgs | \
    $vcfanno
