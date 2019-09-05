#!/bin/bash

#set -x

alignmentUrls=$1
alignmentIndices=$2
region=$3
refFastaFile=$4
useSuggestedVariants=$5
clinvarUrl=$6
genomeBuildName=$7
vepREVELFile=$8
vepAF=$9
isRefSeq=${10}
samplesFileStr=${11}
extraArgs=${12}

contigFile=$(mktemp)
echo -e "$contigStr" > $contigFile

samplesFile=$(mktemp)
echo -e "$samplesFileStr" > $samplesFile


# TODO: using samtools and tabix instead of samtools_od and tabix_od would
# likely be faster

freebayesArgs="-s $samplesFile"

# split alignments by ','
IFS=','
read -ra urls <<< "$alignmentUrls"
read -ra indices <<< "$alignmentIndices"
for i in "${!urls[@]}"; do
    url=${urls[$i]}
    indexUrl=${indices[$i]}
    alignmentFile=$(mktemp)
    samtools_od view -b $url $region $indexUrl > $alignmentFile &
    freebayesArgs="$freebayesArgs -b $alignmentFile"
done
IFS=' '

tabixCommand=''
if [ "$useSuggestedVariants" == "true" ]; then
    suggFile=$(mktemp)
    tabix_od -h $clinvarUrl $region | vt view -f "INFO.CLNSIG=~'5|4'" - > $suggFile
    freebayesArgs="$freebayesArgs -@ $suggFile"
fi

freebayesArgs="$freebayesArgs $extraArgs"

vepArgs="--assembly $genomeBuildName --format vcf --allele_number --hgvs --check_existing --fasta $refFastaFile"

if [ "$vepREVELFile" ]; then
    vepArgs="$vepArgs --dir_plugins ./data/vep-cache/Plugins --plugin REVEL,$vepREVELFile"
fi

if [ "$vepAF" == "true" ]; then
    vepArgs="$vepArgs --af --af_gnomad --af_esp --af_1kg --max_af"
fi

if [ "$isRefSeq" == "true" ]; then
    vepArgs="$vepArgs --refseq"
fi


wait

# TODO: vt filter -d used to be "Variants called by iobio" when invoked by
# minion, but for some reason I can't get it to accept more than a single word
# directly in the bash script
freebayes -f $refFastaFile $freebayesArgs | \
    vt normalize -r $refFastaFile - | \
    vt filter -f 'QUAL>1' -t 'PASS' -d 'iobio' - | \
    bcftools annotate -h $contigFile | \
    vep $vepArgs


rm $suggFile
