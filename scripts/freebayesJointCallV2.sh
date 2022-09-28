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
samplesFileStr=${10}
extraArgs=${11}
vepCacheDir=${12}
vepPluginDir=${13}
decompose=${14}
contigStr=${15}
dataDir=${16}

runDir=$PWD
tempDir=$(mktemp -d)
cd $tempDir

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
    alignmentFile=$(mktemp)

    if [ -n "${indexUrl}" ]; then
        samtools-1.11 view -b -X $url $indexUrl $region > $alignmentFile &
    else
        samtools-1.11 view -b $url $region > $alignmentFile &
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
        toml="$dataDir/gnomad/vcfanno_gnomad_3.1_grch38.toml"
else
    toml="$dataDir/gnomad/vcfanno_gnomad_2.1_grch37.toml"
fi

function gnomadAnnotFunc {
    # Add the gnomAD INFO fields to the input vcf
    vcfanno $toml /dev/stdin
}

gnomadAnnotStage=gnomadAnnotFunc


freebayesArgs="$freebayesArgs $extraArgs"

vepBaseArgs="-i STDIN --format vcf --cache --dir_cache $vepCacheDir --offline --vcf -o STDOUT --no_stats --no_escape --sift b --polyphen b --regulatory --fork 4 --merged --fasta $refFastaFile"

vepArgs="$vepBaseArgs --assembly $genomeBuildName --allele_number --hgvs --check_existing"

if [ "$vepREVELFile" ]; then
    vepArgs="$vepArgs --dir_plugins $vepPluginDir --plugin REVEL,$vepREVELFile"
fi

if [ "$vepAF" == "true" ]; then
    vepArgs="$vepArgs --af --af_gnomad --af_esp --af_1kg --max_af"
fi

wait

# TODO: vt filter -d used to be "Variants called by iobio" when invoked by
# minion, but for some reason I can't get it to accept more than a single word
# directly in the bash script
freebayes -f $refFastaFile $freebayesArgs | \
    $decomposeStage | \
    vt normalize -r $refFastaFile - | \
    vt filter -f 'QUAL>1' -t 'PASS' -d 'iobio' - | \
    bcftools annotate -h $contigFile | \
    vep $vepArgs | \
    $gnomadAnnotStage


rm -rf $tempDir
cd $runDir
