#!/bin/bash
#SJG Jan2020

vcfUrl=$1
selectedSamples=$2
regions=$3
somaticFilterPhrase=$4
genomeBuildName=$5
vepCacheDir=$6

#Make temp dir
runDir=$PWD
tempDir=$(mktemp -d)
cd $tempDir

#Get correct fasta file
fastaPath="/home/ubuntu/data/references/GRCh37/human_g1k_v37_decoy_phix.fasta"
if [ "${genomeBuildName}" = "GRCh38" ]; then
	fastaPath="/home/ubuntu/data/references/GRCh38/human_g1k_v38_decoy_phix.fasta"
fi

#Compose args
vepArgs="--assembly $genomeBuildName --format vcf --allele_number --dir_cache $vepCacheDir"
echo -e "$regions" >> regions.txt

#Do work
bcftools view -s $selectedSamples $vcfUrl | \
    bcftools filter -t $regions - | \
    bcftools norm -m - -w 10000 -f $fastaPath - | \
    bcftools filter -i $somaticFilterPhrase - | \
    vep $vepArgs

rm -rf $tempDir
cd $runDir

