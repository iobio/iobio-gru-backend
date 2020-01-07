#!/bin/bash

vcfUrl=$1
qualCutoff=$2
totalReadCutoff=$3
normalCountCutoff=$4
tumorCountCutoff=$5
normalAfCutoff=$6
tumorAfCutoff=$7
normalSampleIdx=$8
totalSampleNum=$9

qualPhrase="QUAL>${qualCutoff}"

# TODO: this needs to be AO cumulative of all genotypes
# if we have an alt in sample we're looking at, must pass quality threshold (
depthPhrase="AN>${totalReadCutoff}"

normalCountPhrase="AC[${normalSampleIdx}]<=${normalCountCutoff}"
normalAfPhrase="(AF[${normalSampleIdx}]<=${normalAfCutoff}||(AC[${normalSampleIdx}]/AN)<=${normalAfCutoff})"

tumorCountPhrase="("
tumorAfPhrase="("
for ((i=0; i<totalSampleNum; i++)); do
	if ((i != normalSampleIdx)); then
		tumorCountPhrase="${tumorCountPhrase}AC[${i}]>=${tumorCountCutoff}"
		tumorAfPhrase="${tumorAfPhrase}((AF[${i}]>=${tumorAfCutoff})||((AC[${i}]/AN)<=${tumorAfCutoff}))"
		if ((i < totalSampleNum-1)); then
			tumorCountPhrase="${tumorCountPhrase}||"
			tumorAfPhrase="$tumorAfPhrase||"
		fi
	fi
done
tumorCountPhrase="${tumorCountPhrase})"
tumorAfPhrase="${tumorAfPhrase})"

queryPhrase="${qualPhrase}&&${normalCountPhrase}&&${normalAfPhrase}&&${tumorCountPhrase}&&${tumorAfPhrase}"
echo $queryPhrase

runDir=$PWD
tempDir=$(mktemp -d)
cd $tempDir

bcftools query -f '%LINE\n' -i $queryPhrase $vcfUrl | head -10

#echo $tempDir
rm -rf $tempDir
cd $runDir

