#!/bin/bash
#SJG Jan2020

vcfUrl=$1
qualCutoff=$2
depthCutoff=$3 #Note: not using for now - doing this on front end
normalCountCutoff=$4
tumorCountCutoff=$5
normalAfCutoff=$6
tumorAfCutoff=$7
normalSampleIdx=$8
totalSampleNum=$9

qualPhrase="QUAL>${qualCutoff}"
normalCountPhrase="AC[${normalSampleIdx}]<=${normalCountCutoff}"
normalAfPhrase="(AF[${normalSampleIdx}]<=${normalAfCutoff}||(AC[${normalSampleIdx}]/AN)<=${normalAfCutoff})"

#format tumor query pieces
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

#format final query
queryPhrase="${qualPhrase}&&${normalCountPhrase}&&${normalAfPhrase}&&${tumorCountPhrase}&&${tumorAfPhrase}"

runDir=$PWD
tempDir=$(mktemp -d)
cd $tempDir

bcftools query -f '%LINE\n' -i $queryPhrase $vcfUrl

#echo $tempDir
rm -rf $tempDir
cd $runDir

