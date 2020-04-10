#!/bin/bash
#SJG Jan2020

vcfUrl=$1
selectedSamples=$2
qualCutoff=$3
depthCutoff=$4 #Note: not using for now - doing this on front end
normalCountCutoff=$5
tumorCountCutoff=$6
normalAfCutoff=$7
tumorAfCutoff=$8
normalSampleIdx=$9
totalSampleNum=$10

qualPhrase="QUAL>${qualCutoff}"
normalCountPhrase="AC[${normalSampleIdx}]<=${normalCountCutoff}"
normalAfPhrase="(AF[${normalSampleIdx}]<=${normalAfCutoff}||(AC[${normalSampleIdx}]/AN)<=${normalAfCutoff})"

#filter for selected samples if provided
#todo: this phrase needs to be tested
samplePhrase=""
if (! -n "$selectedSamples"); then
   samplePhrase="-s $selectedSamples"
fi

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

bcftools query -f '%LINE\n' $samplePhrase -i $queryPhrase $vcfUrl

#echo $tempDir
rm -rf $tempDir
cd $runDir

