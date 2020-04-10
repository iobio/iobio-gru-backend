#!/bin/bash
#SJG Jan2020

vcfUrl=$1
selectedSamples=$2
regions=$3
qualCutoff=$4
depthCutoff=$5 #Note: not using for now - doing this on front end
normalCountCutoff=$6
tumorCountCutoff=$7
normalAfCutoff=$8
tumorAfCutoff=$9
normalSampleIdx=$10
totalSampleNum=$11

echo -e "$regions" > regions.txt

qualPhrase="QUAL>${qualCutoff}"
normalCountPhrase="AC[${normalSampleIdx}]<=${normalCountCutoff}"
normalAfPhrase="(AF[${normalSampleIdx}]<=${normalAfCutoff}||(AC[${normalSampleIdx}]/AN)<=${normalAfCutoff})"

#filter for selected samples if provided
#todo: this phrase needs to be tested
#todo: regions param also needs to be tested
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

bcftools query -f '%LINE\n' $samplePhrase -i $queryPhrase $vcfUrl -R regions.txt

#echo $tempDir
rm -rf $tempDir
cd $runDir

