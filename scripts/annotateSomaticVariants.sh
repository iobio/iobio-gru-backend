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
normalSampleIdx=${10}
totalSampleNum=${11}
genomeBuildName=${12}
vepCacheDir=${13}

echo -e "$regions" >> regions.txt

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


vepArgs="--assembly $genomeBuildName --format vcf --allele_number --dir_cache $vepCacheDir"

runDir=$PWD
tempDir=$(mktemp -d)
cd $tempDir

bcftools view -s $selectedSamples $vcfUrl | \
    bcftools filter -i $queryPhrase -t $regions - | \
    vep $vepArgs

#todo: need to do actual vep annotation call here dummy!

#echo $tempDir
rm -rf $tempDir
cd $runDir

