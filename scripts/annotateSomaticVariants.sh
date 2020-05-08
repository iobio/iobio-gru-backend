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
normalSampleIdxs=${10} #todo: these must be line-delim single string
tumorSampleIdxs=${11} #todo: update point of call here
totalSampleNum=${11}
genomeBuildName=${12}
vepCacheDir=${13}

echo -e "$regions" >> regions.txt
echo -e "$normalSampleIdxs" >> normals.txt
echo -e "$tumorSampleIdxs" >> tumors.txt

qualPhrase="QUAL>${qualCutoff}"
depthPhrase="DP>${depthCutoff}"

cmpdNormalPhrase="("
for idx in normals.txt do
	currPhrase=""
	normalCountPhrase="FORMAT/AO[${idx}:0]<=${normalCountCutoff}"
	normalFreqPhrase="FORMAT/AO[${idx}:0]/(FORMAT/AO[${idx}:0]+FORMAT/RO[${idx}:0])<=${normalAfCutoff}"
	if (count > 0) then
		currPhrase="|("
	else
		currPhrase="("
	fi
	currPhrase="${currPhrase}${normalCountPhrase}|${normalFreqPhrase})"
	cmpdNormalPhrase="${cmpdNormalPhrase}${currPhrase}"
done
cmpdNormalPhrase="${cmpdNormalPhrase})"

cmpdTumorPhrase="("
for idx in tumors.txt do
        currPhrase=""
        tumorCountPhrase="FORMAT/AO[${idx}:0]>=${tumorCountCutoff}"
        tumorFreqPhrase="FORMAT/AO[${idx}:0]/(FORMAT/AO[${idx}:0]+FORMAT/RO[${idx}:0])>=${tumorAfCutoff}"
        if (count > 0) then
                currPhrase="|("
        else
                currPhrase="("
        fi
        currPhrase="${currPhrase}${tumorCountPhrase}|${tumorFreqPhrase})"
        cmpdTumorPhrase="${cmpdtumorPhrase}${currPhrase}"
done
cmpdTumorPhrase="${cmpdTumorPhrase})"


<<'OLD'
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

OLD

#format final query
queryPhrase="${qualPhrase}&${depthPhrase}&${cmpdNormalPhrase}&${cmpdTumorPhrase}"

#debugging
echo "${queryPhrase}"

vepArgs="--assembly $genomeBuildName --format vcf --allele_number --dir_cache $vepCacheDir"

runDir=$PWD
tempDir=$(mktemp -d)
cd $tempDir

bcftools view -s $selectedSamples $vcfUrl | \
    #bcftools norm -N -m -
    bcftools filter -i $queryPhrase -t $regions -
    #vep $vepArgs

rm -rf $tempDir
cd $runDir

