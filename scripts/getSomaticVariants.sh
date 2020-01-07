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

# have to translate to 1 based indexing
$normalSampleIdx = $normalSampleIdx + 1

qualPhrase="QUAL>${qualCutoff}"
depthPhrase="AN>${totalReadCutoff}"
normalCountPhrase="AC[${normalSampleIdx}]<=${normalCountCutoff}"
normalAfPhrase="AF[${normalSampleIdx}]<=${normalAfCutoff}||(AC[${normalSampleIdx}]/AN)<=${normalAfCutoff}"

tumorCountPhrase="("
tumorAfPhrase="("
for i in $totalSampleNum; do
	if ($i != $normalSampleIdx)
		j=i+1
		$tumorCountPhrase="${tumorCountPhrase}AC[${j}]>=${tumorCountCutoff}"
		$tumorAfPhrase="${tumorAfPhrase}((AF[${j}]>=${tumorAfCutoff})||((AC[${j}]/AN)<=${tumorAfCutoff}))"
		if ($i < $totalSampleNum-1)
			$tumorCountPhrase="${tumorCountPhrase}||"
			$tumorAfPhrase="$tumorAfPhrase||"
		fi
	fi
done
$tumorCountPhrase="${tumorCountPhrase})"
$tumorAfPhrase="${tumorAfPhrase})"
echo $tumorCountPhrase
echo $tumorAfPhrase

queryPhrase="${qualPhrase}&&${depthPhrase}&&${normalCountPhrase}&&${normalAfPhrase}&&${tumorCountPhrase}&&${tumorAfPhrase}"
echo $queryPhrase

runDir=$PWD
tempDir=$(mktemp -d)
cd $tempDir

#bcftools query -f '%CHROM %POS %REF %ALT %QUAL %INFO %FILTER' | awk '{ if ($5 >= int($qualCutoff) print $1 $2 $3 $4 $5 $6 $7 }'

bcftools query -f '%LINE\n' -i $queryPhrase $vcfUrl

#awk '{ if (int($6) > int($qualCutoff)) print $6 }'

#echo $tempDir
rm -rf $tempDir
cd $runDir

