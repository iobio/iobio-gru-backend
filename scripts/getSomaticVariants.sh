#!/bin/bash

vcfUrl=$1
totalReadCutoff=$2

#$qualCutoff
#$normalCountCutoff
#$tumorCountCutoff
#$normalAfCutoff
#$tumorAfCutoff
#$normalSampleIdx


runDir=$PWD
tempDir=$(mktemp -d)
cd $tempDir

#bcftools query -f '%CHROM %POS %REF %ALT %QUAL %INFO %FILTER' | awk '{ if (int($5) >= int($qualCutoff) print $1 $2 $3 $4 $5 $6 $7 }'

echo Calling bcftools query...
echo $vcfUrl
bcftools query -f'%QUAL\n' -i'QUAL>0.01' $vcfUrl
echo Done calling bcftools query


#awk '{ if (int($6) > int($qualCutoff)) print $6 }'


#echo $tempDir
rm -rf $tempDir
cd $runDir

