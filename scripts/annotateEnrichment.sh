#!/bin/bash

vcfUrl=$1
tbiUrl=$2
regionStr=$3
contigStr=$4
refFastaFile=$5
genomeBuildName=$6
filterArgs=$7
experIdString=$8
controlIdString=$9

echo -e "$contigStr" > contigs.txt

# todo: put in conditional arg for filtercmd here
filterCmd=cat

if [ "$filterArgs" ]; then
	echo "Filtering with args"
	echo "$filterArgs"

	function filterFunc {
	    vt filter $filterArgs
	}
	filterCmd=filterFunc
fi

tabix_od -h $vcfUrl $regionStr $tbiUrl | \
    bcftools annotate -h contigs.txt - | \
    bcftools view -s controlIdString | \
    # todo: get rid of this line if above works - vt subset -s $controlNamesFile - | \
    filterCmd | \
    gtEnricher $experIdString

rm -rf $tempDir
cd $runDir
