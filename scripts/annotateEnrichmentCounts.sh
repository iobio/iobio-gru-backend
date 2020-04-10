#!/bin/bash

vcfUrl=$1
tbiUrl=$2
regionStr=$3
contigStr=$4
refFastaFile=$5
filterArgs=$6
experStr=$7
controlStr=$8

echo -e "$contigStr" > contigs.txt
echo -e "$controlStr" > controlNames.txt

filterCmd=cat

if [ -z "$filterArgs" ]; then
	echo "Filtering with args"
	echo "$filterArgs"

	function filterFunc {
	    vt filter $filterArgs
	}
	filterCmd=filterFunc
fi

tabix_od -h $vcfUrl $regionStr $tbiUrl | \
    bcftools annotate -h contigs.txt - | \
    vt subset -s controlNames.txt - | \
    $filterCmd | \
    gtEnricher $experStr

rm -rf $tempDir
cd $runDir
