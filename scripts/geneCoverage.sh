#!/bin/bash

url=$1
index_url=$2
samtools_region=$3
region_file_str=$4

runDir=$PWD
tempDir=$(mktemp -d)
cd $tempDir

# TODO: it doesn't seem to work unless I have a file extension on the end...
printf "$region_file_str" > region_file.txt

samtools_od view -b $url $samtools_region $index_url > alignment.aln

samtools index alignment.aln

geneCoverage -b alignment.aln -r region_file.txt

cd $runDir
rm -rf $tempDir
