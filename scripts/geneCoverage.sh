#!/bin/bash

url=$1
index_url=$2
samtools_region=$3
region_file_str=$4
data_dir=$5

runDir=$PWD
tempDir=$(mktemp -d)
cd $tempDir

# TODO: it doesn't seem to work unless I have a file extension on the end...
printf "$region_file_str" > region_file.txt

data_opts=$url
if [ -n "${index_url}" ]; then
    data_opts="-X $url $index_url"
fi

export REF_CACHE=$data_dir/md5_reference_cache/%2s/%2s/%s

samtools-1.11 view -b $data_opts $samtools_region > alignment.aln

samtools-1.11 index alignment.aln

geneCoverage -b alignment.aln -r region_file.txt

cd $runDir
rm -rf $tempDir
