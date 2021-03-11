#!/bin/bash

url=$1
index_url=$2
samtools_region=$3
max_points=$4
spanning_region=$5
coverage_regions=$6
quality_threshold=$7
data_dir=$8

tempDir=$(mktemp -d)
cd $tempDir

#if quality value provided, filter reads by mapq
# otherwise just add binary flag
view_opt="-b"
if [ "quality_threshold" ]; then
    view_opt="-b -q $quality_threshold"
fi

data_opts=$url
if [ -n "${index_url}" ]; then
    data_opts="-X $url $index_url"
fi

export REF_CACHE=$data_dir/md5_reference_cache/%2s/%2s/%s

samtools view $view_opt $data_opts $samtools_region | \
    samtools mpileup - | \
    summarize_coverage $max_points $spanning_region $coverage_regions

rm -rf $tempDir
