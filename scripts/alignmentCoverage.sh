#!/bin/bash

url=$1
index_url=$2
samtools_region=$3
max_points=$4
spanning_region=$5
coverage_regions=$6
quality_threshold=$7

#if quality value provided, filter reads by mapq
# otherwise just add binary flag
view_opt="-b"
if [ "quality_threshold" ]; then
    view_opt="-b -q $quality_threshold"
fi

samtools_od view $view_opt $url $samtools_region $index_url | \
    samtools mpileup - | \
    coverage $max_points $spanning_region $coverage_regions
