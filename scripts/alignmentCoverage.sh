#!/bin/bash

url=$1
index_url=$2
samtools_region=$3
max_points=$4
spanning_region=$5
coverage_regions=$6

samtools_od view -b $url $samtools_region $index_url | \
    samtools mpileup - | \
    coverage $max_points $spanning_region $coverage_regions
