#!/bin/bash

url=$1
index_url=$2
samtools_region=$3
quality_cutoff=$4

samtools_od view -q $quality_cutoff -c $url $samtools_region $index_url
