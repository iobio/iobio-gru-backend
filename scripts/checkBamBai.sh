#!/bin/bash

url=$1
index_url=$2
samtools_region=$3

url_valid=$(curl -s $url | head -n 1 | tr '\0' '\n')
index_valid=$(curl -s $index_url | head -n 1 | tr '\0' '\n')
if [ -z "$url_valid" ] || [ -z "$index_valid" ]; 
then
	exit 1;
fi

samtools_od view -H $url $samtools_region $index_url | head -n 1
