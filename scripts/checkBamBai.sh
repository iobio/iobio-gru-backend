#!/bin/bash

url=$1
index_url=$2
samtools_region=$3

samtools_od view -b $url $samtools_region $index_url | head -n 1
