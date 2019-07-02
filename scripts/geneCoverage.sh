#!/bin/bash

url=$1
index_url=$2
samtools_region=$3
region_file_str=$4

# TODO: should probably re-enable this, but it was going into an endless loop
#while true; do
#    region_file_path=/tmp/iobio_region_file_$(cat /proc/sys/kernel/random/uuid).txt
#    echo ${region_file_path}
#    if [ ! -f ${region_file_path} ]; then
#        break;
#    fi
#done

region_file_path=/tmp/iobio_region_file_$(cat /proc/sys/kernel/random/uuid)
# TODO: it doesn't seem to work unless I have the .txt on the end...
printf "$region_file_str" > ${region_file_path}.txt

samtools_od view -b $url $samtools_region $index_url > ${region_file_path}.aln

samtools index ${region_file_path}.aln

#geneCoverage -b <(samtools index $url $samtools_region $index_url) -r $region_file_path.txt
geneCoverage -b $region_file_path.aln -r $region_file_path.txt

rm -rf $region_file_path
