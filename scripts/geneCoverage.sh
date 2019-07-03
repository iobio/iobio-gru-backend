#!/bin/bash

url=$1
index_url=$2
samtools_region=$3
region_file_str=$4

while true; do
    tmp_dir=/tmp/iobio_gene_coverage_$(cat /proc/sys/kernel/random/uuid)
    if [ ! -d $tmp_dir ]; then
        break;
    fi
done

mkdir $tmp_dir

# TODO: it doesn't seem to work unless I have a file extension on the end...
printf "$region_file_str" > $tmp_dir/region_file.txt

samtools_od view -b $url $samtools_region $index_url > $tmp_dir/alignment.aln

samtools index $tmp_dir/alignment.aln

geneCoverage -b $tmp_dir/alignment.aln -r $tmp_dir/region_file.txt

rm -rf $tmp_dir
