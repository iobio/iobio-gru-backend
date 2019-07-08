#!/bin/bash


vcfUrl=$1
tbiUrl=$2
refName=$3
regions=$4
contigStr=$5
refFastaFile=$6

while true; do
    tmp_dir=/tmp/iobio_gene_coverage_$(cat /proc/sys/kernel/random/uuid)
    if [ ! -d $tmp_dir ]; then
        break;
    fi
done

mkdir $tmp_dir

# TODO: it doesn't seem to work unless I have a file extension on the end...
printf "$contigStr" > $tmp_dir/contigs.txt

#echo vcfUrl: $vcfUrl tbiUrl: $tbiUrl refName: $refName regions: $regions contigStr: $contigStr refFastaFile: $refFastaFile

tabix_od -h $vcfUrl $regions $tbiUrl | \
    bcftools annotate -h $tmp_dir/contigs.txt - | \
    vt normalize -q -n -r $refFastaFile -

rm -rf $tmp_dir
