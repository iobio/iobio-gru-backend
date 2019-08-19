#!/bin/bash

while true; do
    uuid=$(cat /proc/sys/kernel/random/uuid)
    if [ ! -d /tmp/$uuid ]; then
        break;
    fi
done

DIR=$PWD

mkdir -p /tmp/$uuid
cd /tmp/$uuid

last_arg=${@: -1}

# check if the last argument is an URL for index
# for crai files, we want to download it and
# name it the simple name of basename.crai
#  example:  a-test-alignment.crai
if [[ $last_arg == http*crai* ]]; then
    cram_basename="";
    for i in $@; do
        if [[ $i == *.cram* && $i != *.crai* ]]; then
            cram_basename=${i##*/}
            cram_basename=${cram_basename%%.cram*}
        fi
    done

    # download the crai file and name it the simple name
    curl -o ${cram_basename}.crai -sL $last_arg

    # remove the crai url from the args
    len=$(expr $# - 1)
    samtools_args=${@:1:$len}
# check if the last argument is an URL for the bai index
elif [[ $last_arg == http*bai* ]]; then
    bam_basename="";
    bam_longname="";
    for i in $@; do
        if [[ $i == *.bam* && $i != *.bai* ]]; then
            bam_longname=${i##*/}
            bam_basename=${i##*/}
            bam_basename=${bam_basename%%.bam*}
        fi
    done


    curl -o ${bam_basename}.bam.bai -sL $last_arg
    bai_to_csi ${bam_basename}.bam.bai ${bam_basename}.bam

    # signed URLs need a special from of the .csi filename which is the filename + args, with
    # extension .csi
    if [ ! -f ${bam_longname}.csi ]; then
        cp ${bam_basename}.bam.csi ${bam_longname}.csi
    fi

    len=$(expr $# - 1)
    samtools_args=${@:1:$len}
else
    samtools_args=$@
fi

export REF_PATH=$DIR/data/md5_reference_cache/%2s/%2s/%s:http://www.ebi.ac.uk/ena/cram/md5/%s
export REF_CACHE=$DIR/data/md5_reference_cache/%2s/%2s/%s

samtools $samtools_args
samtoolsRetCode=$?

cd ..
rm -rf /tmp/$uuid
exit $samtoolsRetCode 
