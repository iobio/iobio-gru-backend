#!/bin/bash

while true; do
    uuid=$(cat /proc/sys/kernel/random/uuid)
    if [ ! -d $uuid ]; then
        break;
    fi
done

mkdir -p $uuid
cd $uuid

last_arg=${@: -1}
tbi_url=""

# check if the last argument is an URL for index
if [[ $last_arg == http*tbi* ]]; then
    tbi_url=$last_arg
    vcf_basename=""
    for i in $@; do
        if [[ $i == *.vcf.gz* && $i != *.vcf.gz.tbi* ]]; then
            # strip down to the base filename without any extension
            vcf_basename=${i##*/}
            vcf_basename=${vcf_basename%%.vcf*}
        fi
    done
    # download the .tbi to basename.vcf.tbi
    curl -o ${vcf_basename}.vcf.tbi -sL $tbi_url

    # convert the .tbi to .csi
    tbi_to_csi ${vcf_basename}.vcf.tbi ${vcf_basename}.vcf

    # for DNANexus, we need the csi have the filename basename.vcf.gz.csi
    cp ${vcf_basename}.vcf.csi ${vcf_basename}.vcf.gz.csi

    len=$(expr $# - 1)
    tabix_args=${@:1:$len}
else
    tabix_args=$@
fi

tabix $tabix_args
tabixRetCode=$?
cd ..
rm -rf $uuid
exit $tabixRetCode
