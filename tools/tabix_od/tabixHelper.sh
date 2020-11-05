#!/bin/bash

runDir=$PWD
tempDir=$(mktemp -d)
cd $tempDir

last_arg=${@: -1}
tbi_url=""

# check if the last argument is an URL for index
if [[ $last_arg == http*tbi* ]]; then
    tbi_url=$last_arg
    vcf_basename=""
    args=()
    idx=0
    for i in $@; do
        if [[ $i == *.vcf.gz* && $i != *.vcf.gz.tbi* ]]; then
            # strip down to the base filename without any extension
            vcf_basename=${i##*/}
            vcf_basename=${vcf_basename%%.vcf*}
            # for tabix to work with a presigned URL that has 'content-disposition'
            # (contains '.'), we need to tell tabix the name of the local
            # index file that has been downloaded.  We do this by appending
            # ##idx##[csi file name goes here].
            args[$idx]="${i}##idx##${vcf_basename}.vcf.gz.csi"
    	else
            # Since we are appending ##idx## to the end of the presigned vcf
            # url, we can keep all other args, but drop the tbi URL as
            # it will not be used
            if [[ $i != *.vcf.gz.tbi* ]]; then
                args[$idx]=${i}
            fi
    	fi
        idx=$((idx+1))
    done
    # download the .tbi to basename.vcf.tbi
    curl -o ${vcf_basename}.vcf.tbi -sL $tbi_url

    # convert the .tbi to .csi
    tbi_to_csi ${vcf_basename}.vcf.tbi ${vcf_basename}.vcf

    # for DNANexus, we need the csi have the filename basename.vcf.gz.csi
    cp ${vcf_basename}.vcf.csi ${vcf_basename}.vcf.gz.csi

    tabix_args=${args[@]}
else
    tabix_args=$@
fi
tabix $tabix_args
tabixRetCode=$?

cd $runDir
rm -rf $tempDir

exit $tabixRetCode
