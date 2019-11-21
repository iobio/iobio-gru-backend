#!/bin/bash

vcfUrl=$1
tbiUrl=$2
region=$3
contigStr=$4
vcfSampleNamesStr=$5
refFastaFile=$6
genomeBuildName=$7
vepCacheDir=$8
vepREVELFile=$9
vepAF=${10}
vepPluginDir=${11}
isRefSeq=${12}
hgvsNotation=${13}
getRsId=${14}
gnomadUrl=${15}
gnomadRegionFileStr=${16}
gnomadHeaderFile=${17}
decompose=${18}

# default optional stages to no-op
subsetStage=cat
gnomadAnnotStage=cat
decomposeStage=cat

runDir=$PWD
tempDir=$(mktemp -d)
cd $tempDir

if [ "$vcfSampleNamesStr" ]; then
    echo -e "$vcfSampleNamesStr" > samples.txt

    function subsetFunc {
        vt subset -s samples.txt -
    }

    subsetStage=subsetFunc
fi


if [ "$decompose" == "true" ]; then
    function decomposeFunc {
        vt decompose -s  -
    }

    decomposeStage=decomposeFunc
fi



if [ "$gnomadUrl" ]; then
    printf "$gnomadRegionFileStr" > gnomad_regions.txt
    
    # These are the INFO fields to clear out
    annotsToRemove=AF,AN,AC
    
    # These are the gnomAD INFO fields to add to the input vcf
    annotsToAdd=CHROM,POS,REF,ALT,INFO/AF,INFO/AN,INFO/AC,INFO/nhomalt_raw,INFO/AF_popmax,INFO/AF_fin,INFO/AF_nfe,INFO/AF_oth,INFO/AF_amr,INFO/AF_afr,INFO/AF_asj,INFO/AF_eas,INFO/AF_sas
    

    function gnomadAnnotFunc {
        vt rminfo -t $annotsToRemove - | bgzip -c > gnomad.vcf.gz
        tabix gnomad.vcf.gz
        
        # Add the gnomAD INFO fields to the input vcf
        bcftools annotate -a $gnomadUrl -h $gnomadHeaderFile -c $annotsToAdd -R gnomad_regions.txt gnomad.vcf.gz
    }

    gnomadAnnotStage=gnomadAnnotFunc
fi


echo -e "$contigStr" > contigs.txt


# TODO: remove --dir_cache from vep Dockerfile since we're overriding it here.
# Actually, move all the vep args into this script so they aren't split
# between two locations.
vepArgs="--assembly $genomeBuildName --format vcf --allele_number --dir_cache $vepCacheDir"

if [ "$vepREVELFile" ]; then
    vepArgs="$vepArgs --dir_plugins $vepPluginDir --plugin REVEL,$vepREVELFile"
fi

if [ "$vepAF" == "true" ]; then
    vepArgs="$vepArgs --af --af_gnomad --af_esp --af_1kg --max_af"
fi

if [ "$isRefSeq" == "true" ]; then
    vepArgs="$vepArgs --refseq"
fi

if [ "$hgvsNotation" == "true" ]; then
    vepArgs="$vepArgs --hgvs"
fi

if [ "$getRsId" == "true" ]; then
    vepArgs="$vepArgs --check_existing"
fi

if [ "$hgvsNotation" == "true" ] || [ "$getRsId" == "true" ] || [ "$isRefSeq" == "true" ]; then
    vepArgs="$vepArgs --fasta $refFastaFile"
fi

tabix_od -h $vcfUrl $region $tbiUrl | \
    bcftools annotate -h contigs.txt - | \
    $subsetStage | \
    $decomposeStage | \
    vt normalize -n -r $refFastaFile - | \
    vep $vepArgs | \
    $gnomadAnnotStage

#echo $tempDir
rm -rf $tempDir
cd $runDir
