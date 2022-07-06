#!/bin/bash
set -euo pipefail

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
hgvsNotation=${12}
getRsId=${13}
gnomadUrl=${14}
gnomadRegionFileStr=${15}
gnomadHeaderFile=${16}
decompose=${17}
gnomadRenameChr=${18}

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
    if [ "$genomeBuildName" == "GRCh38" ]; then
        annotsToAdd=CHROM,POS,REF,ALT,INFO/AF,INFO/AN,INFO/AC,INFO/nhomalt_raw,INFO/AF_popmax,INFO/AF-fin,INFO/AF-nfe,INFO/AF-oth,INFO/AF-amr,INFO/AF-afr,INFO/AF-asj,INFO/AF-eas,INFO/AF-sas
    fi 

    function gnomadAnnotFunc {
        if [ "$gnomadRenameChr" ]; then
                echo -e "$gnomadRenameChr" > gnomad_rename_chr.txt
                vt rminfo -t $annotsToRemove - | bcftools annotate --rename-chrs gnomad_rename_chr.txt | bgzip -c > gnomad.vcf.gz
        else
                vt rminfo -t $annotsToRemove - | bgzip -c > gnomad.vcf.gz
        fi        

        tabix gnomad.vcf.gz
        
        # Add the gnomAD INFO fields to the input vcf
        bcftools annotate -a $gnomadUrl -h $gnomadHeaderFile -c $annotsToAdd -R gnomad_regions.txt gnomad.vcf.gz
    }

    gnomadAnnotStage=gnomadAnnotFunc
fi


echo -e "$contigStr" > contigs.txt

vepBaseArgs="-i STDIN --format vcf --cache --dir_cache $vepCacheDir --offline --vcf -o STDOUT --no_stats --no_escape --sift b --polyphen b --regulatory --fork 4 --merged --fasta $refFastaFile"

vepArgs="$vepBaseArgs --assembly $genomeBuildName --allele_number"

if [ "$vepREVELFile" ]; then
    vepArgs="$vepArgs --dir_plugins $vepPluginDir --plugin REVEL,$vepREVELFile"
fi

if [ "$vepAF" == "true" ]; then
    vepArgs="$vepArgs --af --af_gnomad --af_esp --af_1kg --max_af"
fi

if [ "$hgvsNotation" == "true" ]; then
    vepArgs="$vepArgs --hgvs"
fi

if [ "$getRsId" == "true" ]; then
    vepArgs="$vepArgs --check_existing"
fi

tabixVcfArg=$vcfUrl
if [ -n "${tbiUrl}" ]; then
    tabixVcfArg="$vcfUrl##idx##$tbiUrl"
fi

tabix -h $tabixVcfArg $region | \
    bcftools annotate -h contigs.txt - | \
    $subsetStage | \
    $decomposeStage | \
    vt normalize -n -r $refFastaFile - | \
    vep $vepArgs | \
    $gnomadAnnotStage

#echo $tempDir
rm -rf $tempDir
cd $runDir
