#!/bin/bash
set -euo pipefail

vcfUrl=$1
tbiUrl=$2
region=$3
contigStr=$4
refFastaFile=$5
genomeBuildName=$6
gnomadUrl=$7
gnomadRegionFileStr=$8
gnomadHeaderFile=$9
gnomadRenameChr=${10}
pathoFilterPhrase=${11}

echo -e "$contigStr" > contigs.txt

tabixVcfArg=$vcfUrl
if [ -n "${tbiUrl}" ]; then
    tabixVcfArg="$vcfUrl##idx##$tbiUrl"
fi

#optional gnomad stage
gnomadAnnotStage=cat
if [ "$gnomadUrl" ]; then
    printf "$gnomadRegionFileStr" > gnomad_regions.txt

    # These are the INFO fields to clear out
    annotsToRemove=AF,AN,AC

    # These are the gnomAD INFO fields to add to the input vcf
    annotsToAdd=CHROM,POS,REF,ALT,INFO/AF,INFO/AN,INFO/AC,INFO/nhomalt_raw,INFO/AF_popmax,INFO/AF_fin,INFO/AF_nfe,INFO/AF_oth,INFO/AF_amr,INFO/AF_afr,INFO/AF_asj,INFO/AF_eas,INFO/AF_sas


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


tabix -h $tabixVcfArg $region | \
    bcftools annotate -h contigs.txt - | \
    vt normalize -n -r $refFastaFile - | \
    bcftools filter -i $pathoFilterPhrase - | \
    $gnomadAnnotStage
