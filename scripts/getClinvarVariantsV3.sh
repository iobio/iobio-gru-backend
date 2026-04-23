#!/bin/bash
set -euo pipefail
#
# getClinvarVariantsV3.sh
#
# This is a special version of annotateVariantsVx.sh. Since Clinvar variants 
# are already annotated and don't have genotypes, they require different
# steps, more easily accomplished in a separate bash script. The steps are:
#  1) annotate clinvar variants with gnomad allele frequencies 
#  2) normalize the variants (i.e. convert to canonical representation)
#  3) filter the variants via an input arg. For example, return onlyvariants 
#     that are pathogenic or likely pathogenic.  
#
# Example usage:
# /iobio-gru-backend/scripts/getClinvarVariantsV3.sh \
#     "https://backend.iobio.io/static/clinvar/GRCh38/clinvar.vcf.gz" \
#     "https://backend.iobio.io/static/clinvar/GRCh38/clinvar.vcf.gz.tbi" \
#     "9:1979290-2194624" \
#     "1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y" \
#     "/gru_data/references/GRCh38/Homo_sapiens_assembly38.fasta" \
#     "GRCh38" \
#     "true" \
#     'INFO/CLNSIG="Pathogenic,Likely_pathogenic"'


vcfUrl=$1              # vcf url
tbiUrl=$2              # optional tbi url
region=$3              # region to extract (e.g. 9:1979290-2194624)
contigStr=$4           # contigs for header (e.g. chr1,...,chr22,chrX,chrY,chrM  or 1,2,X,Y,M)
refFastaFile=$5        # reference fasta file
genomeBuildName=$6     # genome build name (e.g. GRCh38, GRCh37)
gnomadMergeAnnots=$7   # true to annotate with gnomad, "" to skip
pathoFilterPhrase=${8} # bcftools filter phrase, e.g. "INFO/CLNSIG=\"Pathogenic,Likely_pathogenic\"" 


echo -e "$contigStr" > contigs.txt

tabixVcfArg=$vcfUrl
if [ -n "${tbiUrl}" ]; then
    tabixVcfArg="$vcfUrl##idx##$tbiUrl"
fi

#optional gnomad stage
gnomadAnnotStage=cat
if [ "$gnomadMergeAnnots" ]; then
    
    if [ "$genomeBuildName" == "GRCh38" ]; then
        toml="/gru_data/annotations/GRCh38/vcfanno_annotate_variants.toml"
    else
        toml="/gru_data/annotations/GRCh37/vcfanno_annotate_variants.toml"
    fi
    custom_lua="/gru_data/annotations/vcfanno_custom.lua"
    gnomadAnnotStage="vcfanno --lua $custom_lua $toml /dev/stdin"
fi

tabix -h $tabixVcfArg $region | \
    bcftools annotate -h contigs.txt - | \
    vt normalize -n -r $refFastaFile - | \
    bcftools filter -i $pathoFilterPhrase - | \
    $gnomadAnnotStage
