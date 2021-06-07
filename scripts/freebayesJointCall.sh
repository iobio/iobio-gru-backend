#!/bin/bash

#set -x

alignmentUrls=$1
alignmentIndices=$2
region=$3
refFastaFile=$4
useSuggestedVariants=$5
clinvarUrl=$6
genomeBuildName=$7
vepREVELFile=$8
vepAF=$9
samplesFileStr=${10}
extraArgs=${11}
vepCacheDir=${12}
vepPluginDir=${13}
gnomadUrl=${14}
gnomadRegionFileStr=${15}
gnomadHeaderFile=${16}
decompose=${17}
dataDir=${18}

runDir=$PWD
tempDir=$(mktemp -d)
cd $tempDir

contigFile="contig.txt"
echo -e "$contigStr" > $contigFile

samplesFile="samples.txt"
echo -e "$samplesFileStr" > $samplesFile

export REF_CACHE=$dataDir/md5_reference_cache/%2s/%2s/%s

freebayesArgs="-s $samplesFile"

# split alignments by ','
IFS=','
read -ra urls <<< "$alignmentUrls"
read -ra indices <<< "$alignmentIndices"
for i in "${!urls[@]}"; do
    url=${urls[$i]}
    indexUrl=${indices[$i]}
    alignmentFile=$(mktemp)

    if [ -n "${indexUrl}" ]; then
        samtools-1.11 view -b -X $url $indexUrl $region > $alignmentFile &
    else
        samtools-1.11 view -b $url $region > $alignmentFile &
    fi

    freebayesArgs="$freebayesArgs -b $alignmentFile"
done
IFS=' '

tabixCommand=''
if [ "$useSuggestedVariants" == "true" ]; then
    suggFile="sugg.vcf"
    tabix -h $clinvarUrl $region | vt view -f "INFO.CLNSIG=~'5|4'" - > $suggFile
    freebayesArgs="$freebayesArgs -@ $suggFile"
fi

decomposeStage=''
if [ "$decompose" == "true" ]; then
    function decomposeFunc {
        vt decompose -s  -
    }
    decomposeStage=decomposeFunc
fi


gnomadAnnotStage=''
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

freebayesArgs="$freebayesArgs $extraArgs"

vepBaseArgs="-i STDIN --format vcf --cache --dir_cache $vepCacheDir --offline --vcf -o STDOUT --no_stats --no_escape --sift b --polyphen b --regulatory --fork 4 --merged --fasta $refFastaFile"

vepArgs="$vepBaseArgs --assembly $genomeBuildName --allele_number --hgvs --check_existing"

if [ "$vepREVELFile" ]; then
    vepArgs="$vepArgs --dir_plugins $vepPluginDir --plugin REVEL,$vepREVELFile"
fi

if [ "$vepAF" == "true" ]; then
    vepArgs="$vepArgs --af --af_gnomad --af_esp --af_1kg --max_af"
fi

wait

# TODO: vt filter -d used to be "Variants called by iobio" when invoked by
# minion, but for some reason I can't get it to accept more than a single word
# directly in the bash script
freebayes -f $refFastaFile $freebayesArgs | \
    $decomposeStage | \
    vt normalize -r $refFastaFile - | \
    vt filter -f 'QUAL>1' -t 'PASS' -d 'iobio' - | \
    bcftools annotate -h $contigFile | \
    vep $vepArgs | \
    $gnomadAnnotStage


rm -rf $tempDir
cd $runDir
