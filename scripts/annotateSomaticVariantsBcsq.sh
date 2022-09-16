#!/bin/bash
#SJG updated Sept2022
set -euo pipefail

vcfUrl=$1
selectedSamples=$2
regions=$3
somaticFilterPhrase=$4
refFastaFile=$5
gffFile=$6

runDir=$PWD
tempDir=$(mktemp -d)
cd $tempDir

# Add region filter stage if we have regions
regionFilterStage=cat
if [ ! -z "$regions" ]; then
    function regionFilterFunc {
        bcftools filter -t $regions -
    }
    regionFilterStage=regionFilterFunc
fi

#Add somatic filter stage if we have a filter
somFilterStage=cat
if [ ! -z "$somaticFilterPhrase" ]; then
    function somFilterFunc {
        bcftools filter -i $somaticFilterPhrase -
    }
    somFilterStage=somFilterFunc
fi

echo -e "$regions" >> regions.txt

#remove prefix from gff file
#gffFile="/iobio-gru-backend/static/ensembl/GRCh38/geneSet38.gff3.gz"
#if [ "$genomeBuild" == "GRCh38" ]; then
#    cat $gffFile | gzip -d - | \
#    perl -pe 's/^([0-9]+|[X]|[Y]|[M])/chr$1/' - | \
#    gzip - >> noPrefix.gff3.gz
#    gffFile="/iobio-gru-backend/static/ensembl/GRCh37/geneSet37.gff3.gz"
#else
#    cat $gffFile >> noPrefix.gff3.gz
#fi

#cat $gffFile | gzip -d -

#Do work
bcftools view -s $selectedSamples $vcfUrl | \
    $regionFilterStage | \
    bcftools norm -m - -w 10000 -f $refFastaFile - | \
    $somFilterStage | \
    bcftools csq -f $refFastaFile -g $gffFile - -Ov -l

rm -rf $tempDir
cd $runDir

