#!/bin/bash

VEP_VERSION=101.0

outDir=$1

if [ -z $outDir ];then
    echo "Usage: $0 OUT_DIR"
    exit 1
fi

singularity pull ${outDir}/vep.sif docker://ensemblorg/ensembl-vep:release_${VEP_VERSION}

cp vep.sh ${outDir}/vep
