#!/bin/bash
set -euo pipefail

clinvarUrl=$1
region=$2
binLength=$3
regionParts=$4
annotationMode=$5
requiresVepService=$6
vepArgs=$7

tempDir=$(mktemp -d)
cd $tempDir

binLengthArg=
regionPartsArg=
annotationModeArg=

if [ "$binLength" ]; then
    binLengthArg="-b $binLength"
elif [ "$regionParts" ]; then
    regionPartsArg="-p $regionParts"
fi

if [ "$annotationMode" ]; then
    annotationModeArg="-m $annotationMode"
fi

# Pipe into VEP if we want to return counts by VEP categories but haven't already annotated it
if [ "$requiresVepService" = true ]; then
    tabix -h $clinvarUrl $region | \
        vep $vepArgs | \
            knownVariants -r $region $binLengthArg $regionPartsArg $annotationModeArg
else
    tabix -h $clinvarUrl $region | \
        knownVariants -r $region $binLengthArg $regionPartsArg $annotationModeArg
fi

rm -rf $tempDir
