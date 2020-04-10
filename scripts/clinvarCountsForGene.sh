#!/bin/bash

clinvarUrl=$1
region=$2
binLength=$3
regionParts=$4
annotationMode=$5
requiresVepService=$6
vepArgs=$7

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
    tabix_od -h $clinvarUrl $region | \
        vep $vepArgs | \
            knownVariants_2 -r $region $binLengthArg $regionPartsArg $annotationModeArg
else
    tabix_od -h $clinvarUrl $region | \
        knownVariants_2 -r $region $binLengthArg $regionPartsArg $annotationModeArg
fi
