#!/bin/bash

clinvarUrl=$1
region=$2
binLength=$3
regionParts=$4

if [ "$binLength" ]; then
    binLengthArg="-b $binLength"
elif [ "$regionParts" ]; then
    regionPartsArg="-p $regionParts"
fi

tabix_od -h $clinvarUrl $region | \
    knownVariants -r $region $binLengthArg $regionPartsArg
