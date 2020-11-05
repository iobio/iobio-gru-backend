#!/bin/bash

cacheDir=$1

singularity exec vep.sif vep \
    -i STDIN --format vcf --cache --dir_cache $cacheDir --offline \
    --vcf -o STDOUT --no_stats --no_escape --sift b --polyphen b --regulatory \
    --fork 4 $@
