#!/bin/bash

# Incantation taken from https://stackoverflow.com/a/246128/943814
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

singularity exec $SCRIPT_DIR/vep.sif vep \
    -i STDIN --format vcf --cache --dir_cache $SCRIPT_DIR/../vep-cache \
    --offline --vcf -o STDOUT --no_stats --no_escape --sift b --polyphen b \
    --regulatory --fork 4 $@
