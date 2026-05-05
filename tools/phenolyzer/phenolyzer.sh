#!/bin/bash

searchTerm=$1

perl /tool_bin/lib/phenolyzer/disease_annotation.pl \
    "$searchTerm" \
    -p -ph -logistic \
    -addon_gg DB_MENTHA_GENE_GENE_INTERACTION \
    -addon_gg_weight 0.05 \
    -out phenolyzer.out \
    > /dev/null

cat phenolyzer.out.final_gene_list
