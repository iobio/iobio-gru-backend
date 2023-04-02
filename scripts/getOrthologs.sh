#!/bin/bash
set -euo pipefail

geneName=$1
eutilsKey=$2
speciesIds=$3
taxonLevelId=$4

tempDir=$(mktemp -d)
cd $tempDir

geneNameArg="-g $geneName"
eutilsKeyArg="-k $eutilsKey"
speciesIdArg="-s $speciesIds"
taxonLevelArg="-t $taxonLevelId"

get_orthologs $geneNameArg $eutilsKeyArg $speciesIdArg $taxonLevelArg

rm -rf $tempDir
