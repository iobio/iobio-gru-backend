#!/bin/bash
set -euo pipefail

vcfUrl=$1
region=$2

runDir=$PWD
tempDir=$(mktemp -d)
cd $tempDir

bcftools view -r $region $vcfUrl

#echo $tempDir
rm -rf $tempDir
cd $runDir

