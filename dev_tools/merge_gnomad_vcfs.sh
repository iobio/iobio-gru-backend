#/bin/bash

set -euo pipefail

bcftools="${1}"

gnomad_version="4.1"

command="${bcftools} merge -O z --write-index=tbi --threads 8 -o gnomad.genomes.v${gnomad_version}.iobio.vcf.bgz "
for i in {1..22}; do
    command=$command"gnomad.genomes.v${gnomad_version}.sites.chr"${i}".iobio.vcf.bgz "
done

command=$command"gnomad.genomes.v${gnomad_version}.sites.chrX.iobio.vcf.bgz "
command=$command"gnomad.genomes.v${gnomad_version}.sites.chrY.iobio.vcf.bgz "

echo $command

$command
