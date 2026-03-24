#!/bin/bash
#SBATCH --account=notchpeak-shared-short
#SBATCH --partition=notchpeak-shared-short
#SBATCH --qos=notchpeak-shared-short
#SBATCH --time=8:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --requeue


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
