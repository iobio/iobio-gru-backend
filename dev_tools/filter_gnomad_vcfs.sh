#!/bin/bash
#SBATCH --account=marth
#SBATCH --partition=notchpeak-freecycle
#SBATCH --qos=notchpeak-freecycle
#SBATCH --time=04:00:00
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=4
#SBATCH --mem=4G
#SBATCH --array=1-24
#SBATCH --requeue

set -euo pipefail

gnomad_dir="${1}"
bcftools="${2}"

case "${SLURM_ARRAY_TASK_ID}" in
    1|2|3|4|5|6|7|8|9|10|11|12|13|14|15|16|17|18|19|20|21|22) chr="${SLURM_ARRAY_TASK_ID}" ;;
    23) chr="X" ;;
    24) chr="Y" ;;
esac

in_file="${gnomad_dir}/gnomad.genomes.v4.1.sites.chr${chr}.vcf.bgz"
out_file="gnomad.genomes.v4.1.sites.chr${chr}.iobio.vcf.bgz"
out_file_tmp=tmp.${out_file}

echo "Processing ${in_file}"

if [ ! -f "${out_file}" ]; then
    echo "Generating ${out_file}"

    ${bcftools} annotate --threads 4 -x \
    '^INFO/AF_mid,^INFO/AF_asj,^INFO/AF_ami,^INFO/AF_fin,^INFO/AF_nfe,^INFO/AF_amr,^INFO/AF_afr,^INFO/AF_eas,^INFO/AF_sas,^INFO/AF,^INFO/AF_remaining,^INFO/nhomalt,^INFO/AC_XY,^INF
    O/AC,^INFO/AN,^INFO/AF_grpmax,^INFO/VarDP,^INFO/non_par' \
    -O z -o ${out_file_tmp} \
    ${in_file}

    mv ${out_file_tmp} ${out_file}
fi

if [ ! -f "${out_file}.tbi" ]; then
    echo "Generating ${out_file}.tbi"

    ${bcftools} index --tbi --threads 4 ${out_file} -o tmp.${out_file}.tbi
    mv tmp.${out_file}.tbi ${out_file}.tbi
fi
