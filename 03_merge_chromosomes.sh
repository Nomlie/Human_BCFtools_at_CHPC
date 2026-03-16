#!/bin/bash
#PBS -N merge_chroms
#PBS -q seriallong
#PBS -P CBBI1039
#PBS -l select=1:ncpus=8:mem=32gb
#PBS -l walltime=24:00:00

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
source "${SCRIPT_DIR}/pipeline.conf"
source "${CONFIG_TXT}"
source "${CONDA_ACTIVATE}"

show_help() {
    cat << 'EOF'
03_merge_chromosomes.sh
------------------------

Merges chunked VCFs into one final VCF per chromosome for the chromosome group
defined in pipeline.conf.
EOF
}

if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
    show_help
    exit 0
fi

CHROM_GROUP_TAG=$(IFS=_; echo "${CHROM_LIST[*]}")

VCF_DIR="${vcf%/}"
GROUP_DIR="${VCF_DIR}/group_runs/${CHROM_GROUP_TAG}"
FINAL_DIR="${GROUP_DIR}/final"
mkdir -p "${FINAL_DIR}"

echo "Starting chromosome merges"
echo "Chromosome group: ${CHROM_GROUP_TAG}"
echo "Group directory: ${GROUP_DIR}"
echo "Start=$(date)"

for CHR in "${CHROM_LIST[@]}"; do
    CHUNK_DIR="${GROUP_DIR}/chunks/${CHR}"
    OUT_VCF="${FINAL_DIR}/HCC.${CHR}.bcftools.annotated.filtered.vcf.gz"
    OUT_DONE="${FINAL_DIR}/HCC.${CHR}.merge.done"

    if [ -f "${OUT_DONE}" ] && [ -f "${OUT_VCF}" ] && [ -f "${OUT_VCF}.tbi" ]; then
        echo "Skipping ${CHR}: already merged"
        continue
    fi

    if [ ! -d "${CHUNK_DIR}" ]; then
        echo "Skipping ${CHR}: no chunk directory"
        continue
    fi

    LIST=$(mktemp)
    find "${CHUNK_DIR}" -maxdepth 1 -type f -name "HCC.${CHR}.*.filtered.vcf.gz" | sort -V > "${LIST}"

    if [ ! -s "${LIST}" ]; then
        echo "Skipping ${CHR}: no filtered chunk files"
        rm -f "${LIST}"
        continue
    fi

    rm -f "${OUT_VCF}" "${OUT_VCF}.tbi" "${OUT_DONE}"

    bcftools concat -f "${LIST}" -Oz -o "${OUT_VCF}"
    bcftools index -t "${OUT_VCF}"
    bcftools view -h "${OUT_VCF}" > /dev/null

    date "+%Y-%m-%d %H:%M:%S" > "${OUT_DONE}"
    rm -f "${LIST}"
    echo "Merged ${CHR}"
done

echo "Finished chromosome merges"
echo "Final outputs: ${FINAL_DIR}"
echo "End=$(date)"