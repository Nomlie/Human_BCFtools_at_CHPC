#!/bin/bash
#PBS -N make_regions
#PBS -q serial
#PBS -P CBBI1039
#PBS -l select=1:ncpus=1:mem=4gb
#PBS -l walltime=01:00:00
#PBS -m abe
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
source "${SCRIPT_DIR}/pipeline.conf"
source "${CONFIG_TXT}"

show_help() {
    cat << 'EOF'
01_make_regions.sh
-------------------

Creates a regions.tsv file for the chromosome group defined in pipeline.conf.

This script:
- reads hg38.fa.fai
- keeps only chromosomes listed in CHROM_LIST
- splits them into fixed-size chunks
- writes:
    <VCF_DIR>/group_runs/<chrom_group_tag>/regions.tsv
EOF
}

if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
    show_help
    exit 0
fi

FAI="${fastq%/}/hg38.fa.fai"
VCF_DIR="${vcf%/}"

CHROM_GROUP_TAG=$(IFS=_; echo "${CHROM_LIST[*]}")
GROUP_DIR="${VCF_DIR}/group_runs/${CHROM_GROUP_TAG}"
OUT="${GROUP_DIR}/regions.tsv"

mkdir -p "${GROUP_DIR}"

if [ ! -f "${FAI}" ]; then
    echo "ERROR: FASTA index not found: ${FAI}"
    exit 1
fi

CHROM_REGEX="^($(printf "%s|" "${CHROM_LIST[@]}" | sed 's/|$//'))$"

{
    printf "chrom\tstart\tend\tregion\n"
    awk -v chunk="${CHUNK_SIZE}" -v chrom_re="${CHROM_REGEX}" '
    BEGIN { OFS="\t" }
    $1 ~ chrom_re {
        chrom = $1
        len = $2
        for (start = 1; start <= len; start += chunk) {
            end = start + chunk - 1
            if (end > len) end = len
            print chrom, start, end, chrom ":" start "-" end
        }
    }
    ' "${FAI}"
} > "${OUT}"

TOTAL_REGIONS=$(( $(wc -l < "${OUT}") - 1 ))

echo "Chromosome group: ${CHROM_GROUP_TAG}"
echo "Wrote regions file: ${OUT}"
echo "Total regions: ${TOTAL_REGIONS}"
echo "Chunk size: ${CHUNK_SIZE}"