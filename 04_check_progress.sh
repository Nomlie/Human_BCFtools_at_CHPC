#!/bin/bash
#PBS -N check_progress
#PBS -q serial
#PBS -P CBBI1039
#PBS -l select=1:ncpus=1:mem=4gb
#PBS -l walltime=00:30:00

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
source "${PBS_O_WORKDIR}/pipeline.config"
source "${CONFIG_TXT}"

show_help() {
    cat << 'EOF'
04_check_progress.sh
---------------------

Summarizes chunk progress for the chromosome group defined in pipeline.conf.
EOF
}

if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
    show_help
    exit 0
fi

CHROM_GROUP_TAG=$(IFS=_; echo "${CHROM_LIST[*]}")

VCF_DIR="${vcf%/}"
GROUP_DIR="${VCF_DIR}/group_runs/${CHROM_GROUP_TAG}"
REPORT="${GROUP_DIR}/chunk_progress.tsv"

mkdir -p "${GROUP_DIR}"

{
    echo -e "chromosome_group\tchromosome\tfinished_chunks\tfailed_chunks\trunning_chunks"
    for CHR in "${CHROM_LIST[@]}"; do
        STAMPDIR="${GROUP_DIR}/checkpoints/${CHR}"

        if [ ! -d "${STAMPDIR}" ]; then
            echo -e "${CHROM_GROUP_TAG}\t${CHR}\t0\t0\t0"
            continue
        fi

        DONE=$(find "${STAMPDIR}" -maxdepth 1 -name "*.filter.done" | wc -l)
        FAILED=$(find "${STAMPDIR}" -maxdepth 1 -name "*.failed" | wc -l)
        RUNNING=0

        echo -e "${CHROM_GROUP_TAG}\t${CHR}\t${DONE}\t${FAILED}\t${RUNNING}"
    done
} > "${REPORT}"

echo "Chromosome group: ${CHROM_GROUP_TAG}"
echo "Wrote progress report: ${REPORT}"
cat "${REPORT}"