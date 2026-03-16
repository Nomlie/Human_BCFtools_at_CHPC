#!/bin/bash
#PBS -N bcft_group_parallel
#PBS -q smp
#PBS -P CBBI1039
#PBS -l select=1:ncpus=24:mpiprocs=24:mem=120gb
#PBS -l walltime=96:00:00

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
source "${SCRIPT_DIR}/pipeline.conf"
source "${CONFIG_TXT}"
source "${CONDA_ACTIVATE}"

show_help() {
    cat << 'EOF'
02_run_chunks.sh
-----------------

Runs the chromosome group defined in pipeline.conf using internal parallelism.

Design:
- one PBS smp job requests 24 CPUs
- up to 3 chromosomes run at the same time inside the node
- each chromosome gets 8 threads
- each chromosome is processed chunk by chunk from regions.tsv
- completed chunks are skipped automatically
- outputs go under:
    <VCF_DIR>/group_runs/<chrom_group_tag>/

Special chromosome handling:
- chrM is called as haploid
- chrY is called as haploid
- chrX currently uses the default calling mode
EOF
}

if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
    show_help
    exit 0
fi

VCF_DIR="${vcf%/}"
FASTQ_DIR="${fastq%/}"
DBSNP_DIR="${dbsnp%/}"
REF="${FASTQ_DIR}/hg38.fa"

CHROM_GROUP_TAG=$(IFS=_; echo "${CHROM_LIST[*]}")
GROUP_DIR="${VCF_DIR}/group_runs/${CHROM_GROUP_TAG}"
REGIONS_FILE="${GROUP_DIR}/regions.tsv"

AUTO_CLEAN_INCOMPLETE="${AUTO_CLEAN_INCOMPLETE:-0}"

if [ ! -f "${REGIONS_FILE}" ]; then
    echo "ERROR: Regions file not found: ${REGIONS_FILE}"
    echo "Run setup first."
    exit 1
fi

MASTER_LOG_DIR="${GROUP_DIR}/master_logs"
mkdir -p "${MASTER_LOG_DIR}"
MASTER_LOG="${MASTER_LOG_DIR}/run_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "${MASTER_LOG}") 2>&1

format_duration() {
    local total_seconds=$1
    local days=$((total_seconds / 86400))
    local hours=$(((total_seconds % 86400) / 3600))
    local minutes=$(((total_seconds % 3600) / 60))
    local seconds=$((total_seconds % 60))

    if [ "${days}" -gt 0 ]; then
        printf "%dd %02dh %02dm %02ds" "${days}" "${hours}" "${minutes}" "${seconds}"
    elif [ "${hours}" -gt 0 ]; then
        printf "%02dh %02dm %02ds" "${hours}" "${minutes}" "${seconds}"
    elif [ "${minutes}" -gt 0 ]; then
        printf "%02dm %02ds" "${minutes}" "${seconds}"
    else
        printf "%02ds" "${seconds}"
    fi
}

validate_bcf() {
    local f="$1"
    [ -f "${f}" ] && [ -s "${f}" ] && bcftools view -h "${f}" > /dev/null 2>&1
}

validate_vcfgz() {
    local f="$1"
    [ -f "${f}" ] && [ -s "${f}" ] && bcftools view -h "${f}" > /dev/null 2>&1
}

call_variants_for_chrom() {
    local chr="$1"
    local raw_bcf="$2"

    if [ "${chr}" = "chrY" ] || [ "${chr}" = "chrM" ]; then
        bcftools call \
            --threads "${THREADS_PER_CHROM}" \
            --ploidy 1 \
            -mv \
            -Ob \
            -o "${raw_bcf}"
    else
        bcftools call \
            --threads "${THREADS_PER_CHROM}" \
            -mv \
            -Ob \
            -o "${raw_bcf}"
    fi
}

run_chunk() {
    local chr="$1"
    local start="$2"
    local end="$3"
    local region="$4"

    local outdir="${GROUP_DIR}/chunks/${chr}"
    local logdir="${GROUP_DIR}/logs/${chr}"
    local stampdir="${GROUP_DIR}/checkpoints/${chr}"
    local summarydir="${GROUP_DIR}/summaries/${chr}"

    mkdir -p "${outdir}" "${logdir}" "${stampdir}" "${summarydir}"

    local prefix="${outdir}/HCC.${chr}.${start}_${end}"

    local raw_bcf="${prefix}.bcf"
    local raw_csi="${raw_bcf}.csi"
    local ann_bcf="${prefix}.annotated.bcf"
    local ann_csi="${ann_bcf}.csi"
    local flt_vcf="${prefix}.filtered.vcf.gz"
    local flt_tbi="${flt_vcf}.tbi"

    local raw_done="${stampdir}/HCC.${chr}.${start}_${end}.mpileup_call.done"
    local ann_done="${stampdir}/HCC.${chr}.${start}_${end}.annotate.done"
    local flt_done="${stampdir}/HCC.${chr}.${start}_${end}.filter.done"

    local summary_tsv="${summarydir}/HCC.${chr}.${start}_${end}.timing.tsv"
    local status_tsv="${summarydir}/HCC.${chr}.${start}_${end}.status.tsv"
    local logfile="${logdir}/HCC.${chr}.${start}_${end}.log"

    {
        echo "=================================================="
        echo "Chromosome group: ${CHROM_GROUP_TAG}"
        echo "Chromosome: ${chr}"
        echo "Region: ${region}"
        echo "Threads per chromosome: ${THREADS_PER_CHROM}"
        echo "Start: $(date)"
        echo "=================================================="

        [ -f "${summary_tsv}" ] || echo -e "region\tstep\tstart_epoch\tend_epoch\tseconds\tformatted_duration\tstatus" > "${summary_tsv}"
        [ -f "${status_tsv}" ] || echo -e "region\tstep\tstatus\ttimestamp" > "${status_tsv}"

        append_timing() {
            local step="$1"
            local start_epoch="$2"
            local end_epoch="$3"
            local status="$4"
            local elapsed=$((end_epoch - start_epoch))
            echo -e "${region}\t${step}\t${start_epoch}\t${end_epoch}\t${elapsed}\t$(format_duration "${elapsed}")\t${status}" >> "${summary_tsv}"
        }

        append_status() {
            local step="$1"
            local status="$2"
            echo -e "${region}\t${step}\t${status}\t$(date "+%Y-%m-%d %H:%M:%S")" >> "${status_tsv}"
        }

        region_start_ts=$(date +%s)

        if [ -f "${raw_done}" ] && [ -f "${raw_csi}" ] && validate_bcf "${raw_bcf}"; then
            echo "[SKIP] mpileup+call already completed"
            append_status "mpileup_call" "skipped"
        else
            if { [ -f "${raw_bcf}" ] || [ -f "${raw_csi}" ]; } && [ "${AUTO_CLEAN_INCOMPLETE}" = "1" ]; then
                rm -f "${raw_bcf}" "${raw_csi}" "${raw_done}"
            fi

            if { [ -f "${raw_bcf}" ] || [ -f "${raw_csi}" ]; } && [ ! -f "${raw_done}" ] && [ "${AUTO_CLEAN_INCOMPLETE}" != "1" ]; then
                echo "ERROR: Incomplete mpileup_call output for ${region}"
                echo "Delete ${raw_bcf} ${raw_csi} ${raw_done} and rerun, or use AUTO_CLEAN_INCOMPLETE=1"
                exit 1
            fi

            step_start=$(date +%s)
            append_status "mpileup_call" "started"
            rm -f "${raw_bcf}" "${raw_csi}" "${raw_done}"

            bcftools mpileup \
                --threads "${THREADS_PER_CHROM}" \
                --max-depth "${MAX_DEPTH}" \
                -Ou \
                -r "${region}" \
                -f "${REF}" \
                --bam-list "${BAM_LIST}" \
            | call_variants_for_chrom "${chr}" "${raw_bcf}"

            bcftools index "${raw_bcf}"
            validate_bcf "${raw_bcf}"
            touch "${raw_done}"

            step_end=$(date +%s)
            append_status "mpileup_call" "done"
            append_timing "mpileup_call" "${step_start}" "${step_end}" "done"
            echo "[DONE] mpileup+call :: ${region}"
            echo "[TIME] mpileup+call: $(format_duration "$((step_end - step_start))")"
        fi

        if [ -f "${ann_done}" ] && [ -f "${ann_csi}" ] && validate_bcf "${ann_bcf}"; then
            echo "[SKIP] annotate already completed"
            append_status "annotate" "skipped"
        else
            if { [ -f "${ann_bcf}" ] || [ -f "${ann_csi}" ]; } && [ "${AUTO_CLEAN_INCOMPLETE}" = "1" ]; then
                rm -f "${ann_bcf}" "${ann_csi}" "${ann_done}"
            fi

            if { [ -f "${ann_bcf}" ] || [ -f "${ann_csi}" ]; } && [ ! -f "${ann_done}" ] && [ "${AUTO_CLEAN_INCOMPLETE}" != "1" ]; then
                echo "ERROR: Incomplete annotate output for ${region}"
                echo "Delete ${ann_bcf} ${ann_csi} ${ann_done} and rerun, or use AUTO_CLEAN_INCOMPLETE=1"
                exit 1
            fi

            step_start=$(date +%s)
            append_status "annotate" "started"
            rm -f "${ann_bcf}" "${ann_csi}" "${ann_done}"

            bcftools annotate \
                --threads "${THREADS_PER_CHROM}" \
                -a "${DBSNP_DIR}/dbsnp_146.hg38.vcf.gz" \
                -c ID \
                -Ob \
                -o "${ann_bcf}" \
                "${raw_bcf}"

            bcftools index "${ann_bcf}"
            validate_bcf "${ann_bcf}"
            touch "${ann_done}"

            step_end=$(date +%s)
            append_status "annotate" "done"
            append_timing "annotate" "${step_start}" "${step_end}" "done"
            echo "[DONE] annotate :: ${region}"
            echo "[TIME] annotate: $(format_duration "$((step_end - step_start))")"
        fi

        if [ -f "${flt_done}" ] && [ -f "${flt_tbi}" ] && validate_vcfgz "${flt_vcf}"; then
            echo "[SKIP] filter already completed"
            append_status "filter" "skipped"
        else
            if { [ -f "${flt_vcf}" ] || [ -f "${flt_tbi}" ]; } && [ "${AUTO_CLEAN_INCOMPLETE}" = "1" ]; then
                rm -f "${flt_vcf}" "${flt_tbi}" "${flt_done}"
            fi

            if { [ -f "${flt_vcf}" ] || [ -f "${flt_tbi}" ]; } && [ ! -f "${flt_done}" ] && [ "${AUTO_CLEAN_INCOMPLETE}" != "1" ]; then
                echo "ERROR: Incomplete filter output for ${region}"
                echo "Delete ${flt_vcf} ${flt_tbi} ${flt_done} and rerun, or use AUTO_CLEAN_INCOMPLETE=1"
                exit 1
            fi

            step_start=$(date +%s)
            append_status "filter" "started"
            rm -f "${flt_vcf}" "${flt_tbi}" "${flt_done}"

            bcftools filter \
                --threads "${THREADS_PER_CHROM}" \
                -Oz \
                -o "${flt_vcf}" \
                -s LOWQUAL \
                -i 'QUAL>=20' \
                "${ann_bcf}"

            bcftools index -t "${flt_vcf}"
            validate_vcfgz "${flt_vcf}"
            touch "${flt_done}"

            step_end=$(date +%s)
            append_status "filter" "done"
            append_timing "filter" "${step_start}" "${step_end}" "done"
            echo "[DONE] filter :: ${region}"
            echo "[TIME] filter: $(format_duration "$((step_end - step_start))")"
        fi

        region_end_ts=$(date +%s)
        append_timing "total_region_runtime" "${region_start_ts}" "${region_end_ts}" "done"

        echo "[DONE] Region ${region}"
        echo "[TIME] Total region runtime: $(format_duration "$((region_end_ts - region_start_ts))")"
        echo "End: $(date)"
        echo
    } >> "${logfile}" 2>&1
}

run_chromosome() {
    local chr="$1"
    local chrom_start_ts
    chrom_start_ts=$(date +%s)

    echo ">>> Starting chromosome ${chr}"

    awk -F'\t' -v target="${chr}" 'NR > 1 && $1 == target {print $1 "\t" $2 "\t" $3 "\t" $4}' "${REGIONS_FILE}" \
    | while IFS=$'\t' read -r c s e r; do
        run_chunk "${c}" "${s}" "${e}" "${r}"
    done

    local chrom_end_ts
    chrom_end_ts=$(date +%s)

    echo "<<< Finished chromosome ${chr} in $(format_duration "$((chrom_end_ts - chrom_start_ts))")"
}

MASTER_START_TS=$(date +%s)

echo "=================================================="
echo "PBS Job: ${PBS_JOBID:-NA}"
echo "Chromosome group: ${CHROM_GROUP_TAG}"
echo "Chromosomes: ${CHROM_LIST[*]}"
echo "Node CPUs: ${NODE_CPUS}"
echo "Parallel chromosomes: ${PARALLEL_CHROMS}"
echo "Threads per chromosome: ${THREADS_PER_CHROM}"
echo "AUTO_CLEAN_INCOMPLETE: ${AUTO_CLEAN_INCOMPLETE}"
echo "Start: $(date)"
echo "=================================================="

declare -a PIDS=()

wait_for_slot() {
    while [ "${#PIDS[@]}" -ge "${PARALLEL_CHROMS}" ]; do
        local new_pids=()
        local pid
        for pid in "${PIDS[@]}"; do
            if kill -0 "${pid}" 2>/dev/null; then
                new_pids+=("${pid}")
            fi
        done
        PIDS=("${new_pids[@]}")
        sleep 5
    done
}

for CHR in "${CHROM_LIST[@]}"; do
    wait_for_slot
    run_chromosome "${CHR}" &
    PIDS+=("$!")
done

for pid in "${PIDS[@]}"; do
    wait "${pid}"
done

MASTER_END_TS=$(date +%s)

echo "=================================================="
echo "[DONE] Chromosome group ${CHROM_GROUP_TAG}"
echo "[TIME] Total group runtime: $(format_duration "$((MASTER_END_TS - MASTER_START_TS))")"
echo "End: $(date)"
echo "==================================================