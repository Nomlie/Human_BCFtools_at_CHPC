#!/bin/bash
set -euo pipefail

show_help() {
    cat << 'EOF'
run_pipeline.sh
---------------


This mini pipeline is a small Bash + PBS workflow for human autosomal variant
calling with bcftools, designed so each person handles only their assigned
chromosome group and writes everything into a clearly named group-specific
folder such as chr1_chr14_chr21. The pipeline first creates a regions.tsv file
for that chromosome set by splitting each assigned chromosome into fixed-size
chunks, then submits an smp PBS array job that runs up to 3 chromosomes in parallel, each chromosome handled by one worker using 8 threads.
index using bcftools mpileup, call, annotate, and filter, while writing
checkpoints, logs, timing summaries, and status files for every region. Once
all chunk jobs for that person are done, the merge step concatenates the
chunk-level filtered VCFs back into one final VCF per chromosome inside that
person’s output folder.

What each person should do:
- Edit pipeline.conf once and set only their assigned CHROM_LIST
- Run:
    bash run_pipeline.sh setup

Wait for it to finish, this will produce regions.tsv file that the next script uses to divide the chromosomes into chunks
    bash run_pipeline.sh submit
- When chunk processing is complete, run:
    bash run_pipeline.sh merge

How to check progress:
- Run:
    bash run_pipeline.sh progress
- This submits the progress job and writes a report showing, for each
  chromosome in that group, how many chunks are finished, failed, or running

How to resume if smp stops midway:
- Do not start from scratch
- Resubmit the chunk job with:
    bash run_pipeline.sh submit
- Completed chunks are skipped automatically because the pipeline checks for
  valid outputs and .done checkpoint files

How to handle incomplete chunks:
- If a chunk failed midway and left partial files behind, check the log for the
  exact region and step that failed
- You can delete only the listed partial files for that region and resubmit
- Or submit the chunk job with automatic cleanup enabled:
    qsub -v AUTO_CLEAN_INCOMPLETE=1 02_run_chunks.sh

Commands:
  setup     Submit the region-generation job
  submit    Runs up to 3 chromosomes in parallel, each chromosome handled by one worker using 8 threads.
  progress  Submit the progress-report job
  merge     Submit the chromosome merge job
  help      Show this help message

All person-specific settings are read from pipeline.conf.
EOF
}

PIPELINE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SETUP_PBS="${PIPELINE_DIR}/01_make_regions.sh"
CHUNKS_PBS="${PIPELINE_DIR}/02_run_chunks.sh"
MERGE_PBS="${PIPELINE_DIR}/03_merge_chromosomes.sh"
PROGRESS_PBS="${PIPELINE_DIR}/04_check_progress.sh"

command="${1:-help}"

case "${command}" in
    setup) qsub "${SETUP_PBS}" ;;
    submit) qsub "${CHUNKS_PBS}" ;;
    progress) qsub "${PROGRESS_PBS}" ;;
    merge) qsub "${MERGE_PBS}" ;;
    help|-h|--help) show_help ;;
    *)
        echo "ERROR: Unknown command: ${command}"
        echo
        show_help
        exit 1
        ;;
esac