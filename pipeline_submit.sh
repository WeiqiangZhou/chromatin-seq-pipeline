#!/bin/bash

# Main pipeline submission wrapper.
#
# Usage:
#   bash pipeline_submit.sh <input_folder> <output_folder> [--spike] [--use-min-spike]
#
#   <input_folder>   : top-level folder whose subdirectories each represent one sample.
#                      Each subdirectory must contain paired-end files *_1.fq.gz and *_2.fq.gz.
#   <output_folder>  : root folder for all output. Will be created if absent.
#   --spike          : enable spike-in alignment and generate spike-in-normalised bigWigs.
#   --use-min-spike  : (only with --spike) use the minimum spike-in read count across
#                      all samples as the normalisation constant rather than 1,000,000.
#
# Job dependency chain:
#   [align_S1] [align_S2] ... (one per sample, parallel)
#        └─────────┴── afterok ──► [stats]  (one job; writes alignment_stats.tsv
#                                            and spike_scale.env)
#                                      └── afterok ──► [bam2bw_S1] [bam2bw_S2] ...
#                                                       (one per sample, parallel)
#
# Tool paths, genome indices, and module names are read from pipeline.conf,
# which must live in the same directory as this script.
#
# Output layout:
#   <output_folder>/
#     mouse/               final (name-sorted, deduped) main-genome BAMs
#     spike/               final (name-sorted, deduped) spike-in BAMs  [--spike only]
#     bigwig/              bigWig files (*.libnorm.bw or *.spike_norm.bw)
#     logs/                per-job SLURM stdout/stderr and tool log files
#     alignment_stats.tsv  per-sample alignment metrics
#     spike_scale.env      SCALE_CONST used for bigWig normalisation  [--spike only]

INPUT_FOLDER=$1
OUTPUT_FOLDER=$2
SPIKE=false
USE_MIN_SPIKE=false

for arg in "$@"; do
    [[ "$arg" == "--spike" ]]         && SPIKE=true
    [[ "$arg" == "--use-min-spike" ]] && USE_MIN_SPIKE=true
done

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONF="${SCRIPT_DIR}/pipeline.conf"
ALIGN_WORKER="${SCRIPT_DIR}/pipeline_align.sh"
STATS_WORKER="${SCRIPT_DIR}/pipeline_stats.sh"
BW_WORKER="${SCRIPT_DIR}/pipeline_bam2bw.sh"
CHECK_SCRIPT="${SCRIPT_DIR}/check_tools.sh"

if [[ -z "$INPUT_FOLDER" || -z "$OUTPUT_FOLDER" ]]; then
    echo "ERROR: Usage: bash pipeline_submit.sh <input_folder> <output_folder> [--spike] [--use-min-spike]"
    exit 1
fi

if [[ ! -f "$CONF" ]]; then
    echo "ERROR: Config file not found: $CONF"
    echo "       Edit pipeline.conf in the same directory as this script."
    exit 1
fi
# shellcheck source=pipeline.conf
source "$CONF"

for script in "$ALIGN_WORKER" "$STATS_WORKER" "$BW_WORKER"; do
    if [[ ! -f "$script" ]]; then
        echo "ERROR: Required worker script not found: $script"
        exit 1
    fi
done

# Run tool / path check before submitting anything
if [[ -f "$CHECK_SCRIPT" ]]; then
    bash "$CHECK_SCRIPT" || exit 1
else
    echo "WARNING: check_tools.sh not found — skipping dependency check."
fi

if [[ "$USE_MIN_SPIKE" == true && "$SPIKE" == false ]]; then
    echo "WARNING: --use-min-spike has no effect without --spike; ignoring."
    USE_MIN_SPIKE=false
fi

# ── Validate file/index paths ─────────────────────────────────────────────────
_CONF_ERRORS=0
_check_file() {
    local label=$1 path="${2/#\~/$HOME}"
    [[ -f "$path" ]] || { echo "ERROR [pipeline.conf]: $label not found: $path"; (( _CONF_ERRORS++ )); }
}
_check_index() {
    local label=$1 prefix=$2
    [[ -f "${prefix}.1.bt2" || -f "${prefix}.1.bt2l" ]] || {
        echo "ERROR [pipeline.conf]: $label index not found at: $prefix"
        (( _CONF_ERRORS++ ))
    }
}
_check_file  "TRIMMOMATIC_JAR"      "${TRIMMOMATIC_JAR/#\~/$HOME}"
_check_file  "TRIMMOMATIC_ADAPTERS" "${TRIMMOMATIC_ADAPTERS/#\~/$HOME}"
_check_file  "PICARD_JAR"           "${PICARD_JAR/#\~/$HOME}"
_check_index "BOWTIE2_INDEX_MAIN"   "$BOWTIE2_INDEX_MAIN"
[[ "$SPIKE" == true ]] && _check_index "BOWTIE2_INDEX_SPIKE" "$BOWTIE2_INDEX_SPIKE"
if [[ $_CONF_ERRORS -gt 0 ]]; then
    echo "Aborting: fix the errors above in $CONF before resubmitting."
    exit 1
fi

mkdir -p "${OUTPUT_FOLDER}/logs"

echo "Pipeline configuration"
echo "  Config             : $CONF"
echo "  Input folder       : $INPUT_FOLDER"
echo "  Output folder      : $OUTPUT_FOLDER"
echo "  Spike-in alignment : $SPIKE"
echo "  Use min-spike norm : $USE_MIN_SPIKE"
echo "  Trimmomatic JAR    : ${TRIMMOMATIC_JAR/#\~/$HOME}"
echo "  Adapter file       : ${TRIMMOMATIC_ADAPTERS/#\~/$HOME}"
echo "  Picard JAR         : ${PICARD_JAR/#\~/$HOME}"
echo "  Main genome index  : $BOWTIE2_INDEX_MAIN"
[[ "$SPIKE" == true ]] && echo "  Spike genome index : $BOWTIE2_INDEX_SPIKE"
echo ""

# ── Stage 1: per-sample alignment jobs (parallel) ─────────────────────────────
ALIGN_JOB_IDS=()
SAMPLES=()
SUBMITTED=0
SKIPPED=0

for SUBDIR in "$INPUT_FOLDER"/*/; do
    [[ -d "$SUBDIR" ]] || continue
    SAMPLE=$(basename "$SUBDIR")

    FQ1=$(find "$SUBDIR" -maxdepth 1 -name "*_1.fq.gz" | head -1)
    FQ2=$(find "$SUBDIR" -maxdepth 1 -name "*_2.fq.gz" | head -1)

    if [[ -z "$FQ1" || -z "$FQ2" ]]; then
        echo "WARNING: _1.fq.gz / _2.fq.gz not found in $SUBDIR — skipping '$SAMPLE'."
        (( SKIPPED++ ))
        continue
    fi

    ALIGN_ARGS=(--conf "$CONF")
    [[ "$SPIKE" == true ]] && ALIGN_ARGS+=(--spike)

    JOB_ID=$(sbatch \
        --job-name="align_${SAMPLE}" \
        --output="${OUTPUT_FOLDER}/logs/align_${SAMPLE}_%j.out" \
        --error="${OUTPUT_FOLDER}/logs/align_${SAMPLE}_%j.err" \
        "$ALIGN_WORKER" "$FQ1" "$FQ2" "${OUTPUT_FOLDER}/${SAMPLE}" "${ALIGN_ARGS[@]}" \
        | awk '{print $NF}')

    echo "Submitted alignment job ${JOB_ID} : ${SAMPLE}"
    ALIGN_JOB_IDS+=("$JOB_ID")
    SAMPLES+=("$SAMPLE")
    (( SUBMITTED++ ))
done

echo ""
echo "Alignment jobs : ${SUBMITTED} submitted, ${SKIPPED} skipped"

if [[ ${#ALIGN_JOB_IDS[@]} -eq 0 ]]; then
    echo "ERROR: No alignment jobs submitted. Exiting."
    exit 1
fi

# ── Stage 2: stats job (after all alignments) ─────────────────────────────────
ALIGN_DEP="afterok:$(IFS=:; echo "${ALIGN_JOB_IDS[*]}")"

STATS_ARGS=("$OUTPUT_FOLDER" --conf "$CONF")
[[ "$SPIKE"         == true ]] && STATS_ARGS+=(--spike)
[[ "$USE_MIN_SPIKE" == true ]] && STATS_ARGS+=(--use-min-spike)

STATS_JOB_ID=$(sbatch \
    --job-name="stats" \
    --dependency="$ALIGN_DEP" \
    --output="${OUTPUT_FOLDER}/logs/stats_%j.out" \
    --error="${OUTPUT_FOLDER}/logs/stats_%j.err" \
    "$STATS_WORKER" "${STATS_ARGS[@]}" \
    | awk '{print $NF}')

echo "Submitted stats job   ${STATS_JOB_ID} (after: ${ALIGN_JOB_IDS[*]})"

# ── Stage 3: per-sample bam2bw jobs (after stats job, parallel) ───────────────
BW_DEP="afterok:${STATS_JOB_ID}"
BW_JOB_IDS=()

for SAMPLE in "${SAMPLES[@]}"; do
    BW_ARGS=("$OUTPUT_FOLDER" "$SAMPLE" --conf "$CONF")
    [[ "$SPIKE" == true ]] && BW_ARGS+=(--spike)

    BW_JOB_ID=$(sbatch \
        --job-name="bam2bw_${SAMPLE}" \
        --dependency="$BW_DEP" \
        --output="${OUTPUT_FOLDER}/logs/bam2bw_${SAMPLE}_%j.out" \
        --error="${OUTPUT_FOLDER}/logs/bam2bw_${SAMPLE}_%j.err" \
        "$BW_WORKER" "${BW_ARGS[@]}" \
        | awk '{print $NF}')

    echo "Submitted bam2bw job  ${BW_JOB_ID} : ${SAMPLE}"
    BW_JOB_IDS+=("$BW_JOB_ID")
done

echo ""
echo "All jobs submitted."
echo "  Stage 1 — alignment (${#ALIGN_JOB_IDS[@]} jobs) : ${ALIGN_JOB_IDS[*]}"
echo "  Stage 2 — stats     (1 job)                    : ${STATS_JOB_ID}"
echo "  Stage 3 — bam2bw    (${#BW_JOB_IDS[@]} jobs)  : ${BW_JOB_IDS[*]}"
echo ""
echo "  alignment_stats.tsv : ${OUTPUT_FOLDER}/alignment_stats.tsv"
echo "  BigWig output       : ${OUTPUT_FOLDER}/bigwig/"
echo "  Logs                : ${OUTPUT_FOLDER}/logs/"
