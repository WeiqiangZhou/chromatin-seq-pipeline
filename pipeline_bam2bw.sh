#!/bin/bash

#SBATCH -p shared           # partition name
#SBATCH -t 4:00:00          # run time limit
#SBATCH -c 4                # number of CPU cores
#SBATCH --mem 10G           # memory
#SBATCH -o /dev/null        # stdout handled by submission wrapper (--output flag)
#SBATCH -e /dev/null        # stderr handled by submission wrapper (--error flag)

# Per-sample bigWig generation worker — called by pipeline_submit.sh via sbatch.
# Runs after pipeline_stats.sh completes (SLURM afterok dependency).
# Do NOT submit this script directly.
#
# Arguments:
#   $1            : pipeline root output folder
#   $2            : sample name (must match the BAM filename in <output>/mouse/)
#   --conf <path> : full path to pipeline.conf (required)
#   --spike       : generate spike-in-normalised bigWig using SCALE_CONST
#                   from <output>/spike_scale.env written by pipeline_stats.sh
#
# Without --spike : CPM-normalised  →  <bigwig>/<sample>.libnorm.bw
# With    --spike : spike-normalised →  <bigwig>/<sample>.spike_norm.bw

OUT_DIR=$1
SAMPLE=$2
SPIKE=false
CONF=""

for i in "$@"; do
    [[ "$i"    == "--spike" ]] && SPIKE=true
    [[ "$prev" == "--conf" ]]  && CONF=$i
    prev=$i
done

if [[ -z "$OUT_DIR" || -z "$SAMPLE" ]]; then
    echo "ERROR: Usage: sbatch pipeline_bam2bw.sh <output_folder> <sample> --conf <path> [--spike]"
    exit 1
fi
if [[ -z "$CONF" || ! -f "$CONF" ]]; then
    echo "ERROR: --conf <path> is required and the file must exist. Got: '$CONF'"
    exit 1
fi
# shellcheck source=pipeline.conf
source "$CONF"

MOUSE_BAM="${OUT_DIR}/mouse/${SAMPLE}.bam"
SPIKE_BAM="${OUT_DIR}/spike/${SAMPLE}.bam"
BW_DIR="${OUT_DIR}/bigwig"
SPIKE_SCALE_ENV="${OUT_DIR}/spike_scale.env"

if [[ ! -f "$MOUSE_BAM" ]]; then
    echo "ERROR: Mouse BAM not found: $MOUSE_BAM"
    exit 1
fi

mkdir -p "$BW_DIR"

[[ -n "$MODULE_SAMTOOLS" ]]  && module load "$MODULE_SAMTOOLS"
[[ -n "$MODULE_DEEPTOOLS" ]] && module load "$MODULE_DEEPTOOLS"

echo "=========================================="
echo "BigWig generation"
echo "  Sample     : $SAMPLE"
echo "  Mode       : $( [[ "$SPIKE" == true ]] && echo "spike-in normalised" || echo "CPM (library-size) normalised" )"
echo "  Start time : $(date)"
echo "=========================================="

# Coordinate-sort and index the name-sorted BAM for bamCoverage
SORT_BAM="${MOUSE_BAM%.bam}.coord.bam"
echo "  Coordinate-sorting..."
samtools sort -@ 4 -o "$SORT_BAM" "$MOUSE_BAM"
samtools index -@ 4 "$SORT_BAM"

# ── CPM normalisation ──────────────────────────────────────────────────────────
if [[ "$SPIKE" == false ]]; then
    BW_OUT="${BW_DIR}/${SAMPLE}.libnorm.bw"
    echo "  Running bamCoverage (CPM)..."
    bamCoverage -p 4 -b "$SORT_BAM" -o "$BW_OUT" --normalizeUsing CPM
    echo "  Output: $BW_OUT"

# ── Spike-in normalisation ─────────────────────────────────────────────────────
else
    if [[ ! -f "$SPIKE_SCALE_ENV" ]]; then
        echo "ERROR: spike_scale.env not found: $SPIKE_SCALE_ENV"
        echo "       pipeline_stats.sh must complete successfully before this job runs."
        rm -f "$SORT_BAM" "${SORT_BAM}.bai"
        exit 1
    fi
    # shellcheck source=/dev/null
    source "$SPIKE_SCALE_ENV"   # provides SCALE_CONST

    if [[ ! -f "$SPIKE_BAM" ]]; then
        echo "ERROR: Spike BAM not found: $SPIKE_BAM"
        rm -f "$SORT_BAM" "${SORT_BAM}.bai"
        exit 1
    fi

    SPIKE_READS=$(samtools view -c -F 260 "$SPIKE_BAM" 2>/dev/null)
    if [[ ! "$SPIKE_READS" =~ ^[0-9]+$ || "$SPIKE_READS" -eq 0 ]]; then
        echo "ERROR: Zero or invalid spike-in read count for $SAMPLE ($SPIKE_READS)"
        rm -f "$SORT_BAM" "${SORT_BAM}.bai"
        exit 1
    fi

    SCALE_FACTOR=$(awk "BEGIN {printf \"%.6f\", $SCALE_CONST / $SPIKE_READS}")
    BW_OUT="${BW_DIR}/${SAMPLE}.spike_norm.bw"

    echo "  SCALE_CONST  : $SCALE_CONST"
    echo "  Spike reads  : $SPIKE_READS"
    echo "  Scale factor : $SCALE_FACTOR"
    echo "  Running bamCoverage (spike-normalised)..."
    bamCoverage -p 4 -b "$SORT_BAM" -o "$BW_OUT" --scaleFactor "$SCALE_FACTOR"
    echo "  Output: $BW_OUT"
fi

rm -f "$SORT_BAM" "${SORT_BAM}.bai"

echo ""
echo "=========================================="
echo "Done. End time: $(date)"
echo "=========================================="
