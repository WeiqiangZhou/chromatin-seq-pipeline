#!/bin/bash

#SBATCH -p shared           # partition name
#SBATCH -t 4:00:00          # run time limit
#SBATCH -c 4                # number of CPU cores
#SBATCH --mem 8G            # memory
#SBATCH -o /dev/null        # stdout handled by submission wrapper (--output flag)
#SBATCH -e /dev/null        # stderr handled by submission wrapper (--error flag)

# Alignment statistics and spike-in normalisation prep job.
# Runs after ALL per-sample alignment jobs complete (SLURM afterok dependency).
# Must finish before any per-sample bam2bw job starts.
#
# Arguments:
#   $1              : pipeline root output folder
#   --conf <path>   : full path to pipeline.conf (required)
#   --spike         : also count spike-in reads and write spike_scale.env
#   --use-min-spike : (with --spike) set SCALE_CONST to the minimum spike-in
#                     read count across all samples rather than 1,000,000
#
# Outputs:
#   <output_folder>/alignment_stats.tsv   — per-sample metrics table
#   <output_folder>/spike_scale.env       — SCALE_CONST for bam2bw jobs [--spike only]

OUT_DIR=$1
SPIKE=false
USE_MIN_SPIKE=false
CONF=""

for i in "$@"; do
    [[ "$i"    == "--spike" ]]         && SPIKE=true
    [[ "$i"    == "--use-min-spike" ]] && USE_MIN_SPIKE=true
    [[ "$prev" == "--conf" ]]          && CONF=$i
    prev=$i
done

if [[ -z "$OUT_DIR" ]]; then
    echo "ERROR: Usage: sbatch pipeline_stats.sh <output_folder> --conf <path> [--spike] [--use-min-spike]"
    exit 1
fi
if [[ -z "$CONF" || ! -f "$CONF" ]]; then
    echo "ERROR: --conf <path> is required and the file must exist. Got: '$CONF'"
    exit 1
fi
# shellcheck source=pipeline.conf
source "$CONF"

LOG_DIR="${OUT_DIR}/logs"
MOUSE_DIR="${OUT_DIR}/mouse"
SPIKE_DIR="${OUT_DIR}/spike"
STATS_TSV="${OUT_DIR}/alignment_stats.tsv"
SPIKE_SCALE_ENV="${OUT_DIR}/spike_scale.env"

[[ -n "$MODULE_SAMTOOLS" ]] && module load "$MODULE_SAMTOOLS"

echo "=========================================="
echo "Alignment statistics"
echo "  Output folder : $OUT_DIR"
echo "  Spike mode    : $SPIKE"
echo "  Start time    : $(date)"
echo "=========================================="

# ── Collect per-sample stats into arrays (single pass) ────────────────────────
SAMPLE_NAMES=()
INPUT_READS_ARR=()
ALIGNED_READS_ARR=()
ALIGN_RATE_ARR=()
FINAL_READS_ARR=()
DUP_RATE_ARR=()
SPIKE_READS_ARR=()   # populated only in spike mode

for MOUSE_BAM in "$MOUSE_DIR"/*.bam; do
    [[ "$MOUSE_BAM" == *.coord.bam ]] && continue
    BASENAME=$(basename "$MOUSE_BAM" .bam)
    echo "  Processing: $BASENAME"

    # ── bowtie2 log ────────────────────────────────────────────────────────────
    BT2_LOG="${LOG_DIR}/${BASENAME}_bowtie2.log"
    INPUT_READS="NA"; ALIGNED_READS="NA"; ALIGN_RATE="NA"
    if [[ -f "$BT2_LOG" ]]; then
        INPUT_READS=$(grep -m1 "reads; of these:" "$BT2_LOG" | awk '{print $1}')
        ALIGNED_READS=$(grep -E "aligned concordantly exactly 1 time|aligned concordantly >1 times" \
            "$BT2_LOG" | awk '{s += $1} END {print s+0}')
        ALIGN_RATE=$(grep "overall alignment rate" "$BT2_LOG" | awk '{print $1}')
    else
        echo "    WARNING: bowtie2 log not found: $BT2_LOG"
    fi

    # ── Picard MarkDuplicates log ──────────────────────────────────────────────
    PICARD_LOG="${LOG_DIR}/${BASENAME}_picard.log"
    DUP_RATE="NA"
    if [[ -f "$PICARD_LOG" ]]; then
        DUP_RATE=$(awk '
            /^## METRICS CLASS/ { in_m=1; next }
            in_m && /^LIBRARY/  { hdr=$0; next }
            in_m && hdr && NF>0 {
                n=split(hdr,h,"\t"); split($0,d,"\t")
                for(i=1;i<=n;i++)
                    if(h[i]=="PERCENT_DUPLICATION"){ print d[i]; exit }
            }
            /^## HISTOGRAM/ { in_m=0 }
        ' "$PICARD_LOG")
    else
        echo "    WARNING: Picard log not found: $PICARD_LOG"
    fi

    # ── Final read count ───────────────────────────────────────────────────────
    # -F 2308 excludes unmapped (4) + secondary (256) + supplementary (2048).
    # Divide by 2: each fragment = two reads in a paired-end BAM.
    RAW_COUNT=$(samtools view -c -F 2308 "$MOUSE_BAM" 2>/dev/null)
    FINAL_READS=$(( RAW_COUNT / 2 ))

    # ── Spike-in reads ─────────────────────────────────────────────────────────
    SPIKE_READS="NA"
    if [[ "$SPIKE" == true ]]; then
        SPIKE_BAM="${SPIKE_DIR}/${BASENAME}.bam"
        if [[ -f "$SPIKE_BAM" ]]; then
            SPIKE_READS=$(samtools view -c -F 260 "$SPIKE_BAM" 2>/dev/null)
        else
            echo "    WARNING: Spike BAM not found: $SPIKE_BAM"
        fi
    fi

    SAMPLE_NAMES+=("$BASENAME")
    INPUT_READS_ARR+=("$INPUT_READS")
    ALIGNED_READS_ARR+=("$ALIGNED_READS")
    ALIGN_RATE_ARR+=("$ALIGN_RATE")
    FINAL_READS_ARR+=("$FINAL_READS")
    DUP_RATE_ARR+=("$DUP_RATE")
    SPIKE_READS_ARR+=("$SPIKE_READS")
done

if [[ ${#SAMPLE_NAMES[@]} -eq 0 ]]; then
    echo "ERROR: No BAM files found in $MOUSE_DIR"
    exit 1
fi

# ── Determine SCALE_CONST (spike mode only) ────────────────────────────────────
SCALE_CONST=1000000
if [[ "$SPIKE" == true ]]; then
    if [[ "$USE_MIN_SPIKE" == true ]]; then
        MIN_READS=""
        for i in "${!SAMPLE_NAMES[@]}"; do
            SR="${SPIKE_READS_ARR[$i]}"
            [[ "$SR" =~ ^[0-9]+$ && "$SR" -gt 0 ]] || continue
            if [[ -z "$MIN_READS" || "$SR" -lt "$MIN_READS" ]]; then
                MIN_READS="$SR"
                MIN_SAMPLE="${SAMPLE_NAMES[$i]}"
            fi
        done
        if [[ -n "$MIN_READS" ]]; then
            SCALE_CONST="$MIN_READS"
            echo ""
            echo "--use-min-spike: SCALE_CONST = $SCALE_CONST (sample: $MIN_SAMPLE)"
        else
            echo "WARNING: No valid spike reads found; falling back to SCALE_CONST=1000000"
        fi
    else
        echo ""
        echo "Fixed SCALE_CONST = $SCALE_CONST"
    fi

    # Write SCALE_CONST for downstream per-sample bam2bw jobs to source
    echo "SCALE_CONST=${SCALE_CONST}" > "$SPIKE_SCALE_ENV"
    echo "Wrote: $SPIKE_SCALE_ENV"
fi

# ── Write alignment_stats.tsv ──────────────────────────────────────────────────
echo ""
if [[ "$SPIKE" == true ]]; then
    printf "sample\tinput_reads\taligned_reads\talignment_rate\tfinal_reads\tdup_rate\tspike_reads\tscale_factor\n" \
        > "$STATS_TSV"
else
    printf "sample\tinput_reads\taligned_reads\talignment_rate\tfinal_reads\tdup_rate\n" \
        > "$STATS_TSV"
fi

for i in "${!SAMPLE_NAMES[@]}"; do
    BASENAME="${SAMPLE_NAMES[$i]}"
    INPUT_READS="${INPUT_READS_ARR[$i]}"
    ALIGNED_READS="${ALIGNED_READS_ARR[$i]}"
    ALIGN_RATE="${ALIGN_RATE_ARR[$i]}"
    FINAL_READS="${FINAL_READS_ARR[$i]}"
    DUP_RATE="${DUP_RATE_ARR[$i]}"

    if [[ "$SPIKE" == true ]]; then
        SR="${SPIKE_READS_ARR[$i]}"
        SCALE_FACTOR="NA"
        if [[ "$SR" =~ ^[0-9]+$ && "$SR" -gt 0 ]]; then
            SCALE_FACTOR=$(awk "BEGIN {printf \"%.6f\", $SCALE_CONST / $SR}")
        fi
        printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
            "$BASENAME" "$INPUT_READS" "$ALIGNED_READS" "$ALIGN_RATE" \
            "$FINAL_READS" "$DUP_RATE" "$SR" "$SCALE_FACTOR" >> "$STATS_TSV"
    else
        printf "%s\t%s\t%s\t%s\t%s\t%s\n" \
            "$BASENAME" "$INPUT_READS" "$ALIGNED_READS" "$ALIGN_RATE" \
            "$FINAL_READS" "$DUP_RATE" >> "$STATS_TSV"
    fi
done

echo "Statistics written to: $STATS_TSV"
echo ""
echo "=========================================="
echo "Stats complete. End time: $(date)"
echo "=========================================="
