#!/bin/bash

#SBATCH -p shared           # partition name
#SBATCH -t 2-0:00           # run time limit (days-hours:minutes)
#SBATCH -c 8                # number of CPU cores
#SBATCH --mem 20G           # memory
#SBATCH -o /dev/null        # stdout handled by submission wrapper (--output flag)
#SBATCH -e /dev/null        # stderr handled by submission wrapper (--error flag)

# Per-sample alignment worker — called by pipeline_submit.sh via sbatch.
# Do NOT submit this script directly.
#
# Arguments:
#   $1        : path to R1 input file (*_1.fq.gz)
#   $2        : path to R2 input file (*_2.fq.gz)
#   $3        : output prefix  (e.g. /output_folder/SampleA)
#   --spike   : (optional) also align to spike-in genome (Drosophila BDGP6)
#               and produce a deduplicated spike-in BAM for downstream
#               spike-in normalisation by pipeline_bam2bw.sh.
#
# Output (relative to output prefix):
#   <out_dir>/mouse/<sample>.bam         name-sorted, deduped main-genome BAM
#   <out_dir>/spike/<sample>.bam         name-sorted, deduped spike-in BAM [--spike only]
#   <out_dir>/logs/<sample>_*.log        tool log files

FQ1=$1
FQ2=$2
OUT=$3
SPIKE=false
CONF=""

for i in "$@"; do
    [[ "$i" == "--spike" ]]  && SPIKE=true
    [[ "$prev" == "--conf" ]] && CONF=$i
    prev=$i
done

if [[ -z "$FQ1" || -z "$FQ2" || -z "$OUT" ]]; then
    echo "ERROR: Usage: sbatch pipeline_align.sh <FQ1> <FQ2> <output_prefix> [--conf <path>] [--spike]"
    exit 1
fi

if [[ -z "$CONF" ]]; then
    echo "ERROR: --conf <path> is required (pass the full path to pipeline.conf)"
    exit 1
fi
if [[ ! -f "$CONF" ]]; then
    echo "ERROR: Config file not found: $CONF"
    exit 1
fi
# shellcheck source=pipeline.conf
source "$CONF"
# Resolve tilde in paths that may come from pipeline.conf
TRIMMOMATIC_JAR="${TRIMMOMATIC_JAR/#\~/$HOME}"
TRIMMOMATIC_ADAPTERS="${TRIMMOMATIC_ADAPTERS/#\~/$HOME}"
PICARD_JAR="${PICARD_JAR/#\~/$HOME}"

SAMPLE=$(basename "$OUT")
OUT_DIR=$(dirname "$OUT")
LOG_DIR="${OUT_DIR}/logs"
OUT_MOUSE="${OUT_DIR}/mouse/${SAMPLE}"
mkdir -p "$LOG_DIR" "${OUT_DIR}/mouse"

if [[ "$SPIKE" == true ]]; then
    OUT_SPIKE="${OUT_DIR}/spike/${SAMPLE}"
    mkdir -p "${OUT_DIR}/spike"
fi

echo "=========================================="
echo "Processing sample : $SAMPLE"
echo "  R1             : $FQ1"
echo "  R2             : $FQ2"
echo "  Output prefix  : $OUT"
echo "  Spike-in       : $SPIKE"
echo "  Start time     : $(date)"
echo "=========================================="

[[ -n "$MODULE_BOWTIE2" ]]  && module load "$MODULE_BOWTIE2"
[[ -n "$MODULE_SAMTOOLS" ]] && module load "$MODULE_SAMTOOLS"

# ── Trim adapters ──────────────────────────────────────────────────────────────
java -Xmx10g -jar "$TRIMMOMATIC_JAR" PE \
    "$FQ1" "$FQ2" \
    "${OUT}.R1.paired.gz"   "${OUT}.R1.unpaired.gz" \
    "${OUT}.R2.paired.gz"   "${OUT}.R2.unpaired.gz" \
    ILLUMINACLIP:"${TRIMMOMATIC_ADAPTERS}":2:30:10:2:True \
    LEADING:3 TRAILING:3 MINLEN:36
rm "${OUT}.R1.unpaired.gz"
rm "${OUT}.R2.unpaired.gz"

# ── Align to main genome ──────────────────────────────────────────────────────
bowtie2 -x "$BOWTIE2_INDEX_MAIN" \
    -p 8 -I 10 -X 700 -q \
    --local --very-sensitive-local --no-unal --no-mixed --no-discordant \
    -1 "${OUT}.R1.paired.gz" -2 "${OUT}.R2.paired.gz" \
    -S "${OUT_MOUSE}.sam" 2>"${LOG_DIR}/${SAMPLE}_bowtie2.log"
samtools view -@ 8 -bS -o "${OUT_MOUSE}.bam" "${OUT_MOUSE}.sam"
rm "${OUT_MOUSE}.sam"

# Filter by mapping quality and remove unwanted chromosomes
samtools view -b -q 10 -o "${OUT_MOUSE}.filtered.bam" "${OUT_MOUSE}.bam"
samtools view -h "${OUT_MOUSE}.filtered.bam" \
    | sed '/chrM/d;/random/d;/chrUn/d' \
    | samtools view -Sb - > "${OUT_MOUSE}.filtered.clean.bam"
rm "${OUT_MOUSE}.filtered.bam"
rm "${OUT_MOUSE}.bam"

# Remove PCR duplicates; keep name-sorted BAM as final output
samtools sort -o "${OUT_MOUSE}.sort.bam" "${OUT_MOUSE}.filtered.clean.bam"
java -Xmx16g -jar "$PICARD_JAR" MarkDuplicates \
    INPUT="${OUT_MOUSE}.sort.bam" \
    OUTPUT="${OUT_MOUSE}.filtered.clean.dup.bam" \
    METRICS_FILE="${LOG_DIR}/${SAMPLE}_picard.log" \
    REMOVE_DUPLICATES=true
samtools sort -n "${OUT_MOUSE}.filtered.clean.dup.bam" \
    -o "${OUT_MOUSE}.filtered.clean.dup.namesort.bam"
rm "${OUT_MOUSE}.sort.bam"
rm "${OUT_MOUSE}.filtered.clean.bam"
rm "${OUT_MOUSE}.filtered.clean.dup.bam"
mv "${OUT_MOUSE}.filtered.clean.dup.namesort.bam" "${OUT_MOUSE}.bam"

# ── Align to spike-in genome ──────────────────────────────────────────────────
if [[ "$SPIKE" == true ]]; then
    bowtie2 -x "$BOWTIE2_INDEX_SPIKE" \
        -p 8 -I 10 -X 700 -q \
        --local --very-sensitive-local --no-unal --no-mixed --no-discordant \
        -1 "${OUT}.R1.paired.gz" -2 "${OUT}.R2.paired.gz" \
        -S "${OUT_SPIKE}.sam" 2>"${LOG_DIR}/${SAMPLE}_spike_bowtie2.log"
    samtools view -@ 8 -bS -o "${OUT_SPIKE}.bam" "${OUT_SPIKE}.sam"
    rm "${OUT_SPIKE}.sam"

    # Remove PCR duplicates; keep name-sorted BAM as final output
    samtools sort -o "${OUT_SPIKE}.sort.bam" "${OUT_SPIKE}.bam"
    java -Xmx16g -jar "$PICARD_JAR" MarkDuplicates \
        INPUT="${OUT_SPIKE}.sort.bam" \
        OUTPUT="${OUT_SPIKE}.dup.bam" \
        METRICS_FILE="${LOG_DIR}/${SAMPLE}_spike_picard.log" \
        REMOVE_DUPLICATES=true
    samtools sort -n "${OUT_SPIKE}.dup.bam" -o "${OUT_SPIKE}.dup.namesort.bam"
    rm "${OUT_SPIKE}.bam"
    rm "${OUT_SPIKE}.sort.bam"
    rm "${OUT_SPIKE}.dup.bam"
    mv "${OUT_SPIKE}.dup.namesort.bam" "${OUT_SPIKE}.bam"
fi

# Clean up trimmed reads after all alignments complete
rm "${OUT}.R1.paired.gz"
rm "${OUT}.R2.paired.gz"

echo "=========================================="
echo "Finished sample : $SAMPLE"
echo "  End time      : $(date)"
echo "=========================================="
