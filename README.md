# ChIP-seq / CUT&RUN Analysis Pipeline

A SLURM-based pipeline for paired-end ChIP-seq and CUT&RUN data. Handles adapter
trimming, alignment, deduplication, and bigWig generation. Optionally supports
spike-in normalization using a Drosophila (*BDGP6*) genome.

---

## Table of Contents

1. [Pipeline overview](#1-pipeline-overview)
2. [Prerequisites](#2-prerequisites)
3. [Quick start](#3-quick-start)
4. [Input data structure](#4-input-data-structure)
5. [Configuration — pipeline.conf](#5-configuration--pipelineconf)
6. [Running the pipeline](#6-running-the-pipeline)
7. [Job dependency chain](#7-job-dependency-chain)
8. [Processing steps](#8-processing-steps)
9. [Output files](#9-output-files)
10. [Spike-in normalization details](#10-spike-in-normalization-details)

---

## 1. Pipeline overview

| Script | Role |
|---|---|
| `pipeline_submit.sh` | Entry point — validates config, runs tool check, submits all SLURM jobs |
| `pipeline_align.sh` | Per-sample worker — trims, aligns, filters, deduplicates |
| `pipeline_stats.sh` | Aggregates alignment statistics; computes spike-in scale constants |
| `pipeline_bam2bw.sh` | Per-sample worker — generates one bigWig file per sample |
| `pipeline.conf` | User-editable configuration (tool paths, genome indices, modules) |
| `check_tools.sh` | Standalone tool dependency checker; also called automatically by `pipeline_submit.sh` |

---

## 2. Prerequisites

The following tools must be available either on `PATH` or via `module load`.
Run `bash check_tools.sh` to verify your environment before submitting jobs.

| Tool | Purpose | Tested version |
|---|---|---|
| **Java** (JDK 8+) | Required to run Trimmomatic and Picard | 8, 11, 17 |
| **Trimmomatic** | Adapter trimming | 0.39 |
| **bowtie2** | Read alignment | 2.x |
| **samtools** | BAM processing and read counting | 1.x |
| **Picard** | PCR duplicate removal | 2.x |
| **deeptools** (`bamCoverage`) | bigWig generation | 3.x |

You will also need:

- A **bowtie2 index** for your main genome (e.g. mm10, hg38).
- A **bowtie2 index** for the spike-in genome (e.g. Drosophila BDGP6) if using `--spike`.
- The **Trimmomatic adapter file** matching your library kit (e.g. `TruSeq3-PE-2.fa`).

---

## 3. Quick start

```bash
# 1. Edit pipeline.conf with your tool paths and genome indices
vim pipeline.conf

# 2. (Optional) Verify tools are accessible
bash check_tools.sh

# 3. Run the pipeline — no spike-in
bash pipeline_submit.sh /path/to/input /path/to/output

# 3. Run the pipeline — with spike-in normalization
bash pipeline_submit.sh /path/to/input /path/to/output --spike

# 3. Run with spike-in, using the sample with the fewest spike reads as the reference
bash pipeline_submit.sh /path/to/input /path/to/output --spike --use-min-spike
```

---

## 4. Input data structure

Organize your raw FASTQ files as one subdirectory per sample. Each subdirectory
must contain exactly one pair of gzipped FASTQ files named `*_1.fq.gz` and `*_2.fq.gz`.

```
input_folder/
├── SampleA/
│   ├── SampleA_1.fq.gz      # R1 reads
│   └── SampleA_2.fq.gz      # R2 reads
├── SampleB/
│   ├── SampleB_1.fq.gz
│   └── SampleB_2.fq.gz
└── SampleC/
    ├── SampleC_1.fq.gz
    └── SampleC_2.fq.gz
```

- The subdirectory name becomes the **sample name** used for all output files.
- Subdirectories without a matching `*_1.fq.gz` / `*_2.fq.gz` pair are skipped with a warning.

---

## 5. Configuration — pipeline.conf

Edit `pipeline.conf` before running. All pipeline scripts read this file automatically.

```bash
# ── Tool paths ────────────────────────────────────────────────────────────────

# Full path to the Trimmomatic JAR file
TRIMMOMATIC_JAR=$HOME/Trimmomatic-0.39/trimmomatic-0.39.jar

# Adapter FASTA used by Trimmomatic ILLUMINACLIP
# Choose the file that matches your library preparation kit:
#   TruSeq2-PE.fa       — older GA IIx sequencers
#   TruSeq3-PE-2.fa     — HiSeq / MiSeq (most common)
#   NexteraPE-PE.fa     — Nextera XT libraries
TRIMMOMATIC_ADAPTERS=$HOME/Trimmomatic-0.39/adapters/TruSeq3-PE-2.fa

# Full path to the Picard JAR file
PICARD_JAR=/path/to/picard.jar

# ── Genome index prefixes (bowtie2) ──────────────────────────────────────────
# Do NOT include the .bt2 / .bt2l extension — provide the prefix only.

# Main genome (the organism you are studying, e.g. mm10, hg38, dm6)
BOWTIE2_INDEX_MAIN=/path/to/bowtie2_index/mm10

# Spike-in genome (only used with --spike; default: Drosophila BDGP6)
BOWTIE2_INDEX_SPIKE=/path/to/bowtie2_index/BDGP6

# ── Environment modules ───────────────────────────────────────────────────────
# Set to the exact module name on your HPC, or leave empty ("") to skip
# module loading if the tool is already on PATH (e.g. managed by conda).

MODULE_BOWTIE2=bowtie2      # used during alignment
MODULE_SAMTOOLS=samtools    # used during alignment and bigWig generation
MODULE_DEEPTOOLS=deeptools  # used during bigWig generation (provides bamCoverage)
```

> **Tip — building a bowtie2 index:**
> ```bash
> bowtie2-build genome.fa /path/to/index_prefix
> ```

---

## 6. Running the pipeline

All jobs are submitted from a **login or interactive node**. Do not run
`pipeline_submit.sh` inside a batch job.

```bash
bash pipeline_submit.sh <input_folder> <output_folder> [OPTIONS]
```

| Option | Description |
|---|---|
| `--spike` | Align reads to the spike-in genome; generate spike-in-normalised bigWigs |
| `--use-min-spike` | *(requires `--spike`)* Set the normalisation constant to the minimum spike-in read count across all samples. The sample with the fewest spike reads gets scale factor = 1.0; all others scale down proportionally. Default: constant = 1,000,000 |

`pipeline_submit.sh` will:

1. Verify that `pipeline.conf` exists and all configured paths are valid.
2. Run `check_tools.sh` to confirm all tools are accessible.
3. Submit SLURM jobs and print the job IDs for all three stages.

---

## 7. Job dependency chain

```
Stage 1 — alignment (one job per sample, all run in parallel)
  align_SampleA ──┐
  align_SampleB ──┼── afterok ──► Stage 2
  align_SampleC ──┘

Stage 2 — statistics (one job, runs after all Stage 1 jobs succeed)
  stats ──────────────────────── afterok ──► Stage 3

Stage 3 — bigWig generation (one job per sample, all run in parallel)
  bam2bw_SampleA
  bam2bw_SampleB
  bam2bw_SampleC
```

- If **any** alignment job fails, the stats job and all bigWig jobs are
  automatically cancelled by SLURM (afterok semantics).
- Each stage's jobs can be monitored with `squeue -u $USER`.

---

## 8. Processing steps

### Per sample (Stage 1 — `pipeline_align.sh`)

1. **Adapter trimming** — Trimmomatic PE with ILLUMINACLIP, `LEADING:3 TRAILING:3 MINLEN:36`.
2. **Alignment to main genome** — bowtie2 with `--local --very-sensitive-local --no-unal --no-mixed --no-discordant`, insert size 10–700 bp.
3. **Mapping quality filter** — retain reads with MAPQ ≥ 10.
4. **Chromosome filter** — remove chrM, random, and unplaced (chrUn) contigs.
5. **PCR duplicate removal** — Picard MarkDuplicates (`REMOVE_DUPLICATES=true`).
6. *(--spike only)* **Spike-in alignment** — bowtie2 against the spike-in genome with the same settings as step 2.
7. *(--spike only)* **Spike-in deduplication** — Picard MarkDuplicates on the spike-in BAM.

Final BAMs are **name-sorted** (compatible with MACS2 peak calling).

### Aggregate (Stage 2 — `pipeline_stats.sh`)

- Parses per-sample bowtie2 and Picard log files.
- Counts final reads from each BAM with `samtools view -c`.
- *(--spike)* Counts mapped spike-in reads; determines the normalisation constant (`SCALE_CONST`).
- Writes `alignment_stats.tsv` and *(--spike)* `spike_scale.env`.

### Per sample (Stage 3 — `pipeline_bam2bw.sh`)

- Coordinate-sorts the name-sorted BAM for bamCoverage.
- **Without `--spike`**: `bamCoverage --normalizeUsing CPM`.
- **With `--spike`**: reads `SCALE_CONST` from `spike_scale.env`; computes `scale_factor = SCALE_CONST / spike_reads`; runs `bamCoverage --scaleFactor`.
- Removes the temporary coordinate-sorted BAM.

---

## 9. Output files

```
output_folder/
│
├── mouse/                         Main-genome BAMs (name-sorted, deduped)
│   ├── SampleA.bam
│   ├── SampleB.bam
│   └── SampleC.bam
│
├── spike/                         Spike-in BAMs [--spike only]
│   ├── SampleA.bam
│   ├── SampleB.bam
│   └── SampleC.bam
│
├── bigwig/                        Genome coverage tracks
│   ├── SampleA.libnorm.bw         CPM-normalised  [without --spike]
│   ├── SampleA.spike_norm.bw      Spike-normalised [with --spike]
│   ├── SampleB.libnorm.bw / .spike_norm.bw
│   └── SampleC.libnorm.bw / .spike_norm.bw
│
├── alignment_stats.tsv            Per-sample alignment metrics table (see below)
│
├── spike_scale.env                SCALE_CONST used for bigWig generation [--spike only]
│   (contents: SCALE_CONST=1000000)
│
└── logs/
    ├── align_SampleA_<jobid>.out/.err      Alignment job stdout/stderr
    ├── align_SampleB_<jobid>.out/.err
    ├── stats_<jobid>.out/.err              Stats job stdout/stderr
    ├── bam2bw_SampleA_<jobid>.out/.err     BigWig job stdout/stderr
    ├── SampleA_bowtie2.log                 bowtie2 alignment summary
    ├── SampleA_picard.log                  Picard duplication metrics
    ├── SampleA_spike_bowtie2.log           Spike-in alignment summary [--spike only]
    └── SampleA_spike_picard.log            Spike-in duplication metrics [--spike only]
```

### alignment_stats.tsv columns

| Column | Description |
|---|---|
| `sample` | Sample name |
| `input_reads` | Read pairs fed into bowtie2 (after Trimmomatic) |
| `aligned_reads` | Concordantly aligned read pairs (unique + multi-mapping) |
| `alignment_rate` | Overall alignment rate reported by bowtie2 (e.g. `83.59%`) |
| `final_reads` | Read pairs remaining after MAPQ filter, chromosome filter, and deduplication |
| `dup_rate` | Fraction of duplicate read pairs removed by Picard |
| `spike_reads` | Mapped read pairs in the spike-in BAM *(--spike only)* |
| `scale_factor` | `SCALE_CONST / spike_reads` used for bigWig generation *(--spike only)* |

---

## 10. Spike-in normalization details

Spike-in normalization corrects for differences in ChIP efficiency across samples
by scaling each sample's signal relative to the number of reads that mapped to the
spike-in genome.

```
scale_factor = SCALE_CONST / spike_reads_for_this_sample
```

**Default (`--spike` without `--use-min-spike`)**

`SCALE_CONST = 1,000,000`. All samples are expressed as signal per million
spike-in-mapped reads, analogous to RPM.

**Minimum-spike mode (`--spike --use-min-spike`)**

`SCALE_CONST = min(spike_reads across all samples)`. The sample with the fewest
spike-in reads receives `scale_factor = 1.0` and all other samples are scaled
down proportionally. This keeps signal values close to raw read depth and avoids
very large scale factors when spike read counts vary widely.

`SCALE_CONST` is determined once by `pipeline_stats.sh`, written to
`spike_scale.env`, and read by each per-sample `pipeline_bam2bw.sh` job.
