#!/usr/bin/env bash
set -euo pipefail

########################################
# Usage & argument parsing
########################################

usage() {
  cat >&2 <<EOF
Usage: $0 -r mt_reference.fa -o outdir [-t threads] bam1.bam [bam2.bam ...]

Description:
  For each whole-genome BAM, this script:
    1) Converts BAM to FASTQ (paired-end assumed)
    2) Realigns reads to the provided mitochondrial reference (FASTA)
    3) Calls variants (haploid) and builds consensus mtDNA FASTA
    4) Outputs per-sample mt BAM and FASTA
    5) Merges all per-sample mt FASTA into all_samples_mt.fa

  Temporary VCF files are created per sample and deleted automatically.

Outputs (per sample, in outdir):
  SAMPLE.mt.bam     - reads realigned to mt reference (sorted, indexed)
  SAMPLE_mt.fa      - consensus mtDNA sequence
Additionally:
  all_samples_mt.fa - merged FASTA of all samples

Options:
  -r  Path to mitochondrial reference FASTA (e.g. mt.fa)
  -o  Output directory (will be created if not exists)
  -t  Number of threads [default: 50]
  -h  Show this help and exit

Example:
  $0 -r mt.fa -o mt_result /path/to/*.bam
EOF
  exit 1
}

MT_REF=""
OUT_DIR=""
THREADS=50   # 默认 50 线程

while getopts "r:o:t:h" opt; do
  case "$opt" in
    r) MT_REF="$OPTARG" ;;
    o) OUT_DIR="$OPTARG" ;;
    t) THREADS="$OPTARG" ;;
    h) usage ;;
    *) usage ;;
  esac
done

shift $((OPTIND - 1))

# Remaining args are BAM files
if [ -z "${MT_REF}" ] || [ -z "${OUT_DIR}" ] || [ "$#" -lt 1 ]; then
  usage
fi

BAMS=("$@")


########################################
# Basic checks
########################################

if [ ! -f "$MT_REF" ]; then
  echo "[ERROR] mt reference FASTA not found: $MT_REF" >&2
  exit 1
fi

for tool in samtools bwa bcftools; do
  if ! command -v "$tool" >/dev/null 2>&1; then
    echo "[ERROR] '$tool' not found in PATH. Please install or load the environment." >&2
    exit 1
  fi
done

mkdir -p "$OUT_DIR"

MERGED_FASTA="${OUT_DIR}/all_samples_mt.fa"
> "$MERGED_FASTA"   # truncate / create

echo "[INFO] mt reference: $MT_REF"
echo "[INFO] Output directory: $OUT_DIR"
echo "[INFO] Threads: $THREADS"
echo "[INFO] Merged FASTA will be: $MERGED_FASTA"


########################################
# Index mt reference for samtools & bwa
########################################

# 1) samtools faidx
if [ ! -f "${MT_REF}.fai" ]; then
  echo "[INFO] Indexing mt reference for samtools (faidx)..."
  samtools faidx "$MT_REF"
fi

# 2) bwa index
if [ ! -f "${MT_REF}.bwt" ]; then
  echo "[INFO] Indexing mt reference for bwa..."
  bwa index "$MT_REF"
fi


########################################
# Process each BAM
########################################

for BAM in "${BAMS[@]}"; do
  if [ ! -f "$BAM" ]; then
    echo "[WARN] BAM not found, skip: $BAM" >&2
    continue
  fi

  BASENAME=$(basename "$BAM")
  # 样本名：取第一个点前面的部分，如 1098.sort.fixmate.rmdup.bam -> 1098
  SAMPLE=${BASENAME%%.*}

  echo "==============================================="
  echo "[INFO] Processing sample: ${SAMPLE}"
  echo "[INFO]   Input BAM: $BAM"

  # 输出前缀
  SAMPLE_PREFIX="${OUT_DIR}/${SAMPLE}"

  # 中间 FASTQ
  FQ1="${SAMPLE_PREFIX}.R1.fq"
  FQ2="${SAMPLE_PREFIX}.R2.fq"

  # 输出 BAM & FASTA
  MT_BAM="${SAMPLE_PREFIX}.mt.bam"
  FASTA="${SAMPLE_PREFIX}_mt.fa"

  ######################################
  # 1) BAM -> FASTQ (paired-end)
  ######################################
  echo "[INFO]   Converting BAM to FASTQ (paired-end assumed)..."
  samtools fastq -@ "$THREADS" \
    -1 "$FQ1" \
    -2 "$FQ2" \
    -0 /dev/null \
    -s /dev/null \
    -n "$BAM"

  ######################################
  # 2) Align FASTQ to mt reference
  ######################################
  echo "[INFO]   Aligning reads to mt reference with bwa mem..."
  bwa mem -t "$THREADS" "$MT_REF" "$FQ1" "$FQ2" \
    | samtools sort -@ "$THREADS" -o "$MT_BAM"

  samtools index "$MT_BAM"

  ######################################
  # 3) Call variants -> temporary VCF.gz
  ######################################
  VCF_GZ="${SAMPLE_PREFIX}.mt.tmp.vcf.gz"

  echo "[INFO]   Calling mt variants (haploid) to temporary VCF..."
  bcftools mpileup -Ou -f "$MT_REF" "$MT_BAM" \
    | bcftools call -mv --ploidy 1 -Oz -o "$VCF_GZ"

  bcftools index "$VCF_GZ"

  ######################################
  # 4) Build consensus mt sequence
  ######################################
  echo "[INFO]   Building consensus mtDNA from temporary VCF..."
  bcftools consensus -f "$MT_REF" "$VCF_GZ" > "$FASTA"

  # 把 FASTA header 改成 >SAMPLE
  TMPFA="${FASTA}.tmp"
  {
    read -r first_line
    echo ">${SAMPLE}"
    cat
  } < "$FASTA" > "$TMPFA"
  mv "$TMPFA" "$FASTA"

  # 追加到合并 FASTA
  cat "$FASTA" >> "$MERGED_FASTA"

  ######################################
  # 5) 清理中间文件（只留 bam 和 fasta）
  ######################################
  echo "[INFO]   Cleaning intermediate FASTQ and VCF..."
  rm -f "$FQ1" "$FQ2" "$VCF_GZ" "${VCF_GZ}.tbi" "${VCF_GZ}.csi" 2>/dev/null || true

  echo "[INFO]   Kept mt BAM:    $MT_BAM (and ${MT_BAM}.bai)"
  echo "[INFO]   Kept mt FASTA:  $FASTA"
  echo "[INFO]   Done for ${SAMPLE}"
done

echo "==============================================="
echo "[INFO] All samples processed."
echo "[INFO] Per-sample mt BAM:   ${OUT_DIR}/*.mt.bam"
echo "[INFO] Per-sample mt FASTA: ${OUT_DIR}/*_mt.fa"
echo "[INFO] Merged mt FASTA:     ${MERGED_FASTA}"

