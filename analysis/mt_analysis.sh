touch extract_mt_from_bam.sh

#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat >&2 <<EOF
Usage:
  $0 -r mt_reference.fa -o outdir [-t threads] bam1.bam [bam2.bam ...]

Description:
  Extract mitochondrial consensus genome from WGS BAM files.
  For each BAM:
    1. Convert BAM to paired FASTQ
    2. Align reads to mitochondrial reference
    3. Keep high-quality mt alignments
    4. Call haploid variants
    5. Mask low-depth regions as N
    6. Generate consensus mt FASTA

Options:
  -r  mitochondrial reference FASTA
  -o  output directory
  -t  threads, default 8
EOF
  exit 1
}

MT_REF=""
OUT_DIR=""
THREADS=8

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

if [[ -z "$MT_REF" || -z "$OUT_DIR" || "$#" -lt 1 ]]; then
  usage
fi

BAMS=("$@")

for tool in samtools bwa bcftools bgzip; do
  command -v "$tool" >/dev/null 2>&1 || {
    echo "[ERROR] $tool not found in PATH" >&2
    exit 1
  }
done

[[ -f "$MT_REF" ]] || {
  echo "[ERROR] mt reference not found: $MT_REF" >&2
  exit 1
}

mkdir -p "$OUT_DIR"

MERGED_FASTA="${OUT_DIR}/all_samples_mt.fa"
: > "$MERGED_FASTA"

# index mt reference
[[ -f "${MT_REF}.fai" ]] || samtools faidx "$MT_REF"
[[ -f "${MT_REF}.bwt" ]] || bwa index "$MT_REF"

for BAM in "${BAMS[@]}"; do
  [[ -f "$BAM" ]] || {
    echo "[WARN] BAM not found, skip: $BAM" >&2
    continue
  }

  BASENAME=$(basename "$BAM")
  SAMPLE=${BASENAME%%.*}

  echo "======================================"
  echo "[INFO] Processing $SAMPLE"
  echo "[INFO] BAM: $BAM"

  PREFIX="${OUT_DIR}/${SAMPLE}"

  FQ1="${PREFIX}.R1.fq.gz"
  FQ2="${PREFIX}.R2.fq.gz"

  RAW_MT_BAM="${PREFIX}.mt.raw.bam"
  MT_BAM="${PREFIX}.mt.bam"
  VCF="${PREFIX}.mt.vcf.gz"
  LOWDEPTH_BED="${PREFIX}.mt.lowdepth.bed"
  FASTA="${PREFIX}_mt.fa"

  echo "[INFO] BAM -> FASTQ"

  samtools fastq -@ "$THREADS" \
    -1 >(gzip -c > "$FQ1") \
    -2 >(gzip -c > "$FQ2") \
    -0 /dev/null \
    -s /dev/null \
    -n "$BAM"

  echo "[INFO] Align reads to mitochondrial reference"

  bwa mem -t "$THREADS" "$MT_REF" "$FQ1" "$FQ2" | \
    samtools sort -@ "$THREADS" -o "$RAW_MT_BAM"

  samtools index "$RAW_MT_BAM"

  echo "[INFO] Filter mt alignments: MAPQ >= 30"

  samtools view -@ "$THREADS" -b -q 30 -F 4 "$RAW_MT_BAM" | \
    samtools sort -@ "$THREADS" -o "$MT_BAM"

  samtools index "$MT_BAM"

  echo "[INFO] Call haploid mt variants"

  bcftools mpileup -Ou \
    -f "$MT_REF" \
    -q 30 \
    -Q 20 \
    "$MT_BAM" | \
  bcftools call \
    -mv \
    --ploidy 1 \
    -Oz \
    -o "$VCF"

  bcftools index -t "$VCF"

  echo "[INFO] Build low-depth mask: DP < 5"

  samtools depth -a "$MT_BAM" | \
    awk 'BEGIN{OFS="\t"} $3 < 5 {print $1, $2-1, $2}' \
    > "$LOWDEPTH_BED"

  echo "[INFO] Build consensus FASTA"

  bcftools consensus \
    -f "$MT_REF" \
    -m "$LOWDEPTH_BED" \
    "$VCF" \
    > "$FASTA"

  # rename FASTA header to sample name
  awk -v s="$SAMPLE" '
    /^>/ {print ">"s; next}
    {print}
  ' "$FASTA" > "${FASTA}.tmp"

  mv "${FASTA}.tmp" "$FASTA"

  cat "$FASTA" >> "$MERGED_FASTA"

  echo "[INFO] Clean FASTQ and raw BAM"
  rm -f "$FQ1" "$FQ2" "$RAW_MT_BAM" "${RAW_MT_BAM}.bai"

  echo "[OK] Kept:"
  echo "  $MT_BAM"
  echo "  ${MT_BAM}.bai"
  echo "  $VCF"
  echo "  ${VCF}.tbi"
  echo "  $LOWDEPTH_BED"
  echo "  $FASTA"
done

echo "======================================"
echo "[INFO] All done"
echo "[INFO] Merged FASTA: $MERGED_FASTA"



bash extract_mt_from_bam.sh  -r sequence.fasta  -o mt_result -t 150 ./*.bam > mt.log 2>&1 &
