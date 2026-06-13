#!/usr/bin/env bash
set -euo pipefail

#############################################
# Pairwise FST for 12 lake populations
# 使用 vcftools --weir-fst-pop
#############################################

############################
# 1. 输入路径
############################

VCF="/media/server/Elements/BC2026040063-BN260407NJ01S02N5-134个异尾高原鳅重测序数据变异分析/gvcf_joint/snapp_filter_missing075/fish.02_missing075.DP10.MAF005.vcf.gz"

# 三列格式：
# sample_id    sample_id    population
POPMAP="/media/server/Elements/BC2026040063-BN260407NJ01S02N5-134个异尾高原鳅重测序数据变异分析/gvcf_joint/snapp_filter_missing075/plink_result/population_table.txt"

OUTDIR="/media/server/Elements/BC2026040063-BN260407NJ01S02N5-134个异尾高原鳅重测序数据变异分析/gvcf_joint/FST_pairwise"

# 12 个湖泊群体，不包括 TO 外群
POPS=("BC" "CEND" "CN" "GN" "GRC" "MJC" "PC" "QCG" "SLC" "YC" "YQC" "ZGC")

THREADS=8

############################
# 2. 不需要修改
############################

mkdir -p "$OUTDIR"
mkdir -p "$OUTDIR/pop_keep"
mkdir -p "$OUTDIR/pairwise_logs"
mkdir -p "$OUTDIR/pairwise_fst_sites"

echo
echo "========================================"
echo "Pairwise FST analysis"
echo "VCF: $VCF"
echo "POPMAP: $POPMAP"
echo "OUTDIR: $OUTDIR"
echo "Populations: ${POPS[*]}"
echo "========================================"
echo

#############################################
# Step 0. 检查文件和软件
#############################################

echo "========== Step 0: check files and tools =========="

if [[ ! -s "$VCF" ]]; then
    echo "[ERROR] VCF not found:"
    echo "$VCF"
    exit 1
fi

if [[ ! -s "$POPMAP" ]]; then
    echo "[ERROR] population table not found:"
    echo "$POPMAP"
    exit 1
fi

command -v vcftools >/dev/null 2>&1 || { echo "[ERROR] vcftools not found"; exit 1; }
command -v bcftools >/dev/null 2>&1 || { echo "[ERROR] bcftools not found"; exit 1; }

if [[ ! -s "${VCF}.tbi" && ! -s "${VCF}.csi" ]]; then
    echo "[WARNING] VCF index not found. Building index..."
    tabix -f -p vcf "$VCF"
fi

echo "[OK] files and tools checked."

#############################################
# Step 1. 检查样本名
#############################################

echo
echo "========== Step 1: check sample names =========="

bcftools query -l "$VCF" | sort > "$OUTDIR/vcf.samples.txt"
awk 'NF>=3{print $1}' "$POPMAP" | sort > "$OUTDIR/popmap.samples.txt"

comm -23 "$OUTDIR/popmap.samples.txt" "$OUTDIR/vcf.samples.txt" > "$OUTDIR/samples_in_popmap_not_in_vcf.txt"
comm -13 "$OUTDIR/popmap.samples.txt" "$OUTDIR/vcf.samples.txt" > "$OUTDIR/samples_in_vcf_not_in_popmap.txt"

if [[ -s "$OUTDIR/samples_in_popmap_not_in_vcf.txt" ]]; then
    echo "[ERROR] 以下 popmap 样本不在 VCF 中："
    cat "$OUTDIR/samples_in_popmap_not_in_vcf.txt"
    exit 1
fi

echo "[OK] popmap 样本都在 VCF 中。"

if [[ -s "$OUTDIR/samples_in_vcf_not_in_popmap.txt" ]]; then
    echo "[WARNING] 以下 VCF 样本不在 popmap 中："
    cat "$OUTDIR/samples_in_vcf_not_in_popmap.txt"
fi

#############################################
# Step 2. 为每个群体生成 keep 文件
#############################################

echo
echo "========== Step 2: make population keep files =========="

echo -e "population\tn_samples" > "$OUTDIR/population_sample_counts.tsv"

for pop in "${POPS[@]}"; do
    keepfile="$OUTDIR/pop_keep/${pop}.keep"

    awk -v p="$pop" 'NF>=3 && $3==p {print $1}' "$POPMAP" > "$keepfile"

    n=$(wc -l < "$keepfile")

    echo -e "${pop}\t${n}" >> "$OUTDIR/population_sample_counts.tsv"

    if [[ "$n" -lt 2 ]]; then
        echo "[WARNING] $pop 样本数少于 2，FST 可能不稳定。n=$n"
    fi
done

cat "$OUTDIR/population_sample_counts.tsv"

#############################################
# Step 3. 计算所有两两 FST
#############################################

echo
echo "========== Step 3: calculate pairwise FST =========="

SUMMARY="$OUTDIR/pairwise_fst_summary.tsv"

echo -e "pop1\tpop2\tmean_FST\tweighted_FST\tn_sites" > "$SUMMARY"

N=${#POPS[@]}

for ((i=0; i<N; i++)); do
    for ((j=i+1; j<N; j++)); do

        pop1="${POPS[$i]}"
        pop2="${POPS[$j]}"

        keep1="$OUTDIR/pop_keep/${pop1}.keep"
        keep2="$OUTDIR/pop_keep/${pop2}.keep"

        prefix="$OUTDIR/pairwise_fst_sites/${pop1}_vs_${pop2}"

        echo
        echo "Running FST: ${pop1} vs ${pop2}"

        vcftools \
          --gzvcf "$VCF" \
          --weir-fst-pop "$keep1" \
          --weir-fst-pop "$keep2" \
          --out "$prefix" \
          > "$OUTDIR/pairwise_logs/${pop1}_vs_${pop2}.log" 2>&1

        log="$OUTDIR/pairwise_logs/${pop1}_vs_${pop2}.log"

        mean_fst=$(grep "Weir and Cockerham mean Fst estimate" "$log" | awk -F ":" '{gsub(/ /,"",$2); print $2}')
        weighted_fst=$(grep "Weir and Cockerham weighted Fst estimate" "$log" | awk -F ":" '{gsub(/ /,"",$2); print $2}')

        if [[ -s "${prefix}.weir.fst" ]]; then
            n_sites=$(awk 'NR>1 && $3!="nan" {count++} END{print count+0}' "${prefix}.weir.fst")
        else
            n_sites="0"
        fi

        if [[ -z "${mean_fst}" ]]; then
            mean_fst="NA"
        fi

        if [[ -z "${weighted_fst}" ]]; then
            weighted_fst="NA"
        fi

        echo -e "${pop1}\t${pop2}\t${mean_fst}\t${weighted_fst}\t${n_sites}" >> "$SUMMARY"

    done
done

echo
echo "[OK] Pairwise FST finished."
echo "Summary:"
echo "$SUMMARY"

#############################################
# Step 4. 生成 FST 矩阵
#############################################

echo
echo "========== Step 4: make FST matrix =========="

python3 - <<PY
import pandas as pd
import numpy as np

summary = "${SUMMARY}"
out_matrix_weighted = "${OUTDIR}/pairwise_weighted_fst_matrix.tsv"
out_matrix_mean = "${OUTDIR}/pairwise_mean_fst_matrix.tsv"

df = pd.read_csv(summary, sep="\t")

pops = ${POPS[@]+"[" + ", ".join([repr(p) for p in POPS]) + "]"} if False else None
PY
