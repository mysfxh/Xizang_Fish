##!/usr/bin/env bash
set -euo pipefail

#############################################
# TreeMix workflow from LD-pruned VCF
# 1. Rename LG1/LG2/... to 1/2/...
# 2. Convert VCF to TreeMix input
# 3. Run TreeMix m=0-10 with k=200
#############################################

#############################################
# 需要你修改的参数
#############################################

# 输入：已经 LD pruning 后的 VCF
VCF="/media/server/Elements/BC2026040063-BN260407NJ01S02N5-134个异尾高原鳅重测序数据变异分析/gvcf_joint/snapp_filter_missing075/plink_result/fish.02_missing075.DP10.MAF005.snapp.filtered.unique.100_1_0.8.unlinked.vcf.gz"

# population table，两列格式：
# sample_id    population
CLUST="/media/server/Elements/BC2026040063-BN260407NJ01S02N5-134个异尾高原鳅重测序数据变异分析/gvcf_joint/Treemix/population_table.txt"

# 输出目录
WORKDIR="/media/server/Elements/BC2026040063-BN260407NJ01S02N5-134个异尾高原鳅重测序数据变异分析/gvcf_joint/Treemix/work"
RESULTDIR="/media/server/Elements/BC2026040063-BN260407NJ01S02N5-134个异尾高原鳅重测序数据变异分析/gvcf_joint/Treemix/result"

# TreeMix root，必须是 population_table.txt 第二列里的群体名
# 如果没有外群，可以设置 ROOT=""
ROOT="Rana_chensinensis"

# TreeMix 参数
K_BLOCK=200
M_MIN=0
M_MAX=10
N_RUN=5

# 同时运行多少个 TreeMix 任务
MAX_JOBS=5

#############################################
# 不需要修改
#############################################

mkdir -p "${WORKDIR}" "${RESULTDIR}"

BASE=$(basename "${VCF}" .vcf.gz)

CHR_MAP="${WORKDIR}/${BASE}.chr.rename.txt"
RENAMED_VCF="${WORKDIR}/${BASE}.chrnum.vcf.gz"

TREEMIX_PREFIX="${RENAMED_VCF%.vcf.gz}"
TREEMIX_IN="${TREEMIX_PREFIX}.treemix.frq.gz"

echo
echo "============================================="
echo "TreeMix workflow"
echo "Input VCF: ${VCF}"
echo "Population table: ${CLUST}"
echo "Workdir: ${WORKDIR}"
echo "Resultdir: ${RESULTDIR}"
echo "Root: ${ROOT}"
echo "k block: ${K_BLOCK}"
echo "m range: ${M_MIN}-${M_MAX}"
echo "runs per m: ${N_RUN}"
echo "============================================="
echo

#############################################
# Step 0: check input files
#############################################

echo "========== Step 0: check input files =========="

if [[ ! -s "${VCF}" ]]; then
    echo "[ERROR] VCF not found:"
    echo "${VCF}"
    exit 1
fi

if [[ ! -s "${VCF}.tbi" && ! -s "${VCF}.csi" ]]; then
    echo "[WARNING] VCF index not found. Building index..."
    tabix -f -p vcf "${VCF}"
fi

if [[ ! -s "${CLUST}" ]]; then
    echo "[ERROR] population table not found:"
    echo "${CLUST}"
    exit 1
fi

echo "[OK] input files checked."

#############################################
# Step 1: make chromosome rename table
#############################################

echo
echo "========== Step 1: make chromosome rename table =========="

bcftools query -f '%CHROM\n' "${VCF}" | sort -V | uniq | \
awk '
{
    old=$1;
    new=old;

    # LG1 -> 1, LG2 -> 2
    if (new ~ /^LG[0-9]+$/) {
        sub(/^LG/, "", new);
    }
    # Chr01 -> 1, chr01 -> 1
    else if (new ~ /^[Cc]hr[0-9]+$/) {
        sub(/^[Cc]hr/, "", new);
        new = new + 0;
    }
    # chromosome_1 -> 1
    else if (new ~ /^chromosome_[0-9]+$/) {
        sub(/^chromosome_/, "", new);
    }
    # 其他非数字染色体名保留原名
    print old"\t"new;
}' > "${CHR_MAP}"

echo "Chromosome rename table:"
cat "${CHR_MAP}"

#############################################
# Step 2: rename chromosomes
#############################################

echo
echo "========== Step 2: rename chromosomes =========="

bcftools annotate \
  --rename-chrs "${CHR_MAP}" \
  -Oz -o "${RENAMED_VCF}" \
  "${VCF}"

tabix -f -p vcf "${RENAMED_VCF}"

echo "[OK] renamed VCF:"
echo "${RENAMED_VCF}"

echo "Check chromosome names after renaming:"
bcftools query -f '%CHROM\n' "${RENAMED_VCF}" | sort -V | uniq | head -n 20

#############################################
# Step 3: check sample names
#############################################

echo
echo "========== Step 3: check sample names =========="

bcftools query -l "${RENAMED_VCF}" | sort > "${WORKDIR}/${BASE}.vcf.samples.txt"
awk 'NF{print $1}' "${CLUST}" | sort > "${WORKDIR}/${BASE}.clust.samples.txt"

comm -23 "${WORKDIR}/${BASE}.clust.samples.txt" "${WORKDIR}/${BASE}.vcf.samples.txt" > "${WORKDIR}/${BASE}.samples_not_in_vcf.txt"
comm -13 "${WORKDIR}/${BASE}.clust.samples.txt" "${WORKDIR}/${BASE}.vcf.samples.txt" > "${WORKDIR}/${BASE}.samples_not_in_clust.txt"

if [[ -s "${WORKDIR}/${BASE}.samples_not_in_vcf.txt" ]]; then
    echo "[ERROR] 以下 population_table 样本不在 VCF 中："
    cat "${WORKDIR}/${BASE}.samples_not_in_vcf.txt"
    exit 1
fi

echo "[OK] population_table 样本都在 VCF 中。"

if [[ -s "${WORKDIR}/${BASE}.samples_not_in_clust.txt" ]]; then
    echo "[WARNING] 以下 VCF 样本不在 population_table 中："
    cat "${WORKDIR}/${BASE}.samples_not_in_clust.txt"
    echo "这些样本可能不会被正确分群，建议确认。"
fi

echo "Populations in population_table:"
awk 'NF{print $2}' "${CLUST}" | sort | uniq

#############################################
# Step 4: convert VCF to TreeMix input
#############################################

echo
echo "========== Step 4: convert VCF to TreeMix input =========="

vcf2treemix.sh "${RENAMED_VCF}" "${CLUST}"

if [[ ! -s "${TREEMIX_IN}" ]]; then
    echo "[ERROR] TreeMix input was not found:"
    echo "${TREEMIX_IN}"
    echo
    echo "请检查 vcf2treemix.sh 实际输出文件名。"
    echo "当前目录可能生成的 treemix 文件："
    ls -lh "${WORKDIR}"/*.treemix.frq.gz 2>/dev/null || true
    exit 1
fi

echo "[OK] TreeMix input:"
echo "${TREEMIX_IN}"

echo "TreeMix population header:"
zcat "${TREEMIX_IN}" | head -n 1

#############################################
# Step 5: check root
#############################################

echo
echo "========== Step 5: check root =========="

if [[ -n "${ROOT}" ]]; then
    if zcat "${TREEMIX_IN}" | head -n 1 | tr ' ' '\n' | grep -qx "${ROOT}"; then
        echo "[OK] Root population found: ${ROOT}"
    else
        echo "[ERROR] Root population not found in TreeMix input:"
        echo "${ROOT}"
        echo
        echo "TreeMix populations are:"
        zcat "${TREEMIX_IN}" | head -n 1
        echo
        echo "请确认 ROOT 是否是 population_table.txt 第二列里的群体名。"
        exit 1
    fi
else
    echo "[WARNING] ROOT is empty. TreeMix will run without root."
fi

#############################################
# Step 6: run TreeMix
#############################################

echo
echo "========== Step 6: run TreeMix =========="

for m in $(seq "${M_MIN}" "${M_MAX}"); do
    for run in $(seq 1 "${N_RUN}"); do

        OUT_PREFIX="${RESULTDIR}/${BASE}.chrnum.k${K_BLOCK}.m${m}_run${run}"
        LOG_FILE="${RESULTDIR}/treemix.k${K_BLOCK}.m${m}_run${run}.log"

        echo "Running TreeMix: m=${m}, run=${run}, k=${K_BLOCK}"

        if [[ -n "${ROOT}" ]]; then
            treemix \
              -i "${TREEMIX_IN}" \
              -m "${m}" \
              -o "${OUT_PREFIX}" \
              -root "${ROOT}" \
              -bootstrap \
              -k "${K_BLOCK}" \
              > "${LOG_FILE}" 2>&1 &
        else
            treemix \
              -i "${TREEMIX_IN}" \
              -m "${m}" \
              -o "${OUT_PREFIX}" \
              -bootstrap \
              -k "${K_BLOCK}" \
              > "${LOG_FILE}" 2>&1 &
        fi

        # 控制同时运行的任务数
        while [[ $(jobs -rp | wc -l) -ge "${MAX_JOBS}" ]]; do
            sleep 30
        done

    done
done

wait

echo
echo "========== All TreeMix runs finished =========="
echo "Results:"
echo "${RESULTDIR}"
echo
echo "Check likelihoods:"
echo "grep -H 'likelihood' ${RESULTDIR}/*.log"
