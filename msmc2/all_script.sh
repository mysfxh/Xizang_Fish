
#!/usr/bin/env bash
set -euo pipefail

# ========== 参数 ==========
REF="/home/xiongh/2026/Fish/fitered/new_1/tos/reference/HiC.review.assembly.chr.fa"
WORKDIR="/home/xiongh/2026/Fish/fitered/new_1/msmc2/mask_work"
SEQBILITY="/home/xiongh/2026/Fish/fitered/new_1/msmc2/misc-master/seq/seqbility"
MSMCTOOLS="/home/xiongh/2026/Fish/fitered/new_1/msmc2/msmc-tools-master"

K=35
THREADS=8
PREFIX="ref"

mkdir -p "${WORKDIR}"
cd "${WORKDIR}"

# ========== 1. 准备参考 ==========
cp "${REF}" ./
REFNAME=$(basename "${REF}")

# 建索引
bwa index "${REFNAME}"

# ========== 2. 切分为 overlapping k-mers ==========
"${SEQBILITY}/splitfa" "${REFNAME}" ${K} | split -l 20000000 - "${PREFIX}.split."

# 合并 split 结果
cat ${PREFIX}.split.* > "${PREFIX}.split.${K}"

# ========== 3. 回贴到参考 ==========
bwa aln -t ${THREADS} -R 1000000 -O 3 -E 3 \
  "${REFNAME}" "${PREFIX}.split.${K}" > "${PREFIX}.split.${K}.sai"

bwa samse -f "${PREFIX}.split.${K}.sam" \
  "${REFNAME}" \
  "${PREFIX}.split.${K}.sai" \
  "${PREFIX}.split.${K}"

# ========== 4. 生成 raw mask ==========
"${SEQBILITY}/gen_raw_mask.pl" \
  "${PREFIX}.split.${K}.sam" \
  > "${PREFIX}_rawMask.${K}.fa"

# ========== 5. 生成 final mappability mask ==========
"${SEQBILITY}/gen_mask" \
  -l ${K} -r 0.5 \
  "${PREFIX}_rawMask.${K}.fa" \
  > "${PREFIX}_mask.${K}.50.fa"

# ========== 6. 转成 MSMC 用的 BED ==========
python3 "${MSMCTOOLS}/makeMappabilityMask.py" \
  "${PREFIX}_mask.${K}.50.fa" \
  "${PREFIX}_mask.${K}.50"


#!/usr/bin/env bash
set -euo pipefail

########################################
# 配置
########################################

REF="/home/xiongh/2026/Fish/fitered/new_1/tos/reference/HiC.review.assembly.chr.fa"
MSMCTOOLS="/home/xiongh/2026/Fish/fitered/new_1/msmc2/msmc-tools-master"
MASKDIR="/home/xiongh/2026/Fish/fitered/new_1/msmc2/mask_work"
WORKDIR="/home/xiongh/2026/Fish/fitered/new_1/msmc2"

# 是否删除旧的 MSMC2 中间结果
# 1 = 删除旧文件后重跑
# 0 = 不删除，已有文件可能被跳过
CLEAN_OLD=1

# 参考 mappability mask 的含义：
# 0 = 正向 mask，表示这些区域可以用，推荐默认
# 1 = 负向 mask，只有当 ref_mask.* 明确表示这些区域要排除时才设为 1
USE_NEGATIVE_GENOMASK=0

# samtools / bcftools 过滤阈值
MIN_MAPQ=20
MIN_BASEQ=20
ADJUST_MQ=50

# 并行数
NPROC_STEP1=6
NPROC_STEP2=6
NPROC_STEP3=6

# 染色体列表：LG1 ~ LG20
# 如果你的参考是 LG1~LG25，就改成：CHR_LIST=($(seq 1 25))
CHR_LIST=($(seq 1 25))

########################################
# 样本定义
########################################

LAKE_IDS=(L01 L02 L03 L04 L05 L06 L07 L08 L09 L10 L11 L12)

declare -A LAKES
LAKES["L01"]="/home/xiongh/2026/Fish/bam/bam/BC-11_markdup.bam /home/xiongh/2026/Fish/bam/bam/BC-13_markdup.bam"
LAKES["L02"]="/home/xiongh/2026/Fish/bam/bam/PC-1_markdup.bam /home/xiongh/2026/Fish/bam/bam/PC-2_markdup.bam"
LAKES["L03"]="/home/xiongh/2026/Fish/bam/bam/CEND-5_markdup.bam /home/xiongh/2026/Fish/bam/bam/CEND-8_markdup.bam"
LAKES["L04"]="/home/xiongh/2026/Fish/bam/bam/CN-1_markdup.bam /home/xiongh/2026/Fish/bam/bam/CN-3_markdup.bam"
LAKES["L05"]="/home/xiongh/2026/Fish/bam/bam/GN-1_markdup.bam /home/xiongh/2026/Fish/bam/bam/GN-2_markdup.bam"
LAKES["L06"]="/home/xiongh/2026/Fish/bam/bam/GRC-3_markdup.bam /home/xiongh/2026/Fish/bam/bam/GRC-8_markdup.bam"
LAKES["L07"]="/home/xiongh/2026/Fish/bam/bam/MJC-YWGYQ-2_markdup.bam /home/xiongh/2026/Fish/bam/bam/MJC-YWGYQ-3_markdup.bam"
LAKES["L08"]="/home/xiongh/2026/Fish/bam/bam/QCG-8_markdup.bam /home/xiongh/2026/Fish/bam/bam/QCG-9_markdup.bam"
LAKES["L09"]="/home/xiongh/2026/Fish/bam/bam/SLC-8_markdup.bam /home/xiongh/2026/Fish/bam/bam/SLC-9_markdup.bam"
LAKES["L10"]="/home/xiongh/2026/Fish/bam/bam/YC-YWGYQ-1_markdup.bam /home/xiongh/2026/Fish/bam/bam/YC-YWGYQ-2_markdup.bam"
LAKES["L11"]="/home/xiongh/2026/Fish/bam/bam/YQC-YWGYQ-1_markdup.bam /home/xiongh/2026/Fish/bam/bam/YQC-YWGYQ-2_markdup.bam"
LAKES["L12"]="/home/xiongh/2026/Fish/bam/bam/ZGC-1_markdup.bam /home/xiongh/2026/Fish/bam/bam/ZGC-4_markdup.bam"

########################################
# 输出目录
########################################

mkdir -p \
  "${WORKDIR}/sites" \
  "${WORKDIR}/input" \
  "${WORKDIR}/results" \
  "${WORKDIR}/logs"

########################################
# 工具函数
########################################

log() {
    echo "[$(date '+%F %T')] $*"
}

check_cmd() {
    command -v "$1" >/dev/null 2>&1 || {
        echo "错误: 未找到命令 $1" >&2
        exit 1
    }
}

check_file() {
    [[ -f "$1" ]] || {
        echo "错误: 文件不存在 -> $1" >&2
        exit 1
    }
}

wait_for_slot() {
    local max_jobs="$1"
    while true; do
        local running
        running=$(jobs -rp | wc -l)
        if (( running < max_jobs )); then
            break
        fi
        wait -n
    done
}

make_empty_vcf_and_empty_mask() {
    local sample_name="$1"
    local chr="$2"
    local vcfgz="$3"
    local maskgz="$4"

    # 空 mask：表示该样本该染色体没有任何 callable 区域
    : | gzip -c > "$maskgz"

    {
        echo "##fileformat=VCFv4.2"
        echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t${sample_name}"
    } | gzip -c > "$vcfgz"
}

########################################
# 删除旧文件
########################################

if [[ "$CLEAN_OLD" -eq 1 ]]; then
    log "=== 删除旧的 MSMC2 中间文件 ==="

    rm -f "${WORKDIR}"/sites/*.vcf.gz
    rm -f "${WORKDIR}"/sites/*.mask.bed.gz
    rm -f "${WORKDIR}"/input/*.multihetsep.txt
    rm -f "${WORKDIR}"/input/*.err
    rm -f "${WORKDIR}"/logs/*.log

    log "旧文件删除完成"
fi

########################################
# 预检查
########################################

log "=== 预检查 ==="

check_cmd samtools
check_cmd bcftools
check_cmd python3
check_cmd gzip
check_cmd awk

check_file "$REF"
check_file "${MSMCTOOLS}/bamCaller.py"
check_file "${MSMCTOOLS}/generate_multihetsep.py"

if [[ ! -f "${REF}.fai" ]]; then
    log "参考基因组未建立 faidx，开始建立..."
    samtools faidx "$REF"
fi

for lake in "${LAKE_IDS[@]}"; do
    read -r bam1 bam2 <<< "${LAKES[$lake]}"
    check_file "$bam1"
    check_file "$bam2"
done

for chr in "${CHR_LIST[@]}"; do
    lg="LG${chr}"
    check_file "${MASKDIR}/ref_mask.35.50.${lg}.mask.bed.gz"
done

log "预检查完成"

########################################
# Step 1: 每个样本每条染色体生成 VCF.gz + mask.bed.gz
########################################

build_sample_chr() {
    local bam="$1"
    local sample_name="$2"
    local lg="$3"

    local vcfgz="${WORKDIR}/sites/${sample_name}.${lg}.vcf.gz"
    local maskgz="${WORKDIR}/sites/${sample_name}.${lg}.mask.bed.gz"
    local logfile="${WORKDIR}/logs/${sample_name}.${lg}.step1.log"

    if [[ -s "$vcfgz" && -s "$maskgz" ]]; then
        echo "[SKIP] ${sample_name} ${lg} 已存在" >> "$logfile"
        return 0
    fi

    if [[ ! -f "${bam}.bai" && ! -f "${bam%.bam}.bai" ]]; then
        echo "[INFO] 为 ${bam} 建索引" >> "$logfile"
        samtools index "$bam" >> "$logfile" 2>&1
    fi

    local depth

    # 使用 -a，把 0 覆盖位点也计入平均深度
    depth=$(samtools depth -a -r "$lg" "$bam" | \
        awk '{s+=$3;n++} END {if(n>0) printf "%.6f\n", s/n; else print "0"}')

    echo "[INFO] ${sample_name} ${lg} mean_depth=${depth}" >> "$logfile"

    if awk -v d="$depth" 'BEGIN{exit !(d<=0)}'; then
        make_empty_vcf_and_empty_mask "$sample_name" "$lg" "$vcfgz" "$maskgz"
        echo "[WARN] ${sample_name} ${lg} 深度为 0，已生成 empty mask + empty VCF，不参与该染色体分析" >> "$logfile"
        return 0
    fi

    samtools mpileup \
        -q "$MIN_MAPQ" \
        -Q "$MIN_BASEQ" \
        -C "$ADJUST_MQ" \
        -u \
        -r "$lg" \
        -f "$REF" \
        "$bam" 2>> "$logfile" | \
    bcftools call \
        -c \
        -V indels 2>> "$logfile" | \
    python3 "${MSMCTOOLS}/bamCaller.py" \
        "$depth" \
        "$maskgz" 2>> "$logfile" | \
    gzip -c > "$vcfgz"

    if [[ ! -s "$vcfgz" ]]; then
        echo "错误: 生成失败，VCF 为空 -> $vcfgz" >> "$logfile"
        return 1
    fi

    if [[ ! -s "$maskgz" ]]; then
        echo "错误: 生成失败，mask 为空 -> $maskgz" >> "$logfile"
        return 1
    fi

    echo "[OK] ${sample_name} ${lg} 完成" >> "$logfile"
}

log "=== Step 1: 处理 24 个个体，生成单样本单染色体 VCF + mask ==="

for lake in "${LAKE_IDS[@]}"; do
    read -r bam1 bam2 <<< "${LAKES[$lake]}"

    for idx in 1 2; do
        if [[ "$idx" == "1" ]]; then
            bam="$bam1"
        else
            bam="$bam2"
        fi

        sample_name="${lake}_${idx}"
        log "处理样本 ${sample_name}"

        for chr in "${CHR_LIST[@]}"; do
            lg="LG${chr}"
            wait_for_slot "$NPROC_STEP1"
            build_sample_chr "$bam" "$sample_name" "$lg" &
        done

        wait
    done
done

log "Step 1 完成"

########################################
# Step 2: 生成 12 个湖泊的 single-population multihetsep
########################################

run_singlepop_chr() {
    local lake="$1"
    local lg="$2"

    local genomask="${MASKDIR}/ref_mask.35.50.${lg}.mask.bed.gz"
    local outtxt="${WORKDIR}/input/${lake}.${lg}.multihetsep.txt"
    local errtxt="${WORKDIR}/input/${lake}.${lg}.err"

    local vcf1="${WORKDIR}/sites/${lake}_1.${lg}.vcf.gz"
    local vcf2="${WORKDIR}/sites/${lake}_2.${lg}.vcf.gz"
    local mask1="${WORKDIR}/sites/${lake}_1.${lg}.mask.bed.gz"
    local mask2="${WORKDIR}/sites/${lake}_2.${lg}.mask.bed.gz"

    check_file "$genomask"
    check_file "$vcf1"
    check_file "$vcf2"
    check_file "$mask1"
    check_file "$mask2"

    if [[ "$USE_NEGATIVE_GENOMASK" -eq 1 ]]; then
        python3 "${MSMCTOOLS}/generate_multihetsep.py" \
            --chr "$lg" \
            --negative_mask "$genomask" \
            --mask "$mask1" \
            --mask "$mask2" \
            "$vcf1" "$vcf2" \
            > "$outtxt" 2> "$errtxt"
    else
        python3 "${MSMCTOOLS}/generate_multihetsep.py" \
            --chr "$lg" \
            --mask "$genomask" \
            --mask "$mask1" \
            --mask "$mask2" \
            "$vcf1" "$vcf2" \
            > "$outtxt" 2> "$errtxt"
    fi
}

log "=== Step 2: 生成 12 个湖泊的 single-population 输入 ==="

for lake in "${LAKE_IDS[@]}"; do
    log "生成 ${lake} ..."

    for chr in "${CHR_LIST[@]}"; do
        lg="LG${chr}"
        wait_for_slot "$NPROC_STEP2"
        run_singlepop_chr "$lake" "$lg" &
    done

    wait

    final="${WORKDIR}/input/${lake}.multihetsep.txt"
    : > "$final"

    for chr in "${CHR_LIST[@]}"; do
        lg="LG${chr}"
        part="${WORKDIR}/input/${lake}.${lg}.multihetsep.txt"
        check_file "$part"
        cat "$part" >> "$final"
    done

    log "完成 ${lake}.multihetsep.txt -> $(ls -lh "$final" | awk '{print $5}')"
done

log "Step 2 完成"

########################################
# Step 3: 生成 66 对 cross-population multihetsep
########################################

run_crosspop_chr() {
    local l1="$1"
    local l2="$2"
    local lg="$3"

    local genomask="${MASKDIR}/ref_mask.35.50.${lg}.mask.bed.gz"
    local outtxt="${WORKDIR}/input/${l1}_${l2}.${lg}.multihetsep.txt"
    local errtxt="${WORKDIR}/input/${l1}_${l2}.${lg}.err"

    local vcf1="${WORKDIR}/sites/${l1}_1.${lg}.vcf.gz"
    local vcf2="${WORKDIR}/sites/${l1}_2.${lg}.vcf.gz"
    local vcf3="${WORKDIR}/sites/${l2}_1.${lg}.vcf.gz"
    local vcf4="${WORKDIR}/sites/${l2}_2.${lg}.vcf.gz"

    local mask1="${WORKDIR}/sites/${l1}_1.${lg}.mask.bed.gz"
    local mask2="${WORKDIR}/sites/${l1}_2.${lg}.mask.bed.gz"
    local mask3="${WORKDIR}/sites/${l2}_1.${lg}.mask.bed.gz"
    local mask4="${WORKDIR}/sites/${l2}_2.${lg}.mask.bed.gz"

    check_file "$genomask"
    check_file "$vcf1"
    check_file "$vcf2"
    check_file "$vcf3"
    check_file "$vcf4"
    check_file "$mask1"
    check_file "$mask2"
    check_file "$mask3"
    check_file "$mask4"

    if [[ "$USE_NEGATIVE_GENOMASK" -eq 1 ]]; then
        python3 "${MSMCTOOLS}/generate_multihetsep.py" \
            --chr "$lg" \
            --negative_mask "$genomask" \
            --mask "$mask1" \
            --mask "$mask2" \
            --mask "$mask3" \
            --mask "$mask4" \
            "$vcf1" "$vcf2" "$vcf3" "$vcf4" \
            > "$outtxt" 2> "$errtxt"
    else
        python3 "${MSMCTOOLS}/generate_multihetsep.py" \
            --chr "$lg" \
            --mask "$genomask" \
            --mask "$mask1" \
            --mask "$mask2" \
            --mask "$mask3" \
            --mask "$mask4" \
            "$vcf1" "$vcf2" "$vcf3" "$vcf4" \
            > "$outtxt" 2> "$errtxt"
    fi
}

log "=== Step 3: 生成 66 对 cross-population 输入 ==="

for ((i=0; i<${#LAKE_IDS[@]}; i++)); do
    for ((j=i+1; j<${#LAKE_IDS[@]}; j++)); do
        l1="${LAKE_IDS[$i]}"
        l2="${LAKE_IDS[$j]}"

        log "生成 ${l1}_vs_${l2} ..."

        for chr in "${CHR_LIST[@]}"; do
            lg="LG${chr}"
            wait_for_slot "$NPROC_STEP3"
            run_crosspop_chr "$l1" "$l2" "$lg" &
        done

        wait

        final="${WORKDIR}/input/${l1}_${l2}.multihetsep.txt"
        : > "$final"

        for chr in "${CHR_LIST[@]}"; do
            lg="LG${chr}"
            part="${WORKDIR}/input/${l1}_${l2}.${lg}.multihetsep.txt"
            check_file "$part"
            cat "$part" >> "$final"
        done

        log "完成 ${l1}_${l2}.multihetsep.txt -> $(ls -lh "$final" | awk '{print $5}')"
    done
done

log "Step 3 完成"

########################################
# 汇总
########################################

log ""
log "=== 所有输入文件大小汇总 ==="

log "--- 单群体文件 ---"
for lake in "${LAKE_IDS[@]}"; do
    f="${WORKDIR}/input/${lake}.multihetsep.txt"
    if [[ -f "$f" ]]; then
        echo "$(basename "$f") $(ls -lh "$f" | awk '{print $5}')"
    fi
done

log ""
log "--- 跨群体文件 ---"
for ((i=0; i<${#LAKE_IDS[@]}; i++)); do
    for ((j=i+1; j<${#LAKE_IDS[@]}; j++)); do
        f="${WORKDIR}/input/${LAKE_IDS[$i]}_${LAKE_IDS[$j]}.multihetsep.txt"
        if [[ -f "$f" ]]; then
            echo "$(basename "$f") $(ls -lh "$f" | awk '{print $5}')"
        fi
    done
done

log ""
log "=== 所有 MSMC2 输入文件准备完成 ==="





#!/usr/bin/env bash
set -euo pipefail

########################################
# MSMC2 pairwise cross-population only
# 不合并 multihetsep，直接传入每条染色体文件
########################################

INPUTDIR="/home/xiongh/2026/Fish/fitered/new_1/msmc2/input"
OUTDIR="/home/xiongh/2026/Fish/fitered/new_1/msmc2/results"
MSMCTOOLS="//home/xiongh/2026/Fish/fitered/new_1/msmc2/msmc-tools-master"

MSMC2_BIN="msmc2_Linux"

THREADS=100

# 你的染色体范围
CHR_LIST=($(seq 1 25))

mkdir -p "${OUTDIR}"

POPS=(L01 L02 L03 L04 L05 L06 L07 L08 L09 L10 L11 L12)

WITHIN_POP1_INDICES="0,1,2,3"
WITHIN_POP2_INDICES="4,5,6,7"
CROSS_INDICES="0-4,0-5,0-6,0-7,1-4,1-5,1-6,1-7,2-4,2-5,2-6,2-7,3-4,3-5,3-6,3-7"

########################################
# Check commands and files
########################################

command -v "${MSMC2_BIN}" >/dev/null 2>&1 || {
    echo "[ERROR] 找不到 ${MSMC2_BIN}，请确认 msmc2_Linux 在 PATH 中" >&2
    exit 1
}

command -v python3 >/dev/null 2>&1 || {
    echo "[ERROR] 找不到 python3" >&2
    exit 1
}

if [[ ! -f "${MSMCTOOLS}/combineCrossCoal.py" ]]; then
    echo "[ERROR] 找不到 combineCrossCoal.py: ${MSMCTOOLS}/combineCrossCoal.py" >&2
    exit 1
fi

########################################
# Pairwise cross-population MSMC2
########################################

echo "========== Pairwise cross-population MSMC2 only =========="
echo "Input mode: per-chromosome multihetsep files, not merged"

for ((i=0; i<${#POPS[@]}; i++)); do
    for ((j=i+1; j<${#POPS[@]}; j++)); do

        POP1="${POPS[$i]}"
        POP2="${POPS[$j]}"
        PAIR="${POP1}_${POP2}"

        # 收集每条染色体文件，保持 LG1, LG2, ..., LG25 顺序
        input_files=()
        for chr in "${CHR_LIST[@]}"; do
            f="${INPUTDIR}/${PAIR}.LG${chr}.multihetsep.txt"
            if [[ ! -s "$f" ]]; then
                echo "[ERROR] pair chromosome input file not found or empty: $f" >&2
                exit 1
            fi
            input_files+=("$f")
        done

        prefix_within1="${OUTDIR}/${PAIR}.within_${POP1}"
        prefix_within2="${OUTDIR}/${PAIR}.within_${POP2}"
        prefix_across="${OUTDIR}/${PAIR}.across"

        within1_final="${prefix_within1}.final.txt"
        within2_final="${prefix_within2}.final.txt"
        across_final="${prefix_across}.final.txt"

        combined="${OUTDIR}/${PAIR}.combined.txt"
        donefile="${OUTDIR}/${PAIR}.combined.done"

        if [[ -s "$donefile" && -s "$combined" ]]; then
            echo "[SKIP] ${PAIR} already combined: ${combined}"
            continue
        fi

        echo
        echo "========== Running pair: ${PAIR} =========="
        echo "  input files: ${#input_files[@]} chromosomes"
        echo "  first input: ${input_files[0]}"
        echo "  last input:  ${input_files[-1]}"

        ########################################
        # 1. within POP1
        ########################################

        if [[ ! -s "$within1_final" ]]; then
            echo "Running within ${POP1} for ${PAIR}"
            echo "  -I: ${WITHIN_POP1_INDICES}"
            echo "  output: ${prefix_within1}"

            rm -f \
                "${prefix_within1}.final.txt" \
                "${prefix_within1}.loop.txt" \
                "${prefix_within1}.log" \
                "${prefix_within1}.stdout.log" \
                "${prefix_within1}.stderr.log"

            /usr/bin/time -v "${MSMC2_BIN}" \
                -t "${THREADS}" \
                -I "${WITHIN_POP1_INDICES}" \
                -o "${prefix_within1}" \
                "${input_files[@]}" \
                > "${prefix_within1}.stdout.log" \
                2> "${prefix_within1}.stderr.log"

            if [[ ! -s "$within1_final" ]]; then
                echo "[ERROR] MSMC2 failed: ${within1_final}" >&2
                echo "Check log: ${prefix_within1}.stderr.log" >&2
                exit 1
            fi
        else
            echo "[SKIP] within ${POP1} exists: ${within1_final}"
        fi

        ########################################
        # 2. within POP2
        ########################################

        if [[ ! -s "$within2_final" ]]; then
            echo "Running within ${POP2} for ${PAIR}"
            echo "  -I: ${WITHIN_POP2_INDICES}"
            echo "  output: ${prefix_within2}"

            rm -f \
                "${prefix_within2}.final.txt" \
                "${prefix_within2}.loop.txt" \
                "${prefix_within2}.log" \
                "${prefix_within2}.stdout.log" \
                "${prefix_within2}.stderr.log"

            /usr/bin/time -v "${MSMC2_BIN}" \
                -t "${THREADS}" \
                -I "${WITHIN_POP2_INDICES}" \
                -o "${prefix_within2}" \
                "${input_files[@]}" \
                > "${prefix_within2}.stdout.log" \
                2> "${prefix_within2}.stderr.log"

            if [[ ! -s "$within2_final" ]]; then
                echo "[ERROR] MSMC2 failed: ${within2_final}" >&2
                echo "Check log: ${prefix_within2}.stderr.log" >&2
                exit 1
            fi
        else
            echo "[SKIP] within ${POP2} exists: ${within2_final}"
        fi

        ########################################
        # 3. across POP1-POP2
        ########################################

        if [[ ! -s "$across_final" ]]; then
            echo "Running across ${POP1}-${POP2}"
            echo "  -I: ${CROSS_INDICES}"
            echo "  output: ${prefix_across}"

            rm -f \
                "${prefix_across}.final.txt" \
                "${prefix_across}.loop.txt" \
                "${prefix_across}.log" \
                "${prefix_across}.stdout.log" \
                "${prefix_across}.stderr.log"

            /usr/bin/time -v "${MSMC2_BIN}" \
                -t "${THREADS}" \
                -I "${CROSS_INDICES}" \
                -o "${prefix_across}" \
                "${input_files[@]}" \
                > "${prefix_across}.stdout.log" \
                2> "${prefix_across}.stderr.log"

            if [[ ! -s "$across_final" ]]; then
                echo "[ERROR] MSMC2 failed: ${across_final}" >&2
                echo "Check log: ${prefix_across}.stderr.log" >&2
                exit 1
            fi
        else
            echo "[SKIP] across exists: ${across_final}"
        fi

        ########################################
        # 4. combineCrossCoal
        ########################################

        echo "Combining cross-coalescence for ${PAIR}"

        python3 "${MSMCTOOLS}/combineCrossCoal.py" \
            "${across_final}" \
            "${within1_final}" \
            "${within2_final}" \
            > "${combined}"

        if [[ ! -s "$combined" ]]; then
            echo "[ERROR] combineCrossCoal failed: ${combined}" >&2
            exit 1
        fi

        touch "${donefile}"
        echo "[OK] ${PAIR} combined: ${combined}"

    done
done

echo
echo "All pairwise cross-population MSMC2 analyses completed!"
