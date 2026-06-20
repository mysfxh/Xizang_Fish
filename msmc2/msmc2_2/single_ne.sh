#!/usr/bin/env bash
set -euo pipefail

########################################
# MSMC2 single-population Ne history
# 两个二倍体个体 = 4 haplotypes
########################################

INPUTDIR="/home/xiongh/2026/Fish/fitered/new_1/msmc2/input"
OUTDIR="/home/xiongh/2026/Fish/fitered/new_1/msmc2/results_singlepop"

MSMC2_BIN="msmc2_Linux"

THREADS=100
CHR_LIST=($(seq 1 25))

POPS=(L01 L02 L03 L04 L05 L06 L07 L08 L09 L10 L11 L12)

# 每个湖泊两个二倍体个体 = 4 条 haplotype
INDICES="0,1,2,3"

mkdir -p "${OUTDIR}"

command -v "${MSMC2_BIN}" >/dev/null 2>&1 || {
    echo "[ERROR] 找不到 ${MSMC2_BIN}，请确认 msmc2_Linux 在 PATH 中" >&2
    exit 1
}

for POP in "${POPS[@]}"; do

    echo
    echo "========== Running single-population MSMC2: ${POP} =========="

    input_files=()
    for chr in "${CHR_LIST[@]}"; do
        f="${INPUTDIR}/${POP}.LG${chr}.multihetsep.txt"
        if [[ ! -s "$f" ]]; then
            echo "[ERROR] single-pop input file not found or empty: $f" >&2
            exit 1
        fi
        input_files+=("$f")
    done

    prefix="${OUTDIR}/${POP}"
    final="${prefix}.final.txt"

    if [[ -s "$final" ]]; then
        echo "[SKIP] ${POP} already exists: ${final}"
        continue
    fi

    rm -f \
        "${prefix}.final.txt" \
        "${prefix}.loop.txt" \
        "${prefix}.log" \
        "${prefix}.stdout.log" \
        "${prefix}.stderr.log"

    /usr/bin/time -v "${MSMC2_BIN}" \
        -t "${THREADS}" \
        -I "${INDICES}" \
        -o "${prefix}" \
        "${input_files[@]}" \
        > "${prefix}.stdout.log" \
        2> "${prefix}.stderr.log"

    if [[ ! -s "$final" ]]; then
        echo "[ERROR] MSMC2 failed for ${POP}: ${final}" >&2
        echo "Check log: ${prefix}.stderr.log" >&2
        exit 1
    fi

    echo "[OK] ${POP}: ${final}"

done

echo
echo "All single-population MSMC2 analyses completed!"
