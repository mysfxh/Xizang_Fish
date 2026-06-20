#!/usr/bin/env bash
set -euo pipefail

########################################
# MSMC2 single-population bootstrap
# 1 Mb blocks, artificial 400 Mb genome, 20 bootstraps
########################################

INPUTDIR="/home/xiongh/2026/Fish/fitered/new_1/msmc2/input"
BOOTDIR="/home/xiongh/2026/Fish/fitered/new_1/msmc2/bootstrap_singlepop"
OUTDIR="/home/xiongh/2026/Fish/fitered/new_1/msmc2/results_singlepop_bootstrap"
MSMCTOOLS="/home/xiongh/2026/Fish/fitered/new_1/msmc2/msmc-tools-master"

MSMC2_BIN="msmc2_Linux"

THREADS=100
CHR_LIST=($(seq 1 25))
POPS=(L01 L02 L03 L04 L05 L06 L07 L08 L09 L10 L11 L12)

INDICES="0,1,2,3"

N_BOOT=20
CHUNK_SIZE=1000000
CHUNKS_PER_CHR=20
N_ARTIFICIAL_CHR=20

mkdir -p "${BOOTDIR}" "${OUTDIR}"

command -v "${MSMC2_BIN}" >/dev/null 2>&1 || {
    echo "[ERROR] 找不到 ${MSMC2_BIN}，请确认 msmc2_Linux 在 PATH 中" >&2
    exit 1
}

if [[ ! -f "${MSMCTOOLS}/multihetsep_bootstrap.py" ]]; then
    echo "[ERROR] 找不到 multihetsep_bootstrap.py: ${MSMCTOOLS}/multihetsep_bootstrap.py" >&2
    exit 1
fi

for POP in "${POPS[@]}"; do

    echo
    echo "========== Bootstrap for ${POP} =========="

    input_files=()
    for chr in "${CHR_LIST[@]}"; do
        f="${INPUTDIR}/${POP}.LG${chr}.multihetsep.txt"
        if [[ ! -s "$f" ]]; then
            echo "[ERROR] input file not found or empty: $f" >&2
            exit 1
        fi
        input_files+=("$f")
    done

    pop_boot_prefix="${BOOTDIR}/${POP}/bootstrap"
    mkdir -p "${BOOTDIR}/${POP}"

    echo "[INFO] Generating bootstrap datasets for ${POP}"

    python3 "${MSMCTOOLS}/multihetsep_bootstrap.py" \
        -n "${N_BOOT}" \
        -s "${CHUNK_SIZE}" \
        --chunks_per_chromosome "${CHUNKS_PER_CHR}" \
        --nr_chromosomes "${N_ARTIFICIAL_CHR}" \
        "${pop_boot_prefix}" \
        "${input_files[@]}"

    for b in $(seq 1 "${N_BOOT}"); do

        echo "[INFO] Running MSMC2 bootstrap ${POP} ${b}"

        bdir="${BOOTDIR}/${POP}/bootstrap_${b}"
        out_prefix="${OUTDIR}/${POP}.bootstrap_${b}"
        final="${out_prefix}.final.txt"

        if [[ -s "$final" ]]; then
            echo "[SKIP] ${POP} bootstrap ${b} exists: ${final}"
            continue
        fi

        boot_inputs=()
        for chr in $(seq 1 "${N_ARTIFICIAL_CHR}"); do
            f="${bdir}/bootstrap_multihetsep.chr${chr}.txt"
            if [[ ! -s "$f" ]]; then
                echo "[ERROR] bootstrap input missing: $f" >&2
                exit 1
            fi
            boot_inputs+=("$f")
        done

        /usr/bin/time -v "${MSMC2_BIN}" \
            -t "${THREADS}" \
            -I "${INDICES}" \
            -o "${out_prefix}" \
            "${boot_inputs[@]}" \
            > "${out_prefix}.stdout.log" \
            2> "${out_prefix}.stderr.log"

        if [[ ! -s "$final" ]]; then
            echo "[ERROR] MSMC2 bootstrap failed: ${final}" >&2
            echo "Check log: ${out_prefix}.stderr.log" >&2
            exit 1
        fi

        echo "[OK] ${POP} bootstrap ${b}"

    done

done

echo
echo "All single-population bootstrap MSMC2 analyses completed!"
