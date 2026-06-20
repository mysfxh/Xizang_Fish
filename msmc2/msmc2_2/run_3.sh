#!/usr/bin/env bash
set -euo pipefail

########################################
# MSMC2 pairwise cross-population
# unphased-friendly version
# Output to a new result folder
########################################

INPUTDIR="/home/xiongh/2026/Fish/fitered/new_1/msmc2/input"
OUTDIR="/home/xiongh/2026/Fish/fitered/new_1/msmc2/results_s_unphased"
MSMCTOOLS="/home/xiongh/2026/Fish/fitered/new_1/msmc2/msmc-tools-master"

MSMC2_BIN="msmc2_Linux"

THREADS=100
CHR_LIST=($(seq 1 25))

mkdir -p "${OUTDIR}"

POPS=(L01 L02 L03 L04 L05 L06 L07 L08 L09 L10 L11 L12)

# unphased genomes:
# each diploid individual contributes one usable within-individual pair
# POP1: individual1 = 0-1, individual2 = 2-3
# POP2: individual1 = 4-5, individual2 = 6-7
WITHIN_POP1_INDICES="0-1,2-3"
WITHIN_POP2_INDICES="4-5,6-7"

# cross-population pairs
CROSS_INDICES="0-4,0-5,0-6,0-7,1-4,1-5,1-6,1-7,2-4,2-5,2-6,2-7,3-4,3-5,3-6,3-7"

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

echo "========== Pairwise cross-population MSMC2 =========="
echo "Mode: unphased-friendly + -s skipAmbiguous"
echo "Output directory: ${OUTDIR}"

for ((i=0; i<${#POPS[@]}; i++)); do
    for ((j=i+1; j<${#POPS[@]}; j++)); do

        POP1="${POPS[$i]}"
        POP2="${POPS[$j]}"
        PAIR="${POP1}_${POP2}"

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
        echo "Input chromosomes: ${#input_files[@]}"
        echo "Within ${POP1}: ${WITHIN_POP1_INDICES}"
        echo "Within ${POP2}: ${WITHIN_POP2_INDICES}"
        echo "Across: ${CROSS_INDICES}"

        ########################################
        # 1. within POP1
        ########################################

        if [[ ! -s "$within1_final" ]]; then
            echo "[RUN] within ${POP1} for ${PAIR}"

            rm -f \
                "${prefix_within1}.final.txt" \
                "${prefix_within1}.loop.txt" \
                "${prefix_within1}.log" \
                "${prefix_within1}.stdout.log" \
                "${prefix_within1}.stderr.log"

            /usr/bin/time -v "${MSMC2_BIN}" \
                -t "${THREADS}" \
                -s \
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
            echo "[RUN] within ${POP2} for ${PAIR}"

            rm -f \
                "${prefix_within2}.final.txt" \
                "${prefix_within2}.loop.txt" \
                "${prefix_within2}.log" \
                "${prefix_within2}.stdout.log" \
                "${prefix_within2}.stderr.log"

            /usr/bin/time -v "${MSMC2_BIN}" \
                -t "${THREADS}" \
                -s \
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
            echo "[RUN] across ${POP1}-${POP2}"

            rm -f \
                "${prefix_across}.final.txt" \
                "${prefix_across}.loop.txt" \
                "${prefix_across}.log" \
                "${prefix_across}.stdout.log" \
                "${prefix_across}.stderr.log"

            /usr/bin/time -v "${MSMC2_BIN}" \
                -t "${THREADS}" \
                -s \
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

        echo "[RUN] Combining cross-coalescence for ${PAIR}"

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
echo "New results are in: ${OUTDIR}"
