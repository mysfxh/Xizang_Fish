#!/usr/bin/env bash
set -euo pipefail

########################################
# SNP filtering + LD pruning for PLINK/SNAPP
# Run two missing thresholds: 0.75 and 0.95
########################################

QUAL_MIN=30
QD_MIN=2
FS_MAX=60
MQ_MIN=40
SOR_MAX=4
THREADS=16

MIN_DP=10
MAF=0.05

INPUT_VCF="./all_sample.snp.vcf.gz"

LD_WINDOW=100
LD_STEP=1
LD_R2=0.8

MISSING_LIST=("0.75" "0.95")

SITE_KEEP_EXPR="QUAL>=${QUAL_MIN} && INFO/QD>=${QD_MIN} && INFO/FS<=${FS_MAX} && INFO/MQ>=${MQ_MIN} && INFO/SOR<=${SOR_MAX}"

########################################
# Step 1: site-level filtering
# This step is shared by 0.75 and 0.95
########################################

PREFIX1="fish.01.sitefiltereds"

echo "========== Step 1: site-level filtering =========="

if [[ ! -s "${PREFIX1}.vcf.gz" ]]; then
    bcftools view --threads "${THREADS}" \
      -m2 -M2 -v snps \
      -i "${SITE_KEEP_EXPR}" \
      -Oz -o "${PREFIX1}.vcf.gz" \
      "${INPUT_VCF}"

    tabix -f -p vcf "${PREFIX1}.vcf.gz"
else
    echo "[SKIP] ${PREFIX1}.vcf.gz already exists"
fi


########################################
# Step 2-7: run two missing thresholds
########################################

for MAX_MISSING in "${MISSING_LIST[@]}"; do

    LABEL=$(echo "${MAX_MISSING}" | sed 's/\.//g')
    # 0.75 -> 075
    # 0.95 -> 095

    OUTDIR="snapp_filter_missing${LABEL}"
    PREFIX2="fish.02_missing${LABEL}.DP${MIN_DP}.MAF005"

    echo
    echo "=================================================="
    echo "Running filtering with max-missing = ${MAX_MISSING}"
    echo "Output directory: ${OUTDIR}"
    echo "Prefix: ${PREFIX2}"
    echo "=================================================="

    mkdir -p "${OUTDIR}"

    ########################################
    # Step 2: depth, missing rate, MAF filtering
    ########################################

    echo "========== Step 2: genotype depth, missing rate, MAF filtering =========="

    if [[ ! -s "${OUTDIR}/${PREFIX2}.vcf.gz" ]]; then
        vcftools \
          --gzvcf "${PREFIX1}.vcf.gz" \
          --minDP "${MIN_DP}" \
          --max-missing "${MAX_MISSING}" \
          --maf "${MAF}" \
          --min-alleles 2 \
          --max-alleles 2 \
          --recode \
          --recode-INFO-all \
          --stdout | bgzip -c > "${OUTDIR}/${PREFIX2}.vcf.gz"

        tabix -f -p vcf "${OUTDIR}/${PREFIX2}.vcf.gz"
    else
        echo "[SKIP] ${OUTDIR}/${PREFIX2}.vcf.gz already exists"
    fi


    ########################################
    # Step 3: convert VCF to PLINK bed
    ########################################

    echo "========== Step 3: convert VCF to PLINK bed =========="

    PLINKDIR="${OUTDIR}/plink_result"
    mkdir -p "${PLINKDIR}"

    plink \
      --vcf "${OUTDIR}/${PREFIX2}.vcf.gz" \
      --allow-extra-chr \
      --double-id \
      --make-bed \
      --out "${PLINKDIR}/${PREFIX2}"


    ########################################
    # Step 4: reset SNP IDs
    ########################################

    echo "========== Step 4: reset SNP IDs =========="

    plink2 \
      --bfile "${PLINKDIR}/${PREFIX2}" \
      --allow-extra-chr \
      --set-all-var-ids @:#\$r,\$a \
      --new-id-max-allele-len 1000 truncate \
      --make-bed \
      --out "${PLINKDIR}/${PREFIX2}.unique"


    ########################################
    # Step 5: LD pruning
    ########################################

    echo "========== Step 5: LD pruning =========="

    LD_PREFIX="${PLINKDIR}/${PREFIX2}.unique.${LD_WINDOW}_${LD_STEP}_${LD_R2}"

    plink \
      --bfile "${PLINKDIR}/${PREFIX2}.unique" \
      --allow-extra-chr \
      --indep-pairwise "${LD_WINDOW}" "${LD_STEP}" "${LD_R2}" \
      --out "${LD_PREFIX}"


    ########################################
    # Step 6: extract LD-pruned SNPs
    ########################################

    echo "========== Step 6: extract LD-pruned SNPs =========="

    PRUNED_PREFIX="${PLINKDIR}/${PREFIX2}.snapp.filtered.unique.${LD_WINDOW}_${LD_STEP}_${LD_R2}"

    plink \
      --bfile "${PLINKDIR}/${PREFIX2}.unique" \
      --allow-extra-chr \
      --extract "${LD_PREFIX}.prune.in" \
      --make-bed \
      --out "${PRUNED_PREFIX}"


    ########################################
    # Step 7: export LD-pruned VCF
    ########################################

    echo "========== Step 7: export LD-pruned VCF =========="

    FINAL_PREFIX="${PRUNED_PREFIX}.unlinked"

    plink \
      --bfile "${PRUNED_PREFIX}" \
      --allow-extra-chr \
      --recode vcf-iid bgz \
      --out "${FINAL_PREFIX}"

    tabix -f -p vcf "${FINAL_PREFIX}.vcf.gz"


    ########################################
    # Step 8: count SNPs and missing rate
    ########################################

    echo "========== Step 8: summary =========="

    bcftools view -H "${OUTDIR}/${PREFIX2}.vcf.gz" | wc -l > "${OUTDIR}/${PREFIX2}.sitefiltered.snp.count.txt"
    bcftools view -H "${FINAL_PREFIX}.vcf.gz" | wc -l > "${OUTDIR}/${PREFIX2}.LDpruned.snp.count.txt"

    vcftools \
      --gzvcf "${FINAL_PREFIX}.vcf.gz" \
      --missing-indv \
      --out "${FINAL_PREFIX}.missing_indv"

    vcftools \
      --gzvcf "${FINAL_PREFIX}.vcf.gz" \
      --missing-site \
      --out "${FINAL_PREFIX}.missing_site"

    echo "[OK] max-missing ${MAX_MISSING} finished"
    echo "Filtered VCF before LD pruning:"
    echo "${OUTDIR}/${PREFIX2}.vcf.gz"
    echo "Final LD-pruned VCF:"
    echo "${FINAL_PREFIX}.vcf.gz"

done

echo
echo "========== All analyses finished =========="
echo "Results:"
echo "1. snapp_filter_missing075/"
echo "2. snapp_filter_missing095/"

