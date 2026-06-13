#!/usr/bin/env bash
set -euo pipefail

############################################################
# Resume ADMIXTOOLS / qpDstat workflow
# 目的：
#   1. 尽量复用已经生成的 VCF，避免从头重跑
#   2. 去除非染色体位点，如 Scaffold19
#   3. 只保留 LG 主染色体或数字染色体
#   4. 继续完成 PLINK -> convertf -> qpDstat -> qp3Pop -> 作图
#
# 注意：
#   OUTDIR 虽然叫 Dsuite，但本脚本运行的是 ADMIXTOOLS：
#   convertf / qpDstat / qp3Pop / qpF4ratio
############################################################

############################
# 1. 输入路径
############################

ORIGINAL_VCF="/media/server/Elements/BC2026040063-BN260407NJ01S02N5-134个异尾高原鳅重测序数据变异分析/gvcf_joint/snapp_filter_missing075/fish.02_missing075.DP10.MAF005.vcf.gz"

POPMAP3="/media/server/Elements/BC2026040063-BN260407NJ01S02N5-134个异尾高原鳅重测序数据变异分析/gvcf_joint/snapp_filter_missing075/plink_result/population_table.txt"

OUTDIR="/media/server/Elements/BC2026040063-BN260407NJ01S02N5-134个异尾高原鳅重测序数据变异分析/gvcf_joint/Dsuite"

############################
# 2. 分析参数
############################

# 12 个湖泊群体，不包括外群 TO
LAKE_POPS="BC CEND CN GN GRC MJC PC QCG SLC YC YQC ZGC"

# 外群
OUTGROUP="TO"

# ADMIXTOOLS block size
BLGSIZE=0.01

# 新前缀，避免和之前包含 Scaffold 的 fish_admix 混淆
PREFIX="fish_admix_noScaffold"

# 如果想强制重跑所有步骤：
# FORCE=1 bash run_qpDstat_resume_noScaffold.sh
FORCE="${FORCE:-0}"

mkdir -p "$OUTDIR"
cd "$OUTDIR"

echo
echo "===================================================="
echo "Resume ADMIXTOOLS / qpDstat workflow"
echo "ORIGINAL_VCF: $ORIGINAL_VCF"
echo "POPMAP3:      $POPMAP3"
echo "OUTDIR:       $OUTDIR"
echo "Lake pops:    $LAKE_POPS"
echo "Outgroup:     $OUTGROUP"
echo "Prefix:       $PREFIX"
echo "FORCE:        $FORCE"
echo "===================================================="
echo

############################################################
# Step 0. Check files and tools
############################################################

echo "========== Step 0: check files and tools =========="

if [[ ! -s "$ORIGINAL_VCF" ]]; then
    echo "[ERROR] Original VCF not found:"
    echo "$ORIGINAL_VCF"
    exit 1
fi

if [[ ! -s "$POPMAP3" ]]; then
    echo "[ERROR] population_table.txt not found:"
    echo "$POPMAP3"
    exit 1
fi

command -v bcftools >/dev/null 2>&1 || { echo "[ERROR] bcftools not found"; exit 1; }
command -v tabix >/dev/null 2>&1 || { echo "[ERROR] tabix not found"; exit 1; }
command -v plink >/dev/null 2>&1 || { echo "[ERROR] plink not found"; exit 1; }
command -v convertf >/dev/null 2>&1 || { echo "[ERROR] convertf not found; please activate admixtools env"; exit 1; }
command -v qpDstat >/dev/null 2>&1 || { echo "[ERROR] qpDstat not found; please activate admixtools env"; exit 1; }
command -v qp3Pop >/dev/null 2>&1 || { echo "[ERROR] qp3Pop not found; please activate admixtools env"; exit 1; }

if command -v qpF4ratio >/dev/null 2>&1; then
    HAS_F4RATIO="yes"
else
    HAS_F4RATIO="no"
    echo "[WARNING] qpF4ratio not found. Optional f4-ratio step will be skipped."
fi

if command -v Rscript >/dev/null 2>&1; then
    HAS_R="yes"
else
    HAS_R="no"
    echo "[WARNING] Rscript not found. Plotting step will be skipped."
fi

echo "[OK] files and tools checked."

############################################################
# Step 1. Make two-column popmap and keep list
############################################################

echo
echo "========== Step 1: make popmap and keep list =========="

POPMAP2="${OUTDIR}/fish_12lakes_TO.popmap2.txt"
KEEP="${OUTDIR}/fish_12lakes_TO.keep.txt"

awk -v lakes="$LAKE_POPS" -v outgroup="$OUTGROUP" '
BEGIN {
    split(lakes, a, " ");
    for (i in a) keep[a[i]]=1;
}
NF >= 3 {
    sample=$1;
    pop=$3;
    if (pop == outgroup || pop in keep) {
        print sample, pop;
    }
}
' "$POPMAP3" > "$POPMAP2"

awk '{print $1}' "$POPMAP2" > "$KEEP"

sort "$KEEP" > keep.samples.txt

echo "[OK] two-column popmap:"
echo "$POPMAP2"

echo
echo "[INFO] Population counts:"
awk '{count[$2]++} END{for(p in count) print p, count[p]}' "$POPMAP2" | sort

############################################################
# Step 2. Choose existing VCF to avoid rerunning
############################################################

echo
echo "========== Step 2: choose existing VCF =========="

FINAL_VCF="${OUTDIR}/${PREFIX}.biallelicSNP.vcf.gz"

SOURCE_VCF=""

if [[ -s "$FINAL_VCF" && "$FORCE" != "1" ]]; then
    echo "[SKIP] Final no-Scaffold VCF already exists:"
    echo "$FINAL_VCF"
else
    CANDIDATE_VCFS=(
        "${OUTDIR}/fish_admix.biallelicSNP.vcf.gz"
        "${OUTDIR}/fish.02_missing075.DP10.MAF005.onlyLG.chrnum.vcf.gz"
        "${OUTDIR}/fish.02_missing075.DP10.MAF005.chrnum.vcf.gz"
        "${OUTDIR}/fish.02_missing075.DP10.MAF005.onlyLG.vcf.gz"
        "$ORIGINAL_VCF"
    )

    for f in "${CANDIDATE_VCFS[@]}"; do
        if [[ -s "$f" ]]; then
            echo "[INFO] Found candidate VCF:"
            echo "$f"

            if [[ ! -s "${f}.tbi" && ! -s "${f}.csi" ]]; then
                echo "[INFO] Index not found. Indexing candidate VCF..."
                tabix -f -p vcf "$f"
            fi

            bcftools query -l "$f" | sort > candidate.vcf.samples.txt

            comm -23 keep.samples.txt candidate.vcf.samples.txt > candidate.missing.samples.txt

            if [[ -s candidate.missing.samples.txt ]]; then
                echo "[WARNING] Candidate VCF misses some required samples. Skip this VCF."
                head candidate.missing.samples.txt
            else
                SOURCE_VCF="$f"
                echo "[OK] Selected source VCF:"
                echo "$SOURCE_VCF"
                break
            fi
        fi
    done

    if [[ -z "$SOURCE_VCF" ]]; then
        echo "[ERROR] No suitable VCF found."
        exit 1
    fi
fi

############################################################
# Step 3. Remove non-chromosome contigs and create final VCF
############################################################

echo
echo "========== Step 3: remove non-chromosome contigs =========="

if [[ ! -s "$FINAL_VCF" || "$FORCE" == "1" ]]; then

    bcftools query -f '%CHROM\n' "$SOURCE_VCF" | sort -V | uniq > source.chroms.txt

    echo "[INFO] Chromosomes/contigs in selected source VCF:"
    cat source.chroms.txt

    MODE=""

    if grep -Eq '^LG[0-9]+$' source.chroms.txt; then
        MODE="LG"
        grep -E '^LG[0-9]+$' source.chroms.txt > keep.chroms.txt
    elif grep -Eq '^[0-9]+$' source.chroms.txt; then
        MODE="NUMERIC"
        grep -E '^[0-9]+$' source.chroms.txt > keep.chroms.txt
    else
        echo "[ERROR] No LG or numeric chromosomes found in source VCF."
        echo "Your VCF chromosomes are:"
        cat source.chroms.txt
        exit 1
    fi

    echo
    echo "[INFO] Chromosome mode: $MODE"
    echo "[INFO] Chromosomes retained:"
    cat keep.chroms.txt

    CHRS_CSV=$(paste -sd, keep.chroms.txt)

    TMP_KEEP="${OUTDIR}/${PREFIX}.tmp.keepchrom.vcf.gz"
    TMP_CHRNUM="${OUTDIR}/${PREFIX}.tmp.chrnum.vcf.gz"
    CHR_RENAME="${OUTDIR}/${PREFIX}.chr.rename.txt"

    rm -f "$TMP_KEEP" "$TMP_KEEP.tbi" "$TMP_CHRNUM" "$TMP_CHRNUM.tbi" "$FINAL_VCF" "$FINAL_VCF.tbi"

    if [[ "$MODE" == "LG" ]]; then
        echo "[INFO] Keeping only LG chromosomes..."
        bcftools view \
          -r "$CHRS_CSV" \
          -Oz -o "$TMP_KEEP" \
          "$SOURCE_VCF"

        tabix -f -p vcf "$TMP_KEEP"

        echo "[INFO] Renaming LG chromosomes to numeric..."
        awk '{
            old=$1;
            new=old;
            sub(/^LG/, "", new);
            print old"\t"new;
        }' keep.chroms.txt > "$CHR_RENAME"

        bcftools annotate \
          --rename-chrs "$CHR_RENAME" \
          -Oz -o "$TMP_CHRNUM" \
          "$TMP_KEEP"

        tabix -f -p vcf "$TMP_CHRNUM"

    else
        echo "[INFO] Keeping only numeric chromosomes..."
        bcftools view \
          -r "$CHRS_CSV" \
          -Oz -o "$TMP_CHRNUM" \
          "$SOURCE_VCF"

        tabix -f -p vcf "$TMP_CHRNUM"
    fi

    echo "[INFO] Chromosomes after removing non-chromosome contigs:"
    bcftools query -f '%CHROM\n' "$TMP_CHRNUM" | sort -n | uniq

    echo
    echo "[INFO] Checking sample names after chromosome filtering..."
    bcftools query -l "$TMP_CHRNUM" | sort > tmp.vcf.samples.txt
    comm -23 keep.samples.txt tmp.vcf.samples.txt > samples_in_popmap_not_in_vcf.txt

    if [[ -s samples_in_popmap_not_in_vcf.txt ]]; then
        echo "[ERROR] These popmap samples are not in the selected VCF:"
        cat samples_in_popmap_not_in_vcf.txt
        exit 1
    fi

    echo "[INFO] Creating final biallelic SNP VCF for ADMIXTOOLS..."
    bcftools view \
      -S "$KEEP" \
      -m2 -M2 -v snps \
      -Oz -o "$FINAL_VCF" \
      "$TMP_CHRNUM"

    tabix -f -p vcf "$FINAL_VCF"

else
    echo "[SKIP] Final no-Scaffold VCF exists:"
    echo "$FINAL_VCF"

    if [[ ! -s "${FINAL_VCF}.tbi" && ! -s "${FINAL_VCF}.csi" ]]; then
        tabix -f -p vcf "$FINAL_VCF"
    fi
fi

echo
echo "[OK] Final VCF used for ADMIXTOOLS:"
echo "$FINAL_VCF"

echo
echo "[INFO] Final VCF chromosomes:"
bcftools query -f '%CHROM\n' "$FINAL_VCF" | sort -n | uniq

echo
echo "[INFO] Final VCF sample number:"
bcftools query -l "$FINAL_VCF" | wc -l

echo "[INFO] Final VCF SNP number:"
bcftools view -H "$FINAL_VCF" | wc -l

############################################################
# Step 4. VCF to PLINK binary
############################################################

echo
echo "========== Step 4: VCF to PLINK =========="

if [[ "$FORCE" == "1" || ! -s "${PREFIX}.bed" || ! -s "${PREFIX}.bim" || ! -s "${PREFIX}.fam" ]]; then

    rm -f ${PREFIX}.bed ${PREFIX}.bim ${PREFIX}.fam ${PREFIX}.log ${PREFIX}.nosex

    plink \
      --vcf "$FINAL_VCF" \
      --make-bed \
      --out "$PREFIX" \
      --double-id \
      --allow-extra-chr \
      --snps-only just-acgt \
      --biallelic-only strict \
      --noweb

else
    echo "[SKIP] PLINK files already exist:"
    ls -lh ${PREFIX}.bed ${PREFIX}.bim ${PREFIX}.fam
fi

echo "[OK] PLINK files:"
ls -lh ${PREFIX}.bed ${PREFIX}.bim ${PREFIX}.fam

############################################################
# Step 5. Clean BIM for convertf
############################################################

echo
echo "========== Step 5: clean BIM for convertf =========="

awk 'BEGIN{OFS="\t"} {
    $2 = "snp_" NR;
    $3 = 0;
    print $1,$2,$3,$4,$5,$6
}' "${PREFIX}.bim" > "${PREFIX}.bim.clean"

mv "${PREFIX}.bim.clean" "${PREFIX}.bim"

echo "[INFO] First 10 rows of cleaned BIM:"
head "${PREFIX}.bim"

BAD_CHR_COUNT=$(awk '$1 !~ /^[0-9]+$/ {c++} END{print c+0}' "${PREFIX}.bim")

if [[ "$BAD_CHR_COUNT" -gt 0 ]]; then
    echo "[ERROR] Non-numeric chromosome names still exist in BIM:"
    awk '$1 !~ /^[0-9]+$/ {print; if(++n==20) exit}' "${PREFIX}.bim"
    exit 1
fi

echo "[OK] All chromosomes in BIM are numeric."

echo
echo "[INFO] Chromosome counts in BIM:"
awk '{count[$1]++} END{for(c in count) print c, count[c]}' "${PREFIX}.bim" | sort -n

############################################################
# Step 6. convertf PLINK to EIGENSTRAT
############################################################

echo
echo "========== Step 6: convertf PLINK to EIGENSTRAT =========="

if [[ "$FORCE" == "1" || ! -s "${PREFIX}.geno" || ! -s "${PREFIX}.snp" || ! -s "${PREFIX}.ind" ]]; then

    rm -f ${PREFIX}.geno ${PREFIX}.snp ${PREFIX}.ind convertf.log

    cat > convert.par <<EOF
genotypename:    ${PREFIX}.bed
snpname:         ${PREFIX}.bim
indivname:       ${PREFIX}.fam
inputformat:     PACKEDPED
outputformat:    EIGENSTRAT
genotypeoutname: ${PREFIX}.geno
snpoutname:      ${PREFIX}.snp
indivoutname:    ${PREFIX}.ind
familynames:     NO
numchrom:        100
EOF

    echo "[INFO] convert.par:"
    cat convert.par

    set +e
    convertf -p convert.par > convertf.log 2>&1
    CONVERT_STATUS=$?
    set -e

    if [[ "$CONVERT_STATUS" -ne 0 ]]; then
        echo "[WARNING] convertf EIGENSTRAT failed. Showing convertf.log:"
        tail -n 80 convertf.log

        echo
        echo "[INFO] Retrying with PACKEDANCESTRYMAP..."

        rm -f ${PREFIX}.geno ${PREFIX}.snp ${PREFIX}.ind convertf.log

        cat > convert.par <<EOF
genotypename:    ${PREFIX}.bed
snpname:         ${PREFIX}.bim
indivname:       ${PREFIX}.fam
inputformat:     PACKEDPED
outputformat:    PACKEDANCESTRYMAP
genotypeoutname: ${PREFIX}.geno
snpoutname:      ${PREFIX}.snp
indivoutname:    ${PREFIX}.ind
familynames:     NO
numchrom:        100
EOF

        convertf -p convert.par > convertf.log 2>&1
    fi

else
    echo "[SKIP] ADMIXTOOLS files already exist:"
    ls -lh ${PREFIX}.geno ${PREFIX}.snp ${PREFIX}.ind
fi

echo "[OK] ADMIXTOOLS files:"
ls -lh ${PREFIX}.geno ${PREFIX}.snp ${PREFIX}.ind

############################################################
# Step 7. Fix .ind population column
############################################################

echo
echo "========== Step 7: fix .ind population labels =========="

awk '
NR==FNR {
    pop[$1]=$2;
    next
}
{
    id=$1;
    if (id in pop) {
        print $1, $2, pop[id];
    } else {
        print $1, $2, "Ignore";
    }
}
' "$POPMAP2" "${PREFIX}.ind" > "${PREFIX}.ind.new"

echo "[OK] fixed ind file:"
echo "${PREFIX}.ind.new"

echo
echo "[INFO] Population counts in .ind.new:"
awk '{count[$3]++} END{for(p in count) print p, count[p]}' "${PREFIX}.ind.new" | sort

if grep -q "Ignore" "${PREFIX}.ind.new"; then
    echo "[WARNING] Some samples labelled as Ignore:"
    grep "Ignore" "${PREFIX}.ind.new" | head
fi

############################################################
# Step 8. Generate listD
# D(A, B; C, TO)
############################################################

echo
echo "========== Step 8: generate listD =========="

python3 - <<PY
from itertools import combinations

pops = "${LAKE_POPS}".split()
outgroup = "${OUTGROUP}"

with open("listD", "w") as out:
    for C in pops:
        others = [p for p in pops if p != C]
        for A, B in combinations(others, 2):
            out.write(f"{A} {B} {C} {outgroup}\n")

print("[OK] listD generated")
print("Number of D-stat tests:", sum(1 for _ in open("listD")))
PY

echo "[INFO] First 10 rows of listD:"
head listD

############################################################
# Step 9. Run qpDstat
############################################################

echo
echo "========== Step 9: run qpDstat =========="

cat > parD <<EOF
indivname:    ${PREFIX}.ind.new
snpname:      ${PREFIX}.snp
genotypename: ${PREFIX}.geno
popfilename:  listD
printsd:      YES
blgsize:      ${BLGSIZE}
EOF

if [[ "$FORCE" == "1" || ! -s "Dstat.out" ]]; then
    qpDstat -p parD > Dstat.out
else
    echo "[SKIP] Dstat.out already exists."
fi

echo "[OK] qpDstat finished."

############################################################
# Step 10. Extract Dstat table
############################################################

echo
echo "========== Step 10: extract Dstat table =========="

(echo "popA popB popC popO D stderr Z BABA ABBA NSNP"; \
grep "result" Dstat.out | awk '{print $2,$3,$4,$5,$6,$7,$8,$9,$10,$11}') > Dstat.table

awk '
NR==1 {
    print $0, "significance";
    next
}
NR>1 {
    z=$7;
    absz=(z<0 ? -z : z);
    if (absz >= 4) sig="strong_Z4";
    else if (absz >= 3) sig="significant_Z3";
    else sig="ns";
    print $0, sig;
}
' Dstat.table > Dstat.table.significance

awk '
NR==1 {
    print "absZ", $0;
    next
}
NR>1 {
    z=$7;
    absz=(z<0 ? -z : z);
    print absz, $0;
}
' Dstat.table | sort -k1,1gr > Dstat.table.sorted_by_absZ

echo "[OK] Dstat output tables:"
echo "${OUTDIR}/Dstat.table"
echo "${OUTDIR}/Dstat.table.significance"
echo "${OUTDIR}/Dstat.table.sorted_by_absZ"

echo
echo "[INFO] Top 30 D-stat results:"
head -n 31 Dstat.table.sorted_by_absZ

############################################################
# Step 11. Generate listF3 and run qp3Pop
############################################################

echo
echo "========== Step 11: run qp3Pop =========="

python3 - <<PY
from itertools import combinations

pops = "${LAKE_POPS}".split()

with open("listF3", "w") as out:
    for target in pops:
        sources = [p for p in pops if p != target]
        for s1, s2 in combinations(sources, 2):
            out.write(f"{target} {s1} {s2}\n")

print("[OK] listF3 generated")
print("Number of f3 tests:", sum(1 for _ in open("listF3")))
PY

cat > parF3 <<EOF
indivname:    ${PREFIX}.ind.new
snpname:      ${PREFIX}.snp
genotypename: ${PREFIX}.geno
popfilename:  listF3
printsd:      YES
blgsize:      ${BLGSIZE}
EOF

if [[ "$FORCE" == "1" || ! -s "F3.out" ]]; then
    qp3Pop -p parF3 > F3.out
else
    echo "[SKIP] F3.out already exists."
fi

(echo "target source1 source2 f3 stderr Z NSNP"; \
grep "result" F3.out | awk '{print $2,$3,$4,$5,$6,$7,$8}') > F3.table || true

awk '
NR==1 {
    print $0, "significance";
    next
}
NR>1 {
    z=$6;
    f3=$4;
    if (f3 < 0 && z <= -4) sig="strong_negative_Z4";
    else if (f3 < 0 && z <= -3) sig="significant_negative_Z3";
    else sig="ns";
    print $0, sig;
}
' F3.table > F3.table.significance || true

echo "[OK] qp3Pop output:"
echo "${OUTDIR}/F3.table"
echo "${OUTDIR}/F3.table.significance"

############################################################
# Step 12. Optional qpF4ratio
############################################################

echo
echo "========== Step 12: optional qpF4ratio =========="

if [[ -s "f4_models.tsv" && "$HAS_F4RATIO" == "yes" ]]; then

    echo "[INFO] f4_models.tsv found. Running qpF4ratio..."

    awk '
    NF>=5{
        print $1, $2, ":", $3, $4, "::", $1, $2, ":", $5, $4
    }
    ' f4_models.tsv > listF4

    cat > parF4ratio <<EOF
indivname:    ${PREFIX}.ind.new
snpname:      ${PREFIX}.snp
genotypename: ${PREFIX}.geno
popfilename:  listF4
printsd:      YES
blgsize:      ${BLGSIZE}
EOF

    qpF4ratio -p parF4ratio > F4ratio.out

    (echo "RefA Outgroup Target SourceB SourceA alpha stderr Zscore NSNP"; \
    grep "result" F4ratio.out | awk '{print $2,$3,$4,$5,$9,$11,$12,$13,$14}') > F4ratio.table

    echo "[OK] qpF4ratio output:"
    echo "${OUTDIR}/F4ratio.out"
    echo "${OUTDIR}/F4ratio.table"

else
    echo "[SKIP] No f4_models.tsv found or qpF4ratio unavailable."
fi

############################################################
# Step 13. Plot results
############################################################

echo
echo "========== Step 13: plot results =========="

PLOTDIR="${OUTDIR}/plots"
mkdir -p "$PLOTDIR"

if [[ "$HAS_R" == "yes" ]]; then

cat > plot_admixtools_results.R <<'RSCRIPT'
#!/usr/bin/env Rscript

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("ggplot2 is not installed. Please install.packages('ggplot2').")
}

suppressPackageStartupMessages({
  library(ggplot2)
})

plotdir <- "plots"

if (file.exists("Dstat.table")) {

  d <- read.table("Dstat.table", header = TRUE, stringsAsFactors = FALSE)

  d$D <- as.numeric(d$D)
  d$stderr <- as.numeric(d$stderr)
  d$Z <- as.numeric(d$Z)

  d$signif <- ifelse(abs(d$Z) >= 4, "**",
                     ifelse(abs(d$Z) >= 3, "*", ""))

  d$popB <- factor(d$popB, levels = unique(d$popB[order(d$D)]))

  pd <- position_dodge(width = 0.6)

  p <- ggplot(d, aes(x = popB, y = D, color = popC)) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4) +
    geom_point(position = pd, size = 2, alpha = 0.85) +
    geom_errorbar(aes(ymin = D - stderr, ymax = D + stderr),
                  position = pd,
                  width = 0.15,
                  linewidth = 0.35,
                  alpha = 0.7) +
    geom_text(aes(label = signif),
              position = pd,
              vjust = -0.8,
              color = "black",
              size = 4) +
    theme_bw(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 60, hjust = 1),
      panel.grid.minor = element_blank()
    ) +
    labs(
      x = "Population B",
      y = "D-statistic",
      color = "Population C",
      title = "D-statistics among 12 lake populations",
      subtitle = "* |Z| >= 3; ** |Z| >= 4"
    )

  ggsave(file.path(plotdir, "Dstat_plot.pdf"), p, width = 12, height = 7)
  ggsave(file.path(plotdir, "Dstat_plot.png"), p, width = 12, height = 7, dpi = 300)

  d_sig <- d[abs(d$Z) >= 3, ]
  d_sig <- d_sig[order(-abs(d_sig$Z)), ]

  write.table(
    d_sig,
    file.path(plotdir, "Dstat_significant_Z3.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
}

if (file.exists("F3.table")) {

  f3 <- read.table("F3.table", header = TRUE, stringsAsFactors = FALSE)

  if (nrow(f3) > 0) {

    f3$f3 <- as.numeric(f3$f3)
    f3$stderr <- as.numeric(f3$stderr)
    f3$Z <- as.numeric(f3$Z)

    f3$signif <- ifelse(f3$f3 < 0 & f3$Z <= -4, "**",
                        ifelse(f3$f3 < 0 & f3$Z <= -3, "*", ""))

    f3$target <- factor(f3$target, levels = unique(f3$target))

    p3 <- ggplot(f3, aes(x = target, y = f3)) +
      geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4) +
      geom_point(size = 2, alpha = 0.85) +
      geom_errorbar(aes(ymin = f3 - stderr, ymax = f3 + stderr),
                    width = 0.15,
                    linewidth = 0.35,
                    alpha = 0.7) +
      geom_text(aes(label = signif),
                vjust = -0.8,
                color = "black",
                size = 4) +
      theme_bw(base_size = 12) +
      theme(
        axis.text.x = element_text(angle = 60, hjust = 1),
        panel.grid.minor = element_blank()
      ) +
      labs(
        x = "Target population",
        y = "f3 statistic",
        title = "f3 admixture tests among 12 lake populations",
        subtitle = "* f3 < 0 and Z <= -3; ** f3 < 0 and Z <= -4"
      )

    ggsave(file.path(plotdir, "F3_plot.pdf"), p3, width = 12, height = 7)
    ggsave(file.path(plotdir, "F3_plot.png"), p3, width = 12, height = 7, dpi = 300)

    f3_sig <- f3[f3$f3 < 0 & f3$Z <= -3, ]
    f3_sig <- f3_sig[order(f3_sig$Z), ]

    write.table(
      f3_sig,
      file.path(plotdir, "F3_significant_negative_Z3.tsv"),
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )
  }
}

if (file.exists("F4ratio.table")) {

  f4 <- read.table("F4ratio.table", header = TRUE, stringsAsFactors = FALSE)

  if (nrow(f4) > 0) {

    f4$alpha <- as.numeric(f4$alpha)
    f4$stderr <- as.numeric(f4$stderr)
    f4$Zscore <- as.numeric(f4$Zscore)

    f4$Model <- paste(f4$Target, f4$SourceA, f4$SourceB, sep = "_")
    f4$Model <- factor(f4$Model, levels = f4$Model)

    f4$signif <- ifelse(abs(f4$Zscore) >= 4, "**",
                        ifelse(abs(f4$Zscore) >= 3, "*", ""))

    p4 <- ggplot(f4, aes(x = Model, y = alpha)) +
      geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3) +
      geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.3) +
      geom_point(size = 3, alpha = 0.9) +
      geom_errorbar(aes(ymin = alpha - stderr, ymax = alpha + stderr),
                    width = 0.18,
                    linewidth = 0.4) +
      geom_text(aes(label = signif),
                vjust = -0.8,
                color = "black",
                size = 4) +
      theme_bw(base_size = 12) +
      theme(
        axis.text.x = element_text(angle = 60, hjust = 1),
        panel.grid.minor = element_blank()
      ) +
      labs(
        x = "F4-ratio model",
        y = "Admixture proportion",
        title = "F4-ratio ancestry proportion estimates",
        subtitle = "* |Z| >= 3; ** |Z| >= 4"
      )

    ggsave(file.path(plotdir, "F4ratio_plot.pdf"), p4, width = 12, height = 7)
    ggsave(file.path(plotdir, "F4ratio_plot.png"), p4, width = 12, height = 7, dpi = 300)
  }
}
RSCRIPT

    Rscript plot_admixtools_results.R
    echo "[OK] plots saved in: $PLOTDIR"

else
    echo "[SKIP] Rscript not found. Plotting skipped."
fi

############################################################
# Done
############################################################

echo
echo "===================================================="
echo "All done."
echo "Main outputs:"
echo "$FINAL_VCF"
echo "${OUTDIR}/Dstat.table"
echo "${OUTDIR}/Dstat.table.significance"
echo "${OUTDIR}/Dstat.table.sorted_by_absZ"
echo "${OUTDIR}/F3.table"
echo "${OUTDIR}/F3.table.significance"
echo "${OUTDIR}/plots/Dstat_plot.pdf"
echo "${OUTDIR}/plots/Dstat_plot.png"
echo "${OUTDIR}/plots/F3_plot.pdf"
echo "${OUTDIR}/plots/F3_plot.png"
echo "===================================================="
