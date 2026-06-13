 #!/usr/bin/env bash
set -euo pipefail

#############################################
# ADMIXTOOLS / qpDstat workflow for fish
# 12 lake populations + TO outgroup
#
# Analyses:
# 1. VCF -> PLINK -> EIGENSTRAT
# 2. qpDstat: D(A, B; C, TO)
# 3. qp3Pop: f3(Target; Source1, Source2)
# 4. Optional qpF4ratio if f4_models.tsv is provided
# 5. Plot Dstat and f3 results
#############################################

#############################################
# 1. 你已经修改好的路径
#############################################

# 输入 VCF：不需要 LD pruning
VCF="/media/server/Elements/BC2026040063-BN260407NJ01S02N5-134个异尾高原鳅重测序数据变异分析/gvcf_joint/snapp_filter_missing075/fish.02_missing075.DP10.MAF005.vcf.gz"

# 你的三列表格，格式：
# sample_id    sample_id    population
POPMAP3="/media/server/Elements/BC2026040063-BN260407NJ01S02N5-134个异尾高原鳅重测序数据变异分析/gvcf_joint/snapp_filter_missing075/plink_result/population_table.txt"

# 输出目录
# 注意：目录名字叫 Dsuite 没关系，但这个脚本运行的是 qpDstat / qp3Pop / qpF4ratio，不是 Dsuite
OUTDIR="/media/server/Elements/BC2026040063-BN260407NJ01S02N5-134个异尾高原鳅重测序数据变异分析/gvcf_joint/Dsuite"

#############################################
# 2. 分析参数
#############################################

# 12 个湖泊群体
LAKE_POPS="BC CEND CN GN GRC MJC PC QCG SLC YC YQC ZGC"

# 外群
OUTGROUP="TO"

# ADMIXTOOLS block size
# 练习中是 0.01；WGS 数据先用 0.01 可以，后续可测试 0.05
BLGSIZE=0.01

# 输出前缀
PREFIX="fish_admix"

#############################################
# 3. 不需要修改
#############################################

mkdir -p "$OUTDIR"
cd "$OUTDIR"

echo
echo "============================================="
echo "ADMIXTOOLS / qpDstat workflow"
echo "VCF: $VCF"
echo "POPMAP3: $POPMAP3"
echo "OUTDIR: $OUTDIR"
echo "Lake populations: $LAKE_POPS"
echo "Outgroup: $OUTGROUP"
echo "Block size: $BLGSIZE"
echo "============================================="
echo

#############################################
# Step 0. Check tools
#############################################

echo "========== Step 0: check files and tools =========="

if [[ ! -s "$VCF" ]]; then
    echo "[ERROR] VCF not found:"
    echo "$VCF"
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

echo "[OK] all required tools found."

#############################################
# Step 1. Rename chromosomes
# LG1 -> 1, LG2 -> 2
#############################################

echo
echo "========== Step 1: rename chromosomes =========="

VCF_RENAME="${OUTDIR}/fish.02_missing075.DP10.MAF005.chrnum.vcf.gz"
CHR_MAP="${OUTDIR}/chr.rename.txt"

if [[ ! -s "${VCF}.tbi" && ! -s "${VCF}.csi" ]]; then
    echo "[WARNING] VCF index not found. Building index..."
    tabix -f -p vcf "$VCF"
fi

bcftools query -f '%CHROM\n' "$VCF" | sort -V | uniq | \
awk '
{
    old=$1;
    new=old;

    if (new ~ /^LG[0-9]+$/) {
        sub(/^LG/, "", new);
    } else if (new ~ /^[Cc]hr[0-9]+$/) {
        sub(/^[Cc]hr/, "", new);
        new = new + 0;
    }

    print old"\t"new;
}' > "$CHR_MAP"

echo "Chromosome rename table:"
cat "$CHR_MAP"

if [[ ! -s "$VCF_RENAME" ]]; then
    bcftools annotate \
      --rename-chrs "$CHR_MAP" \
      -Oz -o "$VCF_RENAME" \
      "$VCF"

    tabix -f -p vcf "$VCF_RENAME"
else
    echo "[SKIP] renamed VCF already exists:"
    echo "$VCF_RENAME"
fi

echo "[OK] renamed VCF:"
echo "$VCF_RENAME"

#############################################
# Step 2. Make two-column popmap and keep list
#############################################

echo
echo "========== Step 2: make popmap and keep list =========="

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

echo "[OK] two-column popmap:"
echo "$POPMAP2"

echo
echo "Population counts:"
awk '{count[$2]++} END{for(p in count) print p, count[p]}' "$POPMAP2" | sort

#############################################
# Step 3. Check sample names
#############################################

echo
echo "========== Step 3: check sample names =========="

bcftools query -l "$VCF_RENAME" | sort > vcf.samples.txt
sort "$KEEP" > keep.samples.txt

comm -23 keep.samples.txt vcf.samples.txt > samples_in_popmap_not_in_vcf.txt
comm -13 keep.samples.txt vcf.samples.txt > samples_in_vcf_not_in_popmap.txt

if [[ -s samples_in_popmap_not_in_vcf.txt ]]; then
    echo "[ERROR] 以下 popmap 样本不在 VCF 中："
    cat samples_in_popmap_not_in_vcf.txt
    exit 1
fi

echo "[OK] all popmap samples are present in VCF."

if [[ -s samples_in_vcf_not_in_popmap.txt ]]; then
    echo "[WARNING] 以下 VCF 样本不在本次分析中，将被剔除："
    cat samples_in_vcf_not_in_popmap.txt
fi

#############################################
# Step 4. Prepare VCF
# 不做 LD pruning
#############################################

echo
echo "========== Step 4: prepare VCF for ADMIXTOOLS =========="

VCF_ADMIX="${OUTDIR}/${PREFIX}.biallelicSNP.vcf.gz"

if [[ ! -s "$VCF_ADMIX" ]]; then
    bcftools view \
      -S "$KEEP" \
      -m2 -M2 -v snps \
      -Oz -o "$VCF_ADMIX" \
      "$VCF_RENAME"

    tabix -f -p vcf "$VCF_ADMIX"
else
    echo "[SKIP] ADMIXTOOLS VCF already exists:"
    echo "$VCF_ADMIX"
fi

echo "Sample number:"
bcftools query -l "$VCF_ADMIX" | wc -l

echo "SNP number:"
bcftools view -H "$VCF_ADMIX" | wc -l

#############################################
# Step 5. VCF to PLINK binary
#############################################

echo
echo "========== Step 5: VCF to PLINK =========="

plink \
  --vcf "$VCF_ADMIX" \
  --make-bed \
  --out "$PREFIX" \
  --double-id \
  --allow-extra-chr \
  --snps-only just-acgt \
  --biallelic-only strict \
  --noweb

# 修正 .bim 中 SNP ID，避免 . 或重复 ID 影响 convertf
awk 'BEGIN{OFS="\t"} {$2=$1":"$4":"$5":"$6":"NR; print}' "${PREFIX}.bim" > "${PREFIX}.bim.tmp"
mv "${PREFIX}.bim.tmp" "${PREFIX}.bim"

echo "[OK] PLINK files generated:"
ls -lh ${PREFIX}.bed ${PREFIX}.bim ${PREFIX}.fam

#############################################
# Step 6. PLINK to EIGENSTRAT
#############################################

echo
echo "========== Step 6: convertf PLINK to EIGENSTRAT =========="

cat > convert.par <<EOF
genotypename:    ${PREFIX}.bed
snpname:         ${PREFIX}.bim
indivname:       ${PREFIX}.fam
outputformat:    EIGENSTRAT
genotypeoutname: ${PREFIX}.geno
snpoutname:      ${PREFIX}.snp
indivoutname:    ${PREFIX}.ind
familynames:     NO
EOF

convertf -p convert.par > convertf.log 2>&1

echo "[OK] EIGENSTRAT files:"
ls -lh ${PREFIX}.geno ${PREFIX}.snp ${PREFIX}.ind

#############################################
# Step 7. Fix .ind population column
#############################################

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
echo "Population counts in .ind.new:"
awk '{count[$3]++} END{for(p in count) print p, count[p]}' "${PREFIX}.ind.new" | sort

if grep -q "Ignore" "${PREFIX}.ind.new"; then
    echo "[WARNING] Some samples labelled as Ignore:"
    grep "Ignore" "${PREFIX}.ind.new" | head
fi

#############################################
# Step 8. Generate listD
# D(A, B; C, TO)
#############################################

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

echo "First 10 listD rows:"
head listD

#############################################
# Step 9. Run qpDstat
#############################################

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

qpDstat -p parD > Dstat.out

echo "[OK] qpDstat finished."

#############################################
# Step 10. Extract Dstat table
#############################################

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
echo "Top 30 D-stat results:"
head -n 31 Dstat.table.sorted_by_absZ

#############################################
# Step 11. Generate listF3
# f3(Target; Source1, Source2)
#############################################

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

qp3Pop -p parF3 > F3.out

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

#############################################
# Step 12. Optional qpF4ratio
# 只有存在 f4_models.tsv 时才运行
#
# f4_models.tsv 五列：
# RefA  Outgroup  Target  SourceB  SourceA
#
# 自动转换为：
# RefA Outgroup : Target SourceB :: RefA Outgroup : SourceA SourceB
#############################################

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
    echo "如果后续要做 f4-ratio，请在 OUTDIR 下新建 f4_models.tsv。"
fi

#############################################
# Step 13. Plot Dstat, F3, and optional F4ratio
#############################################

echo
echo "========== Step 13: plot results =========="

PLOTDIR="${OUTDIR}/plots"
mkdir -p "$PLOTDIR"

cat > plot_admixtools_results.R <<'RSCRIPT'
#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
})

plotdir <- "plots"

#############################################
# Plot Dstat
#############################################

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

#############################################
# Plot f3
#############################################

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

#############################################
# Plot optional F4ratio
#############################################

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

echo "[OK] Plots saved in:"
echo "${PLOTDIR}"

echo
echo "========== All done =========="
echo "Main outputs:"
echo "${OUTDIR}/Dstat.table"
echo "${OUTDIR}/Dstat.table.significance"
echo "${OUTDIR}/Dstat.table.sorted_by_absZ"
echo "${OUTDIR}/F3.table"
echo "${OUTDIR}/F3.table.significance"
echo "${PLOTDIR}/Dstat_plot.pdf"
echo "${PLOTDIR}/Dstat_plot.png"
echo "${PLOTDIR}/F3_plot.pdf"
echo "${PLOTDIR}/F3_plot.png"
