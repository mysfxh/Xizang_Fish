#!/usr/bin/env bash
set -euo pipefail

############################################################
# Resume f4-statistic + f4-ratio + plotting
#
# 功能：
# 1. 不重跑已有 VCF / PLINK / convertf / D-statistics
# 2. 基于已有 ADMIXTOOLS 文件跑 f4 statistic
# 3. 汇总 f4 statistic 为 66 个湖泊对
# 4. 如果存在 f4_models.tsv，则跑 qpF4ratio
# 5. 自动绘制 f4-ratio 图和 f4 statistic 图
############################################################

OUTDIR="/media/server/Elements/BC2026040063-BN260407NJ01S02N5-134个异尾高原鳅重测序数据变异分析/gvcf_joint/Dsuite"

PREFIX="fish_admix_noScaffold"

LAKE_POPS="BC CEND CN GN GRC MJC PC QCG SLC YC YQC ZGC"
OUTGROUP="TO"

BLGSIZE=0.01

# 默认不覆盖已有结果
FORCE="${FORCE:-0}"

# 如果你修改了 f4_models.tsv，想重新跑 f4-ratio：
# FORCE_F4RATIO=1 bash run_f4_all_resume.sh
FORCE_F4RATIO="${FORCE_F4RATIO:-0}"

cd "$OUTDIR"

echo
echo "===================================================="
echo "Run f4-statistic + f4-ratio + plotting"
echo "OUTDIR:         $OUTDIR"
echo "PREFIX:         $PREFIX"
echo "FORCE:          $FORCE"
echo "FORCE_F4RATIO:  $FORCE_F4RATIO"
echo "===================================================="
echo

############################################################
# Step 0. Check tools and files
############################################################

echo "========== Step 0: check tools and files =========="

command -v qpDstat >/dev/null 2>&1 || {
    echo "[ERROR] qpDstat not found. Please activate ADMIXTOOLS environment."
    exit 1
}

if command -v qpF4ratio >/dev/null 2>&1; then
    HAS_QPF4RATIO="yes"
else
    HAS_QPF4RATIO="no"
    echo "[WARNING] qpF4ratio not found. f4-ratio step will be skipped."
fi

if command -v Rscript >/dev/null 2>&1; then
    HAS_R="yes"
else
    HAS_R="no"
    echo "[WARNING] Rscript not found. Plotting will be skipped."
fi

for f in "${PREFIX}.geno" "${PREFIX}.snp" "${PREFIX}.ind.new"; do
    if [[ ! -s "$f" ]]; then
        echo "[ERROR] Required ADMIXTOOLS file missing:"
        echo "$OUTDIR/$f"
        exit 1
    fi
done

echo "[OK] Required ADMIXTOOLS files found."

############################################################
# Step 1. Generate listD for f4 statistic
############################################################

echo
echo "========== Step 1: generate listD =========="

if [[ "$FORCE" == "1" || ! -s "listD" ]]; then

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
print("Number of f4 tests:", sum(1 for _ in open("listD")))
PY

else
    echo "[SKIP] listD already exists."
    echo "Number of f4 tests: $(wc -l < listD)"
fi

############################################################
# Step 2. Run f4 statistic using qpDstat f4mode
############################################################

echo
echo "========== Step 2: run f4 statistic =========="

cat > parF4stat <<EOF
indivname:    ${PREFIX}.ind.new
snpname:      ${PREFIX}.snp
genotypename: ${PREFIX}.geno
popfilename:  listD
printsd:      YES
blgsize:      ${BLGSIZE}
f4mode:       YES
EOF

if [[ "$FORCE" == "1" || ! -s "F4stat.out" ]]; then
    qpDstat -p parF4stat > F4stat.out
else
    echo "[SKIP] F4stat.out already exists."
fi

echo "[OK] f4 statistic finished."

############################################################
# Step 3. Extract f4 statistic table
############################################################

echo
echo "========== Step 3: extract f4 statistic table =========="

if [[ "$FORCE" == "1" || ! -s "F4stat.table" ]]; then

python3 - <<'PY'
out = open("F4stat.table", "w")
out.write("popA\tpopB\tpopC\tpopO\tf4\tstderr\tZ\tNSNP\n")

with open("F4stat.out") as f:
    for line in f:
        if not line.startswith("result:"):
            continue

        x = line.strip().split()

        popA = x[1]
        popB = x[2]
        popC = x[3]
        popO = x[4]

        f4 = x[5] if len(x) > 5 else "NA"
        stderr = x[6] if len(x) > 6 else "NA"
        Z = x[7] if len(x) > 7 else "NA"
        NSNP = x[8] if len(x) > 8 else "NA"

        out.write(f"{popA}\t{popB}\t{popC}\t{popO}\t{f4}\t{stderr}\t{Z}\t{NSNP}\n")

out.close()
PY

else
    echo "[SKIP] F4stat.table already exists."
fi

if [[ "$FORCE" == "1" || ! -s "F4stat.table.significance" ]]; then

awk '
BEGIN{OFS="\t"}
NR==1 {
    print $0, "absZ", "significance";
    next
}
NR>1 {
    z=$7;
    absz=(z<0 ? -z : z);

    if (absz >= 4) sig="strong_Z4";
    else if (absz >= 3) sig="significant_Z3";
    else sig="ns";

    print $0, absz, sig;
}
' F4stat.table > F4stat.table.significance

else
    echo "[SKIP] F4stat.table.significance already exists."
fi

if [[ "$FORCE" == "1" || ! -s "F4stat.table.sorted_by_absZ" ]]; then
    {
        head -n 1 F4stat.table.significance
        awk 'NR>1' F4stat.table.significance | sort -k9,9gr
    } > F4stat.table.sorted_by_absZ
else
    echo "[SKIP] F4stat.table.sorted_by_absZ already exists."
fi

############################################################
# Step 4. Summarize f4 statistic into 66 lake pairs
############################################################

echo
echo "========== Step 4: summarize f4 statistic into 66 lake pairs =========="

if [[ "$FORCE" == "1" || ! -s "F4stat_66_lake_pair_summary.tsv" ]]; then

python3 - <<'PY'
import pandas as pd
import numpy as np
from itertools import combinations

infile = "F4stat.table.significance"

pops = ["BC", "CEND", "CN", "GN", "GRC", "MJC", "PC", "QCG", "SLC", "YC", "YQC", "ZGC"]

df = pd.read_csv(infile, sep="\t")

for col in ["f4", "stderr", "Z", "NSNP"]:
    df[col] = pd.to_numeric(df[col], errors="coerce")

df["absF4"] = df["f4"].abs()
df["absZ"] = df["Z"].abs()

sig = df[df["absZ"] >= 3].copy()

# f4(A, B; C, O)
# f4 > 0: A 与 C 共享更多
# f4 < 0: B 与 C 共享更多
sig["shared_pop1"] = np.where(sig["f4"] > 0, sig["popA"], sig["popB"])
sig["shared_pop2"] = sig["popC"]

sig["candidate_pair"] = sig.apply(
    lambda x: "-".join(sorted([x["shared_pop1"], x["shared_pop2"]])),
    axis=1
)

sig_out = sig[
    [
        "popA", "popB", "popC", "popO",
        "f4", "stderr", "Z", "absZ", "absF4",
        "NSNP", "significance",
        "shared_pop1", "shared_pop2", "candidate_pair"
    ]
].copy()

sig_out = sig_out.sort_values(["candidate_pair", "absZ"], ascending=[True, False])
sig_out.to_csv("F4stat_significant_rows_inferred_pairs.tsv", sep="\t", index=False)

summary = (
    sig.groupby("candidate_pair")
    .agg(
        n_supported_tests=("candidate_pair", "size"),
        n_strong_Z4=("significance", lambda x: (x == "strong_Z4").sum()),
        mean_abs_f4=("absF4", "mean"),
        max_abs_f4=("absF4", "max"),
        mean_abs_Z=("absZ", "mean"),
        max_abs_Z=("absZ", "max"),
        mean_NSNP=("NSNP", "mean")
    )
    .reset_index()
)

all_pairs = pd.DataFrame({
    "candidate_pair": ["-".join(sorted(x)) for x in combinations(pops, 2)]
})

summary = all_pairs.merge(summary, on="candidate_pair", how="left")

for col in [
    "n_supported_tests", "n_strong_Z4",
    "mean_abs_f4", "max_abs_f4",
    "mean_abs_Z", "max_abs_Z", "mean_NSNP"
]:
    summary[col] = summary[col].fillna(0)

summary["support_rate"] = summary["n_supported_tests"] / 20

def classify(row):
    n = row["n_supported_tests"]
    if n >= 10:
        return "robust"
    elif n >= 5:
        return "moderate"
    elif n >= 1:
        return "weak"
    else:
        return "no_evidence"

summary["evidence_level"] = summary.apply(classify, axis=1)

order = {"robust": 0, "moderate": 1, "weak": 2, "no_evidence": 3}
summary["rank_order"] = summary["evidence_level"].map(order)

summary = summary.sort_values(
    ["rank_order", "n_supported_tests", "mean_abs_Z", "mean_abs_f4"],
    ascending=[True, False, False, False]
).drop(columns=["rank_order"])

summary.to_csv("F4stat_66_lake_pair_summary.tsv", sep="\t", index=False)

print("Done.")
print("Top 30 f4 candidate pairs:")
print(summary.head(30).to_string(index=False))
PY

else
    echo "[SKIP] F4stat_66_lake_pair_summary.tsv already exists."
fi

############################################################
# Step 5. Prepare f4-ratio model file
############################################################

echo
echo "========== Step 5: prepare f4-ratio model file =========="

if [[ ! -s "f4_models.tsv" ]]; then

cat > f4_models.template.tsv <<EOF
# Format:
# RefA    Outgroup    Target    SourceB    SourceA    group
#
# Example:
# BC      TO          CN        GN         CEND       directE
# BC      TO          GN        CN         CEND       directE
# CEND    TO          BC        PC         CN         directW
# CEND    TO          PC        BC         CN         directW
# BC      TO          YC        YQC        GRC        directW
# BC      TO          YQC       YC         GRC        directW
# BC      TO          CEND      MJC        CN         indirectW
# BC      TO          MJC       CEND       CN         indirectW
# BC      TO          QCG       ZGC        CN         indirectW
# BC      TO          ZGC       QCG        CN         indirectW
# BC      TO          SLC       CN         GN         indirectW
# BC      TO          SLC       GN         CN         indirectW
EOF

    echo "[WARNING] f4_models.tsv not found."
    echo "[INFO] Template created:"
    echo "${OUTDIR}/f4_models.template.tsv"
    echo
    echo "Please edit f4_models.template.tsv and save as f4_models.tsv, then rerun this script."
else
    echo "[OK] f4_models.tsv found."
fi

############################################################
# Step 6. Run qpF4ratio
############################################################

echo
echo "========== Step 6: run qpF4ratio =========="

if [[ "$HAS_QPF4RATIO" != "yes" ]]; then
    echo "[SKIP] qpF4ratio not found."
elif [[ ! -s "f4_models.tsv" ]]; then
    echo "[SKIP] f4_models.tsv not found. f4-ratio not run."
else

    echo "[INFO] Cleaning f4_models.tsv..."

    awk '
    BEGIN{OFS="\t"}
    NF>=5 && $1 !~ /^#/ {
        group="f4ratio";
        if (NF>=6) group=$6;
        print $1,$2,$3,$4,$5,group;
    }
    ' f4_models.tsv > f4_models.clean.tsv

    echo "[INFO] Cleaned models:"
    cat f4_models.clean.tsv

    if [[ ! -s "f4_models.clean.tsv" ]]; then
        echo "[ERROR] f4_models.clean.tsv is empty. Please check f4_models.tsv."
        exit 1
    fi

    awk '
    BEGIN{OFS=" "}
    NF>=6 {
        print $1, $2, ":", $3, $4, "::", $1, $2, ":", $5, $4
    }
    ' f4_models.clean.tsv > listF4

    echo "[INFO] Converted listF4:"
    cat listF4

    cat > parF4ratio <<EOF
indivname:    ${PREFIX}.ind.new
snpname:      ${PREFIX}.snp
genotypename: ${PREFIX}.geno
popfilename:  listF4
printsd:      YES
blgsize:      ${BLGSIZE}
EOF

    if [[ "$FORCE" == "1" || "$FORCE_F4RATIO" == "1" || ! -s "F4ratio.out" ]]; then
        qpF4ratio -p parF4ratio > F4ratio.out
    else
        echo "[SKIP] F4ratio.out already exists."
        echo "       To rerun after editing f4_models.tsv, use:"
        echo "       FORCE_F4RATIO=1 bash run_f4_all_resume.sh"
    fi

    echo "[OK] qpF4ratio finished or skipped."

    ########################################################
    # Extract qpF4ratio table
    ########################################################

    if [[ "$FORCE" == "1" || "$FORCE_F4RATIO" == "1" || ! -s "F4ratio.table" ]]; then

python3 - <<'PY'
import re

models = []
with open("f4_models.clean.tsv") as f:
    for line in f:
        line = line.strip()
        if not line:
            continue
        x = line.split()
        if len(x) >= 6:
            models.append(x[:6])

result_lines = []
with open("F4ratio.out") as f:
    for line in f:
        if line.startswith("result:"):
            result_lines.append(line.strip())

out = open("F4ratio.table", "w")
out.write("RefA\tOutgroup\tTarget\tSourceB\tSourceA\tgroup\talpha\tstderr\tZscore\tNSNP\n")

for i, line in enumerate(result_lines):
    if i < len(models):
        RefA, Outgroup, Target, SourceB, SourceA, group = models[i]
    else:
        RefA = Outgroup = Target = SourceB = SourceA = group = "NA"

    x = line.split()

    numeric = []
    for item in x:
        try:
            float(item)
            numeric.append(item)
        except:
            pass

    alpha = "NA"
    stderr = "NA"
    zscore = "NA"
    nsnp = "NA"

    # 常见 qpF4ratio 输出末尾包含 alpha stderr Zscore NSNP
    if len(numeric) >= 4:
        alpha = numeric[-4]
        stderr = numeric[-3]
        zscore = numeric[-2]
        nsnp = numeric[-1]
    elif len(numeric) >= 3:
        alpha = numeric[-3]
        stderr = numeric[-2]
        zscore = numeric[-1]

    out.write(f"{RefA}\t{Outgroup}\t{Target}\t{SourceB}\t{SourceA}\t{group}\t{alpha}\t{stderr}\t{zscore}\t{nsnp}\n")

out.close()
print(f"[OK] Extracted {len(result_lines)} qpF4ratio result lines.")
PY

    else
        echo "[SKIP] F4ratio.table already exists."
    fi

    if [[ -s "F4ratio.table" ]]; then

awk '
BEGIN{OFS="\t"}
NR==1 {
    print $0, "alpha_valid", "absZ", "significance";
    next
}
NR>1 {
    alpha=$7;
    z=$9;
    absz=(z<0 ? -z : z);

    if (alpha >= 0 && alpha <= 1) valid="yes";
    else valid="no";

    if (absz >= 4) sig="strong_Z4";
    else if (absz >= 3) sig="significant_Z3";
    else sig="ns";

    print $0, valid, absz, sig;
}
' F4ratio.table > F4ratio.table.significance

        echo "[OK] F4ratio.table.significance generated."
    fi
fi

############################################################
# Step 7. Plot f4 statistic and f4-ratio
############################################################

echo
echo "========== Step 7: plotting =========="

PLOTDIR="${OUTDIR}/plots"
mkdir -p "$PLOTDIR"

if [[ "$HAS_R" == "yes" ]]; then

cat > plot_f4_all.R <<'RSCRIPT'
#!/usr/bin/env Rscript

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  message("ggplot2 not installed. Skip plotting.")
  quit(save = "no", status = 0)
}

suppressPackageStartupMessages(library(ggplot2))

plotdir <- "plots"
dir.create(plotdir, showWarnings = FALSE)

############################################################
# 1. Plot f4 statistic 66-pair summary
############################################################

if (file.exists("F4stat_66_lake_pair_summary.tsv")) {

  df <- read.table(
    "F4stat_66_lake_pair_summary.tsv",
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE
  )

  df <- df[order(df$support_rate, df$mean_abs_f4, decreasing = TRUE), ]
  top <- head(df, 30)

  top$candidate_pair <- factor(top$candidate_pair, levels = rev(top$candidate_pair))

  p1 <- ggplot(top, aes(x = candidate_pair, y = support_rate, fill = evidence_level)) +
    geom_col(width = 0.75) +
    coord_flip() +
    theme_bw(base_size = 12) +
    theme(panel.grid.minor = element_blank()) +
    labs(
      x = "Candidate lake pair",
      y = "Support rate",
      fill = "Evidence level",
      title = "Top candidate gene-flow pairs based on f4-statistics",
      subtitle = "Support rate = number of significant supporting tests / 20"
    )

  ggsave(file.path(plotdir, "F4stat_66_pair_support_rate_top30.pdf"),
         p1, width = 9, height = 8)
  ggsave(file.path(plotdir, "F4stat_66_pair_support_rate_top30.png"),
         p1, width = 9, height = 8, dpi = 300)

  p2 <- ggplot(top, aes(x = candidate_pair, y = mean_abs_f4, fill = evidence_level)) +
    geom_col(width = 0.75) +
    coord_flip() +
    theme_bw(base_size = 12) +
    theme(panel.grid.minor = element_blank()) +
    labs(
      x = "Candidate lake pair",
      y = "Mean |f4|",
      fill = "Evidence level",
      title = "Magnitude of f4-statistic",
      subtitle = "Top 30 candidate lake pairs"
    )

  ggsave(file.path(plotdir, "F4stat_66_pair_mean_absF4_top30.pdf"),
         p2, width = 9, height = 8)
  ggsave(file.path(plotdir, "F4stat_66_pair_mean_absF4_top30.png"),
         p2, width = 9, height = 8, dpi = 300)
}

############################################################
# 2. Plot f4-ratio results
############################################################

if (file.exists("F4ratio.table.significance")) {

  d <- read.table(
    "F4ratio.table.significance",
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  if (nrow(d) > 0) {

    d$alpha <- as.numeric(d$alpha)
    d$stderr <- as.numeric(d$stderr)
    d$Zscore <- as.numeric(d$Zscore)
    d$absZ <- as.numeric(d$absZ)

    if (!"group" %in% colnames(d)) {
      d$group <- "f4ratio"
    }

    d$model_id <- seq_len(nrow(d))

    # 你原来的逻辑：
    # indirectW 用 1 - alpha 统一方向
    d$alpha_plot <- ifelse(d$group == "indirectW", 1 - d$alpha, d$alpha)

    d$ymin <- d$alpha_plot - d$stderr
    d$ymax <- d$alpha_plot + d$stderr

    d$source_label <- ifelse(
      d$group == "indirectW",
      paste0("B: ", d$SourceB),
      paste0("A: ", d$SourceA)
    )

    d$valid_label <- ifelse(d$alpha_valid == "yes", "", "invalid")

    d$model_label <- paste0(
      d$Target, "\n",
      d$SourceA, "/", d$SourceB
    )

    d$model_label <- factor(d$model_label, levels = unique(d$model_label))

    p3 <- ggplot(d, aes(x = model_label, y = alpha_plot)) +
      geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.35) +
      geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.35) +
      geom_point(size = 3) +
      geom_errorbar(aes(ymin = ymin, ymax = ymax), width = 0.18, linewidth = 0.4) +
      geom_text(aes(label = source_label), vjust = -0.8, size = 3) +
      geom_text(aes(label = valid_label), vjust = 1.8, size = 3) +
      facet_wrap(~ group, scales = "free_x") +
      theme_bw(base_size = 12) +
      theme(
        axis.text.x = element_text(angle = 60, hjust = 1),
        panel.grid.minor = element_blank()
      ) +
      labs(
        x = "Target and source model",
        y = "Admixture proportion",
        title = "f4-ratio admixture proportion estimates",
        subtitle = "For indirectW models, plotted value is 1 - alpha"
      )

    ggsave(file.path(plotdir, "F4ratio_admixture_proportion.pdf"),
           p3, width = 11, height = 7)
    ggsave(file.path(plotdir, "F4ratio_admixture_proportion.png"),
           p3, width = 11, height = 7, dpi = 300)

    p4 <- ggplot(d, aes(x = model_label, y = Zscore)) +
      geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.35) +
      geom_hline(yintercept = 3, linetype = "dotted", linewidth = 0.35) +
      geom_hline(yintercept = -3, linetype = "dotted", linewidth = 0.35) +
      geom_point(size = 3) +
      facet_wrap(~ group, scales = "free_x") +
      theme_bw(base_size = 12) +
      theme(
        axis.text.x = element_text(angle = 60, hjust = 1),
        panel.grid.minor = element_blank()
      ) +
      labs(
        x = "Target and source model",
        y = "Z-score",
        title = "Z-scores of f4-ratio models"
      )

    ggsave(file.path(plotdir, "F4ratio_Zscore.pdf"),
           p4, width = 11, height = 7)
    ggsave(file.path(plotdir, "F4ratio_Zscore.png"),
           p4, width = 11, height = 7, dpi = 300)
  }
}
RSCRIPT

    Rscript plot_f4_all.R
    echo "[OK] Plots saved in: $PLOTDIR"
else
    echo "[SKIP] Rscript not available."
fi

############################################################
# Done
############################################################

echo
echo "===================================================="
echo "All done."
echo
echo "f4 statistic outputs:"
echo "1. ${OUTDIR}/F4stat.out"
echo "2. ${OUTDIR}/F4stat.table"
echo "3. ${OUTDIR}/F4stat.table.significance"
echo "4. ${OUTDIR}/F4stat_66_lake_pair_summary.tsv"
echo "5. ${OUTDIR}/F4stat_significant_rows_inferred_pairs.tsv"
echo
echo "f4-ratio outputs:"
echo "6. ${OUTDIR}/f4_models.tsv"
echo "7. ${OUTDIR}/listF4"
echo "8. ${OUTDIR}/F4ratio.out"
echo "9. ${OUTDIR}/F4ratio.table"
echo "10. ${OUTDIR}/F4ratio.table.significance"
echo
echo "plots:"
echo "11. ${OUTDIR}/plots/F4stat_66_pair_support_rate_top30.pdf"
echo "12. ${OUTDIR}/plots/F4stat_66_pair_mean_absF4_top30.pdf"
echo "13. ${OUTDIR}/plots/F4ratio_admixture_proportion.pdf"
echo "14. ${OUTDIR}/plots/F4ratio_Zscore.pdf"
echo "===================================================="
