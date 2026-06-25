#!/usr/bin/env bash
set -euo pipefail

############################################################
# Outgroup-f3 workflow for paper-style heatmap
# Compute f3(X, Y; OUTGROUP) using qp3Pop
############################################################

OUTDIR="/media/server/Elements/BC2026040063-BN260407NJ01S02N5-134个异尾高原鳅重测序数据变异分析/gvcf_joint/Dsuite"
PREFIX="fish_admix_noScaffold"

LAKE_POPS="BC CEND CN GN GRC MJC PC QCG SLC YC YQC ZGC"
OUTGROUP="TO"

BLGSIZE=0.01

cd "$OUTDIR"

echo "========== Outgroup-f3 analysis =========="
echo "OUTDIR:   $OUTDIR"
echo "PREFIX:   $PREFIX"
echo "Pops:     $LAKE_POPS"
echo "Outgroup: $OUTGROUP"

for f in "${PREFIX}.geno" "${PREFIX}.snp" "${PREFIX}.ind.new"; do
    if [[ ! -s "$f" ]]; then
        echo "[ERROR] Missing file: $f"
        exit 1
    fi
done

command -v qp3Pop >/dev/null 2>&1 || {
    echo "[ERROR] qp3Pop not found. Please activate ADMIXTOOLS environment."
    exit 1
}

mkdir -p outgroup_f3 plots_outgroup_f3

############################################################
# 1. Generate outgroup-f3 list
# qp3Pop line: OUTGROUP X Y
# This corresponds to outgroup-f3 f3(X, Y; OUTGROUP)
############################################################

python3 - <<PY
from itertools import combinations

pops = "${LAKE_POPS}".split()
outgroup = "${OUTGROUP}"

with open("outgroup_f3/listF3out", "w") as out:
    for x, y in combinations(pops, 2):
        out.write(f"{outgroup} {x} {y}\n")

print("[OK] listF3out generated")
print("Number of outgroup-f3 tests:", sum(1 for _ in open("outgroup_f3/listF3out")))
PY

head outgroup_f3/listF3out

############################################################
# 2. Run qp3Pop
############################################################

cat > outgroup_f3/parF3out <<EOF
indivname:    ${PREFIX}.ind.new
snpname:      ${PREFIX}.snp
genotypename: ${PREFIX}.geno
popfilename:  outgroup_f3/listF3out
printsd:      YES
blgsize:      ${BLGSIZE}
EOF

echo "[INFO] Running qp3Pop for outgroup-f3..."

qp3Pop -p outgroup_f3/parF3out > outgroup_f3/F3out.out

############################################################
# 3. Extract table
############################################################

# qp3Pop output usually:
# result: pop1 pop2 pop3 f3 stderr Z SNPs
# Here pop1 = outgroup, pop2 = X, pop3 = Y
# We rename to X Y outgroup f3 stderr Z NSNP

(echo "X Y outgroup f3 stderr Z NSNP"; \
grep "result" outgroup_f3/F3out.out | awk '{print $3,$4,$2,$5,$6,$7,$8}') \
> outgroup_f3/OutgroupF3.table

echo "[OK] Outgroup-f3 table:"
echo "${OUTDIR}/outgroup_f3/OutgroupF3.table"

head outgroup_f3/OutgroupF3.table

############################################################
# 4. Plot heatmap and clustering tree
############################################################

cat > outgroup_f3/plot_outgroup_f3_heatmap.R <<'RSCRIPT'
#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
})

infile <- "outgroup_f3/OutgroupF3.table"
plotdir <- "plots_outgroup_f3"
dir.create(plotdir, showWarnings = FALSE, recursive = TRUE)

d <- read.table(infile, header = TRUE, stringsAsFactors = FALSE)

d$f3 <- as.numeric(d$f3)
d$stderr <- as.numeric(d$stderr)
d$Z <- as.numeric(d$Z)

pops <- sort(unique(c(d$X, d$Y)))

# Build symmetric matrix
mat <- matrix(NA, nrow = length(pops), ncol = length(pops),
              dimnames = list(pops, pops))

for (i in seq_len(nrow(d))) {
  x <- d$X[i]
  y <- d$Y[i]
  v <- d$f3[i]
  mat[x, y] <- v
  mat[y, x] <- v
}

diag(mat) <- NA

# Cluster order: high f3 = close, so distance = max - f3
dist_mat <- max(mat, na.rm = TRUE) - mat
diag(dist_mat) <- 0
hc <- hclust(as.dist(dist_mat), method = "average")
ord <- hc$labels[hc$order]

# Long format for triangular heatmap
df <- expand.grid(X = pops, Y = pops, stringsAsFactors = FALSE)
df$f3 <- mapply(function(x, y) mat[x, y], df$X, df$Y)

df$X <- factor(df$X, levels = ord)
df$Y <- factor(df$Y, levels = rev(ord))

df$ix <- as.numeric(df$X)
df$iy <- as.numeric(df$Y)

# keep lower triangle in displayed clustered order
n <- length(ord)
df_tri <- df[df$ix <= (n - df$iy + 1), ]

p <- ggplot(df_tri, aes(x = X, y = Y, fill = f3)) +
  geom_tile(color = "white", linewidth = 0.25) +
  scale_fill_viridis_c(option = "inferno", direction = -1, na.value = "grey90") +
  theme_bw(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1),
    axis.title = element_blank(),
    legend.position = "right"
  ) +
  labs(
    title = "Outgroup-f3 shared drift",
    subtitle = paste0("f3(X, Y; ", unique(d$outgroup), ")"),
    fill = "f3"
  )

ggsave(file.path(plotdir, "OutgroupF3_triangle_heatmap.pdf"), p, width = 7.5, height = 6.5)
ggsave(file.path(plotdir, "OutgroupF3_triangle_heatmap.png"), p, width = 7.5, height = 6.5, dpi = 300)

# Save ordered matrix
write.table(
  mat[ord, ord],
  file.path(plotdir, "OutgroupF3_matrix_clustered.tsv"),
  sep = "\t",
  quote = FALSE,
  col.names = NA
)

# Save pairwise table sorted by f3
d_sorted <- d[order(-d$f3), ]
write.table(
  d_sorted,
  file.path(plotdir, "OutgroupF3_pairs_sorted.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

print("[OK] Heatmap saved.")
RSCRIPT

Rscript outgroup_f3/plot_outgroup_f3_heatmap.R

echo
echo "========== Done =========="
echo "Main outputs:"
echo "${OUTDIR}/outgroup_f3/OutgroupF3.table"
echo "${OUTDIR}/plots_outgroup_f3/OutgroupF3_triangle_heatmap.pdf"
echo "${OUTDIR}/plots_outgroup_f3/OutgroupF3_triangle_heatmap.png"
echo "${OUTDIR}/plots_outgroup_f3/OutgroupF3_matrix_clustered.tsv"
echo "${OUTDIR}/plots_outgroup_f3/OutgroupF3_pairs_sorted.tsv"
