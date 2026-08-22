############################################################
## 02_deseq2_IT_vs_CT.R
## Differential expression: Infiltrating Tumor vs Cellular Tumor
##
## Produces: PCA, MA plot, volcano plot, and the DE results table.
############################################################

source("R/00_setup.R")

counts      <- readRDS(file.path(RESULTS_DIR, "counts_clean.rds"))
meta_subset <- readRDS(file.path(RESULTS_DIR, "meta_clean.rds"))
gene_map    <- readRDS(file.path(RESULTS_DIR, "gene_map.rds"))

############################
## STEP 1 — Build the DESeq2 object
############################

meta_subset$structure_name <- as.factor(meta_subset$structure_name)

# Cellular Tumor is the reference level, so a positive log2FoldChange
# means "higher in the infiltrating edge".
meta_subset$structure_name <- relevel(meta_subset$structure_name,
                                      ref = CT_LABEL)

dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData   = meta_subset,
  design    = ~ structure_name
)

# Drop genes with almost no signal — noise only slows the fit down
# and costs multiple-testing power.
keep <- rowSums(counts(dds) >= 10) >= 3
dds  <- dds[keep, ]
cat("Genes after filtering:", nrow(dds), "\n")

dds <- DESeq(dds)

resultsNames(dds)   # confirm the contrast name before using it below

############################
## STEP 2 — PCA (quality control)
############################

vsd <- vst(dds, blind = TRUE)

pca_plot <- plotPCA(vsd, intgroup = "structure_name") +
  ggtitle("Corrected PCA: IT vs CT") +
  theme_minimal() +
  theme(
    legend.position  = "bottom",
    legend.direction = "vertical",
    legend.title     = element_blank(),
    legend.text      = element_text(size = 7),
    plot.title       = element_text(hjust = 0.5, face = "bold")
  ) +
  scale_color_discrete(labels = function(x) str_wrap(x, width = 45))

print(pca_plot)
ggsave(file.path(FIG_DIR, "PCA_IT_vs_CT.png"), pca_plot,
       width = 8, height = 6, dpi = 300)

# Reading: IT and CT do NOT split into two clean clusters. The difference
# is real but subtle globally, and PC1 is partly donor-to-donor variation.
# That is the argument for using a sensitive per-gene test rather than
# relying on global clustering.

############################
## STEP 3 — Pull the IT vs CT contrast
############################

res <- results(
  dds,
  contrast = c("structure_name", IT_LABEL, CT_LABEL)
)

summary(res)

############################
## STEP 4 — MA plot (before and after LFC shrinkage)
############################

png(file.path(FIG_DIR, "MA_plot_raw.png"), width = 1800, height = 1400, res = 200)
plotMA(res, ylim = c(-3, 3),
       main = "Global Expression Changes: IT vs CT (Glioblastoma)")
dev.off()

# Shrinkage pulls in the wild fold changes of low-count genes so that
# the biomarkers we report are not just noise from lightly sequenced genes.
res_shrunk <- lfcShrink(
  dds,
  coef = grep("structure_name", resultsNames(dds), value = TRUE)[1],
  type = "apeglm"
)

png(file.path(FIG_DIR, "MA_plot_shrunk.png"), width = 1800, height = 1400, res = 200)
plotMA(res_shrunk, ylim = c(-3, 3),
       main = "LFC-shrunk: IT vs CT (Glioblastoma)")
dev.off()

# Sanity check: the gray cloud should sit centred on zero. If it does,
# library-size normalisation worked and the differences are biological.

############################
## STEP 5 — Attach gene symbols
############################

res_table <- as.data.frame(res_shrunk)
res_table$gene_id     <- rownames(res_table)
res_table$gene_symbol <- gene_map[as.character(res_table$gene_id)]

res_table <- res_table[order(res_table$padj), ]
head(res_table)

############################
## STEP 6 — Significant genes and top hits
############################

sig_genes <- res_table[which(res_table$padj < 0.1), ]
cat("Significant genes (FDR < 0.1):", nrow(sig_genes), "\n")

# If the FDR cut leaves too little to look at, the deck used nominal p
# for the volcano. Keep the two thresholds clearly labelled — they are
# not interchangeable.
sig_nominal <- res_table[which(res_table$pvalue < 0.05), ]
cat("Genes at nominal p < 0.05:", nrow(sig_nominal), "\n")

cat("\nTOP 10 UP IN INFILTRATING TUMOR:\n")
print(head(sig_nominal[order(-sig_nominal$log2FoldChange),
                       c("gene_symbol", "log2FoldChange", "pvalue", "padj")], 10))

cat("\nTOP 10 UP IN CELLULAR TUMOR:\n")
print(head(sig_nominal[order(sig_nominal$log2FoldChange),
                       c("gene_symbol", "log2FoldChange", "pvalue", "padj")], 10))

############################
## STEP 7 — Volcano plot
############################

res_table$direction <- "Non-significant"
res_table$direction[res_table$pvalue < 0.05 & res_table$log2FoldChange >  0.5] <- "Higher in IT"
res_table$direction[res_table$pvalue < 0.05 & res_table$log2FoldChange < -0.5] <- "Higher in CT"

volcano <- ggplot(res_table,
                  aes(x = log2FoldChange, y = -log10(pvalue), colour = direction)) +
  geom_point(alpha = 0.7, size = 1.4) +
  scale_colour_manual(values = c(
    "Higher in IT"    = "#D81B47",
    "Higher in CT"    = "#17A2B8",
    "Non-significant" = "grey70"
  )) +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dotted") +
  labs(
    title    = "Differential Expression: Infiltrating vs Cellular Tumor",
    subtitle = "Top gene candidates (nominal p < 0.05)",
    x        = "Log2 Fold Change",
    y        = "-log10 (nominal p-value)",
    colour   = NULL
  ) +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold"))

print(volcano)
ggsave(file.path(FIG_DIR, "volcano_IT_vs_CT.png"), volcano,
       width = 9, height = 6, dpi = 300)

############################
## STEP 8 — Export
############################

write.csv(res_table,   file.path(RESULTS_DIR, "DESeq2_IT_vs_CT_full.csv"),  row.names = FALSE)
write.csv(sig_genes,   file.path(RESULTS_DIR, "DESeq2_IT_vs_CT_FDR0.1.csv"), row.names = FALSE)

saveRDS(dds, file.path(RESULTS_DIR, "dds_gbm.rds"))
saveRDS(vsd, file.path(RESULTS_DIR, "vsd_gbm.rds"))
