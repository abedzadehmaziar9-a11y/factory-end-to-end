############################################################
## 03_wgcna.R
## Weighted gene co-expression network analysis on the Ivy GAP samples
##
## Produces: soft-threshold diagnostics, module assignments,
## kME table, hub genes per module, and the module-trait heatmap.
############################################################

source("R/00_setup.R")
library(WGCNA)
library(pheatmap)

enableWGCNAThreads()
options(stringsAsFactors = FALSE)

dds         <- readRDS(file.path(RESULTS_DIR, "dds_gbm.rds"))
meta_subset <- readRDS(file.path(RESULTS_DIR, "meta_clean.rds"))
gene_map    <- readRDS(file.path(RESULTS_DIR, "gene_map.rds"))

############################
## STEP 1 — Variance-stabilised expression
############################

# WGCNA wants a variance-stabilised matrix, not raw counts.
vsd_mat <- getVarianceStabilizedData(dds)
dim(vsd_mat)

############################
## STEP 2 — Keep the most variable genes
############################

# Top 25% by variance. Genes that barely move carry no correlation
# structure, and dropping them makes the network far faster to build.
rv  <- apply(vsd_mat, 1, var)
q75 <- quantile(rv, 0.75)
expr_normalized <- vsd_mat[rv > q75, ]

dim(expr_normalized)

############################
## STEP 3 — Transpose: WGCNA expects samples in rows
############################

input_mat <- t(expr_normalized)
dim(input_mat)          # should be samples x genes

# Flag genes with zero/near-zero variance or too many NAs
gsg <- goodSamplesGenes(input_mat, verbose = 3)
if (!gsg$allOK) {
  input_mat <- input_mat[gsg$goodSamples, gsg$goodGenes]
}
dim(input_mat)

############################
## STEP 4 — Build the trait matrix
############################

# Traits here are the anatomical structure labels, one-hot encoded so
# they can be correlated against module eigengenes.
traitData <- meta_subset[rownames(input_mat), "structure_name", drop = FALSE]
traitData$structure_name <- as.factor(traitData$structure_name)

traitData_numeric <- as.data.frame(
  model.matrix(~ structure_name, data = traitData)[, -1, drop = FALSE]
)
rownames(traitData_numeric) <- rownames(traitData)

# Clean up the very long auto-generated column names for plotting
colnames(traitData_numeric) <- gsub("structure_name", "", colnames(traitData_numeric))
colnames(traitData_numeric) <- gsub(" sampled by reference histology", "",
                                    colnames(traitData_numeric))

head(traitData_numeric)

############################
## STEP 5 — Pick the soft-threshold power
############################

powers <- c(1:10, seq(from = 12, to = 20, by = 2))

sft <- pickSoftThreshold(
  input_mat,
  powerVector = powers,
  networkType = "signed",
  verbose     = 5
)

png(file.path(FIG_DIR, "wgcna_soft_threshold.png"),
    width = 2000, height = 1000, res = 180)
par(mfrow = c(1, 2))
cex1 <- 0.9

plot(sft$fitIndices[, 1],
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft Threshold (power)",
     ylab = "Scale Free Topology Model Fit, signed R^2",
     type = "n", main = "Scale independence")
text(sft$fitIndices[, 1],
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = powers, cex = cex1, col = "red")
abline(h = 0.80, col = "red")

plot(sft$fitIndices[, 1], sft$fitIndices[, 5],
     xlab = "Soft Threshold (power)", ylab = "Mean Connectivity",
     type = "n", main = "Mean connectivity")
text(sft$fitIndices[, 1], sft$fitIndices[, 5],
     labels = powers, cex = cex1, col = "red")
dev.off()

# Choose the lowest power where signed R^2 crosses ~0.8 without
# flattening connectivity. LOOK AT THE PLOT before setting this.
picked_power <- 10
cat("Suggested power from pickSoftThreshold:", sft$powerEstimate, "\n")

############################
## STEP 6 — Build the network
############################

# WGCNA::cor must mask stats::cor or blockwiseModules errors out.
temp_cor <- cor
cor <- WGCNA::cor

netwk <- blockwiseModules(
  input_mat,
  power             = picked_power,
  networkType       = "signed",
  TOMType           = "signed",
  deepSplit         = 2,
  pamRespectsDendro = FALSE,
  minModuleSize     = 30,
  reassignThreshold = 0,
  mergeCutHeight    = 0.25,
  saveTOMs          = FALSE,
  numericLabels     = FALSE,
  verbose           = 3
)

cor <- temp_cor

moduleColors <- netwk$colors
MEs          <- netwk$MEs

# How many genes ended up in each module?
table(moduleColors)

############################
## STEP 7 — Dendrogram
############################

png(file.path(FIG_DIR, "wgcna_dendrogram.png"),
    width = 2000, height = 1200, res = 180)
plotDendroAndColors(
  netwk$dendrograms[[1]],
  moduleColors[netwk$blockGenes[[1]]],
  "Module colors",
  dendroLabels = FALSE,
  hang = 0.03, addGuide = TRUE, guideHang = 0.05,
  main = "Gene dendrogram and module colors"
)
dev.off()

############################
## STEP 8 — Module table with gene symbols
############################

gene_id <- colnames(input_mat)

module_df <- data.frame(
  gene_id     = gene_id,
  gene_symbol = gene_map[as.character(gene_id)],
  module      = moduleColors,
  stringsAsFactors = FALSE
)

write.table(module_df,
            file.path(RESULTS_DIR, "gene_modules_full_table.txt"),
            sep = "\t", row.names = FALSE, quote = FALSE)

############################
## STEP 9 — Module membership (kME)
############################

# kME = correlation between a gene and its module eigengene.
# High kME means the gene sits at the centre of the module — a hub.
kME <- signedKME(input_mat, MEs)

geneInfo <- data.frame(
  gene_id     = colnames(input_mat),
  gene_symbol = gene_map[as.character(colnames(input_mat))],
  module      = moduleColors,
  kME,
  stringsAsFactors = FALSE
)

head(geneInfo[, 1:5])

############################
## STEP 10 — Hub genes per module
############################

genes_by_module <- split(geneInfo, geneInfo$module)

hub_genes_by_module <- lapply(names(genes_by_module), function(mod) {
  df      <- genes_by_module[[mod]]
  kME_col <- paste0("kME", mod)
  if (!kME_col %in% colnames(df)) return(NULL)
  df <- df[order(-abs(df[[kME_col]])), ]
  rownames(df) <- NULL
  df
})
names(hub_genes_by_module) <- names(genes_by_module)
hub_genes_by_module <- hub_genes_by_module[!sapply(hub_genes_by_module, is.null)]

# Top 20 hubs per module, stacked into one table
top_hubs_table <- do.call(rbind, lapply(names(hub_genes_by_module), function(mod) {
  df      <- hub_genes_by_module[[mod]]
  kME_col <- paste0("kME", mod)
  keep_cols <- c("module", "gene_id", "gene_symbol", kME_col)
  keep_cols <- keep_cols[keep_cols %in% colnames(df)]
  df_out <- head(df[, keep_cols, drop = FALSE], 20)
  colnames(df_out)[colnames(df_out) == kME_col] <- "kME"
  df_out$rank_in_module <- seq_len(nrow(df_out))
  df_out
}))

write.csv(top_hubs_table,
          file.path(RESULTS_DIR, "top20_hub_genes_all_modules.csv"),
          row.names = FALSE)

# The two modules the write-up follows:
subset(top_hubs_table, module == "blue")
subset(top_hubs_table, module == "turquoise")

############################
## STEP 11 — Module-trait correlations
############################

moduleTraitCor    <- cor(MEs, traitData_numeric, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nrow(input_mat))

textMatrix <- paste0(signif(moduleTraitCor, 2), "\n(",
                     signif(moduleTraitPvalue, 2), ")")
dim(textMatrix) <- dim(moduleTraitCor)

png(file.path(FIG_DIR, "wgcna_module_trait_heatmap.png"),
    width = 2200, height = 1600, res = 180)
par(mar = c(10, 10, 3, 3))
labeledHeatmap(
  Matrix        = moduleTraitCor,
  xLabels       = colnames(traitData_numeric),
  yLabels       = colnames(MEs),
  ySymbols      = colnames(MEs),
  colorLabels   = FALSE,
  colors        = blueWhiteRed(50),
  textMatrix    = textMatrix,
  setStdMargins = FALSE,
  cex.text      = 0.6,
  zlim          = c(-1, 1),
  main          = paste0("Module-Trait Correlations (N = ", nrow(input_mat), " samples)")
)
dev.off()

write.csv(moduleTraitCor,
          file.path(RESULTS_DIR, "module_trait_correlations.csv"))
write.csv(moduleTraitPvalue,
          file.path(RESULTS_DIR, "module_trait_pvalues.csv"))

############################
## STEP 12 — Heatmap of the top hub genes in one module
############################

plot_module_hubs <- function(mod, n = 20) {
  hubs <- head(hub_genes_by_module[[mod]]$gene_id, n)
  mat  <- expr_normalized[hubs, , drop = FALSE]
  mat  <- t(scale(t(mat)))   # z-score each gene across samples
  rownames(mat) <- gene_map[as.character(hubs)]

  pheatmap(
    mat,
    show_colnames = FALSE,
    main = paste0(str_to_title(mod), " Module: Top ", n,
                  " Hub Genes Transcriptomic Profile"),
    filename = file.path(FIG_DIR, paste0("hubgenes_heatmap_", mod, ".png")),
    width = 12, height = 7
  )
}

plot_module_hubs("turquoise", 20)
plot_module_hubs("blue", 20)

############################
## STEP 13 — Save
############################

saveRDS(list(netwk = netwk, MEs = MEs, moduleColors = moduleColors,
             geneInfo = geneInfo, input_mat = input_mat,
             traitData_numeric = traitData_numeric),
        file.path(RESULTS_DIR, "wgcna_objects.rds"))
