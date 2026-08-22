############################################################
## 04_enrichment.R
## Functional enrichment of the WGCNA modules
##
## GO biological process + KEGG pathways for the Blue and
## Turquoise modules, producing the dot plot and the pathway
## network figure.
##
## NOTE: this step was not part of the course scripts. It is
## written from scratch to reproduce the enrichment figures in
## the presentation.
############################################################

source("R/00_setup.R")
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)

wgcna    <- readRDS(file.path(RESULTS_DIR, "wgcna_objects.rds"))
geneInfo <- wgcna$geneInfo

res_table <- read.csv(file.path(RESULTS_DIR, "DESeq2_IT_vs_CT_full.csv"))

############################
## STEP 1 — Helper: module genes to Entrez IDs
############################

module_entrez <- function(mod) {
  syms <- geneInfo$gene_symbol[geneInfo$module == mod]
  syms <- unique(syms[!is.na(syms) & syms != ""])

  mapped <- bitr(syms,
                 fromType = "SYMBOL",
                 toType   = "ENTREZID",
                 OrgDb    = org.Hs.eg.db)
  mapped
}

# Background: every gene that entered the network, not the whole genome.
# Using the wrong universe inflates enrichment p-values badly.
universe_map <- bitr(unique(na.omit(geneInfo$gene_symbol)),
                     fromType = "SYMBOL",
                     toType   = "ENTREZID",
                     OrgDb    = org.Hs.eg.db)

############################
## STEP 2 — GO enrichment for the Blue module
############################

blue <- module_entrez("blue")
cat("Blue module genes mapped:", nrow(blue), "\n")

go_blue <- enrichGO(
  gene          = blue$ENTREZID,
  universe      = universe_map$ENTREZID,
  OrgDb         = org.Hs.eg.db,
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE
)

head(as.data.frame(go_blue)[, c("Description", "GeneRatio", "p.adjust", "Count")], 15)

p_blue <- dotplot(go_blue, showCategory = 15) +
  ggtitle("Biological Pathways: Module blue")

ggsave(file.path(FIG_DIR, "GO_blue_module.png"), p_blue,
       width = 10, height = 8, dpi = 300)

write.csv(as.data.frame(go_blue),
          file.path(RESULTS_DIR, "GO_BP_blue_module.csv"), row.names = FALSE)

############################
## STEP 3 — GO enrichment for the Turquoise module
############################

turq <- module_entrez("turquoise")
cat("Turquoise module genes mapped:", nrow(turq), "\n")

go_turq <- enrichGO(
  gene          = turq$ENTREZID,
  universe      = universe_map$ENTREZID,
  OrgDb         = org.Hs.eg.db,
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE
)

p_turq <- dotplot(go_turq, showCategory = 15) +
  ggtitle("Biological Pathways: Module turquoise")

ggsave(file.path(FIG_DIR, "GO_turquoise_module.png"), p_turq,
       width = 10, height = 8, dpi = 300)

write.csv(as.data.frame(go_turq),
          file.path(RESULTS_DIR, "GO_BP_turquoise_module.csv"), row.names = FALSE)

############################
## STEP 4 — KEGG pathways for the Turquoise module
############################

kegg_turq <- enrichKEGG(
  gene         = turq$ENTREZID,
  universe     = universe_map$ENTREZID,
  organism     = "hsa",
  pvalueCutoff = 0.05
)

kegg_turq <- setReadable(kegg_turq, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")

head(as.data.frame(kegg_turq)[, c("Description", "GeneRatio", "p.adjust", "Count")], 10)

write.csv(as.data.frame(kegg_turq),
          file.path(RESULTS_DIR, "KEGG_turquoise_module.csv"), row.names = FALSE)

############################
## STEP 5 — KEGG pathway network (cnetplot)
############################

# Colour the nodes by fold change from the DESeq2 run, so the network
# shows direction as well as membership.
fc_vector <- res_table$log2FoldChange
names(fc_vector) <- res_table$gene_symbol
fc_vector <- fc_vector[!is.na(names(fc_vector))]
fc_vector <- fc_vector[!duplicated(names(fc_vector))]

p_net <- cnetplot(
  kegg_turq,
  showCategory = 5,
  foldChange   = fc_vector,
  node_label   = "category",
  cex_label_gene = 0.4
) +
  ggtitle("Turquoise - Top 5 KEGG pathways")

ggsave(file.path(FIG_DIR, "KEGG_network_turquoise.png"), p_net,
       width = 13, height = 9, dpi = 300)

############################
## STEP 6 — KEGG for the Blue module (angiogenesis check)
############################

kegg_blue <- enrichKEGG(
  gene         = blue$ENTREZID,
  universe     = universe_map$ENTREZID,
  organism     = "hsa",
  pvalueCutoff = 0.05
)
kegg_blue <- setReadable(kegg_blue, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")

write.csv(as.data.frame(kegg_blue),
          file.path(RESULTS_DIR, "KEGG_blue_module.csv"), row.names = FALSE)

############################
## STEP 7 — Sanity check on the named hub genes
############################

# The write-up calls out KDR, FLT4, EFNB2, FOXC1 (blue) and
# CHKA, HSPA12A, SYP, GABRA2 (turquoise). Confirm they really are
# in those modules with the kME values reported.
called_out <- c("KDR", "FLT4", "EFNB2", "FOXC1",
                "CHKA", "HSPA12A", "SYP", "GABRA2")

check <- geneInfo[geneInfo$gene_symbol %in% called_out,
                  c("gene_symbol", "module", "kMEblue", "kMEturquoise")]
print(check)
