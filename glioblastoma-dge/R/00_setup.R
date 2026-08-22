############################################################
## 00_setup.R
## Glioblastoma (Ivy GAP) — IT vs CT
## Package installation and shared project configuration
############################################################

## ---- Install (run once) -----------------------------------------------

# install.packages("BiocManager")
# BiocManager::install(c(
#   "DESeq2", "WGCNA", "biomaRt", "apeglm",
#   "clusterProfiler", "org.Hs.eg.db", "enrichplot"
# ))
# install.packages(c("readr", "dplyr", "tidyr", "ggplot2",
#                    "stringr", "pheatmap"))

## ---- Load --------------------------------------------------------------

library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(DESeq2)

## ---- Project paths -----------------------------------------------------
## Edit DATA_DIR to point at your own copy of the Ivy GAP files.
## On Posit Cloud this was: /cloud/project/Spring_2026/Datasets/glioblastoma_data

DATA_DIR    <- "data/glioblastoma_data"
RESULTS_DIR <- "results"
FIG_DIR     <- "figures"

dir.create(RESULTS_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(FIG_DIR,     showWarnings = FALSE, recursive = TRUE)

## Expected input files inside DATA_DIR:
##   columns-samples.csv    sample metadata (structure_name, tumor_id, rna_well_id)
##   rows-genes.csv         gene annotation (gene_id, gene_symbol)
##   raw_counts_table.csv   gene x sample raw count matrix

## ---- The two structures we compare -------------------------------------
## Ivy GAP writes these out in full; keep the exact strings in one place
## so every downstream script uses the same labels.

IT_LABEL <- "Infiltrating Tumor sampled by reference histology"
CT_LABEL <- "Cellular Tumor sampled by reference histology"

## ---- Reproducibility ---------------------------------------------------

set.seed(42)
options(timeout = max(300, getOption("timeout")))  # biomaRt can be slow

sessionInfo()
