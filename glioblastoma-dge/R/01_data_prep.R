############################################################
## 01_data_prep.R
## Glioblastoma (Ivy GAP) — load, subset, and align the data
##
## Goal: end up with a counts matrix and a metadata table whose
## columns and rows line up exactly, restricted to the reference
## histology samples.
############################################################

source("R/00_setup.R")

############################
## STEP 1 — Load the three raw files
############################

columns_samples <- read_csv(file.path(DATA_DIR, "columns-samples.csv"))
row_genes       <- read_csv(file.path(DATA_DIR, "rows-genes.csv"))
counts_raw      <- read_csv(file.path(DATA_DIR, "raw_counts_table.csv"))

dim(columns_samples)
dim(counts_raw)

############################
## STEP 2 — Explore the metadata
############################

colnames(columns_samples)

# What structures exist, and how many samples of each?
table(columns_samples$structure_name)

# How many samples per tumor?
table(columns_samples$tumor_id)

# Cross-tab: does one tumor contribute several samples of the same structure?
# (Yes — laser microdissection of different regions. These are biological
#  replicates within a tumor, not technical replicates.)
table(columns_samples$tumor_id, columns_samples$structure_name)

############################
## STEP 3 — Subset to reference histology samples
############################

# Ivy GAP includes both reference-histology samples and ISH-control samples.
# We keep the reference histology set only.
reference_rows <- grep("reference histology",
                       columns_samples$structure_name,
                       ignore.case = TRUE)

meta_subset <- as.data.frame(columns_samples[reference_rows, ])

dim(meta_subset)
table(meta_subset$structure_name)

############################
## STEP 4 — Format the counts matrix
############################

counts <- as.data.frame(counts_raw)
rownames(counts) <- counts[[1]]   # first column holds the gene IDs
counts <- counts[, -1]            # drop the now-redundant ID column

head(colnames(counts))            # these should be rna_well_id values

############################
## STEP 5 — Align samples between counts and metadata
############################

# rna_well_id is the key that links the two tables.
rownames(meta_subset) <- meta_subset$rna_well_id

common_samples <- intersect(colnames(counts), rownames(meta_subset))
length(common_samples)

counts      <- counts[, common_samples]
meta_subset <- meta_subset[common_samples, ]

# This MUST return TRUE before DESeq2 will give meaningful answers.
stopifnot(all(colnames(counts) == rownames(meta_subset)))

############################
## STEP 6 — Map Ensembl IDs to gene symbols
############################

# Ivy GAP ships its own gene annotation, so we do not need biomaRt here.
colnames(row_genes)

gene_annotation <- row_genes[, c("gene_id", "gene_symbol")]

# Keep the mapping as a lookup table rather than overwriting rownames now —
# DESeq2 is happier with stable unique IDs, and duplicate symbols exist.
gene_map <- setNames(gene_annotation$gene_symbol,
                     as.character(gene_annotation$gene_id))

head(gene_map)

############################
## STEP 7 — Save the cleaned objects
############################

saveRDS(counts,          file.path(RESULTS_DIR, "counts_clean.rds"))
saveRDS(meta_subset,     file.path(RESULTS_DIR, "meta_clean.rds"))
saveRDS(gene_map,        file.path(RESULTS_DIR, "gene_map.rds"))

write.csv(counts,      file.path(RESULTS_DIR, "counts_clean.csv"))
write.csv(meta_subset, file.path(RESULTS_DIR, "meta_clean.csv"), row.names = FALSE)

cat("Samples kept:", ncol(counts), "\n")
cat("Genes kept:  ", nrow(counts), "\n")
