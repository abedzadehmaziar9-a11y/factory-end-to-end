# Differential Gene Expression Analysis of Glioblastoma Microstructures

**Insights from the Ivy Glioblastoma Atlas Project (Ivy GAP)**

Computational transcriptomics · bioinformatics · R

---

## Research question

How does gene expression shift as tumor cells migrate from the dense core into
surrounding healthy tissue?

Glioblastoma (GBM) behaves as a "vascular parasite" — it hijacks the brain's blood
supply to fuel rapid, invasive growth. This analysis compares two anatomical zones
within the same tumors:

- **Cellular Tumor (CT)** — the dense, high-density core
- **Infiltrating Tumor (IT)** — the invasive edge where cells migrate into healthy tissue

## Dataset

| | |
|---|---|
| RNA-seq samples | 122 |
| Unique tumors | 10 |
| Structures compared | 2 (IT vs CT) |
| Genes quantified | ~20,000 |

Source: Ivy GAP — high-resolution RNA-seq from laser-microdissected tumor structures,
allowing gene expression to be read at specific anatomical zones within the same tumor.

## Methods

Two complementary analyses, both in R:

**DESeq2 — individual biomarkers**
Differential expression analysis to find genes significantly up- or down-regulated
between IT and CT.
Metrics: log₂ fold change (magnitude), adjusted p-value (confidence), LFC shrinkage
(stabilizes low-count genes).

**WGCNA — gene co-expression networks**
Groups thousands of genes into color-coded modules based on co-expression, revealing
the biological functions behind the disease.
Metrics: kME (hub score) to identify master-regulator genes, module–trait correlation
to link modules to structures, scale-free topology to keep the network biologically
realistic.

### Quality control

- **PCA** — IT and CT do not form two distinct clusters (PC1 27%, PC2 13% variance).
  Transcriptomic differences are real but subtle globally; the invasive front is a
  specialized variation of the core rather than a separate biological entity. This
  justifies the sensitive statistical approach.
- **MA plot** — the gray cloud sits centered on the zero line, confirming successful
  library-size normalization. Observed differences are biological, not technical.

## Key findings

**1. The Blue module is GBM's angiogenic engine**
560 co-expressed genes coordinating tumor neovascularization. Hub genes (kME > 0.93):
KDR (VEGFR2), FLT4 (VEGFR3), EFNB2, FOXC1. Enriched for angiogenesis and endothelial
cell migration. KDR is already the target of bevacizumab (Avastin) in clinical practice.

**2. The Turquoise module is GBM's infiltration engine**
20 high-connectivity hub genes (kME > 0.92) — CHKA (0.937), HSPA12A (0.936),
SYP (0.928), GABRA2 (0.926). Top KEGG pathways cluster around neurotransmitter
signaling: neuroactive ligand signaling, glutamatergic synapse, GABAergic synapse,
cAMP signaling, nicotine addiction. The tumor mimics and hijacks the brain's own
neural machinery to spread.

**3. DESeq2 pinpoints individual biomarkers at the edge**
Genes upregulated in IT separate the microscopic invasive front from the tumor core —
a molecular map that could help plan resection margins and identify cells beyond the
visible tumor boundary.

### Module–trait correlations (N = 122)

| Module × trait | r | Note |
|---|---|---|
| MEpink × ID1 | 0.47 | strongest positive (p ≈ 4×10⁻⁷) — stemness-associated |
| MEblue × PDPN | 0.24 | podoplanin — vascular invasion, supports angiogenic reading |
| MEturquoise × HIF1A | 0.22 | hypoxia-responsive — low-oxygen microenvironment |

## Conclusion

Glioblastoma's aggressive infiltration is biologically powered by its ability to hijack
the brain's own vascular and neural machinery.

- **Surgical planning** — a molecular map of the invasive edge can help guide resection margins.
- **Therapeutic targets** — KDR, FLT4, CHKA, SYP as actionable leads.
- **Next steps** — validate hub genes in independent cohorts; integrate with spatial transcriptomics.

## Repository contents

