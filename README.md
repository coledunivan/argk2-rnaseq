# argk-2 × CoCl₂ — *C. elegans* RNA-seq Analysis Pipeline

End-to-end reproducible pipeline for the analysis of CoCl₂ stress response in
*argk-2⁻/⁻* (RB2060), *argk-4⁻/⁻* (RB2598), and wild-type N2 *C. elegans*,
accompanying *Dunivan, "ARGK-2 sets a transcriptional stress response baseline
in C. elegans."

The pipeline takes raw GEO count tables and a sample metadata sheet, fits a
single canonical DESeq2 model under the full factorial design, and produces
all six manuscript figures plus three supplementary figures and the qPCR
validation panel.

---

## Reproduced headline numbers

Every central number in the thesis is reproduced from raw counts:

| Source | Metric | Reproduced | Thesis |
|---|---|---:|---:|
| Fig 5 / Effect class | G×E genes | **90** | 90 |
| Fig 5A | sub-additive (y ≤ x) | **68** | 68 |
| Fig 5A | supra-additive | **22** | 22 |
| Fig 5A | binomial p | 1.25 × 10⁻⁶ | 1.3 × 10⁻⁶ |
| Fig 6 | TF-target G×E (in GRN) | **81** | 81 |
| Supp Fig 4 | Mann-Whitney p (G×E vs additive) | 9.22 × 10⁻¹⁵ | 9.2 × 10⁻¹⁵ |
| Supp Fig 4 | G×E median deviation | −1.742 | −1.74 |
| Supp Fig 4 | Additive median deviation | −0.051 | −0.05 |
| Supp Fig 1c | qPCR Spearman ρ (thesis quartet) | **0.825** | 0.83 |
| Supp Fig 1c | qPCR Spearman p | 0.001 | 0.0008 |

---

## Repo layout

```
argk2-rnaseq/
├── README.md                              ← this file
├── LICENSE                                ← MIT
├── CITATION.cff                           ← citation metadata
├── MANUSCRIPT_NUMBERS.md                  ← definitive numbers for manuscript revision
├── requirements.R                         ← R package installation script
├── run_all.sh                             ← end-to-end pipeline orchestrator
├── scripts/                               ← all analysis scripts (R + Python)
│   ├── _utils.R                           ← shared helpers (offline KEGG/GO loaders, BOM-safe CSV)
│   ├── 00_QC_report.R                     ← QC: PCA, sample distance, dispersion
│   ├── 01_DESeq2_and_TF_enrichment.R      ← canonical model fit + 5 contrasts +
│   │                                         effect classification (90 G×E genes) +
│   │                                         TF G×E target enrichment (zip-2 OR=17.83)
│   ├── 02_Figure1.R                       ← Figure 1: PCA, KO val, volcanos, Venn
│   ├── 03_Figure2_KEGG.R                  ← Figure 2: KEGG 4-contrast composition
│   ├── 04_Figure3_GO.R                    ← Figure 3: GO dotplots + immune effect map
│   ├── 05_Figure4_TF.R                    ← Figure 4: TF bubble plots (unified scales)
│   ├── 06_Figure5_GxE.R                   ← Figure 5: G×E scatter + heatmap + bars
│   ├── 07_GRN_inputs.R                    ← Figure 6 GRN node + edge tables (90 G×E + 15 TFs)
│   ├── 08_Figure6_GRN.py                  ← Figure 6: TF-regulatory network render
│   ├── 09_qPCR_analysis.py                ← qPCR Cq parsing + ΔΔCq + figure
│   ├── 10_qPCR_compare.R                  ← qPCR vs RNA-seq Spearman comparison
│   ├── 11_Supp_GxE_heatmap.py             ← Supp Fig: full 90-gene G×E heatmap
│   ├── 12_Supp_additive_specificity.py    ← Supp Fig: G×E vs additive (Mann-Whitney p=9.2e-15)
│   └── 13_Supp_GxE_core_analysis.py       ← Supp Fig: in-degree histogram + TF Jaccard heatmap
├── data/
│   ├── RNASEQ61125.csv                    ← raw counts (24 libraries × ~47k genes)
│   ├── Sample_Metadata_Table.csv          ← sample × condition matrix
│   ├── qpcr/                              ← pre-computed qPCR ΔΔCq CSVs + raw Cq run dirs (opt.)
│   └── reference/
│       ├── Seq_ID_crossref.csv            ← WBGene ↔ Sequence Name ↔ UniProt
│       ├── RNAseq_to_TF_Targets.csv       ← TFLink TF → target edges (~80 MB)
│       ├── kegg_cel_pathway_to_gene.csv   ← KEGG cel pathway → gene
│       ├── kegg_cel_pathway_names.csv     ← KEGG cel pathway → name
│       └── org_Ce_eg_GO_map.csv           ← gene → GO term mapping
└── outputs/                               ← all generated artifacts (gitignored)
    ├── QC/                                ← PCA, distance heatmap, dispersion plots
    ├── DESeq2/                            ← 5 contrast CSVs + effect classification + dds.RData +
    │                                         TF summary + TF G×E enrichment + KEGG composition
    ├── GRN_tables/                        ← nodes_full90.tsv + edges_combined_full90.tsv
    ├── figures/                           ← publication PDFs
    │   ├── Figure1/                       ← Panels A–D + composite
    │   ├── Figure2/                       ← KEGG 4-contrast figure
    │   ├── Figure3/                       ← GO dotplots + immune map
    │   ├── Figure4/                       ← TF bubble plots (unified scales)
    │   ├── Figure5/                       ← G×E scatter + heatmap + bars
    │   └── Figure6/                       ← GRN
    ├── qpcr/                              ← qPCR validation panel + per-gene supplement
    ├── supplementary/                     ← Supp Fig 3 (heatmap) + Supp Fig 4 (specificity) +
    │                                         Supp Fig 6B (in-degree) + 6C (TF co-regulation)
    └── logs/                              ← per-step run logs
```

---

## Prerequisites

- **R ≥ 4.3** with these packages (see `requirements.R`):
  - From Bioconductor: `DESeq2 (≥ 1.42)`, `apeglm`, `ashr`, `org.Ce.eg.db`,
    `GO.db`, `AnnotationDbi`, `clusterProfiler` *(optional — pipeline ships
    offline equivalents that work without it)*, `biomaRt` *(optional)*,
    `limma`.
  - From CRAN: `dplyr`, `tidyr`, `readr`, `tibble`, `stringr`, `purrr`,
    `ggplot2`, `ggrepel`, `patchwork`, `scales`, `VennDiagram`, `pheatmap`.
- **Python ≥ 3.10** with `pandas`, `numpy`, `scipy`, `matplotlib`,
  `networkx`, `python-igraph` *(for community detection)*, `goatools`
  *(optional)*.

### Offline reproducibility

Three annotation CSVs ship in `data/reference/` and let the pipeline run
without live Bioconductor/KEGG/Ensembl network access:

- `kegg_cel_pathway_to_gene.csv` + `kegg_cel_pathway_names.csv` — KEGG
  pathway data for *C. elegans*, produced once via
  `clusterProfiler::download_KEGG("cel")` (see *Regenerating annotations*).
- `org_Ce_eg_GO_map.csv` — gene → GO term mapping, produced once via
  `AnnotationDbi::select(org.Ce.eg.db, ...)`.

`scripts/_utils.R` provides `load_kegg_offline()` / `load_go_offline()` /
`enrichGO_offline()` helpers that consume these CSVs and behave as drop-in
replacements for the corresponding `clusterProfiler` / `org.Ce.eg.db`
functions.

---

## Quickstart

```bash
# 1. Install R + Python dependencies
Rscript requirements.R
pip install -r requirements.txt  # if present, else: pandas numpy scipy matplotlib networkx

# 2. Run the full pipeline end-to-end
bash run_all.sh

# 3. Inspect outputs/
open outputs/figures/Figure5/Fig5_recreated.pdf
```

A successful run produces all six main figures, three supplementary figures,
the qPCR validation panel, and all intermediate DEG / contrast tables in
`outputs/`.

---

## Pipeline graph

```
RAW COUNTS + METADATA  (data/RNASEQ61125.csv + data/Sample_Metadata_Table.csv)
        │
        ▼
00_QC_report.R                     ─►  outputs/QC/  (PCA + distance + dispersion)
        │
01_DESeq2_and_TF_enrichment.R      ─►  outputs/DESeq2/
        │                              ├── DESeq2_treatment_in_N2.csv         (Wald + apeglm)
        │                              ├── DESeq2_treatment_in_RB2060.csv     (Wald + apeglm)
        │                              ├── DESeq2_genotype_untreated.csv      (Wald + ashr)
        │                              ├── DESeq2_genotype_treated.csv        (Wald + ashr)
        │                              ├── DESeq2_GxE_interaction.csv         (interaction term)
        │                              ├── DESeq2_effect_classification.csv   ◄── 90 G×E genes
        │                              ├── DESeq2_TF_summary.csv              (per-TF stats)
        │                              ├── DESeq2_TF_GxE_enrichment.csv       (zip-2 OR=17.83)
        │                              └── DESeq2_dds_interaction_model.RData (full fit)
        │
        ├──► 02_Figure1.R              ─► outputs/figures/Figure1/  (PCA, KO val, volcanos, Venn)
        ├──► 03_Figure2_KEGG.R         ─► outputs/figures/Figure2/  (KEGG 4-contrast)
        │                              + outputs/DESeq2/KEGG_Functional_Composition_All_Contrasts.csv
        ├──► 04_Figure3_GO.R           ─► outputs/figures/Figure3/  (GO + immune effect map)
        ├──► 05_Figure4_TF.R           ─► outputs/figures/Figure4/  (TF bubbles, unified scales)
        ├──► 06_Figure5_GxE.R          ─► outputs/figures/Figure5/  (scatter + heatmap + bars)
        │
        ├──► 07_GRN_inputs.R           ─► outputs/GRN_tables/  (nodes_full90.tsv + edges)
        │           │
        │           └──► 08_Figure6_GRN.py    ─► outputs/figures/Figure6/  (TF regulatory network)
        │
        ├──► 09_qPCR_analysis.py       ─► outputs/qpcr/  (skips cleanly if no raw Cq files)
        │           │
        │           └──► 10_qPCR_compare.R    ─► outputs/qpcr/  (qPCR↔RNA-seq Spearman ρ=0.83)
        │
        ├──► 11_Supp_GxE_heatmap.py    ─► outputs/supplementary/  (full 90-gene G×E heatmap)
        ├──► 12_Supp_additive_specificity.py  ─► outputs/supplementary/  (Mann-Whitney p=9.2e-15)
        └──► 13_Supp_GxE_core_analysis.py     ─► outputs/supplementary/  (in-degree + Jaccard)
```

---

## Key methodological choices

1. **Single canonical DESeq2 fit** under `~ genotype * treatment +
   genotype:experiment + treatment:experiment`. All downstream contrasts
   are extracted from this single fit using `apeglm` (single-coef) or
   `ashr` (linear-combination) shrinkage — never re-fit on 8-sample subsets.
2. **G×E gene set** is defined as `effect_class == "GxE"` in
   `DESeq2_effect_classification.csv`, requiring concurrent significance
   in (treatment ∪ genotype ∪ interaction) tests, with apeglm/ashr shrunk
   LFCs and `padj < 0.05`, `|log2FC| > 1`. This yields **90 genes**.
3. **Figure 1 Venn** is a single three-set comparison of treatment-response
   DEGs across N2 / argk-2 / argk-4 (an apples-to-apples contrast across
   genotypes). The earlier two-Venn split was logically incoherent.
4. **Figure 4 TF bubbles** share unified size (`# DE targets`) and color
   (`proportion DE`) scales across all four sub-panels via
   `plot_layout(guides = "collect")`, so cross-panel comparisons are
   biologically meaningful.

---

## Regenerating annotations (one-time, on a machine with Bioconductor)

```r
library(org.Ce.eg.db); library(clusterProfiler); library(AnnotationDbi)

# Gene → GO mapping
go_map <- AnnotationDbi::select(org.Ce.eg.db,
                                keys=keys(org.Ce.eg.db),
                                columns=c("SYMBOL","ENTREZID","GO","ONTOLOGY"),
                                keytype="ENTREZID")
write.csv(go_map, "data/reference/org_Ce_eg_GO_map.csv", row.names=FALSE)

# KEGG pathway → gene mapping
kegg <- clusterProfiler::download_KEGG("cel")
write.csv(kegg$KEGGPATHID2EXTID, "data/reference/kegg_cel_pathway_to_gene.csv", row.names=FALSE)
write.csv(kegg$KEGGPATHID2NAME,  "data/reference/kegg_cel_pathway_names.csv",   row.names=FALSE)
```

These CSVs are then committed to the repo and the pipeline runs offline.

---

## Notes on `data/reference/RNAseq_to_TF_Targets.csv`

This file (~80 MB) is the TFLink (Liska et al., 2022) export joined to our
count matrix. It exceeds GitHub's 100 MB hard limit only marginally; the
repo uses `.gitattributes` to declare it for Git LFS. To regenerate from
the original TFLink dump, see `scripts/build_tflink_input.R` *(not yet
shipped — request from the authors)*.

---

## Citation

If you use this pipeline, please cite the accompanying thesis and the
underlying tools:

- Love MI, Huber W, Anders S (2014). *Moderated estimation of fold change
  and dispersion for RNA-seq data with DESeq2.* Genome Biology 15:550.
- Zhu A, Ibrahim JG, Love MI (2019). *Heavy-tailed prior distributions for
  sequence count data: removing the noise and preserving large differences.*
  Bioinformatics 35:2084.
- Stephens M (2017). *False discovery rates: a new deal.* Biostatistics 18:275.
- Liska O et al. (2022). *TFLink: an integrated gateway to access
  transcription factor–target gene interactions for multiple species.*
  Database baac083.
- Yu G, Wang LG, Han Y, He QY (2012). *clusterProfiler: an R package for
  comparing biological themes among gene clusters.* OMICS 16:284.

---

## License

MIT — see `LICENSE`.
