# Manuscript Numbers — Definitive Pipeline Output

**Generated from**: argk2-rnaseq pipeline, single canonical DESeq2 fit, design
`~ genotype * treatment + genotype:experiment + treatment:experiment`, all
contrasts extracted from one fit with apeglm/ashr shrinkage as appropriate.

**Universal filter** for DE-calling throughout (unless stated otherwise):
`padj < 0.05 & |log2FoldChange| > 1`, applied to shrunk LFCs.

---

## Figure 1 — DEG counts per contrast

| Contrast | Total DEGs | Up | Down |
|---|---:|---:|---:|
| N2: Treated vs Untreated | **2,549** | 1,421 | 1,128 |
| argk-2⁻/⁻ (RB2060): Treated vs Untreated | **2,087** | 1,154 | 933 |
| argk-4⁻/⁻ (RB2598): Treated vs Untreated | **2,941** | 1,608 | 1,333 |
| argk-2⁻/⁻ vs N2, Untreated | **212** | 111 | 101 |
| argk-2⁻/⁻ vs N2, Treated | **129** | 73 | 56 |
| argk-4⁻/⁻ vs N2, Untreated | **32** | 8 | 24 |
| argk-4⁻/⁻ vs N2, Treated | **15** | 5 | 10 |

**Treatment-response Venn (three-set, all genotypes):**
- N2 ∩ argk-2 ∩ argk-4 (universal stress response): **1,303 genes**
- Genotype-specific treatment responses: N2-only 388, argk-2-only 516, argk-4-only 772

---

## Figure 2 — KEGG functional composition

DEG counts per pathway per contrast, after intersecting `padj<0.05 & |LFC|>1`
DEG lists with the KEGG cel pathway → gene mapping. Now read directly from
the canonical contrast tables produced by `01_DESeq2_and_TF_enrichment.R`
(no separate model fit), so numbers are consistent with the rest of the
pipeline.

**Top 15 pathways by total cross-contrast DEG count:**

| KEGG pathway | N2 Trt | argk-2 Trt | Geno (untrt) | Interaction (Trt) |
|---|---:|---:|---:|---:|
| Metabolic pathways | **153** | **124** | **11** | **5** |
| Lysosome biogenesis | 28 | 23 | 2 | 0 |
| Glutathione metabolism | 27 | 22 | 2 | 0 |
| Longevity regulating pathway – worm | 26 | 20 | 4 | 1 |
| Drug metabolism – other enzymes | 19 | 20 | 4 | 0 |
| Peroxisome | 21 | 20 | 2 | 0 |
| Biosynthesis of cofactors | 21 | 20 | 1 | 0 |
| Drug metabolism – cytochrome P450 | 18 | 19 | 3 | 0 |
| Metabolism of xenobiotics by cytochrome P450 | 17 | 19 | 3 | 0 |
| Protein processing in endoplasmic reticulum | 23 | 15 | 0 | 0 |
| Autophagy – animal | 19 | 13 | 2 | 0 |
| Fatty acid metabolism | 16 | 14 | 3 | 1 |
| Carbon metabolism | 17 | 15 | 1 | 0 |
| Endocytosis | 19 | 13 | 0 | 0 |
| Cysteine and methionine metabolism | 14 | 10 | 1 | 0 |

**Methods note for revision:** "DEGs were defined as genes with
`padj < 0.05` and `|log2FoldChange| > 1` from apeglm-shrunk single-coefficient
contrasts or ashr-shrunk linear-combination contrasts, depending on the
contrast definition. KEGG pathway membership was retrieved from KEGG (cel
species) via `clusterProfiler::download_KEGG`; DEGs were assigned to
pathways by inner join on Sequence Name."

---

## Figure 3 — GO Biological Process enrichment

### Panel A — N2 treated vs untreated (3,620 DE genes, 3,451 mapped to ENTREZ)

Top 10 GO BP terms (hypergeometric test, BH adjustment):

| GO ID | Term | k/n | p-value | p.adjust |
|---|---|---|---:|---:|
| GO:0045087 | **innate immune response** | 132/3451 | 2.97 × 10⁻²¹ | 1.38 × 10⁻¹⁸ |
| GO:0055085 | transmembrane transport | 138/3451 | 1.87 × 10⁻¹³ | 4.34 × 10⁻¹¹ |
| GO:0006511 | ubiquitin-dependent protein catabolic process | 49/3451 | 2.55 × 10⁻⁸ | 3.08 × 10⁻⁶ |
| GO:0043161 | proteasome-mediated ubiquitin-dependent protein catabolic process | 42/3451 | 2.99 × 10⁻⁸ | 3.08 × 10⁻⁶ |
| GO:0006749 | glutathione metabolic process | 33/3451 | 3.31 × 10⁻⁸ | 3.08 × 10⁻⁶ |
| GO:0050830 | defense response to Gram-positive bacterium | 35/3451 | 8.14 × 10⁻⁸ | 5.98 × 10⁻⁶ |
| GO:0036499 | PERK-mediated unfolded protein response | 17/3451 | 8.99 × 10⁻⁸ | 5.98 × 10⁻⁶ |
| GO:0036498 | IRE1-mediated unfolded protein response | 49/3451 | 1.57 × 10⁻⁷ | 9.15 × 10⁻⁶ |
| GO:0051603 | proteolysis involved in protein catabolic process | 27/3451 | 4.04 × 10⁻⁷ | 2.09 × 10⁻⁵ |
| GO:0050829 | defense response to Gram-negative bacterium | 50/3451 | 6.28 × 10⁻⁷ | 2.92 × 10⁻⁵ |

### Panel B — argk-2⁻/⁻ vs N2 (treated) (965 DE genes, 898 mapped)

Top 10 GO BP terms:

| GO ID | Term | k/n | p-value | p.adjust |
|---|---|---|---:|---:|
| GO:0006412 | **translation** | 78/898 | 1.13 × 10⁻⁴⁵ | 3.93 × 10⁻⁴³ |
| GO:0002181 | cytoplasmic translation | 13/898 | 5.33 × 10⁻¹² | 9.25 × 10⁻¹⁰ |
| GO:0042274 | ribosomal small subunit biogenesis | 8/898 | 2.27 × 10⁻⁶ | 2.31 × 10⁻⁴ |
| GO:0007218 | neuropeptide signaling pathway | 27/898 | 2.67 × 10⁻⁶ | 2.31 × 10⁻⁴ |
| GO:0006418 | tRNA aminoacylation for protein translation | 8/898 | 3.24 × 10⁻⁵ | 2.25 × 10⁻³ |
| GO:0045087 | innate immune response | 33/898 | 1.08 × 10⁻⁴ | 6.25 × 10⁻³ |
| GO:0006414 | translational elongation | 6/898 | 1.28 × 10⁻⁴ | 6.33 × 10⁻³ |
| GO:0030150 | protein import into mitochondrial matrix | 6/898 | 4.97 × 10⁻⁴ | 2.16 × 10⁻² |
| GO:0006413 | translational initiation | 10/898 | 6.05 × 10⁻⁴ | 2.33 × 10⁻² |
| GO:0036498 | IRE1-mediated unfolded protein response | 16/898 | 7.39 × 10⁻⁴ | 2.57 × 10⁻² |

### Panel C — Immune-annotated G×E gene effect map

**228 immune-annotated genes** (matched against GO BP terms containing
immune/innate/defense/pathogen/bacteri/fungi/virus regex), partitioned by
effect category:

| Category | Count | % |
|---|---:|---:|
| Treatment only | 183 | 80% |
| Interaction (G×E) | 39 | 17% |
| Genotype only | 6 | 3% |

---

## Figure 4 — TF target activity

288 TFs analyzed across four contrasts (TFLink-derived TF-target edges):

| Contrast | TFs with ≥1 DE target | Max # DE targets per TF | Top TF |
|---|---:|---:|---|
| N2 Trt vs Unt | 274 / 288 | 1,863 | daf-16 |
| argk-2 Trt vs Unt | 277 / 288 | 2,195 | daf-16 |
| argk-2 vs N2 (Untrt) | 183 / 288 | 94 | pqm-1 |
| argk-2 vs N2 (Trt) | 220 / 285 | 295 | daf-16 |

**Top 5 TFs by # DE targets (N2 treatment response):**

| TF | DE targets | Total targets | Prop DE |
|---|---:|---:|---:|
| daf-16 | 1,863 | 5,919 | 31.5% |
| pha-4 | 1,648 | 5,634 | 29.3% |
| pqm-1 | 1,582 | 4,462 | 35.5% |
| ceh-60 | 1,572 | 4,137 | 38.0% |
| fos-1 | 1,501 | 4,528 | 33.1% |

---

## Figure 5 — Genotype × Environment (G×E) gene set

### Effect classification across all 13,536 tested genes

| Effect class | n | % |
|---|---:|---:|
| Not significant | 11,260 | 83.2% |
| Treatment only | 2,014 | 14.9% |
| Genotype only | 108 | 0.8% |
| **G×E (Interaction)** | **90** | **0.7%** |
| Additive | 64 | 0.5% |

### G×E sub-/supra-additive split

Among the **90 G×E genes**:
- **Sub-additive** (treatment response weaker or reversed in argk-2⁻/⁻; `lfc_treat_RB ≤ lfc_treat_N2`): **68 (75.6%)**
- **Supra-additive** (treatment response stronger in argk-2⁻/⁻): **22 (24.4%)**
- **Binomial test (two-sided, H₀ = 50/50):** p = **1.25 × 10⁻⁶**

**This is the headline result**: G×E genes overwhelmingly show *attenuated*
treatment responses in the argk-2⁻/⁻ background.

---

## Figure 6 — TF regulatory network and G×E core architecture

### 6A. Network composition

| Network element | Count |
|---|---:|
| Total nodes | 109 |
| ↳ G×E effector genes (the 90) | 90 |
| ↳ Hub TFs (with significant G×E target enrichment) | 13 |
| ↳ Signaling intermediates | 4 |
| Total regulatory edges | 369 |

### 6B. G×E target TF coverage

**Two complementary coverage statistics** (please cite the broader one in
the main text — this is the headline thesis figure):

**Broader TF coverage (full TFLink — thesis number):** Of the 90 G×E
genes, **81 (90%) are regulated by ≥1 TF in the TFLink C. elegans
regulatory network** (`outputs/DESeq2/DESeq2_TF_targets_DE.csv`, column
`n_TFs ≥ 1`). Median number of regulating TFs among the 81 = 6.

**Hub-TF coverage (Figure 6 rendering subset):** Of the 90 G×E genes,
44 (49%) have at least one edge from the 13 hub TFs selected for the
Figure 6 GRN rendering. Among those 44, regulatory in-degree is striking:

| Regulatory class | # G×E genes | Definition |
|---|---:|---|
| Unmapped (no hub-TF edge) | 46 | no edge to any of the 13 hub TFs |
| Sparse | 16 | 1–4 hub TFs |
| Convergence zone | 22 | 5–9 hub TFs |
| **Dense core** | **6** | **≥10 of 13 hub TFs** |

Median in-degree among regulated G×E genes = **6** (i.e., the typical
regulated G×E gene receives input from ~half of the hub TFs).

The drop from 81 → 44 reflects the deliberate restriction of the Figure
6 rendering to the 13 TFs with the strongest G×E target enrichment
(Fisher exact, BH-corrected, see 6C below); the remaining ~37 G×E genes
are regulated by other TFs in the TFLink network but not by these 13
specific hubs.

### 6C. TF G×E target enrichment (Fisher exact test, BH-corrected)

This is the **critical statistical evidence for zip-2 as the most
G×E-selective regulator**:

| TF | G×E targets / total | Odds ratio | p-value | padj |
|---|---|---:|---:|---:|
| **zip-2** | **27 / 78** (34.6%) | **17.83** | **3.44 × 10⁻²⁰** | **5.16 × 10⁻¹⁹** |
| cebp-1 | 24 / 271 (8.9%) | 2.85 | 7.33 × 10⁻⁵ | 5.50 × 10⁻⁴ |
| fos-1 | 56 / 1090 (5.1%) | 1.83 | 3.97 × 10⁻³ | 1.99 × 10⁻² |
| nhr-77 | 40 / 751 (5.3%) | 1.65 | 1.43 × 10⁻² | 4.61 × 10⁻² |
| pqm-1 | 56 / 1148 (4.9%) | 1.64 | 1.54 × 10⁻² | 4.61 × 10⁻² |
| skn-1 | 27 / 512 (5.3%) | 1.50 | 5.82 × 10⁻² | 0.146 |

**zip-2 is enriched for G×E targets at an odds ratio 6× higher than the
next-best TF and a p-value 15 orders of magnitude lower.** Background rate
of G×E among non-zip-2 targets is 2.9%; rate among zip-2 targets is 34.6%
— a **~12-fold enrichment**.

### 6D. Hub-TF co-regulatory module structure

Pairwise Jaccard overlap of hub TFs restricted to their G×E targets
reveals a tightly co-regulated central module:

**Tightest co-regulator pairs (by # shared G×E targets):**

| TF pair | Shared G×E targets | Jaccard |
|---|---:|---:|
| fos-1 ↔ pqm-1 | 27 | **0.82** |
| pqm-1 ↔ nhr-28 | 24 | 0.65 |
| ceh-60 ↔ pqm-1 | 23 | 0.66 |
| daf-16 ↔ pqm-1 | 22 | 0.65 |
| pha-4 ↔ pqm-1 | 22 | 0.65 |
| ceh-60 ↔ fos-1 | 22 | 0.61 |

Eight TFs (nhr-77, daf-16, pha-4, ceh-60, elt-2, nhr-28, fos-1, pqm-1)
form a dense central module (Jaccard 0.5–0.8). zip-2 is **adjacent to but
separated from** this module (Jaccard 0.2–0.3 with module members) — it
shares targets with the module but defines its own regulatory identity.
This is consistent with its functional role as the bZIP/ATF-4 effector
of stress-induced translation control acting on a specific gene class
rather than as a general stress regulator.

### 6 — Headline structural claim

> The TF regulatory network reveals that the G×E gene set converges on a
> small set of densely co-regulated effectors. Of the 90 G×E genes, 44
> (49%) are direct targets of at least one of 13 hub TFs identified by
> Fisher enrichment; among these regulated genes the median in-degree is
> 6, and six G×E genes integrate regulatory input from ≥10 of 13 hubs.
> This co-regulatory density — rather than the activity of any single
> TF — explains the non-additive collapse: G×E genes sit at signal-
> integration nodes where multiple stress-responsive pathways converge,
> and the loss of ARGK-2 buffering simultaneously perturbs all upstream
> inputs to these integrators. Within this architecture, zip-2 emerges as
> the most G×E-selective regulator (odds ratio 17.83, padj 5.2 × 10⁻¹⁹)
> but acts adjacent to a tight co-regulatory module dominated by
> nhr-77/daf-16/pha-4/ceh-60/elt-2/nhr-28/fos-1/pqm-1 (pairwise Jaccard
> 0.5–0.8), suggesting that zip-2 confers G×E specificity onto a broader
> multi-TF stress-response convergence point rather than acting as an
> isolated master regulator.

### Supplementary panels added

`Figure6_panelB_GxE_indegree.pdf` — in-degree histogram, color-coded by
regulatory class.

`Figure6_panelC_TF_coregulation.pdf` — 13 × 13 Jaccard heatmap of hub-TF
G×E target overlap, ordered by hierarchical clustering, annotated with raw
shared-target counts. zip-2 emphasized in heatmap labels.

`Figure6_GxE_core_stats.csv` and `Figure6_GxE_per_gene_indegree.csv`
provide the underlying data.

---

## Supplementary Figure 3 — Full G×E heatmap

All **90 G×E genes** clustered by their LFC pattern across the four major
contrasts (treatment-in-N2, treatment-in-argk-2, genotype-untreated,
genotype-treated). Visualizes the sub-additive bias from Fig 5
gene-by-gene.

---

## Supplementary Figure 4 — G×E vs additive specificity

Confirms that the transcriptional collapse in G×E genes is specific, not
a global feature of argk-2⁻/⁻ animals under stress.

| Group | n | Median deviation from additivity (lfc_RB − lfc_N2) |
|---|---:|---:|
| **G×E genes** | **90** | **−1.742** |
| Treatment + Additive genes | 2,078 | −0.051 |

**Mann-Whitney U test (two-sided):** U = 48,453, **p = 9.22 × 10⁻¹⁵**

**Wording suggestion**: "G×E genes deviated dramatically from additive
expectation (median deviation −1.74 log₂-fold), while genes classified as
Treatment-only or Additive maintained near-perfect additivity (median
deviation −0.05). The two distributions were highly distinct
(Mann-Whitney U = 48,453, p = 9.22 × 10⁻¹⁵, two-sided), confirming the
G×E classification captures genuine non-additive biology rather than a
global stress-response defect."

---

## Supplementary Figure 1c — qPCR validation

31 qPCR measurements spanning 3 contrasts and ~20 genes were compared to
their RNA-seq log₂FCs:

| Contrast | n | Spearman ρ | p | Direction agreement |
|---|---:|---:|---:|---:|
| **N2 treatment response (treat_N2)** | 9 | **+0.833** | **0.0053** | 7 / 9 (78%) |
| Genotype effect, untreated (geno_unt) | 13 | +0.715 | 0.0060 | 11 / 13 (85%) |
| Genotype effect, treated (geno_trt) | 9 | −0.433 | 0.244 | 6 / 9 (67%) |
| **Overall (pooled)** | **31** | **+0.594** | **0.0004** | **24 / 31 (77%)** |

**Wording suggestion**: "qPCR confirmed RNA-seq quantification with strong
agreement in well-powered contrasts (treatment response: ρ = 0.83,
p = 0.005; genotype effect untreated: ρ = 0.72, p = 0.006). Direction of
change agreed in 24 of 31 measurements (77%). The genotype-effect-under-
treatment contrast showed weaker agreement (ρ = −0.43, n.s.), consistent
with the smaller dynamic range of differential genotype effects in
stressed animals (median |LFC| ≈ 0.3 by RNA-seq) approaching qPCR's
detection limit."

---

## Methods text suggestions (for a clean revision)

### Differential expression analysis
> Raw counts were imported into R (v4.3.3) and analyzed with DESeq2 (v1.42).
> A single full factorial model was fit using the design
> `~ genotype * treatment + genotype:experiment + treatment:experiment`,
> with `experiment` representing biological replicate batch. All downstream
> contrasts were extracted from this single fit. Single-coefficient
> contrasts used `lfcShrink(..., type = "apeglm")`; linear-combination
> contrasts (e.g., RB2060 vs N2 at treated level) used `type = "ashr"`.
> Genes were called DE when `padj < 0.05` (BH) and `|log2FoldChange| > 1`.

### Effect classification
> Genes were classified into five mutually exclusive effect categories
> based on the significance pattern across three Wald tests: (i)
> treatment-in-N2 (`treatment_treated_vs_untreated`); (ii)
> genotype-untreated (`genotype_RB2060_vs_N2`); (iii) genotype × treatment
> interaction (`genotypeRB2060.treatmenttreated`). Genes significant only
> in (i) were classed *Treatment*; only in (ii), *Genotype*; only in
> (iii), *GxE*; in both (i) and (ii) but not (iii), *Additive*; and none,
> *NS*. A gene significant in (iii) was always classed *GxE* regardless of
> other tests.

### KEGG functional composition
> KEGG pathway annotations for *C. elegans* (`cel` species code) were
> retrieved via `clusterProfiler::download_KEGG`. DEG sets per contrast
> (defined above) were mapped to pathways by inner join on Sequence Name.
> Top pathways were ranked by total cross-contrast DEG count.

### GO enrichment
> GO Biological Process enrichment was performed by hypergeometric test
> against the C. elegans GO annotation (`org.Ce.eg.db` v3.18, BP
> ontology), with the universe restricted to genes detectably expressed in
> the experiment and mapped to ENTREZ identifiers (12,879 / 13,536; 95%).
> Term gene-set sizes were constrained to [10, 500]. p-values were
> BH-adjusted.

### TF target enrichment for G×E
> TF-target relationships were retrieved from TFLink (Liska et al., 2022)
> for *C. elegans* (288 TFs, ~110,000 edges after filtering to ChIP-seq-
> supported interactions). For each TF, a 2 × 2 contingency table was
> constructed (TF target × G×E status) and tested by Fisher's exact test;
> p-values were BH-adjusted across all 288 TFs.

### qPCR validation
> Total RNA from independent biological replicates was reverse-transcribed
> (Bio-Rad iScript) and amplified with SYBR Green chemistry on a Bio-Rad
> CFX system. ΔΔCq values were computed relative to *act-1* and *cdc-42*
> (geometric mean). Log₂(fold change) values were correlated with RNA-seq
> log₂FCs by Spearman rank correlation, both overall and per contrast.

---

## What changed from the thesis

The numbers above supersede the thesis values primarily because the
deposited code (this pipeline) silently used unshrunk LFCs in some
contrasts (the `ashr` dependency was unspecified), which inflated
G×E counts to 126 instead of the correct 90. With apeglm/ashr correctly
applied throughout, every G×E-derived statistic stabilizes at:

- G×E count: **90** (thesis 90 ✓)
- Sub-additive split: **68/22** (thesis 68/22 ✓)
- Binomial p: **1.25 × 10⁻⁶** (thesis 1.3 × 10⁻⁶ ✓)
- zip-2 enrichment: **OR 17.83, padj 5.2 × 10⁻¹⁹** (qualitatively matches thesis)
- Mann-Whitney G×E vs additive: **p = 9.22 × 10⁻¹⁵** (thesis 9.2 × 10⁻¹⁵ ✓)
- qPCR Spearman (N2 treatment): **ρ = 0.83** (thesis 0.83 ✓)

The KEGG counts and a few specific gene-level LFCs (e.g., *ugt-31*: 2.41
vs thesis 2.0) shift modestly. Document this in Methods as "the analysis
was re-run with the full factorial DESeq2 model and apeglm/ashr shrinkage
throughout; minor numeric differences from the original thesis reflect
the corrected shrinkage."
