# bioana

Analysis for the Biodes suite. Bioana takes a **Biodes design folder filled
with data** — simulated by Simitall or measured by you — and runs, for every
data set in it, the analysis that data type calls for, with the model the
design calls for. When the data came from Simitall it also scores each result
against the simulated truth.

```
question --Biodes--> design folder --Simitall or the lab--> data + data_manifest.json --Bioana--> results (+ truth scores) --Biopaint--> figures
```

| agent | repo | does |
|---|---|---|
| Biodes | [BioDesignGPT](https://github.com/nirwan1265/BioDesignGPT) (R package `studyforge`) | reads the question, designs the experiment: `study_spec.json`, `field_layout.csv` |
| Simitall | [simitall](https://github.com/nirwan1265/simitall) | simulates the data the design would produce, with known truth |
| **Bioana** | this repo (R package `bioana`) | analyses the data with the design's model |
| Biopaint | (later) | figures |

Bioana never reads the question: Biodes already turned it into a spec. It
never simulates: it treats simulated and measured data the same way.

## Install

```r
# install.packages("remotes")
remotes::install_github("nirwan1265/bioana")
```

Needs `jsonlite`, `Matrix` and `lme4`; `lmerTest`, `MASS`, `ordinal` and
`survival` are used for mixed models and non-Gaussian outcomes, `rrBLUP`,
`BGLR` and `ranger` for GWAS and genomic prediction, `ggplot2` and `patchwork`
for figures, when installed. Tests that simulate data need `simitall`.

## Quick start

Simulated data (Simitall writes `data_manifest.json` for you):

```r
simitall::simitall_from_biodes("~/biodes_runs/my_design")
bioana::bioana_analyze("~/biodes_runs/my_design")
#> bioana: design 'rcbd' (y ~ treatment + block)
#>   [outcome] continuous outcome, lm ; truth: 1/1 coefficients in their CI
#>   [rnaseq] 500 features tested, treatmenttreatment_2 10 significant; truth: treatmenttreatment_2 found 9/50 (FDR 0.1)
#>   [methylation] 4809 features tested, ...; truth: found 80/119 (FDR 0.059)
```

Your own data: put the files in the design folder (matrix columns = the
`sample_id` values of `sample_metadata_template.tsv`), then say what they are:

```r
library(bioana)
bioana_write_manifest("~/biodes_runs/my_design",
  list(id = "outcome", data_type = "outcome", files = list(table = "sample_metadata_template.tsv"),
       column = "primary_outcome", outcome_type = "continuous"),
  list(id = "rnaseq", data_type = "feature_counts", assay = "rnaseq", files = list(counts = "counts.tsv")))
bioana_analyze("~/biodes_runs/my_design")
```

## What it analyses

`bioana_analyses()` prints `inst/rules/analyses.tsv`:

| data_type | files | assays | analysis |
|---|---|---|---|
| `outcome` | table | primary outcome | the design's `analysis.model`: lm / lmer (continuous), glm / glmer binomial (binary), negative binomial (count), cumulative logit (ordinal), Cox (survival, needs an `event` column) |
| `feature_counts` | counts | RNA-seq, 3' tag, small RNA, microbiome, single-cell pseudobulk | log2-CPM, per-feature linear model with the design's fixed terms, empirical-Bayes moderated t (as limma) |
| `peak_counts_long` | counts, libraries | ChIP-seq, DAP-seq, CUT&Tag | IP reads per peak and library, mapped to units; then as above |
| `read_positions` | reads, libraries | ATAC-seq, DNase-seq, FAIRE-seq | reads counted in windows (no peak calling); then as above |
| `bisulfite_counts` | sites, methylated, coverage | WGBS, RRBS, EM-seq | sites pooled in windows per context, M-values, moderated t |
| `methylation_beta` | beta | 450k / EPIC arrays | M-values, moderated t per probe |
| `ip_input_counts` | ip, input | m6A MeRIP-seq, MeDIP, RIP | log2(IP/input) per window, so expression changes cancel |
| `genotypes` | vcf | GWAS, GBS, WGS, QTL mapping | line means (BLUEs) of the primary outcome per line from the design's `analysis.model`, then a kinship mixed-model GWAS (rrBLUP, or native EMMAX without it); scored against `qtn_truth.tsv` by exact marker or LD |
| `single_cell_counts` | counts, cells | scRNA-seq (cell level) | planned in `bioana_analyze()`; `cluster_scrnaseq()` works directly |

Settings per data type: `inst/rules/method_parameters.tsv`, overridden with
`bioana_analyze(dir, params = list(bisulfite_counts = list(window_bp = 500)))`.

## Other analyses (moved from Simitall)

These work on files directly (Simitall's simulator outputs or the same formats
from real data) and are not yet wired into `bioana_analyze()`:

| area | functions |
|---|---|
| GWAS | `analyze_gwas()` (mixed model with kinship and PCs, rrBLUP), `analyze_binary_gwas()`, `benchmark_gwas()` (recovery of causal variants), `plot_gwas_results()` (Manhattan, QQ) |
| Genomic prediction and selection | `calculate_genomic_kinship()`, `fit_genomic_model()` (GBLUP, RR-BLUP, BGLR Bayesian models, random forest), `predict_genomic_values()`, `cross_validate_genomic_prediction()`, `select_parents()`, `design_parent_crosses()`, `run_genomic_selection()` |
| eQTL and RNA-seq | `benchmark_eqtl()`, `benchmark_celltype_eqtls()` (donor pseudobulk per cell type), `summarize_rnaseq_results()`, `plot_rnaseq_eqtl_results()` |
| Single cell | `cluster_scrnaseq()`, `score_scrnaseq_clusters()`, `test_scrnaseq_conditions()`, `score_scrnaseq_conditions()`, `plot_scrnaseq_results()` |
| ChIP / ATAC peaks | `call_peaks_windowed()`, `score_peak_calls()`, `test_differential_peaks()`, `score_differential_peaks()`, `plot_chipseq_results()` |
| Genes and GO | `read_gene_coordinates()` (bundled catalogues: Arabidopsis, barley, rice, wheat, maize, fly), `read_go_gaf()`, `join_gene_go()`, `summarize_go_terms()` |
| Genome assembly | `run_unicycler_grid()`, `run_quast_grid()`, `summarize_quast()` |
| Figures | `export_bioana_plot_panels()` saves each panel of a results figure |

```r
g <- bioana::analyze_gwas("cohort.vcf", phenotypes, "gwas_out", trait = "yield", n_pcs = 3)
bioana::benchmark_gwas(g$paths$results, "cohort.causal.tsv")
fit <- bioana::fit_genomic_model("cohort.geno.tsv", phenotypes, trait = "yield", model = "gblup")
```

## Output

`<folder>/bioana/`:

- `<id>/results.tsv` — outcome: one row per coefficient (`term, estimate, se, df, statistic, p, role, padj`);
  features: one row per feature and treatment contrast (`feature_id, contrast, estimate, se, t, p, padj`).
- `<id>/recovery.tsv` — outcome only, simulated data: true vs estimated coefficients.
- genotypes: `<id>/results.tsv` one row per tested marker (`marker_id, seqname, pos, estimate, se, p, padj, maf, significant`;
  `estimate` = effect per alternative allele, on the outcome scale); `<id>/line_means.tsv` (`line, mean, se, n_units`);
  `<id>/qtn_recovery.tsv`, simulated data: one row per QTN (`marker_id, seqname, pos, type, effect, maf, tested, recovered,
  estimate, se, p, padj`); `<id>/gwas.gwas.tsv` is `analyze_gwas()`'s own table, for `plot_gwas_results()`.
- `bioana_summary.json` — per data set: method, model fitted, counts tested and significant,
  warnings, and the truth scores (`tpr`, `fdr`, `sign_agreement`, `estimate_cor`, CI coverage).

## The link to Simitall: `data_manifest.json`

```json
{
  "manifest_version": "0.1",
  "written_by": "simitall 0.2.0",
  "data": [
    {"id": "outcome", "data_type": "outcome", "files": {"table": "sample_metadata_template.tsv"},
     "column": "primary_outcome", "outcome_type": "continuous",
     "truth": {"kind": "outcome", "file": "sim_truth.json"}},
    {"id": "rnaseq", "data_type": "feature_counts", "assay": "rnaseq", "files": {"counts": "counts.tsv"},
     "truth": {"kind": "features", "file": "truth/rnaseq.tsv", "term": "treatment", "large": 1}},
    {"id": "methylation", "data_type": "bisulfite_counts", "assay": "methylation",
     "files": {"sites": "methylation/methylation.sites.tsv.gz", "methylated": "...", "coverage": "..."},
     "truth": {"kind": "regions", "file": "truth/methylation.tsv", "term": "treatment"}}
  ]
}
```

Each simulated data type in Simitall maps to one `data_type` here; a new
simulator in Simitall adds a manifest entry, and a new analysis here adds a
row to `analyses.tsv`. Truth files: `features` = `feature_id, condition,
effect`; `regions` = `region_id, chr, start, end, condition, effect`; the
tested contrast is `term` + `condition` (e.g. `treatmenttreatment_2`).

## Future analyses to check

Each needs an analysis here, a row in `inst/rules/analyses.tsv`, a test, and
(for simulated data) the matching Simitall manifest entry. "Questions" = how
many of Biodes' 4,695 paper questions need it.

**Wire existing functions into `bioana_analyze()`**

- [x] `genotypes`: GWAS / QTL mapping with `analyze_gwas()` on the design's
      line means, scored with `benchmark_gwas()` against `genetic_truth/qtn_truth.tsv`.
      Next for it: GWAS per environment under GxE, ordinal and survival outcomes,
      QTL mapping with cross-aware genotype probabilities.
- [ ] eQTL designs: `benchmark_eqtl()` on `counts.tsv` + `genotypes.vcf`.
- [ ] Single cell at cell level: `cluster_scrnaseq()`, cell-type proportions,
      `test_scrnaseq_conditions()`; cell-type eQTL.
- [ ] Peaks: `call_peaks_windowed()` + `test_differential_peaks()` as the
      peak-calling alternative to fixed windows for ChIP / ATAC.
- [ ] Genomic prediction for breeding designs (`fit_genomic_model()`,
      `cross_validate_genomic_prediction()`), scored against the true genetic values.

**Analyses for data types Simitall will simulate**

- [ ] Metabolomics / lipidomics (54 questions): log intensities, missing
      values, moderated t; pathway / class summaries.
- [ ] Genetics with no lines (36): BSA / QTL-seq (allele-frequency difference
      between bulks, G' statistic), pool-seq (allele frequencies over
      generations), population scans (Fst, selection).
- [ ] Ribo-seq (42): translational efficiency, the RNA x footprint interaction.
- [ ] scATAC-seq (41): peaks x cells, cluster-level differential accessibility.
- [ ] Hi-C (38): differential contacts, loops and TAD boundaries.
- [ ] MeDIP (31) and RIP/CLIP (14): already `ip_input_counts`; check on their
      simulated data.
- [ ] Degradome (22): cleavage sites and miRNA targets.
- [ ] Spatial (19): spatially variable genes, domains.
- [ ] MNase-seq (16): nucleosome positions and occupancy changes.
- [ ] CRISPR screens: guide-level counts to gene scores (MAGeCK-like).
- [ ] Amplicon genotyping: editing efficiency, allele frequencies per sample.
- [ ] TCR/BCR repertoires: clonotype abundance and diversity.
- [ ] Proteomics: log intensities with missing values (limma-like).
- [ ] STARR-seq / MPRA: RNA over DNA barcode counts per element.
- [ ] Compositional outcomes; survival with random effects (frailty / coxme).

**Methods and checks**

- [ ] Package adapters when installed, with the native methods as fallback:
      DESeq2 / edgeR / limma-voom (counts), DSS / methylKit (methylation),
      ANCOM-BC / ALDEx2 (microbiome, compositional), exomePeak2 (m6A),
      csaw / DiffBind (peaks).
- [ ] Random terms in feature tests (block or subject correlation, as limma
      duplicateCorrelation / dream); the feature engine now uses fixed terms only.
- [ ] Contrasts beyond "each level vs the reference": pairwise and interaction
      contrasts for multi-level factors.
- [ ] Power for every data type: Simitall's power sweep uses Bioana's scores
      for the outcome and counts; add methylation (DMRs), m6A, peaks, GWAS.
- [ ] Real-data readers: featureCounts / Salmon tables, Bismark coverage,
      MACS2 peaks, VCF; QC and normalisation reports.
- [ ] Check low power where the simulation says so: small RNA and microbiome
      templates vary a lot between real samples (0 found with 5 per group);
      m6A at 5 per group (2/40).

## Status

Version 0.3: `bioana_analyze()` for the data types in the first table, and
all of Simitall's former analysis functions (second table), tested on
synthetic folders and on Simitall output (`devtools::test()`; the moved tests
run when Simitall is installed). See `AGENTS.md` for how the code is organised.
