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

Needs `jsonlite` and `lme4`; `lmerTest`, `MASS`, `ordinal` and `survival` are
used for mixed models and non-Gaussian outcomes when installed.

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
| `genotypes` | vcf | GWAS, QTL mapping | planned |
| `single_cell_counts` | counts, cells | scRNA-seq (cell level) | planned |

Settings per data type: `inst/rules/method_parameters.tsv`, overridden with
`bioana_analyze(dir, params = list(bisulfite_counts = list(window_bp = 500)))`.

## Output

`<folder>/bioana/`:

- `<id>/results.tsv` — outcome: one row per coefficient (`term, estimate, se, df, statistic, p, role, padj`);
  features: one row per feature and treatment contrast (`feature_id, contrast, estimate, se, t, p, padj`).
- `<id>/recovery.tsv` — outcome only, simulated data: true vs estimated coefficients.
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

## Status

Version 0.1: the analyses in the table above, tested on synthetic folders and
end to end on Simitall output. See `AGENTS.md` for how the code is organised
and what comes next.
