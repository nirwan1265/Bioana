# Analysis engines: which package Bioana calls, and how to check it

Bioana does not invent statistics. For every data set it calls a published
package (the **engine**), with the model the design asks for, and records what
it did so you can rerun it yourself.

- `inst/rules/engines.tsv`: one row per engine. It gives the package, the
  function called, the model, the test statistic, the multiple-testing method,
  the estimate scale, the assumptions, the reference, and the tutorial page.
  Print it with `bioana_engines()`, which also shows what is installed.
- `inst/rules/engine_selection.tsv`: which engine runs for each data type and
  design, in rank order. **`fixed`** means `analysis.model` has no random terms
  (independent units). **`random`** means it has random terms (the same subject,
  plot or animal measured more than once, or nested units). Biodes decides this
  when it writes the model, so Bioana never has to ask "are these the same
  individuals over time?". `bioana_engine_plan("feature_counts", "random")`
  shows the plan.
- **Primary** = the first installed engine in that order. **Secondary** = the
  next gold-standard engine. `bioana_analyze(dir, secondary = TRUE)` runs both
  and reports how well they agree; with `secondary = NA` (default) Bioana asks
  in an interactive session. **Native Bioana code is fallback only**: the
  report marks it `FALLBACK`, and the summary has `fallback_used: true` and
  says why the gold engines were not used.
- Pick an engine yourself:
  `bioana_analyze(dir, params = list(feature_counts = list(engine = "edger")))`.

## The provenance block

Every data set in `bioana/bioana_summary.json` has an `engine` block with:
`package`, `version`, `call` (the exact call), `formula`, `model`, `input` (what
went in, with dimensions), `assumptions`, `test_statistic`, `multiple_testing`,
`estimate_scale`, `random_terms` (what was done with them), `experimental_unit`,
`n_units`, `fallback_used`, `not_used` (engines skipped and why), `reference`,
`tutorial`.

## Gold-standard checks

`tests/testthat/test-engines.R` fits the same data twice: once through
Bioana's adapter, and once by calling the package directly, as its vignette
does. The test requires the **same numbers** (estimates, statistics,
p values, adjusted p values; tolerance 1e-6 or tighter). The native fallback
is checked against limma: it gives identical estimates, and its t statistics
are within 5%.

| engine | package | data types (rank 1 when) | tutorial | gold check |
|---|---|---|---|---|
| `deseq2` | DESeq2 | counts, fixed design | [deseq2.md](deseq2.md) | `DESeq2 adapter = DESeq2 called directly` |
| `edger` | edgeR | peaks and windows, fixed design; counts (rank 2) | [edger.md](edger.md) | `edgeR adapter = edgeR QL pipeline called directly` |
| `edger_ip_input` | edgeR | m6A / MeDIP / RIP (IP and input), fixed | [edger.md](edger.md) | `edgeR IP x input: enrichment changes are found, expression changes are not` |
| `limma_voom` | limma | counts (rank 3); counts with one random term (rank 2) | [limma.md](limma.md) | `limma-voom adapter = ...`, `limma with one random term = duplicateCorrelation ...` |
| `limma` | limma | methylation arrays, M-values, log ratios | [limma.md](limma.md) | `limma on values = lmFit/eBayes` |
| `dream` | variancePartition | counts, peaks, windows with random terms | [dream.md](dream.md) | `dream adapter = variancePartition::dream called directly` |
| `dss` | DSS | bisulfite counts, fixed design | [dss.md](dss.md) | `DSS adapter = DMLfit.multiFactor / DMLtest.multiFactor called directly` |
| `stats_lm`, `lmer`, `glm_*`, `glmer_*`, `glm_nb`, `glmer_nb`, `clm`, `coxph`, `glmmtmb` | stats, lme4/lmerTest, MASS, ordinal, survival, glmmTMB | primary outcome by type | [outcome.md](outcome.md) | `glmmTMB outcome engine matches glmmTMB called directly, and agrees with lm` |
| `rrblup`, `gapit_mlm`, `bioana_emmax` | rrBLUP, GAPIT, Bioana | genotypes (GWAS on line means) | [gwas.md](gwas.md) | `native and rrBLUP scans agree` |
| `bioana_native` | Bioana | fallback for every feature type | [native.md](native.md) | `the native fallback reproduces limma` |

## Installing the engines

```r
install.packages("BiocManager")
BiocManager::install(c("DESeq2", "edgeR", "limma", "variancePartition", "DSS", "bsseq"))
install.packages(c("lme4", "lmerTest", "glmmTMB", "MASS", "ordinal", "survival", "rrBLUP"))
remotes::install_github("jiabowang/GAPIT")   # GAPIT 3 is distributed on GitHub
```

`bioana_engines()[, c("engine", "package", "installed")]` shows what is missing.

## Coming next (same pattern: a registry row, a selection row, an adapter, a tutorial, a gold check)

qtl2 (QTL mapping in crosses), muscat (single-cell pseudobulk), csaw and
DiffBind (peaks), ANCOM-BC2 and ALDEx2 (microbiome), exomePeak2 (m6A from
BAM), MSstats and proDA (proteomics), clusterProfiler and fgsea (enrichment),
mixOmics (multi-omics). Python engines (Scanpy/scVI-tools, PyDESeq2,
statsmodels, scikit-learn, COBRApy) come after the R engines.
