# bioana 0.4.0

- Gold-standard engines. `inst/rules/engines.tsv` (registry) and
  `inst/rules/engine_selection.tsv` (order per data type and design) decide
  which package fits each data set: DESeq2, edgeR, limma-voom, dream, limma,
  DSS, an edgeR IP x input interaction model, lm/lmer/glm/glmer/glm.nb/clm/
  coxph, glmmTMB, rrBLUP, GAPIT. Native Bioana code is the fallback only and
  is flagged.
- Every data set records an `engine` provenance block (package, version,
  exact call, formula, model, input, assumptions, test statistic, multiple
  testing, experimental unit, fallback used and why).
- `bioana_analyze(secondary = )`: run the secondary engine and report
  concordance with the primary; asks in interactive sessions by default.
- `bioana_engines()`, `bioana_engine_plan()`; `params = list(<data_type> =
  list(engine = "..."))` picks an engine.
- Random terms in feature tests: dream, or limma with duplicateCorrelation.
- Feature results: column `t` renamed `statistic` (Wald z, QL F or moderated t).
- `docs/engines/`: a tutorial per engine; `tests/testthat/test-engines.R`:
  gold checks (adapter = package called directly).
- `genotypes`: the `backend` parameter is replaced by the engine rules;
  `analyze_gwas(backend = "GAPIT")` added.

# bioana 0.3.0

- `genotypes` is ready in `bioana_analyze()`: line means (BLUEs) of the
  primary outcome per level of the genotype term, from the design's
  `analysis.model`, then a kinship mixed-model GWAS on the VCF. Simulated
  folders are scored against `qtn_truth.tsv` (tpr, fdr, sign, estimate r) and
  the line means against the true genetic values. Writes
  `bioana/<id>/results.tsv`, `line_means.tsv`, `qtn_recovery.tsv`.
- `analyze_gwas(backend = )`: `"native"` EMMAX scan (same model as rrBLUP with
  `P3D = TRUE`) when rrBLUP is not installed; it also returns per-marker `beta`
  and `se`. The `.gwas.tsv` table gains `beta` and `se` columns (NA with rrBLUP).
- `benchmark_gwas(ld_r2 = , genotypes = )`: a hit in LD with a causal marker
  counts (for crosses); `ld_r2` is added to the metrics table.
- New rules: `inst/rules/genetics.json` (genotype term pattern, QTN types
  scored, cross mating schemes) and `genotypes` rows in `method_parameters.tsv`.

# bioana 0.2.0

- All of Simitall's analysis functions moved here: GWAS (`analyze_gwas()`,
  `analyze_binary_gwas()`, `benchmark_gwas()`, `plot_gwas_results()`), genomic
  prediction and selection, peak calling and differential peaks, single-cell
  clustering and condition tests, eQTL benchmarks, RNA-seq results summary,
  GO annotation with the gene-coordinate catalogues, and assembly evaluation.
- `bioana_analyze()` returns the truth tables (`recovery`) and fits, and takes
  `quiet`.

# bioana 0.1.0

- `bioana_analyze()`: every data set in a Biodes design folder, by data type
  from `data_manifest.json`, with the design's model; scored against the truth
  when simulated.
