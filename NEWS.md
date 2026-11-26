# bioana 0.7.6

- Per-level truth (Simitall `truth.file`) is scored only when the goal is `treatment_effect`; under interaction,
  time_course or trend the goal is scored by `goal_recovery` (no more "found 0/200" for an untested contrast).
- With an interaction in the model, the report says each treatment coefficient is the effect at the reference level
  of the other factor (`recovery_note`).
- Checked with Biodes 38dce99 level names: two-level factors use the first level as reference
  (`treatmentdrought:treatment2mutant`).

# bioana 0.7.5

- Simitall hidden factors (`sim_truth.json` `counts.hidden_factors`, or under `assays.<assay>`): unplanted calls
  whose hidden part (loadings x scores, fitted with Bioana's design model) moves them along the tested effect are
  reported as `called_hidden_factor`, apart from `called_not_planted`. Rule in `goal_rules.json`.
- Checked with Simitall 399c1ff (hidden factors off by default): dose design, 7 of 109 DESeq2 calls unplanted (6%).

# bioana 0.7.4

- Trend goals on simulated data: the per-level truth check is skipped (it crashed on a numeric term:
  "0 (non-NA) cases"); the slope is scored by `goal_recovery`. Same for features: no per-level "found 0/100".
- `goal_recovery` for features splits calls outside the goal set into `called_other_planted` (Simitall's
  other effects) and `called_not_planted`, with `false_call_rate`.
- Checked against Simitall 279be4b (trend predictor fixed in 847cf84): `predictor_match = TRUE`.

# bioana 0.7.3

- Models with no fixed terms (`y ~ (1 | entry)`, p-rep catalog) no longer fail: every fixed-part formula is
  `reformulate(c("1", fixed_labels))` (design checks, outcome, features, linear test).
- Biodes 0.6.3: `goal_effect.trend.per_unit` reported with the slope; the factor's `scale` decides log, and a
  disagreeing `per_unit` is flagged. Exact level labels ("wild type", "mutant") for within-group time courses.

# bioana 0.7.2

- Biodes 0.6.2: trend on `log_base(x + log_offset)` from the spec factor (defaults in `goal_rules.json`:
  base 10, smallest positive value); trend size may be absent.
- `time_course` within one group (`goal_effect.time_course` `scope`, `group`, `group_factor`): that group is
  the reference level and its time effects are tested jointly. Equals `anova()` and DESeq2's LRT on the same
  matrices (`test-goals.R`).
- Outcome joint tests: reduced model = full model matrix without the goal's columns.
- `goal_recovery` (simulated data), keyed to manifest `truth.goal` and `sim_truth.json` `goal_truth`: terms
  match, detection, estimate vs true effect, trend predictor check; for features the share of goal features
  found (`tpr`). Joint goals on features are now scored.

# bioana 0.7.1

- Reads Biodes 0.6.0 specs: `handoff_version` 0.6; `assumptions` keyed by field path (reason kept);
  `goal_effect` reported with the engine's scale; trend factors with `values` / `unit` / `scale`
  (level labels mapped to values, log scale with a zero-dose offset from `inst/rules/goal_rules.json`).

# bioana 0.7.0

- Analysis goals (contract with Biodes: `docs/contracts/analysis_goals.md`). `engine_selection.tsv`
  is keyed by goal; `engines.tsv` says which tests each engine runs (`wald`, `lrt`).
  `treatment_effect`, `interaction` (Wald or joint), `time_course` (joint test without the time
  terms), `trend` (numeric slope); goals read from `analysis.goals` / `primary_goal` /
  `goal_terms`, assumed goals flagged, default with a note when the model has time or interaction
  terms. `bioana_goals()`; `params$<id>$goal`.
- Joint tests: DESeq2 LRT with a reduced model matrix, edgeR QL F and limma / voom / dream moderated
  F on a coefficient set, outcome `anova(reduced, full)` (written to `goal_test.tsv`). Each equals the
  package called directly (`test-goals.R`).

# bioana 0.6.0

- Design checks before every analysis (`bioana_check_design()`, `inst/rules/design_checks.json`):
  stop on a confounded treatment, no residual df, repeated subjects without a random term,
  technical replicates, unidentifiable random terms, normalised values given as counts, outcome
  values that do not match the type; warn on thin replication and few random-term levels.
- No silent substitution: `design_fit` (full / approx / drops_random) per engine row and
  `max_random` per engine. An engine that cannot hold the random terms is skipped; after a runtime
  failure only an engine that fits the design fully takes over; approximations are marked in the
  report and provenance; dropping random terms needs `allow_fallback = TRUE`.
- Fit diagnostics: warnings raised while fitting, singular fits, DESeq2 non-converged genes,
  glmmTMB Hessian, in `summary$diagnostics` and the report.
- Survival outcomes with random terms: `coxme`, or a gamma frailty Cox model as an approximation;
  `polr` is an engine (fixed designs).

# bioana 0.5.0

- `inst/rules/requirements.tsv`: every package Bioana uses (source, profile, what for);
  `bioana_check_dependencies()`, `bioana_install_dependencies()` (profiles minimal, analysis, full,
  suite), and `install_bioana.R` at the top of the repo.

- Real data: `bioana_example("airway")` (the DESeq2 vignette data set, shipped
  in `inst/extdata/real/airway/`). Bioana's DESeq2 run equals the vignette code
  (log2FC to 5e-9) and finds the known dexamethasone genes
  (`test-real-airway.R`).
- BAM input: `bioana_manifest_from_files(annotation = , libraries = )` counts reads with Rsubread
  featureCounts in genes (GTF), peaks (SAF or peak files) or exon windows (IP/input); paired-end detected
  without an index; single-end ChIP reads extended to 125 bp (`file_formats.json` `chip_fragment_bp`).
  HapMap genotypes (converted to VCF) and DSS count text are read too.
- Sample names also match when a sample_id is the only prefix of a file name (`S01_subset.bam`).
  Appending to a manifest keeps existing ids; reading the same files again replaces their entry.
- Fix: bisulfite windows treat positions as 1-based (window k covers (k-1)*w+1 .. k*w).
- `test-real-data.R`: airway BAMs, DiffBind ChIP, maize 282 GWAS, DSS, bsseq Bismark files, exomePeak2
  MeRIP, from github.com/nirwan1265/bioana-testdata.
- Factor levels follow `study_spec.json` `factors.levels` (reference first)
  instead of alphabetical order; levels the spec does not list come after.

- Only published packages run by default: `bioana_analyze(allow_fallback =
  FALSE)`. A data set without an installed engine stops and names the
  packages to install; Bioana's own code needs `allow_fallback = TRUE` (or
  naming it as the engine) and stays marked `FALLBACK`.
- Real data: `bioana_manifest_from_files()` detects formats
  (`inst/rules/file_formats.tsv`), converts them and writes the manifest with
  the source recorded: featureCounts, Salmon and kallisto through tximport,
  STAR ReadsPerGene (strandedness detected), Bismark coverage and cytosine
  reports, VCF, plain count / beta / IP-input tables. `bioana_detect_format()`,
  `bioana_consensus_peaks()` (MACS peaks -> SAF for featureCounts).
- Fix: DSS (bisulfite) and the edgeR IP x input model were never chosen inside
  `bioana_analyze()` (the chain compared input kinds to input names), so those
  data sets went to limma. Tests now check the engine that runs end to end.
- `docs/contracts/biodes_from_data.md`: the Biodes "describe my data" contract.

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
