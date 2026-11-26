# Bioana: handoff for coding agents (Claude, Codex, ...)

Read this before editing. It is the whole context a new session needs.

## The suite

Nirwan Tandukar (PhD student, plant genomics, NC State) is building four
agents that hand work to each other through files, never through prompts:

1. **Biodes** — github.com/nirwan1265/BioDesignGPT, R package `studyforge`.
   Reads a research question, designs the experiment and writes a **design
   folder**: `study_spec.json`, `field_layout.csv`,
   `sample_metadata_template.tsv` (`write_design_package()`). Schema:
   `inst/registry/study_spec.schema.json` there (handoff 0.5.x).
2. **Simitall** — github.com/nirwan1265/simitall. Fills the design folder
   with simulated data and the truth (`sim_truth.json`, `truth/`), and writes
   `data_manifest.json`. Routes on spec fields via `inst/suite/simulators.tsv`.
3. **Bioana** — this repo, R package `bioana`. Analyses the data in a filled
   design folder, simulated or measured, and scores results against the truth
   when there is one.
4. **Biopaint** — figures (later; some plotting lives in Simitall now).

Rules that hold across the suite:

- Bioana never reads the user's question and never simulates. Its inputs are
  the design folder and the data listed in `data_manifest.json`.
- Decisions live in rule tables (`inst/rules/*.tsv`, `*.json`), not in code:
  adding a data type or a setting should be a table row plus one method.
- Simulated and measured data take the same path. Truth only adds scoring.
- Biodes decides the experimental unit, blocks and random terms
  (`analysis.model`). Bioana fits that model; it refuses to treat
  subsamples as independent replicates.
- Statistics come from published packages (the engines in
  `inst/rules/engines.tsv`, ordered per data type and design in
  `engine_selection.tsv`). Native Bioana code is the last fallback only.
  Every adapter has a gold check in `tests/testthat/test-engines.R`: the
  adapter must give the same numbers as the package called directly.

## Inputs

**Design folder** (from Biodes):

- `study_spec.json`: `handoff_version`, `organism`, `system`, `factors`
  (name, n_levels, levels, applied_to), `design.chosen`, `design.reps`,
  `primary_outcome` (name, `type` = continuous / binary / count / ordinal /
  survival), `assays` (coarse labels) and `assay_kinds` (finer: dnase_seq,
  amplicon_16s, crispr_screen, ...), `population` (mating, generations),
  `analysis.model` (an lme4 formula with generic terms: treatment,
  treatment2, block, rep, row, col, unit, subject, time, env, ...) and
  `analysis.terms` (which layout column each term is).
- `field_layout.csv`: one row per experimental unit, with the factor columns.
- `sample_metadata_template.tsv`: one row per sample (`sample_id`), the layout
  columns, `primary_outcome` (filled with data), `event` for survival.

**`data_manifest.json`** (Simitall writes it; users write it with
`bioana_write_manifest()`): `data` = list of entries with `id`, `data_type`
(a row of `inst/rules/analyses.tsv`), `assay`, `files` (named, relative
paths) and optional `truth` (`kind` features / regions / outcome / qtn,
`file`, `term` = the model term the conditions belong to, `large`).
A `genotypes` entry may carry `genotype_term` (the model term that holds the
lines); without it Bioana reads `sim_truth.json` `genetics.genotype_term`,
then matches `inst/rules/genetics.json` `genotype_pattern`. Without a
manifest, `bioana_manifest()` infers the outcome and `counts.tsv`.

**Data files**: every feature x sample matrix has the feature id in the first
column and one column per `sample_id`. Long tables (ChIP peak counts, ATAC
reads) have a `library_id` column and a libraries table mapping libraries to
`sample_id` (or `sample`).

## Outputs

`<dir>/bioana/<id>/results.tsv`, `<dir>/bioana/outcome/recovery.tsv`,
`<dir>/bioana/bioana_summary.json` (see README). Feature contrasts are
model-matrix column names (`treatmenttreatment_2`).

## Code map

| file | what |
|---|---|
| `R/analyze.R` | `bioana_analyze()`: read design + manifest, run each data set's method from `analyses.tsv`, write the summary, print one line per data set |
| `R/manifest.R` | `bioana_manifest()`, `bioana_write_manifest()`, `bioana_analyses()`, truth loading |
| `R/read_design.R` | reads and checks the Biodes folder, parses `analysis.model` into roles and columns. **A copy of Simitall's `R/suite_handoff.R` (`.suite_` -> `.ba_`): change both together** |
| `R/engines.R` | engine registry and selection (`bioana_engines()`, `bioana_engine_plan()`), feature design, the engine chain, provenance, secondary-engine concordance |
| `R/engine_features.R` | one adapter per feature engine: DESeq2, edgeR, limma-voom, limma, dream, DSS, edgeR IP x input, native |
| `inst/rules/engines.tsv`, `engine_selection.tsv` | the registry (package, calls, model, test, multiple testing, assumptions, reference, tutorial) and the per-data-type, per-design ranking |
| `R/readers.R`, `inst/rules/file_formats.tsv`, `file_formats.json` | real data: format detection (`bioana_detect_format()`), readers (featureCounts, tximport for Salmon/kallisto, STAR, Bismark, VCF, plain tables, MACS peaks), `bioana_manifest_from_files()`, `bioana_consensus_peaks()` |
| `docs/contracts/` | shared contracts proposed to the other agents (`biodes_from_data.md`) |
| `docs/engines/` | one tutorial per engine (the package's example code, then the same through Bioana, then the gold check) |
| `R/linear_test.R` | the shared feature engine: per-feature linear model with the design's fixed terms, empirical-Bayes moderated t (limma's squeezeVar, native), BH; truth scoring by feature or by region overlap |
| `R/method_outcome.R` | primary outcome by type (lm/lmer, glm/glmer, glm.nb, clm/clmm, coxph) and its truth scoring |
| `R/method_features.R` | each feature data type: builds the matrix on a testable scale, then the engine |
| `R/method_genotypes.R` | `genotypes`: genotype term, line means (BLUEs) from the design's model, GWAS, QTN scoring (exact marker or LD) |
| `R/gwas.R`, `R/genomic_selection.R` | GWAS (rrBLUP mixed model, binary traits, truth benchmark, Manhattan/QQ) and genomic prediction / selection (moved from Simitall) |
| `R/rnaseq_eqtl.R`, `R/scrnaseq.R` | eQTL benchmarks, RNA-seq results summary and figure; single-cell clustering, condition tests, cell-type eQTL benchmark and figure (moved) |
| `R/peaks.R` | windowed peak calling, differential peaks, scoring against truth peaks, ChIP results figure (moved) |
| `R/go_annotations.R`, `inst/extdata/gene_coordinates`, `inst/extdata/go` | gene-coordinate catalogues and GO joins (moved) |
| `R/assembly.R`, `R/plot_panels.R`, `R/globals.R` | Unicycler / QUAST evaluation; saving figure panels; ggplot column names (moved) |
| `inst/rules/analyses.tsv` | data_type -> method, status, files needed, assays, note |
| `inst/rules/method_parameters.tsv` | settings per data type (window sizes, filters) |
| `inst/rules/genetics.json` | genotype term pattern (**a copy of Simitall's `inst/suite/genetics.json` pattern: change both**), QTN types scored, cross mating schemes |
| `inst/rules/handoff.json`, `model_terms.json` | Biodes contract (supported handoff versions, model-term vocabulary); copies of Simitall's |
| `tests/testthat/` | `test-analyze.R`: synthetic data in a fixture RCBD folder (no Biodes or Simitall needed); `test-gwas-selection.R`, `test-legacy-analyses.R`: the tests moved with the functions, which simulate with Simitall (skipped without it) |

## Analysis goals (v0.7, 2026-11-25)

Contract with Biodes: [docs/contracts/analysis_goals.md](docs/contracts/analysis_goals.md).
Biodes owns the vocabulary (`studyforge` `inst/registry/analysis_goals.json`) and writes
`analysis.goals` / `primary_goal` / `goal_terms` (+ `assumptions`). Bioana keys
`engine_selection.tsv` `goals` to those ids (no copy of the vocabulary; `bioana_goals()`
compares with the Biodes registry when studyforge is installed). `R/goals.R`: goal
resolution per data set, test plans (Wald coefficients or LRT drop set), outcome LRT.
Implemented: treatment_effect, interaction, time_course, trend, association, survival.
Truth scoring of joint tests is not done yet (Simitall should plant truth for
`primary_goal`; the truth `term` names the goal term).

## Safeguards (2026-11-25)

- `R/design_checks.R` + `inst/rules/design_checks.json`: checks before every data set
  (estimability/confounding, residual df, units per level, repeated subjects without a
  random term, technical replicates, unidentifiable random terms); stop or warn.
  `bioana_check_design(dir)` exposes them.
- `engine_selection.tsv` `design_fit` (full / approx / drops_random) and `engines.tsv`
  `max_random` (0 / 1 / any). `.ba_run_chain()` runs every chain (features, outcome,
  genotypes): after a runtime failure only a `full` engine may take over; warnings
  raised during fitting are collected into `summary$diagnostics`.
- Outcome fits go by engine id (`.ba_fit_outcome()`); new engines `coxme`,
  `coxph_frailty`, `polr`.
- `R/read_design.R` (copy of Simitall's `suite_handoff.R`): random terms now carry
  `intercept_only`. Tell Simitall (same change in `.suite_parse_model`).

## Biodes 0.6.0 (8132e0c), 2026-11-25

- `handoff.json` accepts 0.6. `assumptions` is read as an object keyed by field path.
- `R/read_design.R` (copy of Simitall's `suite_handoff.R`): `parsed$values[[term]]`
  = {levels, values, unit, scale} from `factors[]`; `.ba_frame` maps level labels
  ("10 mg") to values for numeric terms and applies `values[[term]]$transform` (set by
  the trend goal for `scale = "log"`). Tell Simitall (same change in its parser).
- `inst/rules/goal_rules.json`: zero-dose log rule (`offset_min_positive`) and the
  meaning of `goal_effect.scale = "outcome"` (engine scale). Simitall must plant on these.
- Biodes 0.6.2 (4bcb848): trend `log_base`/`log_offset` read from the factor
  (`parsed$values` carries them: same change for Simitall's parser); `time_course`
  `scope = "within_group"` with `group`/`group_factor` runs (group as reference level, time
  main-effect columns tested jointly). Outcome joint tests now fit the reduced model on the
  full model matrix minus the goal's columns (`X0_`), not a reduced formula.
- `goal_recovery` per data set, keyed to manifest `truth.goal` and `sim_truth.json`
  `goal_truth` (contract section 4). Simitall's power sweep reads `goal_recovery$tpr`
  (features); outcome goals give `detected`.
- Open with Simitall (fa37c0c): `.suite_trend_x` puts a zero dose one decade below the
  smallest positive value; the agreed predictor is `log10(x + log_offset)`.
  `goal_recovery$predictor_match` is FALSE until that changes. Fixed in Simitall 847cf84
  (checked at 279be4b: x = 0, 0.301, 1.041, 2.004, predictor_match TRUE).
- Open with Simitall: `simulate_rnaseq_experiment` defaults `n_latent = 2`, `latent_effect_sd = 0.2`,
  not in `simulator_parameters.tsv` nor in `sim_truth.json`. In a 16-sample dose design they made
  132 of 237 DESeq2 calls genes with no planted effect (7 with `n_latent = 0`). Simitall 399c1ff:
  off by default; when on, scores/loadings files are recorded at `counts.hidden_factors` (top level
  of sim_truth.json, not under `assays.<assay>`; Bioana reads both). With n_latent = 2 Bioana attributes
  125 of those 132 to the hidden factors (share 0.25); 7 remain, as with n_latent = 0.
- 2026-11-25, Biodes 38dce99 genotype x drought RNA-seq (4 per cell, 16 libraries, 12 residual df):
  DESeq2 Wald is liberal on null genes even with no goal planted (p<0.05 7.4%, p<0.001 0.67%);
  edgeR QL is calibrated (4.6%, 0.06%). Interaction goal: DESeq2 41 calls (12 unplanted), edgeR QL 0.
  Open: rule for small designs (DESeq2 useT, or edgeR QL first), for Nirwan to decide.
- 2026-11-25 classical ANOVA (0.8.0): `docs/contracts/stat_tests.md`. Asks Biodes for an optional
  `analysis.posthoc` {type, control, adjust}; asks Simitall for optional correlated occasions
  (AR(1) / random slope) recorded in sim_truth.json. Neither blocks: Bioana runs without them.

## Packages

Every package Bioana calls is a row of `inst/rules/requirements.tsv`
(`test-dependencies.R` fails otherwise). Add the row with a new engine or
reader; `bioana_install_dependencies()` and `install_bioana.R` read it.

## Adding an engine

1. A row in `inst/rules/engines.tsv` (package, requires, input kind, calls,
   model, test statistic, multiple testing, random terms, estimate scale,
   assumptions, reference, doc) and its rank in `engine_selection.tsv`.
2. An adapter `.ba_engine_<engine>(inputs, design, alpha, args)` returning
   `.ba_engine_result()` (table: feature_id, contrast, estimate, se,
   statistic, p, padj; contrast = model-matrix column name).
3. A tutorial `docs/engines/<engine>.md` and a gold check in
   `test-engines.R` (adapter = the package called directly, same numbers).

## Adding a data type

1. A row in `inst/rules/analyses.tsv` (data_type, method, status `ready`,
   files, assays, note) and its settings in `method_parameters.tsv`.
2. A method `.ba_method_<name>(entry, design, outdir, alpha, args)` that
   returns `list(summary, paths, results)`; for features, build a features x
   samples matrix and call `.ba_feature_test()`, then `.ba_finish_features()`.
3. In Simitall, the simulator's manifest entry (`R/suite_manifest.R`) with
   the same `data_type` and its truth in the standard format.
4. A test with synthetic data in `tests/testthat/test-analyze.R`.

## Status (v0.2, 2026-11-23)

End to end on one Simitall folder (RCBD, 2 treatments x 5 blocks, human):
outcome CI covered; RNA-seq 9/50 true changes found (FDR 0.10, estimate r
0.96); ChIP 19/20 (FDR 0.17); ATAC 22/30 windows (FDR 0.04); WGBS 80/119
DMR windows (FDR 0.06); m6A 2/40 (5 per group is little power for m6A);
small RNA and microbiome 0 found (the real templates vary a lot between
samples: signs right, r 0.75 and 0.42, but 5 per group is not enough).

## Real test data (2026-11-24)

github.com/nirwan1265/bioana-testdata holds small real data sets (airway
BAMs, DiffBind ChIP chr18 1-20 Mb, exomePeak2 MeRIP BAMs, TASSEL maize 282
panel, rrBLUP source); its README gives sources and citations.
`tests/testthat/test-real-data.R` uses it via `BIOANA_TESTDATA` or a clone
next to `bioana/`, and skips otherwise. Keep data out of this repo:
`inst/extdata/real/` holds only the 380 KB airway counts for
`bioana_example()`.

Known gap: on the DiffBind ChIP data Bioana calls fewer DB sites than
DiffBind (52-96 vs 167-284 of 314). DiffBind subtracts input reads (its folds
are about 1.5x Bioana's; subtracting input reproduces that) and the rest is
not yet explained. Running DiffBind itself (not installable in the build
environment) on the same subset would close it.

## Level order (2026-11-23, tell Simitall)

`R/read_design.R` changed: `.ba_parse_model()` now returns `levels` (the spec
factor's `levels`, reference first) and `.ba_frame()` orders factor levels by
it, then natural order for unlisted levels. Before, levels were always natural
order (so `trt` came before `untrt`). This file is a copy of Simitall's
`R/suite_handoff.R`: the same change is needed there (`.suite_parse_model`,
`.suite_frame`).

## User data: Biodes describes, Bioana analyses (proposed 2026-11-23)

For a user's own data the flow is **Biodes -> Bioana**, no Simitall: Biodes
describes the experiment that was done (unit, blocks, random terms,
`analysis.model`) and writes the usual design folder; Bioana detects data
types from the files and runs the engines. Simitall is optional afterwards
("what effect could this design detect?") and never picks the engine; the
engine rules live only in Bioana's `engine_selection.tsv` (Simitall reaches
them through `bioana_analyze()`).

**Contract for the Biodes window**: [docs/contracts/biodes_from_data.md](docs/contracts/biodes_from_data.md)
(`studyforge::biodes_from_data(description, sample_sheet, ...)`, new
optional `study_spec.origin = "observed"`, measured `primary_outcome` filled
in, `unit_id` for subsamples, random terms for repeated units, confounding in
`unit_audit.json`).

## Engines (v0.4, 2026-11-23)

Since v0.5 `bioana_analyze(allow_fallback = FALSE)` is the default: a data
set with no citable package installed stops and names the packages to
install; Bioana's own code runs only with `allow_fallback = TRUE` or when
named as the engine. Simitall's power sweep gets an error entry (not a crash)
for such a data set; it can pass `allow_fallback = TRUE` if it wants numbers
without the packages.

Gold-standard packages replace the native statistics: DESeq2, edgeR,
limma-voom, limma, dream, DSS, edgeR IP x input, lm/lmer/glm/glmer/glm.nb/
clm/coxph, glmmTMB, rrBLUP, GAPIT. Results tables use `statistic` (not `t`).
`bioana_analyze(secondary = NA)` asks (interactive) or skips; `TRUE` runs the
secondary engine and writes `results_<engine>.tsv` and `concordance`.
Simitall calls `bioana_analyze(quiet = TRUE)`, so it is never asked.
Tested here: DESeq2 1.42, edgeR 4.0, limma 3.58, DSS 2.50, glmmTMB 1.1.8.
Not installable in the build environment (CRAN/Bioconductor blocked):
variancePartition (dream), rrBLUP, GAPIT; their tests skip without them.

## Genotypes (v0.3, 2026-11-23)

`genotypes` runs in `bioana_analyze()`. Line means: the outcome fitted with
`analysis.model` (genotype term fixed; refitted as fixed if random), each
line's fixed-effect prediction averaged over the other fixed factors. Scan:
`analyze_gwas()`, rrBLUP when installed, else the native EMMAX scan (also
gives per-marker effects). `n_pcs = auto`: 0 for crosses, 3 for panels.
Truth: additive and dominance QTN count in tpr; `gxe_slope` rows dropped;
epistatic QTN are not false positives but not in tpr; a hit counts when it is
the QTN marker or has r2 >= `ld_r2` (0.5) with it.

**Contract change for Simitall** (`R/suite_manifest.R`): add
`genotype_term = truth$genetics$genotype_term` to the `genotypes` manifest
entry. Bioana already falls back to `sim_truth.json`, so nothing breaks
before then. Simitall's power sweep can read `analyses$<id>$recovery$tpr`
for genotypes (same name as the feature types).

## Moved from Simitall (2026-11-23)

All of Simitall's analysis code now lives here: GWAS, genomic selection,
peak calling and differential peaks, single-cell clustering and tests, eQTL
benchmarks, RNA-seq results summaries, GO annotation, assembly evaluation and
their results figures, plus the old `bioana_fit_biodes()` (replaced by
`bioana_analyze()`). Internal helpers were renamed `.simitall_*` ->
`.bioana_*`; `export_simitall_plot_panels()` -> `export_bioana_plot_panels()`;
class `simitall_genomic_model` -> `bioana_genomic_model`. Simitall keeps only
simulation (and the helpers its simulators need, e.g. pseudobulk
aggregation), and calls `bioana::bioana_analyze()` for its power sweep and
`simitall_from_question()`.

## Next

The checklist is "Future analyses to check" in README.md. In short: wire the
moved functions into `bioana_analyze()` (eQTL, single cell,
peak calling, genomic prediction), add analyses for the data types Simitall
will simulate (metabolomics, BSA/pool-seq, Ribo-seq, scATAC, Hi-C, degradome,
spatial, CRISPR screens, repertoires, proteomics, reporter assays), package
adapters with native fallbacks, random terms in feature tests, more contrasts,
power per data type and real-data readers.

## Conventions

- Commit as `Nirwan Tandukar <tandukar.nirwan@gmail.com>`, no AI attribution
  lines, push straight to `main`, no branches.
- Run `devtools::test()` before pushing; keep tests synthetic and small.
- Explain results in plain words; Nirwan wants working code and short
  answers.
