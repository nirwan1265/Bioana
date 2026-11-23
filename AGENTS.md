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
`file`, `term` = the model term the conditions belong to, `large`). Without a
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
| `R/linear_test.R` | the shared feature engine: per-feature linear model with the design's fixed terms, empirical-Bayes moderated t (limma's squeezeVar, native), BH; truth scoring by feature or by region overlap |
| `R/method_outcome.R` | primary outcome by type (lm/lmer, glm/glmer, glm.nb, clm/clmm, coxph) and its truth scoring |
| `R/method_features.R` | each feature data type: builds the matrix on a testable scale, then the engine |
| `inst/rules/analyses.tsv` | data_type -> method, status, files needed, assays, note |
| `inst/rules/method_parameters.tsv` | settings per data type (window sizes, filters) |
| `inst/rules/handoff.json`, `model_terms.json` | Biodes contract (supported handoff versions, model-term vocabulary); copies of Simitall's |
| `tests/testthat/` | synthetic data in a fixture RCBD folder (no Biodes or Simitall needed) |

## Adding a data type

1. A row in `inst/rules/analyses.tsv` (data_type, method, status `ready`,
   files, assays, note) and its settings in `method_parameters.tsv`.
2. A method `.ba_method_<name>(entry, design, outdir, alpha, args)` that
   returns `list(summary, paths, results)`; for features, build a features x
   samples matrix and call `.ba_feature_test()`, then `.ba_finish_features()`.
3. In Simitall, the simulator's manifest entry (`R/suite_manifest.R`) with
   the same `data_type` and its truth in the standard format.
4. A test with synthetic data in `tests/testthat/test-analyze.R`.

## Status (v0.1, 2026-11-23)

End to end on one Simitall folder (RCBD, 2 treatments x 5 blocks, human):
outcome CI covered; RNA-seq 9/50 true changes found (FDR 0.10, estimate r
0.96); ChIP 19/20 (FDR 0.17); ATAC 22/30 windows (FDR 0.04); WGBS 80/119
DMR windows (FDR 0.06); m6A 2/40 (5 per group is little power for m6A);
small RNA and microbiome 0 found (the real templates vary a lot between
samples: signs right, r 0.75 and 0.42, but 5 per group is not enough).

## Next

- `genotypes`: GWAS and QTL mapping on the line genotypes, scored against
  `genetic_truth/qtn_truth.tsv`. Simitall still has `analyze_gwas()`,
  genomic selection, ChIP and single-cell result code (~3,000 lines) to move
  here in stages.
- Single-cell at cell level (clusters, proportions).
- Package adapters, used when installed: DESeq2 / edgeR / limma-voom for
  counts, DSS / methylKit for methylation, ANCOM-BC for microbiome
  (compositional), exomePeak2 for m6A. The native methods stay as the
  fallback.
- Analyses for the data types Simitall will add: Ribo-seq, scATAC, Hi-C,
  MeDIP (already `ip_input_counts`), metabolomics/lipidomics, BSA/pool-seq.
- Simitall's power sweep and `simitall_from_question()` still use
  Simitall's own copy of the old `bioana_fit_biodes()`; switch them to
  `bioana::bioana_analyze()` once both are installed together.

## Conventions

- Commit as `Nirwan Tandukar <tandukar.nirwan@gmail.com>`, no AI attribution
  lines, push straight to `main`, no branches.
- Run `devtools::test()` before pushing; keep tests synthetic and small.
- Explain results in plain words; Nirwan wants working code and short
  answers.
