# Contract: Biodes "describe my data" (for the Biodes window)

Status: proposed by Bioana, 2026-11-23. Owner: Biodes (BioDesignGPT /
`studyforge`). Bioana reads the result. Simitall can read it later for "what
could this design detect?".

## Why

Bioana should analyse a user's own data the same way it analyses a simulated
experiment: from a design folder. Today only `write_design_package()` makes
that folder, and it designs a new experiment from a question. For real data
the experiment already exists. Biodes has to **describe** it rather than
design it, but write the same folder. Biodes keeps every design decision: the
experimental unit, the blocks, the random terms and `analysis.model`. Bioana
fits that model as given and picks the package from its
`engine_selection.tsv`. A random term in `analysis.model` makes Bioana use
dream or limma duplicateCorrelation instead of DESeq2, so getting the random
terms right matters most.

## Function

```r
studyforge::biodes_from_data(
  description,           # the user's words: organism, treatments, times, what was measured, how sampled
  sample_sheet,          # data.frame or path: one row per sample (library / measurement), sample_id + any columns
  data_files = NULL,     # optional file names; passed through unread (Bioana detects data types)
  out_dir,               # the design folder to write
  ask = interactive()    # ask the user when the description and sheet do not decide something
)
```

It returns `out_dir` invisibly and writes the same three files as
`write_design_package()`.

## What Bioana needs in the folder

**`study_spec.json`** (`handoff_version` 0.5.x is fine; Bioana ignores fields
it does not know):

| field | content for observed data |
|---|---|
| `origin` | **new, optional**: `"observed"` (from data) or `"designed"` (default when absent) |
| `organism`, `system` | from the description |
| `factors` | one per treatment factor: `name`, `n_levels`, `levels` = the levels **in the sample sheet**, reference (control) first |
| `design.chosen` | what the layout is: `crd`, `rcbd`, `split_plot`, `repeated_measures`, `factorial`, ... |
| `design.reps` | the smallest number of experimental units per treatment cell actually present |
| `experimental_unit` | the unit that received the treatment (plant, plot, pot, animal, flask, ...) |
| `primary_outcome` | `name`, `type` (continuous, binary, count, ordinal, survival) if a phenotype column was measured; omit if the data are omics only |
| `assays`, `assay_kinds` | from the description and file names (`rnaseq`, `wgbs`, `merip_seq`, `gwas`, ...) |
| `population`, `genetics` | for crosses and panels, as today (`mating`, `unit = line`) |
| `analysis.model` | the lme4 formula with generic terms, **random terms for repeated or nested units** |
| `analysis.terms` | which sample-sheet column each generic term is (`treatment = "drought"`, `time = "hours"`, `subject = "plant_id"`) |
| `power` | omit, or `solve_for = "detectable_effect"` (Simitall's post-hoc question) |

**`field_layout.csv`**: one row per **observed** experimental unit (unique
units from the sheet, never invented ones), with the columns `analysis.terms`
names.

**`sample_metadata_template.tsv`**: the sample sheet. `sample_id` must equal
the column names in the data files. It carries the layout columns,
`primary_outcome` **filled with the measured values** (not empty), `event` for
survival, and `unit_id` when several rows belong to one unit.

**`unit_audit.json`** (optional; Bioana prints its warnings): observed
replicates per cell, missing cells, confounding found, subsamples found.

## Rules Biodes applies

1. **Repeated measures.** If the same unit id appears at several times or
   conditions, add `(1 | subject)` (or the right nesting) to `analysis.model`.
   If the sheet has no unit id and the description does not say, ask:
   *"Were the samples at each time point taken from the same plants, or from
   different plants?"* The answer decides the random term.
2. **Subsamples are not replicates.** Several rows from one unit under the
   same condition (technical replicates, several leaves of one plant) get the
   unit's `unit_id`. The model stays at the unit level; Bioana averages or
   nests subsamples and never counts them as replicates.
3. **Blocks.** A batch, block, run or lane column that is crossed with
   treatment becomes `block` (fixed, or `(1 | block)` with many levels, as
   Biodes already decides).
4. **Confounding.** If treatment is confounded with batch (every control in
   run 1, every treated sample in run 2), write it in `unit_audit.json` and
   in the reply to the user. Do not fit a model that cannot separate them.
5. **Reference level.** Control or untreated or 0 h comes first in
   `factors.levels`. Ask if none is obvious.
6. **No more than three questions.** Ask only what changes `analysis.model`
   or the unit, and say in one line what each answer decides.

## Acceptance examples (Bioana will test these once the function exists)

| description | sheet | expected `analysis.model` |
|---|---|---|
| maize leaves, drought vs control, 0/6/12/24/48 h, 3 plants per time, **different plants each time** | 30 rows, no plant id | `y ~ treatment * time` |
| same, **the same 3 plants per treatment sampled at every time** | 30 rows, `plant_id` repeats | `y ~ treatment * time + (1 \| subject)`, `subject = plant_id` |
| sorghum field, 4 N rates, RCBD with 5 blocks, 3 leaves per plot sequenced | 60 rows, `plot` repeats 3x | `y ~ treatment + block`, `unit_id = plot` (subsamples) |
| 200 inbred lines x 2 reps, plant height, genotyped by GBS | 400 rows + VCF | `y ~ treatment + block`, `treatment = line` (genotype term), `genetics.unit = line` |

## What Bioana does with it

`bioana_manifest_from_files(dir, files)` detects each data type and writes
`data_manifest.json`. `bioana_analyze(dir)` then fits `analysis.model` with
the package `engine_selection.tsv` names, and records package, version, call
and citation. With `origin = "observed"` there is no truth, so no recovery
scores.
