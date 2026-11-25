# Contract: analysis goals (Biodes writes, Bioana reads)

Status: agreed with the Biodes window, 2026-11-25. Owner of the vocabulary: Biodes
(BioDesignGPT / `studyforge`). Bioana keys its engine rules to the goal ids and keeps no
copy of the vocabulary. Simitall plants truth for `primary_goal`.

## 1. The registry: `studyforge` `inst/registry/analysis_goals.json`

One entry per goal. Bioana reads it from the installed package
(`system.file("registry", "analysis_goals.json", package = "studyforge")`) to check a
design against the goal. Fields Bioana uses:

| field | meaning |
|---|---|
| `id` | the goal id written into `study_spec.json` |
| `test` | `wald` (coefficients one at a time), `lrt` (a joint test: full vs reduced model), `wald_or_lrt` (Wald for one coefficient, LRT for several), or `other` (not a coefficient test: diversity, prediction, ranking) |
| `required_terms` | what `analysis.model` must contain, as goal-term slots: e.g. `["factor_a", "factor_b"]` for interaction, `["time"]` for time_course, `["numeric"]` for trend, `[]` for treatment_effect |
| `data_types` | data kinds the goal applies to (assay labels or Bioana data types; Bioana ignores unknown ones) |
| `description` | one line, shown to the user |

Proposed entries (Biodes owns the final wording):

| id | test | required_terms | Bioana tests |
|---|---|---|---|
| `treatment_effect` | wald | — | each treatment coefficient vs the reference |
| `interaction` | wald_or_lrt | factor_a, factor_b (time allowed) | the interaction coefficients of the two `goal_terms`: Wald if one, joint (LRT / F) if several |
| `time_course` | lrt | time | any change over time: full model vs the model without every term that contains time |
| `trend` | wald | numeric | the slope of a numeric term (dose, gradient) |
| `transcript_usage` | lrt | — | (Bioana: later; DRIMSeq / satuRn) |
| `differential_abundance` | wald | — | (Bioana: later; ANCOM-BC2 / ALDEx2) |
| `community_composition` | other | — | (later; diversity, PERMANOVA) |
| `cell_composition` | other | — | (later; propeller / miloR) |
| `association` | wald | — | GWAS on line means (ready) |
| `qtl_mapping` | other | — | (later; qtl2) |
| `genomic_prediction` | other | — | (later) |
| `variance_components` | other | — | (later; heritability, BLUPs, G×E stability) |
| `ranking` | other | — | (later) |
| `survival` | wald | — | the survival outcome engines (ready: coxph, coxme) |

## 2. The study spec: `analysis.goals`, `analysis.primary_goal`, `analysis.goal_terms`

```json
"analysis": {
  "model": "y ~ treatment * time + (1 | subject)",
  "terms": {"treatment": "genotype", "time": "hours", "subject": "plant_id"},
  "goals": ["interaction", "treatment_effect", "time_course"],
  "primary_goal": "interaction",
  "goal_terms": {"interaction": ["treatment", "time"], "time_course": ["time"]},
  "goal_effect": {"interaction": {"size": 1.5, "scale": "outcome",
                                  "definition": "difference of differences"}}
},
"assumptions": {"analysis.goals": "assumed: interaction (the question asks whether the response over time differs between genotypes)"}
```

Schema 0.6.x (Biodes 8132e0c; 0.6.2 at 4bcb848; 0.6.3 at 2e747fe). Bioana accepts `handoff_version` 0.6.

- `goals`: every goal the study answers, most important first. `primary_goal` is one of them;
  Biodes sizes the experiment for it (an interaction needs about 4× the replicates of a
  main effect) and Simitall plants truth for it.
- `goal_terms`: per goal id, the **generic model terms** (the names in `analysis.model`, as
  `analysis.terms` maps them to columns) the goal tests. Required for `interaction` (two
  terms), `time_course` (the time term, default `time`) and `trend` (the numeric term).
- The model must contain what the goal tests: `interaction` needs the interaction term
  (`treatment * time` or `treatment:time`); `time_course` needs a time term; `trend` needs the
  numeric term as a main effect. Its column must hold numbers (Bioana fits it as numeric for
  this goal only).
- `assumptions`: an object keyed by field path. A goal Biodes chose is the key
  `analysis.goals` (or `analysis.primary_goal`); the value says why. Bioana reads the key,
  prints "goal assumed by Biodes" and keeps the reason in `summary$...$goal$assumed_why`.
  (An older array of `{field, ...}` entries is still read.)
- `goal_effect`: per goal, `{size, scale, definition}`; `time_course` adds `shape`,
  `peak`/`half_time` and `scope`. Bioana reads `scale = "outcome"` as **the scale the engine
  estimates on** (`engines.tsv` `estimate_scale`): outcome units for continuous outcomes, log2
  fold change for counts, log odds for binary, M-value for methylation. Bioana reports it as
  `effect_scale_used`. A trend size is present only when stated (per log10 unit when
  log-scaled); otherwise it is absent and listed in `assumptions`. `goal_effect.trend.per_unit`
  ("kg N/ha", "log10 uM") is reported with the slope. The factor's `scale` decides log or not;
  a `per_unit` that disagrees with it is flagged in the report.
- `time_course` with `scope = "within_group"`, `group` (the exact level label, e.g. `"mutant"`; Biodes
  drops `group`/`group_factor` when no level matches) and `group_factor`
  (e.g. `"genotype"`, a spec factor or model term): Bioana makes `group` the reference level
  of that term and jointly tests the time main-effect coefficients, i.e. time within that
  group (full model vs the model without that group's time effects). It needs the
  `group_factor:time` interaction in the model; without it, or when the group is not a level,
  time is tested over all groups and the report says why.
- Trend factors: `factors[].values` (numbers aligned with `levels`, e.g. levels `"0 mg"`,
  values `0`), `unit`, `scale` (`linear` / `log`). Bioana maps level labels in the layout to
  `values`, and with `scale = "log"` fits the slope on `log_base(value + log_offset)`, both from
  the factor (Biodes 0.6.2: `log_base: 10`, `log_offset` = smallest positive level when there is
  a zero). Without them, `inst/rules/goal_rules.json` gives base 10 and the smallest positive
  value. The report gives the slope's unit. **Simitall must plant the slope on the
  same scale.**
- The question Biodes asks when the text does not decide, e.g. for "RNA-seq of wild-type and
  mutant roots at 5 times": *Genes that differ between genotypes, genes whose response over
  time differs between genotypes, or genes that change over time at all?* (`treatment_effect`,
  `interaction`, `time_course`).
- Biodes' existing top-level `goal` ("analysis" / "design" intent) is unrelated and stays.

## 3. What Bioana does

For each data set Bioana takes the first goal in `analysis.goals` that
`engine_selection.tsv` has rows for (its `goal` column) for that data type. With no
`analysis.goals`, it uses `treatment_effect` (`association` for genotypes) and says so in
the report when the model has a time or interaction term. For the goal it builds a test
plan (the coefficients for a Wald test, or the reduced model for a joint test), picks the
first installed engine that can run that test (`engines.tsv` `tests`), and records the goal,
test and plan in `summary$analyses$<id>$goal`.

Joint tests per engine: DESeq2 `test = "LRT"` with the reduced model matrix; edgeR
`glmQLFTest` on the coefficient set (QL F); limma, voom and dream `topTable` on the
coefficient set (moderated F); outcome engines `anova(reduced, full)`. Engines without a
joint test (DSS, the IP × input model, Bioana's own code) are skipped for `lrt` goals.

A goal id Bioana does not know yet stops that data set with the list it supports. A goal
id that is not in the Biodes registry (when `studyforge` is installed) is reported as unknown
to Biodes.

## 4. Goal recovery (simulated data)

Keyed to each data set's manifest `truth.goal` and `sim_truth.json` `goal_truth`
(`goal`, `model_terms`, `status`). `summary$analyses$<id>$goal_recovery`:

- every data set: `goal`, `status`, `model_terms` (Simitall), `tested_terms` (Bioana),
  `terms_match` (same variable sets), `scored` (planted and same goal), `within` /
  `group_match` (within-group time course vs `groups_changing`), `predictor_match` (trend:
  Bioana's `log_base(values + log_offset)` against `goal_truth.x`; a mismatch stops the slope
  comparison and is reported).
- outcome: `p`, `detected`; for a one-coefficient goal `estimate`, `true_effect`, `bias`,
  `covered` (continuous) or `sign_agreement`.
- features (`truth.goal.file`: feature_id, term, effect, change): `tpr` (share of goal features
  found by the goal test), `found`, `n_goal_features`, `n_tested`, `n_called`,
  `called_outside_goal`; trend also `estimate_cor`, `sign_agreement`.

## 5. For Simitall

Plant truth for `primary_goal` on the scales in section 2 (outcome scale per data type; log
trend with the zero-dose offset): an interaction effect (the difference of differences) for
`interaction`, a temporal pattern for `time_course`, a slope for `trend`. The truth file's
`term` names the goal term (`treatment:time`, `time`, the numeric term) so Bioana can score it.
