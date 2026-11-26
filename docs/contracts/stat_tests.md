# Classical ANOVA and post-hoc tests: what Bioana reads (2026-11-25)

Bioana runs a classical ANOVA next to the model fit for every continuous outcome
(`docs/engines/anova.md`). It needs nothing new to run; these fields make it better.

## Biodes (optional, new)

```json
"analysis": {
  "posthoc": {"type": "vs_control", "control": "untreated", "adjust": "dunnettx"}
}
```

- `type`: `pairwise` (all pairs; default), `vs_control` (each level against one control), `none`.
  Fill it when the question says it: "compared with the control / untreated / wild type" ->
  `vs_control` with `control` = that level label; "which genotypes differ" -> `pairwise`.
- `control`: the exact level label (default: the first level, the reference).
- `adjust`: an emmeans method (`holm`, `bonferroni`, `tukey`, `dunnettx`, `none`). Leave it out
  to use Bioana's default (`inst/rules/posthoc.json`: Holm; Dunnett for vs_control).
- Already used, keep as is: the subject as a `unit`/`subject` random intercept
  (`y ~ treatment * time + (1 | unit)`), and `repeated_factors` (Bioana checks it against the layout).

## Simitall (optional, new)

- Repeated measures are now simulated with a random intercept only. That gives compound
  symmetry, so sphericity always holds and the Greenhouse-Geisser correction is never exercised.
  Please add an option for correlated occasions, e.g. AR(1) residuals over time or a random slope,
  and record it in `sim_truth.json`, e.g.
  `"within_subject": {"structure": "ar1", "rho": 0.7, "over": "time"}`.
- No new truth is needed for the post-hoc tests: Bioana can derive true contrasts from
  `treatment_effects` / `interaction_effects`.

## Files Bioana writes (for Biopaint)

`bioana/<id>/anova_table.tsv` (term, num_df, den_df, MSE, F, ges, pes, p, correction) and
`bioana/<id>/posthoc.tsv` (family, by, contrast, estimate, se, df, statistic, p, p_adj, adjust).
