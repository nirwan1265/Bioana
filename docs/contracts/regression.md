# Regression, covariates, robust tests and power (2026-11-26, Bioana 0.13.0)

Bioana reads what is below; Biodes writes it, Simitall simulates and records the truth. Methods and numbers:
docs/engines/regression.md and docs/engines/robust.md (datanovia regression, power and robust pages).

## Biodes

1. **Covariates (measured predictors).** Generic model terms `covariate`, `covariate2`, `covariate3`, ...
   (any term listed in `analysis.covariates` is a covariate, so `covariate5` and beyond work). Role covariate:
   always numeric, one value per unit, a column of `sample_metadata_template.tsv` (not the layout).

   ```json
   "analysis": {
     "model": "y ~ covariate + covariate2",
     "covariates": [
       {"term": "covariate",  "name": "glucose", "column": "glucose", "unit": "mg/dL"},
       {"term": "covariate2", "name": "BMI",     "column": "bmi",     "unit": "kg/m2", "transform": "log"}
     ]
   }
   ```

   `transform`: `log`, `log10`, `log2`, `sqrt` (applied before fitting; slopes are per unit of the
   transformed value). The template gets one blank column per covariate.
2. **Designs.**
   - `regression`: no randomized factor, `y ~ covariate (+ covariate2 ...)`, one row per unit.
   - ANCOVA: a randomized factor and a baseline covariate, `y ~ treatment + covariate` (the catalog's
     `ancova_baseline` can now carry its covariate).
   - Interactions: `covariate * covariate2`, `treatment * covariate`.
   - Categorical predictors stay treatment-type terms (`treatment`, `treatment2`, or the column name) with
     `factors[].levels` giving the reference first.
3. **Goal `regression`.** `analysis.goals: ["regression"]`, `goal_terms.regression: ["covariate", ...]`
   (default: every covariate), `goal_effect.regression: {size, scale}` with scale `per_unit` (slope, continuous
   outcome) or `odds_ratio_per_unit` (yes/no outcome). Yes/no outcomes with covariates are logistic regression.
4. **Options** `analysis.regression`:
   - `predict_at: [{"glucose": 150}, ...]`: points to predict at (keys: term, covariate name or column),
     confidence and prediction intervals (probabilities for logistic).
   - `robust: "mm"`: MM regression (robustbase) next to OLS. HC3 robust SEs are always reported.
   - `order_column`: a column giving measurement order, for the Durbin-Watson test.
5. **Robust tests** (trimmed means, WRS2) on request: `analysis.stat_tests: ["robust"]` or an id (`yuen`,
   `yuen_paired`, `t1way`, `t2way`, `rmanova_trimmed`, `bwtrim`). Robust correlation:
   `analysis.correlation.method` `percentage_bend` or `winsorized`.
6. **Power belongs to Biodes** (a-priori, before data). With `pwr` (datanovia power page):
   - two groups `pwr.t.test(d, power)`;
   - k groups `pwr.anova.test(k, f)`;
   - correlation `pwr.r.test(r)`;
   - two proportions `pwr.2p.test(ES.h(p1, p2))`;
   - goodness of fit / table `pwr.chisq.test(ES.w1 / ES.w2, df)`;
   - regression `pwr.f2.test(u = predictors tested, f2 = R2 / (1 - R2))`, n = v + u + 1.

   Write `power: {method, effect_size_name, effect_size, alpha, power, n_per_group | n_total, source}`
   (`source`: pilot, literature, or Cohen's convention: d 0.2/0.5/0.8, f 0.1/0.25/0.4, r 0.1/0.3/0.5,
   w 0.1/0.3/0.5, f2 0.02/0.15/0.35). Do not power on the study's own data. Note per-group vs total n.

## Simitall

1. Simulate covariates per unit: normal by default (`mean`, `sd` from `analysis.covariates` if given), with a
   correlation between covariates (`covariate_cor`) so VIF can be checked. Optional planted problems for the
   diagnostics: heteroscedastic noise, a curved effect, a few outliers or high-leverage units.
2. Outcome = intercept + treatment effects + slope x covariate (+ interactions) + noise; logistic for yes/no
   (log odds = intercept + log OR per unit x covariate).
3. Truth in `goal_truth`: `true_slopes: {term: slope per unit (after transform)}` (continuous) or
   `true_log_odds_ratios: {term: log OR per unit}` (binary), `model_terms`, `status`. Optional
   `regression_truth: {r_squared, residual_sd, covariate_cor, planted_outliers: [rows]}`.
4. Power by simulation (checks Biodes' a-priori n): repeat the folder over seeds, run Bioana, share of seeds
   with the goal detected.

## Bioana

Reads the above (`read_design.R`: covariate terms are looked up in the sample metadata; this file is a copy of
Simitall's `suite_handoff.R`, so the change goes to both). Writes, per outcome:
`regression_coefficients.tsv`, `regression_type2.tsv`, `regression_influence.tsv`,
`regression_influence_refit.tsv`, `regression_vif.tsv`, `regression_simple_slopes.tsv`,
`regression_adjusted_means.tsv`, `regression_predictions.tsv`, `regression_sensitivity.tsv`,
`regression_robust.tsv`; `summary.regression` (fit, checks, coefficients, interactions, ANCOVA, AUC,
Hosmer-Lemeshow, sensitivity). Goal truth: `goal_recovery.slopes` (each slope vs `true_slopes` with its CI).
Bioana never reports observed power. It reports the smallest effect the sample could detect (sensitivity).
