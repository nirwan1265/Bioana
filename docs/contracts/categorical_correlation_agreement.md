# Categorical outcomes, correlation between traits, agreement between raters (2026-11-26)

Bioana 0.11.0. Nothing here blocks: Bioana picks the test from the design and the outcome type.

## Biodes

1. `primary_outcome.type = "nominal"` for unordered classes (flower colour, seed coat class,
   disease type), with `primary_outcome.categories` (exact labels; the first is the reference).
   Yes/no stays `binary` (0/1); scores stay `ordinal`.
2. `analysis.expected_ratio = {label: weight}` when the question states a ratio: segregation
   ("3:1 purple to white" -> {"purple": 3, "white": 1}; dihybrid 9:3:3:1; 1:2:1). Labels must be
   the category labels in the data (`model: "y ~ 1"` when there is no treatment).
3. Paired yes/no designs (same plants before/after) need nothing new: `(1 | unit)` +
   `repeated_factors` gives McNemar / Cochran's Q.
4. Yes/no over a numeric dose: the `trend` goal with the factor's values gives Cochran-Armitage.
5. `analysis.correlation = {variables: [columns], method: "pearson" | "spearman" | "kendall"}` when
   the question asks how traits relate ("is height related to yield"). Please list the traits as
   `secondary_outcomes` [{name, unit, column}] and add their columns to
   `sample_metadata_template.tsv`. `within_design` defaults to true.
6. `analysis.agreement = {raters: [columns], scale: "nominal" | "ordinal" | "continuous", levels}`
   when several raters or methods score the same units ("two people scored rust 1-9 on each
   plot"; "manual vs imaging leaf area"). Long format instead: `rater_column`, `score_column`,
   `unit_column`. Ordinal levels lowest to highest.

## Simitall

1. Nominal outcomes: categories from a multinomial logit with the treatment effects; truth: the
   category probabilities per group. Segregation: classes drawn from `expected_ratio` (status quo)
   or a distorted ratio as the planted effect (truth: the true ratio).
2. Paired yes/no: a unit random effect on the logit scale (truth: the conditional odds ratio).
3. Secondary traits for correlation: traits sharing treatment effects and a within-design
   residual correlation matrix (truth: the within-design correlations; the raw ones follow).
4. Ratings: continuous = true value + rater bias + rater noise (truth: the true ICC2 from the
   variances); categorical = a misclassification matrix per rater (truth: the expected kappa).

## Read since (Bioana 0.12.0; Biodes 0.6.7 89e9b49, Simitall 480f38d)

- `single_group` (`y ~ 1`, `analysis.goals: []`): no effect goal. The test comes from `expected_ratio`
  (goodness of fit), `correlation` or `agreement`. With long-format agreement (`rater_column`,
  `score_column` = the outcome column) the outcome entry is skipped: the column holds the scores.
- Binary with two `primary_outcome.categories` (or `expected_ratio` labels): text values coded
  0 = the first category (reference), 1 = the second; 0/1 values are read that way for the goodness of fit.
- Biodes 0e5dcc9: an outcome tested against `expected_ratio` is nominal (two classes too), categories in ratio order.
  `secondary_outcomes[].type` sets each pair's method when the method is unstated (Pearson for continuous,
  Spearman when a count or ordinal trait is involved; `correlation_rules.json method_by_type`).
- `analysis.correlation.variables`: trait names. `primary_outcome.name` -> column `primary_outcome`,
  `secondary_outcomes[{name, column}]` -> their columns (names are matched case- and punctuation-free).
- Survival with one row per subsample: the unit is `analysis.survival.cluster` (or `plot` when there is a
  `subsample` column); clusters get a shared frailty and the robust score test.
- `stat_test$p` (and `statistic_name`, `statistic`, `df1`, `headline_term`): the headline test, the row of the
  tested factor. Before 0.12.0 only `stat_test$table[[i]]$p` existed, and `$p` partially matched `$package`.
- Small fixed blocks in binary / ordinal / nominal models (< 10 units per block, `stat_rules.json`): binary
  and ordinal get a random block intercept (ORs conditional on the block); nominal drops the block from the
  multinomial model, and the generalized CMH test stratifies on it. Fixed coefficients for pairs double the log OR.

Truth scored (`goal_recovery` / `recovery`):

| Simitall truth | Bioana estimate | Output |
|---|---|---|
| `outcome_model.group_probability` | observed share per treatment cell, Clopper-Pearson CI | `group_probability_recovery.tsv`, `goal_recovery.probabilities` |
| `outcome_model.conditional_odds_ratios` | model OR (glmer with random block = conditional), Wald CI | `goal_recovery.odds_ratios` |
| `outcome_model.true_ratio`, `distorted` | goodness-of-fit p (detected = distorted?), class shares with CIs | `ratio_recovery.tsv`, `goal_recovery` (goal `segregation_ratio`) |
| `secondary_traits.within_design_correlation`; `expected_observed_correlation` for pairs with a count trait (Simitall 22031bb) | within-design (residual) correlation; Spearman/Kendall pairs vs the truth converted from Pearson | `correlation_recovery.tsv` (`truth_source`, `true_for_method`), `recovery` |
| `agreement.icc2` / `icc3` | psych ICC2 / ICC3 with CI | `agreement_recovery.tsv` |
| `agreement.expected_kappa`, `expected_weighted_kappa_quadratic` | Cohen's (2 raters, CI) or Fleiss' and Light's kappa; quadratic weighted kappa (mean of pairs for 3+) | `agreement_recovery.tsv` |

## Files Bioana writes (for Biopaint)

`cells.tsv` (chi-squared residuals), `descriptives.tsv` (counts, proportions, exact CIs),
`correlation.tsv`, `correlation_matrix*.tsv`, `agreement.tsv`, `agreement_table.tsv`, and the recovery
tables above.
