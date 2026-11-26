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

## Files Bioana writes (for Biopaint)

`cells.tsv` (chi-squared residuals), `descriptives.tsv` (counts, proportions, exact CIs),
`correlation.tsv`, `correlation_matrix*.tsv`, `agreement.tsv`, `agreement_table.tsv`.
