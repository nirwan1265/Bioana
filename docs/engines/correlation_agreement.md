# Correlation between traits and agreement between raters

Code: `R/correlation.R` (data_type `traits`), `R/agreement.R` (data_type `ratings`). Rules:
`correlation_rules.json`, `agreement_rules.json`. Both run from the spec (`analysis.correlation`,
`analysis.agreement`) on the sample metadata, or from a manifest entry with its own table.
Checked: `tests/testthat/test-correlation-agreement.R` (mtcars numbers of the correlation tutorial,
Hmisc::rcorr, psych::cohen.kappa, vcd::Kappa, the Shrout-Fleiss ICC formula).

## Correlation (`analysis.correlation = {variables, method, adjust, within_design}`)

Every pair: `cor.test` (pearson default; spearman, kendall), CI (Pearson), Holm-adjusted p,
strength label. Also within the design: the residuals of each trait on the design's fixed terms
(and random block intercepts), because across plots two traits that both respond to the treatment
correlate even with no relation inside a treatment. Files: `correlation.tsv` (long, both scopes),
`correlation_matrix.tsv`, `correlation_matrix_within.tsv` (Biopaint draws the correlogram).
Pearson: Shapiro-Wilk per trait is reported; Spearman is the request for non-normal traits.

## Agreement (`analysis.agreement = {raters | rater_column + score_column + unit_column, scale, levels}`)

| scale | 2 raters | 3+ raters |
|---|---|---|
| nominal | Cohen's kappa with CI (vcd) | Fleiss' kappa (irr), else Light's kappa (mean pairwise Cohen) |
| ordinal | weighted kappa, linear (primary) and quadratic, plus unweighted | Light's / Fleiss |
| continuous | ICC, six Shrout-Fleiss forms (psych); primary ICC2: two-way, absolute agreement, single rater | same |

Interpretation: Landis & Koch for kappa; Koo & Li for the ICC, on the lower 95% bound.
Files: `agreement.tsv`, `agreement_table.tsv` (2 raters; Biopaint draws the agreement chart).
Plant uses: disease or lodging scores by two people; manual vs image-based measurements (ICC);
repeated scoring by one rater (test-retest).

Not here: Bland-Altman limits of agreement, partial correlation beyond the design adjustment.
