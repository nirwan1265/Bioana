# Classical ANOVA next to the model (afex + emmeans)

Rule table: `inst/rules/stat_tests.tsv`. Rules: `inst/rules/stat_rules.json`. Code: `R/stat_tests.R`.
Output per outcome data set: `bioana/<id>/anova_table.tsv`, `bioana/<id>/posthoc.tsv`,
`summary$analyses$<id>$stat_test`.

## Which test

The design picks it (not an assumption test):

| design (from analysis.model and the layout) | test | package |
|---|---|---|
| one random intercept on a `unit`/`subject` term, factors varying within it, every subject at every within cell | repeated-measures ANOVA (`within`) or mixed / split-plot ANOVA (`mixed`, with between-subject factors) | `afex::aov_car(y ~ between + Error(subject/(within)))` |
| no random terms, factors fully crossed | factorial ANOVA, type III (`between`); the between part as the model writes it (`treatment + block` stays additive) | `afex::aov_car(y ~ terms + Error(unit))` |
| subjects missing a within cell, other random terms (blocks, whole plots), a numeric (trend) term | not run, with the reason; the mixed model (lmer) covers these | - |

Repeated-measures ANOVA drops every subject with a missing cell; the mixed model keeps them,
so Bioana does not run the ANOVA then.

## What is reported

- F per term with numerator and denominator df, p, generalized eta squared (`ges`, Bakeman 2005)
  and partial eta squared (`pes`).
- Sphericity: Mauchly's test and the Greenhouse-Geisser and Huynh-Feldt epsilons for every
  within-subject term with more than two levels; the correction is `stat_rules.json`
  `sphericity_correction` (`GG` default; `girden` = HF when the GG epsilon > 0.75).
- Post-hoc (emmeans): families from the goal (`stat_rules.json` `families`); within-subject
  contrasts from the multivariate model (each contrast its own variance: a paired t-test without
  between factors). Holm adjustment within each family; `vs_control` uses Dunnett.
  Biodes can ask: `analysis.posthoc = {type: pairwise | vs_control | none, control, adjust}`.
- Checks (reported, never used to switch tests): extreme outliers per cell (3 IQR),
  Shapiro-Wilk on residuals, Levene across treatment groups (per within cell for mixed designs).

## Checked against

`tests/testthat/test-stat-tests.R`: repeated-measures and mixed ANOVA = `rstatix::anova_test`
(GG) and base `aov(y ~ day + Error(plant/day))`; Mauchly p = rstatix; post-hoc p = paired
t-tests (`rstatix::pairwise_t_test(paired = TRUE)`, Holm); factorial ANOVA =
`car::Anova(type = 3)` with sum-to-zero contrasts.

## References

Singmann H et al. afex (CRAN). Lenth RV. emmeans (CRAN). Greenhouse SW, Geisser S (1959)
Psychometrika 24:95. Huynh H, Feldt LS (1976) J Educ Stat 1:69. Girden ER (1992) ANOVA:
Repeated Measures. Bakeman R (2005) Behav Res Methods 37:379. Olejnik S, Algina J (2003)
Psychol Methods 8:434. Holm S (1979) Scand J Stat 6:65.
