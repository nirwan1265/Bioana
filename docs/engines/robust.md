# Robust tests (R/robust.R, WRS2; stat_tests.tsv family "robust")

Trimmed means (20% from each tail, `stat_rules.json` `trim`) and robust correlations (Wilcox 2017; Mair &
Wilcox 2020, WRS2). Run on request (`analysis.stat_tests: ["robust"]` or an id), and next to the primary test
when its checks find extreme outliers (`robust_when_outliers`). The primary test never changes.

| id | Layout | Call | Post-hoc |
|---|---|---|---|
| yuen | 2 independent groups | `WRS2::yuen` + `yuen.effect.ci` (bootstrap, seed 123) | none |
| yuen_paired | 2 levels within subjects | `WRS2::yuend` on the per-subject means | none |
| t1way | 3+ independent groups | `WRS2::t1way` | `WRS2::lincon` (Hochberg) |
| t2way | two between factors | `WRS2::t2way` | not run |
| rmanova_trimmed | within subjects | `WRS2::rmanova` | `WRS2::rmmcp` |
| bwtrim | between x within | `WRS2::bwtrim` | not run |

Correlation: `analysis.correlation.method` `percentage_bend` (`WRS2::pbcor`, beta 0.2) or `winsorized`
(`WRS2::wincor`, 20%). Regression: HC3 standard errors always (regression.R); MM regression with
`analysis.regression.robust: "mm"` (`robustbase::lmrob`).

Datanovia numbers reproduced: genderweight Yuen t 21.9556 (df 15.18), difference -22.45 [-24.62, -20.27],
effect 0.96; PlantGrowth t1way F 8.2818 (2, 9.64), p 0.00807, lincon trt1 - trt2 p 0.0071; jobsatisfaction
t2way gender 0.1788, education 339.66, interaction 8.7873; mtcars mpg-wt percentage-bend r -0.8962.
