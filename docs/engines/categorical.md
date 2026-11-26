# Categorical outcomes and tests

Rows in `inst/rules/stat_tests.tsv`, adapters in `R/stat_tests.R`, rules in `stat_rules.json`.
Checked against stats, rstatix, epitools, nnet and the categorical tutorials (datanovia):
`tests/testthat/test-categorical.R`.

| design | yes/no outcome (`binary`) | unordered categories (`nominal`) |
|---|---|---|
| 2 independent groups | proportion test (Fisher when an expected count < 5), difference with CI, Cohen's h, odds ratio and relative risk (epitools, Wald) | chi-squared + Cramer's V, cell residuals (`cells.tsv`); Fisher when sparse |
| 3+ independent groups | proportion test + pairwise | as above |
| groups within blocks (RCBD, strata) | Cochran-Mantel-Haenszel, common odds ratio | - |
| numeric dose (`trend` goal) | Cochran-Armitage trend test (scores = dose values) | - |
| 2 levels on the same subjects | McNemar (exact binomial on the discordant pairs when <= 20) | - |
| 3+ levels on the same subjects | Cochran's Q + pairwise McNemar | - |
| expected ratio (`analysis.expected_ratio`) | goodness of fit per group and overall, Cohen's w; exact multinomial when sparse | same |

The model next to the tests: logistic (binary, glm / glmer) and multinomial logistic for nominal
outcomes (`nnet::multinom`, one coefficient per category vs the reference category).
Descriptives: counts and proportions with exact (Clopper-Pearson) CIs per group and category.

Plant uses: F2 segregation (3:1, 9:3:3:1, 1:2:1) with `expected_ratio`; flower or seed colour
classes by genotype; germinated yes/no by treatment, per block (CMH) or over doses (trend);
disease present before and after a treatment on the same plants (McNemar).

References: Pearson 1900; Fisher 1935; McNemar 1947 Psychometrika 12:153; Cochran 1950 Biometrika
37:256; Cochran 1954 Biometrics 10:417; Armitage 1955 Biometrics 11:375; Mantel & Haenszel 1959 J
Natl Cancer Inst 22:719; Cramer 1946; Agresti 2002 Categorical Data Analysis; Aragon epitools (CRAN).
