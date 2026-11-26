# Tests on one treatment factor: t-tests, Wilcoxon, sign, Kruskal-Wallis, Friedman, Welch ANOVA, proportions

Rule table: `inst/rules/stat_tests.tsv`. Rules: `inst/rules/stat_rules.json`. Code: `R/stat_tests.R`.
Output: `bioana/<id>/stat_tests.tsv` (one row per test or ANOVA term), `posthoc.tsv`,
`descriptives.tsv`, `summary$analyses$<id>$stat_test` (`$also` for tests that were asked for).

## Which test runs

The first `auto = yes` row whose conditions match the outcome and the design, in the family of
the outcome type (`stat_rules.json` `default_family`: continuous parametric, ordinal
non-parametric, binary parametric). Conditions: outcome type, layout (between / within subjects),
number of treatment factors and their levels, blocking (none, a block of an RCBD, a subject).

| design | continuous (parametric) | continuous or ordinal (non-parametric) | binary |
|---|---|---|---|
| 2 independent groups | Welch t (`t_test_welch`) | Wilcoxon rank-sum | proportion test (Fisher when an expected count < 5) |
| 2 levels on the same subjects | paired t | Wilcoxon signed-rank | (McNemar: next step) |
| 3+ independent groups | one-way ANOVA | Kruskal-Wallis + Dunn | proportion test + pairwise |
| 3+ levels on the same subjects, or an RCBD | repeated-measures ANOVA / ANOVA with block | Friedman + pairwise signed-rank | - |

On request only (`auto = no`): Student t (`t_test_student`), sign test (`sign_test`), Welch ANOVA
with Games-Howell (`welch_anova`). Ask with `analysis.stat_tests` in study_spec.json or
`bioana_analyze(params = list(outcome = list(stat_tests = c("sign_test", "nonparametric"))))`:
a test id, or a family to add its first matching test. A request that does not fit the design is
reported, not run.

## Effect sizes

Cohen's d (Welch: root-mean-square SD; Student: pooled SD; paired: SD of the differences) and
Hedges' g (rstatix); r = Z / sqrt(N) for Wilcoxon (coin's standardized statistic, rstatix);
eta squared H = (H - k + 1) / (n - k) for Kruskal-Wallis; Kendall's W for Friedman; difference in
proportions with CI and Cohen's h (effectsize) for two proportions; odds ratio for Fisher 2x2.

## Checks (reported, never used to switch tests)

Extreme outliers (3 IQR) per group, Shapiro-Wilk per group, Levene (median centre) and
Fligner-Killeen; for paired designs the differences (outliers, Shapiro-Wilk, skewness: the
signed-rank test assumes symmetric differences, the sign test does not). Ties: normal
approximation, noted. Why not switch on a Shapiro-Wilk result: pre-testing changes the error
rates of the test that follows, small samples rarely fail and large samples fail on trivial
departures. Ask for the non-parametric test (or a transformation in the model) instead.

## Descriptives

Per treatment cell: n, mean, sd, se, 95% CI half-width, median, q1, q3, iqr, min, max, skewness,
and `report` = "median (IQR)" when |skewness| > `stat_rules.json` `skew_report_median`, else
"mean (SD)". Binary: n, successes, proportion. Same numbers as
`rstatix::get_summary_stats(type = "common")`.

## Checked against (`tests/testthat/test-two-groups.R`)

Welch and Student t = `t.test`, rstatix; Cohen's d = `effectsize::cohens_d`; paired t =
`t.test(paired = TRUE)` and its t^2 = the repeated-measures ANOVA F; Wilcoxon = `wilcox.test`,
r from `coin::wilcox_test`; sign test = `rstatix::sign_test`; Kruskal-Wallis = `kruskal.test`;
Welch ANOVA = `oneway.test`; Friedman on subjects and on RCBD blocks = `friedman.test`, post-hoc =
`pairwise.wilcox.test(paired = TRUE)`; proportions = `prop.test`, Cohen's h by its formula;
Fisher = `fisher.test`; descriptives = `rstatix::get_summary_stats`.

## References

Welch 1947 Biometrika 34:28; Student 1908 Biometrika 6:1; Wilcoxon 1945 Biometrics Bull 1:80;
Mann & Whitney 1947 Ann Math Stat 18:50; Kruskal & Wallis 1952 JASA 47:583; Dunn 1964
Technometrics 6:241; Friedman 1937 JASA 32:675; Welch 1951 Biometrika 38:330; Games & Howell 1976
J Educ Stat 1:113; Dixon & Mood 1946 JASA 41:557; Cochran 1954 Biometrics 10:417; Newcombe 1998
Stat Med 17:873; Fisher 1935 J R Stat Soc 98:39; Cohen 1988 Statistical Power Analysis;
Kassambara A. rstatix (CRAN); Ben-Shachar et al. effectsize (CRAN); Hothorn et al. coin (CRAN).
