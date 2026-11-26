# Survival analysis next to the Cox engine

Code: `R/survival_tests.R`, `R/adtte.R`. Rules: `inst/rules/survival_rules.json`; rows in
`inst/rules/stat_tests.tsv`. Examples: `bioana_example("lung" | "colon" | "mgus2" | "rats")`.
Checked against survival, cmprsk, survMisc, survminer, KMsurv and the numbers of the survival
tutorials (datanovia): `tests/testthat/test-survival.R`.

In plant and lab work the same methods cover days to germination, flowering, disease onset,
death of seedlings or insects, shelf life: any time to an event with censoring (the event had not
happened by the end of the experiment, or the plant was lost).

## What runs for a survival outcome

| always | on request (`analysis.survival`) | when the data say so |
|---|---|---|
| Cox model (the engine; frailty for a random cluster) with hazard ratios, CIs, concordance, LR/Wald/score tests | weighted log-rank, one prespecified weight (`logrank_weight`: peto, gehan, tarone_ware, fh_1_0, fh_0_1, fh_1_1) | stratified log-rank when the design has blocks (randomization strata), next to the unstratified one |
| Kaplan-Meier per group: medians with CIs, quartiles, curve (`km_curve.tsv`) | landmarks (`landmarks`: times) -> `landmarks.tsv` | pairwise log-rank for 3+ groups (`analysis.posthoc`; Holm by default) |
| log-rank | AFT models (`aft`: true or distributions) -> time ratios of the lowest AIC; Weibull HR = TR^(-shape) | competing risks when the event column has codes 2+: CIF and Gray's test per cause, Fine-Gray, cause-specific Cox |
| proportional-hazards check (`cox.zph`) | life table (`life_table`: interval width) -> `life_table.tsv` | "median not reached" -> report landmarks or RMST |
| RMST difference and ratio (two groups; tau: `rmst_tau` or the smallest of the groups' largest follow-up times) | extra strata (`strata`) | events per Cox coefficient < 10 flagged |

A failed PH check is reported with what to do (the HR is then an average; RMST does not assume PH;
stratify a nuisance factor or ask for an AFT model); it does not switch the test.

## Files

`km_curve.tsv` (group, time, n_risk, n_event, n_censor, surv, lower, upper), `survival_summary.tsv`
(group, n, events, median and CI, q25, q75, rmst, rmst_se, tau, max_followup), `landmarks.tsv`,
`life_table.tsv`, `cif_curve.tsv` (competing risks), `stat_tests.tsv`, `posthoc.tsv`.

## CDISC ADTTE

`bioana_adtte_design(file, dir, paramcd = "OS", reference = "Placebo")` writes a design folder from
an ADTTE data set: event = 1 - CNSR (CDISC codes 1 = censored; passing CNSR to Surv() runs and is
wrong), one PARAMCD, TRT01P / TRT01A / ARM as the treatment, optional STRATA columns.

## Not here (yet)

Time-varying covariates (counting-process data, `tmerge`), time-varying coefficients (`tt()`),
interval censoring, stepwise/LASSO model selection (exploratory; designed experiments prespecify
the model), MaxCombo.

## References

Kaplan & Meier 1958 JASA 53:457; Mantel 1966 Cancer Chemother Rep 50:163; Peto & Peto 1972 JRSS A
135:185; Fleming & Harrington 1991; Cox 1972 JRSS B 34:187; Grambsch & Therneau 1994 Biometrika
81:515; Royston & Parmar 2013 BMC Med Res Methodol 13:152; Tian et al. 2020 Stat Med 39:4048
(tau); Gray 1988 Ann Stat 16:1141; Fine & Gray 1999 JASA 94:496; Wei 1992 Stat Med 11:1871;
Cutler & Ederer 1958 J Chronic Dis 8:699; Schoenfeld 1983 Biometrics 39:499 (events needed);
Therneau survival (CRAN); Gray cmprsk (CRAN); Kassambara survminer (CRAN).
