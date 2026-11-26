# Regression (R/regression.R, rules inst/rules/regression_rules.json)

Runs when `analysis.model` has covariate terms and the outcome is continuous (lm) or yes/no (logistic glm).
It replaces the group tests (stat_tests.tsv) for that outcome. Checks are reported and never switch the model.
Tutorials: datanovia biostatistics/regression (simple, multiple, diagnostics, categorical, interactions,
VIF, predict, logistic, logistic diagnostics) and the power page. Tests reproduce their numbers
(tests/testthat/test-regression.R).

| Output | Linear | Logistic |
|---|---|---|
| Coefficients | b, SE, t, p, CI (`confint`), standardized beta, HC3 SE / p / CI (`sandwich`) | log odds, OR with profile CI (`confint`), Wald CI |
| Fit | R2, adj R2, F(df1, df2), residual SE, AIC, BIC, f2 | LR chi-squared vs null, McFadden / Cox-Snell / Nagelkerke R2, AIC |
| Terms | Type II F (`car::Anova`) | Type II likelihood-ratio (`car::Anova`) |
| Outliers | standardized residual > 3 | standardized residual > 3 |
| Leverage / influence | hat > 2(p+1)/n; Cook's D > 4/(n-p-1), > 0.5 strong; refit without the largest D (shift in SE and %) | same |
| Assumptions | Shapiro-Wilk on residuals, Breusch-Pagan, RESET, Durbin-Watson (with `order_column`) | Box-Tidwell (x log x, LR) per positive predictor, Hosmer-Lemeshow (10 groups), separation, events per variable (>= 10) |
| Collinearity | VIF (`car::vif`; GVIF^(1/(2df))^2 for factors; with interactions per predictor) > 5 | same |
| Interactions | nested F of each covariate interaction; simple slopes (`emmeans::emtrends`) at mean -1 SD / mean / +1 SD or per level | nested LR; slopes on the log-odds scale |
| ANCOVA | equal-slopes test (treatment x covariate added); adjusted means at the covariate mean (`emmeans`) | same, probabilities |
| Prediction | `predict_at`: fit, confidence and prediction intervals, extrapolation flag | probability with CI from the link scale |
| Other | sensitivity: smallest partial f2 (and slope) detectable at 80% power (`pwr::pwr.f2.test`) | AUC with DeLong CI (`pROC`), in-sample classification at 0.5 |

Datanovia numbers reproduced: marketing `sales ~ youtube` b 0.047537 [0.04223, 0.05284], R2 0.6119, F 312.1,
max hat 0.02021, 13 Cook's D above 0.0202, slope 0.04834 without the most influential row, predictions 13.19
and 17.95; full model F(3, 196) 570.3, adj R2 0.8956; interaction F 429.6; Salaries Type II F 0.6324,
100.2572, 33.8582, 1.5137; Boston VIF tax 9.008554, rad 7.484496; cars prediction intervals; Pima glucose OR
1.043 [1.034, 1.054], probabilities 0.093 / 0.567 / 0.916, pedigree OR 3.13 [1.38, 7.37], age VIF 1.974.

What fails and what to do (from the tutorials): curvature -> `transform` (log, sqrt) or a polynomial term;
unequal variance -> HC3 SEs (in the table) or log the outcome; influential unit -> check for an entry error,
keep and report the refit; VIF > 5 -> drop or combine one predictor, or ridge / PCA; logit curvature ->
log, polynomial or spline; few events per variable or separation -> fewer predictors or Firth's penalized
likelihood.

Power: a-priori power is Biodes' (contracts/regression.md). Observed ("post-hoc") power is circular and is
never reported; the sensitivity line says what the sample could have detected.
