# Regression layer (R/regression.R) against base R / car / lmtest and the datanovia regression tutorials;
# robust tests (R/robust.R) against WRS2 and the datanovia robust page.

reg_dir <- function(data, y, covs, model, type = "continuous", factors = list(), terms = list(), extra = list(), goals = list("regression"), truth = NULL) {
  d <- file.path(tempfile("reg_"), "design"); dir.create(d, recursive = TRUE)
  n <- nrow(data); data$plot <- seq_len(n); data$sample_id <- paste0("S", seq_len(n)); data$primary_outcome <- data[[y]]
  cv <- lapply(seq_along(covs), function(i) list(term = if (i == 1) "covariate" else paste0("covariate", i), name = covs[[i]], column = covs[[i]]))
  spec <- list(handoff_version = "0.6.7", origin = "planned", organism = "test", system = "lab", factors = factors, design = list(chosen = "regression"),
               primary_outcome = list(name = y, type = type),
               analysis = c(list(model = model, terms = terms, covariates = cv, goals = goals, goal_terms = list(regression = list("covariate"))), extra))
  jsonlite::write_json(spec, file.path(d, "study_spec.json"), auto_unbox = TRUE, pretty = TRUE, null = "null", digits = NA)
  utils::write.csv(data[c("plot", intersect(unlist(terms$layout_columns), names(data)))], file.path(d, "field_layout.csv"), row.names = FALSE)
  utils::write.table(data, file.path(d, "sample_metadata_template.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
  if (!is.null(truth)) jsonlite::write_json(truth, file.path(d, "sim_truth.json"), auto_unbox = TRUE, digits = NA)
  bioana_write_manifest(d, list(id = "outcome", data_type = "outcome", files = list(table = "sample_metadata_template.tsv"), column = "primary_outcome",
                                truth = if (!is.null(truth)) list(kind = "outcome", file = "sim_truth.json", goal = list(goal = "regression", field = "goal_truth", terms = list("covariate"), test = "wald"))))
  d
}
reg_run <- function(d) { utils::capture.output(r <- bioana_analyze(d, secondary = FALSE)); list(r = r, rg = r$summary$analyses$outcome$regression, dir = file.path(d, "bioana", "outcome")) }
rtsv <- function(x, f) utils::read.delim(file.path(x$dir, paste0("regression_", f, ".tsv")), stringsAsFactors = FALSE)

test_that("covariates: terms from analysis.covariates read from the sample metadata, numeric, named; the regression goal tests their slopes", {
  set.seed(1); n <- 60; x <- stats::runif(n, 1, 10); y <- 2 + 0.5 * x + stats::rnorm(n)
  d <- reg_dir(data.frame(height = x, yield = y), "yield", list("height"), "y ~ covariate",
               truth = list(goal_truth = list(goal = "regression", model_terms = list("covariate"), status = "planted", true_slopes = list(covariate = 0.5))))
  x_ <- reg_run(d); a <- x_$r$summary$analyses$outcome
  ref <- stats::lm(y ~ x)
  expect_equal(x_$r$results$outcome$estimate[2], unname(stats::coef(ref)[2]))
  expect_equal(a$goal$id, "regression"); expect_equal(unlist(a$goal$tested_terms), "covariate")
  expect_equal(x_$rg$coefficients[[2]]$predictor, "height")
  expect_equal(a$goal_recovery$true_effect, 0.5); expect_true(a$goal_recovery$covered)
  ci <- stats::confint(ref)[2, ]; expect_equal(c(a$goal_recovery$slopes$ci_low, a$goal_recovery$slopes$ci_high), unname(ci))
})

test_that("simple and multiple linear regression = datanovia marketing (coefficients, CI, R2, F, diagnostics, influence refit, predictions)", {
  skip_if_not_installed("datarium"); skip_if_not_installed("car"); skip_if_not_installed("lmtest")
  data("marketing", package = "datarium", envir = environment())
  x <- reg_run(reg_dir(marketing, "sales", list("youtube"), "y ~ covariate", extra = list(regression = list(predict_at = list(list(youtube = 100), list(youtube = 200))))))
  co <- rtsv(x, "coefficients")
  expect_equal(round(co$estimate[2], 6), 0.047537); expect_equal(round(c(co$ci_low[2], co$ci_high[2]), 8), c(0.04223072, 0.05284256))
  expect_equal(round(x$rg$fit$r_squared, 4), 0.6119); expect_equal(round(x$rg$fit$F, 1), 312.1)
  inf <- rtsv(x, "influence")
  expect_equal(round(max(inf$leverage), 5), 0.02021); expect_equal(sum(inf$influential), 13L); expect_equal(round(max(inf$cooks_distance), 4), 0.0605)
  expect_equal(round(rtsv(x, "influence_refit")$without_row, 5), 0.04834)
  pr <- rtsv(x, "predictions"); expect_equal(round(pr$fit, 5), c(13.19278, 17.94644)); expect_equal(round(pr$ci_low, 5), c(12.51317, 17.38703))
  x <- reg_run(reg_dir(marketing, "sales", list("youtube", "facebook", "newspaper"), "y ~ covariate + covariate2 + covariate3"))
  expect_equal(round(x$rg$fit$adj_r_squared, 4), 0.8956); expect_equal(round(x$rg$fit$F, 1), 570.3)
  co <- rtsv(x, "coefficients"); expect_equal(round(co$p[4], 2), 0.86)
  ref <- stats::lm(sales ~ youtube + facebook + newspaper, data = marketing)
  expect_equal(co$std_beta[2], unname(stats::coef(ref)[2] * stats::sd(marketing$youtube) / stats::sd(marketing$sales)))
  expect_equal(co$se_hc3, unname(sqrt(diag(sandwich::vcovHC(ref, type = "HC3")))))
  bp <- Filter(function(z) z$check == "homoscedasticity", x$rg$checks)[[1]]; expect_equal(bp$value, unname(lmtest::bptest(ref)$p.value))
})

test_that("interaction (datanovia): nested F 429.6 and simple slopes = emtrends at mean -1 SD / mean / +1 SD", {
  skip_if_not_installed("datarium"); skip_if_not_installed("emmeans")
  data("marketing", package = "datarium", envir = environment())
  x <- reg_run(reg_dir(marketing, "sales", list("youtube", "facebook"), "y ~ covariate * covariate2"))
  expect_equal(round(x$rg$interactions[[1]]$statistic, 2), 429.59)
  ss <- rtsv(x, "simple_slopes"); ref <- stats::lm(sales ~ youtube * facebook, data = marketing); b <- stats::coef(ref)
  m <- mean(marketing$facebook); s <- stats::sd(marketing$facebook)
  expect_equal(ss$slope, unname(b["youtube"] + b["youtube:facebook"] * c(m - s, m, m + s)))
})

test_that("categorical predictors and ANCOVA (datanovia Salaries): Type II F, coefficients, equal slopes, adjusted means", {
  skip_if_not_installed("carData"); skip_if_not_installed("car")
  data("Salaries", package = "carData", envir = environment()); S <- Salaries; for (v in c("rank", "discipline", "sex")) S[[v]] <- as.character(S[[v]])
  fac <- list(list(name = "rank", n_levels = 3L, levels = list("AsstProf", "AssocProf", "Prof")), list(name = "discipline", n_levels = 2L, levels = list("A", "B")),
              list(name = "sex", n_levels = 2L, levels = list("Female", "Male")))
  x <- reg_run(reg_dir(S, "salary", list("yrs.service"), "y ~ covariate + rank + discipline + sex", factors = fac,
                       terms = list(rank = "rank", discipline = "discipline", sex = "sex", layout_columns = list("rank", "discipline", "sex"))))
  a2 <- rtsv(x, "type2"); expect_equal(round(a2$statistic, 4), c(0.6324, 100.2572, 33.8582, 1.5137))
  co <- rtsv(x, "coefficients"); expect_equal(round(co$estimate[co$predictor == "sexMale"], 2), 4771.25)
  ref <- stats::lm(salary ~ yrs.service + rank + discipline + sex, data = Salaries)
  full <- stats::update(ref, . ~ . + rank:yrs.service + discipline:yrs.service + sex:yrs.service)
  expect_equal(x$rg$ancova$equal_slopes_p, stats::anova(ref, full)[2, "Pr(>F)"])
  am <- rtsv(x, "adjusted_means"); em <- as.data.frame(emmeans::emmeans(ref, ~ rank)); expect_equal(am$adjusted_mean, em$emmean)
})

test_that("VIF (datanovia Boston: tax 9.01, rad 7.48) and prediction vs confidence intervals (datanovia cars)", {
  skip_if_not_installed("car")
  B <- MASS::Boston; nm <- setdiff(names(B), "medv")
  x <- reg_run(reg_dir(B, "medv", as.list(nm), paste("y ~", paste(c("covariate", paste0("covariate", 2:13)), collapse = " + "))))
  v <- rtsv(x, "vif"); expect_equal(round(v$vif[v$predictor == "tax"], 6), 9.008554); expect_equal(round(v$vif[v$predictor == "rad"], 6), 7.484496)
  expect_true(any(vapply(x$rg$checks, function(z) z$check == "collinearity" && z$status == "warn", logical(1))))
  x <- reg_run(reg_dir(datasets::cars, "dist", list("speed"), "y ~ covariate", extra = list(regression = list(predict_at = list(list(speed = 12), list(speed = 19), list(speed = 24), list(speed = 30))))))
  pr <- rtsv(x, "predictions")
  expect_equal(round(pr$pi_low[1:3], 6), c(-1.749529, 25.761756, 44.752478)); expect_equal(round(pr$ci_high[2], 5), 62.44421)
  expect_equal(pr$extrapolation, c(FALSE, FALSE, FALSE, TRUE))
})

test_that("logistic regression (datanovia Pima): profile-CI odds ratios, predicted probabilities, LRT, pseudo R2, VIF, AUC, Hosmer-Lemeshow", {
  skip_if_not_installed("mlbench"); skip_if_not_installed("car")
  data("PimaIndiansDiabetes2", package = "mlbench", envir = environment()); pd <- stats::na.omit(PimaIndiansDiabetes2); pd$dia <- as.integer(pd$diabetes == "pos")
  x <- reg_run(reg_dir(pd, "dia", list("glucose"), "y ~ covariate", type = "binary", extra = list(regression = list(predict_at = list(list(glucose = 90), list(glucose = 150), list(glucose = 200))))))
  co <- rtsv(x, "coefficients")
  expect_equal(round(co$odds_ratio[2], 6), 1.043334); expect_equal(round(c(co$or_ci_low[2], co$or_ci_high[2]), 4), c(1.0340, 1.0536))
  expect_equal(round(rtsv(x, "predictions")$probability, 6), c(0.092992, 0.566510, 0.915960))
  ref <- stats::glm(dia ~ glucose, data = pd, family = stats::binomial())
  expect_equal(x$rg$fit$lr_chisq, ref$null.deviance - ref$deviance)
  expect_equal(x$rg$fit$mcfadden_r2, 1 - as.numeric(stats::logLik(ref)) / as.numeric(stats::logLik(stats::glm(dia ~ 1, data = pd, family = stats::binomial()))))
  r <- rank(stats::fitted(ref)); n1 <- sum(pd$dia); n0 <- nrow(pd) - n1
  expect_equal(x$rg$auc$auc, (sum(r[pd$dia == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0))
  nm <- c("pregnant", "glucose", "pressure", "triceps", "insulin", "mass", "pedigree", "age")
  x <- reg_run(reg_dir(pd, "dia", as.list(nm), paste("y ~", paste(c("covariate", paste0("covariate", 2:8)), collapse = " + ")), type = "binary"))
  co <- rtsv(x, "coefficients"); expect_equal(round(co$odds_ratio[co$predictor == "pedigree"], 2), 3.13)
  expect_equal(round(c(co$or_ci_low[co$predictor == "pedigree"], co$or_ci_high[co$predictor == "pedigree"]), 2), c(1.38, 7.37))
  v <- rtsv(x, "vif"); expect_equal(round(v$vif[v$predictor == "age"], 3), 1.974)
  expect_true(any(vapply(x$rg$checks, function(z) z$check == "linearity_logit_covariate8" && z$status == "warn", logical(1))))   # Box-Tidwell: age bends
})

test_that("robust tests on request = WRS2 and the datanovia robust page (Yuen, t1way + lincon, t2way); outliers add the robust test", {
  skip_if_not_installed("WRS2")
  one <- function(df, y, g, levels, stat_tests = list("robust"), g2 = NULL, levels2 = NULL) {
    d <- file.path(tempfile("rob_"), "design"); dir.create(d, recursive = TRUE)
    df$plot <- seq_len(nrow(df)); df$sample_id <- paste0("S", df$plot); df$primary_outcome <- df[[y]]; df$treatment <- as.character(df[[g]])
    fac <- list(list(name = g, n_levels = length(levels), levels = as.list(levels))); model <- "y ~ treatment"; lc <- list("treatment")
    if (!is.null(g2)) { df$treatment2 <- as.character(df[[g2]]); fac <- c(fac, list(list(name = g2, n_levels = length(levels2), levels = as.list(levels2)))); model <- "y ~ treatment * treatment2"; lc <- list("treatment", "treatment2") }
    spec <- list(handoff_version = "0.6.7", origin = "planned", organism = "t", system = "lab", factors = fac, design = list(chosen = "crd"),
                 primary_outcome = list(name = y, type = "continuous"), analysis = list(model = model, terms = list(treatment = "treatment", treatment2 = if (!is.null(g2)) "treatment2", layout_columns = lc), stat_tests = stat_tests))
    jsonlite::write_json(spec, file.path(d, "study_spec.json"), auto_unbox = TRUE, null = "null")
    utils::write.csv(df[c("plot", unlist(lc))], file.path(d, "field_layout.csv"), row.names = FALSE)
    utils::write.table(df, file.path(d, "sample_metadata_template.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
    bioana_write_manifest(d, list(id = "outcome", data_type = "outcome", files = list(table = "sample_metadata_template.tsv"), column = "primary_outcome"))
    utils::capture.output(r <- bioana_analyze(d, secondary = FALSE)); st <- r$summary$analyses$outcome$stat_test
    c(list(st), st$also)
  }
  rob <- function(sts, id) Filter(function(s) identical(s$test, id), sts)[[1]]
  pg <- one(datasets::PlantGrowth, "weight", "group", c("ctrl", "trt1", "trt2"))
  t1 <- rob(pg, "t1way"); expect_equal(round(t1$table[[1]]$statistic, 4), 8.2818); expect_equal(round(t1$table[[1]]$df2, 2), 9.64)
  expect_equal(round(t1$table[[1]]$p, 5), 0.00807)
  skip_if_not_installed("datarium")
  data("genderweight", package = "datarium", envir = environment())
  gw <- one(as.data.frame(genderweight), "weight", "group", c("F", "M"))
  y <- rob(gw, "yuen")$table[[1]]; ref <- WRS2::yuen(weight ~ group, data = genderweight)
  expect_equal(y$statistic, ref$test); expect_equal(y$df1, ref$df); expect_equal(round(y$statistic, 4), 21.9556)
  data("jobsatisfaction", package = "datarium", envir = environment())
  js <- one(as.data.frame(jobsatisfaction), "score", "gender", c("male", "female"), g2 = "education_level", levels2 = c("school", "college", "university"))
  t2 <- rob(js, "t2way")$table; expect_equal(round(t2[[1]]$statistic, 4), 0.1788); expect_equal(round(t2[[3]]$statistic, 4), 8.7873)
  # extreme outlier, no request: the robust test is reported next to the Welch t-test
  set.seed(9); df <- data.frame(group = rep(c("a", "b"), each = 20), w = c(stats::rnorm(20, 10), stats::rnorm(19, 11), 60))
  sts <- one(df, "w", "group", c("a", "b"), stat_tests = list())
  expect_equal(sts[[1]]$test, "t_test_welch"); expect_true(any(vapply(sts, function(s) identical(s$test, "yuen"), logical(1))))
})

test_that("robust correlation (WRS2 percentage bend = datanovia mtcars r = -0.8962)", {
  skip_if_not_installed("WRS2")
  d <- file.path(tempfile("pb_"), "design"); dir.create(d, recursive = TRUE)
  m <- datasets::mtcars; m$plot <- seq_len(nrow(m)); m$sample_id <- rownames(m); m$primary_outcome <- m$mpg
  spec <- list(handoff_version = "0.6.7", origin = "planned", organism = "t", system = "lab", factors = list(), design = list(chosen = "single_group"),
               primary_outcome = list(name = "mpg", type = "continuous"), secondary_outcomes = list(list(name = "wt", column = "wt")),
               analysis = list(model = "y ~ 1", goals = list(), correlation = list(variables = list("mpg", "wt"), method = "percentage_bend")))
  jsonlite::write_json(spec, file.path(d, "study_spec.json"), auto_unbox = TRUE, null = "null")
  utils::write.csv(m["plot"], file.path(d, "field_layout.csv"), row.names = FALSE)
  utils::write.table(m, file.path(d, "sample_metadata_template.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
  utils::capture.output(r <- bioana_analyze(d, secondary = FALSE, only = "correlation"))
  expect_equal(round(r$results$correlation$r[1], 4), -0.8962); expect_equal(r$summary$analyses$correlation$package, "WRS2")
})
