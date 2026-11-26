# Biodes 0.6.7 / Simitall 480f38d: one-group studies, binary labels, small blocks in categorical models, truth
# scoring for class probabilities, conditional odds ratios, segregation ratios, trait correlations and agreement,
# and one template row per survival subsample.

cat_dir <- function(meta, spec, truth = NULL, manifest = TRUE, column = "primary_outcome", event = NULL) {
  d <- file.path(tempfile("cat_"), "design"); dir.create(d, recursive = TRUE)
  spec <- utils::modifyList(list(handoff_version = "0.6.7", origin = "planned", organism = "test", system = "lab"), spec)
  jsonlite::write_json(spec, file.path(d, "study_spec.json"), auto_unbox = TRUE, pretty = TRUE, null = "null")
  lay <- meta[intersect(c("plot", "unit", "rep", "block", "treatment", "rater"), names(meta))]
  utils::write.csv(lay, file.path(d, "field_layout.csv"), row.names = FALSE)
  utils::write.table(meta, file.path(d, "sample_metadata_template.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
  if (!is.null(truth)) jsonlite::write_json(truth, file.path(d, "sim_truth.json"), auto_unbox = TRUE, digits = NA)
  if (manifest) bioana_write_manifest(d, list(id = "outcome", data_type = "outcome", files = list(table = "sample_metadata_template.tsv"), column = column,
                                              event_column = event, truth = if (!is.null(truth)) list(kind = "outcome", file = "sim_truth.json",
                                                                                                       goal = list(goal = "treatment_effect", field = "goal_truth", terms = list("treatment"), test = "wald"))))
  d
}
run <- function(d, ...) { utils::capture.output(r <- bioana_analyze(d, secondary = FALSE, ...)); r }

test_that("nominal in many small blocks: block dropped from multinom, generalized CMH = mantelhaen.test, stat_test$p, class probabilities scored", {
  set.seed(1); nb <- 60
  m <- data.frame(plot = 1:(2 * nb), block = rep(1:nb, each = 2), treatment = rep(c("A", "B"), nb), stringsAsFactors = FALSE)
  m$sample_id <- paste0("S", m$plot)
  pr <- list(A = c(white = 0.4, pink = 0.3, purple = 0.3), B = c(white = 0.2, pink = 0.3, purple = 0.5))
  m$primary_outcome <- vapply(m$treatment, function(g) sample(names(pr[[g]]), 1, prob = pr[[g]]), "")
  spec <- list(factors = list(list(name = "treatment", n_levels = 2L, levels = list("A", "B"))), design = list(chosen = "rcbd"),
               primary_outcome = list(name = "colour", type = "nominal", categories = list("white", "pink", "purple")),
               analysis = list(model = "y ~ treatment + block", terms = list(treatment = "treatment", block = "block", layout_columns = list("treatment", "block"))))
  tr <- list(outcome_model = list(family = "multinomial_logit", group_probability = lapply(pr, as.list)),
             goal_truth = list(goal = "treatment_effect", model_terms = list("treatment"), status = "planted", true_effect = 1))
  d <- cat_dir(m, spec, tr); r <- run(d); a <- r$summary$analyses$outcome
  expect_match(paste(a$warnings, collapse = " "), "dropped from the model")
  expect_false(any(grepl("^block", r$results$outcome$term)))
  st <- a$stat_test; expect_equal(st$test, "cmh")
  ref <- stats::mantelhaen.test(table(factor(m$treatment), factor(m$primary_outcome, c("white", "pink", "purple")), factor(m$block)))
  expect_equal(st[["p"]], ref$p.value); expect_equal(st$statistic, unname(ref$statistic))
  gp <- a$goal_recovery$probabilities; expect_equal(gp$cells, 6L)
  obs <- mean(m$primary_outcome[m$treatment == "B"] == "purple")
  tb <- utils::read.delim(list.files(d, "group_probability_recovery.tsv", recursive = TRUE, full.names = TRUE)[1])
  expect_equal(tb$estimate[tb$group == "B" & tb$category == "purple"], obs)
})

test_that("paired yes/no in blocks of 2: random block (glmer), conditional OR scored against conditional_odds_ratios", {
  skip_if_not_installed("lme4")
  set.seed(2); nb <- 300; u <- stats::rnorm(nb, 0, 0.6)
  m <- data.frame(plot = 1:(2 * nb), block = rep(1:nb, each = 2), treatment = rep(c("before", "after"), nb), stringsAsFactors = FALSE)
  m$sample_id <- paste0("S", m$plot)
  m$primary_outcome <- stats::rbinom(nrow(m), 1, stats::plogis(-1.2 + u[m$block] + log(3) * (m$treatment == "after")))
  spec <- list(factors = list(list(name = "treatment", n_levels = 2L, levels = list("before", "after"))), design = list(chosen = "rcbd"),
               primary_outcome = list(name = "wilted", type = "binary"),
               analysis = list(model = "y ~ treatment + block", terms = list(treatment = "treatment", block = "block", layout_columns = list("treatment", "block"))))
  tr <- list(outcome_model = list(family = "logit", conditional_odds_ratios = list(treatment = list(before = 1, after = 3)), random_effect_sd_logit = list(block = 0.6)),
             goal_truth = list(goal = "treatment_effect", model_terms = list("treatment"), status = "planted", true_effect = 1))
  r <- run(cat_dir(m, spec, tr)); a <- r$summary$analyses$outcome
  expect_true(a$engine$engine %in% c("glmer_binomial", "glmmtmb"))
  expect_match(paste(a$warnings, collapse = " "), "random intercept \\(1 \\| block\\)")
  orr <- a$goal_recovery$odds_ratios
  expect_equal(orr$estimate, exp(r$results$outcome$estimate[r$results$outcome$term == "treatmentafter"]))
  expect_equal(a$goal_recovery$true_effect, 3)
})

test_that("single_group 3:1 (Biodes y ~ 1, goals []): no effect goal; binary labels coded by the categories; segregation ratio scored", {
  set.seed(3); n <- 400
  m <- data.frame(plot = 1:n, sample_id = paste0("S", 1:n), stringsAsFactors = FALSE)
  m$primary_outcome <- sample(c("purple", "white"), n, TRUE, prob = c(0.65, 0.35))
  spec <- list(factors = list(), design = list(chosen = "single_group"), primary_outcome = list(name = "flower colour", type = "nominal", categories = list("purple", "white")),
               analysis = list(model = "y ~ 1", goals = list(), expected_ratio = list(purple = 3, white = 1)))
  tr <- list(outcome_model = list(family = "multinomial_logit", expected_ratio = list(purple = 3, white = 1), true_ratio = list(purple = 0.65, white = 0.35), distorted = TRUE))
  r <- run(cat_dir(m, spec, tr)); a <- r$summary$analyses$outcome
  g <- a$goal_recovery; expect_equal(g$goal, "segregation_ratio")
  ref <- stats::chisq.test(table(factor(m$primary_outcome, c("purple", "white"))), p = c(0.75, 0.25))
  expect_equal(g$p, ref$p.value); expect_true(g$detected); expect_true(g$correct)
  expect_equal(a$stat_test[["p"]], ref$p.value)
  # binary with two named categories: 0 = purple (first), 1 = white
  spec$primary_outcome$type <- "binary"
  r2 <- run(cat_dir(m, spec, NULL)); a2 <- r2$summary$analyses$outcome
  expect_match(paste(a2$warnings, collapse = " "), "coded 0 = purple")
  expect_equal(a2$stat_test[["p"]], ref$p.value)
  expect_null(a2$goal_recovery)
})

test_that("correlation: trait names map to template columns (primary_outcome, secondary_outcomes); within-design correlations scored", {
  set.seed(4); n <- 150
  z <- matrix(stats::rnorm(n * 3), n) %*% chol(matrix(c(1, 0.5, 0.5, 0.5, 1, 0.5, 0.5, 0.5, 1), 3))
  m <- data.frame(plot = 1:n, sample_id = paste0("S", 1:n), treatment = rep(c("A", "B", "C"), length.out = n), stringsAsFactors = FALSE)
  eff <- c(A = 0, B = 1, C = 2)[m$treatment]
  m$primary_outcome <- 10 + eff + z[, 1]; m$tiller_number <- 5 + eff + z[, 2]; m$grain_weight <- 3 + eff + z[, 3]
  spec <- list(factors = list(list(name = "treatment", n_levels = 3L, levels = list("A", "B", "C"))), design = list(chosen = "crd"),
               primary_outcome = list(name = "plant_height", type = "continuous"),
               secondary_outcomes = list(list(name = "tiller_number", column = "tiller_number"), list(name = "grain_weight", column = "grain_weight")),
               analysis = list(model = "y ~ treatment", terms = list(treatment = "treatment", layout_columns = list("treatment")),
                               correlation = list(variables = list("plant_height", "tiller_number", "grain_weight"), method = "pearson")))
  R <- list(primary_outcome = list(primary_outcome = 1, tiller_number = 0.5, grain_weight = 0.5), tiller_number = list(primary_outcome = 0.5, tiller_number = 1, grain_weight = 0.5),
            grain_weight = list(primary_outcome = 0.5, tiller_number = 0.5, grain_weight = 1))
  spec$assumptions <- list(analysis.correlation.method = "not stated: pearson")
  r <- run(cat_dir(m, spec, list(secondary_traits = list(within_design_correlation = R))), only = "correlation"); a <- r$summary$analyses$correlation
  expect_equal(a$method_assumed, "not stated: pearson")
  expect_equal(unlist(a$variables), c("primary_outcome", "tiller_number", "grain_weight"))
  res <- sapply(c("primary_outcome", "tiller_number"), function(v) stats::resid(stats::lm(m[[v]] ~ m$treatment)))
  expect_equal(r$results$correlation$r[r$results$correlation$scope == "within design"][1], stats::cor(res)[1, 2])
  expect_equal(a$recovery$pairs, 3L)
})

test_that("agreement: long format (rater column) skips the outcome; Simitall truth (ICC2/3, expected kappa) scored", {
  skip_if_not_installed("psych"); skip_if_not_installed("vcd")
  set.seed(5); n <- 120; tv <- stats::rnorm(n, 50, 10)
  m <- data.frame(plot = rep(1:n, each = 3), rater = rep(c("r1", "r2", "r3"), n), stringsAsFactors = FALSE); m$sample_id <- paste0("S", seq_len(nrow(m)))
  m$primary_outcome <- tv[m$plot] + c(r1 = 0, r2 = 2, r3 = -2)[m$rater] + stats::rnorm(nrow(m), 0, 4)
  spec <- list(factors = list(), design = list(chosen = "single_group"), primary_outcome = list(name = "severity", type = "continuous"),
               analysis = list(model = "y ~ 1", goals = list(), agreement = list(rater_column = "rater", score_column = "primary_outcome", unit_column = "plot", scale = "continuous")))
  vt <- 100; icc2 <- vt / (vt + 4 + 16); icc3 <- vt / (vt + 16)
  r <- run(cat_dir(m, spec, list(agreement = list(scale = "continuous", icc2 = icc2, icc3 = icc3))))
  expect_equal(r$summary$not_analysed$outcome$status, "skipped")
  rec <- r$summary$analyses$agreement$recovery
  ic <- psych::ICC(as.matrix(stats::reshape(m[c("plot", "rater", "primary_outcome")], idvar = "plot", timevar = "rater", direction = "wide")[-1]), lmer = FALSE)$results
  expect_equal(rec$measures[[1]]$estimate, ic$ICC[ic$type == "ICC2"]); expect_equal(rec$measures[[1]]$true, icc2)
  # ordinal, 3 raters: mean pairwise quadratic weighted kappa against expected_weighted_kappa_quadratic
  m$primary_outcome <- pmin(5, pmax(1, round((tv[m$plot] - 50) / 8 + 3 + stats::rnorm(nrow(m), 0, 0.5))))
  spec$primary_outcome$type <- "ordinal"; spec$analysis$agreement$scale <- "ordinal"; spec$analysis$agreement$levels <- as.list(1:5)
  r <- run(cat_dir(m, spec, list(agreement = list(scale = "ordinal", expected_kappa = 0.5, expected_weighted_kappa_quadratic = 0.8))))
  tab <- r$results$agreement; i <- grep("^mean pairwise weighted kappa \\(quadratic", tab$measure)
  W <- stats::reshape(m[c("plot", "rater", "primary_outcome")], idvar = "plot", timevar = "rater", direction = "wide")[-1]
  F_ <- lapply(W, factor, levels = 1:5)
  ref <- mean(utils::combn(3, 2, function(ij) vcd::Kappa(table(F_[[ij[1]]], F_[[ij[2]]]), weights = "Fleiss-Cohen")$Weighted[["value"]]))
  expect_equal(tab$estimate[i], ref)
  expect_true(any(vapply(r$summary$analyses$agreement$recovery$measures, function(x) grepl("quadratic", x$measure), logical(1))))
})

test_that("survival, one template row per seed (Biodes subsample column, analysis.survival.cluster = plot): plots are the clusters", {
  set.seed(6); np <- 20
  m <- data.frame(plot = rep(1:np, each = 6), subsample = rep(1:6, np), block = rep(rep(1:10, each = 2), each = 6), treatment = rep(rep(c("no", "yes"), 10), each = 6), stringsAsFactors = FALSE)
  m$sample_id <- paste0("S", seq_len(nrow(m)))
  t <- -log(stats::runif(nrow(m))) / (0.1 * exp(stats::rnorm(np, 0, 0.3)[m$plot]) * ifelse(m$treatment == "yes", 0.6, 1))
  m$primary_outcome <- round(pmin(t, 21), 2); m$event <- as.integer(t <= 21)
  spec <- list(factors = list(list(name = "treatment", n_levels = 2L, levels = list("no", "yes"))), design = list(chosen = "crd"),
               primary_outcome = list(name = "days", type = "survival", event_column = "event"),
               analysis = list(model = "y ~ treatment", terms = list(treatment = "treatment", layout_columns = list("treatment")), survival = list(cluster = "plot")))
  r <- run(cat_dir(m, spec, NULL, event = "event")); a <- r$summary$analyses$outcome
  expect_match(paste(a$warnings, collapse = " "), "20 units kept as clusters")
  expect_true(any(grepl("robust score", vapply(a$stat_test$table, `[[`, "", "statistic_name"))))
  rc <- survival::coxph(survival::Surv(primary_outcome, event) ~ treatment + cluster(plot), data = transform(m, treatment = factor(treatment, c("no", "yes"))))
  rs <- utils::tail(Filter(function(x) grepl("robust score", x$statistic_name), a$stat_test$table), 1)[[1]]
  expect_equal(rs$statistic, unname(summary(rc)$robscore["test"]), tolerance = 1e-6)
})

test_that("correlation method per trait type when unstated (Biodes 0e5dcc9): count or ordinal pairs Spearman, continuous pairs Pearson; a stated method applies to all", {
  set.seed(7); n <- 80
  m <- data.frame(plot = 1:n, sample_id = paste0("S", 1:n), stringsAsFactors = FALSE)
  m$primary_outcome <- stats::rnorm(n, 50, 5); m$leaf_area <- m$primary_outcome + stats::rnorm(n, 0, 5); m$tillers <- stats::rpois(n, exp(1 + 0.02 * m$primary_outcome))
  spec <- list(factors = list(), design = list(chosen = "single_group"), primary_outcome = list(name = "height", type = "continuous"),
               secondary_outcomes = list(list(name = "leaf_area", column = "leaf_area", type = "continuous"), list(name = "tillers", column = "tillers", type = "count")),
               assumptions = list(analysis.correlation.method = "not stated: pearson"),
               analysis = list(model = "y ~ 1", goals = list(), correlation = list(variables = list("height", "leaf_area", "tillers"), method = "pearson")))
  r <- run(cat_dir(m, spec, NULL), only = "correlation"); tb <- r$results$correlation
  expect_equal(tb$method, c("pearson", "spearman", "spearman"))
  expect_equal(tb$r[2], unname(suppressWarnings(stats::cor.test(m$primary_outcome, m$tillers, method = "spearman"))$estimate))
  spec$assumptions <- NULL
  r <- run(cat_dir(m, spec, NULL), only = "correlation"); expect_true(all(r$results$correlation$method == "pearson"))
})

test_that("count trait pairs scored against Simitall's expected_observed_correlation, on Spearman's scale (6/pi asin(r/2))", {
  set.seed(8); n <- 200
  m <- data.frame(plot = 1:n, sample_id = paste0("S", 1:n), stringsAsFactors = FALSE)
  z <- matrix(stats::rnorm(n * 2), n) %*% chol(matrix(c(1, 0.5, 0.5, 1), 2))
  m$primary_outcome <- 50 + 5 * z[, 1]; m$tillers <- stats::rpois(n, pmax(0.1, 6 + 2 * z[, 2]))
  spec <- list(factors = list(), design = list(chosen = "single_group"), primary_outcome = list(name = "height", type = "continuous"),
               secondary_outcomes = list(list(name = "tillers", column = "tillers", type = "count")),
               analysis = list(model = "y ~ 1", goals = list(), correlation = list(variables = list("height", "tillers"))))
  R <- list(primary_outcome = list(primary_outcome = 1, tillers = 0.5), tillers = list(primary_outcome = 0.5, tillers = 1))
  Ro <- list(primary_outcome = list(primary_outcome = 1, tillers = 0.4), tillers = list(primary_outcome = 0.4, tillers = 1))
  d <- cat_dir(m, spec, list(secondary_traits = list(within_design_correlation = R, expected_observed_correlation = Ro))); r <- run(d, only = "correlation")
  rec <- utils::read.delim(list.files(d, "correlation_recovery.tsv", recursive = TRUE, full.names = TRUE)[1])
  expect_equal(rec$truth_source, "expected_observed_correlation"); expect_equal(rec$method, "spearman")
  expect_equal(rec$true_for_method, 6 / pi * asin(0.4 / 2))
})
