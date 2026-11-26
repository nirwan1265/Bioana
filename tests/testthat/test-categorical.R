# Categorical tests (R/stat_tests.R): equality with stats, rstatix, epitools, nnet and the tutorials' numbers.

cat_folder <- function(rows, model, terms, type = "binary", factors, an = list(), unit = NULL) {
  d <- file.path(tempfile("cat_"), "design"); dir.create(d, recursive = TRUE)
  rows$sample_id <- sprintf("S%04d", seq_len(nrow(rows)))
  spec <- list(handoff_version = "0.6.6", origin = "observed", organism = "test", system = "lab", factors = factors,
               design = list(chosen = "crd"), primary_outcome = list(name = "y", type = type),
               analysis = c(list(model = model, terms = terms), an))
  jsonlite::write_json(spec, file.path(d, "study_spec.json"), auto_unbox = TRUE, pretty = TRUE)
  utils::write.csv(rows[, setdiff(names(rows), "primary_outcome"), drop = FALSE], file.path(d, "field_layout.csv"), row.names = FALSE)
  utils::write.table(rows, file.path(d, "sample_metadata_template.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
  bioana_write_manifest(d, list(id = "outcome", data_type = "outcome", files = list(table = "sample_metadata_template.tsv"), column = "primary_outcome"))
  d
}
run_cat <- function(d) { utils::capture.output(r <- bioana_analyze(d, secondary = FALSE)); r }
fac <- function(name, lev) list(list(name = name, n_levels = length(lev), levels = as.list(lev)))
row1 <- function(st) st$table[[1]]

test_that("McNemar on paired yes/no = mcnemar.test and the tutorial (chi2 7.26, exact p 0.006)", {
  pairs <- rbind(matrix(c(0, 0), 25, 2, byrow = TRUE), matrix(c(0, 1), 6, 2, byrow = TRUE), matrix(c(1, 0), 21, 2, byrow = TRUE), matrix(c(1, 1), 10, 2, byrow = TRUE))
  rows <- data.frame(unit = rep(sprintf("P%02d", 1:62), 2), time = rep(c("before", "after"), each = 62), primary_outcome = c(pairs[, 1], pairs[, 2]), stringsAsFactors = FALSE)
  d <- cat_folder(rows, "y ~ treatment + (1 | unit)", list(treatment = "time", unit = "unit", layout_columns = list("time", "unit")), factors = fac("time", c("before", "after")))
  st <- run_cat(d)$summary$analyses$outcome$stat_test
  expect_equal(st$test, "mcnemar")
  tb <- table(factor(pairs[, 1], 0:1), factor(pairs[, 2], 0:1))
  expect_equal(row1(st)$statistic, unname(stats::mcnemar.test(tb)$statistic)); expect_equal(round(row1(st)$statistic, 2), 7.26)
  expect_equal(round(st$exact_p, 3), 0.006)
})

test_that("Cochran's Q = rstatix::cochran_qtest; pairwise McNemar (Holm)", {
  set.seed(3); n <- 40
  rows <- data.frame(unit = rep(sprintf("P%02d", 1:n), 3), task = rep(c("t1", "t2", "t3"), each = n), stringsAsFactors = FALSE)
  rows$primary_outcome <- as.integer(stats::runif(3 * n) < c(0.25, 0.45, 0.65)[match(rows$task, c("t1", "t2", "t3"))])
  d <- cat_folder(rows, "y ~ treatment + (1 | unit)", list(treatment = "task", unit = "unit", layout_columns = list("task", "unit")), factors = fac("task", c("t1", "t2", "t3")))
  st <- run_cat(d)$summary$analyses$outcome$stat_test
  expect_equal(st$test, "cochran_q")
  rows$task <- factor(rows$task); rows$unit <- factor(rows$unit)
  q <- rstatix::cochran_qtest(rows, primary_outcome ~ task | unit)
  expect_equal(row1(st)$statistic, unname(q$statistic)); expect_equal(st$posthoc$n_comparisons, 3)
})

test_that("CMH for yes/no in blocks = mantelhaen.test; common OR is treatment vs reference", {
  set.seed(4); g <- expand.grid(rep = 1:10, trt = c("ctrl", "trt"), block = paste0("B", 1:5), stringsAsFactors = FALSE)
  g$primary_outcome <- as.integer(stats::runif(nrow(g)) < stats::plogis(-0.5 + 1.2 * (g$trt == "trt") + 0.3 * as.integer(factor(g$block))))
  g$unit <- sprintf("U%03d", seq_len(nrow(g)))
  d <- cat_folder(g[, c("unit", "trt", "block", "primary_outcome")], "y ~ treatment + block", list(treatment = "trt", block = "block", layout_columns = list("trt", "block")),
                  factors = fac("trt", c("ctrl", "trt")))
  st <- run_cat(d)$summary$analyses$outcome$stat_test
  expect_equal(st$test, "cmh")
  tb <- table(factor(g$trt, c("ctrl", "trt")), factor(g$primary_outcome, 0:1), g$block)
  mh <- stats::mantelhaen.test(tb)
  expect_equal(row1(st)$p, mh$p.value)
  or_mh <- sum(tb["trt", "1", ] * tb["ctrl", "0", ] / apply(tb, 3, sum)) / sum(tb["ctrl", "1", ] * tb["trt", "0", ] / apply(tb, 3, sum))
  expect_equal(row1(st)$estimate, or_mh, tolerance = 1e-10)
  expect_gt(row1(st)$estimate, 1)
})

test_that("Cochran-Armitage on a numeric dose = prop.trend.test and rstatix (tutorial practice table)", {
  ev <- c(5, 12, 20, 28); no <- c(45, 38, 30, 22); dose <- c(0, 1, 2, 3)
  rows <- do.call(rbind, lapply(1:4, function(i) data.frame(dose = dose[i], primary_outcome = rep(c(1L, 0L), c(ev[i], no[i])))))
  rows$unit <- sprintf("U%03d", seq_len(nrow(rows)))
  d <- cat_folder(rows, "y ~ treatment", list(treatment = "dose", layout_columns = list("dose")), factors = fac("dose", as.character(dose)),
                  an = list(goals = list("trend"), primary_goal = "trend", goal_terms = list(trend = list("treatment"))))
  st <- run_cat(d)$summary$analyses$outcome$stat_test
  expect_equal(st$test, "cochran_armitage")
  ct <- stats::prop.trend.test(ev, ev + no, score = dose)
  expect_equal(row1(st)$statistic, unname(ct$statistic)); expect_lt(row1(st)$p, 1e-4)
  expect_equal(row1(st)$statistic, unname(rstatix::prop_trend_test(rbind(ev, no), score = dose)$statistic), tolerance = 1e-6)
})

test_that("two proportions: odds ratio and relative risk = epitools (Wald) and the cross-product; exact CIs in descriptives", {
  set.seed(6); rows <- data.frame(trt = rep(c("placebo", "drug"), each = 60), stringsAsFactors = FALSE)
  rows$primary_outcome <- as.integer(stats::runif(120) < ifelse(rows$trt == "drug", 0.55, 0.3)); rows$unit <- sprintf("U%03d", 1:120)
  d <- cat_folder(rows, "y ~ treatment", list(treatment = "trt", layout_columns = list("trt")), factors = fac("trt", c("placebo", "drug")))
  st <- run_cat(d)$summary$analyses$outcome$stat_test
  x <- tapply(rows$primary_outcome, factor(rows$trt, c("placebo", "drug")), sum); n <- c(60, 60)
  orr <- Filter(function(r) grepl("^odds ratio", r$statistic_name), st$table)[[1]]
  rr <- Filter(function(r) grepl("^relative risk", r$statistic_name), st$table)[[1]]
  expect_equal(orr$estimate, unname((x[2] / (n[2] - x[2])) / (x[1] / (n[1] - x[1]))), tolerance = 1e-10)
  expect_equal(rr$estimate, unname((x[2] / n[2]) / (x[1] / n[1])), tolerance = 1e-10)
  se <- sqrt(sum(1 / c(x, n - x))); expect_equal(orr$ci_low, exp(log(orr$estimate) - stats::qnorm(0.975) * se), tolerance = 1e-8)
  ds <- utils::read.delim(file.path(d, "bioana", "outcome", "descriptives.tsv"))
  expect_equal(ds$ci_low[ds$treatment == "drug"], stats::binom.test(x[["drug"]], 60)$conf.int[1])
})

test_that("nominal outcome: multinom engine = nnet; chi-squared and Cramer's V = chisq.test and rstatix; cells; Fisher when sparse", {
  skip_if_not_installed("nnet")
  set.seed(7); rows <- data.frame(line = rep(c("A", "B", "C"), each = 60), stringsAsFactors = FALSE)
  pr <- list(A = c(.6, .3, .1), B = c(.3, .4, .3), C = c(.2, .3, .5))
  rows$primary_outcome <- unlist(lapply(c("A", "B", "C"), function(l) sample(c("green", "yellow", "brown"), 60, TRUE, pr[[l]])))
  rows$unit <- sprintf("U%03d", seq_len(nrow(rows)))
  d <- cat_folder(rows, "y ~ treatment", list(treatment = "line", layout_columns = list("line")), type = "nominal", factors = fac("line", c("A", "B", "C")))
  r <- run_cat(d); st <- r$summary$analyses$outcome$stat_test
  expect_equal(st$test, "chisq_independence")
  tab <- table(factor(rows$line), factor(rows$primary_outcome, sort(unique(rows$primary_outcome))))
  cs <- stats::chisq.test(tab)
  expect_equal(row1(st)$statistic, unname(cs$statistic)); expect_equal(row1(st)$effect, rstatix::cramer_v(tab))
  cells <- utils::read.delim(file.path(d, "bioana", "outcome", "cells.tsv")); expect_equal(sum(cells$contribution_pct), 100)
  m <- nnet::multinom(primary_outcome ~ line, data = transform(rows, primary_outcome = factor(primary_outcome, sort(unique(primary_outcome)))), trace = FALSE)
  res <- r$results$outcome
  expect_equal(res$estimate[res$term == "treatmentB" & res$category == "green"], stats::coef(m)["green", "lineB"], tolerance = 1e-3)
  small <- rows[c(1:5, 61:65, 121:125), ]
  d2 <- cat_folder(small, "y ~ treatment", list(treatment = "line", layout_columns = list("line")), type = "nominal", factors = fac("line", c("A", "B", "C")))
  s2 <- run_cat(d2)$summary$analyses$outcome$stat_test
  expect_equal(s2$test, "fisher_exact"); expect_match(s2$note, "Fisher's exact test instead of chi-squared")
})

test_that("goodness of fit to an expected ratio = chisq.test(p) and the tutorial tulips (3:2:1: chi2 0.20, p 0.90)", {
  rows <- data.frame(primary_outcome = rep(c("red", "yellow", "white"), c(81, 50, 27)), stringsAsFactors = FALSE)
  rows$unit <- sprintf("U%03d", seq_len(nrow(rows)))
  d <- cat_folder(rows, "y ~ 1", list(layout_columns = list("unit")), type = "nominal", factors = list(),
                  an = list(expected_ratio = list(red = 3, yellow = 2, white = 1)))
  st <- run_cat(d)$summary$analyses$outcome$stat_test
  expect_equal(st$test, "goodness_of_fit")
  g <- row1(st); ref <- stats::chisq.test(c(81, 50, 27), p = c(3, 2, 1) / 6)
  expect_equal(g$statistic, unname(ref$statistic)); expect_equal(round(g$statistic, 2), 0.20); expect_equal(round(g$p, 2), 0.90)
})
