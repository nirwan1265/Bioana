# Correlation between traits (R/correlation.R) and agreement between raters (R/agreement.R).

meta_folder <- function(rows, an, model = "y ~ 1", terms = list(layout_columns = list("unit")), factors = list()) {
  d <- file.path(tempfile("ca_"), "design"); dir.create(d, recursive = TRUE)
  rows$sample_id <- sprintf("S%04d", seq_len(nrow(rows)))
  spec <- list(handoff_version = "0.6.6", origin = "observed", organism = "test", system = "lab", factors = factors, design = list(chosen = "crd"),
               primary_outcome = list(name = "y", type = "continuous"), analysis = c(list(model = model, terms = terms), an))
  jsonlite::write_json(spec, file.path(d, "study_spec.json"), auto_unbox = TRUE, pretty = TRUE)
  utils::write.csv(rows[intersect(c("unit", unlist(terms$layout_columns)), names(rows))], file.path(d, "field_layout.csv"), row.names = FALSE)
  utils::write.table(rows, file.path(d, "sample_metadata_template.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
  d
}
run_ca <- function(d, id) { utils::capture.output(r <- bioana_analyze(d, secondary = FALSE)); list(s = r$summary$analyses[[id]], res = r$results[[id]], r = r) }

test_that("correlation of mtcars traits = cor.test, Hmisc::rcorr and the tutorial (mpg-wt r -0.87, t -9.559, CI -0.934 to -0.744)", {
  vars <- c("mpg", "disp", "hp", "drat", "wt", "qsec")
  rows <- cbind(unit = rownames(mtcars), mtcars[vars]); rownames(rows) <- NULL
  x <- run_ca(meta_folder(rows, list(correlation = list(variables = as.list(vars)))), "correlation")
  t <- x$res[x$res$scope == "raw", ]
  mw <- t[t$var1 == "mpg" & t$var2 == "wt", ]
  expect_equal(round(mw$r, 2), -0.87); expect_equal(round(mw$statistic, 3), -9.559); expect_equal(round(c(mw$ci_low, mw$ci_high), 3), c(-0.934, -0.744))
  hr <- Hmisc::rcorr(as.matrix(mtcars[vars]))
  for (i in seq_len(nrow(t))) { expect_equal(t$r[i], hr$r[t$var1[i], t$var2[i]], tolerance = 1e-8); expect_equal(t$p[i], hr$P[t$var1[i], t$var2[i]], tolerance = 1e-6) }
  expect_equal(t$p_adj, stats::p.adjust(t$p, "holm"))
  xs <- run_ca(meta_folder(rows, list(correlation = list(variables = list("mpg", "wt"), method = "spearman"))), "correlation")
  expect_equal(round(xs$res$r[1], 2), -0.89)
  xk <- run_ca(meta_folder(rows, list(correlation = list(variables = list("mpg", "wt"), method = "kendall"))), "correlation")
  expect_equal(round(xk$res$r[1], 2), -0.73)
})

test_that("within the design: traits that both respond to a treatment correlate raw but not within treatments", {
  set.seed(9); n <- 30; trt <- rep(c("low N", "high N"), each = n)
  rows <- data.frame(unit = sprintf("P%02d", 1:(2 * n)), nitrogen = trt, height = 100 + 20 * (trt == "high N") + stats::rnorm(2 * n, 0, 5),
                     chlorophyll = 30 + 8 * (trt == "high N") + stats::rnorm(2 * n, 0, 3), stringsAsFactors = FALSE)
  x <- run_ca(meta_folder(rows, list(correlation = list(variables = list("height", "chlorophyll"))), model = "y ~ treatment",
                          terms = list(treatment = "nitrogen", layout_columns = list("unit", "nitrogen")),
                          factors = list(list(name = "nitrogen", n_levels = 2L, levels = list("low N", "high N")))), "correlation")
  raw <- x$res[x$res$scope == "raw", ]; wd <- x$res[x$res$scope == "within design", ]
  expect_gt(raw$r, 0.6); expect_lt(abs(wd$r), 0.3)
  rh <- stats::residuals(stats::lm(height ~ nitrogen, rows)); rc <- stats::residuals(stats::lm(chlorophyll ~ nitrogen, rows))
  expect_equal(wd$r, stats::cor(rh, rc), tolerance = 1e-10)
  expect_match(x$s$within_design, "residuals of each trait on treatment")
})

test_that("two raters: Cohen's kappa = vcd and psych; ordinal weighted kappa (quadratic) = psych::cohen.kappa; agreement table", {
  set.seed(2); n <- 60; truth <- sample(1:4, n, TRUE)
  r1 <- pmin(4, pmax(1, truth + sample(c(-1, 0, 0, 0, 1), n, TRUE))); r2 <- pmin(4, pmax(1, truth + sample(c(-1, 0, 0, 1), n, TRUE)))
  rows <- data.frame(unit = sprintf("L%02d", 1:n), rater_a = r1, rater_b = r2)
  x <- run_ca(meta_folder(rows, list(agreement = list(raters = list("rater_a", "rater_b"), scale = "ordinal"))), "agreement")
  m <- x$res
  ck <- psych::cohen.kappa(cbind(r1, r2))
  expect_equal(m$estimate[m$measure == "Cohen's kappa (unweighted)"], ck$kappa, tolerance = 1e-8)
  expect_equal(m$estimate[grepl("quadratic", m$measure)], ck$weighted.kappa, tolerance = 1e-8)
  vk <- vcd::Kappa(table(factor(r1, 1:4), factor(r2, 1:4)))
  expect_equal(m$estimate[grepl("linear", m$measure)], vk$Weighted[["value"]], tolerance = 1e-10)
  expect_equal(x$s$primary$measure, "weighted kappa (linear, Cicchetti-Allison)")
  expect_true(file.exists(file.path(dirname(x$r$paths[["agreement/agreement"]]), "agreement_table.tsv")))
})

test_that("three raters (long format): Light's kappa = mean pairwise Cohen without irr; continuous ratings: ICC = psych and the Shrout-Fleiss ICC2 formula", {
  set.seed(5); n <- 25; truth <- sample(c("healthy", "mild", "severe"), n, TRUE)
  flip <- function() ifelse(stats::runif(n) < 0.2, sample(c("healthy", "mild", "severe"), n, TRUE), truth)
  long <- data.frame(unit = rep(sprintf("P%02d", 1:n), 3), rater = rep(c("A", "B", "C"), each = n), score = c(flip(), flip(), flip()), stringsAsFactors = FALSE)
  x <- run_ca(meta_folder(long, list(agreement = list(rater_column = "rater", score_column = "score", unit_column = "unit", scale = "nominal"))), "agreement")
  W <- stats::reshape(long, idvar = "unit", timevar = "rater", direction = "wide")[-1]
  pk <- utils::combn(3, 2, function(ij) vcd::Kappa(table(factor(W[[ij[1]]], c("healthy", "mild", "severe")), factor(W[[ij[2]]], c("healthy", "mild", "severe"))))$Unweighted[["value"]])
  expect_equal(x$res$estimate[grepl("Light", x$res$measure)], mean(pk), tolerance = 1e-10)
  # continuous: 3 raters measure leaf length
  set.seed(8); L <- stats::rnorm(20, 50, 8)
  rows <- data.frame(unit = sprintf("L%02d", 1:20), r1 = L + stats::rnorm(20, 0, 2), r2 = L + 1 + stats::rnorm(20, 0, 2), r3 = L - 1 + stats::rnorm(20, 0, 2))
  y <- run_ca(meta_folder(rows, list(agreement = list(raters = list("r1", "r2", "r3"), scale = "continuous"))), "agreement")
  M <- as.matrix(rows[c("r1", "r2", "r3")]); n <- 20; k <- 3
  ms <- function(ss, df) ss / df; gm <- mean(M)
  MSR <- ms(k * sum((rowMeans(M) - gm)^2), n - 1); MSC <- ms(n * sum((colMeans(M) - gm)^2), k - 1)
  MSE <- ms(sum((M - outer(rowMeans(M), colMeans(M), "+") + gm)^2), (n - 1) * (k - 1))
  icc2 <- (MSR - MSE) / (MSR + (k - 1) * MSE + k * (MSC - MSE) / n)
  expect_equal(y$s$primary$estimate, icc2, tolerance = 1e-10)
  expect_equal(nrow(y$res), 6)
})
