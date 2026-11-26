# Tests on one treatment factor (R/stat_tests.R .ba_st_simple): equality with stats, rstatix, effectsize, coin.

st_folder <- function(kind = "between", k = 2, n = 10, type = "continuous", seed = 1, an = list(), shift = 1) {
  d <- file.path(tempfile("st_"), "design"); dir.create(d, recursive = TRUE)
  lev <- paste0("L", seq_len(k)); set.seed(seed)
  rows <- switch(kind,
    between = data.frame(cond = rep(lev, each = n), unit = sprintf("U%03d", seq_len(k * n))),
    within  = data.frame(cond = rep(lev, times = n), unit = rep(sprintf("P%02d", seq_len(n)), each = k)),
    rcbd    = data.frame(cond = rep(lev, times = n), block = rep(sprintf("B%02d", seq_len(n)), each = k)))
  rows$cond <- as.character(rows$cond)
  g <- match(rows$cond, lev) - 1
  base <- if (kind == "between") 0 else stats::rnorm(n, 0, 1)[rep(seq_len(n), each = k)]
  z <- 10 + shift * g + base + stats::rexp(nrow(rows), 1)                   # skewed noise
  rows$primary_outcome <- switch(type, continuous = z, ordinal = as.integer(cut(z, stats::quantile(z, 0:5 / 5), include.lowest = TRUE)),
                                 binary = as.integer(stats::runif(nrow(rows)) < stats::plogis(-1 + shift * g)))
  rows$sample_id <- sprintf("S%03d", seq_len(nrow(rows)))
  model <- switch(kind, between = "y ~ treatment", within = "y ~ treatment + (1 | unit)", rcbd = "y ~ treatment + block")
  terms <- switch(kind, between = list(treatment = "cond", layout_columns = list("cond")),
                  within = list(treatment = "cond", unit = "unit", layout_columns = list("cond", "unit")),
                  rcbd = list(treatment = "cond", block = "block", layout_columns = list("cond", "block")))
  spec <- list(handoff_version = "0.6.4", origin = "observed", organism = "test", system = "lab",
               factors = list(list(name = "cond", n_levels = k, levels = as.list(lev))), design = list(chosen = kind, reps = n),
               primary_outcome = list(name = "y", type = type), analysis = c(list(model = model, terms = terms), an))
  jsonlite::write_json(spec, file.path(d, "study_spec.json"), auto_unbox = TRUE, pretty = TRUE)
  utils::write.csv(rows[, setdiff(names(rows), "primary_outcome")], file.path(d, "field_layout.csv"), row.names = FALSE)
  utils::write.table(rows, file.path(d, "sample_metadata_template.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
  bioana_write_manifest(d, list(id = "outcome", data_type = "outcome", files = list(table = "sample_metadata_template.tsv"), column = "primary_outcome"))
  rows$cond <- factor(rows$cond, lev)
  list(dir = d, data = rows, lev = lev)
}
run_one <- function(f) { out <- utils::capture.output(r <- bioana_analyze(f$dir, secondary = FALSE)); list(st = r$summary$analyses$outcome$stat_test, out = out) }
first <- function(st) st$table[[1]]
also <- function(st, id) Filter(function(x) identical(x$test, id), st$also)[[1]]

test_that("two independent groups: Welch t = t.test and rstatix; Cohen's d (Welch SD) = effectsize; Student, Wilcoxon on request", {
  skip_if_not_installed("rstatix"); skip_if_not_installed("effectsize")
  f <- st_folder("between", an = list(stat_tests = list("t_test_student", "nonparametric")))
  r <- run_one(f); st <- r$st
  expect_equal(st$test, "t_test_welch")
  tt <- stats::t.test(primary_outcome ~ cond, data = f$data)
  expect_equal(first(st)$statistic, unname(tt$statistic)); expect_equal(first(st)$df1, unname(tt$parameter)); expect_equal(first(st)$p, tt$p.value)
  ef <- effectsize::cohens_d(primary_outcome ~ cond, data = f$data, pooled_sd = FALSE)
  expect_equal(first(st)$effect, ef$Cohens_d, tolerance = 1e-6)
  rt <- rstatix::t_test(f$data, primary_outcome ~ cond); expect_equal(first(st)$p, rt$p, tolerance = 1e-3)
  stu <- also(st, "t_test_student"); ts <- stats::t.test(primary_outcome ~ cond, data = f$data, var.equal = TRUE)
  expect_equal(first(stu)$p, ts$p.value)
  expect_equal(first(stu)$effect, effectsize::cohens_d(primary_outcome ~ cond, data = f$data)$Cohens_d, tolerance = 1e-6)
  wx <- also(st, "wilcoxon_rank_sum"); w <- stats::wilcox.test(primary_outcome ~ cond, data = f$data)
  expect_equal(first(wx)$statistic, unname(w$statistic)); expect_equal(first(wx)$p, w$p.value)
  z <- coin::statistic(coin::wilcox_test(primary_outcome ~ cond, data = f$data), type = "standardized")
  expect_equal(first(wx)$effect, abs(as.numeric(z)) / sqrt(nrow(f$data)), tolerance = 1e-6)
  expect_true(file.exists(file.path(f$dir, "bioana", "outcome", "stat_tests.tsv")))
  expect_true(any(grepl("t_test_welch \\[stats", r$out)))
})

test_that("two levels on the same subjects: paired t = t.test(paired) and F of the RM ANOVA = t^2; d_z; signed-rank and sign test on request", {
  skip_if_not_installed("rstatix"); skip_if_not_installed("afex")
  f <- st_folder("within", an = list(stat_tests = list("anova_within", "wilcoxon_signed_rank", "sign_test")))
  st <- run_one(f)$st
  expect_equal(st$test, "t_test_paired")
  w <- reshape(f$data[, c("unit", "cond", "primary_outcome")], idvar = "unit", timevar = "cond", direction = "wide")
  x <- w$primary_outcome.L1; y <- w$primary_outcome.L2
  tt <- stats::t.test(x, y, paired = TRUE)
  expect_equal(first(st)$statistic, unname(tt$statistic)); expect_equal(first(st)$p, tt$p.value)
  expect_equal(first(st)$effect, mean(x - y) / stats::sd(x - y), tolerance = 1e-8)
  expect_equal(also(st, "anova_within")$table[[1]]$F, unname(tt$statistic)^2, tolerance = 1e-8)
  expect_equal(first(also(st, "wilcoxon_signed_rank"))$p, stats::wilcox.test(x, y, paired = TRUE)$p.value)
  sgn <- rstatix::sign_test(f$data[order(f$data$unit), ], primary_outcome ~ cond)
  expect_lt(abs(first(also(st, "sign_test"))$p - sgn$p), 5e-4)   # rstatix rounds p to 3 significant digits
})

test_that("ordinal outcome, three groups: Kruskal-Wallis = kruskal.test, eta2[H] by formula, Dunn post-hoc; Welch ANOVA on request = oneway.test", {
  skip_if_not_installed("rstatix")
  f <- st_folder("between", k = 3, n = 12, type = "ordinal")
  st <- run_one(f)$st
  expect_equal(st$test, "kruskal_wallis")
  kt <- stats::kruskal.test(primary_outcome ~ cond, data = f$data)
  expect_equal(first(st)$statistic, unname(kt$statistic)); expect_equal(first(st)$p, kt$p.value)
  expect_equal(first(st)$effect, (unname(kt$statistic) - 3 + 1) / (nrow(f$data) - 3), tolerance = 1e-8)
  expect_equal(st$posthoc$n_comparisons, 3)
  f2 <- st_folder("between", k = 3, n = 12, an = list(stat_tests = list("welch_anova")))
  st2 <- run_one(f2)$st
  expect_equal(st2$test, "anova_between")
  ow <- stats::oneway.test(primary_outcome ~ cond, data = f2$data, var.equal = FALSE)
  expect_equal(first(also(st2, "welch_anova"))$p, ow$p.value)
})

test_that("Friedman on subjects and on RCBD blocks = friedman.test; Kendall's W; post-hoc = pairwise.wilcox.test(paired)", {
  skip_if_not_installed("rstatix")
  for (kind in c("within", "rcbd")) {
    f <- st_folder(kind, k = 3, n = 9, type = "ordinal")
    st <- run_one(f)$st
    expect_equal(st$test, "friedman")
    b <- if (kind == "within") "unit" else "block"
    ft <- stats::friedman.test(stats::as.formula(paste("primary_outcome ~ cond |", b)), data = f$data)
    expect_equal(first(st)$statistic, unname(ft$statistic)); expect_equal(first(st)$p, ft$p.value)
    expect_equal(first(st)$effect, unname(ft$statistic) / (9 * (3 - 1)), tolerance = 1e-8)
    o <- f$data[order(f$data[[b]]), ]
    pw <- suppressWarnings(stats::pairwise.wilcox.test(o$primary_outcome, o$cond, paired = TRUE, p.adjust.method = "holm"))
    ph <- utils::read.delim(file.path(f$dir, "bioana", "outcome", "posthoc.tsv"))
    expect_equal(sort(ph$p_adj), sort(stats::na.omit(c(pw$p.value))), tolerance = 1e-8)
  }
})

test_that("binary outcome: proportion test = prop.test with Cohen's h; Fisher when expected counts are small; pairwise for three groups", {
  skip_if_not_installed("effectsize")
  f <- st_folder("between", n = 40, type = "binary", shift = 1.5)
  st <- run_one(f)$st
  expect_equal(st$test, "prop_test")
  x <- tapply(f$data$primary_outcome, f$data$cond, sum); n <- tapply(f$data$primary_outcome, f$data$cond, length)
  pt <- stats::prop.test(x, n)
  expect_equal(first(st)$p, pt$p.value)
  p1 <- x[1] / n[1]; p2 <- x[2] / n[2]
  expect_equal(abs(first(st)$effect), abs(unname(2 * asin(sqrt(p1)) - 2 * asin(sqrt(p2)))), tolerance = 1e-8)
  small <- st_folder("between", n = 6, type = "binary", seed = 4)
  s2 <- run_one(small)$st
  expect_equal(s2$test, "fisher_exact"); expect_match(s2$note, "Fisher's exact test instead")
  expect_equal(first(s2)$p, stats::fisher.test(table(small$data$cond, small$data$primary_outcome))$p.value)
  f3 <- st_folder("between", k = 3, n = 40, type = "binary", seed = 2)
  s3 <- run_one(f3)$st
  expect_equal(s3$posthoc$n_comparisons, 3)
})

test_that("descriptives = rstatix get_summary_stats('common'); a requested test that does not fit is reported", {
  skip_if_not_installed("rstatix"); skip_if_not_installed("dplyr")
  f <- st_folder("between", k = 3, an = list(stat_tests = list("t_test_paired", "made_up")))
  r <- run_one(f)
  ds <- utils::read.delim(file.path(f$dir, "bioana", "outcome", "descriptives.tsv"))
  ref <- as.data.frame(rstatix::get_summary_stats(dplyr::group_by(f$data, cond), primary_outcome, type = "common"))
  m <- match(ref$cond, ds$treatment)
  for (k in c("n", "mean", "sd", "se", "ci", "median", "iqr", "min", "max")) expect_equal(ds[[k]][m], ref[[k]], tolerance = 1e-3)
  expect_match(r$st$note, "'t_test_paired' does not fit")
  expect_match(r$st$note, "'made_up' is not in stat_tests.tsv")
})
