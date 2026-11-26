# Classical ANOVA next to the model (R/stat_tests.R): equality with rstatix, base aov, car, paired t-tests.

rm_folder <- function(mixed = TRUE, drop_one = FALSE, posthoc = NULL, seed = 3) {
  d <- file.path(tempfile("rm_"), "design"); dir.create(d, recursive = TRUE)
  days <- c("d0", "d7", "d14", "d21"); n <- 8
  plants <- sprintf("P%02d", seq_len(if (mixed) 2 * n else n))
  water <- if (mixed) rep(c("well-watered", "drought"), each = n) else rep("well-watered", n)
  set.seed(seed)
  rows <- do.call(rbind, lapply(seq_along(plants), function(i) {
    e <- as.numeric(stats::arima.sim(list(ar = 0.8), length(days), sd = 1))     # correlated days: sphericity fails
    data.frame(plant = plants[i], water = water[i], day = days, primary_outcome = 20 + stats::rnorm(1, 0, 2) +
                 (water[i] == "drought") * c(0, -1, -3, -5) + c(0, 2, 4, 5) + e, stringsAsFactors = FALSE)
  }))
  if (drop_one) rows <- rows[-nrow(rows), ]
  rows$sample_id <- sprintf("S%03d", seq_len(nrow(rows)))
  model <- if (mixed) "y ~ treatment * time + (1 | unit)" else "y ~ time + (1 | unit)"
  terms <- if (mixed) list(treatment = "water", time = "day", unit = "plant", layout_columns = list("water", "day", "plant"))
           else list(time = "day", unit = "plant", layout_columns = list("day", "plant"))
  an <- list(model = model, terms = terms)
  if (!is.null(posthoc)) an$posthoc <- posthoc
  spec <- list(handoff_version = "0.6.4", origin = "observed", organism = "test", system = "greenhouse", experimental_unit = "plant",
               repeated_factors = list("day"),
               factors = c(if (mixed) list(list(name = "water", n_levels = 2L, levels = list("well-watered", "drought"))),
                           list(list(name = "day", n_levels = 4L, levels = as.list(days)))),
               design = list(chosen = "repeated_measures", reps = n), primary_outcome = list(name = "height", type = "continuous"), analysis = an)
  jsonlite::write_json(spec, file.path(d, "study_spec.json"), auto_unbox = TRUE, pretty = TRUE)
  utils::write.csv(rows[, c("sample_id", "plant", "water", "day")], file.path(d, "field_layout.csv"), row.names = FALSE)
  utils::write.table(rows, file.path(d, "sample_metadata_template.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
  bioana_write_manifest(d, list(id = "outcome", data_type = "outcome", files = list(table = "sample_metadata_template.tsv"), column = "primary_outcome"))
  rows$plant <- factor(rows$plant); rows$day <- factor(rows$day, levels = days); rows$water <- factor(rows$water, levels = c("well-watered", "drought"))
  list(dir = d, data = rows)
}

run_st <- function(f) { utils::capture.output(r <- bioana_analyze(f$dir, secondary = FALSE)); r$summary$analyses$outcome$stat_test }
row_of <- function(st, term) st$table[[which(vapply(st$table, `[[`, "", "term") == term)]]

test_that("repeated-measures ANOVA = rstatix::anova_test (GG) and base aov(Error()); Mauchly and post-hoc = rstatix", {
  skip_if_not_installed("afex"); skip_if_not_installed("emmeans"); skip_if_not_installed("rstatix")
  f <- rm_folder(mixed = FALSE)
  st <- run_st(f)
  expect_equal(st$status, "run"); expect_equal(st$layout, "within"); expect_equal(st$correction, "GG")
  ref <- rstatix::anova_test(f$data, dv = primary_outcome, wid = plant, within = day)
  gg <- rstatix::get_anova_table(ref, correction = "GG")
  b <- row_of(st, "time")
  # rstatix rounds: F, ges to 3 decimals, p to 3 decimals (or 3 significant digits when small), df to 2
  expect_equal(b$F, gg$F, tolerance = 1e-3); expect_lt(abs(b$num_df - gg$DFn), 0.006)
  expect_lt(abs(b$den_df - gg$DFd), 0.06); expect_lt(abs(b$p - gg$p), 5e-4); expect_lt(abs(b$ges - gg$ges), 5e-4)
  expect_equal(st$sphericity[[1]]$mauchly_p, ref$`Mauchly's Test for Sphericity`$p, tolerance = 1e-3)
  # uncorrected F = the classical Error() ANOVA
  a <- summary(stats::aov(primary_outcome ~ day + Error(plant / day), data = f$data))
  expect_equal(b$F, a[["Error: plant:day"]][[1]]["day", "F value"], tolerance = 1e-8)
  # post-hoc (multivariate model, no between factors) = paired t-tests, Holm
  pt <- rstatix::pairwise_t_test(f$data, primary_outcome ~ day, paired = TRUE, p.adjust.method = "holm")
  ph <- utils::read.delim(file.path(f$dir, "bioana", "outcome", "posthoc.tsv"))
  key <- function(a, b) paste(pmin(a, b), pmax(a, b))
  m <- match(key(pt$group1, pt$group2), key(sub(" - .*", "", ph$contrast), sub(".* - ", "", ph$contrast)))
  expect_false(anyNA(m))
  expect_equal(ph$p[m], pt$p, tolerance = 1e-3); expect_equal(ph$p_adj[m], pt$p.adj, tolerance = 1e-3)
})

test_that("mixed ANOVA = rstatix::anova_test(between, within) with GG; simple effects of day within water for time_course", {
  skip_if_not_installed("afex"); skip_if_not_installed("emmeans"); skip_if_not_installed("rstatix")
  f <- rm_folder(mixed = TRUE)
  st <- run_st(f)
  expect_equal(st$layout, "mixed")
  gg <- rstatix::get_anova_table(rstatix::anova_test(f$data, dv = primary_outcome, wid = plant, between = water, within = day), correction = "GG")
  for (t in c("treatment", "time", "treatment:time")) {
    rs <- c(treatment = "water", time = "day", `treatment:time` = "water:day")[[t]]
    expect_equal(row_of(st, t)$F, gg$F[gg$Effect == rs], tolerance = 1e-3)
    expect_lt(abs(row_of(st, t)$p - gg$p[gg$Effect == rs]), 5e-4)
    expect_lt(abs(row_of(st, t)$ges - gg$ges[gg$Effect == rs]), 5e-4)
  }
  expect_true(is.numeric(st$checks$levene_min_p))
  # goal time_course -> pairs of days within each water level
  s <- jsonlite::read_json(file.path(f$dir, "study_spec.json"))
  s$analysis$goals <- list("time_course"); s$analysis$primary_goal <- "time_course"; s$analysis$goal_terms <- list(time_course = list("time"))
  jsonlite::write_json(s, file.path(f$dir, "study_spec.json"), auto_unbox = TRUE, pretty = TRUE)
  st <- run_st(f)
  expect_equal(unlist(st$posthoc$families), "time | treatment")
  expect_equal(st$posthoc$n_comparisons, 12)
})

test_that("factorial ANOVA keeps the model's additive block term and equals car::Anova type III", {
  skip_if_not_installed("afex"); skip_if_not_installed("car")
  d <- file.path(tempfile("rcbd_"), "design"); dir.create(d, recursive = TRUE)
  g <- expand.grid(block = paste0("B", 1:4), trt = c("ctrl", "N1", "N2"), stringsAsFactors = FALSE)
  set.seed(5); g$primary_outcome <- 10 + (g$trt == "N2") * 2 + as.integer(factor(g$block)) * 0.5 + stats::rnorm(nrow(g))
  g$sample_id <- sprintf("S%02d", seq_len(nrow(g)))
  spec <- list(handoff_version = "0.6.4", origin = "observed", organism = "test", system = "field",
               factors = list(list(name = "trt", n_levels = 3L, levels = list("ctrl", "N1", "N2"))), design = list(chosen = "rcbd", reps = 4L),
               primary_outcome = list(name = "yield", type = "continuous"),
               analysis = list(model = "y ~ treatment + block", terms = list(treatment = "trt", block = "block", layout_columns = list("trt", "block")),
                               posthoc = list(type = "vs_control", control = "ctrl")))
  jsonlite::write_json(spec, file.path(d, "study_spec.json"), auto_unbox = TRUE, pretty = TRUE)
  utils::write.csv(g[, c("sample_id", "block", "trt")], file.path(d, "field_layout.csv"), row.names = FALSE)
  utils::write.table(g, file.path(d, "sample_metadata_template.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
  bioana_write_manifest(d, list(id = "outcome", data_type = "outcome", files = list(table = "sample_metadata_template.tsv"), column = "primary_outcome"))
  st <- run_st(list(dir = d))
  expect_equal(st$layout, "between"); expect_equal(length(st$table), 2)
  op <- options(contrasts = c("contr.sum", "contr.poly")); on.exit(options(op))
  ca <- car::Anova(stats::lm(primary_outcome ~ trt + block, data = g), type = 3)
  expect_equal(row_of(st, "treatment")$F, ca["trt", "F value"], tolerance = 1e-8)
  expect_equal(row_of(st, "block")$p, ca["block", "Pr(>F)"], tolerance = 1e-8)
  expect_equal(st$posthoc$type, "vs_control"); expect_equal(st$posthoc$adjust, "dunnettx"); expect_equal(st$posthoc$n_comparisons, 2)
})

test_that("a subject with a missing day: ANOVA not run (would drop it), the mixed model still fits; Girden rule switches to HF", {
  skip_if_not_installed("afex")
  f <- rm_folder(mixed = FALSE, drop_one = TRUE)
  utils::capture.output(r <- bioana_analyze(f$dir, secondary = FALSE))
  st <- r$summary$analyses$outcome$stat_test
  expect_equal(st$status, "not run"); expect_match(st$note, "1 of 8 unit levels miss some time cells")
  expect_equal(r$summary$analyses$outcome$engine$engine, "lmer")
  # girden: HF when the GG epsilon > 0.75
  f2 <- rm_folder(mixed = FALSE, seed = 9)
  fit <- suppressMessages(afex::aov_car(primary_outcome ~ Error(plant / day), data = f2$data))
  tab <- bioana:::.ba_anova_corrected(fit, "girden", 0.75)
  eps <- suppressWarnings(summary(fit$Anova, multivariate = FALSE))$pval.adjustments["day", "GG eps"]
  expect_equal(tab$correction, if (eps > 0.75) "HF" else "GG")
})
