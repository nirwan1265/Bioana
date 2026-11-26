# Analysis goals (docs/contracts/analysis_goals.md): the goal picks the test; joint tests must equal
# the package called directly.

tc_folder <- function(goals = NULL, primary = NULL, goal_terms = NULL, assumptions = NULL, reps = 3, seed = 7, G = 200) {
  d <- file.path(tempfile("goal_"), "design"); dir.create(d, recursive = TRUE)
  grid <- expand.grid(rep = seq_len(reps), time = c("0", "6", "24"), treatment = c("control", "drought"), stringsAsFactors = FALSE)
  grid$sample_id <- sprintf("S%02d", seq_len(nrow(grid)))
  an <- list(model = "y ~ treatment * time", terms = list(treatment = "treatment", time = "time", layout_columns = list("treatment", "time")))
  if (!is.null(goals)) an$goals <- as.list(goals)
  if (!is.null(primary)) an$primary_goal <- primary
  if (!is.null(goal_terms)) an$goal_terms <- goal_terms
  spec <- list(handoff_version = "0.5.1", origin = "observed", organism = "Zea mays", system = "growth chamber",
               factors = list(list(name = "treatment", n_levels = 2L, levels = list("control", "drought")),
                              list(name = "time", n_levels = 3L, levels = list("0", "6", "24"))),
               design = list(chosen = "factorial", reps = reps), assays = list("rnaseq"),
               primary_outcome = list(name = "root length", type = "continuous"), analysis = an)
  if (!is.null(assumptions)) spec$assumptions <- assumptions
  jsonlite::write_json(spec, file.path(d, "study_spec.json"), auto_unbox = TRUE, pretty = TRUE)
  set.seed(seed)
  trt <- grid$treatment == "drought"; tm <- match(grid$time, c("0", "6", "24")) - 1
  grid$primary_outcome <- 10 + 1 * trt + 0.5 * tm + 1.2 * trt * tm + stats::rnorm(nrow(grid), 0, 0.6)
  utils::write.csv(grid[, c("sample_id", "treatment", "time", "rep")], file.path(d, "field_layout.csv"), row.names = FALSE)
  utils::write.table(grid, file.path(d, "sample_metadata_template.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
  mu <- 2^stats::rnorm(G, 7, 1)
  inter <- c(rep(1.5, 20), rep(0, G - 20)); timeonly <- c(rep(0, 20), rep(1.2, 20), rep(0, G - 40))
  cnt <- sapply(seq_len(nrow(grid)), function(j) stats::rnbinom(G, mu = mu * 2^(inter * trt[j] * (tm[j] > 0) + timeonly * (tm[j] > 0)), size = 15))
  dimnames(cnt) <- list(sprintf("g%03d", 1:G), grid$sample_id)
  write_mat(cnt, file.path(d, "counts.tsv"))
  bioana_write_manifest(d, list(id = "outcome", data_type = "outcome", files = list(table = "sample_metadata_template.tsv"), column = "primary_outcome"),
                        list(id = "rna", data_type = "feature_counts", assay = "rnaseq", files = list(counts = "counts.tsv")))
  list(dir = d, counts = cnt, meta = grid)
}
cd_of <- function(x) data.frame(treatment = factor(x$meta$treatment, levels = c("control", "drought")),
                                time = factor(x$meta$time, levels = c("0", "6", "24")), row.names = x$meta$sample_id)

test_that("time_course with DESeq2 = DESeq2's LRT called directly (reduced = ~ treatment)", {
  skip_if_not_installed("DESeq2")
  x <- tc_folder(goals = "time_course", goal_terms = list(time_course = list("time")))
  utils::capture.output(r <- bioana_analyze(x$dir, only = "rna", secondary = FALSE))
  a <- r$summary$analyses$rna
  expect_equal(a$goal$id, "time_course"); expect_equal(a$goal$test, "lrt"); expect_equal(a$engine$package, "DESeq2")
  dds <- DESeq2::DESeqDataSetFromMatrix(x$counts, cd_of(x), ~ treatment * time)
  dds <- suppressMessages(DESeq2::DESeq(dds, test = "LRT", reduced = ~ treatment, quiet = TRUE))
  g <- DESeq2::results(dds, alpha = 0.05)
  t <- r$results$rna[match(rownames(g), r$results$rna$feature_id), ]
  expect_equal(t$statistic, g$stat, tolerance = 1e-6)
  expect_equal(t$p, g$pvalue, tolerance = 1e-6)
  expect_identical(is.na(t$padj), is.na(g$padj))
  expect_gt(mean(t$padj[1:40] < 0.05, na.rm = TRUE), 0.8)        # genes that change over time are found
})

test_that("interaction: DESeq2 LRT, edgeR QL F and limma-voom F equal each package called directly", {
  skip_if_not_installed("DESeq2"); skip_if_not_installed("edgeR"); skip_if_not_installed("limma")
  x <- tc_folder(goals = c("interaction", "treatment_effect"), primary = "interaction",
                 goal_terms = list(interaction = list("treatment", "time")),
                 assumptions = list(analysis.goals = "assumed: interaction (the question does not say which comparison)"))
  out <- utils::capture.output(r <- bioana_analyze(x$dir, only = "rna", secondary = FALSE))
  a <- r$summary$analyses$rna
  expect_equal(a$goal$id, "interaction"); expect_true(a$goal$assumed_by_biodes)
  expect_equal(unlist(a$goal$coefficients), c("treatmentdrought:time6", "treatmentdrought:time24"))
  expect_true(any(grepl("\\{goal interaction, joint test, assumed by Biodes\\}", out)))
  cd <- cd_of(x)
  dds <- suppressMessages(DESeq2::DESeq(DESeq2::DESeqDataSetFromMatrix(x$counts, cd, ~ treatment * time), test = "LRT", reduced = ~ treatment + time, quiet = TRUE))
  g <- DESeq2::results(dds, alpha = 0.05)
  t <- r$results$rna[match(rownames(g), r$results$rna$feature_id), ]
  expect_equal(t$statistic, g$stat, tolerance = 1e-6)
  expect_gt(mean(t$padj[1:20] < 0.05, na.rm = TRUE), 0.6)       # interaction genes found (3 reps per cell)
  expect_lt(mean(t$padj[21:40] < 0.05, na.rm = TRUE), 0.25)      # time-only genes are not interaction
  X <- stats::model.matrix(~ treatment * time, cd)
  utils::capture.output(r2 <- bioana_analyze(x$dir, only = "rna", secondary = FALSE, params = list(feature_counts = list(engine = "edger"))))
  y <- bioana:::.ba_edger_norm(edgeR::DGEList(x$counts))
  fit <- edgeR::glmQLFit(edgeR::estimateDisp(y, X), X)
  ge <- edgeR::topTags(edgeR::glmQLFTest(fit, coef = 5:6), n = Inf, sort.by = "none")$table
  t2 <- r2$results$rna[match(rownames(ge), r2$results$rna$feature_id), ]
  expect_equal(t2$statistic, ge$F, tolerance = 1e-8); expect_equal(t2$padj, ge$FDR, tolerance = 1e-8)
  utils::capture.output(r3 <- bioana_analyze(x$dir, only = "rna", secondary = FALSE, params = list(feature_counts = list(engine = "limma_voom"))))
  v <- limma::voom(y, X)
  gl <- limma::topTable(limma::eBayes(limma::lmFit(v, X)), coef = 5:6, number = Inf, sort.by = "none")
  t3 <- r3$results$rna[match(rownames(gl), r3$results$rna$feature_id), ]
  expect_equal(t3$statistic, gl$F, tolerance = 1e-8)
})

test_that("outcome goals: interaction LRT = anova(lm); trend = the slope of a numeric term", {
  x <- tc_folder(goals = "interaction", goal_terms = list(interaction = list("treatment", "time")))
  out <- utils::capture.output(r <- bioana_analyze(x$dir, only = "outcome", secondary = FALSE))
  a <- r$summary$analyses$outcome
  fr <- data.frame(y = x$meta$primary_outcome, cd_of(x))
  av <- stats::anova(stats::lm(y ~ treatment + time, fr), stats::lm(y ~ treatment * time, fr))
  expect_equal(a$goal$result$statistic, av$F[2], tolerance = 1e-8)
  expect_equal(a$goal$result$p, av$`Pr(>F)`[2], tolerance = 1e-8)
  expect_true(file.exists(file.path(x$dir, "bioana", "outcome", "goal_test.tsv")))
  expect_true(any(grepl("interaction treatment:time \\(2 coefficients\\): statistic", out)))
  # trend: time as a numeric term
  s <- jsonlite::read_json(file.path(x$dir, "study_spec.json"))
  s$analysis$model <- "y ~ treatment + time"; s$analysis$goals <- list("trend"); s$analysis$goal_terms <- list(trend = list("time"))
  jsonlite::write_json(s, file.path(x$dir, "study_spec.json"), auto_unbox = TRUE)
  utils::capture.output(r <- bioana_analyze(x$dir, only = "outcome", secondary = FALSE))
  a <- r$summary$analyses$outcome
  expect_equal(a$goal$id, "trend"); expect_equal(unlist(a$goal$coefficients), "time")
  fr$time_num <- as.numeric(as.character(x$meta$time))
  co <- stats::coef(summary(stats::lm(y ~ treatment + time_num, fr)))["time_num", ]
  res <- r$results$outcome
  expect_equal(res$estimate[res$term == "time"], unname(co["Estimate"]), tolerance = 1e-8)
  expect_false(is.na(res$padj[res$term == "time"]))
})

test_that("goal resolution: default with a note, goals that do not apply, explicit goals, joint tests skip Wald-only engines", {
  x <- tc_folder()
  out <- utils::capture.output(r <- bioana_analyze(x$dir, only = "outcome", secondary = FALSE))
  g <- r$summary$analyses$outcome$goal
  expect_equal(g$id, "treatment_effect"); expect_equal(g$source, "default")
  expect_match(g$note, "no analysis.goals")
  expect_true(any(grepl("goal: no analysis.goals", out)))
  y <- tc_folder(goals = c("cell_composition"))
  utils::capture.output(r <- bioana_analyze(y$dir, only = "outcome", secondary = FALSE))
  expect_match(r$summary$analyses$outcome$goal$note, "do not apply to this data set")
  utils::capture.output(r <- bioana_analyze(y$dir, only = "outcome", secondary = FALSE, params = list(outcome = list(goal = "ranking"))))
  expect_match(r$summary$not_analysed$outcome$note, "not available for outcome")
  design <- bioana:::.ba_read_biodes(x$dir)
  design$goal <- list(id = "time_course")
  ch <- tryCatch(bioana:::.ba_engine_chain("bisulfite_counts", design, c("values", "bisulfite"), test = "lrt"), error = function(e) e)
  if (!inherits(ch, "error")) {
    expect_false("dss" %in% ch$chain)
    expect_match(ch$skipped$dss, "no lrt test")
  }
  gl <- bioana_goals()
  expect_true(all(c("treatment_effect", "interaction", "time_course", "trend", "association") %in% gl$goal))
})

trend_folder <- function(factor_extra = list(), effect = NULL, truth_x = NULL, seed = 4) {
  d <- file.path(tempfile("trend_"), "design"); dir.create(d, recursive = TRUE)
  lev <- c("0 mg", "1 mg", "10 mg", "100 mg"); val <- c(0, 1, 10, 100)
  grid <- expand.grid(rep = 1:4, dose = lev, stringsAsFactors = FALSE); grid$sample_id <- sprintf("S%02d", seq_len(nrow(grid)))
  x <- val[match(grid$dose, lev)]
  set.seed(seed); grid$primary_outcome <- 5 + 0.8 * log10(x + 1) + stats::rnorm(nrow(grid), 0, 0.3)
  an <- list(model = "y ~ treatment", terms = list(treatment = "dose", layout_columns = list("dose")),
             goals = list("trend"), primary_goal = "trend", goal_terms = list(trend = list("treatment")))
  if (!is.null(effect)) an$goal_effect <- list(trend = effect)
  spec <- list(handoff_version = "0.6.2", origin = "observed", organism = "test", system = "lab",
               factors = list(utils::modifyList(list(name = "dose", n_levels = 4L, levels = as.list(lev), values = as.list(val), unit = "mg", scale = "log"), factor_extra)),
               design = list(chosen = "crd", reps = 4L), primary_outcome = list(name = "growth", type = "continuous"), analysis = an,
               assumptions = if (is.null(effect$size)) list(analysis.goal_effect.trend.size = "not stated"))
  jsonlite::write_json(spec, file.path(d, "study_spec.json"), auto_unbox = TRUE, pretty = TRUE)
  utils::write.csv(grid[, c("sample_id", "dose", "rep")], file.path(d, "field_layout.csv"), row.names = FALSE)
  utils::write.table(grid, file.path(d, "sample_metadata_template.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
  e <- list(id = "outcome", data_type = "outcome", files = list(table = "sample_metadata_template.tsv"), column = "primary_outcome")
  if (!is.null(truth_x)) {
    # the shape Simitall writes for a trend: per-level effects on the predictor plus goal_truth
    jsonlite::write_json(list(intercept = 5, treatment_effects = list(treatment = stats::setNames(as.list(0.8 * truth_x), lev)),
                              numeric_effects = list(), block_effects = list(), random_sd = list(), residual_sd = 0.3, goal_truth = list(goal = "trend", model_terms = list("treatment"), status = "planted", true_effect = 0.8,
                                                               values = as.list(val), x = as.list(truth_x), scale = "log")),
                         file.path(d, "sim_truth.json"), auto_unbox = TRUE, pretty = TRUE)
    e$truth <- list(kind = "outcome", file = "sim_truth.json", goal = list(goal = "trend", field = "goal_truth", terms = list("treatment"), test = "wald"))
  }
  bioana_write_manifest(d, e)
  list(dir = d, grid = grid, x = x)
}

test_that("trend: level labels mapped to values, log10(x + log_offset) from the spec, size may be absent", {
  tf <- trend_folder(factor_extra = list(log_base = 10, log_offset = 1), truth_x = log10(c(0, 1, 10, 100) + 1))
  utils::capture.output(r <- bioana_analyze(tf$dir, secondary = FALSE))
  g <- r$summary$analyses$outcome$goal
  expect_match(g$scale_note, "log10 unit of treatment \\+ 1 mg \\(offset: factors\\[\\]\\.log_offset\\)")
  expect_null(g$effect$size)
  co <- stats::coef(summary(stats::lm(tf$grid$primary_outcome ~ log10(tf$x + 1))))[2, ]
  res <- r$results$outcome
  expect_equal(res$estimate[res$term == "treatment"], unname(co["Estimate"]), tolerance = 1e-8)
  gr <- r$summary$analyses$outcome$goal_recovery
  expect_true(gr$scored); expect_true(gr$terms_match); expect_true(gr$predictor_match)
  expect_true(gr$detected); expect_true(gr$covered); expect_equal(gr$true_effect, 0.8)
  expect_null(r$summary$analyses$outcome$recovery)   # no per-level truth under a numeric term
})

test_that("trend without log_base / log_offset uses goal_rules.json; a different planted predictor is flagged", {
  # Simitall fa37c0c: zero one decade below the smallest positive value, not log10(x + 1)
  tf <- trend_folder(effect = list(size = 0.5, scale = "outcome", definition = "slope per log10 unit"),
                     truth_x = c(-1, 0, 1, 2))
  utils::capture.output(r <- bioana_analyze(tf$dir, secondary = FALSE))
  g <- r$summary$analyses$outcome$goal
  expect_match(g$scale_note, "log10 unit of treatment \\+ 1 mg \\(offset: the smallest positive value")
  expect_equal(g$effect$size, 0.5); expect_equal(g$effect_scale_used, "outcome units")
  gr <- r$summary$analyses$outcome$goal_recovery
  expect_false(gr$predictor_match); expect_null(gr$covered); expect_match(gr$note, "not comparable")
})

within_folder <- function(seed = 11, G = 150) {
  d <- file.path(tempfile("within_"), "design"); dir.create(d, recursive = TRUE)
  grid <- expand.grid(rep = 1:4, time = c("0", "6", "24"), genotype = c("wt", "mutant"), stringsAsFactors = FALSE)
  grid$sample_id <- sprintf("S%02d", seq_len(nrow(grid)))
  mut <- grid$genotype == "mutant"; tm <- match(grid$time, c("0", "6", "24")) - 1
  set.seed(seed); grid$primary_outcome <- 10 + 0.3 * mut + 1.0 * mut * tm + stats::rnorm(nrow(grid), 0, 0.5)
  spec <- list(handoff_version = "0.6.2", origin = "observed", organism = "test", system = "lab",
               factors = list(list(name = "genotype", n_levels = 2L, levels = list("wt", "mutant")),
                              list(name = "time", n_levels = 3L, levels = list("0", "6", "24"), values = list(0, 6, 24), unit = "h")),
               design = list(chosen = "factorial", reps = 4L), primary_outcome = list(name = "growth", type = "continuous"),
               analysis = list(model = "y ~ treatment * time", terms = list(treatment = "genotype", time = "time", layout_columns = list("genotype", "time")),
                               goals = list("time_course"), primary_goal = "time_course", goal_terms = list(time_course = list("time")),
                               goal_effect = list(time_course = list(size = 2, scale = "outcome", shape = "monotone", scope = "within_group",
                                                                     group = "mutant", group_factor = "genotype"))))
  jsonlite::write_json(spec, file.path(d, "study_spec.json"), auto_unbox = TRUE, pretty = TRUE)
  utils::write.csv(grid[, c("sample_id", "genotype", "time", "rep")], file.path(d, "field_layout.csv"), row.names = FALSE)
  utils::write.table(grid, file.path(d, "sample_metadata_template.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
  mu <- 2^stats::rnorm(G, 7, 1); ch <- c(rep(1.5, 20), rep(0, G - 20))
  cnt <- sapply(seq_len(nrow(grid)), function(j) stats::rnbinom(G, mu = mu * 2^(ch * mut[j] * (tm[j] > 0)), size = 20))
  dimnames(cnt) <- list(sprintf("g%03d", 1:G), grid$sample_id)
  write_mat(cnt, file.path(d, "counts.tsv"))
  dir.create(file.path(d, "truth"))
  utils::write.table(data.frame(feature_id = sprintf("g%03d", 1:20), term = "treatment:time", effect = 1.5, change = 1.5),
                     file.path(d, "truth", "rna_goal.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
  utils::write.table(data.frame(feature_id = sprintf("g%03d", 21:30), condition = "mutant", effect = 0),
                     file.path(d, "truth", "rna.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
  gt <- list(goal = "time_course", model_terms = list("treatment:time"), status = "planted", scope = "within_group",
             groups_changing = list("mutant"), group_reference = "wt", true_effect = 2)
  jsonlite::write_json(list(intercept = 10, goal_truth = gt), file.path(d, "sim_truth.json"), auto_unbox = TRUE, pretty = TRUE)
  bioana_write_manifest(d,
    list(id = "outcome", data_type = "outcome", files = list(table = "sample_metadata_template.tsv"), column = "primary_outcome",
         truth = list(kind = "outcome", file = "sim_truth.json", goal = list(goal = "time_course", field = "goal_truth", terms = list("treatment:time"), test = "lrt"))),
    list(id = "rna", data_type = "feature_counts", assay = "rnaseq", files = list(counts = "counts.tsv"),
         truth = list(kind = "features", file = "truth/rna.tsv", term = "treatment", large = 1,
                      goal = list(goal = "time_course", file = "truth/rna_goal.tsv", term = "treatment:time", effect_scale = "log2 largest change"))))
  list(dir = d, grid = grid, counts = cnt)
}

test_that("time_course within one group = anova of the model without that group's time effects", {
  w <- within_folder()
  utils::capture.output(r <- bioana_analyze(w$dir, only = "outcome", secondary = FALSE))
  a <- r$summary$analyses$outcome
  expect_equal(a$goal$within$level, "mutant")
  expect_equal(unlist(a$goal$coefficients), c("time6", "time24"))
  g2 <- factor(w$grid$genotype, levels = c("mutant", "wt")); tt <- factor(w$grid$time, levels = c("0", "6", "24"))
  X <- stats::model.matrix(~ g2 * tt); X0 <- X[, !colnames(X) %in% c("tt6", "tt24")]
  ref <- stats::anova(stats::lm(w$grid$primary_outcome ~ X0 - 1), stats::lm(w$grid$primary_outcome ~ X - 1))
  expect_equal(a$goal$result$p, ref[2, "Pr(>F)"], tolerance = 1e-10)
  expect_equal(a$goal$result$statistic, ref[2, "F"], tolerance = 1e-10)
  gr <- a$goal_recovery
  expect_true(gr$scored); expect_true(gr$terms_match); expect_true(gr$group_match); expect_true(gr$detected)
})

test_that("goal recovery for count features: share of Simitall's goal genes found by the within-group test", {
  skip_if_not_installed("DESeq2")
  w <- within_folder()
  utils::capture.output(r <- bioana_analyze(w$dir, only = "rna", secondary = FALSE))
  a <- r$summary$analyses$rna
  expect_match(a$goal$label, "time_course within treatment = mutant")
  gr <- a$goal_recovery
  expect_true(gr$terms_match); expect_equal(gr$n_goal_features, 20)
  expect_gt(gr$tpr, 0.8); expect_lt(gr$called_outside_goal, 10)
  expect_equal(gr$called_other_planted + gr$called_not_planted, gr$called_outside_goal)
  # the reference level is mutant, so DESeq2's LRT drops the time columns of mutant
  cd <- data.frame(treatment = factor(w$grid$genotype, levels = c("mutant", "wt")), time = factor(w$grid$time, levels = c("0", "6", "24")))
  X <- stats::model.matrix(~ treatment * time, cd)
  dds <- DESeq2::DESeqDataSetFromMatrix(w$counts, cd, ~ treatment * time)
  dds <- suppressMessages(DESeq2::DESeq(dds, test = "LRT", full = X, reduced = X[, !colnames(X) %in% c("time6", "time24")], quiet = TRUE))
  g <- DESeq2::results(dds)
  t <- r$results$rna[match(rownames(g), r$results$rna$feature_id), ]
  expect_equal(t$p, g$pvalue, tolerance = 1e-6)
})

test_that("trend per_unit: linear factor stays linear (levels that double), per_unit reported; a log per_unit on a linear factor is flagged", {
  tf <- trend_folder(factor_extra = list(scale = "linear"), effect = list(size = 0.01, scale = "outcome", per_unit = "mg"))
  utils::capture.output(r <- bioana_analyze(tf$dir, secondary = FALSE))
  g <- r$summary$analyses$outcome$goal
  expect_equal(g$predictor$scale, "linear"); expect_match(g$scale_note, "linear scale.*size per mg")
  res <- r$results$outcome
  expect_equal(res$estimate[res$term == "treatment"], unname(stats::coef(stats::lm(tf$grid$primary_outcome ~ tf$x))[2]), tolerance = 1e-8)
  tf2 <- trend_folder(factor_extra = list(scale = "linear"), effect = list(size = 0.5, scale = "outcome", per_unit = "log10 mg"))
  utils::capture.output(r2 <- bioana_analyze(tf2$dir, secondary = FALSE))
  expect_match(r2$summary$analyses$outcome$goal$note, "does not match factor scale 'linear'")
})

test_that("within_group: Biodes' exact level label switches the reference; without group the test is over all groups", {
  w <- within_folder()
  s <- jsonlite::read_json(file.path(w$dir, "study_spec.json"))
  s$factors[[1]]$levels <- list("wild type", "mutant")
  jsonlite::write_json(s, file.path(w$dir, "study_spec.json"), auto_unbox = TRUE, pretty = TRUE)
  for (f in c("field_layout.csv", "sample_metadata_template.tsv")) {
    x <- readLines(file.path(w$dir, f)); writeLines(gsub("\\bwt\\b", "wild type", x), file.path(w$dir, f))
  }
  utils::capture.output(r <- bioana_analyze(w$dir, only = "outcome", secondary = FALSE))
  expect_equal(r$summary$analyses$outcome$goal$within$level, "mutant")
  # Biodes 0.6.3: no level matched -> group and group_factor dropped, scope kept
  s$analysis$goal_effect$time_course$group <- NULL; s$analysis$goal_effect$time_course$group_factor <- NULL
  jsonlite::write_json(s, file.path(w$dir, "study_spec.json"), auto_unbox = TRUE, pretty = TRUE)
  utils::capture.output(r <- bioana_analyze(w$dir, only = "outcome", secondary = FALSE))
  g <- r$summary$analyses$outcome$goal
  expect_null(g$within); expect_match(g$note, "names no group"); expect_equal(length(g$coefficients), 4)
})

test_that("a model with no fixed terms (p-rep catalog: y ~ (1 | entry)) runs", {
  skip_if_not_installed("lme4")
  d <- file.path(tempfile("prep_"), "design"); dir.create(d, recursive = TRUE)
  grid <- data.frame(sample_id = sprintf("S%02d", 1:30), entry = sprintf("E%02d", rep(1:15, 2)), stringsAsFactors = FALSE)
  set.seed(2); grid$primary_outcome <- 5 + stats::rnorm(15)[match(grid$entry, sprintf("E%02d", 1:15))] + stats::rnorm(30, 0, 0.5)
  spec <- list(handoff_version = "0.6.3", origin = "observed", organism = "test", system = "field",
               factors = list(list(name = "entry", n_levels = 15L)), design = list(chosen = "prep", reps = 2L),
               primary_outcome = list(name = "yield", type = "continuous"),
               analysis = list(model = "y ~ (1 | entry)", terms = list(entry = "entry", layout_columns = list("entry"))))
  jsonlite::write_json(spec, file.path(d, "study_spec.json"), auto_unbox = TRUE, pretty = TRUE)
  utils::write.csv(grid[, c("sample_id", "entry")], file.path(d, "field_layout.csv"), row.names = FALSE)
  utils::write.table(grid, file.path(d, "sample_metadata_template.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
  bioana_write_manifest(d, list(id = "outcome", data_type = "outcome", files = list(table = "sample_metadata_template.tsv"), column = "primary_outcome"))
  expect_no_error(bioana_check_design(d))
  utils::capture.output(r <- bioana_analyze(d, secondary = FALSE))
  expect_equal(r$summary$analyses$outcome$engine$engine, "lmer")
})
