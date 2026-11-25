# Safeguards: no silent substitution of a method that does not fit the design, and design checks before fitting.

counts_for <- function(meta, G = 120, seed = 1) {
  set.seed(seed)
  m <- sapply(seq_len(nrow(meta)), function(j) stats::rnbinom(G, mu = 2^stats::rnorm(G, 6, 1), size = 10))
  dimnames(m) <- list(sprintf("g%03d", 1:G), meta$sample_id); m
}
# sample metadata and field layout together (model terms are looked up in the layout)
write_both <- function(m, d) {
  write_meta(m, d)
  utils::write.csv(m[, setdiff(names(m), c("primary_outcome", "event")), drop = FALSE], file.path(d, "field_layout.csv"), row.names = FALSE)
}
with_model <- function(d, model, terms = NULL) {
  s <- jsonlite::read_json(file.path(d, "study_spec.json"))
  s$analysis$model <- model
  if (!is.null(terms)) s$analysis$terms <- utils::modifyList(s$analysis$terms, terms)
  jsonlite::write_json(s, file.path(d, "study_spec.json"), auto_unbox = TRUE)
}

test_that("a failed engine is replaced only by one that fits the design fully", {
  run <- function(eng) if (eng == "dream") stop("boom") else "ok"
  ch <- list(chain = c("dream", "limma_voom"), fit = c(dream = "full", limma_voom = "approx"))
  expect_error(bioana:::.ba_run_chain(ch, run, "feature_counts"), "dream failed \\(boom\\).*only approximates this design")
  ch$fit[["limma_voom"]] <- "full"
  r <- bioana:::.ba_run_chain(ch, run, "feature_counts")
  expect_equal(r$engine, "limma_voom")
  expect_equal(names(r$failed), "dream")
  w <- bioana:::.ba_run_chain(list(chain = "x", fit = c(x = "full")), function(eng) { warning("odd fit"); 1 }, "outcome")
  expect_equal(w$fit_warnings, "odd fit")
})

test_that("engines that cannot hold the design's random terms are skipped, never used", {
  d <- rcbd_folder()
  with_model(d, "y ~ treatment + (1 | block) + (1 | rep)", list(rep = "block"))
  design <- bioana:::.ba_read_biodes(d)
  ch <- tryCatch(bioana:::.ba_engine_chain("feature_counts", design, c("counts", "values")), error = function(e) e)
  if (inherits(ch, "error")) {
    expect_match(conditionMessage(ch), "limma_voom: models one random intercept")
  } else {
    expect_false("limma_voom" %in% ch$chain)          # two random terms: only dream
  }
})

test_that("a treatment confounded with another term stops the analysis", {
  d <- rcbd_folder(); m <- read_meta(d)
  m$block <- ifelse(m$treatment == "treatment_1", "1", "2")       # every control in block 1, every treated in block 2
  write_both(m, d)
  write_mat(counts_for(m), file.path(d, "counts.tsv"))
  bioana_write_manifest(d, list(id = "rna", data_type = "feature_counts", files = list(counts = "counts.tsv")))
  utils::capture.output(r <- bioana_analyze(d, secondary = FALSE))
  expect_match(r$summary$not_analysed$rna$note, "confounded")
  ck <- bioana_check_design(d)
  expect_true(any(ck$check == "estimability" & ck$status == "stop"))
})

test_that("a subject measured repeatedly without a random term stops; with the random term it runs", {
  skip_if_not_installed("limma"); skip_if_not_installed("edgeR")
  d <- rcbd_folder(); m <- read_meta(d)
  m$plant_id <- paste0("P", m$block)                               # both treatments measured on the same plant
  write_both(m, d)
  write_mat(counts_for(m), file.path(d, "counts.tsv"))
  bioana_write_manifest(d, list(id = "rna", data_type = "feature_counts", files = list(counts = "counts.tsv")))
  with_model(d, "y ~ treatment")
  utils::capture.output(r <- bioana_analyze(d, secondary = FALSE))
  expect_match(r$summary$not_analysed$rna$note, "repeated measures of one unit")
  with_model(d, "y ~ treatment + (1 | subject)", list(subject = "plant_id"))
  utils::capture.output(r <- bioana_analyze(d, secondary = FALSE))
  e <- r$summary$analyses$rna$engine
  expect_true(e$engine %in% c("dream", "limma_voom"))
  if (e$engine == "limma_voom") {
    expect_equal(e$design_fit, "approx")
    expect_match(e$random_terms, "duplicateCorrelation")
  }
})

test_that("technical replicates are not treated as independent units", {
  d <- rcbd_folder(); m <- read_meta(d)
  m$plant_id <- paste0("P", m$treatment, "_", ceiling(as.integer(m$block) / 2))   # two libraries per plant, same treatment
  write_both(m, d)
  write_mat(counts_for(m), file.path(d, "counts.tsv"))
  bioana_write_manifest(d, list(id = "rna", data_type = "feature_counts", files = list(counts = "counts.tsv")))
  with_model(d, "y ~ treatment")
  utils::capture.output(r <- bioana_analyze(d, secondary = FALSE))
  expect_match(r$summary$not_analysed$rna$note, "technical replicates")
})

test_that("normalised values are refused as counts; outcome values must match the type", {
  d <- rcbd_folder(); m <- read_meta(d)
  write_mat(counts_for(m) / 3.7, file.path(d, "counts.tsv"))
  bioana_write_manifest(d, list(id = "rna", data_type = "feature_counts", files = list(counts = "counts.tsv")))
  utils::capture.output(r <- bioana_analyze(d, secondary = FALSE))
  expect_match(r$summary$not_analysed$rna$note, "not normalised values")
  m$primary_outcome <- c(0, 1, 2, rep(c(0, 1), length.out = nrow(m) - 3))
  write_meta(m, d)
  s <- jsonlite::read_json(file.path(d, "study_spec.json")); s$primary_outcome$type <- "binary"
  jsonlite::write_json(s, file.path(d, "study_spec.json"), auto_unbox = TRUE)
  bioana_write_manifest(d, list(id = "outcome", data_type = "outcome", files = list(table = "sample_metadata_template.tsv"), column = "primary_outcome"))
  utils::capture.output(r <- bioana_analyze(d, secondary = FALSE))
  expect_match(r$summary$not_analysed$outcome$note, "binary outcome must be 0/1")
})

test_that("survival with a random term: coxme, or a gamma frailty; the Cox model never silently drops it", {
  skip_if_not_installed("survival")
  d <- rcbd_folder(); m <- read_meta(d)
  set.seed(3)
  m$primary_outcome <- round(stats::rexp(nrow(m), ifelse(m$treatment == "treatment_2", 0.5, 1)), 3) + 0.01
  m$event <- 1
  write_meta(m, d)
  s <- jsonlite::read_json(file.path(d, "study_spec.json")); s$primary_outcome$type <- "survival"
  s$analysis$model <- "y ~ treatment + (1 | block)"
  jsonlite::write_json(s, file.path(d, "study_spec.json"), auto_unbox = TRUE)
  bioana_write_manifest(d, list(id = "outcome", data_type = "outcome", files = list(table = "sample_metadata_template.tsv"), column = "primary_outcome"))
  utils::capture.output(r <- bioana_analyze(d, secondary = FALSE))
  e <- r$summary$analyses$outcome$engine
  expect_true(e$engine %in% c("coxme", "coxph_frailty"))
  if (e$engine == "coxph_frailty") expect_equal(e$design_fit, "approx")
  expect_true("treatmenttreatment_2" %in% r$results$outcome$term)
  expect_false(any(grepl("frailty", r$results$outcome$term)))
})

test_that("bioana_check_design reports the checks; a random term with one observation per level stops", {
  d <- rcbd_folder()
  ck <- bioana_check_design(d)
  expect_true(all(c("estimability", "residual_df", "replication") %in% ck$check))
  expect_false(any(ck$status == "stop"))
  m <- read_meta(d); m$plant_id <- m$sample_id; write_both(m, d)
  with_model(d, "y ~ treatment + (1 | subject)", list(subject = "plant_id"))
  ck <- bioana_check_design(d)
  expect_true(any(ck$check == "random_terms" & ck$status == "stop"))
})
