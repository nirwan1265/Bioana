test_that("the continuous outcome is fitted with the design's model and scored against the truth", {
  d <- rcbd_folder(); m <- read_meta(d)
  set.seed(1)
  trt <- m$treatment == "treatment_2"
  m$primary_outcome <- 10 + 2 * trt + as.numeric(m$block) * 0.3 + stats::rnorm(nrow(m), 0, 1)
  write_meta(m, d)
  jsonlite::write_json(list(intercept = 10, treatment_effects = list(treatment = list(treatment_1 = 0, treatment_2 = 2))),
                       file.path(d, "sim_truth.json"), auto_unbox = TRUE)
  out <- utils::capture.output(r <- bioana_analyze(d))
  a <- r$summary$analyses$outcome
  expect_equal(a$outcome_type, "continuous")
  expect_match(a$model_fitted, "treatment \\+ block")
  res <- r$results$outcome
  expect_lt(abs(res$estimate[res$term == "treatmenttreatment_2"] - 2), 1.5)
  expect_equal(a$recovery$n_terms, 1)
  expect_true(file.exists(file.path(d, "bioana", "outcome", "recovery.tsv")))
  expect_true(file.exists(file.path(d, "bioana", "bioana_summary.json")))
})

test_that("binary, count, ordinal and survival outcomes get their own families", {
  set.seed(2)
  for (ty in c("binary", "count", "ordinal", "survival")) {
    d <- rcbd_folder(); m <- read_meta(d); trt <- m$treatment == "treatment_2"
    s <- jsonlite::read_json(file.path(d, "study_spec.json")); s$primary_outcome$type <- ty
    jsonlite::write_json(s, file.path(d, "study_spec.json"), auto_unbox = TRUE)
    m$primary_outcome <- switch(ty, binary = rep(c(0, 1), 6)[order(trt)], count = stats::rpois(12, ifelse(trt, 30, 5)),
                                ordinal = ifelse(trt, sample(3:4, 12, TRUE), sample(1:2, 12, TRUE)),
                                survival = ifelse(trt, 50, 10) * stats::rexp(12))
    if (ty == "survival") m$event <- 1
    write_meta(m, d)
    utils::capture.output(r <- suppressWarnings(bioana_analyze(d)))
    a <- r$summary$analyses$outcome
    expect_equal(a$outcome_type, ty)
    est <- r$results$outcome$estimate[r$results$outcome$term == "treatmenttreatment_2"]
    expect_true(if (ty == "survival") est < 0 else est > 0, info = ty)   # survival: longer times, lower hazard
  }
})

test_that("feature counts: planted changes are found, the rest mostly not", {
  d <- rcbd_folder(); m <- read_meta(d); trt <- m$treatment == "treatment_2"
  set.seed(3)
  G <- 400; mu <- 2^stats::rnorm(G, 7, 1.5)
  lfc <- c(rep(c(2, -2), 20), rep(0, G - 40))
  cnt <- sapply(seq_len(12), function(j) stats::rnbinom(G, mu = mu * 2^(lfc * trt[j]), size = 10))
  dimnames(cnt) <- list(sprintf("g%03d", 1:G), m$sample_id)
  write_mat(cnt, file.path(d, "counts.tsv"))
  dir.create(file.path(d, "truth"))
  utils::write.table(data.frame(feature_id = sprintf("g%03d", 1:40), condition = "treatment_2", effect = lfc[1:40]),
                     file.path(d, "truth", "rna.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
  bioana_write_manifest(d, list(id = "rna", data_type = "feature_counts", assay = "rnaseq", files = list(counts = "counts.tsv"),
                                truth = list(kind = "features", file = "truth/rna.tsv", term = "treatment", large = 1)))
  utils::capture.output(r <- bioana_analyze(d))
  rec <- r$summary$analyses$rna$recovery$treatmenttreatment_2
  expect_gt(rec$tpr, 0.8)
  expect_lt(rec$fdr, 0.15)
  expect_gt(rec$estimate_cor, 0.9)
})

test_that("bisulfite counts: windows over a planted DMR change, scored by overlap", {
  d <- rcbd_folder(); m <- read_meta(d); trt <- m$treatment == "treatment_2"
  set.seed(4)
  pos <- seq(1, 200000, by = 50); S <- length(pos)
  base <- stats::runif(S, 0.1, 0.9)
  dmr <- pos >= 50001 & pos <= 53000
  cov <- matrix(stats::rpois(S * 12, 20), S, 12)
  p <- sapply(seq_len(12), function(j) pmin(0.99, pmax(0.01, base + ifelse(dmr & trt[j], 0.3, 0) + stats::rnorm(S, 0, 0.03))))
  me <- matrix(stats::rbinom(S * 12, cov, p), S, 12)
  ids <- paste0("chr1:", pos)
  dimnames(cov) <- dimnames(me) <- list(ids, m$sample_id)
  write_mat(cov, file.path(d, "meth", "coverage.tsv")); write_mat(me, file.path(d, "meth", "methylated.tsv"))
  utils::write.table(data.frame(feature_id = ids, chr = "chr1", pos = pos, context = "CG"), file.path(d, "meth", "sites.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
  dir.create(file.path(d, "truth"))
  utils::write.table(data.frame(region_id = "dmr1", chr = "chr1", start = 50001, end = 53000, condition = "treatment_2", effect = 0.3),
                     file.path(d, "truth", "meth.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
  bioana_write_manifest(d, list(id = "meth", data_type = "bisulfite_counts", files = list(sites = "meth/sites.tsv", methylated = "meth/methylated.tsv", coverage = "meth/coverage.tsv"),
                                truth = list(kind = "regions", file = "truth/meth.tsv", term = "treatment")))
  utils::capture.output(r <- bioana_analyze(d))
  rec <- r$summary$analyses$meth$recovery$treatmenttreatment_2
  if (requireNamespace("DSS", quietly = TRUE)) expect_equal(r$summary$analyses$meth$engine$engine, "dss")
  expect_equal(rec$true_changes, 3)          # the DMR covers three 1 kb windows
  expect_equal(rec$true_positives, 3)
  expect_lt(rec$fdr, 0.3)
})

test_that("IP over input: an expression change alone is not an enrichment change", {
  d <- rcbd_folder(); m <- read_meta(d); trt <- m$treatment == "treatment_2"
  set.seed(5)
  W <- 300; expr <- 2^stats::rnorm(W, 5, 1); enr <- ifelse(seq_len(W) <= 150, 6, 1)
  expr_change <- ifelse(seq_len(W) %in% 151:200, 4, 1)     # background windows whose gene goes up 4-fold
  enr_change <- ifelse(seq_len(W) %in% 1:20, 3, 1)          # peaks whose enrichment triples
  mk <- function(f) sapply(seq_len(12), function(j) stats::rnbinom(W, mu = f(j), size = 20))
  inp <- mk(function(j) expr * if (trt[j]) expr_change else 1)
  ip <- mk(function(j) expr * (if (trt[j]) expr_change else 1) * enr * (if (trt[j]) enr_change else 1))
  dimnames(ip) <- dimnames(inp) <- list(sprintf("w%03d", 1:W), m$sample_id)
  write_mat(ip, file.path(d, "m6a", "ip.tsv")); write_mat(inp, file.path(d, "m6a", "input.tsv"))
  bioana_write_manifest(d, list(id = "m6a", data_type = "ip_input_counts", files = list(ip = "m6a/ip.tsv", input = "m6a/input.tsv")))
  utils::capture.output(r <- bioana_analyze(d))
  if (requireNamespace("edgeR", quietly = TRUE)) expect_equal(r$summary$analyses$m6a$engine$engine, "edger_ip_input")
  x <- r$results$m6a
  sig <- x$feature_id[x$padj < 0.05]
  expect_gt(mean(sprintf("w%03d", 1:20) %in% sig), 0.8)
  expect_lt(mean(sprintf("w%03d", 151:200) %in% sig), 0.1)
})

test_that("a folder without a manifest is inferred; unknown and planned data types are reported, not run", {
  d <- rcbd_folder(); m <- read_meta(d)
  m$primary_outcome <- seq_len(nrow(m)); write_meta(m, d)
  man <- bioana_manifest(d)
  expect_equal(man$data[[1]]$data_type, "outcome")
  expect_error(bioana_write_manifest(d, list(id = "x", data_type = "nope", files = list())), "Unknown data_type")
  expect_error(bioana_write_manifest(d, list(id = "x", data_type = "feature_counts", files = list())), "needs files: counts")
  jsonlite::write_json(list(manifest_version = "0.1", data = list(
    list(id = "outcome", data_type = "outcome", files = list(table = "sample_metadata_template.tsv"), column = "primary_outcome"),
    list(id = "cells", data_type = "single_cell_counts", files = list(counts = "c.tsv", cells = "cells.tsv")))), file.path(d, "data_manifest.json"), auto_unbox = TRUE)
  utils::capture.output(r <- bioana_analyze(d))
  expect_equal(r$summary$not_analysed$cells$status, "planned")
  expect_true("outcome" %in% names(r$summary$analyses))
})
