# Gold-standard checks: each engine adapter must give the same numbers as calling the package
# directly (the way its vignette does). Tutorials: docs/engines/<engine>.md.

engine_counts <- function(seed = 21, G = 300, lfc_n = 30) {
  d <- rcbd_folder(); m <- read_meta(d); trt <- m$treatment == "treatment_2"
  set.seed(seed)
  mu <- 2^stats::rnorm(G, 7, 1.5)
  lfc <- c(rep(c(1.5, -1.5), lfc_n / 2), rep(0, G - lfc_n))
  blk <- stats::rnorm(6, 0, 0.3)[as.integer(m$block)]
  cnt <- sapply(seq_len(nrow(m)), function(j) stats::rnbinom(G, mu = mu * 2^(lfc * trt[j] + blk[j]), size = 8))
  dimnames(cnt) <- list(sprintf("g%03d", 1:G), m$sample_id)
  design <- bioana:::.ba_read_biodes(d)
  list(dir = d, counts = cnt, design = design, frame = bioana:::.ba_frame(design$meta, design$parsed), lfc = lfc)
}

test_that("registry and selection rules are consistent", {
  e <- bioana_engines()
  s <- bioana:::.ba_tsv("engine_selection")
  expect_true(all(s$engine %in% e$engine))
  expect_true(all(c("package", "calls", "model", "test_statistic", "multiple_testing", "assumptions", "reference", "doc") %in% names(e)))
  root <- testthat::test_path("..", "..")
  if (dir.exists(file.path(root, "docs", "engines"))) expect_true(all(file.exists(file.path(root, unique(e$doc)))))
  plan <- bioana_engine_plan("feature_counts", "fixed")
  expect_equal(plan$engine[1], "deseq2")
  expect_equal(plan$role[plan$engine == "bioana_native"], "fallback")
})

test_that("DESeq2 adapter = DESeq2 called directly", {
  skip_if_not_installed("DESeq2")
  x <- engine_counts()
  r <- bioana:::.ba_engine_deseq2(list(counts = x$counts), x$design, 0.05, list())
  dds <- DESeq2::DESeqDataSetFromMatrix(x$counts, colData = x$frame, design = ~ treatment + block)
  dds <- suppressMessages(DESeq2::DESeq(dds, quiet = TRUE))
  g <- DESeq2::results(dds, name = "treatment_treatment_2_vs_treatment_1", alpha = 0.05)
  t <- r$table[match(rownames(g), r$table$feature_id), ]
  expect_equal(t$estimate, g$log2FoldChange, tolerance = 1e-6)
  expect_equal(t$se, g$lfcSE, tolerance = 1e-6)
  expect_equal(t$p, g$pvalue, tolerance = 1e-6)
  expect_equal(t$padj, g$padj, tolerance = 1e-6)
  expect_equal(unique(t$contrast), "treatmenttreatment_2")
})

test_that("edgeR adapter = edgeR QL pipeline called directly", {
  skip_if_not_installed("edgeR")
  x <- engine_counts()
  r <- bioana:::.ba_engine_edger(list(counts = x$counts), x$design, 0.05, list())
  X <- stats::model.matrix(~ treatment + block, x$frame)
  y <- bioana:::.ba_edger_norm(edgeR::DGEList(x$counts))
  fit <- edgeR::glmQLFit(edgeR::estimateDisp(y, X), X)
  g <- edgeR::topTags(edgeR::glmQLFTest(fit, coef = 2), n = Inf, sort.by = "none")$table
  t <- r$table[match(rownames(g), r$table$feature_id), ]
  expect_equal(t$estimate, g$logFC, tolerance = 1e-8)
  expect_equal(t$p, g$PValue, tolerance = 1e-8)
  expect_equal(t$padj, g$FDR, tolerance = 1e-8)
})

test_that("limma-voom adapter = voom/lmFit/eBayes called directly", {
  skip_if_not_installed("limma"); skip_if_not_installed("edgeR")
  x <- engine_counts()
  r <- bioana:::.ba_engine_limma_voom(list(counts = x$counts), x$design, 0.05, list())
  X <- stats::model.matrix(~ treatment + block, x$frame)
  v <- limma::voom(bioana:::.ba_edger_norm(edgeR::DGEList(x$counts)), X)
  fit <- limma::eBayes(limma::lmFit(v, X))
  g <- limma::topTable(fit, coef = "treatmenttreatment_2", number = Inf, sort.by = "none")
  t <- r$table[match(rownames(g), r$table$feature_id), ]
  expect_equal(t$estimate, g$logFC, tolerance = 1e-8)
  expect_equal(t$statistic, g$t, tolerance = 1e-8)
  expect_equal(t$padj, g$adj.P.Val, tolerance = 1e-8)
})

test_that("limma with one random term = duplicateCorrelation called directly", {
  skip_if_not_installed("limma"); skip_if_not_installed("edgeR")
  x <- engine_counts()
  s <- jsonlite::read_json(file.path(x$dir, "study_spec.json"))
  s$analysis$model <- "y ~ treatment + (1 | block)"
  jsonlite::write_json(s, file.path(x$dir, "study_spec.json"), auto_unbox = TRUE)
  design <- bioana:::.ba_read_biodes(x$dir)
  ch <- bioana:::.ba_engine_chain("feature_counts", design, c("counts", "values"))
  if (!requireNamespace("variancePartition", quietly = TRUE)) expect_equal(ch$chain[1], "limma_voom")
  r <- bioana:::.ba_engine_limma_voom(list(counts = x$counts), design, 0.05, list())
  X <- stats::model.matrix(~ treatment, x$frame); blk <- x$frame$block
  y <- bioana:::.ba_edger_norm(edgeR::DGEList(x$counts))
  v <- limma::voom(y, X)
  cr <- limma::duplicateCorrelation(v, X, block = blk)$consensus.correlation
  v <- limma::voom(y, X, block = blk, correlation = cr)
  cr <- limma::duplicateCorrelation(v, X, block = blk)$consensus.correlation
  g <- limma::topTable(limma::eBayes(limma::lmFit(v, X, block = blk, correlation = cr)), coef = 2, number = Inf, sort.by = "none")
  t <- r$table[match(rownames(g), r$table$feature_id), ]
  expect_equal(t$statistic, g$t, tolerance = 1e-8)
  expect_match(r$random, "duplicateCorrelation")
})

test_that("dream adapter = variancePartition::dream called directly", {
  skip_if_not_installed("variancePartition")
  x <- engine_counts()
  s <- jsonlite::read_json(file.path(x$dir, "study_spec.json"))
  s$analysis$model <- "y ~ treatment + (1 | block)"
  jsonlite::write_json(s, file.path(x$dir, "study_spec.json"), auto_unbox = TRUE)
  design <- bioana:::.ba_read_biodes(x$dir)
  r <- bioana:::.ba_engine_dream(list(counts = x$counts), design, 0.05, list())
  form <- ~ treatment + (1 | block)
  y <- bioana:::.ba_edger_norm(edgeR::DGEList(x$counts))
  v <- variancePartition::voomWithDreamWeights(y, form, x$frame, quiet = TRUE)
  fit <- variancePartition::eBayes(variancePartition::dream(v, form, x$frame, quiet = TRUE))
  g <- variancePartition::topTable(fit, coef = "treatmenttreatment_2", number = Inf, sort.by = "none")
  t <- r$table[match(rownames(g), r$table$feature_id), ]
  expect_equal(t$statistic, g$t, tolerance = 1e-6)
})

test_that("limma on values = lmFit/eBayes; the native fallback reproduces limma", {
  skip_if_not_installed("limma")
  x <- engine_counts()
  Y <- log2(x$counts + 1)
  r <- bioana:::.ba_engine_limma(list(values = Y), x$design, 0.05, list())
  X <- stats::model.matrix(~ treatment + block, x$frame)
  fit <- limma::eBayes(limma::lmFit(Y, X))
  g <- limma::topTable(fit, coef = "treatmenttreatment_2", number = Inf, sort.by = "none")
  t <- r$table[match(rownames(g), r$table$feature_id), ]
  expect_equal(t$statistic, g$t, tolerance = 1e-8)
  nat <- bioana:::.ba_engine_bioana_native(list(values = Y), x$design, 0.05, list())$table
  nat <- nat[match(rownames(g), nat$feature_id), ]
  expect_equal(nat$estimate, g$logFC, tolerance = 1e-8)
  expect_gt(stats::cor(nat$statistic, g$t), 0.999)
  expect_lt(max(abs(nat$statistic - g$t)), 0.05 * max(abs(g$t)))
})

test_that("DSS adapter = DMLfit.multiFactor / DMLtest.multiFactor called directly", {
  skip_if_not_installed("DSS")
  x <- engine_counts(); m <- x$design$meta; trt <- m$treatment == "treatment_2"
  set.seed(9); W <- 80
  C <- matrix(stats::rpois(W * 12, 40), W, 12)
  p <- outer(c(rep(0.3, 10), rep(0.5, W - 10)), ifelse(trt, 1, 0)) * 0.6 + 0.3
  M <- matrix(stats::rbinom(W * 12, C, pmin(p, 0.95)), W, 12)
  ids <- paste0("chr1:", seq_len(W) * 1000 + 1)
  dimnames(C) <- dimnames(M) <- list(ids, m$sample_id)
  feats <- data.frame(feature_id = ids, chr = "chr1", start = seq_len(W) * 1000 + 1)
  r <- bioana:::.ba_engine_dss(list(methylated = M, coverage = C, features = feats), x$design, 0.05, list())
  bs <- bsseq::BSseq(chr = feats$chr, pos = feats$start, M = M, Cov = C, sampleNames = m$sample_id)
  utils::capture.output(fit <- DSS::DMLfit.multiFactor(bs, design = x$frame[, c("treatment", "block")], formula = ~ treatment + block))
  g <- DSS::DMLtest.multiFactor(fit, coef = "treatmenttreatment_2")
  expect_equal(r$table$statistic, g$stat, tolerance = 1e-8)
  expect_equal(r$table$padj, g$fdrs, tolerance = 1e-8)
  expect_equal(r$table$feature_id, ids)
  expect_true(all(r$table$padj[1:10] < 0.05))
})

test_that("edgeR IP x input: enrichment changes are found, expression changes are not", {
  skip_if_not_installed("edgeR")
  x <- engine_counts(); m <- x$design$meta; trt <- m$treatment == "treatment_2"
  set.seed(4); W <- 200
  base <- 2^stats::rnorm(W, 6, 1)
  expr <- c(rep(0, 20), rep(1.5, 20), rep(0, W - 40))      # windows 21-40: expression up (IP and input)
  enr  <- c(rep(1.5, 20), rep(0, W - 20))                  # windows 1-20: enrichment up (IP only)
  inp <- sapply(seq_len(12), function(j) stats::rnbinom(W, mu = base * 2^(expr * trt[j]), size = 20))
  ip  <- sapply(seq_len(12), function(j) stats::rnbinom(W, mu = 3 * base * 2^((expr + enr) * trt[j]), size = 20))
  dimnames(inp) <- dimnames(ip) <- list(sprintf("w%03d", 1:W), m$sample_id)
  r <- bioana:::.ba_engine_edger_ip_input(list(ip = ip, input = inp), x$design, 0.05, list())
  sig <- r$table$feature_id[r$table$padj < 0.05]
  expect_gt(mean(sprintf("w%03d", 1:20) %in% sig), 0.8)
  expect_lt(mean(sprintf("w%03d", 21:40) %in% sig), 0.15)
})

test_that("bioana_analyze records provenance and runs the secondary engine on request", {
  skip_if_not_installed("DESeq2"); skip_if_not_installed("edgeR")
  x <- engine_counts()
  write_mat(x$counts, file.path(x$dir, "counts.tsv"))
  bioana_write_manifest(x$dir, list(id = "rna", data_type = "feature_counts", assay = "rnaseq", files = list(counts = "counts.tsv")))
  out <- utils::capture.output(r <- bioana_analyze(x$dir, secondary = FALSE))
  e <- r$summary$analyses$rna$engine
  expect_equal(e$package, "DESeq2")
  expect_equal(e$version, as.character(utils::packageVersion("DESeq2")))
  expect_false(e$fallback_used)
  expect_true(all(c("call", "model", "input", "assumptions", "test_statistic", "multiple_testing", "experimental_unit", "reference") %in% names(e)))
  expect_equal(r$summary$analyses$rna$secondary_engine, "edger")
  expect_true(any(grepl("\\[rna\\] \\[DESeq2", out)))
  expect_true(any(grepl("secondary engines not run", out)))
  out <- utils::capture.output(r <- bioana_analyze(x$dir, secondary = TRUE))
  cc <- r$summary$analyses$rna$secondary$concordance
  expect_gt(cc$estimate_cor, 0.95)
  expect_true(file.exists(file.path(x$dir, "bioana", "rna", "results_edger.tsv")))
  expect_true(any(grepl("secondary edgeR", out)))
  # engine override from params
  utils::capture.output(r <- bioana_analyze(x$dir, params = list(feature_counts = list(engine = "bioana_native")), secondary = FALSE))
  expect_true(r$summary$analyses$rna$engine$fallback_used)
})

test_that("glmmTMB outcome engine matches glmmTMB called directly, and agrees with lm", {
  skip_if_not_installed("glmmTMB")
  d <- rcbd_folder(); m <- read_meta(d)
  set.seed(1)
  m$primary_outcome <- 10 + 2 * (m$treatment == "treatment_2") + stats::rnorm(nrow(m))
  write_meta(m, d)
  utils::capture.output(r <- bioana_analyze(d, secondary = TRUE))
  a <- r$summary$analyses$outcome
  expect_equal(a$engine$engine, "stats_lm")
  expect_equal(a$secondary$engine$engine, "glmmtmb")
  expect_gt(a$secondary$concordance$estimate_cor, 0.999)
  design <- bioana:::.ba_read_biodes(d); fr <- bioana:::.ba_frame(design$meta, design$parsed)
  fr$y_ <- as.numeric(design$meta$primary_outcome)
  g <- summary(glmmTMB::glmmTMB(y_ ~ treatment + block, data = fr))$coefficients$cond
  t2 <- utils::read.delim(file.path(d, "bioana", "outcome", "results_glmmtmb.tsv"))
  expect_equal(t2$estimate[t2$term == "treatmenttreatment_2"], g["treatmenttreatment_2", "Estimate"], tolerance = 1e-6)
})
