# Feature engines: one adapter per package. Each takes the inputs a data type
# gives (counts, values, methylated + coverage, ip + input), the design, alpha
# and args, and returns the standard table (feature_id, contrast, estimate,
# se, statistic, p, padj; contrast = model-matrix column, as the truth files
# use) plus the exact call it made. Tutorials: docs/engines/<engine>.md.

.ba_int_counts <- function(m, ids) {
  cnt <- round(m[, ids, drop = FALSE])
  storage.mode(cnt) <- "integer"
  cnt
}

.ba_edger_norm <- function(y) {
  ns <- asNamespace("edgeR")
  if (exists("normLibSizes", envir = ns)) edgeR::normLibSizes(y) else edgeR::calcNormFactors(y)
}

.ba_engine_result <- function(tab, d, method, alpha, call, input, warnings = character(), random = NULL, diagnostics = NULL) {
  rownames(tab) <- NULL
  list(table = tab, summary = .ba_engine_summary(tab, d, method, alpha), call = call, input = input,
       n_units = length(d$ids), warnings = warnings, random = random, diagnostics = diagnostics)
}

.ba_no_random <- function(d, engine) {
  if (d$n_random) stop(sprintf("%s fits fixed terms only; the model has random terms (%s)", engine, paste(d$random_labels, collapse = ", ")), call. = FALSE)
}

# DESeq2: NB GLM, the design's fixed terms as a model matrix, Wald test per treatment coefficient.
.ba_engine_deseq2 <- function(inputs, design, alpha, args) {
  d <- .ba_feature_design(colnames(inputs$counts), design$meta, design$parsed)
  .ba_need_coefs(d); .ba_no_random(d, "DESeq2")
  cnt <- .ba_int_counts(inputs$counts, d$ids)
  X <- d$X
  key <- stats::setNames(colnames(X), make.names(colnames(X), unique = TRUE))
  colnames(X) <- names(key)
  dds <- DESeq2::DESeqDataSetFromMatrix(cnt, colData = d$frame, design = ~1)
  dds <- suppressMessages(DESeq2::DESeq(dds, full = X, quiet = TRUE))
  tab <- do.call(rbind, lapply(d$coefs, function(cf) {
    r <- DESeq2::results(dds, name = names(key)[key == cf], alpha = alpha)
    data.frame(feature_id = rownames(r), contrast = cf, estimate = r$log2FoldChange, se = r$lfcSE,
               statistic = r$stat, p = r$pvalue, padj = r$padj, stringsAsFactors = FALSE)
  }))
  nc <- sum(!S4Vectors::mcols(dds)$betaConv, na.rm = TRUE)
  .ba_engine_result(tab, d, "DESeq2 negative-binomial GLM, Wald test", alpha,
                    sprintf("DESeq2::DESeq(DESeqDataSetFromMatrix(counts, colData, design = ~1), full = model.matrix(%s)); DESeq2::results(dds, name = <coefficient>, alpha = %g)",
                            paste(deparse(d$rhs), collapse = ""), alpha),
                    sprintf("raw integer counts, %d features x %d units", nrow(cnt), ncol(cnt)),
                    warnings = if (nc) sprintf("DESeq2: %d features did not converge (betaConv FALSE)", nc) else character(),
                    diagnostics = list(non_converged = nc))
}

# edgeR: TMM, NB dispersions, quasi-likelihood F test per treatment coefficient.
.ba_engine_edger <- function(inputs, design, alpha, args) {
  d <- .ba_feature_design(colnames(inputs$counts), design$meta, design$parsed)
  .ba_need_coefs(d); .ba_no_random(d, "edgeR")
  cnt <- .ba_int_counts(inputs$counts, d$ids)
  y <- .ba_edger_norm(edgeR::DGEList(cnt))
  y <- edgeR::estimateDisp(y, d$X)
  fit <- edgeR::glmQLFit(y, d$X)
  tab <- do.call(rbind, lapply(d$coefs, function(cf) {
    q <- edgeR::glmQLFTest(fit, coef = match(cf, colnames(d$X)))
    tt <- edgeR::topTags(q, n = Inf, sort.by = "none")$table
    data.frame(feature_id = rownames(tt), contrast = cf, estimate = tt$logFC, se = NA_real_,
               statistic = tt$F, p = tt$PValue, padj = tt$FDR, stringsAsFactors = FALSE)
  }))
  .ba_engine_result(tab, d, "edgeR quasi-likelihood NB GLM, QL F test", alpha,
                    sprintf("y <- normLibSizes(DGEList(counts)); estimateDisp(y, X); fit <- glmQLFit(y, X); glmQLFTest(fit, coef = <coefficient>); topTags(n = Inf), X = model.matrix(%s)",
                            paste(deparse(d$rhs), collapse = "")),
                    sprintf("raw integer counts, %d features x %d units", nrow(cnt), ncol(cnt)))
}

# limma fit on an EList or matrix, with duplicateCorrelation when the design has one random grouping term.
.ba_limma_fit <- function(obj, d, counts_y = NULL) {
  if (d$n_random > 1) stop("limma models one random grouping term; the design has ", d$n_random, call. = FALSE)
  if (is.null(d$block)) {
    v <- if (!is.null(counts_y)) limma::voom(counts_y, d$X) else obj
    return(list(fit = limma::eBayes(limma::lmFit(v, d$X)), correlation = NULL))
  }
  v <- if (!is.null(counts_y)) limma::voom(counts_y, d$X) else obj
  cor <- limma::duplicateCorrelation(v, d$X, block = d$block)$consensus.correlation
  if (!is.null(counts_y)) {
    v <- limma::voom(counts_y, d$X, block = d$block, correlation = cor)
    cor <- limma::duplicateCorrelation(v, d$X, block = d$block)$consensus.correlation
  }
  list(fit = limma::eBayes(limma::lmFit(v, d$X, block = d$block, correlation = cor)), correlation = cor)
}

.ba_limma_table <- function(fit, coefs) {
  do.call(rbind, lapply(coefs, function(cf) {
    tt <- limma::topTable(fit, coef = cf, number = Inf, sort.by = "none")
    data.frame(feature_id = rownames(tt), contrast = cf, estimate = tt$logFC,
               se = (sqrt(fit$s2.post) * fit$stdev.unscaled[, cf])[rownames(tt)],
               statistic = tt$t, p = tt$P.Value, padj = tt$adj.P.Val, stringsAsFactors = FALSE)
  }))
}

# limma-voom on counts.
.ba_engine_limma_voom <- function(inputs, design, alpha, args) {
  d <- .ba_feature_design(colnames(inputs$counts), design$meta, design$parsed)
  .ba_need_coefs(d)
  cnt <- .ba_int_counts(inputs$counts, d$ids)
  y <- .ba_edger_norm(edgeR::DGEList(cnt))
  lf <- .ba_limma_fit(NULL, d, counts_y = y)
  rnd <- if (is.null(lf$correlation)) NULL else sprintf("%s as a within-group correlation (duplicateCorrelation, consensus %.3f)", d$random_labels, lf$correlation)
  .ba_engine_result(.ba_limma_table(lf$fit, d$coefs), d, "limma-voom precision-weighted linear model, moderated t", alpha,
                    sprintf("v <- voom(normLibSizes(DGEList(counts)), X%s); fit <- eBayes(lmFit(v, X%s)); topTable(fit, coef = <coefficient>, number = Inf), X = model.matrix(%s)",
                            if (is.null(rnd)) "" else ", block, correlation", if (is.null(rnd)) "" else ", block, correlation",
                            paste(deparse(d$rhs), collapse = "")),
                    sprintf("raw integer counts, %d features x %d units", nrow(cnt), ncol(cnt)), random = rnd)
}

# limma on values already on a linear scale (M-values, log ratios).
.ba_engine_limma <- function(inputs, design, alpha, args) {
  d <- .ba_feature_design(colnames(inputs$values), design$meta, design$parsed)
  .ba_need_coefs(d)
  Y <- inputs$values[, d$ids, drop = FALSE]
  Y <- Y[stats::complete.cases(Y) & apply(Y, 1, stats::var) > 0, , drop = FALSE]
  lf <- .ba_limma_fit(Y, d)
  rnd <- if (is.null(lf$correlation)) NULL else sprintf("%s as a within-group correlation (duplicateCorrelation, consensus %.3f)", d$random_labels, lf$correlation)
  .ba_engine_result(.ba_limma_table(lf$fit, d$coefs), d, "limma linear model, moderated t", alpha,
                    sprintf("fit <- eBayes(lmFit(Y, X%s)); topTable(fit, coef = <coefficient>, number = Inf), X = model.matrix(%s)",
                            if (is.null(rnd)) "" else ", block, correlation", paste(deparse(d$rhs), collapse = "")),
                    sprintf("%s, %d features x %d units", inputs$values_scale %||% "values", nrow(Y), ncol(Y)), random = rnd)
}

# dream: precision-weighted linear mixed model with the design's random terms.
.ba_engine_dream <- function(inputs, design, alpha, args) {
  d <- .ba_feature_design(colnames(inputs$counts), design$meta, design$parsed)
  .ba_need_coefs(d)
  cnt <- .ba_int_counts(inputs$counts, d$ids)
  y <- .ba_edger_norm(edgeR::DGEList(cnt))
  form <- d$full_rhs
  vobj <- variancePartition::voomWithDreamWeights(y, form, d$frame, quiet = TRUE)
  fit <- variancePartition::dream(vobj, form, d$frame, quiet = TRUE)
  fit <- variancePartition::eBayes(fit)
  tab <- do.call(rbind, lapply(d$coefs, function(cf) {
    tt <- variancePartition::topTable(fit, coef = cf, number = Inf, sort.by = "none")
    data.frame(feature_id = rownames(tt), contrast = cf, estimate = tt$logFC,
               se = (sqrt(fit$s2.post) * fit$stdev.unscaled[, cf])[rownames(tt)],
               statistic = tt$t, p = tt$P.Value, padj = tt$adj.P.Val, stringsAsFactors = FALSE)
  }))
  .ba_engine_result(tab, d, "dream precision-weighted linear mixed model, moderated t (Satterthwaite df)", alpha,
                    sprintf("vobj <- voomWithDreamWeights(normLibSizes(DGEList(counts)), %s, metadata); fit <- eBayes(dream(vobj, %s, metadata)); topTable(fit, coef = <coefficient>, number = Inf)",
                            paste(deparse(form), collapse = ""), paste(deparse(form), collapse = "")),
                    sprintf("raw integer counts, %d features x %d units", nrow(cnt), ncol(cnt)),
                    random = sprintf("%s fitted as random terms (lme4)", paste(d$random_labels, collapse = ", ")))
}

# DSS: beta-binomial regression on methylated and total reads per window, general design.
.ba_engine_dss <- function(inputs, design, alpha, args) {
  d <- .ba_feature_design(colnames(inputs$coverage), design$meta, design$parsed)
  .ba_need_coefs(d); .ba_no_random(d, "DSS")
  M <- inputs$methylated[, d$ids, drop = FALSE]; C <- inputs$coverage[rownames(M), d$ids, drop = FALSE]
  f <- inputs$features[match(rownames(M), inputs$features$feature_id), , drop = FALSE]
  # one DSS chromosome per chromosome x context, so windows of different contexts at one position stay apart
  chr <- if (!is.null(f$context)) paste(f$chr, f$context, sep = "|") else f$chr
  bs <- bsseq::BSseq(chr = chr, pos = f$start, M = M, Cov = C, sampleNames = d$ids)
  vars <- all.vars(d$rhs)
  utils::capture.output(fit <- DSS::DMLfit.multiFactor(bs, design = d$frame[, vars, drop = FALSE], formula = d$rhs))
  if (!identical(colnames(fit$X), colnames(stats::model.matrix(d$rhs, d$frame)))) {
    stop("DSS built a different model matrix from Bioana's", call. = FALSE)
  }
  k <- ncol(fit$X)
  key <- paste(as.character(GenomeInfoDb::seqnames(fit$gr)), BiocGenerics::start(fit$gr))
  fid <- rownames(M)[match(key, paste(chr, f$start))]
  tab <- do.call(rbind, lapply(d$coefs, function(cf) {
    j <- match(cf, colnames(fit$X))
    tt <- DSS::DMLtest.multiFactor(fit, coef = cf)
    data.frame(feature_id = fid, contrast = cf, estimate = fit$fit$beta[, j], se = sqrt(fit$fit$var.beta[, (j - 1) * k + j]),
               statistic = tt$stat, p = tt$pvals, padj = tt$fdrs, stringsAsFactors = FALSE)
  }))
  .ba_engine_result(tab, d, "DSS beta-binomial regression (arcsine link), Wald test", alpha,
                    sprintf("fit <- DMLfit.multiFactor(BSseq(chr, pos, M, Cov), design, formula = %s); DMLtest.multiFactor(fit, coef = <coefficient>)",
                            paste(deparse(d$rhs), collapse = "")),
                    sprintf("methylated and total reads pooled per window, %d windows x %d units", nrow(M), ncol(M)))
}

# edgeR on IP and input libraries together: unit + library type + library type x the design's fixed terms.
.ba_engine_edger_ip_input <- function(inputs, design, alpha, args) {
  d <- .ba_feature_design(colnames(inputs$ip), design$meta, design$parsed)
  .ba_need_coefs(d); .ba_no_random(d, "edgeR (IP x input)")
  ip <- .ba_int_counts(inputs$ip, d$ids); inp <- .ba_int_counts(inputs$input[rownames(ip), , drop = FALSE], d$ids)
  n <- length(d$ids)
  U <- diag(n); colnames(U) <- paste0("unit_", d$ids)
  Xf <- d$X[, colnames(d$X) != "(Intercept)", drop = FALSE]
  X2 <- rbind(cbind(U, ip = 1, Xf), cbind(U, ip = 0, Xf * 0))
  colnames(X2) <- c(colnames(U), "ip", paste0("ip:", colnames(Xf)))
  both <- cbind(ip, inp)
  colnames(both) <- c(paste0("ip_", d$ids), paste0("input_", d$ids))
  y <- edgeR::DGEList(both)
  # offsets as the native method: input by library size, IP by its median ratio to input (IP efficiency per unit)
  sf_in <- colSums(inp) / mean(colSums(inp))
  sf_ip <- apply((ip + 1) / sweep(inp + 1, 2, sf_in, "/"), 2, stats::median)
  y$samples$lib.size <- 1e6 * c(sf_ip, sf_in)
  y$samples$norm.factors <- 1
  y <- edgeR::estimateDisp(y, X2)
  fit <- edgeR::glmQLFit(y, X2)
  tab <- do.call(rbind, lapply(d$coefs, function(cf) {
    q <- edgeR::glmQLFTest(fit, coef = match(paste0("ip:", cf), colnames(X2)))
    tt <- edgeR::topTags(q, n = Inf, sort.by = "none")$table
    data.frame(feature_id = rownames(tt), contrast = cf, estimate = tt$logFC, se = NA_real_,
               statistic = tt$F, p = tt$PValue, padj = tt$FDR, stringsAsFactors = FALSE)
  }))
  .ba_engine_result(tab, d, "edgeR QL NB GLM on IP and input, IP x treatment interaction", alpha,
                    sprintf("y <- DGEList(cbind(ip, input)) with offsets (input: library size; IP: median IP/input ratio); fit <- glmQLFit(estimateDisp(y, X2), X2); glmQLFTest(fit, coef = 'ip:<coefficient>'), X2 = unit + ip + ip:(model.matrix(%s) without intercept)",
                            paste(deparse(d$rhs), collapse = "")),
                    sprintf("IP and input counts, %d windows x %d units", nrow(ip), n))
}

# Native fallback: per-feature linear model with limma-style variance moderation (Bioana code).
.ba_engine_bioana_native <- function(inputs, design, alpha, args) {
  ft <- .ba_feature_test(inputs$values, design$meta, design$parsed, alpha)
  if (!nrow(ft$table)) stop(ft$summary$method, call. = FALSE)
  ids <- intersect(design$meta$sample_id, colnames(inputs$values))
  list(table = ft$table, summary = ft$summary, warnings = ft$warnings,
       call = sprintf("bioana:::.ba_feature_test(Y, meta, parsed_model); fixed terms %s", paste(deparse(stats::reformulate(design$parsed$fixed_labels)), collapse = "")),
       input = sprintf("%s, %d features", inputs$values_scale %||% "values", nrow(inputs$values)), n_units = ft$summary$n_samples,
       random = if (length(design$parsed$random)) sprintf("dropped: %s", paste(vapply(design$parsed$random, `[[`, "", "label"), collapse = ", ")) else NULL)
}
