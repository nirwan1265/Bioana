# One engine for every feature-level data type (genes, miRNAs, taxa, peaks,
# methylation windows, m6A windows): a linear model per feature with the fixed
# terms of the design's analysis.model, and empirical-Bayes moderated t-tests
# (Smyth 2004, as limma eBayes) on every treatment coefficient. The data type's
# method only decides what goes in (log-CPM, M-values, log IP/input, ...).

# Y: features x samples on the scale to test (column names = sample ids).
.ba_feature_test <- function(Y, meta, p, alpha = 0.05, feature_label = "feature") {
  ids <- intersect(meta$sample_id, colnames(Y))
  if (length(ids) < 3) stop(sprintf("Only %d samples of the data match sample_metadata sample_id.", length(ids)), call. = FALSE)
  frame <- .ba_frame(meta[match(ids, meta$sample_id), , drop = FALSE], p)
  rhs <- stats::reformulate(p$fixed_labels)
  ok_rows <- stats::complete.cases(frame[, intersect(unlist(lapply(p$fixed_labels, function(l) all.vars(str2lang(l)))), names(frame)), drop = FALSE])
  frame <- frame[ok_rows, , drop = FALSE]; ids <- ids[ok_rows]
  X <- stats::model.matrix(rhs, data = frame)
  Y <- Y[, ids, drop = FALSE]
  Y <- Y[stats::complete.cases(Y) & apply(Y, 1, function(z) stats::var(z) > 0), , drop = FALSE]
  qrx <- qr(X)
  df <- nrow(X) - qrx$rank
  rows <- which(vapply(colnames(X), function(t) .ba_coef_role(t, p), "") == "treatment")
  rows <- rows[rows %in% qrx$pivot[seq_len(qrx$rank)]]
  empty <- data.frame(feature_id = character(), contrast = character(), estimate = numeric(), se = numeric(),
                      statistic = numeric(), p = numeric(), padj = numeric())
  if (!length(rows) || df < 1 || !nrow(Y)) {
    why <- if (!length(rows)) "the model has no fixed treatment term" else if (df < 1) "no residual degrees of freedom" else "no feature varies"
    return(list(table = empty, warnings = paste("features not tested:", why),
                summary = list(method = paste("not run:", why), design = paste(deparse(rhs), collapse = ""))))
  }
  Yt <- t(Y)
  beta <- qr.coef(qrx, Yt)
  resid <- Yt - X %*% beta
  s2 <- colSums(resid^2) / df
  # (X'X)^-1 of the estimable columns, in qr pivot order
  rk <- seq_len(qrx$rank)
  xtxi <- chol2inv(qr.R(qrx)[rk, rk, drop = FALSE])
  piv <- qrx$pivot[rk]
  sq <- .ba_squeeze_var(s2, df)
  tab <- do.call(rbind, lapply(rows, function(j) {
    k <- match(j, piv)
    se <- sqrt(sq$var_post * xtxi[k, k])
    tt <- beta[j, ] / se
    data.frame(feature_id = colnames(Yt), contrast = colnames(X)[j], estimate = beta[j, ], se = se, statistic = tt,
               p = 2 * stats::pt(-abs(tt), df + sq$df_prior), stringsAsFactors = FALSE, row.names = NULL)
  }))
  tab$padj <- stats::ave(tab$p, tab$contrast, FUN = function(x) stats::p.adjust(x, "BH"))
  list(table = tab, warnings = character(), summary = list(
    method = "per-feature linear model with empirical-Bayes moderated t-tests (as limma eBayes)",
    design = paste(deparse(rhs), collapse = ""), residual_df = df, prior_df = round(sq$df_prior, 2),
    random_terms_dropped = vapply(p$random, `[[`, "", "label"),
    n_features_tested = ncol(Yt), n_samples = nrow(X), contrasts = colnames(X)[rows],
    significant = as.list(tapply(tab$padj < alpha, tab$contrast, sum, na.rm = TRUE))))
}

# Empirical Bayes variance moderation (Smyth 2004, limma squeezeVar/fitFDist): a scaled
# F prior on the feature variances, fitted by moments of log(s2); posterior = weighted mean.
.ba_squeeze_var <- function(s2, df) {
  ok <- is.finite(s2) & s2 > 0
  if (df <= 0 || sum(ok) < 3) return(list(var_post = s2, df_prior = 0, var_prior = NA))
  z <- log(s2[ok])
  e <- z - digamma(df / 2) + log(df / 2)
  emean <- mean(e)
  evar <- stats::var(e) - trigamma(df / 2)
  if (evar > 0) {
    d0 <- 2 * .ba_trigamma_inverse(evar)
    s0 <- exp(emean + digamma(d0 / 2) - log(d0 / 2))
  } else {
    d0 <- Inf; s0 <- exp(emean)
  }
  post <- if (is.finite(d0)) (d0 * s0 + df * s2) / (d0 + df) else rep(s0, length(s2))
  post[!is.finite(post)] <- s0
  list(var_post = post, df_prior = if (is.finite(d0)) d0 else 1e6, var_prior = s0)
}

# Inverse of trigamma by Newton's method (as limma::trigammaInverse).
.ba_trigamma_inverse <- function(x) {
  if (x > 1e7) return(1 / sqrt(x))
  if (x < 1e-6) return(1 / x)
  y <- 0.5 + 1 / x
  for (i in 1:50) {
    tri <- trigamma(y)
    dif <- tri * (1 - tri / x) / psigamma(y, 2)
    y <- y + dif
    if (-dif / y < 1e-8) break
  }
  y
}

# Which model variable does a coefficient name belong to?
.ba_coef_role <- function(term, p) {
  if (term == "(Intercept)") return("intercept")
  parts <- strsplit(term, ":", fixed = TRUE)[[1]]
  roles <- vapply(parts, function(x) {
    v <- if (startsWith(x, "I(")) intersect(all.vars(str2lang(x)), p$vars) else p$vars[startsWith(x, p$vars)]
    if (length(v)) p$role[[v[which.max(nchar(v))]]] else "other"
  }, "")
  if (all(roles == "treatment")) "treatment" else if (any(roles == "block")) "block" else roles[1]
}

# Score feature calls against a truth table (feature_id, contrast, effect): per contrast, how many
# true changes were found (tpr), how many calls were wrong (fdr), and how well the estimates track
# the truth. `large` = |effect| from which a change counts as large (tpr_large).
.ba_feature_recovery <- function(tab, truth, alpha = 0.05, large = NULL) {
  out <- list()
  for (ct in unique(truth$contrast)) {
    t <- tab[tab$contrast == ct, , drop = FALSE]
    tr <- truth[truth$contrast == ct & truth$effect != 0, , drop = FALSE]
    called <- t$feature_id[!is.na(t$padj) & t$padj < alpha]
    tested <- intersect(tr$feature_id, t$feature_id)
    est <- t$estimate[match(tr$feature_id, t$feature_id)]
    big <- if (!is.null(large)) abs(tr$effect) >= large else rep(FALSE, nrow(tr))
    out[[ct]] <- list(
      contrast = ct, in_contrast = nrow(t) > 0, true_changes = nrow(tr), tested = length(tested), called = length(called),
      true_positives = sum(tr$feature_id %in% called),
      tpr = round(sum(tr$feature_id %in% called) / max(nrow(tr), 1), 3),
      tpr_large = if (any(big)) round(mean(tr$feature_id[big] %in% called), 3) else NA,
      large_threshold = large,
      fdr = round(length(setdiff(called, tr$feature_id)) / max(length(called), 1), 3),
      sign_agreement = if (sum(!is.na(est))) round(mean(sign(est) == sign(tr$effect), na.rm = TRUE), 3) else NA,
      estimate_cor = if (sum(!is.na(est)) > 2 && stats::sd(tr$effect) > 0) round(stats::cor(est, tr$effect, use = "complete.obs"), 3) else NA)
  }
  out
}

# Features (with chr/start/end) overlapping truth regions (chr/start/end + contrast, effect) take the region's truth.
.ba_regions_to_features <- function(features, regions) {
  out <- list()
  for (i in seq_len(nrow(regions))) {
    hit <- features$feature_id[features$chr == regions$chr[i] & features$start <= regions$end[i] & features$end >= regions$start[i]]
    if (length(hit)) out[[length(out) + 1]] <- data.frame(feature_id = hit, contrast = regions$contrast[i], effect = regions$effect[i],
                                                          region_id = regions$region_id[i] %||% i, stringsAsFactors = FALSE)
  }
  if (!length(out)) return(data.frame(feature_id = character(), contrast = character(), effect = numeric(), region_id = character()))
  r <- do.call(rbind, out)
  r[!duplicated(r[c("feature_id", "contrast")]), , drop = FALSE]
}
