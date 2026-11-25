# Genotypes: GWAS on the line means of the primary outcome.
#
# 1. The genotype term: the manifest entry's genotype_term, else sim_truth.json
#    genetics.genotype_term (Simitall), else the treatment term that matches
#    inst/rules/genetics.json genotype_pattern.
# 2. Line means (BLUEs): the design's analysis.model fitted as for the outcome,
#    with the genotype term fixed (refitted as fixed if the model has it random);
#    each line's fixed-effect prediction averaged over the levels of the other
#    fixed factors (block, environment), as emmeans would. Random terms drop out.
# 3. analyze_gwas() on the VCF with the line means (rrBLUP or native EMMAX).
# 4. Truth (simulated): additive and dominance QTN count in tpr; a hit counts
#    when it is the QTN marker or in LD with it (ld_r2) or within window_bp.

.ba_method_genotypes <- function(entry, design, outdir, alpha = 0.05, args = list()) {
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  warn <- character()
  geno <- .ba_genotype_term(entry, design)
  vcf <- file.path(design$dir, entry$files$vcf %||% stop("genotypes entry has no files$vcf", call. = FALSE))
  lm_out <- .ba_line_means(design, geno, entry, args)
  warn <- c(warn, lm_out$warnings)
  means <- lm_out$means
  gt <- .bioana_read_genotypes(vcf)
  in_vcf <- means$line %in% colnames(gt$genotype)
  if (!any(in_vcf)) {
    stop(sprintf("No level of '%s' (layout column '%s') is a sample of %s: the VCF columns must be the line names.",
                 geno, design$parsed$col[[geno]], entry$files$vcf), call. = FALSE)
  }
  if (!all(in_vcf)) warn <- c(warn, sprintf("%d lines have no genotypes and are left out: %s", sum(!in_vcf),
                                            paste(utils::head(means$line[!in_vcf], 5), collapse = ", ")))
  n_pcs <- .ba_gwas_n_pcs(args$n_pcs %||% "auto", design)
  goal <- .ba_resolve_goal("genotypes", design, args)
  design$goal <- goal
  ch <- .ba_engine_chain("genotypes", design, NULL, args$engine, args$allow_fallback)
  backends <- c(rrblup = "rrBLUP", gapit_mlm = "GAPIT", bioana_emmax = "native")
  rc <- .ba_run_chain(ch, function(eng) analyze_gwas(vcf, means[in_vcf & is.finite(means$mean), c("line", "mean")], file.path(outdir, "gwas"),
                                                     trait = "mean", sample_id = "line", n_pcs = n_pcs, min_maf = args$min_maf %||% 0.05,
                                                     fdr_threshold = alpha, backend = backends[[eng]]), "genotypes")
  gw <- rc$result; eng <- rc$engine
  ch$skipped <- c(ch$skipped, lapply(rc$failed, function(m) paste("failed:", m)))
  g <- gw$results
  res <- data.frame(marker_id = g$marker_id, seqname = g$seqname, pos = g$pos, estimate = g$beta, se = g$se,
                    p = g$p_value, padj = g$q_value, maf = g$maf, significant = g$significant, stringsAsFactors = FALSE)
  paths <- c(results = .ba_write_tsv(res, file.path(outdir, "results.tsv")),
             line_means = .ba_write_tsv(means, file.path(outdir, "line_means.tsv")))
  summary <- list(data_type = "genotypes", method = paste("line means (BLUEs) +", gw$summary$backend),
                  model_fitted = lm_out$model, genotype_term = geno, genotype_column = design$parsed$col[[geno]],
                  n_lines = sum(in_vcf & is.finite(means$mean)), n_markers = nrow(gt$genotype),
                  n_features_tested = nrow(res), n_pcs = n_pcs, min_maf = args$min_maf %||% 0.05,
                  null_model_h2 = gw$summary$null_model_h2,
                  significant = list(markers = sum(res$significant, na.rm = TRUE)), warnings = warn,
                  engine = .ba_provenance(eng, design,
                                          sprintf("analyze_gwas(vcf, line_means, backend = '%s', n_pcs = %d, min_maf = %g) after line means from %s",
                                                  backends[[eng]], n_pcs, args$min_maf %||% 0.05, lm_out$model),
                                          sprintf("line means of %d lines; %d markers in %s", sum(in_vcf & is.finite(means$mean)), nrow(gt$genotype), entry$files$vcf),
                                          sum(in_vcf & is.finite(means$mean)), ch,
                                          random = sprintf("line means are BLUEs from %s; kinship random effect in the scan", lm_out$model)),
                  secondary_engine = .ba_secondary_engine(ch$chain, eng),
                  goal = .ba_goal_summary(goal, list(test = "wald", coefs = "each marker", label = NULL)),
                  diagnostics = .ba_diagnostics(NULL, rc$fit_warnings, list(line_means = lm_out$diagnostics)))
  truth <- .ba_qtn_truth(entry, design$dir)
  rec <- NULL
  if (!is.null(truth)) {
    if (all(is.na(res$estimate))) {
      # rrBLUP gives scores only: take marker effects from the native scan (same model) for the truth scores
      nat <- analyze_gwas(vcf, means[in_vcf & is.finite(means$mean), c("line", "mean")], file.path(tempdir(), "bioana_native_effects"),
                          trait = "mean", sample_id = "line", n_pcs = n_pcs, min_maf = args$min_maf %||% 0.05,
                          fdr_threshold = alpha, backend = "native")$results
      k <- match(res$marker_id, nat$marker_id)
      res$estimate <- nat$beta[k]; res$se <- nat$se[k]
    }
    rec <- .ba_qtn_recovery(res, truth, gt$genotype, outdir, alpha, args)
    paths <- c(paths, qtn_recovery = .ba_write_tsv(rec$table, file.path(outdir, "qtn_recovery.tsv")))
    summary$recovery <- rec$summary
    gv <- .ba_true_genetic_values(design$dir)
    if (!is.null(gv)) {
      m <- match(means$line, gv$line)
      ok <- !is.na(m) & is.finite(means$mean)
      if (sum(ok) > 2) summary$line_mean_cor <- round(stats::cor(means$mean[ok], gv$genetic_value[m[ok]]), 3)
    }
  }
  list(summary = summary, paths = paths, fit = lm_out$fit, results = res, recovery = rec$table)
}

# Which model term holds the lines.
.ba_genotype_term <- function(entry, design) {
  p <- design$parsed
  check <- function(v, from) {
    if (!v %in% p$vars) stop(sprintf("%s names genotype term '%s', which is not in analysis.model (%s).", from, v, design$model), call. = FALSE)
    v
  }
  if (!is.null(entry$genotype_term)) return(check(entry$genotype_term, "The manifest entry"))
  st <- file.path(design$dir, "sim_truth.json")
  if (file.exists(st)) {
    v <- jsonlite::read_json(st)$genetics$genotype_term
    if (!is.null(v)) return(check(v, "sim_truth.json"))
  }
  rules <- .ba_rules("genetics")
  norm <- function(x) gsub("[^a-z0-9]+", "_", tolower(trimws(x)))
  aterms <- design$spec$analysis$terms %||% list()
  treat <- p$vars[p$role[p$vars] == "treatment"]
  for (v in treat) {
    cands <- c(if (is.character(aterms[[v]])) aterms[[v]], v, p$col[[v]])
    if (v == "treatment") cands <- c(cands, design$spec$factors[[1]]$name)
    if (any(grepl(rules$genotype_pattern, norm(cands)))) return(v)
  }
  stop(sprintf("No genotype term in analysis.model (%s): no treatment term matches genetics.json genotype_pattern. Give genotype_term in the manifest entry.",
               design$model), call. = FALSE)
}

.ba_gwas_n_pcs <- function(n_pcs, design) {
  if (!identical(as.character(n_pcs), "auto")) return(as.integer(n_pcs))
  rules <- .ba_rules("genetics")
  mating <- tolower(unlist(design$spec$population$mating %||% character()))
  st <- file.path(design$dir, "sim_truth.json")
  src <- if (file.exists(st)) jsonlite::read_json(st)$genetics$genotype_source$source else NULL
  if (identical(src, "cross") || any(mating %in% unlist(rules$cross_mating))) 0L else 3L
}

# Line BLUEs from the design's model.
.ba_line_means <- function(design, geno, entry, args) {
  warn <- character()
  p <- design$parsed
  fixed_vars <- unique(unlist(lapply(p$fixed_labels, function(l) all.vars(str2lang(l)))))
  if (!geno %in% fixed_vars) {
    bars <- .ba_findbars(p$formula)
    keep <- Filter(function(b) !geno %in% all.vars(b[[3]]), bars)
    rhs <- c(p$fixed_labels, geno, vapply(keep, function(b) paste0("(", paste(deparse(b), collapse = ""), ")"), ""))
    p$formula <- stats::reformulate(rhs, response = p$response)
    p$fixed_labels <- c(p$fixed_labels, geno)
    warn <- c(warn, sprintf("'%s' is random in analysis.model; refitted as fixed to get line means (BLUEs) for the GWAS", geno))
    design$parsed <- p
  }
  type <- entry$outcome_type %||% design$spec$primary_outcome$type %||% "continuous"
  if (!type %in% c("continuous", "binary", "count")) {
    stop(sprintf("Line means for %s outcomes are not supported yet (continuous, binary, count).", type), call. = FALSE)
  }
  if (type != "continuous") warn <- c(warn, sprintf("%s outcome: line means on the link scale", type))
  oe <- list(id = "outcome", data_type = "outcome", column = entry$outcome_column %||% "primary_outcome", outcome_type = type)
  out <- .ba_method_outcome(oe, design, file.path(tempdir(), "bioana_line_means"), args = list(subsamples = "aggregate"))
  fit <- out$fit
  warn <- c(warn, out$summary$warnings)
  frame <- stats::model.frame(fit)
  beta <- if (inherits(fit, "merMod")) lme4::fixef(fit) else stats::coef(fit)
  V <- as.matrix(stats::vcov(fit))
  est <- !is.na(beta)
  fixed_f <- .ba_nobars(p$formula)
  tt <- stats::delete.response(stats::terms(fixed_f))
  other <- setdiff(intersect(all.vars(tt), names(frame)), geno)
  levs <- lapply(c(geno, other), function(v) {
    x <- frame[[v]]
    if (is.factor(x)) levels(x) else mean(x, na.rm = TRUE)
  })
  names(levs) <- c(geno, other)
  grid <- expand.grid(levs, stringsAsFactors = FALSE)
  for (v in names(levs)) if (is.factor(frame[[v]])) grid[[v]] <- factor(grid[[v]], levels = levels(frame[[v]]))
  X <- stats::model.matrix(tt, grid)
  X <- X[, names(beta)[est], drop = FALSE]
  L <- rowsum(X, as.character(grid[[geno]]), reorder = FALSE) / as.vector(table(factor(grid[[geno]], levels = unique(grid[[geno]]))))
  b <- beta[est]
  Vb <- V[names(b), names(b), drop = FALSE]
  n_units <- table(frame[[geno]])
  means <- data.frame(line = rownames(L), mean = drop(L %*% b), se = sqrt(pmax(0, rowSums((L %*% Vb) * L))),
                      n_units = as.integer(n_units[rownames(L)]), stringsAsFactors = FALSE)
  means <- means[means$n_units > 0, , drop = FALSE]
  rownames(means) <- NULL
  list(means = means, fit = fit, model = out$summary$model_fitted, warnings = warn, diagnostics = out$summary$diagnostics)
}

# QTN truth (kind qtn) or NULL.
.ba_qtn_truth <- function(entry, dir) {
  tr <- entry$truth
  if (is.null(tr) || is.null(tr$file) || !file.exists(file.path(dir, tr$file))) return(NULL)
  x <- utils::read.delim(file.path(dir, tr$file), stringsAsFactors = FALSE, check.names = FALSE)
  if (!nrow(x)) return(NULL)
  if (is.null(x$type)) x$type <- "additive"
  x
}

.ba_true_genetic_values <- function(dir) {
  st <- file.path(dir, "sim_truth.json")
  if (!file.exists(st)) return(NULL)
  f <- jsonlite::read_json(st)$genetics$genetic_values
  if (is.null(f) || !file.exists(file.path(dir, f))) return(NULL)
  utils::read.delim(file.path(dir, f), stringsAsFactors = FALSE, colClasses = c(line = "character"))
}

# Score the scan against the QTN.
.ba_qtn_recovery <- function(res, truth, dosage, outdir, alpha, args) {
  rules <- .ba_rules("genetics")
  truth <- truth[truth$type != "gxe_slope", , drop = FALSE]
  truth <- truth[!duplicated(truth$marker_id), , drop = FALSE]
  target <- truth$type %in% unlist(rules$target_qtn_types)
  gw <- data.frame(marker_id = res$marker_id, seqname = res$seqname, pos = res$pos, q_value = res$padj, stringsAsFactors = FALSE)
  # every QTN (epistatic too) decides what a false positive is; only targets count in tpr
  bm <- benchmark_gwas(gw, truth[, c("marker_id", "seqname", "pos")], file.path(outdir, "gwas"),
                       fdr_threshold = alpha, window_bp = args$window_bp %||% 0L,
                       ld_r2 = args$ld_r2 %||% 0.5, genotypes = dosage)
  m <- match(truth$marker_id, res$marker_id)
  tab <- data.frame(marker_id = truth$marker_id, seqname = truth$seqname, pos = truth$pos, type = truth$type,
                    effect = truth$effect, maf = truth$maf %||% NA_real_, tested = !is.na(m),
                    recovered = bm$truth$recovered, estimate = res$estimate[m], se = res$se[m],
                    p = res$p[m], padj = res$padj[m], stringsAsFactors = FALSE)
  tgt <- tab[target & tab$tested, , drop = FALSE]
  est_ok <- is.finite(tgt$estimate)
  mt <- bm$metrics
  summary <- list(
    qtn = sum(target), testable = nrow(tgt), untestable = sum(target & !tab$tested), epistatic = sum(truth$type == "epistatic"),
    true_positives = sum(tgt$recovered), tpr = if (nrow(tgt)) round(mean(tgt$recovered), 3) else NA,
    discoveries = mt$discoveries, false_positives = mt$false_positive_markers,
    fdr = if (mt$discoveries) round(mt$false_positive_markers / mt$discoveries, 3) else NA,
    sign_agreement = if (any(est_ok)) round(mean(sign(tgt$estimate[est_ok]) == sign(tgt$effect[est_ok])), 3) else NA,
    estimate_cor = if (sum(est_ok) > 2) round(stats::cor(tgt$estimate[est_ok], tgt$effect[est_ok]), 3) else NA,
    ld_r2 = args$ld_r2 %||% 0.5, window_bp = args$window_bp %||% 0L)
  list(table = tab, summary = summary)
}
