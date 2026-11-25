# Pre-analysis design checks: run before every data set is fitted. A failed check stops the analysis
# with the reason; warnings go into summary$diagnostics. Thresholds: inst/rules/design_checks.json.
#
# 1. estimability: no treatment coefficient may be confounded with other terms (rank-deficient design)
# 2. residual df: at least 1 (else stop), at least min_residual_df (else warn)
# 3. replication: experimental units per treatment level (warn below min_units_per_level)
# 4. repeated units: a subject column that repeats across samples must be a random term, or the
#    samples are treated as independent when they are not (pseudo-replication)
# 5. random terms: each needs more than one observation per level and at least 2 levels

#' Check a design folder before analysis
#'
#' The same checks [bioana_analyze()] runs before fitting each data set, on
#' the whole sample metadata: estimability (no treatment coefficient
#' confounded with another term), residual degrees of freedom, experimental
#' units per treatment level, repeated units that the model does not treat as
#' repeated, and random terms that cannot be estimated.
#'
#' @param dir A design folder.
#' @return A data frame: `check`, `status` (`ok`, `warn`, `stop`), `detail`.
#' @export
bioana_check_design <- function(dir) {
  design <- .ba_read_biodes(dir)
  .ba_check_design(design, stop_on_fail = FALSE)
}

.ba_check_design <- function(design, sample_ids = NULL, what = c("features", "outcome"), meta = NULL, stop_on_fail = TRUE) {
  what <- match.arg(what)
  p <- design$parsed
  meta <- meta %||% design$meta
  if (!is.null(sample_ids)) meta <- meta[meta$sample_id %in% sample_ids, , drop = FALSE]
  rules <- .ba_rules("design_checks")
  rows <- list()
  add <- function(check, status, detail) rows[[length(rows) + 1L]] <<- data.frame(check = check, status = status, detail = detail, stringsAsFactors = FALSE)
  frame <- .ba_frame(meta, p)
  fvars <- intersect(unlist(lapply(p$fixed_labels, function(l) all.vars(str2lang(l)))), names(frame))
  ok <- stats::complete.cases(frame[, fvars, drop = FALSE])
  fr <- frame[ok, , drop = FALSE]; mt <- meta[ok, , drop = FALSE]

  # 1. estimability
  X <- stats::model.matrix(stats::reformulate(p$fixed_labels), data = fr)
  q <- qr(X)
  if (q$rank < ncol(X)) {
    keep <- colnames(X)[q$pivot[seq_len(q$rank)]]
    for (a in colnames(X)[q$pivot[(q$rank + 1L):ncol(X)]]) {
      b <- qr.coef(qr(X[, keep, drop = FALSE]), X[, a])
      with <- names(b)[!is.na(b) & abs(b) > 1e-8]
      msg <- sprintf("%s cannot be estimated: in these samples it equals a combination of %s", a, paste(with, collapse = ", "))
      # stop when a treatment coefficient is part of the dependency, whichever column the QR dropped
      involved <- setdiff(c(a, with), "(Intercept)")
      if (any(vapply(involved, function(t) .ba_coef_role(t, p), "") == "treatment")) {
        add("estimability", "stop", paste0(msg, ". The treatment is confounded with another term: the design cannot separate them."))
      } else add("estimability", "warn", paste0(msg, " (a nuisance term; dropped)"))
    }
  } else add("estimability", "ok", sprintf("%d fixed coefficients, all estimable", ncol(X)))

  # 2. residual degrees of freedom (fixed terms; random terms take more)
  n_obs <- length(unique(.ba_unit_id(mt)))
  df <- n_obs - q$rank
  if (df < 1) add("residual_df", "stop", sprintf("%d experimental units for %d fixed coefficients: no residual degrees of freedom", n_obs, q$rank))
  else if (df < rules$min_residual_df) add("residual_df", "warn", sprintf("only %d residual degrees of freedom (%d units, %d coefficients)", df, n_obs, q$rank))
  else add("residual_df", "ok", sprintf("%d residual degrees of freedom", df))

  # 3. experimental units per treatment level
  unit <- .ba_unit_id(mt)
  for (v in p$vars[p$role[p$vars] == "treatment" & !p$vars %in% p$numeric]) {
    if (!v %in% names(fr) || !is.factor(fr[[v]])) next
    n <- tapply(unit, droplevels(fr[[v]]), function(u) length(unique(u)))
    low <- names(n)[n < rules$min_units_per_level]
    if (length(low)) add("replication", "warn", sprintf("%s: %s with %s unit(s): its effect rests on the variance of the other levels",
                                                         v, paste(low, collapse = ", "), paste(n[low], collapse = ", ")))
    else add("replication", "ok", sprintf("%s: %d to %d units per level", v, min(n), max(n)))
  }

  # 4. repeated units the model does not treat as repeated
  rvars <- unique(unlist(lapply(p$random, `[[`, "vars")))
  rcols <- unique(unlist(p$col[intersect(rvars, names(p$col))]))
  model_cols <- unique(c(unlist(p$col), rcols))
  subj <- grep(rules$subject_pattern, names(mt), value = TRUE, ignore.case = TRUE)
  subj <- setdiff(subj, c(setdiff(model_cols, rcols), "sample_id"))
  for (s in subj) {
    if (s %in% rcols) next
    g <- as.character(mt[[s]])
    if (!anyDuplicated(g[!is.na(g) & nzchar(g)])) next
    cell <- do.call(paste, c(lapply(fvars, function(v) as.character(fr[[v]])), sep = "\r"))
    varies <- tapply(cell, g, function(x) length(unique(x)) > 1)
    same <- tapply(cell, g, function(x) length(x) > 1 && length(unique(x)) == 1)
    if (any(varies, na.rm = TRUE)) {
      add("repeated_units", "stop", sprintf("'%s' repeats across samples that differ in %s (e.g. %s): these are repeated measures of one unit, but analysis.model has no random term for '%s'. Add (1 | subject) with subject = '%s' to the model (Biodes), or the samples are wrongly treated as independent.",
                                            s, paste(fvars, collapse = "/"), names(which(varies))[1], s, s))
    } else if (any(same, na.rm = TRUE) && !(what == "outcome" && identical(s, "unit_id"))) {
      add("technical_replicates", "stop", sprintf("'%s' repeats across samples with the same treatment (e.g. %s): technical replicates or subsamples of one unit, not independent replicates. Sum (counts) or average them per unit first, or give each its own unit.",
                                                  s, names(which(same))[1]))
    }
  }

  # 5. random terms
  for (r in p$random) {
    vars <- intersect(r$vars, names(fr))
    if (!length(vars)) next
    g <- droplevels(interaction(fr[vars], drop = TRUE))
    tab <- table(g)
    if (nlevels(g) < 2) add("random_terms", "stop", sprintf("(%s) has %d level: a random term needs at least 2", r$label, nlevels(g)))
    else if (all(tab == 1)) add("random_terms", "stop", sprintf("(%s) has one observation per level: it cannot be separated from the residual", r$label))
    else if (nlevels(g) < rules$min_random_levels) add("random_terms", "warn", sprintf("(%s) has %d levels: its variance is poorly estimated", r$label, nlevels(g)))
    else add("random_terms", "ok", sprintf("(%s): %d levels", r$label, nlevels(g)))
  }

  out <- do.call(rbind, rows)
  if (stop_on_fail && any(out$status == "stop")) {
    stop(paste0("Design check failed: ", paste(out$detail[out$status == "stop"], collapse = " | ")), call. = FALSE)
  }
  out
}
