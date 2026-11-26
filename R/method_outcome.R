# The primary outcome: the model Biodes wrote (analysis.model), fitted with the
# family the outcome type calls for. Continuous: lm / lmer. Binary: glm /
# glmer (binomial). Count: negative binomial glm (MASS) / glmer.nb. Ordinal:
# cumulative logit, ordinal::clm / clmm (MASS::polr without random terms when
# ordinal is missing, fixed designs only). Survival: Cox model (survival::coxph);
# with random terms coxme, or one random intercept as a gamma frailty. The engine
# chain (engine_selection.tsv) picks the fit; nothing drops random terms unless
# allow_fallback = TRUE.

.ba_method_outcome <- function(entry, design, outdir, alpha = 0.05, args = list()) {
  spec <- design$spec; meta <- design$meta
  subsamples <- match.arg(args$subsamples %||% "aggregate", c("aggregate", "nest", "as_replicates"))
  col <- entry$column %||% "primary_outcome"
  type <- entry$outcome_type %||% spec$primary_outcome$type %||% "continuous"
  goal <- .ba_resolve_goal(paste0("outcome:", type), design, args)
  design <- .ba_goal_design(design, goal)
  goal <- design$goal
  p <- design$parsed
  warn <- character()
  y <- if (type == "nominal") { z <- trimws(as.character(meta[[col]])); z[z == ""] <- NA; z } else suppressWarnings(as.numeric(meta[[col]]))
  if (type == "binary") {
    # two named categories (Biodes 0.6.7 primary_outcome.categories, else expected_ratio labels): 0 = the first
    # (reference), 1 = the second
    z <- trimws(as.character(meta[[col]])); cats <- as.character(unlist(design$spec$primary_outcome$categories %||% names(design$spec$analysis$expected_ratio)))
    lab <- !is.na(z) & z != "" & is.na(suppressWarnings(as.numeric(z)))
    if (any(lab) && length(cats) == 2L && all(z[lab] %in% cats)) {
      y[lab] <- match(z[lab], cats) - 1
      warn <- c(warn, sprintf("binary outcome coded 0 = %s (reference), 1 = %s (primary_outcome.categories)", cats[1], cats[2]))
    }
  }
  if (!any(!is.na(y))) stop(sprintf("Column '%s' of the sample metadata is empty: nothing to fit.", col), call. = FALSE)
  meta$y_ <- y
  if (type == "survival") {
    ev <- suppressWarnings(as.numeric(meta[[entry$event_column %||% "event"]]))
    if (all(is.na(ev))) { ev <- rep(1, nrow(meta)); warn <- c(warn, "no event column: every time treated as an event") }
    meta$event_raw_ <- ev
    if (any(ev > 1, na.rm = TRUE)) {
      # competing events (codes 2+): the Cox model is cause-specific for code 1; survival_tests.R uses the codes
      warn <- c(warn, sprintf("competing events (codes %s): the Cox model is cause-specific for event code 1", paste(sort(setdiff(unique(ev), c(0, 1, NA))), collapse = ", ")))
      ev <- as.numeric(ev == 1)
    }
    meta$event_ <- ev
  }
  unit <- .ba_outcome_unit(meta, spec)
  cell <- do.call(paste, c(list(unit), lapply(p$vars, function(v) meta[[p$col[[v]]]]), sep = "\r"))
  has_subs <- anyDuplicated(cell[!is.na(unit)]) > 0
  if (has_subs && subsamples == "as_replicates") {
    stop(sprintf("These rows are subsamples of %d %s units. Biodes set the experimental unit to '%s', so subsamples are not independent replicates. Use subsamples = 'aggregate' or 'nest'.",
                 length(unique(unit)), spec$experimental_unit %||% "experimental", spec$experimental_unit %||% "unit"), call. = FALSE)
  }
  dat <- meta[!is.na(y), , drop = FALSE]
  dcell <- cell[!is.na(y)]
  if (has_subs && type == "survival") {
    # subsamples of a survival outcome (seeds on a tray): averaging censored times is wrong; every subsample is
    # kept and the unit is a cluster (shared frailty below, robust score test in survival_tests.R)
    if (subsamples != "nest") warn <- c(warn, sprintf("survival subsamples are not averaged: %d units kept as clusters of their subsamples (shared frailty)", length(unique(unit))))
    subsamples <- "nest"
  }
  if (has_subs && type == "nominal") {
    warn <- c(warn, "subsamples of a nominal outcome are counted one by one (categories cannot be averaged); the unit cluster is ignored")
    has_subs <- FALSE
  }
  if (has_subs && type != "survival" && (subsamples == "aggregate" || type != "continuous")) {
    if (type != "continuous" && subsamples == "nest") warn <- c(warn, "subsamples averaged per unit: nesting is fitted for continuous outcomes only")
    first <- dat[!duplicated(dcell), , drop = FALSE]
    first$y_ <- as.numeric(tapply(dat$y_, dcell, mean)[dcell[!duplicated(dcell)]])
    if (type %in% c("binary", "count", "ordinal")) first$y_ <- round(first$y_)
    dat <- first
  }
  f <- p$formula; f[[2]] <- as.name("y_")
  if (has_subs && subsamples == "nest" && type %in% c("continuous", "survival")) {
    dat$unit_id_ <- dcell
    f <- stats::update(f, . ~ . + (1 | unit_id_))
  }
  frame <- .ba_frame(dat, p)
  if (type == "survival") { frame$event_ <- dat$event_; frame$event_raw_ <- dat$event_raw_ }
  .ba_check_outcome_values(frame$y_, frame$event_, type, col)
  if (type == "ordinal") frame$y_ <- factor(frame$y_, ordered = TRUE)
  if (type == "nominal") {
    cats <- as.character(unlist(spec$primary_outcome$categories))
    lv <- unique(frame$y_[!is.na(frame$y_)]); lv <- c(intersect(cats, lv), sort(setdiff(lv, cats)))
    frame$y_ <- factor(frame$y_, levels = lv)
  }
  if (type == "survival") {
    # fixed block terms in a Cox model (survival_rules.json): many blocks -> shared frailty, few -> strata (surv_rhs)
    sr_ <- .ba_rules("survival_rules")
    p$surv_strata <- .ba_survival_strata(c(design$spec$analysis$survival$strata, args$survival$strata), p, frame)
    in_model_ <- unique(unlist(lapply(p$fixed_labels, function(l) all.vars(str2lang(l)))))
    if (length(s_ <- setdiff(p$surv_strata, in_model_))) warn <- c(warn, sprintf("Cox model stratified on %s (analysis.survival.strata)", paste(s_, collapse = ", ")))
    for (l in p$fixed_labels) {
      v <- all.vars(str2lang(l))
      if (length(v) != 1L || !identical(p$role[[v]], "block")) next
      nb <- length(unique(stats::na.omit(frame[[v]])))
      if (nb >= (sr_$frailty_min_blocks %||% 10)) {
        f <- stats::update(f, stats::as.formula(sprintf(". ~ . - %s + (1 | %s)", v, v)))
        p$fixed_labels <- setdiff(p$fixed_labels, l)
        warn <- c(warn, sprintf("block term %s (%d blocks) fitted as a shared frailty (1 | %s): fixed coefficients for many small blocks diverge (survival_rules.json)", v, nb, v))
      } else warn <- c(warn, sprintf("block term %s (%d blocks) fitted as strata(%s) in the Cox model", v, nb, v))
    }
  }
  p_stat <- p   # the stat layer sees the design's blocks (stratified tests) even when the model drops them
  if (type %in% c("binary", "ordinal", "nominal")) {
    # small fixed blocks bias odds ratios (stat_rules.json categorical_block_terms): random (binary, ordinal) or dropped (nominal)
    cb <- .ba_rules("stat_rules")$categorical_block_terms
    for (l in p$fixed_labels) {
      v <- all.vars(str2lang(l))
      if (length(v) != 1L || !identical(p$role[[v]], "block")) next
      per <- sum(!is.na(frame$y_)) / max(1, length(unique(stats::na.omit(frame[[v]]))))
      if (per >= (cb$min_units_per_block %||% 10)) next
      act <- cb[[type]] %||% "random"
      f <- stats::update(f, stats::as.formula(if (act == "random") sprintf(". ~ . - %s + (1 | %s)", v, v) else sprintf(". ~ . - %s", v)))
      p$fixed_labels <- setdiff(p$fixed_labels, l)
      warn <- c(warn, sprintf("block term %s (%.1f units per block) %s: fixed coefficients for small blocks bias the odds ratios away from 1 (stat_rules.json categorical_block_terms)",
                              v, per, if (act == "random") sprintf("fitted as a random intercept (1 | %s); odds ratios are conditional on the block", v)
                                      else "dropped from the model (marginal odds ratios); the test stratifies on it"))
    }
  }
  has_random <- length(.ba_findbars(f)) > 0
  fixed_f <- .ba_nobars(f)
  # the chain sees the random terms of the formula actually fitted (nested subsamples add one)
  d2 <- design
  d2$parsed$random <- lapply(.ba_findbars(f), function(b) list(label = paste(deparse(b[[3]]), collapse = ""), vars = all.vars(b[[3]]),
                                                               intercept_only = identical(paste(deparse(b[[2]]), collapse = ""), "1")))
  checks <- .ba_check_design(d2, meta = dat, what = "outcome")
  col_term <- .ba_column_terms(stats::reformulate(c("1", p$fixed_labels)), frame)
  plan <- .ba_goal_plan(goal, names(col_term), col_term, p)
  ch <- .ba_engine_chain(paste0("outcome:", type), d2, NULL, args$engine, args$allow_fallback, test = plan$test)
  rc <- .ba_run_chain(ch, function(eng) .ba_fit_outcome(eng, f, fixed_f, frame, type, p), "outcome")
  fit <- rc$result; eng <- rc$engine
  goal_result <- NULL
  if (identical(plan$test, "lrt")) {
    # the reduced model: the same engine on the full model matrix without the goal's columns (a formula
    # without the terms would re-code the rest, e.g. treatment + treatment:time spans the full model)
    blk_lab <- if (type == "survival") p$fixed_labels[vapply(p$fixed_labels, function(l) { v <- all.vars(str2lang(l)); length(v) == 1L && identical(p$role[[v]], "block") }, logical(1))] else character()
    rhs <- stats::reformulate(c("1", setdiff(p$fixed_labels, blk_lab)))
    X <- stats::model.matrix(rhs, stats::model.frame(rhs, frame, na.action = stats::na.pass))
    X0 <- X[, !colnames(X) %in% c(plan$drop, "(Intercept)"), drop = FALSE]
    fr0 <- frame
    fixed0 <- c(if (ncol(X0)) { fr0$X0_ <- X0; "X0_" }, blk_lab)
    bars <- vapply(.ba_findbars(f), function(b) paste0("(", paste(deparse(b), collapse = ""), ")"), "")
    f0 <- stats::reformulate(c(if (length(fixed0)) fixed0 else "1", bars), response = "y_")
    p0 <- p; p0$fixed_labels <- fixed0
    goal_result <- tryCatch({
      fit0 <- .ba_fit_outcome(eng, f0, .ba_nobars(f0), fr0, type, p0)
      .ba_lrt(fit0, fit)
    }, error = function(e) list(error = conditionMessage(e)))
    if (!is.null(goal_result$error)) warn <- c(warn, paste("joint test failed:", goal_result$error))
  }
  if (length(rc$failed)) warn <- c(warn, sprintf("engine %s failed (%s); used %s, which fits the design fully", names(rc$failed), unlist(rc$failed), eng))
  if (identical(ch$fit[[eng]], "drops_random")) warn <- c(warn, sprintf("random terms dropped by %s (allow_fallback = TRUE): %s", eng, paste(.ba_random_labels(d2$parsed), collapse = ", ")))
  diag_engine <- list()
  if (has_random && inherits(fit, "merMod") && lme4::isSingular(fit)) {
    warn <- c(warn, "singular fit: a random-effect variance was estimated as ~0")
    diag_engine$singular <- TRUE
  }
  if (inherits(fit, "glmmTMB") && !isTRUE(fit$sdr$pdHess)) {
    warn <- c(warn, "glmmTMB: the Hessian is not positive definite (convergence problem)")
    diag_engine$pdHess <- FALSE
  }
  if (type == "survival") {
    pat <- paste(unlist(.ba_rules("survival_rules")$nonconvergence_patterns), collapse = "|")
    if (nzchar(pat) && any(grepl(pat, rc$fit_warnings %||% character(), ignore.case = TRUE))) {
      warn <- c(warn, "Cox fit did not converge (a coefficient may be infinite): estimates unreliable, goal truth not scored")
      diag_engine$converged <- FALSE
    }
  }
  res <- .ba_coef_table(fit, type, frame, p)
  res$role <- vapply(res$term, function(t) .ba_coef_role(t, p), "")
  res$padj <- rep(NA_real_, nrow(res))
  test <- if (identical(plan$test, "wald")) res$term %in% plan$coefs else res$role == "treatment"
  res$padj[test] <- stats::p.adjust(res$p[test], "BH")
  paths <- c(results = .ba_write_tsv(res, file.path(outdir, "results.tsv")))
  if (!is.null(goal_result) && is.null(goal_result$error)) {
    paths <- c(paths, goal_test = .ba_write_tsv(data.frame(goal = goal$id, tested = plan$label, statistic = goal_result$statistic,
                                                         df = goal_result$df, p = goal_result$p, method = goal_result$method),
                                              file.path(outdir, "goal_test.tsv")))
  }
  # classical ANOVA next to the model (stat_tests.tsv): repeated-measures / mixed / factorial, post-hoc, effect sizes
  st <- NULL
  reg <- NULL
  if (.ba_has_covariates(p) && type %in% c("continuous", "binary") && !isFALSE(args$regression)) {
    # measured predictors (covariates): the regression layer instead of the group tests (R/regression.R)
    reg <- tryCatch(.ba_regression(fit, frame, p, design, alpha, type, args), error = function(e) list(summary = list(status = "error", note = conditionMessage(e))))
    for (k in c("coefficients", "type2", "influence", "influence_refit", "vif", "simple_slopes", "adjusted_means", "predictions", "robust"))
      if (!is.null(reg[[k]])) paths <- c(paths, stats::setNames(.ba_write_tsv(reg[[k]], file.path(outdir, paste0("regression_", k, ".tsv"))), paste0("regression_", k)))
    if (!is.null(reg$summary$sensitivity$slopes)) paths <- c(paths, regression_sensitivity = .ba_write_tsv(reg$summary$sensitivity$slopes, file.path(outdir, "regression_sensitivity.tsv")))
  }
  if (type %in% c("continuous", "ordinal", "binary", "nominal", "survival") && !isFALSE(args$anova) && !isFALSE(args$stat_tests) && is.null(reg)) {
    st <- tryCatch(if (type == "survival") { fr <- frame; fr$event_ <- fr$event_raw_ %||% fr$event_; .ba_survival_tests(fr, design$parsed, design, goal, alpha, args, fit) }
                   else .ba_stat_tests(frame, p_stat, design, goal, alpha, type, args),
                   error = function(e) list(summary = list(status = "error", note = conditionMessage(e))))
    for (k in c("km_curve", "survival_summary", "landmarks", "life_table", "cif_curve", "competing_table", "cells")) if (!is.null(st[[k]])) paths <- c(paths, stats::setNames(.ba_write_tsv(st[[k]], file.path(outdir, paste0(k, ".tsv"))), k))
    if (!is.null(st$descriptives)) paths <- c(paths, descriptives = .ba_write_tsv(st$descriptives, file.path(outdir, "descriptives.tsv")))
    if (!is.null(st$stat_table)) paths <- c(paths, stat_tests = .ba_write_tsv(st$stat_table, file.path(outdir, "stat_tests.tsv")))
    if (!is.null(st$anova_table)) paths <- c(paths, anova_table = .ba_write_tsv(st$anova_table, file.path(outdir, "anova_table.tsv")))
    if (!is.null(st$posthoc)) paths <- c(paths, posthoc = .ba_write_tsv(st$posthoc, file.path(outdir, "posthoc.tsv")))
  }
  aov <- if (type == "continuous") tryCatch(as.data.frame(stats::anova(fit)), error = function(e) NULL) else NULL
  if (!is.null(aov) && has_random && !"Pr(>F)" %in% names(aov) && "F value" %in% names(aov)) {
    # lme4 without lmerTest gives F only: denominator df = rows - fixed coefficients
    aov$DenDF <- max(1, nrow(frame) - ncol(stats::model.matrix(stats::reformulate(c("1", p$fixed_labels)), data = frame)))
    aov[["Pr(>F)"]] <- stats::pf(aov[["F value"]], aov$npar, aov$DenDF, lower.tail = FALSE)
  }
  summary <- list(data_type = "outcome", outcome_type = type, response = spec$primary_outcome$name %||% col,
                  method = paste(class(fit)[1], if (type == "survival") "(Cox proportional hazards)" else ""),
                  model_fitted = .ba_formula_text(attr(fit, "ba_formula") %||% stats::formula(fit)),
                  n_rows = nrow(meta), n_units_fitted = nrow(frame), subsamples = if (has_subs) subsamples else "none",
                  dropped_samples = meta$sample_id[is.na(y)],
                  engine = .ba_provenance(eng, d2, tryCatch(paste(deparse(stats::getCall(fit)), collapse = " "), error = function(e) class(fit)[1]),
                                          sprintf("%s outcome '%s', %d units", type, col, nrow(frame)), nrow(frame), ch,
                                          formula = paste(deparse(f), collapse = "")),
                  diagnostics = .ba_diagnostics(checks, rc$fit_warnings, if (length(diag_engine)) diag_engine),
                  goal = c(.ba_goal_summary(goal, plan), if (!is.null(goal_result)) list(result = goal_result),
                           if (!is.null(goal$effect)) list(effect_scale_used = .ba_tsv("engines")$estimate_scale[.ba_tsv("engines")$engine == eng])),
                  secondary_engine = .ba_secondary_engine(ch$chain, eng),
                  stat_test = st$summary, regression = reg$summary,
                  anova = if (is.null(aov)) NULL else lapply(split(cbind(term = rownames(aov), aov), seq_len(nrow(aov))), as.list),
                  warnings = warn)
  truth <- .ba_entry_truth(entry, design$dir)
  rec <- NULL
  if (!is.null(truth) && identical(st$summary$status, "run")) {
    if (!is.null(st$posthoc) && identical(st$summary$posthoc$package, "emmeans") && type == "continuous") {
      pr <- tryCatch(.ba_posthoc_recovery(st$posthoc, st$summary$posthoc, truth, frame, alpha), error = function(e) NULL)
      if (!is.null(pr)) {
        summary$stat_test$posthoc$recovery <- pr$summary
        paths[["posthoc"]] <- .ba_write_tsv(pr$table, file.path(outdir, "posthoc.tsv"))
      }
    }
    if (!is.null(truth$within_subject)) summary$stat_test$truth_within_subject <- truth$within_subject[intersect(c("structure", "rho", "slope_sd", "over"), names(truth$within_subject))]
  }
  if (!is.null(truth) && isTRUE(goal$none) && length(truth$outcome_model$true_ratio)) {
    summary$goal_recovery <- .ba_ratio_recovery(truth$outcome_model, st$summary, frame, alpha)
    if (!is.null(summary$goal_recovery$table)) { paths[["ratio"]] <- .ba_write_tsv(summary$goal_recovery$table, file.path(outdir, "ratio_recovery.tsv")); summary$goal_recovery$table <- NULL }
  }
  if (!is.null(truth) && !is.null(entry$truth$goal) && !isTRUE(goal$none)) {
    gt <- truth[[entry$truth$goal$field %||% "goal_truth"]]
    if (!is.null(gt)) {
      # Simitall's survival truth is a time ratio (AFT): score it against the AFT model of the truth's family
      aft <- if (type == "survival" && grepl("time ratio", gt$effect_scale %||% "")) .ba_aft_for_truth(frame, p, truth$outcome_model$family %||% "lognormal_aft") else NULL
      summary$goal_recovery <- .ba_goal_recovery_outcome(gt, summary$goal, if (!is.null(aft)) aft$res else res, alpha, if (!is.null(aft)) "aft" else type)
      if (!is.null(aft)) summary$goal_recovery$scale <- aft$scale
      if (type == "survival" && (isFALSE(gt$proportional_hazards) || isFALSE(goal$effect$proportional_hazards)))
        summary$goal_recovery <- .ba_goal_recovery_nonph(summary$goal_recovery, gt, truth$outcome_model, goal, design, frame, p, res, alpha)
      if (type %in% c("nominal", "binary") && length(truth$outcome_model$group_probability)) {
        gpr <- .ba_group_probability_recovery(truth$outcome_model$group_probability, frame, p, alpha)
        if (isTRUE(gpr$scored)) {
          paths[["group_probability"]] <- .ba_write_tsv(gpr$table, file.path(outdir, "group_probability_recovery.tsv"))
          summary$goal_recovery[c("estimate", "true_effect", "bias", "sign_agreement")] <- NULL
          summary$goal_recovery$scale <- "class probability per group (Simitall group_probability)"
          summary$goal_recovery$probabilities <- gpr[c("max_abs_diff", "covered", "cells", "note")]
          summary$goal_recovery$covered <- gpr$covered == gpr$cells
        } else summary$goal_recovery$probabilities <- gpr
      }
      ae <- truth$adjustment_effects %||% gt$adjustment_effects
      if (length(ae) && type == "continuous") summary$goal_recovery$adjustment <- .ba_level_effect_recovery(ae, res, alpha)
      if (identical(goal$id, "regression") && (length(gt$true_slopes) || length(gt$true_log_odds_ratios))) {
        sl <- .ba_slope_recovery(gt, res, alpha, type, p)
        if (!is.null(sl)) {
          summary$goal_recovery$slopes <- sl
          if (nrow(sl) == 1L) summary$goal_recovery[c("estimate", "true_effect", "covered")] <- list(sl$estimate, sl$true, sl$covered)
          else summary$goal_recovery$covered <- all(sl$covered)
          summary$goal_recovery$scale <- if (type == "binary") "odds ratio per unit of each covariate (compared on the log scale)" else "slope per unit of each covariate (as fitted, after any transform)"
        }
      }
      if (type == "binary" && length(truth$outcome_model$conditional_odds_ratios)) {
        orr <- .ba_odds_ratio_recovery(truth$outcome_model$conditional_odds_ratios, res, alpha, has_random)
        if (!is.null(orr)) {
          summary$goal_recovery$odds_ratios <- orr$table
          if (nrow(orr$table) == 1L) summary$goal_recovery[c("estimate", "true_effect", "covered")] <- list(orr$table$estimate, orr$table$true_or, orr$table$covered)
          summary$goal_recovery$scale <- orr$scale
        }
      }
      if (type == "survival" && length(gt$groups %||% truth$outcome_model$groups)) {
        sgr <- .ba_survival_group_recovery(gt$groups %||% truth$outcome_model$groups, frame, p, alpha, tau = suppressWarnings(as.numeric(truth$outcome_model$rmst_tau %||% NA)) %||% NULL)
        if (!is.null(sgr)) {
          paths[["survival_groups"]] <- .ba_write_tsv(sgr, file.path(outdir, "survival_group_recovery.tsv"))
          summary$goal_recovery$groups <- list(covered = sum(sgr$covered), measures = nrow(sgr),
                                               by_measure = lapply(split(sgr, trimws(sub(" (at|to) .*$| \\(.*$", "", sgr$measure))), function(x) list(covered = sum(x$covered), n = nrow(x))),
                                               note = "per treatment cell vs Simitall outcome_model.groups; intervals treat subsamples as independent")
        }
      }
      if (isFALSE(diag_engine$converged)) { summary$goal_recovery$scored <- FALSE
        summary$goal_recovery$note <- paste(c(summary$goal_recovery$note, "Cox fit did not converge: not scored"), collapse = "; ") }
    }
  }
  # a numeric (trend) term has no per-level truth: its slope is scored by goal_recovery
  if (!is.null(truth) && length(intersect(p$numeric, p$vars)) && !is.null(summary$goal_recovery)) truth <- NULL
  if (!is.null(truth)) {
    rec <- .ba_outcome_recovery(res, truth, frame, p, alpha, type)
    if (!is.null(rec)) {
      paths <- c(paths, recovery = .ba_write_tsv(rec$table, file.path(outdir, "recovery.tsv")))
      summary$recovery <- rec$summary
    }
  }
  list(summary = summary, paths = paths, fit = fit, results = res, recovery = if (!is.null(truth)) rec$table)
}

# Coefficient table (term, estimate, se, df, statistic, p) for every supported fit.
.ba_coef_table <- function(fit, type, frame, p) {
  if (type == "survival") {
    if (inherits(fit, "coxme")) {
      b <- coxme::fixef(fit); se <- sqrt(diag(as.matrix(stats::vcov(fit))))
    } else {
      co <- summary(fit)$coefficients
      co <- co[!grepl("frailty\\(", rownames(co)), , drop = FALSE]      # the frailty term is not a fixed coefficient
      b <- stats::setNames(co[, "coef"], rownames(co)); se <- co[, "se(coef)"]
    }
    z <- b / se
    return(data.frame(term = names(b), estimate = unname(b), se = unname(se), df = Inf, statistic = unname(z),
                      p = 2 * stats::pnorm(-abs(unname(z))), stringsAsFactors = FALSE, row.names = NULL))
  }
  if (inherits(fit, "multinom")) {
    # one row per category (vs the reference category) and coefficient
    b <- stats::coef(fit); if (is.null(dim(b))) b <- t(as.matrix(b))
    se <- sqrt(diag(stats::vcov(fit))); se <- matrix(se, nrow = nrow(b), byrow = TRUE)
    if (nrow(b) == 1L && is.null(rownames(b))) rownames(b) <- fit$lev[2]
    out <- do.call(rbind, lapply(seq_len(nrow(b)), function(i) {
      z <- b[i, ] / se[i, ]
      data.frame(term = colnames(b), category = rownames(b)[i], estimate = unname(b[i, ]), se = unname(se[i, ]), df = Inf, statistic = unname(z),
                 p = 2 * stats::pnorm(-abs(unname(z))), stringsAsFactors = FALSE)
    }))
    return(out[out$term != "(Intercept)", , drop = FALSE])
  }
  if (inherits(fit, c("clm", "clmm", "polr"))) {
    co <- if (inherits(fit, "polr")) {
      cf <- stats::coef(summary(fit)); cf <- cf[names(stats::coef(fit)), , drop = FALSE]
      cbind(cf[, 1:2, drop = FALSE], z = cf[, 3], p = 2 * stats::pnorm(-abs(cf[, 3])))
    } else {
      cf <- summary(fit)$coefficients; cf[!rownames(cf) %in% names(fit$alpha), , drop = FALSE]
    }
    if (!length(co) || !nrow(co)) return(data.frame(term = character(), estimate = numeric(), se = numeric(), df = numeric(), statistic = numeric(), p = numeric()))   # y ~ 1
    return(data.frame(term = rownames(co), estimate = co[, 1], se = co[, 2], df = Inf, statistic = co[, 3], p = co[, 4],
                      stringsAsFactors = FALSE, row.names = NULL))
  }
  co <- summary(fit)$coefficients
  if (inherits(fit, "glmmTMB")) co <- co$cond
  pcol <- grep("^Pr\\(", colnames(co), value = TRUE)
  dfv <- if ("df" %in% colnames(co)) co[, "df"] else if (inherits(fit, "lm") && !inherits(fit, "glm")) rep(stats::df.residual(fit), nrow(co)) else
    if (inherits(fit, "glm")) rep(Inf, nrow(co)) else rep(max(1, nrow(frame) - nrow(co)), nrow(co))
  tcol <- grep("^[tz] value$", colnames(co), value = TRUE)[1]
  stat <- if (!is.na(tcol)) co[, tcol] else co[, 1] / co[, 2]
  pv <- if (length(pcol)) co[, pcol] else 2 * stats::pt(-abs(stat), dfv)
  data.frame(term = rownames(co), estimate = co[, 1], se = co[, 2], df = dfv, statistic = stat, p = pv,
             stringsAsFactors = FALSE, row.names = NULL)
}

# True treatment coefficients (Gaussian scale) under the model's parametrisation: fit the fixed
# treatment part of the model to the noise-free treatment means written in sim_truth.json.
# True cell means from Simitall's truth (intercept, treatment, numeric and interaction effects) on a grid.
.ba_true_cell_means <- function(truth, grid) {
  fv <- names(truth$treatment_effects); nv <- names(truth$numeric_effects)
  tv <- intersect(c(fv, nv), names(grid))
  mu <- rep(as.numeric(truth$intercept %||% 0), nrow(grid))
  for (v in intersect(fv, tv)) mu <- mu + unlist(truth$treatment_effects[[v]])[as.character(grid[[v]])]
  for (v in intersect(nv, tv)) mu <- mu + as.numeric(truth$numeric_effects[[v]]) * grid[[v]]
  for (it in names(truth$interaction_effects)) {
    ie <- truth$interaction_effects[[it]]
    vars <- all.vars(str2lang(it))
    if (!is.list(ie) || !all(vars %in% names(grid))) next
    key <- do.call(paste, c(lapply(vars, function(v) as.character(grid[[v]])), sep = ":"))
    v <- unlist(ie$cells)[key]
    mu <- mu + ifelse(is.na(v), 0, v)
  }
  unname(mu)
}

.ba_true_coefs <- function(truth, frame, p) {
  fv <- names(truth$treatment_effects); nv <- names(truth$numeric_effects)
  tv <- intersect(c(fv, nv), names(frame))
  if (!length(tv)) return(NULL)
  grid <- unique(frame[tv])
  mu <- rep(as.numeric(truth$intercept %||% 0), nrow(grid))
  for (v in intersect(fv, tv)) mu <- mu + unlist(truth$treatment_effects[[v]])[as.character(grid[[v]])]
  for (v in intersect(nv, tv)) mu <- mu + as.numeric(truth$numeric_effects[[v]]) * grid[[v]]
  for (it in names(truth$interaction_effects)) {
    ie <- truth$interaction_effects[[it]]
    vars <- all.vars(str2lang(it))
    if (!is.list(ie) || !all(vars %in% names(grid))) next
    key <- do.call(paste, c(lapply(vars, function(v) as.character(grid[[v]])), sep = ":"))
    v <- unlist(ie$cells)[key]
    mu <- mu + ifelse(is.na(v), 0, v)
  }
  labels <- p$fixed_labels[vapply(p$fixed_labels, function(l) all(all.vars(str2lang(l)) %in% tv), logical(1))]
  if (!length(labels)) return(NULL)
  grid$true_mean_ <- mu
  cf <- stats::coef(stats::lm(stats::reformulate(labels, response = "true_mean_"), data = grid))
  cf <- cf[names(cf) != "(Intercept)"]
  cf[abs(cf) < 1e-12] <- 0
  cf
}

# Continuous: interval coverage, bias, RMSE. Other types (the truth is on the latent Gaussian
# scale): sign agreement, detection of the non-zero effects and rank correlation.
.ba_outcome_recovery <- function(res, truth, frame, p, alpha, type) {
  tc <- .ba_true_coefs(truth, frame, p)
  if (is.null(tc)) return(NULL)
  keep <- intersect(names(tc), res$term)
  if (!length(keep)) return(NULL)
  i <- match(keep, res$term)
  tab <- data.frame(term = keep, true = as.numeric(tc[keep]), estimate = res$estimate[i], se = res$se[i], df = res$df[i], p = res$p[i])
  if (type == "continuous") {
    tab$error <- tab$estimate - tab$true
    tab$covered <- abs(tab$error) <= stats::qt(1 - alpha / 2, tab$df) * tab$se
    return(list(table = tab, summary = list(
      scale = "same as the outcome", n_terms = nrow(tab), covered = sum(tab$covered), coverage = round(mean(tab$covered), 3),
      bias = signif(mean(tab$error), 4), rmse = signif(sqrt(mean(tab$error^2)), 4),
      cor = if (nrow(tab) > 2 && stats::sd(tab$true) > 0) round(stats::cor(tab$true, tab$estimate), 4) else NA, nominal = 1 - alpha)))
  }
  # survival: a longer time (higher latent outcome) is a lower hazard
  s <- if (type == "survival") -1 else 1
  nz <- tab$true != 0
  tab$sign_agrees <- ifelse(nz, sign(s * tab$estimate) == sign(tab$true), NA)
  tab$detected <- ifelse(nz, tab$p < alpha, NA)
  list(table = tab, summary = list(
    scale = sprintf("%s fit vs Gaussian (latent) truth: compared by sign and detection%s", type, if (s < 0) "; hazard sign reversed" else ""),
    n_terms = nrow(tab), n_nonzero = sum(nz), sign_agreement = if (any(nz)) round(mean(tab$sign_agrees[nz]), 3) else NA,
    detected = if (any(nz)) round(mean(tab$detected[nz]), 3) else NA,
    false_positives = sum(!nz & tab$p < alpha),
    rank_cor = if (nrow(tab) > 2 && stats::sd(tab$true) > 0) round(stats::cor(s * tab$estimate, tab$true, method = "spearman"), 3) else NA))
}

# One outcome fit per engine id (inst/rules/engines.tsv); the chain decides which.
.ba_fit_outcome <- function(eng, f, fixed_f, frame, type, p) {
  bars <- vapply(.ba_findbars(f), function(b) paste0("(", paste(deparse(b), collapse = ""), ")"), "")
  surv_rhs <- function(extra = character()) {
    # fixed block terms -> strata(block) (survival_rules.json cox_block_terms): coefficients for many small blocks diverge
    lab <- p$fixed_labels
    if (.ba_rules("survival_rules")$cox_block_terms %in% c("strata", "frailty_or_strata")) {   # blocks still fixed here are the few-block case
      blk <- vapply(lab, function(l) { v <- all.vars(str2lang(l)); length(v) == 1L && identical(p$role[[v]], "block") }, logical(1))
      lab[blk] <- sprintf("strata(%s)", lab[blk])
    }
    # analysis.survival.strata not already in the model: stratify the Cox baseline on them (Biodes 01d41ef)
    in_model <- unique(unlist(lapply(lab, function(l) all.vars(str2lang(l)))))
    sv_strata <- setdiff(p$surv_strata, in_model)
    rhs <- c(lab, if (length(sv_strata)) sprintf("strata(%s)", sv_strata), extra); if (!length(rhs)) rhs <- "1"
    stats::as.formula(paste("survival::Surv(y_, event_) ~", paste(rhs, collapse = " + ")),
                      env = list2env(list(strata = survival::strata, cluster = survival::cluster), parent = parent.frame()))   # cox.zph re-evaluates frame there
  }
  switch(eng,
    stats_lm = stats::lm(fixed_f, data = frame),
    lmer = {
      n_fixed <- ncol(stats::model.matrix(stats::reformulate(c("1", p$fixed_labels)), data = frame))
      max_sat <- .ba_rules("handoff")$satterthwaite_max_fixed_coefficients %||% Inf
      if (n_fixed <= max_sat && requireNamespace("lmerTest", quietly = TRUE)) lmerTest::lmer(f, data = frame) else lme4::lmer(f, data = frame)
    },
    glm_binomial = stats::glm(fixed_f, data = frame, family = stats::binomial()),
    glmer_binomial = lme4::glmer(f, data = frame, family = stats::binomial()),
    glm_nb = MASS::glm.nb(fixed_f, data = frame),
    glmer_nb = lme4::glmer.nb(f, data = frame),
    clm = if (length(bars)) ordinal::clmm(f, data = frame) else ordinal::clm(fixed_f, data = frame),
    polr = MASS::polr(fixed_f, data = frame, Hess = TRUE),
    multinom = nnet::multinom(fixed_f, data = frame, trace = FALSE, Hess = TRUE),
    coxph = survival::coxph(surv_rhs(), data = frame),
    coxph_frailty = {
      g <- all.vars(.ba_findbars(f)[[1]][[3]])
      frame$frailty_g_ <- interaction(frame[g], drop = TRUE)
      survival::coxph(surv_rhs("survival::frailty(frailty_g_)"), data = frame)
    },
    coxme = {
      # the formula object goes into the call (do.call), so formula(fit) and printing show it, not "surv_rhs(bars)"
      ff <- surv_rhs(bars); fit_ <- do.call(coxme::coxme, list(formula = ff, data = quote(frame)))
      attr(fit_, "ba_formula") <- ff; fit_
    },
    glmmtmb = {
      if (!type %in% c("continuous", "binary", "count")) stop(sprintf("glmmTMB is set up for continuous, binary and count outcomes, not %s.", type), call. = FALSE)
      glmmTMB::glmmTMB(f, data = frame, family = switch(type, continuous = stats::gaussian(), binary = stats::binomial(), count = glmmTMB::nbinom2()))
    },
    stop(sprintf("No outcome fit for engine '%s'.", eng), call. = FALSE))
}

# The values must be what the outcome type says they are.
.ba_check_outcome_values <- function(y, event, type, col) {
  v <- y[!is.na(y)]
  bad <- switch(type,
    binary = if (!all(v %in% c(0, 1))) sprintf("a binary outcome must be 0/1; '%s' has %s", col, paste(utils::head(setdiff(unique(v), c(0, 1)), 3), collapse = ", ")),
    count = if (any(v < 0 | abs(v - round(v)) > 1e-8)) sprintf("a count outcome must be non-negative whole numbers; '%s' has e.g. %s", col, v[v < 0 | abs(v - round(v)) > 1e-8][1]),
    ordinal = if (length(unique(v)) < 2) sprintf("an ordinal outcome needs at least 2 categories; '%s' has %d", col, length(unique(v))),
    nominal = if (length(unique(v)) < 2) sprintf("a nominal outcome needs at least 2 categories; '%s' has %d", col, length(unique(v))),
    survival = {
      e <- event[!is.na(event)]
      if (any(v < 0)) sprintf("survival times must be >= 0; '%s' has %s", col, min(v))
      else if (!all(e %in% c(0, 1))) "the event column must be 0 (censored) or 1 (event)"
    },
    NULL)
  if (!is.null(bad)) stop("Outcome check failed: ", bad, call. = FALSE)
  invisible(TRUE)
}

# One-line formula text for the report (no deparse line breaks).
.ba_formula_text <- function(f) gsub("\\s+", " ", paste(deparse(f, width.cutoff = 500L), collapse = " "))

# Experimental unit of each metadata row: unit_id (.ba_unit_id, shared with Simitall); else Biodes 0.6.7's one row
# per subsample -- the analysis.survival.cluster column, or plot when there is a subsample column; else sample_id.
.ba_outcome_unit <- function(meta, spec) {
  if ("unit_id" %in% names(meta)) return(meta$unit_id)
  cl <- spec$analysis$survival$cluster %||% spec$analysis$cluster
  if (!is.null(cl) && cl %in% names(meta)) return(as.character(meta[[cl]]))
  if ("subsample" %in% names(meta) && "plot" %in% names(meta)) return(as.character(meta$plot))
  .ba_unit_id(meta)
}
