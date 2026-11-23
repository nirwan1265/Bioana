# The primary outcome: the model Biodes wrote (analysis.model), fitted with the
# family the outcome type calls for. Continuous: lm / lmer. Binary: glm /
# glmer (binomial). Count: negative binomial glm (MASS) / glmer.nb. Ordinal:
# cumulative logit, ordinal::clm / clmm (MASS::polr without random terms when
# ordinal is missing). Survival: Cox model (survival::coxph) on time + event,
# random terms as frailty-free strata are dropped with a warning.

.ba_method_outcome <- function(entry, design, outdir, alpha = 0.05, args = list()) {
  p <- design$parsed; spec <- design$spec; meta <- design$meta
  subsamples <- match.arg(args$subsamples %||% "aggregate", c("aggregate", "nest", "as_replicates"))
  col <- entry$column %||% "primary_outcome"
  type <- entry$outcome_type %||% spec$primary_outcome$type %||% "continuous"
  warn <- character()
  y <- suppressWarnings(as.numeric(meta[[col]]))
  if (!any(!is.na(y))) stop(sprintf("Column '%s' of the sample metadata is empty: nothing to fit.", col), call. = FALSE)
  meta$y_ <- y
  if (type == "survival") {
    ev <- suppressWarnings(as.numeric(meta[[entry$event_column %||% "event"]]))
    if (all(is.na(ev))) { ev <- rep(1, nrow(meta)); warn <- c(warn, "no event column: every time treated as an event") }
    meta$event_ <- ev
  }
  unit <- .ba_unit_id(meta)
  cell <- do.call(paste, c(list(unit), lapply(p$vars, function(v) meta[[p$col[[v]]]]), sep = "\r"))
  has_subs <- anyDuplicated(cell[!is.na(unit)]) > 0
  if (has_subs && subsamples == "as_replicates") {
    stop(sprintf("These rows are subsamples of %d %s units. Biodes set the experimental unit to '%s', so subsamples are not independent replicates. Use subsamples = 'aggregate' or 'nest'.",
                 length(unique(unit)), spec$experimental_unit %||% "experimental", spec$experimental_unit %||% "unit"), call. = FALSE)
  }
  dat <- meta[!is.na(y), , drop = FALSE]
  dcell <- cell[!is.na(y)]
  if (has_subs && (subsamples == "aggregate" || type != "continuous")) {
    if (type != "continuous" && subsamples == "nest") warn <- c(warn, "subsamples averaged per unit: nesting is fitted for continuous outcomes only")
    first <- dat[!duplicated(dcell), , drop = FALSE]
    first$y_ <- as.numeric(tapply(dat$y_, dcell, mean)[dcell[!duplicated(dcell)]])
    if (type %in% c("binary", "count", "ordinal")) first$y_ <- round(first$y_)
    dat <- first
  }
  f <- p$formula; f[[2]] <- as.name("y_")
  if (has_subs && subsamples == "nest" && type == "continuous") {
    dat$unit_id_ <- dcell
    f <- stats::update(f, . ~ . + (1 | unit_id_))
  }
  frame <- .ba_frame(dat, p)
  if (type == "survival") frame$event_ <- dat$event_
  has_random <- length(.ba_findbars(f)) > 0
  fixed_f <- .ba_nobars(f)
  ch <- .ba_engine_chain(paste0("outcome:", type), design, NULL, args$engine)
  if (!length(ch$chain)) stop("No outcome engine installed: ", paste(names(ch$skipped), unlist(ch$skipped), collapse = "; "), call. = FALSE)
  eng <- ch$chain[1]
  fit <- if (identical(eng, "glmmtmb")) {
    if (!type %in% c("continuous", "binary", "count")) stop(sprintf("glmmTMB is set up for continuous, binary and count outcomes, not %s.", type), call. = FALSE)
    glmmTMB::glmmTMB(f, data = frame, family = switch(type, continuous = stats::gaussian(), binary = stats::binomial(), count = glmmTMB::nbinom2()))
  } else switch(type,
    continuous = {
      n_fixed <- ncol(stats::model.matrix(stats::reformulate(c("1", p$fixed_labels)), data = frame))
      max_sat <- .ba_rules("handoff")$satterthwaite_max_fixed_coefficients %||% Inf
      if (has_random && n_fixed > max_sat) warn <- c(warn, sprintf("%d fixed coefficients: lme4 without Satterthwaite df", n_fixed))
      if (has_random) {
        if (n_fixed <= max_sat && requireNamespace("lmerTest", quietly = TRUE)) lmerTest::lmer(f, data = frame) else lme4::lmer(f, data = frame)
      } else stats::lm(f, data = frame)
    },
    binary = if (has_random) lme4::glmer(f, data = frame, family = stats::binomial()) else stats::glm(fixed_f, data = frame, family = stats::binomial()),
    count = {
      if (!requireNamespace("MASS", quietly = TRUE)) stop("count outcomes need MASS", call. = FALSE)
      if (has_random) lme4::glmer.nb(f, data = frame) else MASS::glm.nb(fixed_f, data = frame)
    },
    ordinal = {
      frame$y_ <- factor(frame$y_, ordered = TRUE)
      if (requireNamespace("ordinal", quietly = TRUE)) {
        if (has_random) ordinal::clmm(f, data = frame) else ordinal::clm(fixed_f, data = frame)
      } else {
        if (has_random) warn <- c(warn, "random terms dropped: install 'ordinal' to fit them (clmm)")
        MASS::polr(fixed_f, data = frame, Hess = TRUE)
      }
    },
    survival = {
      if (!requireNamespace("survival", quietly = TRUE)) stop("survival outcomes need the survival package", call. = FALSE)
      if (has_random) warn <- c(warn, "random terms dropped for the Cox model (no frailty fitted)")
      sf <- stats::as.formula(paste("survival::Surv(y_, event_) ~", paste(p$fixed_labels, collapse = " + ")))
      survival::coxph(sf, data = frame)
    },
    stop(sprintf("No outcome analysis for type '%s'.", type), call. = FALSE))
  if (!identical(eng, "glmmtmb")) {
    # the engine the rules name for this type and design is the one the switch above fitted
    eng <- .ba_selection_rows(paste0("outcome:", type), if (has_random) "random" else "fixed")$engine[1]
    if (inherits(fit, "polr")) warn <- c(warn, "ordinal not installed: MASS::polr fitted instead of ordinal::clm")
  }
  if (has_random && inherits(fit, "merMod") && lme4::isSingular(fit)) warn <- c(warn, "singular fit: a random-effect variance was estimated as ~0")
  res <- .ba_coef_table(fit, type, frame, p)
  res$role <- vapply(res$term, function(t) .ba_coef_role(t, p), "")
  res$padj <- NA_real_
  test <- res$role == "treatment"
  res$padj[test] <- stats::p.adjust(res$p[test], "BH")
  paths <- c(results = .ba_write_tsv(res, file.path(outdir, "results.tsv")))
  aov <- if (type == "continuous") tryCatch(as.data.frame(stats::anova(fit)), error = function(e) NULL) else NULL
  if (!is.null(aov) && has_random && !"Pr(>F)" %in% names(aov) && "F value" %in% names(aov)) {
    # lme4 without lmerTest gives F only: denominator df = rows - fixed coefficients
    aov$DenDF <- max(1, nrow(frame) - ncol(stats::model.matrix(stats::reformulate(c("1", p$fixed_labels)), data = frame)))
    aov[["Pr(>F)"]] <- stats::pf(aov[["F value"]], aov$npar, aov$DenDF, lower.tail = FALSE)
  }
  summary <- list(data_type = "outcome", outcome_type = type, response = spec$primary_outcome$name %||% col,
                  method = paste(class(fit)[1], if (type == "survival") "(Cox proportional hazards)" else ""),
                  model_fitted = if (type == "survival") paste(deparse(stats::formula(fit)), collapse = "") else paste(deparse(stats::formula(fit)), collapse = ""),
                  n_rows = nrow(meta), n_units_fitted = nrow(frame), subsamples = if (has_subs) subsamples else "none",
                  dropped_samples = meta$sample_id[is.na(y)],
                  engine = .ba_provenance(eng, design, tryCatch(paste(deparse(stats::getCall(fit)), collapse = " "), error = function(e) class(fit)[1]),
                                          sprintf("%s outcome '%s', %d units", type, col, nrow(frame)), nrow(frame), ch,
                                          formula = paste(deparse(f), collapse = "")),
                  secondary_engine = .ba_secondary_engine(ch$chain, eng),
                  anova = if (is.null(aov)) NULL else lapply(split(cbind(term = rownames(aov), aov), seq_len(nrow(aov))), as.list),
                  warnings = warn)
  truth <- .ba_entry_truth(entry, design$dir)
  rec <- NULL
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
    co <- summary(fit)$coefficients
    return(data.frame(term = rownames(co), estimate = co[, "coef"], se = co[, "se(coef)"], df = Inf,
                      statistic = co[, "z"], p = co[, "Pr(>|z|)"], stringsAsFactors = FALSE, row.names = NULL))
  }
  if (inherits(fit, c("clm", "clmm", "polr"))) {
    co <- if (inherits(fit, "polr")) {
      cf <- stats::coef(summary(fit)); cf <- cf[names(stats::coef(fit)), , drop = FALSE]
      cbind(cf[, 1:2, drop = FALSE], z = cf[, 3], p = 2 * stats::pnorm(-abs(cf[, 3])))
    } else {
      cf <- summary(fit)$coefficients; cf[!rownames(cf) %in% names(fit$alpha), , drop = FALSE]
    }
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
