# Classical ANOVA on the outcome, next to the model fit (rule table inst/rules/stat_tests.tsv,
# rules inst/rules/posthoc.json):
#   repeated-measures / mixed ANOVA with sphericity corrections, or factorial ANOVA (afex::aov_car),
#   post-hoc comparisons (emmeans), generalized and partial eta squared (afex), assumption checks.
# The design picks the test. Checks are reported; they do not switch the test.

# Which ANOVA the design allows: "within", "mixed", "between", or "none" with the reason.
.ba_anova_layout <- function(frame, p) {
  subj <- character(); other <- character()
  for (r in p$random) {
    if (isTRUE(r$intercept_only) && length(r$vars) == 1L && identical(p$role[[r$vars]], "unit")) subj <- c(subj, r$vars)
    else other <- c(other, r$label)
  }
  fixed <- intersect(unique(unlist(lapply(p$fixed_labels, function(l) all.vars(str2lang(l))))), names(frame))
  num <- fixed[vapply(fixed, function(v) is.numeric(frame[[v]]), logical(1))]
  if (length(num)) return(list(layout = "none", reason = sprintf("numeric term %s (trend): the slope is the test", paste(num, collapse = ", "))))
  if (length(other)) return(list(layout = "none", reason = sprintf("random terms other than one subject (%s): see the mixed model's ANOVA table", paste(other, collapse = ", "))))
  if (length(subj) > 1L) return(list(layout = "none", reason = "more than one subject term"))
  if (!length(fixed)) return(list(layout = "none", reason = "no fixed factors"))
  if (!length(subj)) {
    n_cell <- table(interaction(frame[fixed], drop = FALSE))
    if (any(n_cell == 0)) return(list(layout = "none", reason = sprintf("%d empty cells of %s: factors are not fully crossed", sum(n_cell == 0), paste(fixed, collapse = " x "))))
    return(list(layout = "between", subject = NULL, within = character(), between = fixed))
  }
  s <- as.character(frame[[subj]])
  varies <- vapply(fixed, function(v) any(tapply(as.character(frame[[v]]), s, function(z) length(unique(z)) > 1L)), logical(1))
  within <- fixed[varies]; between <- fixed[!varies]
  if (!length(within)) return(list(layout = "none", reason = sprintf("no factor varies within %s", subj)))
  tab <- table(s, interaction(frame[within], drop = FALSE))
  miss <- rowSums(tab == 0) > 0
  if (any(miss)) return(list(layout = "none", reason = sprintf("%d of %d %s levels miss some %s cells: repeated-measures ANOVA would drop them; the mixed model keeps them",
                                                            sum(miss), length(miss), subj, paste(within, collapse = " x "))))
  list(layout = if (length(between)) "mixed" else "within", subject = subj, within = within, between = between,
       replicated = any(tab > 1))
}

.ba_stat_tests <- function(frame, p, design, goal, alpha = 0.05) {
  lay <- .ba_anova_layout(frame, p)
  if (identical(lay$layout, "none")) return(list(summary = list(status = "not run", note = lay$reason)))
  rules <- .ba_tsv("stat_tests"); row <- rules[rules$layout == lay$layout, , drop = FALSE]
  for (pk in strsplit(row$requires, ",", fixed = TRUE)[[1]]) {
    if (!requireNamespace(pk, quietly = TRUE)) return(list(summary = list(status = "not run", note = sprintf("%s not installed", pk))))
  }
  ph <- .ba_rules("posthoc")
  d <- frame
  if (is.null(lay$subject)) { d$id_ <- factor(seq_len(nrow(d))); idv <- "id_" } else idv <- lay$subject
  # between part as the design's model writes it (treatment + block stays additive); afex crosses the
  # within-subject factors with everything
  bl <- p$fixed_labels[vapply(p$fixed_labels, function(l) all(all.vars(str2lang(l)) %in% lay$between), logical(1))]
  bt <- if (length(bl)) paste(bl, collapse = " + ") else "1"
  err <- if (length(lay$within)) sprintf("Error(%s/(%s))", idv, paste(lay$within, collapse = " * ")) else sprintf("Error(%s)", idv)
  f <- stats::as.formula(sprintf("y_ ~ %s + %s", bt, err))
  fit <- suppressMessages(afex::aov_car(f, data = d, type = 3, factorize = TRUE,
                                        anova_table = list(es = c("ges", "pes"), correction = "none")))
  corr <- ph$sphericity_correction %||% "GG"
  tab <- .ba_anova_corrected(fit, corr, ph$girden_threshold %||% 0.75)
  sph <- .ba_sphericity(fit)
  # the design's model may have fewer interactions than the crossed ANOVA
  model_terms <- vapply(p$fixed_labels, function(l) paste(sort(all.vars(str2lang(l))), collapse = ":"), "")
  anova_terms <- vapply(rownames(tab), function(t) paste(sort(strsplit(t, ":", fixed = TRUE)[[1]]), collapse = ":"), "")
  extra <- rownames(tab)[grepl(":", rownames(tab), fixed = TRUE) & !anova_terms %in% model_terms]
  notes <- c(if (length(extra)) sprintf("ANOVA crosses all factors; analysis.model has no %s", paste(extra, collapse = ", ")),
             if (isTRUE(lay$replicated)) sprintf("several observations per %s cell: averaged by afex", idv))
  ph_res <- .ba_posthoc(fit, lay, p, goal, design$spec$analysis$posthoc, ph)
  chk <- .ba_anova_checks(fit, d, lay, ph$outlier_iqr %||% 3, p)
  rows <- lapply(seq_len(nrow(tab)), function(i) c(list(term = rownames(tab)[i]), as.list(tab[i, ])))
  list(summary = list(status = "run", test = row$test, layout = lay$layout, package = "afex", version = as.character(utils::packageVersion("afex")),
                      call = sprintf("afex::aov_car(%s, type = 3)", paste(deparse(f), collapse = "")),
                      subject = lay$subject, within = as.list(lay$within), between = as.list(lay$between),
                      correction = if (length(lay$within)) corr else "none", table = rows, sphericity = sph, checks = chk,
                      posthoc = ph_res$summary, reference = row$reference, tutorial = row$doc, alternative = row$alternative,
                      note = if (length(notes)) paste(notes, collapse = "; ")),
       anova_table = data.frame(term = rownames(tab), tab, row.names = NULL, check.names = FALSE),
       posthoc = ph_res$table)
}

# ANOVA table with the sphericity correction from posthoc.json: GG, HF, girden (HF when GG eps > threshold), none.
.ba_anova_corrected <- function(fit, corr, thr) {
  raw <- function(cc) {
    a <- withCallingHandlers(stats::anova(fit, correction = cc, es = c("ges", "pes")),
                             warning = function(w) if (grepl("HF eps > 1", conditionMessage(w))) invokeRestart("muffleWarning"))
    data.frame(num_df = a[["num Df"]], den_df = a[["den Df"]], MSE = a$MSE, F = a$F, ges = a$ges, pes = a$pes, p = a[["Pr(>F)"]],
               row.names = rownames(a), check.names = FALSE)
  }
  if (identical(corr, "girden")) {
    gg <- raw("GG"); hf <- raw("HF"); none <- raw("none")
    eps <- gg$num_df / none$num_df
    use_hf <- !is.na(eps) & eps > thr & eps < 1
    out <- gg; out[use_hf, ] <- hf[use_hf, ]
    out$correction <- ifelse(eps < 1, ifelse(use_hf, "HF", "GG"), "none")
    return(out)
  }
  cc <- if (corr %in% c("GG", "HF")) corr else "none"
  out <- raw(cc); none <- raw("none")
  out$correction <- ifelse(out$num_df < none$num_df - 1e-9, cc, "none")
  out
}

.ba_sphericity <- function(fit) {
  if (!inherits(fit$Anova, "Anova.mlm")) return(NULL)   # no within-subject factors
  s <- tryCatch(suppressWarnings(summary(fit$Anova, multivariate = FALSE)), error = function(e) NULL)
  if (!is.list(s) || is.null(s$sphericity.tests) || !length(s$sphericity.tests)) return(NULL)
  m <- s$sphericity.tests; a <- s$pval.adjustments
  lapply(rownames(m), function(t) list(term = t, mauchly_W = unname(m[t, "Test statistic"]), mauchly_p = unname(m[t, "p-value"]),
                                       gg_epsilon = unname(a[t, "GG eps"]), hf_epsilon = unname(a[t, "HF eps"])))
}

# Post-hoc comparisons: which families from the goal (posthoc.json families), adjustment from posthoc.json or
# analysis.posthoc. emmeans on the afex fit; multivariate model for within-subject contrasts.
.ba_posthoc <- function(fit, lay, p, goal, req, ph) {
  type <- req$type %||% ph$posthoc_type %||% "pairwise"
  if (identical(type, "none")) return(list(summary = list(status = "not run", note = "analysis.posthoc.type = none")))
  facs <- c(lay$between, lay$within)
  treat <- facs[vapply(facs, function(v) identical(p$role[[v]] %||% "treatment", "treatment"), logical(1))]
  gid <- goal$id %||% "treatment_effect"
  fam <- list()
  if (identical(gid, "trend")) return(list(summary = list(status = "not run", note = "goal trend: the slope is the test")))
  if (identical(gid, "time_course") && "time" %in% facs) {
    grp <- setdiff(treat, "time")
    fam <- if (length(grp)) list(list(specs = "time", by = grp[1])) else list(list(specs = "time", by = NULL))
  } else if (identical(gid, "interaction")) {
    gv <- intersect(goal$terms %||% treat, facs); if (length(gv) < 2) gv <- treat
    if (length(gv) >= 2) fam <- list(list(specs = gv[1], by = gv[2]), list(specs = gv[2], by = gv[1]))
  }
  if (!length(fam)) {
    fam <- lapply(treat, function(v) list(specs = v, by = NULL))
    if (length(treat) >= 2 && any(vapply(p$fixed_labels, function(l) all(treat[1:2] %in% all.vars(str2lang(l))), logical(1))))
      fam <- c(fam, list(list(specs = treat[1], by = treat[2])))
  }
  if (!length(fam)) return(list(summary = list(status = "not run", note = "no treatment factor to compare")))
  adj <- req$adjust %||% if (identical(type, "vs_control")) ph$adjust_vs_control %||% "dunnettx" else ph$adjust %||% "holm"
  emm_model <- ph$emmeans_model %||% "multivariate"
  out <- list()
  for (fm in fam) {
    sp <- if (is.null(fm$by)) stats::as.formula(paste("~", fm$specs)) else stats::as.formula(sprintf("~ %s | %s", fm$specs, fm$by))
    em <- suppressMessages(emmeans::emmeans(fit, sp, model = emm_model))
    ct <- if (identical(type, "vs_control")) {
      lev <- levels(fit$data$long[[fm$specs]]); ref <- match(req$control %||% lev[1], lev); if (is.na(ref)) ref <- 1L
      emmeans::contrast(em, "trt.vs.ctrl", ref = ref)
    } else graphics::pairs(em)
    s_adj <- as.data.frame(summary(ct, adjust = adj)); s_raw <- as.data.frame(summary(ct, adjust = "none"))
    stat <- intersect(c("t.ratio", "z.ratio"), names(s_adj))[1]
    out[[length(out) + 1L]] <- data.frame(
      family = if (is.null(fm$by)) fm$specs else sprintf("%s | %s", fm$specs, fm$by),
      by = if (is.null(fm$by)) NA_character_ else as.character(s_adj[[fm$by]]),
      contrast = as.character(s_adj$contrast), estimate = s_adj$estimate, se = s_adj$SE, df = s_adj$df,
      statistic = s_adj[[stat]], p = s_raw$p.value, p_adj = s_adj$p.value, adjust = adj, stringsAsFactors = FALSE)
  }
  tab <- do.call(rbind, out)
  list(summary = list(status = "run", package = "emmeans", version = as.character(utils::packageVersion("emmeans")), type = type,
                      adjust = adj, emmeans_model = emm_model, families = as.list(unique(tab$family)),
                      n_comparisons = nrow(tab), n_significant = sum(tab$p_adj < 0.05, na.rm = TRUE)),
       table = tab)
}

# Reported, not used to switch tests.
.ba_anova_checks <- function(fit, d, lay, iqr_k, p) {
  facs <- c(lay$between, lay$within)
  cell <- interaction(d[facs], drop = TRUE)
  ext <- unlist(tapply(d$y_, cell, function(y) { q <- stats::quantile(y, c(.25, .75), na.rm = TRUE); i <- diff(q)
    sum(y < q[1] - iqr_k * i | y > q[2] + iqr_k * i, na.rm = TRUE) }))
  r <- tryCatch(stats::residuals(fit), error = function(e) NULL)
  sw <- if (!is.null(r) && length(r) >= 3 && length(r) <= 5000 && stats::sd(r) > 0) stats::shapiro.test(r) else NULL
  lev <- NULL
  btreat <- lay$between[vapply(lay$between, function(v) identical(p$role[[v]] %||% "treatment", "treatment"), logical(1))]
  if (length(btreat) && requireNamespace("car", quietly = TRUE)) {
    bcell <- interaction(d[btreat], drop = TRUE)
    ps <- if (length(lay$within)) {
      wcell <- interaction(d[lay$within], drop = TRUE)
      vapply(split(seq_len(nrow(d)), wcell), function(i) tryCatch(car::leveneTest(d$y_[i], droplevels(bcell[i]))[1, "Pr(>F)"], error = function(e) NA_real_), 0)
    } else tryCatch(car::leveneTest(d$y_, bcell)[1, "Pr(>F)"], error = function(e) NA_real_)
    lev <- if (any(is.finite(ps))) min(ps, na.rm = TRUE) else NULL
  }
  list(extreme_outliers = sum(ext), shapiro_p = if (!is.null(sw)) unname(sw$p.value), levene_min_p = lev,
       warnings = c(if (sum(ext)) sprintf("%d extreme outlier(s) (beyond %g IQR within a cell)", sum(ext), iqr_k),
                    if (!is.null(sw) && sw$p.value < 0.05) sprintf("residuals depart from normality (Shapiro-Wilk p = %.3g)", sw$p.value),
                    if (!is.null(lev) && is.finite(lev) && lev < 0.05) sprintf("unequal variances between groups (Levene p = %.3g)", lev)))
}
