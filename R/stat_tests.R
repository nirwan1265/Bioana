# Classical tests on the outcome, next to the model fit (rule table inst/rules/stat_tests.tsv,
# rules inst/rules/stat_rules.json): t-tests, Wilcoxon, sign test, Kruskal-Wallis, Friedman, Welch ANOVA,
# proportion and Fisher tests (stats, rstatix), repeated-measures / mixed / factorial ANOVA (afex),
# post-hoc comparisons (emmeans, rstatix, stats), effect sizes, assumption checks, descriptives.
# The design and the outcome type pick the test (first auto row that matches); more run when asked.
# Checks are reported; they do not switch the test.

# Which ANOVA the design allows: "within", "mixed", "between", or "none" with the reason.
.ba_anova_layout <- function(frame, p) {
  subj <- character(); rblock <- character(); other <- character()
  for (r in p$random) {
    one <- isTRUE(r$intercept_only) && length(r$vars) == 1L
    if (one && identical(p$role[[r$vars]], "unit")) subj <- c(subj, r$vars)
    else if (one && identical(p$role[[r$vars]], "block")) rblock <- c(rblock, r$vars)   # random block intercept: a block
    else other <- c(other, r$label)
  }
  fixed <- intersect(unique(c(unlist(lapply(p$fixed_labels, function(l) all.vars(str2lang(l)))), rblock)), names(frame))
  num <- fixed[vapply(fixed, function(v) is.numeric(frame[[v]]), logical(1))]
  if (length(num)) return(list(layout = "none", reason = sprintf("numeric term %s (trend): the slope is the test", paste(num, collapse = ", "))))
  if (length(other)) return(list(layout = "none", reason = sprintf("random terms other than one subject or block (%s): see the mixed model's ANOVA table", paste(other, collapse = ", "))))
  if (length(subj) > 1L) return(list(layout = "none", reason = "more than one subject term"))
  if (!length(fixed)) return(list(layout = "none", reason = "no fixed factors"))
  role <- function(v) p$role[[v]] %||% "treatment"
  out <- if (!length(subj)) {
    n_cell <- table(interaction(frame[fixed], drop = FALSE))
    if (any(n_cell == 0)) return(list(layout = "none", reason = sprintf("%d empty cells of %s: factors are not fully crossed", sum(n_cell == 0), paste(fixed, collapse = " x "))))
    list(layout = "between", subject = NULL, within = character(), between = fixed)
  } else {
    s <- as.character(frame[[subj]])
    varies <- vapply(fixed, function(v) any(tapply(as.character(frame[[v]]), s, function(z) length(unique(z)) > 1L)), logical(1))
    within <- fixed[varies]; between <- fixed[!varies]
    if (!length(within)) return(list(layout = "none", reason = sprintf("no factor varies within %s", subj)))
    tab <- table(s, interaction(frame[within], drop = FALSE))
    miss <- rowSums(tab == 0) > 0
    if (any(miss)) return(list(layout = "none", reason = sprintf("%d of %d %s levels miss some %s cells: repeated-measures ANOVA would drop them; the mixed model keeps them",
                                                              sum(miss), length(miss), subj, paste(within, collapse = " x "))))
    list(layout = if (length(between)) "mixed" else "within", subject = subj, within = within, between = between, replicated = any(tab > 1))
  }
  facs <- c(out$between, out$within)
  out$treat <- facs[vapply(facs, function(v) role(v) == "treatment", logical(1))]
  out$blocks <- out$between[vapply(out$between, function(v) role(v) == "block", logical(1))]
  out$random_blocks <- rblock
  out$blocking <- if (!is.null(out$subject)) "subject" else if (!length(out$blocks)) "none" else if (length(out$blocks) == 1L && length(out$treat) == 1L &&
    all(table(frame[[out$treat]], frame[[out$blocks]]) == 1L)) "block" else "blocks"
  out
}

# Does a stat_tests.tsv row fit this outcome and design?
.ba_stat_matches <- function(row, lay, type, frame) {
  sp <- function(x) trimws(strsplit(x, ",", fixed = TRUE)[[1]])
  if (!paste0("outcome:", type) %in% sp(row$applies_to) || !lay$layout %in% sp(row$layout)) return(FALSE)
  nf <- length(lay$treat)
  if (!(row$factors == "any" || (row$factors == "1" && nf == 1L) || (row$factors == "2+" && nf >= 2L))) return(FALSE)
  if (row$levels != "any") {
    if (nf != 1L) return(FALSE)
    nl <- length(unique(stats::na.omit(frame[[lay$treat]])))
    if (!((row$levels == "2" && nl == 2L) || (row$levels == "3+" && nl >= 3L))) return(FALSE)
  }
  row$blocking == "any" || lay$blocking %in% sp(row$blocking)
}

.ba_stat_tests <- function(frame, p, design, goal, alpha = 0.05, type = "continuous", args = list()) {
  lay <- .ba_anova_layout(frame, p)
  sr <- .ba_rules("stat_rules")
  desc <- .ba_descriptives(frame, p, type, sr)
  if (identical(lay$layout, "none")) return(list(summary = list(status = "not run", note = lay$reason), descriptives = desc))
  rules <- .ba_tsv("stat_tests")
  fits <- rules[vapply(seq_len(nrow(rules)), function(i) .ba_stat_matches(rules[i, ], lay, type, frame), logical(1)), , drop = FALSE]
  fam <- sr$default_family[[type]] %||% "parametric"
  ids <- utils::head(fits$test[fits$auto == "yes" & fits$family == fam], 1)
  notes <- character()
  for (r in unique(unlist(c(design$spec$analysis$stat_tests, args$stat_tests)))) {
    if (r %in% c("parametric", "nonparametric")) {
      m <- fits$test[fits$auto == "yes" & fits$family == r]
      if (length(m)) ids <- c(ids, m[1]) else notes <- c(notes, sprintf("no %s test fits this design", r))
    } else if (r %in% fits$test) ids <- c(ids, r)
    else notes <- c(notes, sprintf("requested test '%s' %s", r, if (r %in% rules$test) "does not fit this outcome and design" else "is not in stat_tests.tsv"))
  }
  ids <- unique(ids)
  if (!length(ids)) return(list(summary = list(status = "not run", note = paste(c("no test in stat_tests.tsv fits this outcome and design", notes), collapse = "; ")), descriptives = desc))
  res <- lapply(ids, function(id) {
    row <- rules[rules$test == id, , drop = FALSE]
    miss <- Filter(function(pk) !requireNamespace(pk, quietly = TRUE), strsplit(row$requires, ",", fixed = TRUE)[[1]])
    if (length(miss)) return(list(summary = list(status = "not run", test = id, note = sprintf("%s not installed", paste(miss, collapse = ", ")))))
    tryCatch(if (startsWith(id, "anova_")) .ba_st_anova(row, frame, lay, p, design, goal, sr) else .ba_st_simple(row, frame, lay, p, design, sr, rules),
             error = function(e) list(summary = list(status = "error", test = id, note = conditionMessage(e))))
  })
  main <- res[[1]]$summary
  if (length(notes)) main$note <- paste(c(main$note, notes), collapse = "; ")
  if (length(res) > 1L) main$also <- lapply(res[-1], `[[`, "summary")
  bind <- function(k) { x <- Filter(Negate(is.null), lapply(res, function(r) if (!is.null(r[[k]])) cbind(test = r$summary$test, r[[k]]))); if (length(x)) do.call(rbind, x) }
  at <- Filter(Negate(is.null), lapply(res, `[[`, "anova_table"))
  list(summary = main, anova_table = if (length(at)) at[[1]],
       posthoc = bind("posthoc"), stat_table = bind("flat"), descriptives = desc)
}

.ba_st_anova <- function(row, frame, lay, p, design, goal, ph) {
  d <- frame
  if (is.null(lay$subject)) { d$id_ <- factor(seq_len(nrow(d))); idv <- "id_" } else idv <- lay$subject
  # between part as the design's model writes it (treatment + block stays additive); afex crosses the
  # within-subject factors with everything
  bl <- p$fixed_labels[vapply(p$fixed_labels, function(l) all(all.vars(str2lang(l)) %in% lay$between), logical(1))]
  bl <- c(bl, lay$random_blocks)
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
             if (length(lay$random_blocks)) sprintf("random block intercept %s entered as an additive block (same F for balanced data)", paste(lay$random_blocks, collapse = ", ")),
             if (isTRUE(lay$replicated)) sprintf("several observations per %s cell: averaged by afex", idv))
  orig <- lapply(stats::setNames(c(lay$between, lay$within), c(lay$between, lay$within)), function(v) levels(droplevels(factor(frame[[v]]))))
  ph_res <- .ba_posthoc(fit, lay, p, goal, design$spec$analysis$posthoc, ph, orig)
  chk <- .ba_anova_checks(fit, d, lay, ph$outlier_iqr %||% 3, p)
  rows <- lapply(seq_len(nrow(tab)), function(i) c(list(term = rownames(tab)[i]), as.list(tab[i, ])))
  list(summary = list(status = "run", test = row$test, layout = lay$layout, package = "afex", version = as.character(utils::packageVersion("afex")),
                      call = sprintf("afex::aov_car(%s, type = 3)", paste(deparse(f), collapse = "")),
                      subject = lay$subject, within = as.list(lay$within), between = as.list(lay$between),
                      correction = if (length(lay$within)) corr else "none", table = rows, sphericity = sph, checks = chk,
                      posthoc = ph_res$summary, reference = row$reference, tutorial = row$doc, alternative = row$alternative,
                      note = if (length(notes)) paste(notes, collapse = "; ")),
       anova_table = data.frame(term = rownames(tab), tab, row.names = NULL, check.names = FALSE),
       posthoc = ph_res$table,
       flat = data.frame(term = rownames(tab), contrast = NA_character_, statistic_name = "F", statistic = tab$F, df1 = tab$num_df, df2 = tab$den_df,
                         p = tab$p, estimate = NA_real_, ci_low = NA_real_, ci_high = NA_real_, effect_name = "ges", effect = tab$ges, n = nrow(frame),
                         stringsAsFactors = FALSE))
}

# ANOVA table with the sphericity correction from stat_rules.json: GG, HF, girden (HF when GG eps > threshold), none.
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

# Post-hoc comparisons: which families from the goal (stat_rules.json families), adjustment from stat_rules.json or
# analysis.posthoc. emmeans on the afex fit; multivariate model for within-subject contrasts.
.ba_posthoc <- function(fit, lay, p, goal, req, ph, orig = list()) {
  type <- req$type %||% ph$posthoc_type %||% "pairwise"
  if (identical(type, "none")) return(list(summary = list(status = "not run", note = "analysis.posthoc.type = none")))
  facs <- c(lay$between, lay$within)
  treat <- facs[vapply(facs, function(v) identical(p$role[[v]] %||% "treatment", "treatment"), logical(1))]
  gid <- goal$id %||% "treatment_effect"
  fam <- list(); pnote <- character()
  if (identical(gid, "trend")) return(list(summary = list(status = "not run", note = "goal trend: the slope is the test")))
  fv <- .ba_posthoc_factor(req, p, treat)
  if (is.character(fv) && nzchar(fv) && !fv %in% treat) { pnote <- c(pnote, fv); fv <- NULL }
  if (!is.null(fv) && fv %in% treat) {
    # Biodes analysis.posthoc.factor: that factor, and within each level of a treatment factor it interacts with
    fam <- list(list(specs = fv, by = NULL))
    for (w in setdiff(treat, fv)) if (any(vapply(p$fixed_labels, function(l) all(c(fv, w) %in% all.vars(str2lang(l))), logical(1)))) fam <- c(fam, list(list(specs = fv, by = w)))
  } else if (identical(gid, "time_course") && "time" %in% facs) {
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
  out <- list(); ctrl <- list()
  for (fm in fam) {
    sp <- if (is.null(fm$by)) stats::as.formula(paste("~", fm$specs)) else stats::as.formula(sprintf("~ %s | %s", fm$specs, fm$by))
    em <- suppressMessages(emmeans::emmeans(fit, sp, model = emm_model))
    # afex passes within-subject level labels through make.names ("week 1" -> "week.1"): back to the design's labels
    lv <- levels(em)
    for (v in intersect(names(lv), names(orig))) if (length(orig[[v]]) == length(lv[[v]]) && all(make.names(orig[[v]]) == lv[[v]] | orig[[v]] == lv[[v]])) lv[[v]] <- orig[[v]]
    byv <- em@misc$by.vars
    levels(em) <- lv                                   # drops the by grouping: set it again
    if (length(byv)) em <- stats::update(em, by = byv)
    ct <- if (identical(type, "vs_control")) {
      lev <- levels(em)[[fm$specs]]; ref <- match(req$control %||% lev[1], lev)
      if (is.na(ref)) { pnote <- c(pnote, sprintf("control '%s' is not a level of %s: first level '%s' used", req$control, fm$specs, lev[1])); ref <- 1L }
      ctrl[[fm$specs]] <- lev[ref]
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
                      adjust = adj, emmeans_model = emm_model, families = as.list(unique(tab$family)), control = if (length(ctrl)) ctrl,
                      specs = lapply(fam, function(f) list(specs = f$specs, by = f$by)), requested_by = if (!is.null(req$type)) "analysis.posthoc",
                      n_comparisons = nrow(tab), n_significant = sum(tab$p_adj < 0.05, na.rm = TRUE),
                      note = if (length(pnote)) paste(unique(pnote), collapse = "; ")),
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

# ---- tests on one treatment factor (two groups, several groups, paired, blocked) and binary outcomes ----

.ba_st_simple <- function(row, frame, lay, p, design, sr, rules) {
  id <- row$test; tv <- lay$treat[1]
  d <- frame; d$g_ <- droplevels(factor(d[[tv]])); lev <- levels(d$g_)
  d$yn_ <- if (is.factor(d$y_)) as.integer(d$y_) else as.numeric(d$y_)
  d <- d[!is.na(d$yn_) & !is.na(d$g_), , drop = FALSE]
  req <- design$spec$analysis$posthoc
  adj <- req$adjust %||% sr$adjust_nonparametric %||% sr$adjust %||% "holm"
  notes <- character(); ph <- NULL; chk <- NULL
  quiet <- function(expr) withCallingHandlers(expr, warning = function(w) {
    if (grepl("exact p-value|ties|zeroes|Chi-squared approximation", conditionMessage(w))) { notes <<- unique(c(notes, "ties or zeros: normal approximation")); invokeRestart("muffleWarning") } })
  blk <- lay$subject %||% lay$blocks[1]
  wide <- function() {                       # one value per subject (or block) and level, matched
    w <- tapply(d$yn_, list(as.character(d[[blk]]), d$g_), mean); w[stats::complete.cases(w), , drop = FALSE]
  }
  long_of <- function(w) data.frame(id_ = factor(rep(rownames(w), ncol(w))), g_ = factor(rep(colnames(w), each = nrow(w)), colnames(w)), yn_ = c(w))
  fr <- function(term = tv, contrast = NA_character_, stat_name, stat, df1 = NA_real_, df2 = NA_real_, pv, est = NA_real_, lo = NA_real_, hi = NA_real_,
                 eff_name = NA_character_, eff = NA_real_, n = nrow(d))
    data.frame(term = term, contrast = contrast, statistic_name = stat_name, statistic = unname(stat), df1 = unname(df1), df2 = unname(df2), p = unname(pv),
               estimate = unname(est), ci_low = unname(lo), ci_high = unname(hi), effect_name = eff_name, effect = unname(eff), n = n, stringsAsFactors = FALSE)
  ctr <- sprintf("%s - %s", lev[1], lev[2])
  extra <- list()
  if (id %in% c("t_test_welch", "t_test_student")) {
    ve <- id == "t_test_student"
    tt <- stats::t.test(yn_ ~ g_, data = d, var.equal = ve)
    cd <- rstatix::cohens_d(d, yn_ ~ g_, var.equal = ve); hg <- rstatix::cohens_d(d, yn_ ~ g_, var.equal = ve, hedges.correction = TRUE)
    flat <- fr(contrast = ctr, stat_name = "t", stat = tt$statistic, df1 = tt$parameter, pv = tt$p.value, est = tt$estimate[1] - tt$estimate[2],
               lo = tt$conf.int[1], hi = tt$conf.int[2], eff_name = "Cohen's d", eff = cd$effsize)
    extra <- list(hedges_g = hg$effsize, magnitude = as.character(cd$magnitude))
    chk <- .ba_group_checks(d, sr)
  } else if (id == "t_test_paired") {
    w <- wide(); tt <- stats::t.test(w[, 1], w[, 2], paired = TRUE)
    lw <- long_of(w); cd <- rstatix::cohens_d(lw, yn_ ~ g_, paired = TRUE); hg <- rstatix::cohens_d(lw, yn_ ~ g_, paired = TRUE, hedges.correction = TRUE)
    flat <- fr(contrast = ctr, stat_name = "t", stat = tt$statistic, df1 = tt$parameter, pv = tt$p.value, est = tt$estimate, lo = tt$conf.int[1], hi = tt$conf.int[2],
               eff_name = "Cohen's d (differences)", eff = cd$effsize, n = nrow(w))
    extra <- list(hedges_g = hg$effsize, magnitude = as.character(cd$magnitude))
    chk <- .ba_diff_checks(w[, 1] - w[, 2], sr, "differences")
  } else if (id == "wilcoxon_rank_sum") {
    wt <- quiet(stats::wilcox.test(yn_ ~ g_, data = d)); es <- quiet(rstatix::wilcox_effsize(d, yn_ ~ g_))
    flat <- fr(contrast = ctr, stat_name = "W", stat = wt$statistic, pv = wt$p.value, eff_name = "r", eff = es$effsize)
    extra <- list(magnitude = as.character(es$magnitude)); chk <- .ba_group_checks(d, sr, normality = FALSE)
  } else if (id == "wilcoxon_signed_rank") {
    w <- wide(); wt <- quiet(stats::wilcox.test(w[, 1], w[, 2], paired = TRUE)); es <- quiet(rstatix::wilcox_effsize(long_of(w), yn_ ~ g_, paired = TRUE))
    flat <- fr(contrast = ctr, stat_name = "V", stat = wt$statistic, pv = wt$p.value, eff_name = "r", eff = es$effsize, n = nrow(w))
    extra <- list(magnitude = as.character(es$magnitude)); chk <- .ba_diff_checks(w[, 1] - w[, 2], sr, "differences", normality = FALSE)
  } else if (id == "sign_test") {
    w <- wide(); dif <- w[, 1] - w[, 2]; nz <- dif[dif != 0]
    bt <- stats::binom.test(sum(nz > 0), length(nz), 0.5)
    flat <- fr(contrast = ctr, stat_name = "positive differences", stat = sum(nz > 0), pv = bt$p.value, est = stats::median(dif), n = length(nz))
    if (sum(dif == 0)) notes <- c(notes, sprintf("%d tied pairs dropped", sum(dif == 0)))
  } else if (id == "kruskal_wallis") {
    kt <- stats::kruskal.test(yn_ ~ g_, data = d); es <- rstatix::kruskal_effsize(d, yn_ ~ g_)
    flat <- fr(stat_name = "H", stat = kt$statistic, df1 = kt$parameter, pv = kt$p.value, eff_name = "eta2[H]", eff = es$effsize)
    extra <- list(magnitude = as.character(es$magnitude)); chk <- .ba_group_checks(d, sr, normality = FALSE)
    dn <- rstatix::dunn_test(d, yn_ ~ g_, p.adjust.method = adj)
    ph <- data.frame(family = tv, by = NA_character_, contrast = sprintf("%s - %s", dn$group1, dn$group2), estimate = if (!is.null(dn[["estimate"]])) dn[["estimate"]] else NA_real_, se = NA_real_, df = NA_real_,
                     statistic = dn$statistic, p = dn$p, p_adj = dn$p.adj, adjust = adj, stringsAsFactors = FALSE)
  } else if (id == "friedman") {
    w <- wide(); ft <- stats::friedman.test(w); lw <- long_of(w)
    es <- rstatix::friedman_effsize(lw, yn_ ~ g_ | id_)
    flat <- fr(stat_name = "Friedman chi-squared", stat = ft$statistic, df1 = ft$parameter, pv = ft$p.value, eff_name = "Kendall's W", eff = es$effsize, n = nrow(w))
    extra <- list(magnitude = as.character(es$magnitude), blocked_by = blk)
    if (ncol(w) >= 3) {
      pa <- quiet(stats::pairwise.wilcox.test(lw$yn_, lw$g_, paired = TRUE, p.adjust.method = adj))
      pr <- quiet(stats::pairwise.wilcox.test(lw$yn_, lw$g_, paired = TRUE, p.adjust.method = "none"))
      ph <- .ba_pairwise_table(pa, pr, tv, adj)
    }
  } else if (id == "welch_anova") {
    ow <- stats::oneway.test(yn_ ~ g_, data = d, var.equal = FALSE); gh <- rstatix::games_howell_test(d, yn_ ~ g_)
    flat <- fr(stat_name = "F", stat = ow$statistic, df1 = ow$parameter[1], df2 = ow$parameter[2], pv = ow$p.value)
    chk <- .ba_group_checks(d, sr)
    ph <- data.frame(family = tv, by = NA_character_, contrast = sprintf("%s - %s", gh$group2, gh$group1), estimate = gh$estimate, se = NA_real_, df = NA_real_,
                     statistic = NA_real_, p = NA_real_, p_adj = gh$p.adj, adjust = "games-howell", stringsAsFactors = FALSE)
  } else if (id == "prop_test") {
    x <- tapply(d$yn_, d$g_, sum); n <- tapply(d$yn_, d$g_, length)
    e <- outer(n, c(sum(x), sum(n - x)) / sum(n))
    if (min(e) < (sr$min_expected %||% 5)) {
      r2 <- .ba_st_simple(rules[rules$test == "fisher_exact", , drop = FALSE], frame, lay, p, design, sr, rules)
      r2$summary$note <- paste(c(sprintf("smallest expected count %.2g < %g: Fisher's exact test instead of the proportion test", min(e), sr$min_expected %||% 5), r2$summary$note), collapse = "; ")
      return(r2)
    }
    pt <- stats::prop.test(x, n)
    if (length(x) == 2L) {
      h <- effectsize::cohens_h(t(matrix(c(x, n - x), ncol = 2)))   # effectsize: groups in columns
      flat <- fr(contrast = ctr, stat_name = "X-squared", stat = pt$statistic, df1 = pt$parameter, pv = pt$p.value, est = pt$estimate[1] - pt$estimate[2],
                 lo = pt$conf.int[1], hi = pt$conf.int[2], eff_name = "Cohen's h", eff = h$Cohens_h)
    } else {
      flat <- fr(stat_name = "X-squared", stat = pt$statistic, df1 = pt$parameter, pv = pt$p.value)
      pa <- stats::pairwise.prop.test(x, n, p.adjust.method = adj); pr <- stats::pairwise.prop.test(x, n, p.adjust.method = "none")
      ph <- .ba_pairwise_table(pa, pr, tv, adj)
    }
    extra <- list(proportions = as.list(stats::setNames(as.numeric(x / n), names(x))), yates = length(x) == 2L)
  } else if (id == "fisher_exact") {
    tab <- table(d$g_, factor(d$yn_, levels = c(0, 1)))
    ft <- stats::fisher.test(tab)
    flat <- if (nrow(tab) == 2L) fr(contrast = ctr, stat_name = "odds ratio", stat = ft$estimate, pv = ft$p.value, est = ft$estimate, lo = ft$conf.int[1], hi = ft$conf.int[2],
                                    eff_name = "odds ratio", eff = ft$estimate) else fr(stat_name = "exact", stat = NA_real_, pv = ft$p.value)
    if (nrow(tab) > 2L) {
      pf <- rstatix::pairwise_fisher_test(tab, p.adjust.method = adj)
      ph <- data.frame(family = tv, by = NA_character_, contrast = sprintf("%s - %s", pf$group1, pf$group2), estimate = NA_real_, se = NA_real_, df = NA_real_,
                       statistic = NA_real_, p = pf$p, p_adj = pf$p.adj, adjust = adj, stringsAsFactors = FALSE)
    }
  } else stop(sprintf("no adapter for stat test '%s'", id))
  if (!is.null(ph)) {
    f <- .ba_posthoc_select(ph, req, lev, adj, sr)
    ph <- f$table; if (length(f$note)) notes <- c(notes, f$note)
  }
  pk <- strsplit(row$requires, ",", fixed = TRUE)[[1]][1]
  list(summary = c(list(status = "run", test = id, name = row$model, layout = lay$layout, package = pk, version = as.character(utils::packageVersion(pk)),
                        call = row$calls, factor = tv, levels = as.list(lev), table = lapply(seq_len(nrow(flat)), function(i) as.list(flat[i, ])),
                        checks = chk, posthoc = if (!is.null(ph)) list(status = "run", adjust = adj, n_comparisons = nrow(ph), n_significant = sum(ph$p_adj < 0.05, na.rm = TRUE),
                                                                        families = list(tv)),
                        reference = row$reference, tutorial = row$doc, alternative = row$alternative,
                        note = if (length(notes)) paste(notes, collapse = "; ")), extra),
       posthoc = ph, flat = flat)
}

# stats::pairwise.* result (adjusted and raw) as a post-hoc table.
.ba_pairwise_table <- function(pa, pr, tv, adj) {
  m <- pa$p.value; r <- pr$p.value
  ij <- which(!is.na(m), arr.ind = TRUE)
  data.frame(family = tv, by = NA_character_, contrast = sprintf("%s - %s", colnames(m)[ij[, 2]], rownames(m)[ij[, 1]]), estimate = NA_real_, se = NA_real_,
             df = NA_real_, statistic = NA_real_, p = r[ij], p_adj = m[ij], adjust = adj, stringsAsFactors = FALSE)
}

.ba_skew <- function(x) { x <- x[!is.na(x)]; m <- mean(x); s2 <- mean((x - m)^2); if (s2 > 0) mean((x - m)^3) / s2^1.5 else NA_real_ }

# Reported per group, never used to switch tests: extreme outliers, Shapiro-Wilk, Levene (median centre), Fligner-Killeen.
.ba_group_checks <- function(d, sr, normality = TRUE) {
  k <- sr$outlier_iqr %||% 3
  ext <- sum(unlist(tapply(d$yn_, d$g_, function(y) { q <- stats::quantile(y, c(.25, .75)); i <- diff(q); sum(y < q[1] - k * i | y > q[2] + k * i) })))
  sw <- if (normality) tapply(d$yn_, d$g_, function(y) if (length(y) >= 3 && length(y) <= 5000 && stats::sd(y) > 0) stats::shapiro.test(y)$p.value else NA_real_)
  lev <- if (requireNamespace("car", quietly = TRUE)) tryCatch(car::leveneTest(d$yn_, d$g_)[1, "Pr(>F)"], error = function(e) NA_real_)
  fk <- tryCatch(stats::fligner.test(d$yn_, d$g_)$p.value, error = function(e) NA_real_)
  list(extreme_outliers = ext, shapiro_p = if (!is.null(sw)) as.list(sw), levene_p = lev, fligner_p = fk,
       warnings = c(if (ext) sprintf("%d extreme outlier(s) (beyond %g IQR within a group)", ext, k),
                    if (!is.null(sw) && any(sw < 0.05, na.rm = TRUE)) sprintf("Shapiro-Wilk p < 0.05 in %s", paste(names(sw)[which(sw < 0.05)], collapse = ", ")),
                    if (isTRUE(lev < 0.05)) sprintf("unequal variances (Levene p = %.3g)", lev)))
}

.ba_diff_checks <- function(dif, sr, what, normality = TRUE) {
  k <- sr$outlier_iqr %||% 3; q <- stats::quantile(dif, c(.25, .75)); i <- diff(q)
  ext <- sum(dif < q[1] - k * i | dif > q[2] + k * i)
  sw <- if (normality && length(dif) >= 3 && stats::sd(dif) > 0) stats::shapiro.test(dif)$p.value
  sk <- .ba_skew(dif)
  list(extreme_outliers = ext, shapiro_p = sw, skewness = sk,
       warnings = c(if (ext) sprintf("%d extreme outlier(s) among the %s", ext, what),
                    if (!is.null(sw) && sw < 0.05) sprintf("%s depart from normality (Shapiro-Wilk p = %.3g)", what, sw),
                    if (!normality && is.finite(sk) && abs(sk) > (sr$skew_note %||% 1)) sprintf("%s are skewed (%.2f): the signed-rank test assumes symmetry; the sign test does not", what, sk)))
}

# Descriptives per treatment cell (rstatix get_summary_stats "common" plus quartiles and skewness); binary: n, successes, proportion.
.ba_descriptives <- function(frame, p, type, sr) {
  if (!type %in% c("continuous", "ordinal", "binary")) return(NULL)
  fx <- intersect(unique(unlist(lapply(p$fixed_labels, function(l) all.vars(str2lang(l))))), names(frame))
  gv <- fx[vapply(fx, function(v) identical(p$role[[v]] %||% "treatment", "treatment") && !is.numeric(frame[[v]]), logical(1))]
  y <- if (is.factor(frame$y_)) as.integer(frame$y_) else as.numeric(frame$y_)
  cell <- if (length(gv)) interaction(frame[gv], drop = TRUE, sep = " / ") else factor(rep("all", nrow(frame)))
  out <- lapply(levels(cell), function(cl) {
    v <- y[cell == cl & !is.na(y)]; n <- length(v)
    key <- if (length(gv)) as.list(stats::setNames(strsplit(cl, " / ", fixed = TRUE)[[1]], gv)) else list(group = "all")
    if (type == "binary") return(as.data.frame(c(key, list(n = n, successes = sum(v), proportion = mean(v))), stringsAsFactors = FALSE))
    q <- stats::quantile(v, c(.25, .5, .75)); se <- stats::sd(v) / sqrt(n); sk <- .ba_skew(v)
    as.data.frame(c(key, list(n = n, mean = mean(v), sd = stats::sd(v), se = se, ci = stats::qt(0.975, n - 1) * se, median = unname(q[2]),
                              q1 = unname(q[1]), q3 = unname(q[3]), iqr = unname(q[3] - q[1]), min = min(v), max = max(v), skewness = sk,
                              report = if (is.finite(sk) && abs(sk) > (sr$skew_report_median %||% 1)) "median (IQR)" else "mean (SD)")), stringsAsFactors = FALSE)
  })
  do.call(rbind, out)
}

# Biodes analysis.posthoc.factor (a spec factor name or model term) -> model term; NULL when not given;
# a character note when it names nothing in the model.
.ba_posthoc_factor <- function(req, p, treat) {
  f <- req$factor
  if (is.null(f)) return(NULL)
  hit <- treat[vapply(treat, function(v) identical(v, f) || identical(as.character(p$col[[v]]), f), logical(1))]
  if (length(hit)) hit[1] else sprintf("analysis.posthoc.factor '%s' is not a treatment factor of the model: default comparisons", f)
}

# Post-hoc after rank and proportion tests: analysis.posthoc type none drops it; vs_control keeps the
# comparisons with the control and adjusts again over those (Holm; dunnettx is an emmeans method).
.ba_posthoc_select <- function(ph, req, lev, adj, sr) {
  type <- req$type %||% sr$posthoc_type %||% "pairwise"
  if (identical(type, "none")) return(list(table = NULL, note = "analysis.posthoc type none: no post-hoc"))
  if (!identical(type, "vs_control")) return(list(table = ph))
  ctrl <- req$control %||% lev[1]; note <- character()
  if (!ctrl %in% lev) { note <- sprintf("control '%s' is not a level: '%s' used", ctrl, lev[1]); ctrl <- lev[1] }
  parts <- strsplit(ph$contrast, " - ", fixed = TRUE)
  keep <- vapply(parts, function(x) ctrl %in% x, logical(1))
  ph <- ph[keep, , drop = FALSE]
  a <- if (adj %in% stats::p.adjust.methods) adj else "holm"
  if (nrow(ph) && !all(is.na(ph$p))) { ph$p_adj <- stats::p.adjust(ph$p, a); ph$adjust <- a }
  list(table = ph, note = c(note, sprintf("comparisons with the control '%s' only", ctrl)))
}

# Post-hoc contrasts (emmeans) against Simitall's truth: the same emmeans contrasts on a saturated lm of the
# true cell means of the treatment factors (equal weights over the other factors, as for the estimates).
.ba_posthoc_recovery <- function(ph, phs, truth, frame, alpha) {
  facs <- unique(unlist(lapply(phs$specs, function(s) c(s$specs, s$by))))
  if (!length(facs) || !all(facs %in% names(frame))) return(NULL)
  grid <- expand.grid(lapply(frame[facs], function(x) levels(droplevels(factor(x)))), stringsAsFactors = TRUE)
  grid$true_ <- .ba_true_cell_means(truth, grid)
  tl <- stats::lm(stats::reformulate(paste(facs, collapse = " * "), response = "true_"), data = grid)
  tru <- do.call(rbind, lapply(phs$specs, function(s) {
    sp <- if (is.null(s$by)) stats::as.formula(paste("~", s$specs)) else stats::as.formula(sprintf("~ %s | %s", s$specs, s$by))
    em <- suppressMessages(suppressWarnings(emmeans::emmeans(tl, sp)))
    ct <- if (identical(phs$type, "vs_control")) emmeans::contrast(em, "trt.vs.ctrl", ref = match(phs$control[[s$specs]] %||% levels(grid[[s$specs]])[1], levels(grid[[s$specs]])))
          else graphics::pairs(em)
    x <- suppressWarnings(as.data.frame(summary(ct, adjust = "none")))
    data.frame(family = if (is.null(s$by)) s$specs else sprintf("%s | %s", s$specs, s$by), by = if (is.null(s$by)) NA_character_ else as.character(x[[s$by]]),
               contrast = as.character(x$contrast), true = x$estimate, stringsAsFactors = FALSE)
  }))
  key <- function(d) paste(d$family, ifelse(is.na(d$by), "", d$by), d$contrast, sep = "\r")
  ph$true <- tru$true[match(key(ph), key(tru))]
  ok <- !is.na(ph$true)
  if (!any(ok)) return(NULL)
  q <- stats::qt(1 - alpha / 2, ph$df[ok])
  covered <- abs(ph$estimate[ok] - ph$true[ok]) <= q * ph$se[ok]
  nz <- abs(ph$true[ok]) > 1e-9; sig <- ph$p_adj[ok] < alpha
  list(table = ph, summary = list(n = sum(ok), covered = sum(covered), coverage = round(mean(covered), 3), nominal = 1 - alpha,
                                  true_nonzero = sum(nz), detected = sum(sig & nz), true_zero = sum(!nz), false_calls = sum(sig & !nz)))
}
