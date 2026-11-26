# Robust tests on trimmed means (WRS2; Wilcox 2017 Introduction to Robust Estimation and Hypothesis Testing):
# rows of stat_tests.tsv with family "robust". Not chosen by default; run when asked (analysis.stat_tests:
# ["robust"] or a test id) and, by stat_rules.json robust_when_outliers, next to the primary test when its
# checks find extreme outliers. The primary test does not change.

.ba_robust_ids <- c("yuen", "yuen_paired", "t1way", "t2way", "rmanova_trimmed", "bwtrim")

.ba_st_robust <- function(row, frame, lay, p, design, sr) {
  if (!requireNamespace("WRS2", quietly = TRUE)) stop("robust tests need WRS2", call. = FALSE)
  id <- row$test; tr <- sr$trim %||% 0.2
  d <- frame[!is.na(frame$y_), , drop = FALSE]; d$yn_ <- as.numeric(d$y_)
  tv <- lay$treat
  fr <- function(term, contrast = NA_character_, stat_name, stat, df1 = NA_real_, df2 = NA_real_, pv, est = NA_real_, lo = NA_real_, hi = NA_real_, eff_name = NA_character_, eff = NA_real_)
    data.frame(term = term, contrast = contrast, statistic_name = stat_name, statistic = unname(stat), df1 = unname(df1), df2 = unname(df2), p = unname(pv),
               estimate = unname(est), ci_low = unname(lo), ci_high = unname(hi), effect_name = eff_name, effect = unname(eff), n = nrow(d), stringsAsFactors = FALSE)
  ph <- NULL; extra <- list(trim = tr); adj <- "hochberg"
  ph_tab <- function(comp, lev, fam) {
    data.frame(family = fam, by = NA_character_, contrast = sprintf("%s - %s", lev[comp[, 1]], lev[comp[, 2]]), estimate = comp[, "psihat"], se = NA_real_, df = NA_real_,
               statistic = NA_real_, p = comp[, "p.value"], p_adj = comp[, "p.value"], adjust = "hochberg (lincon)", stringsAsFactors = FALSE)
  }
  if (id == "yuen") {
    d$g_ <- droplevels(factor(d[[tv[1]]])); lev <- levels(d$g_)
    y <- WRS2::yuen(yn_ ~ g_, data = d, tr = tr)
    ef <- tryCatch({ set.seed(sr$seed %||% 123); WRS2::yuen.effect.ci(yn_ ~ g_, data = d, tr = tr) }, error = function(e) NULL)
    flat <- fr(tv[1], sprintf("%s - %s", lev[1], lev[2]), "Yuen's t (trimmed means)", y$test, y$df, NA, y$p.value, y$diff, y$conf.int[1], y$conf.int[2],
               "explanatory measure of effect size (Wilcox & Tian)", y$effsize)
    if (!is.null(ef)) extra$effsize_ci <- as.numeric(ef$CI)
  } else if (id == "yuen_paired") {
    d$g_ <- droplevels(factor(d[[tv[1]]])); lev <- levels(d$g_)
    w <- tapply(d$yn_, list(as.character(d[[lay$subject]]), d$g_), mean); w <- w[stats::complete.cases(w), , drop = FALSE]
    y <- WRS2::yuend(w[, 1], w[, 2], tr = tr)
    flat <- fr(tv[1], sprintf("%s - %s", lev[1], lev[2]), "paired Yuen's t (trimmed means of differences)", y$test, y$df, NA, y$p.value, y$diff, y$conf.int[1], y$conf.int[2],
               "explanatory measure of effect size", y$effsize)
  } else if (id == "t1way") {
    d$g_ <- droplevels(factor(d[[tv[1]]])); lev <- levels(d$g_)
    set.seed(sr$seed %||% 123)
    t1 <- WRS2::t1way(yn_ ~ g_, data = d, tr = tr)
    flat <- fr(tv[1], NA, "robust one-way F (Welch on trimmed means)", t1$test, t1$df1, t1$df2, t1$p.value, eff_name = "explanatory measure of effect size", eff = t1$effsize)
    extra$effsize_ci <- as.numeric(t1$effsize_ci)
    lc <- WRS2::lincon(yn_ ~ g_, data = d, tr = tr); ph <- ph_tab(lc$comp, lev, tv[1])
  } else if (id == "t2way") {
    if (length(tv) != 2L) stop("t2way needs exactly two treatment factors", call. = FALSE)
    d$a_ <- droplevels(factor(d[[tv[1]]])); d$b_ <- droplevels(factor(d[[tv[2]]]))
    t2 <- WRS2::t2way(yn_ ~ a_ * b_, data = d, tr = tr)
    flat <- rbind(fr(tv[1], NA, "robust two-way Q (trimmed means)", t2$Qa, pv = t2$A.p.value), fr(tv[2], NA, "robust two-way Q (trimmed means)", t2$Qb, pv = t2$B.p.value),
                  fr(paste(tv, collapse = ":"), NA, "robust two-way Q (trimmed means)", t2$Qab, pv = t2$AB.p.value))
  } else if (id == "rmanova_trimmed") {
    d$g_ <- droplevels(factor(d[[tv[1]]])); lev <- levels(d$g_)
    w <- tapply(d$yn_, list(as.character(d[[lay$subject]]), d$g_), mean); w <- w[stats::complete.cases(w), , drop = FALSE]
    lw <- data.frame(id_ = factor(rep(rownames(w), ncol(w))), g_ = factor(rep(colnames(w), each = nrow(w)), colnames(w)), y = c(w))
    rm <- WRS2::rmanova(lw$y, lw$g_, lw$id_, tr = tr)
    flat <- fr(tv[1], NA, "robust repeated-measures F (trimmed means)", rm$test, rm$df1, rm$df2, rm$p.value)
    mc <- tryCatch(WRS2::rmmcp(lw$y, lw$g_, lw$id_, tr = tr), error = function(e) NULL)
    if (!is.null(mc)) { cm <- mc$comp
      ph <- data.frame(family = tv[1], by = NA_character_, contrast = sprintf("%s - %s", lev[cm[, 1]], lev[cm[, 2]]), estimate = cm[, "psihat"], se = NA_real_, df = NA_real_,
                       statistic = cm[, "test"], p = cm[, "p.value"], p_adj = cm[, "p.crit"] * 0 + cm[, "p.value"], adjust = "Hochberg critical values (rmmcp)", stringsAsFactors = FALSE) }
  } else if (id == "bwtrim") {
    bv <- lay$between[1]; wv <- lay$within[1]
    w <- stats::aggregate(d$yn_, list(id_ = d[[lay$subject]], b_ = d[[bv]], w_ = d[[wv]]), mean)
    bw <- WRS2::bwtrim(x ~ b_ * w_, id = id_, data = w, tr = tr)
    flat <- rbind(fr(bv, NA, "robust mixed Q (trimmed means)", bw$Qa, pv = bw$A.p.value), fr(wv, NA, "robust mixed Q (trimmed means)", bw$Qb, pv = bw$B.p.value),
                  fr(paste(bv, wv, sep = ":"), NA, "robust mixed Q (trimmed means)", bw$Qab, pv = bw$AB.p.value))
  } else stop(sprintf("no robust adapter for '%s'", id), call. = FALSE)
  list(summary = c(list(status = "run", test = id, name = row$model, layout = lay$layout, package = "WRS2", version = as.character(utils::packageVersion("WRS2")),
                        call = row$calls, factor = paste(tv, collapse = " x "), table = lapply(seq_len(nrow(flat)), function(i) as.list(flat[i, ])),
                        posthoc = if (!is.null(ph)) list(status = "run", adjust = ph$adjust[1], n_comparisons = nrow(ph), n_significant = sum(ph$p_adj < 0.05, na.rm = TRUE), families = list(tv[1])),
                        reference = row$reference, tutorial = row$doc, alternative = row$alternative), extra),
       posthoc = ph, flat = flat)
}
