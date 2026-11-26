# Survival tests next to the Cox engine (rules inst/rules/survival_rules.json, rows in stat_tests.tsv):
# Kaplan-Meier curves, medians and quantiles with CIs, landmark survival, log-rank (stratified by the
# design's blocks), weighted log-rank on request, pairwise log-rank, hazard ratios and the proportional-
# hazards check from the Cox fit, restricted mean survival time, AFT models and life tables on request.
# All numbers come from survival (survfit, survdiff, cox.zph, survreg), survMisc, KMsurv.

.ba_survival_tests <- function(frame, p, design, goal, alpha = 0.05, args = list(), fit = NULL) {
  sr <- .ba_rules("survival_rules")
  sv <- utils::modifyList(design$spec$analysis$survival %||% list(), args$survival %||% list())
  role <- function(v) p$role[[v]] %||% "treatment"
  fx <- intersect(unique(unlist(lapply(p$fixed_labels, function(l) all.vars(str2lang(l))))), names(frame))
  treat <- fx[vapply(fx, function(v) role(v) == "treatment" && !is.numeric(frame[[v]]), logical(1))]
  rblock <- unlist(lapply(p$random, function(r) if (isTRUE(r$intercept_only) && length(r$vars) == 1L && role(r$vars) == "block") r$vars))
  blocks <- intersect(c(fx[vapply(fx, function(v) role(v) == "block", logical(1))], rblock), names(frame))
  extra_strata <- .ba_survival_strata(sv$strata, p, frame)
  strata_v <- unique(c(if (isTRUE(sr$stratify_by_blocks)) blocks, extra_strata))
  d <- frame[!is.na(frame$y_) & !is.na(frame$event_), , drop = FALSE]
  comp <- NULL; notes0 <- character()
  if (any(!d$event_ %in% c(0, 1))) {
    # competing events: the CIF / Gray / Fine-Gray analysis next to the cause-specific one (competing events
    # censored: Kaplan-Meier is then survival net of the competing risk, the log-rank and Cox are cause-specific,
    # as the engine's Cox model)
    comp <- .ba_competing_tests(d, treat, strata_v, p, design, alpha, sv, sr)
    d$event_ <- as.integer(d$event_ == 1)
    notes0 <- "competing events censored here (cause-specific: Kaplan-Meier is survival net of the competing risk); cumulative incidence in the competing-risks block"
  }
  d$grp_ <- if (length(treat)) interaction(d[treat], drop = TRUE, sep = " / ") else factor(rep("all", nrow(d)))
  lev <- levels(d$grp_); k <- length(lev)
  senv <- list2env(list(strata = survival::strata), parent = baseenv())   # strata() must be visible as a special
  S <- function(rhs) stats::as.formula(paste("survival::Surv(y_, event_) ~", rhs), env = senv)
  notes <- character(); rows <- list(); ph <- NULL
  row <- function(term, name, stat = NA_real_, df1 = NA_real_, pv = NA_real_, est = NA_real_, lo = NA_real_, hi = NA_real_, eff_name = NA_character_, eff = NA_real_, contrast = NA_character_)
    data.frame(term = term, contrast = contrast, statistic_name = name, statistic = unname(stat), df1 = unname(df1), df2 = NA_real_, p = unname(pv),
               estimate = unname(est), ci_low = unname(lo), ci_high = unname(hi), effect_name = eff_name, effect = unname(eff), n = nrow(d), stringsAsFactors = FALSE)

  # Kaplan-Meier: curve, medians and quantiles, landmarks, RMST
  km <- survival::survfit(S(if (k > 1) "grp_" else "1"), data = d, conf.type = sr$conf_type %||% "log")
  sm <- summary(km, censored = TRUE)
  curve <- data.frame(group = if (k > 1) sub("^grp_=", "", as.character(sm$strata)) else "all", time = sm$time, n_risk = sm$n.risk, n_event = sm$n.event,
                      n_censor = sm$n.censor, surv = sm$surv, lower = sm$lower, upper = sm$upper, stringsAsFactors = FALSE)
  maxt <- tapply(d$y_, d$grp_, max)
  tau <- if (!is.null(sv$rmst_tau)) as.numeric(sv$rmst_tau) else min(maxt)
  tau_src <- if (!is.null(sv$rmst_tau)) "analysis.survival.rmst_tau" else "smallest of the groups' largest follow-up times"
  tb <- summary(km, rmean = tau)$table
  if (is.null(dim(tb))) tb <- t(as.matrix(tb))
  qs <- stats::quantile(km, probs = unlist(sr$quantiles %||% c(.25, .5, .75)))
  qmat <- function(x) { x <- as.matrix(x); if (nrow(x) != k) x <- t(x); x }
  qe <- qmat(qs$quantile); ql <- qmat(qs$lower); qu <- qmat(qs$upper)
  summ <- data.frame(group = lev, n = tb[, "records"], events = tb[, "events"], median = tb[, "median"], median_lcl = tb[, "0.95LCL"], median_ucl = tb[, "0.95UCL"],
                     q25 = qe[, 1], q75 = qe[, ncol(qe)], rmst = tb[, "rmean"], rmst_se = tb[, "se(rmean)"], tau = tau, max_followup = as.numeric(maxt[lev]),
                     row.names = NULL, stringsAsFactors = FALSE)
  if (any(is.na(summ$median))) notes <- c(notes, sprintf("median not reached in %s: report a landmark survival or the RMST", paste(summ$group[is.na(summ$median)], collapse = ", ")))
  land <- NULL
  if (length(sv$landmarks)) {
    lt <- summary(km, times = as.numeric(unlist(sv$landmarks)), extend = TRUE)
    land <- data.frame(group = if (k > 1) sub("^grp_=", "", as.character(lt$strata)) else "all", time = lt$time, n_risk = lt$n.risk, surv = lt$surv,
                       lower = lt$lower, upper = lt$upper, stringsAsFactors = FALSE)
  }

  clustered <- "unit_id_" %in% names(d) && anyDuplicated(d$unit_id_) > 0
  if (k > 1) {
    # log-rank (stratified by blocks when the design has them), weighted on request
    lr <- survival::survdiff(S("grp_"), data = d)
    rows[[length(rows) + 1L]] <- row("treatment", "log-rank chi-squared", lr$chisq, k - 1, stats::pchisq(lr$chisq, k - 1, lower.tail = FALSE))
    oe <- data.frame(group = lev, observed = as.numeric(lr$obs), expected = as.numeric(lr$exp), stringsAsFactors = FALSE)
    primary <- "log-rank"
    if (clustered) {
      # subsamples clustered in units: robust (Lin-Wei) score test, the log-rank with a cluster-robust variance
      rc <- survival::coxph(stats::as.formula("survival::Surv(y_, event_) ~ grp_ + cluster(unit_id_)", env = list2env(list(cluster = survival::cluster), parent = baseenv())), data = d)
      rs <- summary(rc)$robscore
      rows[[length(rows) + 1L]] <- row("treatment", "robust score chi-squared (Lin-Wei, clustered by unit)", rs[["test"]], rs[["df"]], rs[["pvalue"]])
      primary <- "robust score (clustered)"
      notes <- c(notes, sprintf("%d subsamples in %d units: the robust score test is primary; the plain log-rank, KM and RMST intervals treat subsamples as independent",
                                nrow(d), length(unique(d$unit_id_))))
    }
    if (length(strata_v)) {
      ls <- survival::survdiff(S(sprintf("grp_ + strata(%s)", paste(strata_v, collapse = ", "))), data = d)
      rows[[length(rows) + 1L]] <- row("treatment", sprintf("stratified log-rank chi-squared (strata: %s)", paste(strata_v, collapse = ", ")), ls$chisq, k - 1,
                                       stats::pchisq(ls$chisq, k - 1, lower.tail = FALSE))
      primary <- "stratified log-rank"
    }
    w <- sv$logrank_weight %||% sr$logrank_weight %||% "logrank"
    if (!identical(w, "logrank")) {
      wr <- .ba_weighted_logrank(d, w)
      if (is.character(wr)) notes <- c(notes, wr) else rows[[length(rows) + 1L]] <- row("treatment", sprintf("weighted log-rank chi-squared (%s, prespecified)", w), wr$chisq, wr$df, wr$p)
    }
    # pairwise log-rank (3+ groups): survdiff on each pair, adjusted (analysis.posthoc type / control)
    if (k > 2) ph <- .ba_pairwise_logrank(d, lev, design$spec$analysis$posthoc, sr)
  } else { oe <- NULL; primary <- NULL; notes <- c(notes, "one group: no comparison") }

  # Cox fit from the engine: hazard ratios, global tests, concordance, proportional-hazards check
  cox <- NULL; zph <- NULL
  if (!is.null(fit) && inherits(fit, "coxme")) {
    # coxme (random blocks / units as shared frailties): fixed-effect hazard ratios with Wald intervals; cox.zph
    # does not take coxme fits, so the PH check runs on the fixed part (coxph without the random terms)
    b <- coxme::fixef(fit); se <- sqrt(diag(as.matrix(stats::vcov(fit)))); z <- b / se; q <- stats::qnorm(0.975)
    for (t in names(b)) rows[[length(rows) + 1L]] <- row(t, "hazard ratio (Cox, Wald z)", z[[t]], NA, 2 * stats::pnorm(-abs(z[[t]])),
                                                         exp(b[[t]]), exp(b[[t]] - q * se[[t]]), exp(b[[t]] + q * se[[t]]), "hazard ratio", exp(b[[t]]), t)
    cox <- list(n = fit$n[2], events = fit$n[1], random = lapply(coxme::VarCorr(fit), function(v) unname(v)))
    fx <- tryCatch(survival::coxph(stats::as.formula(paste("survival::Surv(y_, event_) ~", if (length(p$fixed_labels)) paste(p$fixed_labels, collapse = " + ") else "1")), data = d), error = function(e) NULL)
    z2 <- if (!is.null(fx) && length(stats::coef(fx))) tryCatch(survival::cox.zph(fx), error = function(e) NULL)
    if (!is.null(z2)) { zt <- z2$table; gp <- unname(zt["GLOBAL", "p"])
      zph <- lapply(rownames(zt), function(t) list(term = t, chisq = unname(zt[t, "chisq"]), df = unname(zt[t, "df"]), p = unname(zt[t, "p"]), on = "fixed part (coxph without random terms)"))
      if (gp < (sr$ph_alpha %||% 0.05)) notes <- c(notes, sprintf("proportional hazards rejected on the fixed part (cox.zph global p = %.3g)", gp)) }
  }
  if (!is.null(fit) && inherits(fit, "coxph")) {
    s <- summary(fit)
    ci <- s$conf.int; co <- s$coefficients
    keep <- !grepl("frailty\\(", rownames(ci)) & rownames(ci) %in% rownames(co)
    # coxph: z and Pr(>|z|); a penalized (frailty) fit: coef, se(coef), ..., p
    zc <- if ("z" %in% colnames(co)) co[, "z"] else co[, "coef"] / co[, "se(coef)"]
    pc <- if ("Pr(>|z|)" %in% colnames(co)) co[, "Pr(>|z|)"] else co[, "p"]
    names(zc) <- names(pc) <- rownames(co)
    for (t in rownames(ci)[keep]) rows[[length(rows) + 1L]] <- row(t, "hazard ratio (Cox, Wald z)", zc[[t]], NA, pc[[t]],
                                                                 ci[t, "exp(coef)"], ci[t, "lower .95"], ci[t, "upper .95"], "hazard ratio", ci[t, "exp(coef)"], t)
    cox <- list(n = s$n, events = s$nevent, concordance = unname(s$concordance[1]), concordance_se = unname(s$concordance[2]),
                lr_test = unname(s$logtest[c("test", "df", "pvalue")]), wald_test = unname(s$waldtest[c("test", "df", "pvalue")]),
                score_test = unname(s$sctest[c("test", "df", "pvalue")]))
    epv <- s$nevent / max(1, sum(keep))
    if (epv < (sr$min_events_per_term %||% 10)) notes <- c(notes, sprintf("%.1f events per Cox coefficient (< %g): estimates may be overfitted", epv, sr$min_events_per_term %||% 10))
    z <- tryCatch(survival::cox.zph(fit), error = function(e) NULL)
    if (!is.null(z)) {
      zt <- z$table
      zph <- lapply(rownames(zt), function(t) list(term = t, chisq = unname(zt[t, "chisq"]), df = unname(zt[t, "df"]), p = unname(zt[t, "p"])))
      gp <- zt["GLOBAL", "p"]
      if (gp < (sr$ph_alpha %||% 0.05)) notes <- c(notes, sprintf("proportional hazards rejected (cox.zph global p = %.3g): read the hazard ratio as an average; the RMST%s does not assume PH; stratify a nuisance factor (analysis.survival.strata) or ask for an AFT model",
                                                         gp, if (k > 1) " difference" else ""))
    }
  }

  # RMST difference and ratio (two groups; independent groups, survfit's restricted means and SEs)
  rmst <- NULL
  if (k == 2) {
    r1 <- summ$rmst[2]; r0 <- summ$rmst[1]; s1 <- summ$rmst_se[2]; s0 <- summ$rmst_se[1]
    dlt <- r1 - r0; se <- sqrt(s1^2 + s0^2); zq <- stats::qnorm(1 - alpha / 2)
    lr_ <- log(r1 / r0); se_l <- sqrt((s1 / r1)^2 + (s0 / r0)^2)
    rmst <- list(tau = tau, tau_source = tau_src, contrast = sprintf("%s - %s", lev[2], lev[1]), difference = dlt, ci = c(dlt - zq * se, dlt + zq * se),
                 p = 2 * stats::pnorm(-abs(dlt / se)), ratio = exp(lr_), ratio_ci = exp(c(lr_ - zq * se_l, lr_ + zq * se_l)), ratio_p = 2 * stats::pnorm(-abs(lr_ / se_l)))
    rows[[length(rows) + 1L]] <- row("treatment", sprintf("RMST difference to tau = %g", tau), dlt / se, NA, rmst$p, dlt, rmst$ci[1], rmst$ci[2], "RMST ratio", rmst$ratio, rmst$contrast)
  }

  # AFT models on request
  aft <- NULL
  if (!is.null(sv$aft) && !isFALSE(sv$aft) && k > 1) aft <- .ba_aft(d, sv$aft, sr, alpha)
  if (!is.null(aft$rows)) rows <- c(rows, aft$rows)

  # life table on request (actuarial, KMsurv::lifetab)
  life <- NULL
  if (!is.null(sv$life_table) && !isFALSE(sv$life_table)) life <- .ba_life_table(d, sv$life_table)

  # events against the design's event-based power (Biodes power.events_required)
  pw <- design$spec$power
  ev_req <- suppressWarnings(as.numeric(pw$events_required %||% NA))
  if (is.finite(ev_req) && sum(d$event_ == 1) < ev_req) notes <- c(notes, sprintf("%d events observed, %g required by the design (power.events_required): power is below the plan", sum(d$event_ == 1), ev_req))
  flat <- if (length(rows)) do.call(rbind, rows)
  notes <- c(notes0, notes)
  list(summary = list(status = "run", test = "survival", package = "survival", version = as.character(utils::packageVersion("survival")),
                      groups = as.list(lev), factor = if (length(treat)) paste(treat, collapse = " x ") else NULL, strata = as.list(strata_v),
                      events = sum(d$event_ == 1), events_required = if (is.finite(ev_req)) ev_req, clustered = clustered,
                      primary = primary, table = lapply(seq_len(nrow(flat %||% data.frame())), function(i) as.list(flat[i, ])),
                      km = lapply(seq_len(nrow(summ)), function(i) as.list(summ[i, ])), observed_expected = if (!is.null(oe)) lapply(seq_len(nrow(oe)), function(i) as.list(oe[i, ])),
                      cox = cox, ph_check = zph, rmst = rmst, aft = aft$summary, life_table = !is.null(life),
                      posthoc = if (!is.null(ph)) list(status = "run", adjust = attr(ph, "adjust"), n_comparisons = nrow(ph), n_significant = sum(ph$p_adj < 0.05, na.rm = TRUE),
                                                      families = list(paste(treat, collapse = " x ")), type = attr(ph, "type")),
                      reference = "Kaplan & Meier 1958 JASA 53:457; Mantel 1966 Cancer Chemother Rep 50:163; Cox 1972 JRSS B 34:187; Grambsch & Therneau 1994 Biometrika 81:515; Royston & Parmar 2013 BMC Med Res Methodol 13:152; Therneau survival (CRAN)",
                      tutorial = "docs/engines/survival.md", note = if (length(notes)) paste(notes, collapse = "; "),
                      also = if (!is.null(comp)) list(comp$summary)),
       posthoc = ph, flat = flat, stat_table = flat, km_curve = curve, survival_summary = summ, landmarks = land, life_table = life,
       cif_curve = comp$cif_curve, competing_table = comp$flat)
}

# Weighted log-rank by name (survMisc::comp rows; peto also = survdiff rho = 1).
.ba_weighted_logrank <- function(d, w, strata = NULL) {
  # k-sample weighted log-rank in base R (no survMisc): U_g = sum w (d_g - n_g d / n), V from the hypergeometric
  # variance; weights at each event time: logrank 1, gehan n, tarone_ware sqrt(n), fh_p_q S(t-)^p (1 - S(t-))^q
  # with S the pooled Kaplan-Meier just before t (Fleming & Harrington 1991; Klein & Moeschberger 7.3).
  # strata: U and V summed over strata, each with its own pooled KM (stratified test). Equal to survMisc::comp
  # and survival::survdiff (tests). peto keeps survdiff(rho = 1).
  if (identical(w, "peto") && is.null(strata)) { x <- survival::survdiff(survival::Surv(y_, event_) ~ grp_, data = d, rho = 1); df <- length(x$n) - 1
    return(list(chisq = x$chisq, df = df, p = stats::pchisq(x$chisq, df, lower.tail = FALSE))) }
  fh <- regmatches(w, regexec("^fh_([0-9.]+)_([0-9.]+)$", w))[[1]]
  if (identical(w, "peto")) fh <- c("fh_1_0", "1", "0")   # stratified: Peto-Peto ~ FH(1, 0)
  if (!length(fh) && !w %in% c("logrank", "gehan", "tarone_ware"))
    return(sprintf("weighted log-rank '%s' is not a known weight (survival_rules.json logrank_weight)", w))
  g <- droplevels(as.factor(d$grp_)); k <- nlevels(g); if (k < 2) return("weighted log-rank needs 2 or more groups")
  one <- function(dd, gg) {
    tt <- sort(unique(dd$y_[dd$event_ == 1])); if (!length(tt)) return(NULL)
    nr <- vapply(tt, function(t) dd$y_ >= t, logical(nrow(dd))); ev <- vapply(tt, function(t) dd$y_ == t & dd$event_ == 1, logical(nrow(dd)))
    if (is.null(dim(nr))) { nr <- matrix(nr, nrow = nrow(dd)); ev <- matrix(ev, nrow = nrow(dd)) }
    G <- stats::model.matrix(~ gg - 1)
    ng <- t(G) %*% nr; dg <- t(G) %*% ev; n <- colSums(nr); dd_ <- colSums(ev)
    s_minus <- c(1, cumprod(1 - dd_ / n))[seq_along(tt)]
    wt <- if (length(fh)) s_minus^as.numeric(fh[2]) * (1 - s_minus)^as.numeric(fh[3]) else switch(w, logrank = rep(1, length(tt)), gehan = n, tarone_ware = sqrt(n))
    idx <- seq_len(k - 1)
    U <- as.vector((dg - ng * rep(dd_ / n, each = k)) %*% wt)[idx]
    vf <- wt^2 * dd_ * ifelse(n > 1, (n - dd_) / (n - 1), 0) / n
    V <- matrix(0, k - 1, k - 1)
    for (a in idx) for (b in idx) V[a, b] <- sum(vf * ng[a, ] * ((a == b) - ng[b, ] / n))
    list(U = U, V = V)
  }
  parts <- if (is.null(strata)) list(one(d, g)) else lapply(split(seq_len(nrow(d)), droplevels(interaction(d[strata], drop = TRUE))), function(i) one(d[i, , drop = FALSE], g[i]))
  parts <- Filter(Negate(is.null), parts)
  U <- Reduce(`+`, lapply(parts, `[[`, "U")); V <- Reduce(`+`, lapply(parts, `[[`, "V"))
  chisq <- as.numeric(t(U) %*% MASS::ginv(V) %*% U); df <- k - 1
  list(chisq = chisq, df = df, p = stats::pchisq(chisq, df, lower.tail = FALSE))
}

.ba_pairwise_logrank <- function(d, lev, req, sr) {
  type <- req$type %||% "pairwise"
  if (identical(type, "none")) return(NULL)
  pr <- utils::combn(lev, 2, simplify = FALSE)
  ctrl <- req$control %||% lev[1]
  if (identical(type, "vs_control")) pr <- Filter(function(x) ctrl %in% x, pr)
  pv <- vapply(pr, function(x) { s <- d[d$grp_ %in% x, , drop = FALSE]; s$grp_ <- droplevels(s$grp_)
    x2 <- survival::survdiff(survival::Surv(y_, event_) ~ grp_, data = s); stats::pchisq(x2$chisq, 1, lower.tail = FALSE) }, 0)
  adj <- req$adjust %||% sr$pairwise_adjust %||% "holm"; if (!adj %in% stats::p.adjust.methods) adj <- "holm"
  out <- data.frame(family = "log-rank", by = NA_character_, contrast = vapply(pr, function(x) sprintf("%s - %s", x[2], x[1]), ""), estimate = NA_real_, se = NA_real_,
                    df = 1, statistic = NA_real_, p = pv, p_adj = stats::p.adjust(pv, adj), adjust = adj, stringsAsFactors = FALSE)
  attr(out, "adjust") <- adj; attr(out, "type") <- type
  out
}

# AFT (survreg) over the requested or all distributions; time ratios of the chosen one.
.ba_aft <- function(d, req, sr, alpha) {
  dists <- if (isTRUE(req)) unlist(sr$aft_distributions) else unlist(req)
  fits <- lapply(stats::setNames(dists, dists), function(ds) tryCatch(survival::survreg(survival::Surv(y_, event_) ~ grp_, data = d, dist = ds), error = function(e) NULL))
  fits <- Filter(Negate(is.null), fits)
  if (!length(fits)) return(list(summary = list(status = "error", note = "no AFT model converged")))
  aic <- vapply(fits, stats::AIC, 0); best <- names(which.min(aic)); f <- fits[[best]]
  co <- summary(f)$table; tr <- rownames(co)[startsWith(rownames(co), "grp_")]
  zq <- stats::qnorm(1 - alpha / 2)
  shape <- if (best %in% c("weibull", "exponential")) 1 / f$scale else NA_real_
  rows <- lapply(tr, function(t) {
    b <- co[t, "Value"]; se <- co[t, "Std. Error"]
    data.frame(term = "treatment", contrast = sub("^grp_", "", t), statistic_name = sprintf("time ratio (AFT %s, Wald z)", best), statistic = co[t, "z"], df1 = NA_real_, df2 = NA_real_,
               p = co[t, "p"], estimate = exp(b), ci_low = exp(b - zq * se), ci_high = exp(b + zq * se), effect_name = if (is.finite(shape)) "hazard ratio (Weibull)" else "time ratio",
               effect = if (is.finite(shape)) exp(b)^(-shape) else exp(b), n = nrow(d), stringsAsFactors = FALSE)
  })
  close <- names(aic)[aic - min(aic) < 2 & names(aic) != best]
  list(rows = rows, summary = list(status = "run", chosen = best, aic = as.list(aic), shape = shape,
                                   note = if (length(close)) sprintf("AIC within 2 of %s: %s (indistinguishable)", best, paste(close, collapse = ", "))))
}

# Actuarial life table per group (KMsurv::lifetab) on intervals of the requested width (or TRUE: 10 intervals).
.ba_life_table <- function(d, width) {
  if (!requireNamespace("KMsurv", quietly = TRUE)) return(NULL)
  mx <- max(d$y_)
  w <- if (isTRUE(width)) mx / 10 else as.numeric(width)
  br <- seq(0, ceiling(mx / w) * w + w, by = w)
  do.call(rbind, lapply(levels(d$grp_), function(g) {
    s <- d[d$grp_ == g, , drop = FALSE]
    bin <- cut(s$y_, br, right = FALSE, labels = FALSE)
    ne <- tabulate(bin[s$event_ == 1], length(br) - 1); nl <- tabulate(bin[s$event_ == 0], length(br) - 1)
    lt <- KMsurv::lifetab(br, nrow(s), nl, ne)
    data.frame(group = g, start = br[-length(br)], end = br[-1], lt, row.names = NULL, check.names = FALSE, stringsAsFactors = FALSE)
  }))
}

# Competing risks (event codes: 0 censored, 1 the event of interest, 2+ competing; names from
# primary_outcome.event_codes): cumulative incidence (Aalen-Johansen) and Gray's test per cause
# (cmprsk::cuminc), Fine-Gray subdistribution hazard ratio for the event of interest (cmprsk::crr), and
# 1 - Kaplan-Meier next to the CIF (treating competing events as censoring overstates the risk).
# The engine's Cox model is cause-specific (competing events censored: method_outcome.R).
.ba_competing_tests <- function(d, treat, strata_v, p, design, alpha, sv, sr) {
  if (!requireNamespace("cmprsk", quietly = TRUE)) return(list(summary = list(status = "not run", test = "competing_risks", note = "cmprsk not installed")))
  codes <- design$spec$primary_outcome$event_codes
  nm <- function(c) { x <- codes[[as.character(c)]]; if (is.null(x)) as.character(c) else as.character(x) }
  d$grp_ <- if (length(treat)) interaction(d[treat], drop = TRUE, sep = " / ") else factor(rep("all", nrow(d)))
  lev <- levels(d$grp_); k <- length(lev)
  ci <- if (k > 1) cmprsk::cuminc(d$y_, d$event_, d$grp_, cencode = 0) else cmprsk::cuminc(d$y_, d$event_, cencode = 0)
  causes <- sort(setdiff(unique(d$event_), 0))
  rows <- list(); notes <- character()
  if (k > 1 && !is.null(ci$Tests)) for (cz in rownames(ci$Tests))
    rows[[length(rows) + 1L]] <- data.frame(term = "treatment", contrast = nm(cz), statistic_name = sprintf("Gray's test chi-squared (%s)", nm(cz)), statistic = ci$Tests[cz, "stat"],
                                            df1 = ci$Tests[cz, "df"], df2 = NA_real_, p = ci$Tests[cz, "pv"], estimate = NA_real_, ci_low = NA_real_, ci_high = NA_real_,
                                            effect_name = NA_character_, effect = NA_real_, n = nrow(d), stringsAsFactors = FALSE)
  curves <- ci[setdiff(names(ci), "Tests")]
  curve <- do.call(rbind, lapply(names(curves), function(n) { x <- curves[[n]]
    parts <- strsplit(n, " ", fixed = TRUE)[[1]]; cz <- parts[length(parts)]; g <- if (k > 1) paste(parts[-length(parts)], collapse = " ") else "all"
    data.frame(group = g, cause = nm(cz), time = x$time, cif = x$est, var = x$var, stringsAsFactors = FALSE) }))
  # CIF of the event of interest vs the naive 1 - KM at the landmarks (else the median follow-up), per group and overall
  times <- if (length(sv$landmarks)) as.numeric(unlist(sv$landmarks)) else stats::median(d$y_)
  one <- function(dd, g) {
    cif <- cmprsk::timepoints(cmprsk::cuminc(dd$y_, dd$event_, cencode = 0), times)$est
    nv <- summary(survival::survfit(survival::Surv(y_, event_ == 1) ~ 1, data = dd), times = times, extend = TRUE)
    data.frame(group = g, time = times, cif = as.numeric(cif[grep("^1 1$|^1$", rownames(cif))[1], ]), one_minus_km = 1 - nv$surv, stringsAsFactors = FALSE)
  }
  cmp <- do.call(rbind, c(if (k > 1) lapply(lev, function(g) one(d[d$grp_ == g, ], g)), list(one(d, "all"))))
  over <- cmp[cmp$one_minus_km > cmp$cif + 0.01, , drop = FALSE]
  if (nrow(over)) { o <- over[nrow(over), ]; notes <- c(notes, sprintf("1 - Kaplan-Meier overstates the risk of %s (%s at %g: %.3f vs CIF %.3f)", nm(1), o$group, o$time, o$one_minus_km, o$cif)) }
  fg <- NULL
  if (k > 1) {
    X <- stats::model.matrix(~ grp_, d)[, -1, drop = FALSE]
    f <- tryCatch(cmprsk::crr(d$y_, d$event_, cov1 = X, failcode = 1, cencode = 0), error = function(e) NULL)
    if (!is.null(f)) {
      s <- summary(f)$coef; zq <- stats::qnorm(1 - alpha / 2)
      for (i in seq_len(nrow(s))) rows[[length(rows) + 1L]] <- data.frame(term = "treatment", contrast = sub("^grp_", "", colnames(X)[i]),
        statistic_name = sprintf("subdistribution hazard ratio (Fine-Gray, %s)", nm(1)), statistic = s[i, "z"], df1 = NA_real_, df2 = NA_real_, p = s[i, "p-value"],
        estimate = s[i, "exp(coef)"], ci_low = exp(s[i, "coef"] - zq * s[i, "se(coef)"]), ci_high = exp(s[i, "coef"] + zq * s[i, "se(coef)"]),
        effect_name = "subdistribution hazard ratio", effect = s[i, "exp(coef)"], n = nrow(d), stringsAsFactors = FALSE)
      fg <- list(event = nm(1), converged = isTRUE(f$converged))
    }
  }
  flat <- if (length(rows)) do.call(rbind, rows)
  list(summary = list(status = "run", test = "competing_risks", package = "cmprsk", version = as.character(utils::packageVersion("cmprsk")),
                      groups = as.list(lev), causes = as.list(vapply(causes, nm, "")), event_of_interest = nm(1), fine_gray = fg,
                      table = lapply(seq_len(nrow(flat %||% data.frame())), function(i) as.list(flat[i, ])),
                      cif_vs_naive = lapply(seq_len(nrow(cmp)), function(i) as.list(cmp[i, ])),
                      reference = "Gray 1988 Ann Stat 16:1141; Fine & Gray 1999 JASA 94:496; Putter et al. 2007 Stat Med 26:2389; Gray cmprsk (CRAN)",
                      tutorial = "docs/engines/survival.md", note = if (length(notes)) paste(notes, collapse = "; ")),
       flat = flat, stat_table = flat, cif_curve = curve)
}

# AFT fit matching Simitall's survival truth (outcome_model.family), for goal recovery: coefficient table with the
# Cox model's term names (same design matrix), on the log time-ratio scale.
.ba_aft_for_truth <- function(frame, p, family) {
  dist <- c(lognormal_aft = "lognormal", weibull_aft = "weibull", weibull_ph = "weibull", weibull = "weibull", exponential = "exponential",
            loglogistic_aft = "loglogistic")[[family]] %||% "lognormal"
  f <- stats::as.formula(paste("survival::Surv(y_, event_) ~", paste(c(p$fixed_labels, if (!length(p$fixed_labels)) "1"), collapse = " + ")))
  fit <- tryCatch(survival::survreg(f, data = frame, dist = dist), error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  co <- summary(fit)$table; co <- co[!rownames(co) %in% c("(Intercept)", "Log(scale)"), , drop = FALSE]
  list(res = data.frame(term = rownames(co), estimate = co[, "Value"], se = co[, "Std. Error"], df = Inf, statistic = co[, "z"], p = co[, "p"],
                        padj = co[, "p"], stringsAsFactors = FALSE, row.names = NULL),
       scale = sprintf("time ratio, estimate and p from survreg(dist = '%s'), matching Simitall's %s truth (the Cox hazard ratio is in the report)", dist, family))
}

# Goal truth with non-proportional hazards (Simitall delayed_effect / early_effect / crossing): one Cox HR averages
# the phases. Detection: the prespecified weighted log-rank; estimate: piecewise Cox split at change_time, the HR of
# the affected phase against the true HR (log scale), the other phase reported.
.ba_goal_recovery_nonph <- function(out, gt, om, goal, design, frame, p, res, alpha) {
  # Non-proportional hazards (Simitall delayed_effect / early_effect / crossing): detection by the prespecified
  # weighted log-rank (joint over the tested factors' groups; per factor stratified on the others), estimate by
  # a piecewise Cox fit split at change_time with a before / after term per tested two-level factor (robust SEs
  # by the unit when subsamples share one), each against Simitall's hazard_ratios (or true_effect).
  w <- design$spec$analysis$survival$logrank_weight %||% om$logrank_weight %||% sub("^logrank_", "", goal$effect$test %||% "")
  fam <- om$family %||% "delayed_effect"; t0 <- suppressWarnings(as.numeric(om$change_time %||% goal$effect$change_time))
  treat <- names(p$col)[vapply(names(p$col), function(v) identical(p$role[[v]] %||% "treatment", "treatment") && v %in% names(frame) && !is.numeric(frame[[v]]), logical(1))]
  tested <- intersect(unique(unlist(lapply(unlist(goal$tested_terms %||% out$tested_terms %||% treat), function(t) all.vars(str2lang(t))))), treat)
  if (!length(tested)) tested <- treat[1]
  if (!length(tested)) return(out)
  d <- frame[!is.na(frame$y_) & !is.na(frame$event_), , drop = FALSE]
  d$event_ <- as.integer(d$event_ == 1)
  for (v in tested) d[[v]] <- droplevels(factor(d[[v]]))
  two <- tested[vapply(tested, function(v) nlevels(d[[v]]) == 2L, logical(1))]
  hr_true <- function(v) { h <- gt$hazard_ratios[[v]]; if (length(h)) as.numeric(h[[levels(d[[v]])[2]]] %||% NA) else as.numeric(unlist(gt$true_effect))[1] }
  # Cox HR from the engine (one average over the phases), per tested term
  cox_avg <- vapply(two, function(v) { r <- res[res$term == paste0(v, levels(d[[v]])[2]), , drop = FALSE]; if (nrow(r)) exp(r$estimate[1]) else NA_real_ }, 0)
  out$proportional_hazards <- FALSE; out$cox_hr_average <- as.list(cox_avg)
  notes <- sprintf("hazards not proportional (%s): one Cox HR (%s) averages the phases", fam,
                   if (all(is.na(cox_avg))) "NA" else paste(sprintf("%s %s", names(cox_avg), format(signif(cox_avg, 3))), collapse = ", "))
  if (nzchar(w)) {
    d$grp_ <- interaction(d[tested], drop = TRUE, sep = " / ")
    wr <- .ba_weighted_logrank(d, w)
    if (!is.character(wr)) {
      out$p <- wr$p; out$detected <- isTRUE(wr$p < alpha)
      out$test <- sprintf("weighted log-rank %s (prespecified)%s", w, if (length(tested) > 1) sprintf(", joint over %d groups", nlevels(d$grp_)) else "")
      if (length(tested) > 1) out$weighted_logrank_by_factor <- lapply(tested, function(v) {
        dv <- d; dv$grp_ <- dv[[v]]; r <- .ba_weighted_logrank(dv, w, strata = setdiff(tested, v))
        if (is.character(r)) list(term = v, note = r) else list(term = v, chisq = r$chisq, df = r$df, p = r$p, stratified_on = paste(setdiff(tested, v), collapse = ", ")) })
    } else notes <- c(notes, wr)
  }
  if (is.finite(t0) && length(two)) {
    sf <- stats::as.formula("Surv(y_, event_) ~ .", env = list2env(list(Surv = survival::Surv), parent = baseenv()))   # survSplit wants a literal Surv() on the left
    unit <- intersect(c("unit_id_", design$spec$analysis$survival$cluster), names(d))[1]
    keep <- unique(c("y_", "event_", two, if (!is.na(unit)) unit))
    sp <- survival::survSplit(sf, data = d[keep], cut = t0, episode = "phase_", start = "tstart_")
    cv <- character()
    for (v in two) { x <- as.numeric(sp[[v]] == levels(d[[v]])[2]); sp[[paste0(v, "__before")]] <- x * (sp$phase_ == 1); sp[[paste0(v, "__after")]] <- x * (sp$phase_ == 2); cv <- c(cv, paste0(v, c("__before", "__after"))) }
    # bare Surv() / cluster() in an environment that has them: coxph only recognizes the cluster() special unqualified
    fml <- stats::as.formula(paste("Surv(tstart_, y_, event_) ~", paste(cv, collapse = " + "), if (!is.na(unit)) sprintf("+ cluster(%s)", unit) else ""),
                             env = list2env(list(Surv = survival::Surv, cluster = survival::cluster), parent = baseenv()))
    pc <- tryCatch(survival::coxph(fml, data = sp), error = function(e) NULL)
    if (!is.null(pc)) {
      b <- stats::coef(pc); se <- sqrt(diag(stats::vcov(pc))); q <- stats::qnorm(1 - alpha / 2)
      aff <- if (fam %in% c("early_effect", "crossing")) "__before" else "__after"
      tab <- do.call(rbind, lapply(two, function(v) { k <- paste0(v, aff); tv <- hr_true(v)
        data.frame(term = v, level = levels(d[[v]])[2], true_hr = tv, estimate = unname(exp(b[k])), ci_low = unname(exp(b[k] - q * se[k])), ci_high = unname(exp(b[k] + q * se[k])),
                   covered = unname(abs(b[k] - log(tv)) <= q * se[k]), hr_other_phase = unname(exp(b[paste0(v, if (aff == "__after") "__before" else "__after")])),
                   cox_hr_average = unname(cox_avg[v]), stringsAsFactors = FALSE) }))
      out$estimate <- tab$estimate[1]; out$true_effect <- tab$true_hr[1]; out$covered <- all(tab$covered)
      out$bias <- log(tab$estimate[1]) - log(tab$true_hr[1]); out$sign_agreement <- all(sign(log(tab$estimate)) == sign(log(tab$true_hr)))
      out$piecewise <- list(change_time = t0, affected_phase = if (aff == "__after") "after" else "before", robust_by = if (!is.na(unit)) unit,
                            terms = lapply(seq_len(nrow(tab)), function(i) as.list(tab[i, ])),
                            hr_before = if (aff == "__after") tab$hr_other_phase[1] else tab$estimate[1], hr_after = if (aff == "__after") tab$estimate[1] else tab$hr_other_phase[1])
      out$scale <- sprintf("hazard ratio %s t = %g (piecewise Cox split at change_time; log scale)%s", if (aff == "__after") "after" else "before", t0,
                           if (nrow(tab) > 1) sprintf("; estimate / true_effect = %s, all %d terms in piecewise$terms", tab$term[1], nrow(tab)) else "")
    }
  } else if (!is.finite(t0)) notes <- c(notes, "no change_time in the truth: no piecewise estimate")
  out$note <- paste(c(out$note, notes), collapse = "; ")
  out
}

# analysis.survival.strata -> frame columns (a model term or a metadata column). Biodes 01d41ef: blocks named here
# are not also "+ block" in the model, so the Cox model stratifies on them (method_outcome.R surv_rhs) as the
# log-rank does.
.ba_survival_strata <- function(strata, p, frame) {
  intersect(unlist(lapply(unlist(strata), function(s) {
    v <- names(p$col)[vapply(names(p$col), function(v) identical(v, s) || identical(as.character(p$col[[v]]), s), logical(1))]
    if (length(v)) v else s })), names(frame))
}

# Per-group survival truth (Simitall outcome_model.groups / goal_truth.groups, keyed "a:b" by the treatment
# cells): median, survival net of competing events at the landmark days (Kaplan-Meier with competing events
# censored), cumulative incidence of the event (cmprsk::cuminc), RMST to the truth's horizon. Intervals: KM
# log-log / Greenwood, cuminc variance; they treat subsamples as independent (too narrow with shared frailties).
.ba_survival_group_recovery <- function(groups, frame, p, alpha = 0.05, tau = NULL) {
  if (!length(groups)) return(NULL)
  treat <- names(p$col)[vapply(names(p$col), function(v) identical(p$role[[v]] %||% "treatment", "treatment") && v %in% names(frame) && !is.numeric(frame[[v]]), logical(1))]
  d <- frame[!is.na(frame$y_) & !is.na(frame$event_), , drop = FALSE]
  ev_raw <- d$event_raw_ %||% d$event_
  d$e1_ <- as.integer(ev_raw == 1)
  d$cell_ <- if (length(treat)) do.call(paste, c(lapply(d[treat], as.character), sep = ":")) else "all"
  q <- stats::qnorm(1 - alpha / 2); rows <- list()
  add <- function(g, m, tv, est, lo, hi) rows[[length(rows) + 1L]] <<- data.frame(group = g, measure = m, true = as.numeric(tv), estimate = as.numeric(est), ci_low = as.numeric(lo),
                                                                                 ci_high = as.numeric(hi), stringsAsFactors = FALSE)
  for (g in intersect(names(groups), unique(d$cell_))) {
    gt <- groups[[g]]; dg <- d[d$cell_ == g, , drop = FALSE]
    # Simitall: median = the event alone (net of competing events: KM with them censored); survival_at and rmst
    # = free of any event (KM of the first event of any kind), as its truth computes them
    km <- survival::survfit(survival::Surv(y_, e1_) ~ 1, data = dg, conf.type = "log-log", conf.int = 1 - alpha)
    ka <- survival::survfit(survival::Surv(y_, as.integer(ev_raw[d$cell_ == g] > 0)) ~ 1, data = dg, conf.type = "log-log", conf.int = 1 - alpha)
    tb <- summary(km)$table
    if (!is.null(gt$median)) add(g, "median (event alone)", gt$median, tb[["median"]], tb[[grep("LCL", names(tb), value = TRUE)[1]]], tb[[grep("UCL", names(tb), value = TRUE)[1]]])
    for (t in names(gt$survival_at)) { s <- summary(ka, times = as.numeric(t), extend = TRUE)
      add(g, sprintf("event-free survival at %s", t), gt$survival_at[[t]], s$surv, s$lower, s$upper) }
    if (!is.null(gt$rmst)) { ta <- if (length(tau) && is.finite(tau)) tau else max(dg$y_)
      r <- summary(ka, rmean = ta)$table; add(g, sprintf("event-free RMST to %g", ta), gt$rmst, r[["rmean"]], r[["rmean"]] - q * r[["se(rmean)"]], r[["rmean"]] + q * r[["se(rmean)"]]) }
    if (length(gt$cif_event_at) && requireNamespace("cmprsk", quietly = TRUE) && any(ev_raw[d$cell_ == g] > 1)) {
      ci <- cmprsk::cuminc(dg$y_, ev_raw[d$cell_ == g], cencode = 0)
      tp <- cmprsk::timepoints(ci, as.numeric(names(gt$cif_event_at)))
      k <- grep(" 1$", rownames(tp$est))[1]
      for (j in seq_along(gt$cif_event_at)) { e <- tp$est[k, j]; v <- tp$var[k, j]
        add(g, sprintf("cumulative incidence at %s", names(gt$cif_event_at)[j]), gt$cif_event_at[[j]], e, e - q * sqrt(v), e + q * sqrt(v)) }
    }
  }
  if (!length(rows)) return(NULL)
  out <- do.call(rbind, rows)
  # a median not reached has no estimate but an open interval [lcl, Inf): covered when the truth lies in it
  out$covered <- !is.na(out$ci_low) & out$true >= out$ci_low & (is.na(out$ci_high) | out$true <= out$ci_high) & (!is.na(out$estimate) | is.na(out$ci_high))
  out
}
