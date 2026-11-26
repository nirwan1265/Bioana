# Regression next to the outcome model (rules inst/rules/regression_rules.json): when analysis.model has measured
# predictors (role covariate: covariate, covariate2, ...). Linear (lm): coefficients with CIs and standardized
# betas, R2 / adjusted R2 / F, Type II tests, residual diagnostics, VIF, influence and a refit without the most
# influential unit, HC3 robust SEs, interactions (nested F, simple slopes), ANCOVA (equal slopes, adjusted
# means), confidence and prediction intervals, and the slope the design could detect (sensitivity, not observed
# power). Logistic (glm binomial): odds ratios with profile CIs, likelihood-ratio tests, pseudo R2,
# Hosmer-Lemeshow, AUC, Box-Tidwell, events per variable, separation, predicted probabilities.
# Checks are reported, never used to switch the model. Numbers match base R, car, lmtest, sandwich, emmeans and
# the datanovia regression tutorials (tests/testthat/test-regression.R).

.ba_has_covariates <- function(p) any(vapply(p$vars, function(v) identical(p$role[[v]], "covariate"), logical(1)))

.ba_regression <- function(fit, frame, p, design, alpha = 0.05, type = "continuous", args = list()) {
  rr <- .ba_rules("regression_rules")
  ra <- utils::modifyList(design$spec$analysis$regression %||% list(), args$regression %||% list())
  if (!inherits(fit, c("lm", "glm")) || inherits(fit, "negbin")) {
    return(list(summary = list(status = "not run", note = sprintf("regression diagnostics are set up for lm and logistic glm fits, not %s", class(fit)[1]))))
  }
  logistic <- inherits(fit, "glm") && identical(stats::family(fit)$family, "binomial")
  if (inherits(fit, "glm") && !logistic) return(list(summary = list(status = "not run", note = "regression diagnostics: glm families other than binomial are not set up")))
  label <- .ba_reg_labels(p)
  mf <- stats::model.frame(fit); n <- nrow(mf)
  X <- stats::model.matrix(fit); k <- ncol(X) - ("(Intercept)" %in% colnames(X))
  notes <- character(); checks <- list()
  chk <- function(check, status, detail, value = NA_real_) checks[[length(checks) + 1L]] <<- list(check = check, status = status, detail = detail, value = unname(value))
  coefs <- if (logistic) .ba_reg_coef_logistic(fit, alpha) else .ba_reg_coef_linear(fit, alpha)
  coefs$predictor <- .ba_reg_term_label(coefs$term, label)
  bl <- p$vars[vapply(p$vars, function(v) identical(p$role[[v]], "block"), logical(1))]
  coefs$block <- vapply(coefs$term, function(t) any(startsWith(t, bl)) && !grepl(":", t), logical(1))
  # fit statistics
  fitstats <- if (logistic) .ba_reg_fit_logistic(fit) else .ba_reg_fit_linear(fit)
  # Type II tests per model term (car::Anova; F for lm, likelihood ratio for glm)
  a2 <- if (requireNamespace("car", quietly = TRUE)) tryCatch({
    a <- car::Anova(fit, type = 2)
    data.frame(term = rownames(a), predictor = .ba_reg_term_label(rownames(a), label), df = a$Df,
               statistic = if (logistic) a[["LR Chisq"]] else a[["F value"]], statistic_name = if (logistic) "LR chi-squared" else "F value",
               p = a[[ncol(a)]], stringsAsFactors = FALSE, row.names = NULL)[rownames(a) != "Residuals", , drop = FALSE]
  }, error = function(e) NULL)
  # influence: standardized residuals, leverage, Cook's distance (thresholds from the rules)
  rs <- stats::rstandard(fit); h <- stats::hatvalues(fit); cd <- stats::cooks.distance(fit)
  rows <- as.integer(names(rs) %||% seq_along(rs))
  lev_cut <- 2 * (k + 1) / n; cook_cut <- (rr$cooks_numerator %||% 4) / (n - k - 1)
  infl <- data.frame(row = rows, std_residual = unname(rs), leverage = unname(h), cooks_distance = unname(cd), stringsAsFactors = FALSE)
  infl$outlier <- abs(infl$std_residual) > (rr$std_residual_cut %||% 3)
  infl$high_leverage <- infl$leverage > lev_cut
  infl$influential <- infl$cooks_distance > cook_cut
  no <- sum(infl$outlier); nl <- sum(infl$high_leverage); ni <- sum(infl$influential)
  chk("outliers", if (no) "warn" else "ok", sprintf("%d units with |standardized residual| > %g%s", no, rr$std_residual_cut %||% 3,
      if (no) paste0(" (rows ", paste(utils::head(infl$row[infl$outlier], 5), collapse = ", "), ")") else ""), no)
  chk("leverage", if (nl) "note" else "ok", sprintf("%d units with leverage > 2(p+1)/n = %.3g (max %.3g)", nl, lev_cut, max(h)), max(h))
  chk("influence", if (any(cd > (rr$cooks_large %||% 0.5))) "warn" else if (ni) "note" else "ok",
      sprintf("%d units with Cook's distance > 4/(n-p-1) = %.3g; max %.3g%s", ni, cook_cut, max(cd), if (any(cd > (rr$cooks_large %||% 0.5))) sprintf("; %d above %g", sum(cd > (rr$cooks_large %||% 0.5)), rr$cooks_large %||% 0.5) else ""), max(cd))
  # refit without the most influential unit: change of each slope (the page's sensitivity check)
  sens <- NULL
  if (n > k + 3) {
    i <- which.max(cd)
    f2 <- tryCatch(suppressWarnings(.ba_reg_refit(fit, stats::formula(fit), mf[-i, , drop = FALSE])), error = function(e) NULL)
    if (!is.null(f2)) {
      b1 <- stats::coef(fit); b2 <- stats::coef(f2)[names(b1)]; tt <- setdiff(names(b1), "(Intercept)")
      se1 <- sqrt(diag(stats::vcov(fit)))[tt]
      sens <- data.frame(term = tt, predictor = .ba_reg_term_label(tt, label), estimate = unname(b1[tt]), without_row = unname(b2[tt]),
                         change_pct = 100 * unname((b2[tt] - b1[tt]) / b1[tt]), change_se = unname((b2[tt] - b1[tt]) / se1), dropped_row = rows[i], stringsAsFactors = FALSE)
      # a shift is judged in standard errors (a % change of a coefficient near 0 means nothing)
      big <- sens[abs(sens$change_se) > (rr$influence_change_se %||% 0.5), , drop = FALSE]
      chk("influence_refit", if (nrow(big)) "warn" else "ok", sprintf("without row %d (largest Cook's distance): %s", rows[i],
          if (nrow(big)) paste(sprintf("%s %+.1f%% (%+.2f SE)", big$predictor, big$change_pct, big$change_se), collapse = ", ") else
            sprintf("every coefficient moves less than %g SE (largest %s %+.1f%%)", rr$influence_change_se %||% 0.5, sens$predictor[which.max(abs(sens$change_se))], sens$change_pct[which.max(abs(sens$change_se))])),
          max(abs(sens$change_se)))
    }
  }
  # collinearity (car::vif; GVIF^(1/(2 df)) squared for terms with several coefficients)
  vif <- .ba_reg_vif(fit, label)
  if (!is.null(vif)) {
    hi <- vif[vif$vif_comparable > (rr$vif_cut %||% 5), , drop = FALSE]
    chk("collinearity", if (any(vif$vif_comparable > (rr$vif_severe %||% 10))) "warn" else if (nrow(hi)) "warn" else "ok",
        sprintf("largest VIF %.2f (%s)%s", max(vif$vif_comparable), vif$predictor[which.max(vif$vif_comparable)],
                if (nrow(hi)) sprintf("; above %g: %s (standard errors inflated about %.1fx; drop or combine one, or use ridge / PCA)", rr$vif_cut %||% 5,
                                      paste(hi$predictor, collapse = ", "), sqrt(max(hi$vif_comparable))) else ""), max(vif$vif_comparable))
  }
  extra <- list()
  if (!logistic) {
    e <- stats::residuals(fit)
    if (length(e) >= 3 && length(e) <= 5000) { sw <- stats::shapiro.test(e)
      chk("normality", if (sw$p.value < alpha) "note" else "ok", sprintf("Shapiro-Wilk on residuals W = %.3f, p = %.3g%s", sw$statistic, sw$p.value,
          if (sw$p.value < alpha) if (n >= (rr$large_n %||% 100)) " (large n: the tests are fairly robust)" else "; transform the outcome (log, sqrt)" else ""), sw$p.value) }
    if (requireNamespace("lmtest", quietly = TRUE) && k >= 1) {
      bp <- lmtest::bptest(fit)
      chk("homoscedasticity", if (bp$p.value < alpha) "warn" else "ok", sprintf("Breusch-Pagan BP = %.3g, df = %d, p = %.3g%s", bp$statistic, as.integer(bp$parameter), bp$p.value,
          if (bp$p.value < alpha) "; HC3 robust standard errors are in the coefficient table (or log the outcome)" else ""), bp$p.value)
      cv_terms <- intersect(colnames(X), p$vars[vapply(p$vars, function(v) identical(p$role[[v]], "covariate"), logical(1))])
      if (length(cv_terms)) { rt <- tryCatch(lmtest::resettest(fit, power = 2:3, type = "fitted"), error = function(e) NULL)
        if (!is.null(rt)) chk("linearity", if (rt$p.value < alpha) "warn" else "ok", sprintf("RESET (fitted^2, fitted^3) F = %.3g, p = %.3g%s", rt$statistic, rt$p.value,
            if (rt$p.value < alpha) "; curvature: try log(x) or poly(x, 2) (analysis.covariates transform)" else ""), rt$p.value) }
      oc <- ra$order_column
      if (!is.null(oc) && oc %in% names(frame)) { o <- order(suppressWarnings(as.numeric(frame[[oc]][!is.na(frame$y_)])))
        dw <- lmtest::dwtest(stats::lm(stats::formula(fit), data = mf[o, , drop = FALSE]))
        chk("independence", if (dw$p.value < alpha) "warn" else "ok", sprintf("Durbin-Watson DW = %.3f, p = %.3g (order: %s)", dw$statistic, dw$p.value, oc), dw$p.value) }
    }
    extra$sensitivity <- .ba_reg_sensitivity(fit, label, alpha, rr)
  } else {
    lg <- .ba_reg_logistic_checks(fit, mf, p, label, alpha, rr)
    for (cc in lg$checks) chk(cc$check, cc$status, cc$detail, cc$value)
    extra$classification <- lg$classification; extra$hosmer_lemeshow <- lg$hl; extra$auc <- lg$auc
  }
  inter <- .ba_reg_interactions(fit, mf, p, label, alpha, logistic)
  anc <- .ba_reg_ancova(fit, mf, p, label, alpha, logistic)
  if (!is.null(anc$check)) chk(anc$check$check, anc$check$status, anc$check$detail, anc$check$value)
  pred <- .ba_reg_predict(fit, mf, p, ra$predict_at, label, alpha, logistic)
  robust <- NULL
  if (!logistic && identical(ra$robust, "mm") && requireNamespace("robustbase", quietly = TRUE)) {
    rb <- tryCatch(robustbase::lmrob(stats::formula(fit), data = mf), error = function(e) NULL)
    if (!is.null(rb)) { cr <- summary(rb)$coefficients
      robust <- data.frame(term = rownames(cr), predictor = .ba_reg_term_label(rownames(cr), label), estimate = cr[, 1], se = cr[, 2], statistic = cr[, 3], p = cr[, 4],
                           method = "MM (robustbase::lmrob)", stringsAsFactors = FALSE, row.names = NULL) }
  }
  st <- vapply(checks, `[[`, "", "status")
  list(summary = c(list(status = "run", model = if (logistic) "logistic regression (glm binomial)" else "linear regression (lm)",
                        package = "stats", version = as.character(utils::packageVersion("stats")), n = n, n_coefficients = k, fit = fitstats,
                        predictors = as.list(label), checks = checks, n_warn = sum(st == "warn"),
                        coefficients = lapply(seq_len(nrow(coefs)), function(i) as.list(coefs[i, intersect(c("term", "predictor", "estimate", "se", "ci_low", "ci_high", "p", "std_beta", "p_hc3", "odds_ratio", "or_ci_low", "or_ci_high", "block"), names(coefs))])),
                        type2_tests = if (!is.null(a2)) lapply(seq_len(nrow(a2)), function(i) as.list(a2[i, ])),
                        interactions = inter$summary, ancova = anc$summary, predictions = if (!is.null(pred)) nrow(pred),
                        reference = rr$reference, tutorial = "docs/engines/regression.md"), extra),
       coefficients = coefs, type2 = a2, influence = infl, influence_refit = sens, vif = vif, simple_slopes = inter$table,
       adjusted_means = anc$table, predictions = pred, robust = robust)
}

# term -> readable name: covariate terms by analysis.covariates name, other variables by their column
.ba_reg_labels <- function(p) {
  stats::setNames(vapply(p$vars, function(v) as.character(p$values[[v]]$name %||% p$col[[v]] %||% v), ""), p$vars)
}
.ba_reg_term_label <- function(terms, label) {
  vapply(terms, function(t) {
    if (t %in% c("(Intercept)", "Residuals")) return(t)
    parts <- strsplit(t, ":", fixed = TRUE)[[1]]
    paste(vapply(parts, function(x) { v <- names(label)[startsWith(x, names(label))]
      if (!length(v)) return(x); v <- v[which.max(nchar(v))]; paste0(label[[v]], substring(x, nchar(v) + 1)) }, ""), collapse = ":")
  }, "", USE.NAMES = FALSE)
}

.ba_reg_coef_linear <- function(fit, alpha) {
  s <- summary(fit)$coefficients; ci <- stats::confint(fit, level = 1 - alpha)
  out <- data.frame(term = rownames(s), estimate = s[, 1], se = s[, 2], statistic = s[, 3], df = fit$df.residual, p = s[, 4],
                    ci_low = ci[rownames(s), 1], ci_high = ci[rownames(s), 2], stringsAsFactors = FALSE, row.names = NULL)
  # standardized coefficient: b x sd(column of X) / sd(y) (numeric predictors; dummies are "partially standardized")
  X <- stats::model.matrix(fit); y <- stats::model.response(stats::model.frame(fit))
  out$std_beta <- ifelse(out$term == "(Intercept)", NA_real_, out$estimate * apply(X[, out$term, drop = FALSE], 2, stats::sd) / stats::sd(y))
  if (requireNamespace("sandwich", quietly = TRUE)) {
    v <- sandwich::vcovHC(fit, type = "HC3"); se <- sqrt(diag(v))[out$term]; q <- stats::qt(1 - alpha / 2, fit$df.residual)
    out$se_hc3 <- unname(se); out$p_hc3 <- 2 * stats::pt(-abs(out$estimate / se), fit$df.residual)
    out$ci_low_hc3 <- out$estimate - q * se; out$ci_high_hc3 <- out$estimate + q * se
  }
  out
}

.ba_reg_coef_logistic <- function(fit, alpha) {
  s <- summary(fit)$coefficients
  # profile-likelihood intervals (stats / MASS confint for glm); Wald when profiling fails
  ci <- tryCatch(suppressMessages(stats::confint(fit, level = 1 - alpha)), error = function(e) NULL)
  w <- stats::confint.default(fit, level = 1 - alpha)
  if (is.null(dim(ci))) ci <- matrix(ci, nrow = 1, dimnames = list(rownames(s), NULL))
  data.frame(term = rownames(s), estimate = s[, 1], se = s[, 2], statistic = s[, 3], p = s[, 4],
             odds_ratio = exp(s[, 1]), or_ci_low = exp(if (!is.null(ci)) ci[rownames(s), 1] else w[rownames(s), 1]),
             or_ci_high = exp(if (!is.null(ci)) ci[rownames(s), 2] else w[rownames(s), 2]),
             ci_method = if (!is.null(ci)) "profile likelihood" else "Wald", or_wald_low = exp(w[rownames(s), 1]), or_wald_high = exp(w[rownames(s), 2]),
             stringsAsFactors = FALSE, row.names = NULL)
}

.ba_reg_fit_linear <- function(fit) {
  s <- summary(fit); f <- s$fstatistic
  list(r_squared = s$r.squared, adj_r_squared = s$adj.r.squared, sigma = s$sigma, df_residual = fit$df.residual,
       F = if (!is.null(f)) unname(f[1]) else NA_real_, df1 = if (!is.null(f)) unname(f[2]) else NA_real_, df2 = if (!is.null(f)) unname(f[3]) else NA_real_,
       p = if (!is.null(f)) unname(stats::pf(f[1], f[2], f[3], lower.tail = FALSE)) else NA_real_, AIC = stats::AIC(fit), BIC = stats::BIC(fit),
       f2 = s$r.squared / (1 - s$r.squared))
}

.ba_reg_fit_logistic <- function(fit) {
  n <- length(fit$y); ll <- as.numeric(stats::logLik(fit))
  ll0 <- as.numeric(stats::logLik(.ba_reg_refit(fit, stats::update(stats::formula(fit), . ~ 1), stats::model.frame(fit))))
  lr <- fit$null.deviance - fit$deviance; df <- fit$df.null - fit$df.residual
  cs <- 1 - exp(2 * (ll0 - ll) / n)
  list(n = n, events = sum(fit$y == 1), null_deviance = fit$null.deviance, deviance = fit$deviance, lr_chisq = lr, df = df,
       p = stats::pchisq(lr, df, lower.tail = FALSE), AIC = stats::AIC(fit), BIC = stats::BIC(fit),
       mcfadden_r2 = 1 - ll / ll0, cox_snell_r2 = cs, nagelkerke_r2 = cs / (1 - exp(2 * ll0 / n)))
}

.ba_reg_vif <- function(fit, label) {
  if (!requireNamespace("car", quietly = TRUE)) return(NULL)
  nterm <- length(attr(stats::terms(fit), "term.labels"))
  if (nterm < 2) return(NULL)
  # with interactions, VIF per predictor (car type "predictor": a predictor with its interactions); else per term
  ty <- if (any(grepl(":", attr(stats::terms(fit), "term.labels")))) "predictor" else "terms"
  v <- tryCatch(suppressMessages(suppressWarnings(car::vif(fit, type = ty))), error = function(e) tryCatch(suppressMessages(suppressWarnings(car::vif(fit))), error = function(e2) NULL))
  if (is.data.frame(v)) v <- as.matrix(v[, intersect(c("GVIF", "Df", "GVIF^(1/(2*Df))"), names(v)), drop = FALSE])
  if (is.null(v)) return(NULL)
  if (is.null(dim(v))) return(data.frame(term = names(v), predictor = .ba_reg_term_label(names(v), label), vif = unname(v), df = 1, gvif_adj = sqrt(unname(v)),
                                         vif_comparable = unname(v), stringsAsFactors = FALSE))
  ga <- v[, ncol(v)]
  data.frame(term = rownames(v), predictor = .ba_reg_term_label(rownames(v), label), vif = v[, 1], df = v[, 2], gvif_adj = ga,
             vif_comparable = ga^2, stringsAsFactors = FALSE, row.names = NULL)
}

# Sensitivity (datanovia power page: report what the design could detect, never observed power): the smallest
# partial f2 detectable with power `power` at alpha given the residual df (pwr::pwr.f2.test, u = 1), and the
# slope it means for each numeric predictor: b = sqrt(f2 x RSS / SS(x | other predictors)).
.ba_reg_sensitivity <- function(fit, label, alpha, rr) {
  if (!requireNamespace("pwr", quietly = TRUE) || fit$df.residual < 2) return(NULL)
  pw <- rr$sensitivity_power %||% 0.8
  f2 <- tryCatch(pwr::pwr.f2.test(u = 1, v = fit$df.residual, sig.level = alpha, power = pw)$f2, error = function(e) NA_real_)
  X <- stats::model.matrix(fit); rss <- sum(stats::residuals(fit)^2)
  num <- setdiff(colnames(X), "(Intercept)")
  rows <- lapply(num, function(j) {
    o <- setdiff(colnames(X), j)
    ssx <- if (length(o)) sum(stats::lm.fit(X[, o, drop = FALSE], X[, j])$residuals^2) else sum((X[, j] - mean(X[, j]))^2)
    data.frame(term = j, predictor = .ba_reg_term_label(j, label), detectable_abs_slope = sqrt(f2 * rss / ssx), stringsAsFactors = FALSE)
  })
  list(power = pw, alpha = alpha, detectable_partial_f2 = f2, cohen = as.character(cut(f2, c(-Inf, 0.02, 0.15, 0.35, Inf), c("below small", "small", "medium", "large"), right = FALSE)),
       slopes = do.call(rbind, rows), note = "sensitivity: the effect this sample could detect with the stated power, not the power of the observed effect")
}

# Interactions with a covariate: nested F (or LR) test of each interaction term, and simple slopes of the
# covariate at the other variable's mean -1 SD / mean / +1 SD (numeric) or at each level (factor), emmeans::emtrends.
.ba_reg_interactions <- function(fit, mf, p, label, alpha, logistic) {
  tl <- attr(stats::terms(fit), "term.labels")
  it <- tl[grepl(":", tl)]
  cvars <- p$vars[vapply(p$vars, function(v) identical(p$role[[v]], "covariate"), logical(1))]
  it <- it[vapply(it, function(t) any(strsplit(t, ":", fixed = TRUE)[[1]] %in% cvars), logical(1))]
  if (!length(it)) return(list())
  tests <- list(); slopes <- list()
  for (t in it) {
    red <- .ba_reg_refit(fit, stats::update(stats::formula(fit), stats::as.formula(paste(". ~ . -", t))), mf)
    a <- stats::anova(red, fit, test = if (logistic) "LRT" else "F")
    pv <- a[2, ncol(a)]
    tests[[t]] <- list(term = t, predictor = .ba_reg_term_label(t, label), statistic = if (logistic) a[2, "Deviance"] else a[2, "F"],
                       statistic_name = if (logistic) "LR chi-squared" else "F", df = a[2, "Df"], p = pv)
    vs <- strsplit(t, ":", fixed = TRUE)[[1]]; if (length(vs) != 2L || !requireNamespace("emmeans", quietly = TRUE)) next
    x <- intersect(vs, cvars)[1]; mod <- setdiff(vs, x)[1]
    at <- if (is.numeric(mf[[mod]])) { m <- mean(mf[[mod]]); s <- stats::sd(mf[[mod]]); stats::setNames(list(c(m - s, m, m + s)), mod) } else list()
    et <- tryCatch(as.data.frame(summary(emmeans::emtrends(fit, stats::as.formula(paste("~", mod)), var = x, at = at), infer = c(TRUE, TRUE), level = 1 - alpha)),
                   error = function(e) NULL)
    if (is.null(et)) next
    tr <- grep("\\.trend$", names(et), value = TRUE)[1]
    slopes[[t]] <- data.frame(interaction = .ba_reg_term_label(t, label), slope_of = label[[x]], at_variable = label[[mod]],
                              at = if (is.numeric(mf[[mod]])) c("mean - 1 SD", "mean", "mean + 1 SD") else as.character(et[[mod]]),
                              at_value = if (is.numeric(mf[[mod]])) et[[mod]] else NA_real_, slope = et[[tr]], se = et$SE,
                              ci_low = et[[grep("lower|asymp.LCL", names(et), value = TRUE)[1]]], ci_high = et[[grep("upper|asymp.UCL", names(et), value = TRUE)[1]]],
                              p = et$p.value, scale = if (logistic) "log odds per unit" else "outcome per unit", stringsAsFactors = FALSE)
  }
  list(summary = unname(tests), table = if (length(slopes)) do.call(rbind, slopes))
}

# ANCOVA (a treatment factor and a covariate, no interaction between them): homogeneity of regression slopes
# (the treatment x covariate term added, nested test) and covariate-adjusted means (emmeans at the covariate mean).
.ba_reg_ancova <- function(fit, mf, p, label, alpha, logistic) {
  tl <- attr(stats::terms(fit), "term.labels")
  tv <- intersect(tl, p$vars[vapply(p$vars, function(v) identical(p$role[[v]], "treatment"), logical(1))])
  tv <- tv[vapply(tv, function(v) is.factor(mf[[v]]), logical(1))]
  cv <- intersect(tl, p$vars[vapply(p$vars, function(v) identical(p$role[[v]], "covariate"), logical(1))])
  if (!length(tv) || !length(cv)) return(list())
  has_int <- any(vapply(tl[grepl(":", tl)], function(t) any(strsplit(t, ":", fixed = TRUE)[[1]] %in% tv) && any(strsplit(t, ":", fixed = TRUE)[[1]] %in% cv), logical(1)))
  out <- list()
  if (!has_int) {
    full <- tryCatch(.ba_reg_refit(fit, stats::update(stats::formula(fit), stats::as.formula(paste(". ~ . +", paste(outer(tv, cv, paste, sep = ":"), collapse = " + ")))), mf), error = function(e) NULL)
    if (!is.null(full)) {
      a <- stats::anova(fit, full, test = if (logistic) "LRT" else "F"); pv <- a[2, ncol(a)]
      out$check <- list(check = "equal_slopes", status = if (isTRUE(pv < alpha)) "warn" else "ok", value = pv,
                        detail = sprintf("homogeneity of regression slopes (%s x %s): p = %.3g%s", paste(label[tv], collapse = ", "), paste(label[cv], collapse = ", "), pv,
                                         if (isTRUE(pv < alpha)) "; slopes differ between groups: adjusted means depend on the covariate value (fit the interaction)" else ""))
    }
  }
  if (requireNamespace("emmeans", quietly = TRUE)) {
    em <- tryCatch(suppressMessages(emmeans::emmeans(fit, stats::as.formula(paste("~", tv[1])), type = if (logistic) "response" else NULL)), error = function(e) NULL)
    if (!is.null(em)) {
      d <- as.data.frame(summary(em, level = 1 - alpha))
      est <- intersect(c("emmean", "prob", "response"), names(d))[1]
      out$table <- data.frame(factor = label[[tv[1]]], level = as.character(d[[tv[1]]]), adjusted_mean = d[[est]], se = d$SE,
                              ci_low = d[[grep("lower|LCL", names(d), value = TRUE)[1]]], ci_high = d[[grep("upper|UCL", names(d), value = TRUE)[1]]],
                              at = paste(sprintf("%s = %.4g", label[cv], vapply(cv, function(v) mean(mf[[v]]), 0)), collapse = ", "), stringsAsFactors = FALSE)
      out$summary <- list(factor = label[[tv[1]]], covariates = as.list(unname(label[cv])), interaction_in_model = has_int,
                          equal_slopes_p = out$check$value)
    }
  }
  out
}

# Predictions at analysis.regression.predict_at (a list of {term or covariate name or column: value}):
# lm: fit with confidence and prediction intervals; logistic: probability with a CI from the link scale.
.ba_reg_predict <- function(fit, mf, p, at, label, alpha, logistic) {
  if (!length(at)) return(NULL)
  if (!is.null(names(at))) at <- list(at)
  key <- function(nm) { v <- p$vars[p$vars == nm | label[p$vars] == nm | vapply(p$vars, function(z) identical(as.character(p$col[[z]]), nm), logical(1))]; if (length(v)) v[1] else nm }
  nd <- do.call(rbind, lapply(at, function(a) { a <- stats::setNames(a, vapply(names(a), key, "")); as.data.frame(a, stringsAsFactors = FALSE) }))
  vars <- intersect(names(mf)[-1], all.vars(stats::formula(fit)))
  # values are given on the measured scale; a transformed covariate (analysis.covariates transform) is fitted on
  # the transformed scale, so the value is transformed before predicting
  raw <- nd
  for (v in intersect(names(nd), vars)) { tf <- p$values[[v]]$transform; if (is.function(tf)) nd[[v]] <- tf(as.numeric(nd[[v]])) }
  for (v in vars) {
    if (!v %in% names(nd)) nd[[v]] <- if (is.numeric(mf[[v]])) mean(mf[[v]]) else levels(mf[[v]])[1]
    if (is.factor(mf[[v]])) nd[[v]] <- factor(as.character(nd[[v]]), levels = levels(mf[[v]]))
    else nd[[v]] <- as.numeric(nd[[v]])
  }
  out <- nd[vars]
  sc <- vapply(vars, function(v) if (is.function(p$values[[v]]$transform)) sprintf(" (%s scale)", p$values[[v]]$scale) else "", "")
  names(out) <- paste0(label[vars] %||% vars, sc)
  for (v in vars[nzchar(sc)]) out[[label[[v]]]] <- if (v %in% names(raw)) as.numeric(raw[[v]]) else NA_real_   # the measured value, when given
  if (!logistic) {
    cf <- stats::predict(fit, nd, interval = "confidence", level = 1 - alpha); pr <- stats::predict(fit, nd, interval = "prediction", level = 1 - alpha)
    out$fit <- cf[, "fit"]; out$ci_low <- cf[, "lwr"]; out$ci_high <- cf[, "upr"]; out$pi_low <- pr[, "lwr"]; out$pi_high <- pr[, "upr"]
  } else {
    lp <- stats::predict(fit, nd, type = "link", se.fit = TRUE); q <- stats::qnorm(1 - alpha / 2)
    out$probability <- stats::plogis(lp$fit); out$ci_low <- stats::plogis(lp$fit - q * lp$se.fit); out$ci_high <- stats::plogis(lp$fit + q * lp$se.fit)
  }
  num <- vars[vapply(vars, function(v) is.numeric(mf[[v]]), logical(1))]
  out$extrapolation <- if (length(num)) Reduce(`|`, lapply(num, function(v) nd[[v]] < min(mf[[v]]) | nd[[v]] > max(mf[[v]]))) else FALSE
  out
}

.ba_reg_logistic_checks <- function(fit, mf, p, label, alpha, rr) {
  checks <- list(); add <- function(check, status, detail, value = NA_real_) checks[[length(checks) + 1L]] <<- list(check = check, status = status, detail = detail, value = value)
  y <- fit$y; pr <- stats::fitted(fit); k <- length(stats::coef(fit)) - 1
  ev <- min(sum(y == 1), sum(y == 0)); epv <- ev / max(1, k)
  add("events_per_variable", if (epv < (rr$min_epv %||% 10)) "warn" else "ok", sprintf("%d events in the smaller outcome class for %d coefficients: %.1f per coefficient%s", ev, k, epv,
      if (epv < (rr$min_epv %||% 10)) sprintf(" (< %g: estimates may be biased; fewer predictors or Firth's penalized likelihood)", rr$min_epv %||% 10) else ""), epv)
  sep <- any(pr < 1e-8 | pr > 1 - 1e-8) || any(abs(stats::coef(fit)[-1]) > (rr$separation_abs_coef %||% 10), na.rm = TRUE)
  add("separation", if (sep) "warn" else "ok", if (sep) "fitted probabilities of 0 or 1 or a huge coefficient: (quasi-)complete separation; the ML estimates and Wald tests are unreliable (Firth's penalized likelihood)" else "no sign of separation", as.numeric(sep))
  # Hosmer-Lemeshow (g groups by fitted deciles)
  g <- rr$hosmer_lemeshow_groups %||% 10
  br <- unique(stats::quantile(pr, seq(0, 1, length.out = g + 1)))
  hl <- NULL
  if (length(br) > 3) {
    gr <- cut(pr, br, include.lowest = TRUE)
    o1 <- tapply(y, gr, sum); e1 <- tapply(pr, gr, sum); nn <- tapply(y, gr, length)
    x2 <- sum((o1 - e1)^2 / (e1 * (1 - e1 / nn))); dfh <- length(nn) - 2
    hl <- list(statistic = x2, df = dfh, p = stats::pchisq(x2, dfh, lower.tail = FALSE), groups = length(nn))
    add("calibration", if (hl$p < alpha) "warn" else "ok", sprintf("Hosmer-Lemeshow chi-squared(%d) = %.3g, p = %.3g", dfh, x2, hl$p), hl$p)
  }
  # discrimination: AUC = P(fitted of a random event > fitted of a random non-event) (Mann-Whitney); DeLong CI via pROC
  r <- rank(pr); n1 <- sum(y == 1); n0 <- sum(y == 0)
  auc <- (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
  auc_ci <- if (requireNamespace("pROC", quietly = TRUE)) tryCatch(as.numeric(pROC::ci.auc(pROC::roc(y, pr, quiet = TRUE, direction = "<", levels = c(0, 1))))[c(1, 3)], error = function(e) c(NA, NA)) else c(NA, NA)
  add("discrimination", "info", sprintf("AUC %.3f%s", auc, if (!is.na(auc_ci[1])) sprintf(" [%.3f, %.3f] (DeLong)", auc_ci[1], auc_ci[2]) else ""), auc)
  cl <- table(factor(y, c(0, 1)), factor(as.integer(pr >= 0.5), c(0, 1)))
  classification <- list(threshold = 0.5, accuracy = sum(diag(cl)) / sum(cl), sensitivity = cl[2, 2] / sum(cl[2, ]), specificity = cl[1, 1] / sum(cl[1, ]),
                         note = "in-sample, threshold 0.5: optimistic; prediction accuracy belongs to cross-validation")
  # linearity of the logit: Box-Tidwell, x * log(x) added for each positive numeric predictor (likelihood ratio)
  num <- names(mf)[-1][vapply(names(mf)[-1], function(v) is.numeric(mf[[v]]) && all(mf[[v]] > 0), logical(1))]
  for (v in num) {
    mf2 <- mf; mf2$bt_ <- mf[[v]] * log(mf[[v]])
    f2 <- tryCatch(suppressWarnings(stats::glm(stats::update(stats::formula(fit), . ~ . + bt_), data = mf2, family = stats::binomial())), error = function(e) NULL)
    if (is.null(f2)) next
    pv <- stats::pchisq(fit$deviance - f2$deviance, 1, lower.tail = FALSE)
    add(paste0("linearity_logit_", v), if (pv < alpha) "warn" else "ok", sprintf("Box-Tidwell %s: p = %.3g%s", label[[v]] %||% v, pv,
        if (pv < alpha) "; the log odds bend with this predictor: log(x), poly(x, 2) or a spline" else ""), pv)
  }
  list(checks = checks, hl = hl, auc = list(auc = auc, ci_low = auc_ci[1], ci_high = auc_ci[2]), classification = classification)
}

# Refit with another formula or data: lm, or glm with the fit's family (update() would re-evaluate the engine's
# call, whose formula and data live inside the engine function).
.ba_reg_refit <- function(fit, formula, data) {
  if (inherits(fit, "glm")) stats::glm(formula, data = data, family = stats::family(fit)) else stats::lm(formula, data = data)
}
