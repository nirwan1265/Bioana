# Agreement between raters or methods on the same units (data_type "ratings"; rules inst/rules/agreement_rules.json).
# Nominal: Cohen's kappa (2 raters, vcd::Kappa), Fleiss' kappa (3+, irr::kappam.fleiss; Light's kappa = mean pairwise
# Cohen when irr is missing). Ordinal: weighted kappa (linear and quadratic, vcd::Kappa). Continuous: ICC
# (psych::ICC, all six Shrout-Fleiss forms; the primary form from the rules).

.ba_method_ratings <- function(entry, design, outdir, alpha = 0.05, args = list()) {
  rules <- .ba_rules("agreement_rules")
  ag <- utils::modifyList(design$spec$analysis$agreement %||% list(), args)
  dat <- if (!is.null(entry$files$table)) .ba_read_table_any(file.path(design$dir, entry$files$table)) else design$meta
  W <- .ba_ratings_wide(dat, entry, ag)
  scale <- match.arg(entry$scale %||% ag$scale %||% "nominal", c("nominal", "ordinal", "continuous"))
  k <- ncol(W); notes <- character()
  W <- W[stats::complete.cases(W), , drop = FALSE]
  if (nrow(W) < 2) stop("Agreement needs at least two units rated by every rater.", call. = FALSE)
  rows <- list(); tables <- NULL
  add <- function(name, est, lo = NA_real_, hi = NA_real_, z = NA_real_, pv = NA_real_, scale_name = NULL)
    rows[[length(rows) + 1L]] <<- data.frame(measure = name, estimate = est, ci_low = lo, ci_high = hi, z = z, p = pv,
                                             interpretation = .ba_agree_label(if (is.na(lo) || !identical(scale_name, "icc")) est else lo, rules, scale_name %||% "kappa"),
                                             stringsAsFactors = FALSE)
  if (scale == "continuous") {
    if (!requireNamespace("psych", quietly = TRUE)) stop("ICC needs the psych package.", call. = FALSE)
    ic <- suppressMessages(psych::ICC(as.matrix(as.data.frame(lapply(W, as.numeric))), lmer = FALSE))$results
    for (i in seq_len(nrow(ic))) add(sprintf("%s (%s)", ic$type[i], ic$Description[i] %||% ""), ic$ICC[i], ic$`lower bound`[i], ic$`upper bound`[i], NA, ic$p[i], "icc")
    prim <- rules$icc_primary %||% "ICC2"
    primary <- rows[[match(prim, ic$type)]]
    notes <- c(notes, sprintf("primary form %s: %s (stat_rules: agreement_rules.json icc_primary); interpretation from the lower 95%% bound (Koo & Li 2016)", prim, rules$icc_forms[[prim]] %||% ""))
  } else {
    lev <- .ba_rating_levels(W, ag, scale)
    F_ <- as.data.frame(lapply(W, function(z) factor(as.character(z), levels = lev)))
    if (k == 2) {
      tb <- table(F_[[1]], F_[[2]], dnn = names(W))
      tables <- data.frame(rater1 = rep(rownames(tb), ncol(tb)), rater2 = rep(colnames(tb), each = nrow(tb)), count = c(tb), stringsAsFactors = FALSE)
      names(tables)[1:2] <- names(W)
      kp <- .ba_vcd_kappa(tb, "Equal-Spacing"); kq <- .ba_vcd_kappa(tb, "Fleiss-Cohen")
      add("Cohen's kappa (unweighted)", kp$u[1], kp$u[2], kp$u[3], kp$u[4], kp$u[5])
      if (scale == "ordinal") {
        add("weighted kappa (linear, Cicchetti-Allison)", kp$w[1], kp$w[2], kp$w[3], kp$w[4], kp$w[5])
        add("weighted kappa (quadratic, Fleiss-Cohen)", kq$w[1], kq$w[2], kq$w[3], kq$w[4], kq$w[5])
      }
      po <- sum(diag(tb)) / sum(tb)
      primary <- rows[[if (scale == "ordinal") match(rules$ordinal_primary %||% "linear", c("unweighted", "linear", "quadratic")) else 1L]]
      notes <- c(notes, sprintf("observed agreement %.1f%% of %d units", 100 * po, sum(tb)))
    } else {
      R_ <- as.matrix(as.data.frame(lapply(F_, as.character)))
      if (requireNamespace("irr", quietly = TRUE)) {
        # Fleiss (1971) for 3+ raters, overall and per category (irr::kappam.fleiss(detail = TRUE)); Light (1971) via irr
        fk <- irr::kappam.fleiss(R_, detail = TRUE)
        add(sprintf("Fleiss' kappa (%d raters)", k), fk$value, NA, NA, fk$statistic, fk$p.value)
        dt <- fk$detail
        if (!is.null(dt)) for (cg in rownames(dt)) add(sprintf("Fleiss' kappa, category %s", cg), dt[cg, "Kappa"], NA, NA, dt[cg, "z"], dt[cg, "p.value"])
        # irr's variance for Light's kappa can be negative (NaN z and p): kappa kept, test dropped with a note
        lk <- withCallingHandlers(irr::kappam.light(R_), warning = function(w) if (grepl("NaN", conditionMessage(w))) invokeRestart("muffleWarning"))
        z_ok <- is.finite(lk$statistic)
        if (!z_ok) notes <- c(notes, "Light's kappa: irr's variance estimate is negative here, so no z or p")
        add(sprintf("Light's kappa (mean of %d pairwise Cohen's kappas)", choose(k, 2)), lk$value, NA, NA, if (z_ok) lk$statistic else NA, if (z_ok) lk$p.value else NA)
      } else {
        notes <- c(notes, "Fleiss' kappa needs irr (not installed): Light's kappa (mean of pairwise Cohen's kappas, vcd) instead")
        pk <- utils::combn(k, 2, function(ij) .ba_vcd_kappa(table(F_[[ij[1]]], F_[[ij[2]]]), "Equal-Spacing")$u[1])
        add(sprintf("Light's kappa (mean of %d pairwise Cohen's kappas)", length(pk)), mean(pk))
      }
      primary <- rows[[1]]
    }
  }
  tab <- do.call(rbind, rows)
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  paths <- c(agreement = .ba_write_tsv(tab, file.path(outdir, "agreement.tsv")))
  if (!is.null(tables)) paths <- c(paths, agreement_table = .ba_write_tsv(tables, file.path(outdir, "agreement_table.tsv")))
  pk_ <- if (scale == "continuous") "psych" else if (k > 2 && requireNamespace("irr", quietly = TRUE)) "irr" else "vcd"
  summary <- list(data_type = "ratings", scale = scale, raters = as.list(names(W)), n_units = nrow(W), package = pk_, version = as.character(utils::packageVersion(pk_)),
                  primary = as.list(primary), measures = lapply(seq_len(nrow(tab)), function(i) as.list(tab[i, ])),
                  reference = rules$reference, tutorial = "docs/engines/correlation_agreement.md", warnings = notes)
  list(summary = summary, paths = paths, results = tab)
}

# Ratings as a units x raters table: wide columns (entry$columns / analysis.agreement.raters) or long
# (entry$rater_column + entry$score_column + entry$unit_column).
.ba_ratings_wide <- function(dat, entry, ag) {
  rc <- entry$rater_column %||% ag$rater_column
  if (!is.null(rc)) {
    sc <- entry$score_column %||% ag$score_column %||% "score"; uc <- entry$unit_column %||% ag$unit_column %||% "unit_id"
    if (!all(c(rc, sc, uc) %in% names(dat))) stop(sprintf("Long ratings need columns %s, %s, %s.", uc, rc, sc), call. = FALSE)
    w <- stats::reshape(dat[c(uc, rc, sc)], idvar = uc, timevar = rc, direction = "wide")
    out <- w[-1]; names(out) <- sub(paste0("^", sc, "\\."), "", names(out)); return(out)
  }
  cols <- as.character(unlist(entry$columns %||% ag$raters))
  if (length(cols) < 2) stop("Agreement needs at least two rater columns (analysis.agreement.raters).", call. = FALSE)
  miss <- setdiff(cols, names(dat)); if (length(miss)) stop(sprintf("Rater columns not in the table: %s.", paste(miss, collapse = ", ")), call. = FALSE)
  dat[cols]
}

.ba_rating_levels <- function(W, ag, scale) {
  given <- as.character(unlist(ag$levels))
  if (length(given)) return(given)
  v <- unique(as.character(unlist(W))); v <- v[!is.na(v)]
  num <- suppressWarnings(as.numeric(v))
  if (!anyNA(num)) v[order(num)] else if (scale == "ordinal") stop("Ordinal ratings that are not numbers need analysis.agreement.levels (lowest to highest).", call. = FALSE) else sort(v)
}

# vcd::Kappa: estimate, CI (confint), z, p for the unweighted (u) and weighted (w) kappa.
.ba_vcd_kappa <- function(tb, weights) {
  if (!requireNamespace("vcd", quietly = TRUE)) stop("Cohen's kappa needs the vcd package.", call. = FALSE)
  k <- vcd::Kappa(tb, weights = weights); ci <- stats::confint(k)
  f <- function(x, r) { z <- x[["value"]] / x[["ASE"]]; c(x[["value"]], ci[r, 1], ci[r, 2], z, 2 * stats::pnorm(-abs(z))) }
  list(u = f(k$Unweighted, "Unweighted"), w = f(k$Weighted, "Weighted"))
}

.ba_agree_label <- function(x, rules, which) {
  b <- if (which == "icc") rules$icc_scale else rules$kappa_scale
  if (is.null(b) || !is.finite(x)) return(NA_character_)
  as.character(cut(x, c(-Inf, unlist(b$breaks), Inf), unlist(b$labels), right = FALSE))
}
