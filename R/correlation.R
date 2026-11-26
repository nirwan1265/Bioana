# Correlation between measured traits (data_type "traits"; rules inst/rules/correlation_rules.json).
# Every pair: stats::cor.test (pearson / spearman / kendall) with Holm-adjusted p; the same on the residuals of
# the design's fixed part (treatments, blocks), because across plots a correlation also carries the treatment
# differences (two traits that both respond to nitrogen correlate even with no link inside a treatment).

.ba_method_traits <- function(entry, design, outdir, alpha = 0.05, args = list()) {
  rules <- .ba_rules("correlation_rules")
  cs <- utils::modifyList(design$spec$analysis$correlation %||% list(), args)
  vars <- unique(as.character(unlist(entry$columns %||% cs$variables)))
  dat <- if (!is.null(entry$files$table)) .ba_read_table_any(file.path(design$dir, entry$files$table)) else design$meta
  miss <- setdiff(vars, names(dat))
  if (length(vars) < 2) stop("Correlation needs at least two variables (analysis.correlation.variables).", call. = FALSE)
  if (length(miss)) stop(sprintf("Correlation variables not in the table: %s.", paste(miss, collapse = ", ")), call. = FALSE)
  method <- match.arg(cs$method %||% rules$method %||% "pearson", c("pearson", "spearman", "kendall"))
  adj <- cs$adjust %||% rules$adjust %||% "holm"
  X <- as.data.frame(lapply(dat[vars], function(z) suppressWarnings(as.numeric(z))))
  notes <- character()
  one <- function(x, y) withCallingHandlers(stats::cor.test(x, y, method = method),
    warning = function(w) { if (grepl("ties|exact", conditionMessage(w))) { notes <<- unique(c(notes, "ties: normal approximation for the rank correlation p-value")); invokeRestart("muffleWarning") } })
  pairs <- utils::combn(vars, 2, simplify = FALSE)
  long <- function(Z, scope) {
    out <- do.call(rbind, lapply(pairs, function(pr) {
      ok <- stats::complete.cases(Z[pr]); x <- Z[[pr[1]]][ok]; y <- Z[[pr[2]]][ok]
      ct <- one(x, y)
      data.frame(var1 = pr[1], var2 = pr[2], scope = scope, method = method, r = unname(ct$estimate), statistic = unname(ct$statistic),
                 df = if (!is.null(ct$parameter)) unname(ct$parameter) else NA_real_, p = ct$p.value,
                 ci_low = if (!is.null(ct$conf.int)) ct$conf.int[1] else NA_real_, ci_high = if (!is.null(ct$conf.int)) ct$conf.int[2] else NA_real_,
                 n = sum(ok), stringsAsFactors = FALSE)
    }))
    out$p_adj <- stats::p.adjust(out$p, adj); out$adjust <- adj
    out$strength <- as.character(cut(abs(out$r), unlist(rules$strength_breaks %||% c(0, 0.1, 0.3, 0.5, 1.01)),
                                     unlist(rules$strength_labels %||% c("negligible", "weak", "moderate", "strong")), right = FALSE))
    out
  }
  tab <- long(X, "raw")
  # within the design: residuals of each trait on the design's fixed terms (and random block intercepts as factors)
  p <- design$parsed
  role <- function(v) p$role[[v]] %||% "treatment"
  rhs <- c(p$fixed_labels, unlist(lapply(p$random, function(r) if (isTRUE(r$intercept_only) && length(r$vars) == 1L && role(r$vars) == "block") r$vars)))
  rhs <- rhs[vapply(rhs, function(l) all(all.vars(str2lang(l)) %in% p$vars), logical(1))]
  res_note <- NULL
  if (length(rhs) && !isFALSE(cs$within_design %||% rules$within_design)) {
    fr <- tryCatch(.ba_frame(dat, p), error = function(e) NULL)
    if (!is.null(fr)) {
      R <- as.data.frame(lapply(vars, function(v) { fr$v_ <- X[[v]]
        m <- stats::lm(stats::reformulate(rhs, response = "v_"), data = fr, na.action = stats::na.exclude); stats::residuals(m) }))
      names(R) <- vars
      tab <- rbind(tab, long(R, "within design"))
      res_note <- sprintf("within design: residuals of each trait on %s (p-values treat residuals as independent; their df is slightly smaller)", paste(rhs, collapse = " + "))
    }
  }
  # checks for Pearson (reported, not used to switch): normality per trait
  if (method == "pearson") {
    sw <- vapply(vars, function(v) { x <- stats::na.omit(X[[v]]); if (length(x) >= 3 && length(x) <= 5000 && stats::sd(x) > 0) stats::shapiro.test(x)$p.value else NA_real_ }, 0)
    if (any(sw < 0.05, na.rm = TRUE)) notes <- c(notes, sprintf("Shapiro-Wilk p < 0.05 for %s: Spearman (analysis.correlation.method) does not assume normality", paste(names(sw)[which(sw < 0.05)], collapse = ", ")))
  }
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  mat <- function(sc) { M <- diag(length(vars)); dimnames(M) <- list(vars, vars)
    s <- tab[tab$scope == sc, ]; for (i in seq_len(nrow(s))) M[s$var1[i], s$var2[i]] <- M[s$var2[i], s$var1[i]] <- s$r[i]; data.frame(variable = vars, M, check.names = FALSE) }
  paths <- c(correlation = .ba_write_tsv(tab, file.path(outdir, "correlation.tsv")),
             correlation_matrix = .ba_write_tsv(mat("raw"), file.path(outdir, "correlation_matrix.tsv")))
  if (any(tab$scope == "within design")) paths <- c(paths, correlation_matrix_within = .ba_write_tsv(mat("within design"), file.path(outdir, "correlation_matrix_within.tsv")))
  raw <- tab[tab$scope == "raw", ]
  summary <- list(data_type = "traits", method = method, package = "stats", version = as.character(utils::packageVersion("stats")),
                  call = sprintf("stats::cor.test(x, y, method = '%s') per pair; p.adjust('%s')", method, adj), variables = as.list(vars),
                  n_pairs = nrow(raw), n_significant = sum(raw$p_adj < alpha), within_design = res_note,
                  strongest = if (nrow(raw)) { i <- which.max(abs(raw$r)); list(var1 = raw$var1[i], var2 = raw$var2[i], r = raw$r[i], p_adj = raw$p_adj[i]) },
                  reference = rules$reference, tutorial = "docs/engines/correlation_agreement.md", warnings = notes)
  list(summary = summary, paths = paths, results = tab)
}
