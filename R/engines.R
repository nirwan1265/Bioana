# Analysis engines: which package fits each data set.
#
# inst/rules/engines.tsv lists every engine (package, function called, model,
# test statistic, multiple testing, assumptions, reference). inst/rules/
# engine_selection.tsv ranks the engines per data type and per design: "fixed"
# when analysis.model has no random terms (independent units), "random" when
# it has (repeated or nested units), "any" for both. The first installed
# engine that fits the inputs is the primary; native Bioana code is only the
# last fallback. The next gold-standard engine is offered as the secondary.
# Every result carries a provenance block: package, version, exact call,
# model, input, assumptions, test statistic, multiple testing, experimental
# unit, and whether a fallback was used.

#' The analysis engines Bioana can call
#'
#' @return `inst/rules/engines.tsv` with an `installed` column: one row per
#'   engine (package, function, model, test statistic, multiple-testing method,
#'   assumptions, reference, tutorial).
#' @export
bioana_engines <- function() {
  e <- .ba_tsv("engines")
  e$installed <- vapply(seq_len(nrow(e)), function(i) .ba_engine_installed(e[i, ]), logical(1))
  e
}

#' Which engine Bioana would run for a data type and design
#'
#' @param data_type A row of `analyses.tsv`, or `"outcome:<type>"` for the primary outcome.
#' @param design `"fixed"` (no random terms in analysis.model) or `"random"`.
#' @return `inst/rules/engine_selection.tsv` rows for that data type and design, in rank
#'   order, with `installed` and `role` (`primary`, `secondary`, `fallback`, `not installed`).
#' @export
bioana_engine_plan <- function(data_type, design = c("fixed", "random")) {
  design <- match.arg(design)
  s <- .ba_selection_rows(data_type, design)
  e <- .ba_tsv("engines")
  s$package <- e$package[match(s$engine, e$engine)]
  s$installed <- vapply(s$engine, function(x) .ba_engine_installed(e[e$engine == x, ]), logical(1))
  s$role <- "not installed"
  ok <- which(s$installed)
  if (length(ok)) s$role[ok] <- "fallback"
  gold <- ok[s$engine[ok] != "bioana_native" & s$package[ok] != "bioana"]
  if (length(gold)) s$role[gold[1]] <- "primary"
  if (length(gold) > 1) s$role[gold[2]] <- "secondary"
  if (!length(gold) && length(ok)) s$role[ok[1]] <- "primary (fallback)"
  s
}

.ba_engine_installed <- function(row) {
  req <- strsplit(row$requires, ",", fixed = TRUE)[[1]]
  req <- trimws(req[nzchar(req)])
  all(vapply(req, requireNamespace, logical(1), quietly = TRUE))
}

.ba_selection_rows <- function(key, kind) {
  s <- .ba_tsv("engine_selection")
  s <- s[s$data_type == key & s$design %in% c(kind, "any"), , drop = FALSE]
  s[order(as.integer(s$rank)), , drop = FALSE]
}

.ba_design_kind <- function(p) if (length(p$random)) "random" else "fixed"

# Ordered engines to try. `inputs` = the input kinds the data type can give (counts, values, ...).
.ba_engine_chain <- function(key, design, inputs = NULL, override = NULL, allow_fallback = FALSE) {
  e <- .ba_tsv("engines")
  p <- design$parsed
  kind <- .ba_design_kind(p)
  s <- .ba_selection_rows(key, kind)
  if (!nrow(s)) stop(sprintf("No engine rows for '%s' in engine_selection.tsv.", key), call. = FALSE)
  chain <- s$engine
  fit <- stats::setNames(s$design_fit, s$engine)
  if (!is.null(override)) {
    if (!override %in% e$engine) stop(sprintf("Unknown engine '%s'; see bioana_engines().", override), call. = FALSE)
    chain <- unique(c(override, chain))
    if (is.na(fit[override])) fit[override] <- .ba_design_fit_of(e[e$engine == override, ], p)
  }
  skipped <- list()
  keep <- character()
  for (x in chain) {
    row <- e[e$engine == x, ]
    if (!is.null(inputs) && !row$input %in% inputs) { skipped[[x]] <- sprintf("needs %s input", row$input); next }
    if (identical(fit[[x]], "drops_random") && !isTRUE(allow_fallback) && !identical(x, override)) {
      skipped[[x]] <- sprintf("would drop the random terms (%s): only with allow_fallback = TRUE", paste(.ba_random_labels(p), collapse = ", "))
      next
    }
    if (!identical(fit[[x]], "drops_random") && !.ba_random_capacity(row$max_random, p)) {
      skipped[[x]] <- sprintf("models %s random terms; the design has %s", c("0" = "no", "1" = "one random intercept,", any = "any")[[row$max_random]],
                              paste(.ba_random_labels(p), collapse = ", "))
      next
    }
    if (row$package == "bioana" && !isTRUE(allow_fallback) && !identical(x, override)) {
      skipped[[x]] <- "Bioana's own code: off unless allow_fallback = TRUE"; next
    }
    if (!.ba_engine_installed(row)) { skipped[[x]] <- sprintf("not installed (%s)", row$requires); next }
    keep <- c(keep, x)
  }
  if (!length(keep)) {
    ok <- s$engine[s$design_fit != "drops_random" & e$package[match(s$engine, e$engine)] != "bioana"]
    pk <- unique(e$package[match(ok, e$engine)])
    stop(sprintf("No citable engine can run %s for this design (%s; random terms: %s). %s %s",
                 key, kind, if (length(p$random)) paste(.ba_random_labels(p), collapse = ", ") else "none",
                 if (length(pk)) paste0("Install one of: ", paste(pk, collapse = ", "), ".") else "",
                 paste("Engines not used:", paste(names(skipped), unlist(skipped), sep = ": ", collapse = "; "))),
         call. = FALSE)
  }
  list(chain = keep, preferred = chain[1], skipped = skipped, fit = fit[keep])
}

.ba_random_labels <- function(p) vapply(p$random, `[[`, "", "label")

# Can an engine with this max_random hold the design's random terms? "1" = one random intercept.
.ba_random_capacity <- function(max_random, p) {
  n <- length(p$random)
  if (n == 0L || identical(max_random, "any")) return(TRUE)
  if (identical(max_random, "0")) return(FALSE)
  n == 1L && !identical(p$random[[1]]$intercept_only, FALSE)
}

# design_fit of an engine named by the user for a design it has no selection row for.
.ba_design_fit_of <- function(row, p) {
  if (!length(p$random) || identical(row$max_random, "any")) return("full")
  if (.ba_random_capacity(row$max_random, p)) return("approx")
  "drops_random"
}

# Run the chain: the first engine that works. After a failure, only an engine that fits the design
# fully may take over; an approximation or one that drops random terms stops the analysis instead.
# Warnings raised during the fit are collected, not lost.
.ba_run_chain <- function(ch, run, what) {
  failed <- list()
  for (eng in ch$chain) {
    if (length(failed) && !identical(ch$fit[[eng]], "full")) {
      stop(sprintf("%s failed (%s). The next engine, %s, %s, so it was not run in its place. Fix the error, or choose an engine yourself: params = list(%s = list(engine = '%s')).",
                   names(failed)[1], failed[[1]], eng,
                   if (identical(ch$fit[[eng]], "approx")) "only approximates this design" else "drops the design's random terms",
                   what, eng), call. = FALSE)
    }
    w <- character()
    r <- tryCatch(withCallingHandlers(run(eng), warning = function(x) { w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning") }),
                  error = function(err) err)
    if (!inherits(r, "error")) return(list(result = r, engine = eng, failed = failed, fit_warnings = unique(w)))
    failed[[eng]] <- conditionMessage(r)
  }
  stop("Every engine failed: ", paste(names(failed), unlist(failed), sep = ": ", collapse = "; "), call. = FALSE)
}

# The next gold-standard engine after the one used (never native code).
.ba_secondary_engine <- function(chain, used) {
  e <- .ba_tsv("engines")
  rest <- chain[seq_along(chain) > match(used, chain)]
  rest <- rest[e$package[match(rest, e$engine)] != "bioana"]
  if (length(rest)) rest[1] else NULL
}

# Provenance block for one fitted data set.
.ba_provenance <- function(engine, design, call, input, n_units, chain = NULL, random = NULL, formula = NULL) {
  e <- .ba_tsv("engines")
  row <- e[e$engine == engine, ]
  pkg <- row$package
  fit_kind <- if (!is.null(chain$fit) && engine %in% names(chain$fit)) chain$fit[[engine]] else "full"
  list(engine = engine, package = pkg,
       version = if (pkg == "bioana") as.character(utils::packageVersion("bioana")) else
         tryCatch(as.character(utils::packageVersion(pkg)), error = function(err) NA_character_),
       language = row$language, call = call, formula = formula %||% design$model,
       model = row$model, input = input, assumptions = row$assumptions,
       test_statistic = row$test_statistic, multiple_testing = row$multiple_testing,
       estimate_scale = row$estimate_scale, random_terms = random %||% row$random_terms,
       design_fit = fit_kind, design_fit_note = if (!identical(fit_kind, "full")) .ba_fit_note(engine, design) else NULL,
       experimental_unit = design$spec$experimental_unit %||% "unit (one row of sample_metadata per experimental unit)",
       n_units = n_units,
       fallback_used = !is.null(chain) && (row$package == "bioana" || !identical(engine, chain$preferred)),
       not_used = chain$skipped %||% list(), reference = row$reference, tutorial = row$doc)
}

# Design pieces shared by the feature engines.
.ba_feature_design <- function(sample_ids, meta, p) {
  ids <- intersect(meta$sample_id, sample_ids)
  if (length(ids) < 3) stop(sprintf("Only %d samples of the data match sample_metadata sample_id.", length(ids)), call. = FALSE)
  frame <- .ba_frame(meta[match(ids, meta$sample_id), , drop = FALSE], p)
  fvars <- intersect(unlist(lapply(p$fixed_labels, function(l) all.vars(str2lang(l)))), names(frame))
  ok <- stats::complete.cases(frame[, fvars, drop = FALSE])
  frame <- frame[ok, , drop = FALSE]; ids <- ids[ok]
  rhs <- stats::reformulate(p$fixed_labels)
  X <- stats::model.matrix(rhs, data = frame)
  qx <- qr(X)
  X <- X[, sort(qx$pivot[seq_len(qx$rank)]), drop = FALSE]
  coefs <- colnames(X)[vapply(colnames(X), function(t) .ba_coef_role(t, p), "") == "treatment"]
  block <- NULL
  if (length(p$random) == 1L) block <- interaction(frame[p$random[[1]]$vars], drop = TRUE)
  rownames(frame) <- ids
  list(ids = ids, frame = frame, X = X, coefs = coefs, block = block, n_random = length(p$random),
       random_labels = vapply(p$random, `[[`, "", "label"), rhs = rhs,
       full_rhs = stats::as.formula(paste("~", paste(deparse(p$formula[[3]]), collapse = ""))))
}

.ba_need_coefs <- function(d) {
  if (!length(d$coefs)) stop("the model has no estimable fixed treatment coefficient", call. = FALSE)
  if (nrow(d$X) - ncol(d$X) < 1) stop("no residual degrees of freedom", call. = FALSE)
}

# Standard summary for a feature engine result.
.ba_engine_summary <- function(tab, d, method, alpha) {
  list(method = method, design = paste(deparse(d$rhs), collapse = ""), n_features_tested = length(unique(tab$feature_id)),
       n_samples = nrow(d$X), contrasts = unique(tab$contrast),
       significant = as.list(tapply(!is.na(tab$padj) & tab$padj < alpha, tab$contrast, sum)))
}

# Input kinds a method's inputs can feed (engines.tsv `input`).
.ba_input_kinds <- function(inputs) {
  c(if (!is.null(inputs$counts)) "counts", if (!is.null(inputs$values)) "values",
    if (!is.null(inputs$methylated) && !is.null(inputs$coverage)) "bisulfite",
    if (!is.null(inputs$ip) && !is.null(inputs$input)) "ip_input")
}

# Run the chain for a feature data type: first engine that works.
.ba_feature_engines <- function(inputs, entry, design, alpha, args) {
  ch <- .ba_engine_chain(entry$data_type, design, .ba_input_kinds(inputs), args$engine, args$allow_fallback)
  checks <- .ba_check_design(design, .ba_input_samples(inputs), what = "features")
  rc <- .ba_run_chain(ch, function(eng) get(paste0(".ba_engine_", eng), envir = asNamespace("bioana"))(inputs, design, alpha, args),
                      entry$data_type)
  r <- rc$result; eng <- rc$engine; failed <- rc$failed
  ch$skipped <- c(ch$skipped, lapply(failed, function(m) paste("failed:", m)))
  warn <- c(r$warnings, if (length(failed)) sprintf("engine %s failed (%s); used %s, which fits the design fully", names(failed), unlist(failed), eng))
  r$summary$engine <- .ba_provenance(eng, design, r$call, r$input, r$n_units, ch, random = r$random)
  r$summary$diagnostics <- .ba_diagnostics(checks, rc$fit_warnings, r$diagnostics)
  r$summary$secondary_engine <- .ba_secondary_engine(ch$chain, eng)
  list(table = r$table, warnings = warn, summary = r$summary)
}

# Agreement between the primary and secondary results of one data set.
.ba_concordance <- function(a, b, alpha) {
  key <- intersect(c("feature_id", "marker_id", "term", "contrast"), intersect(names(a), names(b)))
  m <- merge(a, b, by = key, suffixes = c("_1", "_2"))
  if (!nrow(m)) return(list(n_shared = 0L))
  s1 <- !is.na(m$padj_1) & m$padj_1 < alpha; s2 <- !is.na(m$padj_2) & m$padj_2 < alpha
  ok <- is.finite(m$estimate_1) & is.finite(m$estimate_2)
  list(n_shared = nrow(m),
       estimate_cor = if (sum(ok) > 2) round(stats::cor(m$estimate_1[ok], m$estimate_2[ok]), 3) else NA,
       sign_agreement = if (any(ok)) round(mean(sign(m$estimate_1[ok]) == sign(m$estimate_2[ok])), 3) else NA,
       log10p_cor = if (sum(is.finite(m$p_1) & is.finite(m$p_2)) > 2) round(stats::cor(-log10(m$p_1), -log10(m$p_2), use = "complete.obs"), 3) else NA,
       significant_both = sum(s1 & s2), significant_primary_only = sum(s1 & !s2), significant_secondary_only = sum(!s1 & s2),
       jaccard = if (sum(s1 | s2)) round(sum(s1 & s2) / sum(s1 | s2), 3) else NA)
}

.ba_fit_note <- function(engine, design) {
  s <- .ba_tsv("engine_selection")
  n <- s$note[s$engine == engine & s$design == .ba_design_kind(design$parsed)]
  n <- n[nzchar(n)]
  if (length(n)) n[1] else NULL
}

# Sample ids an input list covers (the columns of its first matrix).
.ba_input_samples <- function(inputs) {
  m <- inputs$counts %||% inputs$values %||% inputs$coverage %||% inputs$ip
  colnames(m)
}

# Diagnostics block: the design checks, warnings raised while fitting, engine-specific notes.
.ba_diagnostics <- function(checks, fit_warnings = character(), engine_notes = NULL) {
  ck <- if (is.null(checks) || !nrow(checks)) list() else lapply(split(checks, seq_len(nrow(checks))), as.list)
  list(checks = unname(ck), n_check_warnings = sum(checks$status == "warn"),
       fit_warnings = as.list(fit_warnings), engine = engine_notes)
}
