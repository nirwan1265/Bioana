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
.ba_engine_chain <- function(key, design, inputs = NULL, override = NULL) {
  e <- .ba_tsv("engines")
  s <- .ba_selection_rows(key, .ba_design_kind(design$parsed))
  if (!nrow(s)) stop(sprintf("No engine rows for '%s' in engine_selection.tsv.", key), call. = FALSE)
  chain <- s$engine
  if (!is.null(override)) {
    if (!override %in% e$engine) stop(sprintf("Unknown engine '%s'; see bioana_engines().", override), call. = FALSE)
    chain <- unique(c(override, chain))
  }
  skipped <- list()
  keep <- character()
  for (x in chain) {
    row <- e[e$engine == x, ]
    if (!is.null(inputs) && !row$input %in% inputs) { skipped[[x]] <- sprintf("needs %s input", row$input); next }
    if (!.ba_engine_installed(row)) { skipped[[x]] <- sprintf("not installed (%s)", row$requires); next }
    keep <- c(keep, x)
  }
  list(chain = keep, preferred = chain[1], skipped = skipped)
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
  list(engine = engine, package = pkg,
       version = if (pkg == "bioana") as.character(utils::packageVersion("bioana")) else
         tryCatch(as.character(utils::packageVersion(pkg)), error = function(err) NA_character_),
       language = row$language, call = call, formula = formula %||% design$model,
       model = row$model, input = input, assumptions = row$assumptions,
       test_statistic = row$test_statistic, multiple_testing = row$multiple_testing,
       estimate_scale = row$estimate_scale, random_terms = random %||% row$random_terms,
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

# Run the chain for a feature data type: first engine that works.
.ba_feature_engines <- function(inputs, entry, design, alpha, args) {
  ch <- .ba_engine_chain(entry$data_type, design, names(inputs), args$engine)
  if (!length(ch$chain)) stop("No engine available: ", paste(names(ch$skipped), unlist(ch$skipped), sep = " ", collapse = "; "), call. = FALSE)
  failed <- list()
  for (eng in ch$chain) {
    r <- tryCatch(get(paste0(".ba_engine_", eng), envir = asNamespace("bioana"))(inputs, design, alpha, args),
                  error = function(err) err)
    if (!inherits(r, "error")) break
    failed[[eng]] <- conditionMessage(r)
  }
  if (inherits(r, "error")) stop("Every engine failed: ", paste(names(failed), unlist(failed), sep = ": ", collapse = "; "), call. = FALSE)
  ch$skipped <- c(ch$skipped, lapply(failed, function(m) paste("failed:", m)))
  warn <- c(r$warnings, if (length(failed)) sprintf("engine %s failed (%s); used %s", names(failed), unlist(failed), eng))
  r$summary$engine <- .ba_provenance(eng, design, r$call, r$input, r$n_units, ch, random = r$random)
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
