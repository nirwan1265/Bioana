#' Analyse a Biodes experiment: every data set in the folder, with the design's model
#'
#' Reads the Biodes design folder (`study_spec.json`, `field_layout.csv`,
#' `sample_metadata_template.tsv`) and its `data_manifest.json`, and for each
#' data set runs the analysis `inst/rules/analyses.tsv` names for its data
#' type, using the fixed terms of the spec's `analysis.model` (Biodes already
#' decided the experimental unit, blocks and random terms). Works the same on
#' measured data and on data simulated by Simitall; when a data set carries
#' its truth (simulated), the results are scored against it.
#'
#' Results go to `<outdir>/<id>/results.tsv` (and `recovery.tsv` for the
#' outcome) plus `<outdir>/bioana_summary.json`.
#'
#' @param dir Biodes design folder with data.
#' @param outdir Where to write results (default `<dir>/bioana`).
#' @param alpha Significance level (BH-adjusted p for features).
#' @param params Optional per data type settings, e.g.
#'   `list(bisulfite_counts = list(window_bp = 500), outcome = list(subsamples = "nest"))`;
#'   see `inst/rules/method_parameters.tsv`.
#' @param only Optional ids or data types to run.
#' @param quiet `TRUE` skips the one-line-per-data-set report.
#' @param secondary Run the secondary engine (the next gold-standard package in
#'   `inst/rules/engine_selection.tsv`) as well, and report how well it agrees
#'   with the primary: `TRUE`, `FALSE`, or `NA` (default) to ask when the
#'   session is interactive and `quiet = FALSE` (no otherwise). Its results go
#'   to `<outdir>/<id>/results_<engine>.tsv`. Pick an engine per data type with
#'   `params`, e.g. `list(feature_counts = list(engine = "edger"))`; see
#'   [bioana_engines()] and [bioana_engine_plan()].
#' @param allow_fallback `FALSE` (default): every analysis is run by a
#'   published package (citable); a data set whose packages are not installed
#'   stops with the packages to install. `TRUE` lets Bioana's own code run
#'   instead (not citable; marked `FALLBACK` in the report and
#'   `fallback_used` in the summary). Naming a Bioana engine in `params`
#'   (`engine = "bioana_native"`) also allows it for that data set.
#' @return Invisibly, a list: `summary` (per data set, as in `bioana_summary.json`),
#'   `results` (tables by id), `recovery` (truth tables by id, simulated data), `fits`, `paths`.
#' @examples
#' \dontrun{
#' simitall::simitall_from_biodes("~/biodes_runs/my_design")   # or fill it with your data
#' bioana_analyze("~/biodes_runs/my_design")
#' }
#' @export
bioana_analyze <- function(dir, outdir = file.path(dir, "bioana"), alpha = 0.05, params = NULL, only = NULL, quiet = FALSE,
                           secondary = NA, allow_fallback = FALSE) {
  design <- .ba_read_biodes(dir)
  man <- bioana_manifest(dir)
  if (!length(man$data)) {
    stop(sprintf("Nothing to analyse in %s: no data_manifest.json, an empty primary_outcome column and no counts.tsv. Fill the folder with data (simitall::simitall_from_biodes() or your measurements) and list it with bioana_write_manifest().", dir), call. = FALSE)
  }
  tab <- .ba_tsv("analyses")
  defaults <- .ba_tsv("method_parameters")
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  summ <- list(); results <- list(); recovery <- list(); fits <- list(); paths <- character(); skipped <- list(); ran <- list()
  for (e in man$data) {
    if (!is.null(only) && !e$id %in% only && !e$data_type %in% only) next
    row <- tab[tab$data_type == e$data_type, , drop = FALSE]
    if (!nrow(row) || row$status != "ready") {
      skipped[[e$id]] <- list(data_type = e$data_type, status = if (nrow(row)) row$status else "unknown data_type",
                              note = if (nrow(row)) row$note else "not in inst/rules/analyses.tsv")
      next
    }
    d <- defaults[defaults$data_type == e$data_type, , drop = FALSE]
    args <- stats::setNames(lapply(seq_len(nrow(d)), function(i) .ba_parse_value(d$value[i], d$type[i])), d$argument)
    args <- utils::modifyList(args, params[[e$data_type]] %||% list())
    args <- utils::modifyList(args, params[[e$id]] %||% list())
    args$allow_fallback <- isTRUE(args$allow_fallback %||% allow_fallback)
    fn <- get(row$method, envir = asNamespace("bioana"))
    r <- tryCatch(fn(e, design, file.path(outdir, e$id), alpha, args), error = function(err) err)
    if (inherits(r, "error")) {
      skipped[[e$id]] <- list(data_type = e$data_type, status = "error", note = conditionMessage(r))
      next
    }
    summ[[e$id]] <- r$summary
    if (!is.null(e$source)) summ[[e$id]]$source <- e$source
    ran[[e$id]] <- list(entry = e, fn = fn, args = args)
    results[[e$id]] <- r$results
    if (!is.null(r$recovery)) recovery[[e$id]] <- r$recovery
    if (!is.null(r$fit)) fits[[e$id]] <- r$fit
    paths <- c(paths, stats::setNames(r$paths, paste0(e$id, "/", names(r$paths))))
  }
  cands <- Filter(Negate(is.null), lapply(summ, `[[`, "secondary_engine"))
  if (length(cands) && is.na(secondary)) {
    secondary <- FALSE
    if (interactive() && !quiet) {
      e_tab <- .ba_tsv("engines")
      cat(sprintf("Secondary engines available: %s\n", paste(sprintf("%s -> %s", names(cands), e_tab$package[match(unlist(cands), e_tab$engine)]), collapse = ", ")))
      secondary <- tolower(substr(readline("Run them too and compare with the primary? [y/N] "), 1, 1)) == "y"
    }
  }
  if (isTRUE(secondary)) {
    for (id in names(cands)) {
      x <- ran[[id]]; sec <- cands[[id]]
      tmp <- file.path(tempfile("bioana_secondary_"), id)
      r2 <- tryCatch(x$fn(x$entry, design, tmp, alpha, utils::modifyList(x$args, list(engine = sec))), error = function(err) err)
      if (inherits(r2, "error")) {
        summ[[id]]$secondary <- list(engine = sec, status = "error", note = conditionMessage(r2))
        next
      }
      p2 <- .ba_write_tsv(r2$results, file.path(outdir, id, paste0("results_", sec, ".tsv")))
      paths <- c(paths, stats::setNames(p2, paste0(id, "/results_", sec)))
      summ[[id]]$secondary <- list(engine = r2$summary$engine, concordance = .ba_concordance(results[[id]], r2$results, alpha),
                                   recovery = r2$summary$recovery)
    }
  }
  out <- list(bioana_version = as.character(utils::packageVersion("bioana")), design = design$spec$design$chosen,
              model = design$model, handoff_version = design$spec$handoff_version,
              manifest = man$written_by %||% "", alpha = alpha, analyses = summ, not_analysed = skipped)
  paths <- c(paths, summary = .ba_write_json(out, file.path(outdir, "bioana_summary.json")))
  if (!quiet) .ba_report(out)
  invisible(list(summary = out, results = results, recovery = recovery, fits = fits, paths = paths))
}

.ba_parse_value <- function(value, type) {
  if (!nzchar(value)) return(NULL)
  switch(type, integer = as.integer(value), numeric = as.numeric(value), logical = as.logical(value), value)
}

# One line per data set.
.ba_report <- function(s) {
  cat(sprintf("bioana: design '%s' (%s)\n", s$design %||% "?", s$model))
  for (id in names(s$analyses)) {
    a <- s$analyses[[id]]
    line <- if (a$data_type == "outcome") {
      sprintf("%s outcome, %s", a$outcome_type, a$method)
    } else if (a$data_type == "genotypes") {
      sprintf("%d lines x %d markers tested, %d significant", a$n_lines, a$n_features_tested, a$significant$markers)
    } else {
      sig <- unlist(a$significant)
      sprintf("%d features tested, %s significant", a$n_features_tested %||% 0L,
              if (length(sig)) paste(paste0(names(sig), " ", sig), collapse = ", ") else "0")
    }
    rec <- a$recovery
    rtxt <- if (is.null(rec)) "" else if (a$data_type == "genotypes") {
      sprintf("; truth: %d/%d QTN found (FDR %s), sign %s, line means r %s", rec$true_positives, rec$testable,
              rec$fdr, rec$sign_agreement, a$line_mean_cor %||% "NA")
    } else if (a$data_type == "outcome") {
      if (!is.null(rec$coverage)) sprintf("; truth: %d/%d coefficients in their CI", rec$covered, rec$n_terms) else
        sprintf("; truth: sign %s, detected %s", rec$sign_agreement, rec$detected)
    } else {
      paste0("; truth: ", paste(vapply(rec, function(r) sprintf("%s found %d/%d (FDR %s)", r$contrast, r$true_positives, r$true_changes, r$fdr), ""), collapse = ", "))
    }
    eng <- a$engine
    etxt <- if (is.null(eng)) "" else sprintf(" [%s %s%s%s]", eng$package, eng$version, if (isTRUE(eng$fallback_used)) ", FALLBACK" else "",
                                              switch(eng$design_fit %||% "full", approx = ", approximates the design", drops_random = ", RANDOM TERMS DROPPED", ""))
    g <- a$goal
    gtxt <- if (is.null(g) || (identical(g$id, "treatment_effect") && identical(g$source, "default")) || identical(g$id, "association")) "" else
      sprintf(" {goal %s, %s%s}", g$id, if (identical(g$test, "lrt")) "joint test" else "Wald", if (isTRUE(g$assumed_by_biodes)) ", assumed by Biodes" else "")
    cat(sprintf("  [%s]%s%s %s%s\n", id, etxt, gtxt, line, rtxt))
    if (!is.null(g$result$p)) cat(sprintf("      %s: statistic %.3g on %g df, p = %.3g\n", g$label, g$result$statistic, g$result$df, g$result$p))
    if (!is.null(g$note)) cat(sprintf("      goal: %s\n", g$note))
    if (!is.null(g$scale_note)) cat(sprintf("      trend: %s\n", g$scale_note))
    if (!is.null(a$recovery_note)) cat(sprintf("      truth: %s\n", a$recovery_note))
    st <- a$stat_test
    if (!is.null(st)) {
      if (identical(st$status, "run")) {
        cat(sprintf("      ANOVA [%s %s] %s%s: %s\n", st$package, st$version, c(within = "repeated measures", mixed = "mixed (split-plot)", between = "factorial")[[st$layout]],
                    if (!identical(st$correction, "none")) sprintf(", %s-corrected", st$correction) else "",
                    paste(vapply(st$table, function(r) sprintf("%s F(%.3g, %.3g) = %.3g, p = %.3g, ges %.2f", r$term, r$num_df, r$den_df, r$F, r$p, r$ges), ""), collapse = "; ")))
        if (identical(st$posthoc$status, "run")) cat(sprintf("      post-hoc [emmeans %s] %s: %d comparisons (%s), %d with adjusted p < 0.05\n",
          st$posthoc$version, paste(unlist(st$posthoc$families), collapse = ", "), st$posthoc$n_comparisons, st$posthoc$adjust, st$posthoc$n_significant))
        w <- c(st$checks$warnings, st$note)
        if (length(w)) cat(sprintf("      ANOVA checks: %s\n", paste(w, collapse = "; ")))
      } else cat(sprintf("      ANOVA not run: %s\n", st$note))
    }
    gr <- a$goal_recovery
    if (!is.null(gr)) cat(sprintf("      goal truth (%s, terms %s%s): %s\n", gr$goal, paste(unlist(gr$model_terms), collapse = ", "),
      if (isFALSE(gr$terms_match)) ", NOT the terms Bioana tested" else "",
      if (!isTRUE(gr$scored)) gr$note %||% "not scored" else if (!is.null(gr$tpr)) paste0(sprintf("%d/%d goal features found (tpr %.2f)", gr$found, gr$n_tested, gr$tpr),
        if (!is.null(gr$called_not_planted)) sprintf("; of %d called, %d with another planted effect,%s %d with none (%.0f%%)", gr$n_called,
          gr$called_other_planted, if (isTRUE(gr$called_hidden_factor > 0)) sprintf(" %d moved by hidden factors,", gr$called_hidden_factor) else "",
          gr$called_not_planted, 100 * gr$false_call_rate) else sprintf(", %d called outside the goal set", gr$called_outside_goal)) else
        paste0(if (isTRUE(gr$detected)) "detected" else "not detected", sprintf(" (p = %.3g)", gr$p),
               if (!is.null(gr$covered)) sprintf(", estimate %.3g vs true %.3g, %s", gr$estimate, gr$true_effect, if (isTRUE(gr$covered)) "in CI" else "outside CI") else "",
               if (!is.null(gr$note)) paste0("; ", gr$note) else "")))
    dg <- a$diagnostics
    if (!is.null(dg)) {
      cw <- Filter(function(x) identical(x$status, "warn"), dg$checks)
      if (length(cw) || length(dg$fit_warnings)) {
        cat(sprintf("      check: %s\n", vapply(cw, `[[`, "", "detail")), sep = "")
        if (length(dg$fit_warnings)) cat(sprintf("      fit: %d warning(s), e.g. %s\n", length(dg$fit_warnings), substr(dg$fit_warnings[[1]], 1, 120)))
      }
    }
    sec <- a[["secondary"]]
    if (!is.null(sec$concordance)) {
      cc <- sec$concordance
      cat(sprintf("      secondary %s %s: estimates r %s, significant in both %d, primary only %d, secondary only %d\n",
                  sec$engine$package, sec$engine$version, cc$estimate_cor, cc$significant_both, cc$significant_primary_only, cc$significant_secondary_only))
    } else if (!is.null(sec$status)) {
      cat(sprintf("      secondary %s failed: %s\n", sec$engine, sec$note))
    }
  }
  offer <- Filter(Negate(is.null), lapply(s$analyses, function(a) if (is.null(a[["secondary"]])) a[["secondary_engine"]]))
  if (length(offer)) cat(sprintf("  secondary engines not run (%s): bioana_analyze(dir, secondary = TRUE)\n",
                                 paste(sprintf("%s -> %s", names(offer), unlist(offer)), collapse = ", ")))
  for (id in names(s$not_analysed)) cat(sprintf("  [%s] not analysed (%s): %s\n", id, s$not_analysed[[id]]$status, s$not_analysed[[id]]$note))
}
