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
  # analyses the spec asks for on the sample metadata: correlation between traits, agreement between raters
  have <- vapply(man$data, function(e) e$data_type %||% "", "")
  if (!is.null(design$spec$analysis$correlation) && !"traits" %in% have) man$data <- c(man$data, list(list(id = "correlation", data_type = "traits")))
  if (!is.null(design$spec$analysis$agreement) && !"ratings" %in% have) man$data <- c(man$data, list(list(id = "agreement", data_type = "ratings")))
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
    } else if (a$data_type == "traits") {
      sprintf("%s correlation of %d traits: %d/%d pairs with adjusted p < 0.05%s", a$method, length(a$variables), a$n_significant, a$n_pairs,
              if (!is.null(a$strongest)) sprintf("; strongest %s-%s r = %.2f", a$strongest$var1, a$strongest$var2, a$strongest$r) else "")
    } else if (a$data_type == "ratings") {
      pm <- a$primary
      sprintf("agreement of %d raters on %d units (%s): %s %.3f%s%s", length(a$raters), a$n_units, a$scale, pm$measure, pm$estimate,
              if (!is.na(pm$ci_low)) sprintf(" [%.3f, %.3f]", pm$ci_low, pm$ci_high) else "", if (!is.na(pm$interpretation)) paste0(", ", pm$interpretation) else "")
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
    if (!is.null(st)) for (x in c(list(st), st$also)) .ba_print_stat_test(x)
    if (a$data_type %in% c("traits", "ratings") && length(a$warnings)) cat(sprintf("      checks: %s\n", paste(a$warnings, collapse = "; ")))
    if (!is.null(a$within_design)) cat(sprintf("      %s\n", a$within_design))

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

.ba_print_stat_test <- function(st) {
  if (identical(st$status, "run") && identical(st$test, "survival")) return(.ba_print_survival(st))
  if (identical(st$status, "run") && identical(st$test, "competing_risks")) return(.ba_print_competing(st))
  if (!identical(st$status, "run")) { cat(sprintf("      test not run%s: %s\n", if (!is.null(st$test)) paste0(" (", st$test, ")") else "", st$note %||% "")); return(invisible()) }
  if (startsWith(st$test, "anova_")) {
    cat(sprintf("      ANOVA [%s %s] %s%s: %s\n", st$package, st$version, c(within = "repeated measures", mixed = "mixed (split-plot)", between = "factorial")[[st$layout]],
                if (!identical(st$correction, "none")) sprintf(", %s-corrected", st$correction) else "",
                paste(vapply(st$table, function(r) sprintf("%s F(%.3g, %.3g) = %.3g, p = %.3g, ges %.2f", r$term, r$num_df, r$den_df, r$F, r$p, r$ges), ""), collapse = "; ")))
  } else {
    r <- st$table[[1]]
    df <- if (!is.na(r$df1) && !is.na(r$df2)) sprintf("(%.3g, %.3g)", r$df1, r$df2) else if (!is.na(r$df1)) sprintf("(%.3g)", r$df1) else ""
    cat(sprintf("      %s [%s %s] %s: %s%s = %.3g, p = %.3g%s%s\n", st$test, st$package, st$version, st$factor, r$statistic_name, df, r$statistic, r$p,
                if (!is.na(r$estimate)) sprintf(", %s %.3g [%.3g, %.3g]", if (is.na(r$contrast)) "estimate" else r$contrast, r$estimate, r$ci_low, r$ci_high) else "",
                if (!is.na(r$effect)) sprintf(", %s %.2f%s", r$effect_name, r$effect, if (!is.null(st$magnitude)) paste0(" (", st$magnitude, ")") else "") else ""))
    for (r in st$table[-1]) {
      if (grepl("odds ratio|relative risk", r$statistic_name)) cat(sprintf("        %s %s %.3g [%.3g, %.3g]\n", r$statistic_name, r$contrast, r$estimate, r$ci_low, r$ci_high))
      else if (grepl("goodness-of-fit|multinomial", r$statistic_name)) cat(sprintf("        %s (%s): chi-squared(%g) = %.3g, p = %.3g, Cohen's w %.2f\n", r$statistic_name, r$contrast, r$df1, r$statistic, r$p, r$effect))
    }
  }
  ph <- st$posthoc
  if (identical(ph$status, "run")) cat(sprintf("      post-hoc%s %s%s: %d comparisons (%s), %d with adjusted p < 0.05\n",
    if (!is.null(ph$version)) sprintf(" [emmeans %s]", ph$version) else "", paste(unlist(ph$families), collapse = ", "),
    if (identical(ph$type, "vs_control")) sprintf(" vs %s", paste(unlist(ph$control), collapse = ", ")) else "", ph$n_comparisons, ph$adjust, ph$n_significant))
  if (!is.null(ph$recovery)) cat(sprintf("      post-hoc truth: %d/%d true differences in their 95%% CI; %d/%d non-zero differences found, %d/%d zero differences called\n",
    ph$recovery$covered, ph$recovery$n, ph$recovery$detected, ph$recovery$true_nonzero, ph$recovery$false_calls, ph$recovery$true_zero))
  ws <- st$truth_within_subject
  if (!is.null(ws) && length(st$sphericity)) cat(sprintf("      sphericity truth: Simitall planted %s%s; Mauchly p = %s, GG epsilon = %s\n", ws$structure,
    if (!is.null(ws$rho)) sprintf(" (rho %g)", ws$rho) else if (!is.null(ws$slope_sd)) sprintf(" (slope sd %g)", ws$slope_sd) else "",
    paste(signif(vapply(st$sphericity, `[[`, 0, "mauchly_p"), 3), collapse = ", "), paste(signif(vapply(st$sphericity, `[[`, 0, "gg_epsilon"), 3), collapse = ", ")))
  w <- c(st$checks$warnings, st$note, ph$note)
  if (length(w)) cat(sprintf("      checks: %s\n", paste(w, collapse = "; ")))
  invisible()
}

.ba_print_survival <- function(st) {
  f <- function(x) if (is.na(x)) "NA" else format(signif(x, 4))
  km <- vapply(st$km, function(r) sprintf("%s %s [%s, %s] (%d/%d events)", r$group, f(r$median), f(r$median_lcl), f(r$median_ucl), as.integer(r$events), as.integer(r$n)), "")
  cat(sprintf("      Kaplan-Meier [survival %s] median: %s\n", st$version, paste(km, collapse = "; ")))
  for (r in st$table) {
    if (grepl("log-rank", r$statistic_name)) cat(sprintf("      %s(%g) = %.3g, p = %.3g\n", r$statistic_name, r$df1, r$statistic, r$p))
    else if (grepl("hazard ratio \\(Cox", r$statistic_name)) cat(sprintf("      Cox HR %s %.3g [%.3g, %.3g], p = %.3g\n", r$contrast, r$estimate, r$ci_low, r$ci_high, r$p))
    else if (grepl("^RMST", r$statistic_name)) cat(sprintf("      %s (%s): %s %.4g [%.4g, %.4g], p = %.3g; ratio %.3g\n", r$statistic_name, st$rmst$tau_source, r$contrast, r$estimate, r$ci_low, r$ci_high, r$p, r$effect))
    else if (grepl("time ratio", r$statistic_name)) cat(sprintf("      %s %s %.3g [%.3g, %.3g], p = %.3g\n", r$statistic_name, r$contrast, r$estimate, r$ci_low, r$ci_high, r$p))
  }
  if (!is.null(st$cox)) cat(sprintf("      Cox: concordance %.3f, LR test p = %.3g%s\n", st$cox$concordance, st$cox$lr_test[[3]],
    if (length(st$ph_check)) sprintf("; PH check (cox.zph) global p = %.3g", st$ph_check[[length(st$ph_check)]]$p) else ""))
  ph <- st$posthoc
  if (identical(ph$status, "run")) cat(sprintf("      pairwise log-rank (%s%s): %d comparisons, %d with adjusted p < 0.05\n", ph$adjust,
    if (identical(ph$type, "vs_control")) ", vs control" else "", ph$n_comparisons, ph$n_significant))
  if (!is.null(st$note)) cat(sprintf("      checks: %s\n", st$note))
  invisible()
}

.ba_print_competing <- function(st) {
  cat(sprintf("      competing risks [cmprsk %s], event of interest %s, competing: %s\n", st$version, st$event_of_interest, paste(setdiff(unlist(st$causes), st$event_of_interest), collapse = ", ")))
  for (r in st$table) cat(if (!is.na(r$estimate)) sprintf("      %s: %s %.3g [%.3g, %.3g], p = %.3g\n", r$statistic_name, r$contrast, r$estimate, r$ci_low, r$ci_high, r$p)
                          else sprintf("      %s(%g) = %.3g, p = %.3g\n", r$statistic_name, r$df1, r$statistic, r$p))
  cat(sprintf("      cumulative incidence vs 1 - KM: %s\n",
              paste(vapply(st$cif_vs_naive, function(r) sprintf("%s at %g %.3f vs %.3f", r$group, r$time, r$cif, r$one_minus_km), ""), collapse = "; ")))
  if (!is.null(st$note)) cat(sprintf("      checks: %s\n", st$note))
  invisible()
}
