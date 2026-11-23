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
#' @return Invisibly, a list: `summary` (per data set), `results`, `paths`.
#' @examples
#' \dontrun{
#' simitall::simitall_from_biodes("~/biodes_runs/my_design")   # or fill it with your data
#' bioana_analyze("~/biodes_runs/my_design")
#' }
#' @export
bioana_analyze <- function(dir, outdir = file.path(dir, "bioana"), alpha = 0.05, params = NULL, only = NULL) {
  design <- .ba_read_biodes(dir)
  man <- bioana_manifest(dir)
  if (!length(man$data)) {
    stop(sprintf("Nothing to analyse in %s: no data_manifest.json, an empty primary_outcome column and no counts.tsv. Fill the folder with data (simitall::simitall_from_biodes() or your measurements) and list it with bioana_write_manifest().", dir), call. = FALSE)
  }
  tab <- .ba_tsv("analyses")
  defaults <- .ba_tsv("method_parameters")
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  summ <- list(); results <- list(); paths <- character(); skipped <- list()
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
    fn <- get(row$method, envir = asNamespace("bioana"))
    r <- tryCatch(fn(e, design, file.path(outdir, e$id), alpha, args), error = function(err) err)
    if (inherits(r, "error")) {
      skipped[[e$id]] <- list(data_type = e$data_type, status = "error", note = conditionMessage(r))
      next
    }
    summ[[e$id]] <- r$summary
    results[[e$id]] <- r$results
    paths <- c(paths, stats::setNames(r$paths, paste0(e$id, "/", names(r$paths))))
  }
  out <- list(bioana_version = as.character(utils::packageVersion("bioana")), design = design$spec$design$chosen,
              model = design$model, handoff_version = design$spec$handoff_version,
              manifest = man$written_by %||% "", alpha = alpha, analyses = summ, not_analysed = skipped)
  paths <- c(paths, summary = .ba_write_json(out, file.path(outdir, "bioana_summary.json")))
  .ba_report(out)
  invisible(list(summary = out, results = results, paths = paths))
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
    } else {
      sig <- unlist(a$significant)
      sprintf("%d features tested, %s significant", a$n_features_tested %||% 0L,
              if (length(sig)) paste(paste0(names(sig), " ", sig), collapse = ", ") else "0")
    }
    rec <- a$recovery
    rtxt <- if (is.null(rec)) "" else if (a$data_type == "outcome") {
      if (!is.null(rec$coverage)) sprintf("; truth: %d/%d coefficients in their CI", rec$covered, rec$n_terms) else
        sprintf("; truth: sign %s, detected %s", rec$sign_agreement, rec$detected)
    } else {
      paste0("; truth: ", paste(vapply(rec, function(r) sprintf("%s found %d/%d (FDR %s)", r$contrast, r$true_positives, r$true_changes, r$fdr), ""), collapse = ", "))
    }
    cat(sprintf("  [%s] %s%s\n", id, line, rtxt))
  }
  for (id in names(s$not_analysed)) cat(sprintf("  [%s] not analysed (%s): %s\n", id, s$not_analysed[[id]]$status, s$not_analysed[[id]]$note))
}
