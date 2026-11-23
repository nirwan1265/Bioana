# data_manifest.json: what data a design folder holds, so Bioana knows what to run.
#
# Simitall writes it next to the data it simulates; an experimenter writes it
# (bioana_write_manifest()) for measured data. One entry per data set:
#   id         short name, also the results subfolder (bioana/<id>/)
#   data_type  a row of inst/rules/analyses.tsv (outcome, feature_counts, ...)
#   assay      the Biodes assay it came from (rnaseq, chip_seq, ...), for reports
#   files      named paths relative to the folder (which names: analyses.tsv `files`)
#   truth      optional, simulated data only: {kind, file, term, large}
#              kind features: file with feature_id, condition, effect
#              kind regions:  file with region_id, chr, start, end, condition, effect
#              kind outcome:  file = sim_truth.json
# Columns of every matrix are the sample_id values of sample_metadata_template.tsv.

#' Read (or infer) the data manifest of a design folder
#'
#' Reads `data_manifest.json`. A folder without one (an older Simitall run,
#' or measured data filled in by hand) gets an inferred manifest: the primary
#' outcome when its column has values, and `counts.tsv` as gene counts.
#'
#' @param dir Biodes design folder.
#' @return A list with `manifest_version`, `written_by` and `data` (one list per data set).
#' @export
bioana_manifest <- function(dir) {
  p <- file.path(dir, "data_manifest.json")
  if (file.exists(p)) {
    m <- jsonlite::read_json(p, simplifyVector = FALSE)
    m$data <- m$data %||% list()
    return(m)
  }
  data <- list()
  meta_p <- file.path(dir, "sample_metadata_template.tsv")
  if (file.exists(meta_p)) {
    meta <- utils::read.delim(meta_p, colClasses = "character", na.strings = c("", "NA"), check.names = FALSE)
    if (!is.null(meta$primary_outcome) && any(!is.na(meta$primary_outcome))) {
      spec <- tryCatch(jsonlite::read_json(file.path(dir, "study_spec.json")), error = function(e) list())
      e <- list(id = "outcome", data_type = "outcome", files = list(table = "sample_metadata_template.tsv"),
                column = "primary_outcome", outcome_type = spec$primary_outcome$type %||% "continuous")
      if (!is.null(meta$event)) e$event_column <- "event"
      if (file.exists(file.path(dir, "sim_truth.json"))) e$truth <- list(kind = "outcome", file = "sim_truth.json")
      data[[length(data) + 1]] <- e
    }
  }
  if (file.exists(file.path(dir, "counts.tsv"))) {
    data[[length(data) + 1]] <- list(id = "counts", data_type = "feature_counts", assay = "rnaseq", files = list(counts = "counts.tsv"))
  }
  list(manifest_version = "0.1", written_by = "inferred by bioana (no data_manifest.json)", data = data)
}

#' Write a data manifest for measured data
#'
#' For a design folder you filled with your own measurements: list what is
#' there so [bioana_analyze()] runs the right analyses.
#'
#' @param dir Biodes design folder.
#' @param ... One list per data set, e.g.
#'   `list(id = "rnaseq", data_type = "feature_counts", assay = "rnaseq", files = list(counts = "counts.tsv"))`.
#'   See `bioana_analyses()` for the data types and the files each needs.
#' @return The path written, invisibly.
#' @examples
#' d <- tempfile(); dir.create(d)
#' bioana_write_manifest(d,
#'   list(id = "outcome", data_type = "outcome", files = list(table = "sample_metadata_template.tsv"),
#'        column = "primary_outcome", outcome_type = "continuous"),
#'   list(id = "rnaseq", data_type = "feature_counts", assay = "rnaseq", files = list(counts = "counts.tsv")))
#' @export
bioana_write_manifest <- function(dir, ...) {
  data <- list(...)
  tab <- .ba_tsv("analyses")
  for (e in data) {
    if (is.null(e$id) || is.null(e$data_type)) stop("Each data set needs an id and a data_type.", call. = FALSE)
    if (!e$data_type %in% tab$data_type) stop(sprintf("Unknown data_type '%s'; see bioana_analyses().", e$data_type), call. = FALSE)
    need <- strsplit(tab$files[tab$data_type == e$data_type], ",", fixed = TRUE)[[1]]
    miss <- setdiff(need[nzchar(need)], names(e$files))
    if (length(miss)) stop(sprintf("Data set '%s' (%s) needs files: %s", e$id, e$data_type, paste(miss, collapse = ", ")), call. = FALSE)
  }
  p <- file.path(dir, "data_manifest.json")
  jsonlite::write_json(list(manifest_version = "0.1", written_by = "user", data = data), p, auto_unbox = TRUE, pretty = TRUE)
  invisible(p)
}

#' The analyses Bioana runs, by data type
#'
#' @return `inst/rules/analyses.tsv`: one row per data type with the files it
#'   needs, the method, its status and what it writes.
#' @export
bioana_analyses <- function() .ba_tsv("analyses")

# Truth of one manifest entry, ready for scoring (NULL for measured data).
.ba_entry_truth <- function(entry, dir) {
  tr <- entry$truth
  if (is.null(tr) || is.null(tr$file) || !file.exists(file.path(dir, tr$file))) return(NULL)
  f <- file.path(dir, tr$file)
  if (identical(tr$kind, "outcome")) return(jsonlite::read_json(f, simplifyVector = FALSE))
  x <- utils::read.delim(f, stringsAsFactors = FALSE, check.names = FALSE)
  if (!nrow(x)) x <- data.frame(feature_id = character(), condition = character(), effect = numeric())
  x$contrast <- paste0(tr$term %||% "", x$condition)
  list(kind = tr$kind %||% "features", table = x, large = tr$large)
}
