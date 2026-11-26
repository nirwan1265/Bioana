# CDISC ADaM ADTTE (time-to-event analysis data set) -> design folder.
# ADTTE: one row per subject and parameter; AVAL = time from the origin, CNSR = 1 censored, 0 event
# (the reverse of Surv(time, event)), PARAMCD = endpoint code (OS, PFS, ...), TRT01P / ARM = arm.

#' Build a design folder from a CDISC ADTTE data set
#'
#' Reads an ADTTE table (`.csv`, `.tsv`, or `.xpt`/`.sas7bdat` with the haven
#' package), keeps one endpoint (`PARAMCD`), turns the CDISC censoring flag into
#' an event indicator (`event = 1 - CNSR`: CDISC codes 1 = censored, `Surv()`
#' wants 1 = event) and writes a design folder for [bioana_analyze()]: the arm
#' is the treatment, the time is `AVAL` in `AVALU` units.
#'
#' @param file Path to the ADTTE file.
#' @param dir Folder to write.
#' @param paramcd Endpoint to keep (e.g. `"OS"`); required when the file has several.
#' @param arm Arm column: `"TRT01P"` (planned), `"TRT01A"` (actual) or `"ARM"`; the first present by default.
#' @param reference Reference (control) arm; default the first arm in sorted order.
#' @param strata Optional stratification columns (e.g. `"STRATA1"`) for a stratified log-rank.
#' @return The folder path, invisibly.
#' @export
bioana_adtte_design <- function(file, dir, paramcd = NULL, arm = NULL, reference = NULL, strata = NULL) {
  x <- .ba_read_table_any(file)
  names(x) <- toupper(names(x))
  need <- c("USUBJID", "PARAMCD", "AVAL", "CNSR")
  miss <- setdiff(need, names(x))
  if (length(miss)) stop(sprintf("Not an ADTTE data set: missing %s.", paste(miss, collapse = ", ")), call. = FALSE)
  pc <- unique(as.character(x$PARAMCD))
  if (is.null(paramcd)) {
    if (length(pc) > 1) stop(sprintf("ADTTE has several endpoints (PARAMCD %s): choose one with paramcd.", paste(pc, collapse = ", ")), call. = FALSE)
    paramcd <- pc
  }
  x <- x[as.character(x$PARAMCD) == paramcd, , drop = FALSE]
  if (!nrow(x)) stop(sprintf("No rows with PARAMCD '%s' (present: %s).", paramcd, paste(pc, collapse = ", ")), call. = FALSE)
  if (anyDuplicated(x$USUBJID)) stop(sprintf("Several rows per subject for PARAMCD '%s': ADTTE should have one.", paramcd), call. = FALSE)
  cn <- suppressWarnings(as.numeric(x$CNSR))
  if (any(!cn %in% c(0, 1))) stop("CNSR must be 0 (event) or 1 (censored).", call. = FALSE)
  ev <- 1L - as.integer(cn)
  stopifnot(sum(ev) == sum(cn == 0))                       # the check the CDISC conversion needs
  armc <- toupper(arm %||% intersect(c("TRT01P", "TRT01A", "ARM", "ACTARM"), names(x))[1])
  if (is.na(armc) || !armc %in% names(x)) stop("No arm column (TRT01P, TRT01A, ARM): give arm.", call. = FALSE)
  arms <- sort(unique(as.character(x[[armc]])))
  if (!is.null(reference)) { if (!reference %in% arms) stop(sprintf("Reference arm '%s' is not an arm (%s).", reference, paste(arms, collapse = ", ")), call. = FALSE)
    arms <- c(reference, setdiff(arms, reference)) }
  unit <- if ("AVALU" %in% names(x)) as.character(x$AVALU[1]) else "days"
  param <- if ("PARAM" %in% names(x)) as.character(x$PARAM[1]) else paramcd
  d <- data.frame(unit = as.character(x$USUBJID), arm = as.character(x[[armc]]), time = as.numeric(x$AVAL), event = ev, sample_id = as.character(x$USUBJID),
                  stringsAsFactors = FALSE)
  st <- toupper(unlist(strata))
  for (s in st) { if (!s %in% names(x)) stop(sprintf("Strata column '%s' is not in the ADTTE.", s), call. = FALSE); d[[s]] <- as.character(x[[s]]) }
  if (file.exists(file.path(dir, "study_spec.json"))) stop("A design folder already exists at ", dir, call. = FALSE)
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  an <- list(model = "y ~ treatment", terms = list(treatment = "arm", layout_columns = list("arm")))
  if (length(st)) an$survival <- list(strata = as.list(st))
  spec <- list(handoff_version = "0.6.5", origin = "observed", organism = "Homo sapiens", system = "clinical trial (CDISC ADTTE)",
               question = sprintf("%s by arm", param), experimental_unit = "subject",
               factors = list(list(name = "arm", n_levels = length(arms), levels = as.list(arms))), design = list(chosen = "parallel"),
               primary_outcome = list(name = param, type = "survival", unit = unit, event_column = "event", paramcd = paramcd,
                                      note = "event = 1 - CNSR (CDISC: CNSR 1 = censored)"),
               analysis = an, source = list(file = basename(file), format = "CDISC ADaM ADTTE"))
  jsonlite::write_json(spec, file.path(dir, "study_spec.json"), auto_unbox = TRUE, pretty = TRUE)
  utils::write.csv(d[, c("unit", "arm", st), drop = FALSE], file.path(dir, "field_layout.csv"), row.names = FALSE)
  utils::write.table(d, file.path(dir, "sample_metadata_template.tsv"), sep = "\t", row.names = FALSE, quote = FALSE, na = "")
  bioana_write_manifest(dir, list(id = "outcome", data_type = "outcome", files = list(table = "sample_metadata_template.tsv"),
                                  column = "time", event_column = "event"))
  invisible(dir)
}

.ba_read_table_any <- function(file) {
  ext <- tolower(tools::file_ext(sub("\\.gz$", "", file)))
  if (ext %in% c("xpt", "sas7bdat")) {
    if (!requireNamespace("haven", quietly = TRUE)) stop("Reading .xpt / .sas7bdat needs the haven package.", call. = FALSE)
    return(as.data.frame(if (ext == "xpt") haven::read_xpt(file) else haven::read_sas(file)))
  }
  utils::read.table(file, header = TRUE, sep = if (ext == "csv") "," else "\t", quote = "\"", stringsAsFactors = FALSE, check.names = FALSE, comment.char = "")
}
