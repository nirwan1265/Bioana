# Packages Bioana uses: inst/rules/requirements.tsv (package, source, profile, used_for).
# Profiles: core (installed with bioana), engines (the published packages that fit the data),
# readers (BAM, Salmon), genetics (GWAS, genomic prediction), plots, dev (tests), suite (Simitall).

.ba_profiles <- list(
  minimal = "core",
  analysis = c("core", "engines", "readers", "genetics"),
  full = c("core", "engines", "readers", "genetics", "plots", "dev"),
  suite = c("core", "engines", "readers", "genetics", "plots", "dev", "suite"))

.ba_requirements <- function(profile = "full") {
  r <- .ba_tsv("requirements")
  groups <- if (profile %in% names(.ba_profiles)) .ba_profiles[[profile]] else profile
  bad <- setdiff(groups, unique(r$profile))
  if (length(bad)) stop("Unknown profile: ", paste(bad, collapse = ", "), ". Use ", paste(names(.ba_profiles), collapse = ", "),
                        " or groups ", paste(unique(r$profile), collapse = ", "), ".", call. = FALSE)
  r[r$profile %in% groups, , drop = FALSE]
}

#' Check which packages Bioana uses are installed
#'
#' @param profile `"minimal"`, `"analysis"` (default: everything that fits or
#'   reads data), `"full"` (plus plots and tests), `"suite"` (plus Simitall), or
#'   one or more groups of `inst/rules/requirements.tsv`.
#' @return The requirements table with `installed` and `version`.
#' @export
bioana_check_dependencies <- function(profile = "analysis") {
  r <- .ba_requirements(profile)
  r$installed <- vapply(r$package, requireNamespace, logical(1), quietly = TRUE)
  r$version <- vapply(r$package, function(p) tryCatch(as.character(utils::packageVersion(p)), error = function(e) NA_character_), "")
  r
}

#' Install the packages Bioana uses
#'
#' Installs what is missing from CRAN (`install.packages`), Bioconductor
#' (`BiocManager::install`) and GitHub (`remotes::install_github`).
#'
#' @inheritParams bioana_check_dependencies
#' @param update Also update packages that are installed.
#' @return Invisibly, [bioana_check_dependencies()] after installing.
#' @export
bioana_install_dependencies <- function(profile = "analysis", update = FALSE) {
  r <- bioana_check_dependencies(profile)
  todo <- if (update) r else r[!r$installed, , drop = FALSE]
  if (!nrow(todo)) { message("All ", nrow(r), " packages for profile '", profile, "' are installed."); return(invisible(r)) }
  repos <- getOption("repos"); if (is.null(repos) || identical(unname(repos["CRAN"]), "@CRAN@")) options(repos = c(CRAN = "https://cloud.r-project.org"))
  cran <- todo$package[todo$source == "CRAN"]
  bioc <- todo$package[todo$source == "Bioconductor"]
  gh <- sub("^GitHub:", "", todo$source[startsWith(todo$source, "GitHub:")])
  if (length(cran)) utils::install.packages(cran)
  if (length(bioc)) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) utils::install.packages("BiocManager")
    BiocManager::install(bioc, update = FALSE, ask = FALSE)
  }
  if (length(gh)) {
    if (!requireNamespace("remotes", quietly = TRUE)) utils::install.packages("remotes")
    for (g in gh) remotes::install_github(g, upgrade = "never")
  }
  out <- bioana_check_dependencies(profile)
  if (any(!out$installed)) warning("Still missing: ", paste(out$package[!out$installed], collapse = ", "), call. = FALSE)
  invisible(out)
}
