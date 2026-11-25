# Real test data live in github.com/nirwan1265/bioana-testdata (not shipped with the package).
# Point BIOANA_TESTDATA at a clone, or clone it next to the bioana folder; tests skip otherwise.
testdata_dir <- function() {
  cands <- c(Sys.getenv("BIOANA_TESTDATA"), file.path(testthat::test_path("..", "..", ".."), "bioana-testdata"))
  hit <- cands[nzchar(cands) & dir.exists(cands)]
  if (length(hit)) normalizePath(hit[1]) else ""
}
skip_without_testdata <- function(sub) {
  td <- testdata_dir()
  if (!nzchar(td) || !dir.exists(file.path(td, sub))) testthat::skip(paste("bioana-testdata not found:", sub))
  file.path(td, sub)
}

# A design folder written the way Biodes' describe-my-data mode should (docs/contracts/biodes_from_data.md).
real_design <- function(meta, factor_name, levels, assay, model, terms, extra = list()) {
  d <- file.path(tempfile("real_"), "design"); dir.create(d, recursive = TRUE)
  spec <- c(list(handoff_version = "0.5.1", origin = "observed", organism = "test", system = "lab",
                 factors = list(list(name = factor_name, n_levels = length(levels), levels = as.list(levels))),
                 design = list(chosen = "crd", reps = 2L), assays = list(assay),
                 analysis = list(model = model, terms = terms)), extra)
  jsonlite::write_json(spec, file.path(d, "study_spec.json"), auto_unbox = TRUE, pretty = TRUE)
  lay <- meta[, setdiff(names(meta), "primary_outcome"), drop = FALSE]
  utils::write.csv(lay, file.path(d, "field_layout.csv"), row.names = FALSE)
  if (is.null(meta$primary_outcome)) meta$primary_outcome <- ""
  utils::write.table(meta, file.path(d, "sample_metadata_template.tsv"), sep = "\t", row.names = FALSE, quote = FALSE, na = "")
  d
}
