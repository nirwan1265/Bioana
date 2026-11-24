# Real example data sets, shipped in inst/extdata/real/<name>/, each built into a design folder
# the way Biodes' "describe my data" mode should write it (docs/contracts/biodes_from_data.md).

#' Build a design folder from a real example data set
#'
#' `"airway"`: human airway smooth-muscle cells. Four cell lines (donors), each
#' left untreated or treated with dexamethasone (1 uM, 18 h), RNA-seq gene
#' counts (Himes et al. 2014, PLoS One 9:e99625; GEO GSE52778; Bioconductor
#' package airway, LGPL; genes detected in at least 2 samples). It is the data
#' set of the DESeq2 vignette. The design is a randomised complete block: cell
#' line = block, dexamethasone = treatment, untreated = reference, so
#' `analysis.model` is `y ~ treatment + block` (the vignette's `~ cell + dex`).
#'
#' @param name Example name: `"airway"`.
#' @param dir Folder to write (created; must not hold a design already).
#' @return The folder path, invisibly. It holds `study_spec.json`
#'   (`origin = "observed"`), `field_layout.csv`, `sample_metadata_template.tsv`,
#'   the data under `data/`, and `data_manifest.json`.
#' @examples
#' \dontrun{
#' d <- bioana_example("airway")
#' r <- bioana_analyze(d, secondary = TRUE)
#' }
#' @export
bioana_example <- function(name = "airway", dir = file.path(tempfile("bioana_example_"), name)) {
  name <- match.arg(name, c("airway"))
  src <- system.file("extdata", "real", name, package = "bioana")
  if (!nzchar(src)) stop("Example data not found in the installed package: ", name, call. = FALSE)
  if (file.exists(file.path(dir, "study_spec.json"))) stop("A design folder already exists at ", dir, call. = FALSE)
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  s <- utils::read.delim(file.path(src, "samples.tsv"), stringsAsFactors = FALSE)
  spec <- list(
    handoff_version = "0.5.1", origin = "observed", organism = "Homo sapiens", system = "cell culture",
    question = "Which genes respond to dexamethasone in airway smooth-muscle cells?",
    factors = list(list(name = "dex", n_levels = 2L, levels = list("untrt", "trt"))),
    design = list(chosen = "rcbd", reps = 4L),
    experimental_unit = "cell culture (one per cell line and treatment)",
    primary_outcome = list(name = "gene expression", type = "continuous"),
    assays = list("rnaseq"), assay_kinds = list("rnaseq"),
    analysis = list(model = "y ~ treatment + block",
                    terms = list(treatment = "dex", block = "cell", layout_columns = list("dex", "cell"))),
    source = list(citation = "Himes BE et al. (2014) PLoS One 9:e99625", geo = "GSE52778",
                  package = "Bioconductor airway (LGPL)", note = "the DESeq2 vignette data set"))
  jsonlite::write_json(spec, file.path(dir, "study_spec.json"), auto_unbox = TRUE, pretty = TRUE)
  layout <- data.frame(unit = s$sample, cell = s$cell_line, dex = s$group, stringsAsFactors = FALSE)
  utils::write.csv(layout, file.path(dir, "field_layout.csv"), row.names = FALSE)
  meta <- data.frame(sample_id = s$sample, unit = s$sample, cell = s$cell_line, dex = s$group, run = s$run,
                     primary_outcome = NA, stringsAsFactors = FALSE)
  utils::write.table(meta, file.path(dir, "sample_metadata_template.tsv"), sep = "\t", row.names = FALSE, quote = FALSE, na = "")
  raw <- file.path(dir, "raw"); dir.create(raw)
  file.copy(file.path(src, "counts.tsv.gz"), file.path(raw, "airway_counts.tsv.gz"))
  utils::capture.output(bioana_manifest_from_files(dir, file.path(raw, "airway_counts.tsv.gz"), assay = "rnaseq", append = FALSE))
  invisible(dir)
}
