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
#' Survival examples (data from the survival package, Therneau, LGPL; the data sets of the survival
#' tutorials, e.g. datanovia's survival-analysis course):
#' `"lung"`: NCCTG advanced lung cancer, 228 patients, survival by sex (Loprinzi et al. 1994,
#' J Clin Oncol 12:601); landmarks at 180 and 360 days.
#' `"colon"`: adjuvant colon cancer trial, death records, three arms (Obs, Lev, Lev+5FU; Moertel et
#' al. 1990, NEJM 322:352); pairwise log-rank (BH).
#' `"mgus2"`: 1384 patients with MGUS, progression to plasma-cell malignancy with death as a competing
#' event (Kyle et al. 2002, NEJM 346:564); event codes 1 progression, 2 death.
#' `"rats"`: 300 rats in 100 litters, drug vs control, tumour time (Mantel et al. 1977, Cancer Res
#' 37:3863); litter as a random intercept (frailty).
#'
#' @param name Example name: `"airway"`, `"lung"`, `"colon"`, `"mgus2"` or `"rats"`.
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
  name <- match.arg(name, c("airway", "lung", "colon", "mgus2", "rats"))
  if (name != "airway") return(.ba_example_survival(name, dir))
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

# Survival examples from the survival package's data sets, written as design folders.
.ba_example_survival <- function(name, dir) {
  if (!requireNamespace("survival", quietly = TRUE)) stop("The survival examples need the survival package.", call. = FALSE)
  if (file.exists(file.path(dir, "study_spec.json"))) stop("A design folder already exists at ", dir, call. = FALSE)
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  get <- function(x) { e <- new.env(); suppressWarnings(utils::data(list = "cancer", package = "survival", envir = e)); if (is.null(e[[x]])) utils::data(list = x, package = "survival", envir = e); e[[x]] }
  ev_codes <- NULL; surv <- list(); posthoc <- NULL; model <- "y ~ treatment"; terms <- NULL; random_note <- NULL
  if (name == "lung") {
    x <- get("lung")
    d <- data.frame(unit = sprintf("P%03d", seq_len(nrow(x))), sex = c("male", "female")[x$sex], ph.ecog = x$ph.ecog, age = x$age,
                    time = x$time, event = x$status - 1L, stringsAsFactors = FALSE)
    fac <- list(list(name = "sex", n_levels = 2L, levels = list("male", "female")))
    terms <- list(treatment = "sex", layout_columns = list("sex"))
    surv <- list(landmarks = list(180, 360)); unit <- "days"
    q <- "Does survival differ between men and women with advanced lung cancer?"
    cite <- "Loprinzi CL et al. (1994) J Clin Oncol 12:601; survival::lung"
  } else if (name == "colon") {
    x <- get("colon"); x <- x[x$etype == 2, ]
    d <- data.frame(unit = sprintf("P%04d", x$id), rx = as.character(x$rx), time = x$time, event = x$status, stringsAsFactors = FALSE)
    fac <- list(list(name = "rx", n_levels = 3L, levels = list("Obs", "Lev", "Lev+5FU")))
    terms <- list(treatment = "rx", layout_columns = list("rx"))
    posthoc <- list(type = "pairwise", factor = "rx", adjust = "BH"); unit <- "days"
    surv <- list(landmarks = list(365, 1825))
    q <- "Does adjuvant levamisole, alone or with fluorouracil, prolong survival after colon cancer surgery?"
    cite <- "Moertel CG et al. (1990) NEJM 322:352; survival::colon (death records, etype 2)"
  } else if (name == "mgus2") {
    x <- get("mgus2")
    et <- ifelse(x$pstat == 1, x$ptime, x$futime); ev <- ifelse(x$pstat == 1, 1L, 2L * x$death)
    d <- data.frame(unit = sprintf("P%04d", x$id), sex = ifelse(x$sex == "F", "female", "male"), time = et, event = ev, stringsAsFactors = FALSE)
    fac <- list(list(name = "sex", n_levels = 2L, levels = list("female", "male")))
    terms <- list(treatment = "sex", layout_columns = list("sex"))
    ev_codes <- list(`0` = "censored", `1` = "progression", `2` = "death"); unit <- "months"
    surv <- list(landmarks = list(120, 240))
    q <- "Is the risk of progression to plasma-cell malignancy different between women and men with MGUS, with death as a competing event?"
    cite <- "Kyle RA et al. (2002) NEJM 346:564; survival::mgus2"
  } else {
    x <- get("rats")
    d <- data.frame(unit = sprintf("R%03d", seq_len(nrow(x))), rx = c("control", "drug")[x$rx + 1L], litter = sprintf("L%03d", x$litter), sex = x$sex,
                    time = x$time, event = x$status, stringsAsFactors = FALSE)
    fac <- list(list(name = "rx", n_levels = 2L, levels = list("control", "drug")))
    model <- "y ~ treatment + (1 | block)"
    terms <- list(treatment = "rx", block = "litter", layout_columns = list("rx", "litter")); unit <- "days"
    q <- "Does the drug delay tumour appearance in rats, with litter mates sharing a frailty?"
    cite <- "Mantel N et al. (1977) Cancer Res 37:3863; survival::rats"
  }
  d <- d[!is.na(d$time) & !is.na(d$event), , drop = FALSE]
  d$sample_id <- d$unit
  po <- list(name = "time to event", type = "survival", unit = unit, event_column = "event")
  if (!is.null(ev_codes)) po$event_codes <- ev_codes
  an <- list(model = model, terms = terms)
  if (length(surv)) an$survival <- surv
  if (!is.null(posthoc)) an$posthoc <- posthoc
  spec <- list(handoff_version = "0.6.5", origin = "observed", organism = if (name == "rats") "Rattus norvegicus" else "Homo sapiens",
               system = if (name == "rats") "animal study" else "clinical study", question = q, experimental_unit = if (name == "rats") "rat" else "patient",
               factors = fac, design = list(chosen = "observational"), primary_outcome = po, analysis = an,
               source = list(citation = cite, package = "survival (Therneau, LGPL)"))
  jsonlite::write_json(spec, file.path(dir, "study_spec.json"), auto_unbox = TRUE, pretty = TRUE, digits = NA)
  lay <- unique(c("unit", unlist(terms$layout_columns)))
  utils::write.csv(d[, lay, drop = FALSE], file.path(dir, "field_layout.csv"), row.names = FALSE)
  utils::write.table(d, file.path(dir, "sample_metadata_template.tsv"), sep = "\t", row.names = FALSE, quote = FALSE, na = "")
  bioana_write_manifest(dir, list(id = "outcome", data_type = "outcome", files = list(table = "sample_metadata_template.tsv"),
                                  column = "time", event_column = "event"))
  invisible(dir)
}
