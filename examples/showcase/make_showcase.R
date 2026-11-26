# Bioana showcase: analyse Simitall's three showcase folders and keep the main tables here.
#
#   Rscript examples/showcase/make_showcase.R          # from the Bioana repository root
#   SIMITALL_SHOWCASE=/path/to/simitall/examples/showcase Rscript examples/showcase/make_showcase.R
#
# Simitall's folders are read from $SIMITALL_SHOWCASE, else ../simitall/examples/showcase, else GitHub. Each is
# copied to a temporary folder (nothing is written into Simitall's copy), analysed with bioana_analyze() and
# the default engines (rrBLUP, DESeq2, coxme / survival; edgeR as the second engine for RNA-seq), and the main
# tables, the printed report and the truth scoring are copied to examples/showcase/<name>/results/.

if (file.exists("DESCRIPTION") && grepl("^Package: bioana", readLines("DESCRIPTION", 1))) {
  suppressMessages(pkgload::load_all(".", quiet = TRUE, export_all = FALSE))
} else suppressMessages(library(bioana))

`%||%` <- function(x, y) if (is.null(x)) y else x
here <- normalizePath(if (dir.exists("examples/showcase")) "examples/showcase" else ".", mustWork = TRUE)
folders <- c("01_gwas_sorghum", "02_rnaseq_drought", "03_germination_survival")
files <- list(
  "01_gwas_sorghum" = c("study_spec.json", "field_layout.csv", "sample_metadata_template.tsv", "data_manifest.json", "sim_truth.json",
                        "genotypes.vcf.gz", "genetic_truth/qtn_truth.tsv", "genetic_truth/genetic_values.tsv", "genetic_truth/population_structure.tsv"),
  "02_rnaseq_drought" = c("study_spec.json", "field_layout.csv", "sample_metadata_template.tsv", "data_manifest.json", "sim_truth.json",
                          "counts.tsv", "truth/rnaseq.tsv", "truth/rnaseq_goal.tsv"),
  "03_germination_survival" = c("study_spec.json", "field_layout.csv", "sample_metadata_template.tsv", "data_manifest.json", "sim_truth.json"))
keep <- list(
  "01_gwas_sorghum" = c("genetics/results.tsv", "genetics/qtn_recovery.tsv", "genetics/gwas.benchmark_metrics.tsv", "genetics/line_means.tsv"),
  "02_rnaseq_drought" = c("rnaseq/results.tsv"),   # edgeR: in interaction_calls.tsv
  "03_germination_survival" = c("outcome/results.tsv", "outcome/stat_tests.tsv", "outcome/survival_summary.tsv", "outcome/landmarks.tsv",
                                "outcome/km_curve.tsv", "outcome/cif_curve.tsv", "outcome/competing_table.tsv", "outcome/posthoc.tsv",
                                "outcome/survival_group_recovery.tsv"))

src <- Sys.getenv("SIMITALL_SHOWCASE", "")
if (!nzchar(src)) src <- file.path(dirname(dirname(dirname(here))), "simitall", "examples", "showcase")
from_github <- !dir.exists(src)
raw <- "https://raw.githubusercontent.com/nirwan1265/simitall/main/examples/showcase"

for (f in folders) {
  work <- file.path(tempfile("showcase_"), f); dir.create(work, recursive = TRUE)
  for (x in files[[f]]) {
    dir.create(dirname(file.path(work, x)), recursive = TRUE, showWarnings = FALSE)
    if (from_github) utils::download.file(file.path(raw, f, x), file.path(work, x), quiet = TRUE, mode = "wb")
    else file.copy(file.path(src, f, x), file.path(work, x))
  }
  report <- utils::capture.output(r <- bioana_analyze(work, secondary = f == "02_rnaseq_drought"))
  out <- file.path(here, f, "results"); dir.create(out, recursive = TRUE, showWarnings = FALSE)
  for (x in keep[[f]]) if (file.exists(file.path(work, "bioana", x))) file.copy(file.path(work, "bioana", x), file.path(out, basename(x)), overwrite = TRUE)
  writeLines(report, file.path(out, "report.txt"))
  file.copy(file.path(work, "bioana", "bioana_summary.json"), file.path(out, "bioana_summary.json"), overwrite = TRUE)
  rec <- lapply(r$summary$analyses, function(a) a$goal_recovery %||% a$recovery)
  jsonlite::write_json(Filter(Negate(is.null), rec), file.path(out, "truth_recovery.json"), auto_unbox = TRUE, pretty = TRUE, digits = 6, null = "null")
  if (f == "02_rnaseq_drought") {
    # gene-level calls against the truth, for both engines (for Biopaint)
    g <- utils::read.delim(file.path(work, "truth", "rnaseq_goal.tsv")); dr <- utils::read.delim(file.path(work, "truth", "rnaseq.tsv"))
    de <- utils::read.delim(file.path(out, "results.tsv")); ed <- utils::read.delim(file.path(work, "bioana", "rnaseq", "results_edger.tsv"))
    calls <- data.frame(feature_id = de$feature_id,
                        planted = ifelse(de$feature_id %in% g$feature_id, "genotype x drought (goal)", ifelse(de$feature_id %in% dr$feature_id, "drought only", "none")),
                        true_interaction_log2 = g$effect[match(de$feature_id, g$feature_id)],
                        deseq2_estimate = de$estimate, deseq2_padj = de$padj, edger_estimate = ed$estimate[match(de$feature_id, ed$feature_id)],
                        edger_padj = ed$padj[match(de$feature_id, ed$feature_id)], stringsAsFactors = FALSE)
    utils::write.table(calls, file.path(out, "interaction_calls.tsv"), sep = "\t", row.names = FALSE, quote = FALSE, na = "")
  }
  cat(sprintf("%s: %s\n", f, out)); cat(report[-1], sep = "\n")
}
