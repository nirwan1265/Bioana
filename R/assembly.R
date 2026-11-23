# Genome assembly evaluation (moved from Simitall): Unicycler assemblies and QUAST evaluation over read-coverage grids.

#' Summarize QUAST reports
#'
#' Combine all `report.tsv` files beneath a QUAST result directory into one
#' tidy CSV containing coverage levels and key assembly metrics.
#'
#' @param quast_root Directory containing QUAST result folders.
#' @param out_csv Output CSV path.
#'
#' @return The summary data frame, invisibly.
#' @examples
#' \dontrun{
#' summarize_quast("results/quast", "results/quast_summary.csv")
#' }
#' @export
summarize_quast <- function(quast_root, out_csv) {
  reports <- list.files(
    quast_root,
    pattern = "report\\.tsv$",
    recursive = TRUE,
    full.names = TRUE
  )
  if (!length(reports)) {
    stop("No report.tsv files found under: ", quast_root)
  }

  rows <- lapply(reports, function(path) {
    report <- utils::read.delim(
      path,
      header = TRUE,
      sep = "\t",
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    values <- stats::setNames(report[[2L]], report[[1L]])
    get_metric <- function(name) {
      if (name %in% names(values)) as.numeric(values[[name]]) else NA_real_
    }

    folder <- basename(dirname(path))
    match <- regexec("^ill([0-9]+)_pb([0-9]+)$", folder)
    coverage <- regmatches(folder, match)[[1L]]

    data.frame(
      ill_cov = if (length(coverage)) as.integer(coverage[2L]) else NA_integer_,
      pb_cov = if (length(coverage)) as.integer(coverage[3L]) else NA_integer_,
      total_length = get_metric("Total length"),
      n_contigs = get_metric("# contigs"),
      n50 = get_metric("N50"),
      largest_contig = get_metric("Largest contig"),
      gc = get_metric("GC (%)"),
      genome_fraction = get_metric("Genome fraction (%)"),
      mismatches_per_100kbp = get_metric("# mismatches per 100 kbp"),
      indels_per_100kbp = get_metric("# indels per 100 kbp"),
      misassemblies = get_metric("# misassemblies"),
      stringsAsFactors = FALSE
    )
  })

  result <- do.call(rbind, rows)
  result <- result[order(result$ill_cov, result$pb_cov), , drop = FALSE]
  dir.create(dirname(out_csv), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(result, out_csv, row.names = FALSE)
  invisible(result)
}

#' Run QUAST evaluation across assembly grid
#'
#' @param tag Character. Tag used for input/output folders.
#' @param ref_fa Truth reference FASTA.
#' @param assemblies_dir Root assembly directory.
#' @param output_dir Root QUAST output directory.
#' @param threads Number of QUAST threads.
#' @param ill_covs Illumina coverage values.
#' @param pb_covs PacBio coverage values.
#'
#' @return TRUE on success.
#' @examples
#' \dontrun{
#' run_quast_grid("ecoli_demo", ref_fa = "genome.fa")
#' }
#' @export
run_quast_grid <- function(
    tag,
    ref_fa,
    assemblies_dir = "assemblies",
    output_dir = "evaluation",
    threads = 4,
    ill_covs = c(10, 20, 30, 40, 60, 80),
    pb_covs = c(5, 10, 15, 20, 30, 40)) {
  quast <- Sys.which("quast.py")
  if (!nzchar(quast)) {
    quast <- Sys.which("quast")
  }
  if (!nzchar(quast)) {
    stop("QUAST was not found in PATH.")
  }
  if (!file.exists(ref_fa)) {
    stop("Reference FASTA not found: ", ref_fa)
  }

  for (mode in c("CLR", "HIFI")) {
    outroot <- file.path(output_dir, tag, mode)
    dir.create(outroot, showWarnings = FALSE, recursive = TRUE)

    for (ic in ill_covs) {
      for (pc in pb_covs) {
        asm <- file.path(assemblies_dir, tag, mode, paste0("ill", ic, "_pb", pc), "assembly.fasta")
        if (!file.exists(asm)) next
        outdir <- file.path(outroot, paste0("ill", ic, "_pb", pc))
        if (file.exists(file.path(outdir, "report.tsv"))) next

        cmd <- sprintf("%s -r %s -o %s --threads %s %s",
                       shQuote(quast), shQuote(ref_fa), shQuote(outdir),
                       threads, shQuote(asm))
        if (system(cmd) != 0L) stop("QUAST failed for ", basename(outdir))
      }
    }
  }
  TRUE
}

#' Run Unicycler hybrid assemblies across a grid
#'
#' @param tag Character. Tag used for input/output folders.
#' @param threads Integer. Number of threads.
#' @param reads_dir Root directory containing simulated reads.
#' @param output_dir Root assembly output directory.
#' @param ill_covs Illumina coverage values.
#' @param pb_covs PacBio coverage values.
#'
#' @return TRUE on success.
#' @examples
#' \dontrun{
#' run_unicycler_grid("ecoli_repMed", threads = 8)
#' }
#' @export
run_unicycler_grid <- function(
    tag,
    threads = 8,
    reads_dir = "reads",
    output_dir = "assemblies",
    ill_covs = c(10, 20, 30, 40, 60, 80),
    pb_covs = c(5, 10, 15, 20, 30, 40)) {
  if (!nzchar(Sys.which("unicycler"))) {
    stop("Unicycler was not found in PATH.")
  }

  for (ic in ill_covs) {
    r1 <- file.path(reads_dir, tag, "illumina", paste0(tag, ".ill_cov", ic, "_1.fq.gz"))
    r2 <- file.path(reads_dir, tag, "illumina", paste0(tag, ".ill_cov", ic, "_2.fq.gz"))

    for (pc in pb_covs) {
      pb_dir <- file.path(reads_dir, tag, "pacbio_CLR", paste0("cov", pc))
      pb_fastq <- list.files(pb_dir, pattern = paste0("pbCLR_cov", pc, ".*fastq.gz"), full.names = TRUE)
      if (length(pb_fastq) > 0) {
        out_dir <- file.path(output_dir, tag, "CLR", paste0("ill", ic, "_pb", pc))
        if (!file.exists(file.path(out_dir, "assembly.fasta"))) {
          dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
          cmd <- sprintf("unicycler -1 %s -2 %s -l %s -o %s -t %s --mode normal",
                         shQuote(r1), shQuote(r2), shQuote(pb_fastq[1]), shQuote(out_dir), threads)
          if (system(cmd) != 0L) stop("Unicycler failed for ", basename(out_dir))
        }
      }

      pb_dir <- file.path(reads_dir, tag, "pacbio_HIFI", paste0("cov", pc))
      pb_fastq <- list.files(pb_dir, pattern = paste0("pbHIFI_cov", pc, ".*fastq.gz"), full.names = TRUE)
      if (length(pb_fastq) > 0) {
        out_dir <- file.path(output_dir, tag, "HIFI", paste0("ill", ic, "_pb", pc))
        if (!file.exists(file.path(out_dir, "assembly.fasta"))) {
          dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
          cmd <- sprintf("unicycler -1 %s -2 %s -l %s -o %s -t %s --mode normal",
                         shQuote(r1), shQuote(r2), shQuote(pb_fastq[1]), shQuote(out_dir), threads)
          if (system(cmd) != 0L) stop("Unicycler failed for ", basename(out_dir))
        }
      }
    }
  }
  TRUE
}
