# Readers: a user's files -> Bioana data sets.
#
# inst/rules/file_formats.tsv says how each format is recognised (file name
# pattern, a header line, or required columns), whether it is one file per
# sample, what it gives (counts, bisulfite, genotypes, ...) and which reader
# turns it into Bioana's own files. inst/rules/file_formats.json holds the
# rules for sample ids (what to strip from file and column names), which
# data type a count matrix becomes for each assay, and the IP / input
# column suffixes. bioana_manifest_from_files() detects, converts, and writes
# data_manifest.json; every entry records the source files, format and tool.

#' Detect the format of data files
#'
#' @param files Paths.
#' @return A data frame: `file`, `format`, `gives`, `tool`, `per_sample` (one row per file;
#'   `format = NA` when nothing in `inst/rules/file_formats.tsv` matches).
#' @export
bioana_detect_format <- function(files) {
  ff <- .ba_formats()
  out <- lapply(files, function(f) {
    if (!file.exists(f)) stop("File not found: ", f, call. = FALSE)
    head <- if (grepl("\\.(bam|cram)$", f)) character() else .ba_head_lines(f)
    first <- head[!grepl("^#", head)][1] %||% ""
    cols <- strsplit(first, if (grepl("\\.csv(\\.gz)?$", f)) "," else "\t", fixed = TRUE)[[1]]
    hit <- NA_integer_
    for (i in seq_len(nrow(ff))) {
      if (!grepl(ff$file_pattern[i], basename(f), perl = TRUE)) next
      hp <- ff$header_pattern[i]; rc <- ff$required_columns[i]
      ok_h <- nzchar(hp) && any(grepl(hp, head, perl = TRUE))
      ok_c <- nzchar(rc) && all(strsplit(rc, ",", fixed = TRUE)[[1]] %in% cols)
      if (ok_h || ok_c || (!nzchar(hp) && !nzchar(rc))) { hit <- i; break }
    }
    if (is.na(hit)) return(data.frame(file = f, format = NA_character_, gives = NA_character_, tool = NA_character_, per_sample = NA))
    data.frame(file = f, format = ff$format[hit], gives = ff$gives[hit], tool = ff$tool[hit], per_sample = ff$per_sample[hit] == "yes",
               stringsAsFactors = FALSE)
  })
  do.call(rbind, out)
}

#' Turn a user's data files into Bioana data sets and write data_manifest.json
#'
#' Detects each file's format (featureCounts, Salmon, kallisto, STAR
#' ReadsPerGene, Bismark coverage or cytosine report, VCF, MACS peak calls, a
#' plain feature x sample table), converts it to Bioana's files under
#' `<dir>/data/<id>/`, matches sample names to `sample_metadata_template.tsv`,
#' and writes the manifest. Files of a per-sample format (one Salmon folder,
#' one Bismark file per sample) form one data set together.
#'
#' @param dir A design folder (written by Biodes, `write_design_package()` or
#'   `biodes_from_data()`).
#' @param files Data file paths.
#' @param assay Optional assay label for count tables (`"rnaseq"`, `"chipseq"`,
#'   `"atacseq"`, `"merip_seq"`, ...); otherwise taken from `study_spec.json`.
#'   It decides whether a count table is analysed as features or peaks.
#' @param tx2gene Optional two-column table (transcript, gene) or path, for
#'   Salmon and kallisto: transcripts are summed to genes by tximport.
#' @param strandedness For STAR ReadsPerGene: `"auto"`, `"unstranded"`,
#'   `"forward"` or `"reverse"`.
#' @param append Keep the data sets already in the manifest.
#' @param annotation For BAM files: a GTF (genes; exon windows for IP/input
#'   assays), a SAF table, or peak BED / narrowPeak files (merged into consensus
#'   peaks present in at least 2 files).
#' @param libraries For BAM files: a data frame with `file` (BAM file name or
#'   path), `sample_id`, and `library` (`rna`, `chip`, `control`, `ip`,
#'   `input`). Without it, each BAM is one sample named after its file.
#' @param window_bp Window size for IP/input assays (exon windows from the GTF).
#' @return Invisibly, the manifest entries; prints one line per data set.
#' @export
bioana_manifest_from_files <- function(dir, files, assay = NULL, tx2gene = NULL, strandedness = "auto", append = TRUE,
                                       annotation = NULL, libraries = NULL, window_bp = 200L) {
  meta <- utils::read.delim(file.path(dir, "sample_metadata_template.tsv"), stringsAsFactors = FALSE, colClasses = "character")
  spec <- jsonlite::read_json(file.path(dir, "study_spec.json"))
  assay <- assay %||% .ba_spec_assay(spec)
  det <- bioana_detect_format(files)
  if (anyNA(det$format)) {
    stop("Format not recognised (see inst/rules/file_formats.tsv): ", paste(det$file[is.na(det$format)], collapse = ", "), call. = FALSE)
  }
  groups <- split(det, ifelse(det$per_sample, det$format, paste0(det$format, "::", det$file)))
  old <- if (append && file.exists(file.path(dir, "data_manifest.json"))) jsonlite::read_json(file.path(dir, "data_manifest.json"))$data else list()
  old_ids <- vapply(old, `[[`, "", "id")
  old_src <- lapply(old, function(o) sort(unlist(o$source$files)))
  entries <- list()
  for (g in groups) {
    ff <- .ba_formats(); row <- ff[ff$format == g$format[1], ]
    reader <- get(row$reader, envir = asNamespace("bioana"))
    x <- reader(g$file, tx2gene = tx2gene, strandedness = strandedness, annotation = annotation,
                libraries = libraries, assay = assay, window_bp = window_bp)
    # same source files as an existing entry: replace it; otherwise a new id
    same <- which(vapply(old_src, identical, logical(1), sort(normalizePath(g$file))))
    taken <- c(stats::setNames(vector("list", length(old_ids)), old_ids)[setdiff(seq_along(old_ids), same)], entries)
    e <- .ba_entry_from_read(x, g, row, dir, meta, assay, taken, id = if (length(same)) old_ids[same[1]] else NULL)
    entries[[e$id]] <- e
  }
  old <- Filter(function(o) !o$id %in% names(entries), old)
  all <- c(old, unname(entries))
  do.call(bioana_write_manifest, c(list(dir), all))
  for (e in entries) {
    cat(sprintf("  %s: %s (%s, %s) -> %s, %d samples%s\n", e$id, e$source$format, e$source$tool, e$source$reader,
                e$data_type, length(e$source$samples), if (length(e$source$warnings)) paste0("; ", paste(e$source$warnings, collapse = "; ")) else ""))
  }
  invisible(entries)
}

#' Consensus peaks from MACS peak calls, as a SAF file for featureCounts
#'
#' @param files narrowPeak or broadPeak files, one per sample.
#' @param out SAF path to write.
#' @param min_samples A region is kept when it overlaps peaks of at least this many files.
#' @return The SAF path. Count reads in it with `featureCounts -F SAF -a <out>`
#'   and give the result to [bioana_manifest_from_files()].
#' @export
bioana_consensus_peaks <- function(files, out, min_samples = 2L) {
  p <- .ba_read_peaks(files)$peaks
  p <- p[order(p$chr, p$start), ]
  merged <- do.call(rbind, lapply(split(p, p$chr), function(x) {
    grp <- cumsum(c(TRUE, x$start[-1] > cummax(x$end)[-nrow(x)]))
    data.frame(chr = x$chr[1], start = tapply(x$start, grp, min), end = tapply(x$end, grp, max),
               n = tapply(x$sample, grp, function(s) length(unique(s))))
  }))
  merged <- merged[merged$n >= min_samples, ]
  saf <- data.frame(GeneID = sprintf("peak_%s_%d_%d", merged$chr, merged$start, merged$end), Chr = merged$chr,
                    Start = merged$start, End = merged$end, Strand = "+")
  utils::write.table(saf, out, sep = "\t", row.names = FALSE, quote = FALSE)
  invisible(out)
}

# ---- helpers ----------------------------------------------------------------

.ba_formats <- function() utils::read.delim(system.file("rules", "file_formats.tsv", package = "bioana"), quote = "",
                                            colClasses = "character", stringsAsFactors = FALSE)
.ba_format_rules <- function() .ba_rules("file_formats")

.ba_head_lines <- function(f, n = 60L) {
  con <- if (grepl("\\.gz$", f)) gzfile(f, "rt") else file(f, "rt")
  on.exit(close(con))
  readLines(con, n = n, warn = FALSE)
}

.ba_sample_name <- function(x) {
  for (r in unlist(.ba_format_rules()$sample_name_strip)) x <- sub(r, "", x, perl = TRUE)
  x
}

# Sample id of a per-sample file: the file name, or its folder for formats whose file name is fixed (quant.sf).
.ba_sample_from_path <- function(f, format) {
  if (format %in% c("salmon", "kallisto")) basename(dirname(f)) else .ba_sample_name(basename(f))
}

.ba_spec_assay <- function(spec) {
  a <- tolower(unlist(c(spec$assay_kinds, spec$assays)))
  m <- names(.ba_format_rules()$assay_data_type)
  hit <- a[a %in% m]
  if (length(hit)) hit[1] else NULL
}

# Match data sample names to sample_metadata: exact sample_id, then any metadata column
# (library_id, file, bam, ...) holding the name, then names with the strip rules applied.
.ba_match_samples <- function(names_in, meta) {
  ids <- rep(NA_character_, length(names_in))
  try_cols <- c("sample_id", setdiff(names(meta), "sample_id"))
  for (col in try_cols) {
    v <- as.character(meta[[col]])
    for (cand in list(names_in, .ba_sample_name(names_in))) {
      m <- match(cand, v); ok <- is.na(ids) & !is.na(m)
      ids[ok] <- meta$sample_id[m[ok]]
      m <- match(cand, .ba_sample_name(v)); ok <- is.na(ids) & !is.na(m)
      ids[ok] <- meta$sample_id[m[ok]]
    }
  }
  # last pass: a sample_id that is the only prefix of the name before _ . or - (S01_subset, S01.R1)
  for (i in which(is.na(ids))) {
    nm <- .ba_sample_name(names_in[i])
    hit <- meta$sample_id[startsWith(nm, meta$sample_id) & substr(nm, nchar(meta$sample_id) + 1, nchar(meta$sample_id) + 1) %in% c("_", ".", "-")]
    if (length(hit) == 1L) ids[i] <- hit
  }
  ids
}

.ba_write_matrix_tsv <- function(m, path, id = "feature_id") {
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  df <- data.frame(rownames(m), m, check.names = FALSE, stringsAsFactors = FALSE)
  names(df)[1] <- id
  utils::write.table(df, path, sep = "\t", row.names = FALSE, quote = FALSE)
  path
}

# Build the manifest entry from what a reader returned, renaming samples to sample_id.
.ba_entry_from_read <- function(x, g, row, dir, meta, assay, taken, id = NULL) {
  rules <- .ba_format_rules()
  warn <- x$warnings %||% character()
  dt <- switch(x$gives,
    counts = rules$assay_data_type[[assay %||% ""]] %||% "feature_counts",
    bisulfite = "bisulfite_counts", genotypes = "genotypes", beta = "methylation_beta", ip_input = "ip_input_counts",
    peak_calls = stop("Peak calls (", paste(basename(g$file), collapse = ", "), ") have no read counts. Write consensus peaks with bioana_consensus_peaks(), count them with featureCounts -F SAF, and give that table instead.", call. = FALSE),
    stop("Reader gave unknown kind: ", x$gives, call. = FALSE))
  if (identical(dt, "ip_input_counts") && is.null(x$ip)) {
    x <- .ba_split_ip_input(x$counts)
    if (is.null(x)) stop("Assay ", assay, " needs IP and input columns named <sample>_IP and <sample>_input.", call. = FALSE)
  }
  base <- switch(dt, feature_counts = "counts", peak_counts_long = "peaks", bisulfite_counts = "methylation",
                 genotypes = "genotypes", methylation_beta = "methylation_array", ip_input_counts = "ip_input")
  if (is.null(id)) {
    id <- base; k <- 2L
    while (id %in% names(taken)) { id <- paste0(base, k); k <- k + 1L }
  }
  ddir <- file.path("data", id)
  rename <- function(m) {
    s <- .ba_match_samples(colnames(m), meta)
    if (all(is.na(s))) stop(sprintf("None of the sample names in %s (%s ...) match sample_metadata_template.tsv.",
                                    paste(basename(g$file), collapse = ", "), paste(utils::head(colnames(m), 3), collapse = ", ")), call. = FALSE)
    if (anyNA(s)) warn <<- c(warn, sprintf("%d samples not in sample_metadata left out: %s", sum(is.na(s)), paste(utils::head(colnames(m)[is.na(s)], 4), collapse = ", ")))
    m <- m[, !is.na(s), drop = FALSE]; colnames(m) <- s[!is.na(s)]
    miss <- setdiff(meta$sample_id, colnames(m))
    if (length(miss)) warn <<- c(warn, sprintf("%d samples in sample_metadata have no data: %s", length(miss), paste(utils::head(miss, 4), collapse = ", ")))
    m
  }
  files <- switch(dt,
    feature_counts = {
      m <- rename(x$counts)
      list(counts = .ba_rel(.ba_write_matrix_tsv(m, file.path(dir, ddir, "counts.tsv")), dir))
    },
    peak_counts_long = {
      m <- rename(x$counts)
      long <- data.frame(peak_id = rep(rownames(m), ncol(m)), library_id = rep(colnames(m), each = nrow(m)), count = as.vector(m))
      dir.create(file.path(dir, ddir), recursive = TRUE, showWarnings = FALSE)
      utils::write.table(long, file.path(dir, ddir, "peak_counts.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
      utils::write.table(data.frame(library_id = colnames(m), sample_id = colnames(m)), file.path(dir, ddir, "libraries.tsv"),
                         sep = "\t", row.names = FALSE, quote = FALSE)
      list(counts = file.path(ddir, "peak_counts.tsv"), libraries = file.path(ddir, "libraries.tsv"))
    },
    bisulfite_counts = {
      me <- rename(x$methylated); cv <- rename(x$coverage)
      dir.create(file.path(dir, ddir), recursive = TRUE, showWarnings = FALSE)
      utils::write.table(x$sites, file.path(dir, ddir, "sites.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
      list(sites = file.path(ddir, "sites.tsv"),
           methylated = .ba_rel(.ba_write_matrix_tsv(me, file.path(dir, ddir, "methylated.tsv")), dir),
           coverage = .ba_rel(.ba_write_matrix_tsv(cv, file.path(dir, ddir, "coverage.tsv")), dir))
    },
    methylation_beta = list(beta = .ba_rel(.ba_write_matrix_tsv(rename(x$beta), file.path(dir, ddir, "beta.tsv")), dir)),
    ip_input_counts = list(ip = .ba_rel(.ba_write_matrix_tsv(rename(x$ip), file.path(dir, ddir, "ip.tsv")), dir),
                           input = .ba_rel(.ba_write_matrix_tsv(rename(x$input), file.path(dir, ddir, "input.tsv")), dir)),
    genotypes = {
      dir.create(file.path(dir, ddir), recursive = TRUE, showWarnings = FALSE)
      vcf_src <- x$vcf %||% g$file
      to <- file.path(dir, ddir, x$vcf_name %||% basename(vcf_src))
      if (normalizePath(vcf_src) != normalizePath(to, mustWork = FALSE)) file.copy(vcf_src, to, overwrite = TRUE)
      list(vcf = file.path(ddir, basename(to)))
    })
  samples <- switch(dt, genotypes = x$samples, bisulfite_counts = colnames(x$coverage), ip_input_counts = colnames(x$ip),
                    methylation_beta = colnames(x$beta), colnames(x$counts))
  list(id = id, data_type = dt, assay = assay %||% x$gives, files = files,
       source = list(format = row$format, tool = row$tool, reader = x$reader, files = normalizePath(g$file),
                     samples = samples, warnings = warn))
}

.ba_rel <- function(path, dir) sub(paste0("^", normalizePath(dir), "/?"), "", normalizePath(path))

.ba_split_ip_input <- function(m) {
  r <- .ba_format_rules()
  ipc <- grep(r$ip_suffix, colnames(m), perl = TRUE); inc <- grep(r$input_suffix, colnames(m), perl = TRUE)
  if (!length(ipc) || !length(inc)) return(NULL)
  s_ip <- sub(r$ip_suffix, "", colnames(m)[ipc], perl = TRUE); s_in <- sub(r$input_suffix, "", colnames(m)[inc], perl = TRUE)
  s <- intersect(s_ip, s_in)
  if (!length(s)) return(NULL)
  ip <- m[, ipc[match(s, s_ip)], drop = FALSE]; inp <- m[, inc[match(s, s_in)], drop = FALSE]
  colnames(ip) <- colnames(inp) <- s
  list(gives = "ip_input", ip = ip, input = inp, reader = "IP and input columns split by suffix")
}

# ---- readers: each returns list(gives, <matrices>, reader, warnings) ---------

.ba_read_featurecounts <- function(files, ...) {
  x <- utils::read.delim(files, comment.char = "#", check.names = FALSE, stringsAsFactors = FALSE)
  m <- as.matrix(x[, -(1:6), drop = FALSE]); rownames(m) <- x$Geneid
  colnames(m) <- .ba_sample_name(colnames(m))
  list(gives = "counts", counts = m, reader = "featureCounts table: Geneid + one column per BAM")
}

.ba_read_star <- function(files, strandedness = "auto", ...) {
  tabs <- lapply(files, function(f) utils::read.delim(f, header = FALSE, stringsAsFactors = FALSE))
  tabs <- lapply(tabs, function(t) t[!grepl("^N_", t$V1), ])
  s3 <- sum(vapply(tabs, function(t) sum(t$V3), 0)); s4 <- sum(vapply(tabs, function(t) sum(t$V4), 0))
  col <- switch(strandedness, unstranded = 2L, forward = 3L, reverse = 4L,
                auto = if (s3 > 5 * s4) 3L else if (s4 > 5 * s3) 4L else 2L)
  m <- do.call(cbind, lapply(tabs, function(t) t[[col]])); rownames(m) <- tabs[[1]]$V1
  colnames(m) <- vapply(files, .ba_sample_from_path, "", format = "star")
  list(gives = "counts", counts = m, reader = sprintf("STAR ReadsPerGene column %d (%s)", col, c("", "unstranded", "forward", "reverse")[col]),
       warnings = if (strandedness == "auto") sprintf("STAR strandedness detected as %s", c("", "unstranded", "forward", "reverse")[col]) else character())
}

# Salmon / kallisto through tximport (gene counts with tx2gene), the route DESeq2's vignette describes.
.ba_read_tximport <- function(files, tx2gene = NULL, ...) {
  type <- if (all(grepl("quant\\.(genes\\.)?sf", basename(files)))) "salmon" else "kallisto"
  names(files) <- vapply(files, .ba_sample_from_path, "", format = type)
  if (is.character(tx2gene) && length(tx2gene) == 1L) tx2gene <- utils::read.delim(tx2gene, stringsAsFactors = FALSE)
  warn <- character()
  if (requireNamespace("tximport", quietly = TRUE)) {
    txi <- if (is.null(tx2gene)) tximport::tximport(files, type = type, txOut = TRUE)
           else tximport::tximport(files, type = type, tx2gene = tx2gene, countsFromAbundance = "lengthScaledTPM")
    m <- round(txi$counts)
    rd <- sprintf("tximport %s (type = '%s'%s)", utils::packageVersion("tximport"), type,
                  if (is.null(tx2gene)) ", txOut = TRUE: transcript level" else ", tx2gene, countsFromAbundance = 'lengthScaledTPM'")
  } else {
    col <- if (type == "salmon") "NumReads" else "est_counts"; idc <- if (type == "salmon") "Name" else "target_id"
    tabs <- lapply(files, function(f) utils::read.delim(f, stringsAsFactors = FALSE))
    m <- round(do.call(cbind, lapply(tabs, `[[`, col))); rownames(m) <- tabs[[1]][[idc]]
    if (!is.null(tx2gene)) m <- rowsum(m, tx2gene[[2]][match(rownames(m), tx2gene[[1]])])
    rd <- sprintf("%s column read directly (tximport not installed)", col)
    warn <- "tximport not installed: estimated counts read directly (install tximport for the citable route)"
  }
  colnames(m) <- names(files)
  if (is.null(tx2gene)) warn <- c(warn, "no tx2gene: analysed at transcript level")
  list(gives = "counts", counts = m, reader = rd, warnings = warn)
}

.ba_read_bismark <- function(files, ...) {
  rows <- lapply(files, function(f) {
    x <- utils::read.delim(f, header = FALSE, stringsAsFactors = FALSE)
    if (ncol(x) >= 7 && grepl("^(CG|CHG|CHH)$", x$V6[1])) {
      data.frame(chr = x$V1, pos = x$V2, context = x$V6, m = x$V4, u = x$V5)
    } else {
      data.frame(chr = x$V1, pos = x$V2, context = NA_character_, m = x$V5, u = x$V6)
    }
  })
  samples <- vapply(files, .ba_sample_from_path, "", format = "bismark")
  key <- unique(do.call(rbind, lapply(rows, function(r) r[, c("chr", "pos", "context")])))
  sid <- function(chr, pos, ctx) ifelse(is.na(ctx), paste(chr, pos, sep = ":"), paste(chr, pos, ctx, sep = ":"))
  key$feature_id <- sid(key$chr, key$pos, key$context)
  me <- cv <- matrix(0, nrow(key), length(files), dimnames = list(key$feature_id, samples))
  for (j in seq_along(rows)) {
    r <- rows[[j]]; k <- match(sid(r$chr, r$pos, r$context), key$feature_id)
    me[k, j] <- r$m; cv[k, j] <- r$m + r$u
  }
  sites <- data.frame(feature_id = key$feature_id, chr = key$chr, pos = key$pos)
  if (!all(is.na(key$context))) sites$context <- key$context
  list(gives = "bisulfite", sites = sites, methylated = me, coverage = cv,
       reader = "Bismark methylated / unmethylated counts per cytosine, one file per sample")
}

.ba_read_peaks <- function(files, ...) {
  p <- do.call(rbind, lapply(files, function(f) {
    x <- utils::read.delim(f, header = FALSE, stringsAsFactors = FALSE)
    data.frame(chr = x$V1, start = x$V2 + 1L, end = x$V3, sample = .ba_sample_name(sub("_peaks$", "", sub("\\.(narrowPeak|broadPeak)(\\.gz)?$", "", basename(f)))))
  }))
  list(gives = "peak_calls", peaks = p, reader = "MACS peak calls (BED, 0-based start)")
}

.ba_read_vcf_entry <- function(files, ...) {
  hdr <- .ba_head_lines(files, 5000)
  cols <- strsplit(grep("^#CHROM", hdr, value = TRUE)[1], "\t", fixed = TRUE)[[1]]
  list(gives = "genotypes", samples = cols[-(1:9)], reader = "VCF (sample columns = line or sample ids)")
}

.ba_read_plain_matrix <- function(files, ...) {
  sep <- if (grepl("\\.csv(\\.gz)?$", files)) "," else "\t"
  x <- utils::read.table(files, sep = sep, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE, comment.char = "", quote = "\"")
  m <- as.matrix(x[, -1, drop = FALSE]); rownames(m) <- as.character(x[[1]])
  if (!is.numeric(m)) stop("Table ", files, ": every column after the first must be numeric.", call. = FALSE)
  v <- m[is.finite(m)]
  if (all(v >= 0 & v == round(v))) {
    sp <- .ba_split_ip_input(m)
    if (!is.null(sp)) return(sp)
    return(list(gives = "counts", counts = m, reader = "feature x sample count table"))
  }
  if (all(v >= 0 & v <= 1)) return(list(gives = "beta", beta = m, reader = "feature x sample table of methylation beta values"))
  stop("Table ", files, ": values are neither counts (non-negative integers) nor beta values (0-1). Give raw counts, not normalised values: the engines need them.", call. = FALSE)
}

# ---- BAM: Rsubread featureCounts --------------------------------------------

# Library table for BAM files: file, sample_id, library. Without one, each file is a sample (library from the assay).
.ba_bam_libraries <- function(files, libraries, assay) {
  if (is.null(libraries)) {
    lib <- if (!is.null(assay) && grepl("chip|cut|atac|dnase", assay)) "chip" else "rna"
    return(data.frame(file = files, sample_id = .ba_sample_name(files), library = lib, stringsAsFactors = FALSE))
  }
  libraries <- as.data.frame(libraries, stringsAsFactors = FALSE)
  m <- match(basename(files), basename(libraries$file))
  if (anyNA(m)) stop("BAM files missing from libraries$file: ", paste(basename(files)[is.na(m)], collapse = ", "), call. = FALSE)
  data.frame(file = files, sample_id = libraries$sample_id[m], library = tolower(libraries$library[m]), stringsAsFactors = FALSE)
}

# SAF from the annotation: GTF genes (exons per gene_id), GTF exon windows, a SAF file, or peak files.
.ba_saf <- function(annotation, mode, window_bp = 200L) {
  if (is.null(annotation)) stop("BAM files need annotation = a GTF, a SAF file, or peak files", call. = FALSE)
  if (length(annotation) == 1L && grepl("\\.saf(\\.gz)?$", annotation, ignore.case = TRUE)) {
    saf <- utils::read.delim(annotation, stringsAsFactors = FALSE)
    return(list(saf = saf, what = sprintf("SAF %s (%d regions)", basename(annotation), nrow(saf))))
  }
  if (length(annotation) == 1L && grepl("\\.gtf(\\.gz)?$", annotation)) {
    g <- utils::read.delim(annotation, header = FALSE, comment.char = "#", stringsAsFactors = FALSE, quote = "")
    g <- g[g$V3 == "exon", ]
    gid <- sub('.*gene_id "([^"]+)".*', "\\1", g$V9)
    if (mode != "windows") {
      return(list(saf = data.frame(GeneID = gid, Chr = g$V1, Start = g$V4, End = g$V5, Strand = g$V7),
                  what = sprintf("GTF %s exons by gene_id (%d genes)", basename(annotation), length(unique(gid)))))
    }
    # exon windows: merge each gene's exons, tile in window_bp steps
    w <- do.call(rbind, lapply(split(data.frame(chr = g$V1, s = g$V4, e = g$V5, strand = g$V7, gid = gid), gid), function(x) {
      x <- x[order(x$s), ]; grp <- cumsum(c(TRUE, x$s[-1] > cummax(x$e)[-nrow(x)]))
      ex <- data.frame(s = tapply(x$s, grp, min), e = tapply(x$e, grp, max))
      do.call(rbind, lapply(seq_len(nrow(ex)), function(i) {
        st <- seq(ex$s[i], ex$e[i], by = window_bp)
        data.frame(Chr = x$chr[1], Start = st, End = pmin(st + window_bp - 1L, ex$e[i]), Strand = x$strand[1], gene = x$gid[1])
      }))
    }))
    w$GeneID <- sprintf("%s:%s:%d-%d", w$gene, w$Chr, w$Start, w$End)
    return(list(saf = w[, c("GeneID", "Chr", "Start", "End", "Strand")],
                what = sprintf("GTF %s exon windows of %d bp (%d windows in %d genes)", basename(annotation), window_bp, nrow(w), length(unique(w$gene)))))
  }
  # peak files: consensus peaks present in >= 2 files (or the one file)
  tmp <- tempfile(fileext = ".saf")
  bioana_consensus_peaks(annotation, tmp, min_samples = min(2L, length(annotation)))
  saf <- utils::read.delim(tmp, stringsAsFactors = FALSE)
  list(saf = saf, what = sprintf("consensus of %d peak files, regions in >= %d files (%d regions)", length(annotation), min(2L, length(annotation)), nrow(saf)))
}

.ba_read_bam <- function(files, annotation = NULL, libraries = NULL, assay = NULL, window_bp = 200L, strandedness = "auto", ...) {
  if (!requireNamespace("Rsubread", quietly = TRUE)) stop("Counting reads in BAM files needs Rsubread (Bioconductor).", call. = FALSE)
  lib <- .ba_bam_libraries(files, libraries, assay)
  ipin <- any(lib$library %in% c("ip", "input"))
  sa <- .ba_saf(annotation, if (ipin) "windows" else "genes", window_bp)
  paired <- .ba_bam_paired(files[1])
  strand <- switch(strandedness, forward = 1L, reverse = 2L, 0L)
  # ChIP-type single-end reads are extended to the fragment length before counting, as DiffBind
  # (fragmentSize = 125) and csaw do; rules: file_formats.json chip_fragment_bp
  chip <- any(lib$library %in% c("chip", "control")) && !paired
  ext <- 0L
  if (chip) ext <- max(0L, as.integer(.ba_format_rules()$chip_fragment_bp %||% 125L) - .ba_bam_read_length(files[1]))
  utils::capture.output(fc <- Rsubread::featureCounts(files, annot.ext = sa$saf, isGTFAnnotationFile = FALSE, isPairedEnd = paired,
                                                     countReadPairs = paired, strandSpecific = strand, allowMultiOverlap = ipin || chip,
                                                     readExtension3 = ext,
                                                     nthreads = 2L, verbose = FALSE))
  m <- fc$counts
  colnames(m) <- lib$sample_id[match(colnames(m), make.names(basename(files)))]
  if (anyNA(colnames(m))) colnames(m) <- lib$sample_id
  st <- fc$stat; assigned <- colSums(st[st$Status == "Assigned", -1, drop = FALSE])
  reader <- sprintf("Rsubread::featureCounts %s (%s; isPairedEnd = %s, strandSpecific = %d%s%s); %s of reads assigned",
                    utils::packageVersion("Rsubread"), sa$what, paired, strand, if (ipin || chip) ", allowMultiOverlap = TRUE" else "",
                    if (ext > 0) sprintf(", readExtension3 = %d (reads extended to %d bp fragments)", ext, ext + .ba_bam_read_length(files[1])) else "",
                    paste0(round(100 * sum(assigned) / sum(colSums(st[, -1, drop = FALSE]))), "%"))
  if (ipin) {
    ip <- m[, lib$library == "ip", drop = FALSE]; inp <- m[, lib$library == "input", drop = FALSE]
    s <- intersect(colnames(ip), colnames(inp))
    warn <- if (length(setdiff(union(colnames(ip), colnames(inp)), s))) sprintf("IP or input without its pair left out: %s", paste(setdiff(union(colnames(ip), colnames(inp)), s), collapse = ", ")) else character()
    return(list(gives = "ip_input", ip = ip[, s, drop = FALSE], input = inp[, s, drop = FALSE], reader = reader, warnings = warn))
  }
  keep <- !lib$library %in% c("control", "input")
  list(gives = "counts", counts = m[, keep, drop = FALSE], reader = reader,
       warnings = if (any(!keep)) sprintf("%d control / input libraries not used in the test", sum(!keep)) else character())
}

# ---- HapMap -> VCF -----------------------------------------------------------

.ba_read_hapmap <- function(files, ...) {
  x <- utils::read.delim(files, check.names = FALSE, stringsAsFactors = FALSE, colClasses = "character")
  al <- strsplit(x$alleles, "/", fixed = TRUE)
  ref <- vapply(al, `[`, "", 1); alt <- vapply(al, function(a) if (length(a) > 1) a[2] else ".", "")
  calls <- as.matrix(x[, -(1:11), drop = FALSE])
  a1 <- substr(calls, 1, 1); a2 <- substr(calls, nchar(calls), nchar(calls))
  code <- function(a) ifelse(a == ref, "0", ifelse(a == alt, "1", "."))
  gt <- matrix(paste0(code(a1), "/", code(a2)), nrow(calls))
  gt[gt == "1/0"] <- "0/1"; gt[grepl("\\.", gt)] <- "./."
  keep <- alt != "." & grepl("^[ACGT]$", ref) & grepl("^[ACGT]$", alt)
  out <- tempfile(fileext = ".vcf")
  con <- file(out, "w")
  writeLines(c("##fileformat=VCFv4.2", sprintf("##source=bioana HapMap conversion of %s", basename(files)),
               paste(c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", colnames(calls)), collapse = "\t")), con)
  body <- paste(x$chrom[keep], x$pos[keep], x[[1]][keep], ref[keep], alt[keep], ".", "PASS", ".", "GT",
                apply(gt[keep, , drop = FALSE], 1, paste, collapse = "\t"), sep = "\t")
  writeLines(body, con); close(con)
  list(gives = "genotypes", samples = colnames(calls), vcf = out, vcf_name = "genotypes.vcf",
       reader = sprintf("HapMap -> VCF (%d of %d biallelic SNPs; two-letter calls)", sum(keep), nrow(x)),
       warnings = if (any(!keep)) sprintf("%d non-biallelic or indel markers dropped", sum(!keep)) else character())
}

# ---- DSS text: chr, pos, N, X per sample -------------------------------------

.ba_read_dss_txt <- function(files, ...) {
  rows <- lapply(files, function(f) {
    x <- utils::read.delim(f, stringsAsFactors = FALSE)
    data.frame(chr = x$chr, pos = x$pos, context = NA_character_, m = x$X, u = x$N - x$X)
  })
  samples <- vapply(files, function(f) .ba_sample_name(sub("\\.(txt|tsv)(\\.gz)?$", "", basename(f))), "")
  key <- unique(do.call(rbind, lapply(rows, function(r) r[, c("chr", "pos")])))
  key$feature_id <- paste(key$chr, key$pos, sep = ":")
  me <- cv <- matrix(0, nrow(key), length(files), dimnames = list(key$feature_id, samples))
  for (j in seq_along(rows)) { r <- rows[[j]]; k <- match(paste(r$chr, r$pos, sep = ":"), key$feature_id); me[k, j] <- r$m; cv[k, j] <- r$m + r$u }
  list(gives = "bisulfite", sites = data.frame(feature_id = key$feature_id, chr = key$chr, pos = key$pos),
       methylated = me, coverage = cv, reader = "DSS/BSmooth count text: chr, pos, N (reads), X (methylated), one file per sample")
}

# Paired-end if any of the first 1000 reads has the paired flag (no index needed).
.ba_bam_paired <- function(f) {
  if (!requireNamespace("Rsamtools", quietly = TRUE)) return(FALSE)
  bf <- Rsamtools::BamFile(f, index = character(0), yieldSize = 1000L)
  open(bf); on.exit(close(bf))
  fl <- Rsamtools::scanBam(bf, param = Rsamtools::ScanBamParam(what = "flag"))[[1]]$flag
  any(bitwAnd(fl, 1L) > 0)
}

# Median read length of the first 1000 reads.
.ba_bam_read_length <- function(f) {
  bf <- Rsamtools::BamFile(f, index = character(0), yieldSize = 1000L)
  open(bf); on.exit(close(bf))
  as.integer(stats::median(Rsamtools::scanBam(bf, param = Rsamtools::ScanBamParam(what = "qwidth"))[[1]]$qwidth, na.rm = TRUE))
}
