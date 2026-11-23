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
    head <- .ba_head_lines(f)
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
#' @return Invisibly, the manifest entries; prints one line per data set.
#' @export
bioana_manifest_from_files <- function(dir, files, assay = NULL, tx2gene = NULL, strandedness = "auto", append = TRUE) {
  meta <- utils::read.delim(file.path(dir, "sample_metadata_template.tsv"), stringsAsFactors = FALSE, colClasses = "character")
  spec <- jsonlite::read_json(file.path(dir, "study_spec.json"))
  assay <- assay %||% .ba_spec_assay(spec)
  det <- bioana_detect_format(files)
  if (anyNA(det$format)) {
    stop("Format not recognised (see inst/rules/file_formats.tsv): ", paste(det$file[is.na(det$format)], collapse = ", "), call. = FALSE)
  }
  groups <- split(det, ifelse(det$per_sample, det$format, paste0(det$format, "::", det$file)))
  entries <- list()
  for (g in groups) {
    ff <- .ba_formats(); row <- ff[ff$format == g$format[1], ]
    reader <- get(row$reader, envir = asNamespace("bioana"))
    x <- reader(g$file, tx2gene = tx2gene, strandedness = strandedness)
    e <- .ba_entry_from_read(x, g, row, dir, meta, assay, entries)
    entries[[e$id]] <- e
  }
  old <- if (append && file.exists(file.path(dir, "data_manifest.json"))) jsonlite::read_json(file.path(dir, "data_manifest.json"))$data else list()
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
.ba_entry_from_read <- function(x, g, row, dir, meta, assay, taken) {
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
  id <- base; k <- 2L
  while (id %in% names(taken)) { id <- paste0(base, k); k <- k + 1L }
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
      to <- file.path(dir, ddir, basename(g$file))
      if (normalizePath(g$file) != normalizePath(to, mustWork = FALSE)) file.copy(g$file, to, overwrite = TRUE)
      list(vcf = file.path(ddir, basename(g$file)))
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
