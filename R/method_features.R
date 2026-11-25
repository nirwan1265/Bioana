# Feature-level data types. Each turns its files into a features x samples
# matrix on a scale where a linear model fits (log-CPM, M-values, log
# IP/input), then calls .ba_feature_test() with the design's fixed terms.

# Counts per feature (genes, miRNAs, taxa, peaks): log2-CPM, features below min_cpm in most samples dropped.
.ba_logcpm <- function(m, min_cpm = 1) {
  lib <- colSums(m)
  cpm <- sweep(m, 2, lib, "/") * 1e6
  keep <- rowSums(cpm >= min_cpm) >= ceiling(ncol(m) / 2)
  list(y = log2(sweep(m[keep, , drop = FALSE] + 0.5, 2, lib + 1, "/") * 1e6), dropped = sum(!keep))
}

.ba_finish_features <- function(ft, entry, design, outdir, alpha, extra_summary = list(), features = NULL, large = NULL) {
  ft$summary$data_type <- entry$data_type
  ft$summary$assay <- entry$assay
  ft$summary <- c(ft$summary, extra_summary)
  eng <- ft$summary$engine
  if (!is.null(eng) && !identical(.ba_tsv("engines")$input[.ba_tsv("engines")$engine == eng$engine], "values")) {
    ft$summary$scale <- eng$estimate_scale
  }
  paths <- c(results = .ba_write_tsv(ft$table, file.path(outdir, "results.tsv")))
  truth <- .ba_entry_truth(entry, design$dir)
  if (!is.null(truth) && identical(ft$summary$goal$test, "lrt")) {
    ft$warnings <- c(ft$warnings, "truth scoring of joint tests (time_course, interaction) comes later: not scored")
    truth <- NULL
  }
  if (!is.null(truth)) {
    if (identical(truth$kind, "regions")) {
      if (is.null(features)) stop("region truth needs feature coordinates", call. = FALSE)
      truth$table <- .ba_regions_to_features(features, truth$table)
    }
    ft$summary$recovery <- .ba_feature_recovery(ft$table, truth$table, alpha, large = truth$large %||% large)
  }
  ft$summary$warnings <- ft$warnings
  list(summary = ft$summary, paths = paths, results = ft$table)
}

# Counts per feature: the engine chain (DESeq2 / edgeR / limma-voom / dream, native log-CPM fallback).
.ba_count_inputs <- function(m, args, what = "log2 CPM") {
  v <- m[is.finite(m)]
  if (any(v < 0) || any(abs(v - round(v)) > 1e-6)) {
    stop("Count check failed: the count matrix has negative or non-whole values (e.g. ", v[v < 0 | abs(v - round(v)) > 1e-6][1],
         "). The count engines need raw read counts, not normalised values (CPM, TPM, FPKM).", call. = FALSE)
  }
  lc <- .ba_logcpm(m, args$min_cpm %||% 1)
  list(inputs = list(counts = m[rownames(lc$y), , drop = FALSE], values = lc$y, values_scale = what), dropped = lc$dropped)
}

# feature_counts: files$counts = feature x sample count matrix (first column the feature id).
.ba_method_feature_counts <- function(entry, design, outdir, alpha = 0.05, args = list()) {
  m <- .ba_read_matrix(file.path(design$dir, entry$files$counts))
  ci <- .ba_count_inputs(m, args)
  ft <- .ba_feature_engines(ci$inputs, entry, design, alpha, args)
  if (ci$dropped) ft$warnings <- c(ft$warnings, sprintf("%d low-count features not tested (CPM < %g in most samples)", ci$dropped, args$min_cpm %||% 1))
  .ba_finish_features(ft, entry, design, outdir, alpha, list(scale = "log2 CPM"), large = 1)
}

# peak_counts_long: files$counts = long table (peak_id, library_id, count[, assay]); files$libraries maps
# library_id to the unit (column `sample` or `sample_id`). Input/control libraries are left out.
.ba_method_peak_counts_long <- function(entry, design, outdir, alpha = 0.05, args = list()) {
  x <- utils::read.delim(file.path(design$dir, entry$files$counts), stringsAsFactors = FALSE)
  if (!is.null(x$assay)) x <- x[!tolower(x$assay) %in% c("input", "control"), , drop = FALSE]
  lib <- utils::read.delim(file.path(design$dir, entry$files$libraries), stringsAsFactors = FALSE)
  scol <- intersect(c("sample_id", "sample"), names(lib))[1]
  x$sample_id <- lib[[scol]][match(x$library_id, lib$library_id)]
  m <- tapply(x$count, list(x$peak_id, x$sample_id), sum)
  m[is.na(m)] <- 0
  m <- unclass(m); storage.mode(m) <- "double"
  ci <- .ba_count_inputs(m, args, "log2 CPM of IP reads per peak")
  ft <- .ba_feature_engines(ci$inputs, entry, design, alpha, args)
  .ba_finish_features(ft, entry, design, outdir, alpha, list(scale = "log2 CPM of IP reads per peak"), large = 1)
}

# read_positions: reads (seqname, read_position, library_id) counted in fixed windows; no peak calling.
.ba_method_read_positions <- function(entry, design, outdir, alpha = 0.05, args = list()) {
  bin <- as.integer(args$bin_bp %||% 500)
  x <- utils::read.delim(file.path(design$dir, entry$files$reads), stringsAsFactors = FALSE)
  if (!is.null(x$is_duplicate)) x <- x[!as.logical(x$is_duplicate), , drop = FALSE]
  lib <- utils::read.delim(file.path(design$dir, entry$files$libraries), stringsAsFactors = FALSE)
  scol <- intersect(c("sample_id", "sample"), names(lib))[1]
  x$sample_id <- lib[[scol]][match(x$library_id, lib$library_id)]
  x$bin <- paste0(x$seqname, ":", (x$read_position %/% bin) * bin + 1L)
  m <- unclass(table(x$bin, x$sample_id)); storage.mode(m) <- "double"
  m <- m[rowSums(m) >= (args$min_reads %||% 10), , drop = FALSE]
  sc <- sprintf("log2 CPM of reads in %d bp windows", bin)
  ci <- .ba_count_inputs(m, args, sc)
  ft <- .ba_feature_engines(ci$inputs, entry, design, alpha, args)
  feats <- data.frame(feature_id = rownames(m), chr = sub(":.*", "", rownames(m)),
                      start = as.integer(sub(".*:", "", rownames(m))), stringsAsFactors = FALSE)
  feats$end <- feats$start + bin - 1L
  .ba_finish_features(ft, entry, design, outdir, alpha, list(scale = sc, window_bp = bin), features = feats, large = 1)
}

# bisulfite_counts: sites (feature_id, chr, pos[, context]) with methylated and coverage matrices.
# Sites are pooled into windows (per context). DSS fits the pooled reads; limma and the native
# fallback test each window's M-value (log2 methylated/unmethylated).
.ba_method_bisulfite_counts <- function(entry, design, outdir, alpha = 0.05, args = list()) {
  win <- as.integer(args$window_bp %||% 1000); min_cov <- as.numeric(args$min_coverage %||% 10)
  sites <- utils::read.delim(file.path(design$dir, entry$files$sites), stringsAsFactors = FALSE)
  me <- .ba_read_matrix(file.path(design$dir, entry$files$methylated))
  cv <- .ba_read_matrix(file.path(design$dir, entry$files$coverage))
  sites <- sites[match(rownames(cv), sites$feature_id), , drop = FALSE]
  ctx <- if (!is.null(sites$context) && !isTRUE(args$pool_contexts)) sites$context else "all"
  # positions are 1-based (Bismark, DSS): window k covers (k-1)*win+1 .. k*win
  w <- paste0(sites$chr, ":", ((sites$pos - 1L) %/% win) * win + 1L, if (length(ctx) > 1) paste0(":", ctx) else "")
  M <- rowsum(me[, colnames(cv), drop = FALSE], w); C <- rowsum(cv, w)
  Y <- log2((M + 0.5) / (C - M + 0.5))
  Y[C < min_cov] <- NA
  keep <- stats::complete.cases(Y)
  feats <- data.frame(feature_id = rownames(Y), chr = sub(":.*", "", rownames(Y)),
                      start = as.integer(sub("^[^:]*:([0-9]+).*$", "\\1", rownames(Y))), stringsAsFactors = FALSE)
  feats$end <- feats$start + win - 1L
  if (length(ctx) > 1) feats$context <- sub("^[^:]*:[0-9]+:", "", rownames(Y))
  sc <- sprintf("M-value (log2 methylated/unmethylated reads) per %d bp window%s", win, if (length(ctx) > 1) " and context" else "")
  inputs <- list(values = Y[keep, , drop = FALSE], values_scale = sc, methylated = M[keep, , drop = FALSE],
                 coverage = C[keep, , drop = FALSE], features = feats)
  ft <- .ba_feature_engines(inputs, entry, design, alpha, args)
  if (sum(!keep)) ft$warnings <- c(ft$warnings, sprintf("%d windows below %g reads in some sample not tested", sum(!keep), min_cov))
  prop <- C; prop[] <- M / pmax(C, 1)
  ft$table$mean_methylation <- round(rowMeans(prop)[ft$table$feature_id], 4)
  .ba_finish_features(ft, entry, design, outdir, alpha, list(scale = sc, window_bp = win, min_coverage = min_cov), features = feats)
}

# methylation_beta: array beta values (probe x sample) tested as M-values.
.ba_method_methylation_beta <- function(entry, design, outdir, alpha = 0.05, args = list()) {
  b <- .ba_read_matrix(file.path(design$dir, entry$files$beta))
  b <- pmin(pmax(b, 0.001), 0.999)
  sc <- "M-value (log2 beta/(1-beta))"
  ft <- .ba_feature_engines(list(values = log2(b / (1 - b)), values_scale = sc), entry, design, alpha, args)
  .ba_finish_features(ft, entry, design, outdir, alpha, list(scale = sc))
}

# ip_input_counts (MeRIP-seq, MeDIP, RIP): IP and input counts per window. edgeR fits both libraries
# with an IP x treatment interaction; limma and the native fallback test log2(IP / input) per window,
# IP scaled by its median ratio to input, so expression changes (both libraries) cancel.
.ba_method_ip_input <- function(entry, design, outdir, alpha = 0.05, args = list()) {
  ip <- .ba_read_matrix(file.path(design$dir, entry$files$ip))
  inp <- .ba_read_matrix(file.path(design$dir, entry$files$input))
  ids <- intersect(colnames(ip), colnames(inp))
  ip <- ip[, ids, drop = FALSE]; inp <- inp[rownames(ip), ids, drop = FALSE]
  keep <- rowMeans(inp) >= (args$min_input %||% 5)
  sf_in <- colSums(inp) / mean(colSums(inp))
  r <- (ip[keep, , drop = FALSE] + 1) / sweep(inp[keep, , drop = FALSE] + 1, 2, sf_in, "/")
  sf_ip <- apply(r, 2, stats::median)
  Y <- log2(sweep(ip[keep, , drop = FALSE] + 1, 2, sf_ip, "/") / sweep(inp[keep, , drop = FALSE] + 1, 2, sf_in, "/"))
  sc <- "log2(IP / input) per window, IP scaled by its median ratio to input"
  inputs <- list(values = Y, values_scale = sc, ip = ip[keep, , drop = FALSE], input = inp[keep, , drop = FALSE])
  ft <- .ba_feature_engines(inputs, entry, design, alpha, args)
  ft$table$mean_log2_enrichment <- round(rowMeans(Y)[ft$table$feature_id], 3)
  .ba_finish_features(ft, entry, design, outdir, alpha, list(scale = sc, windows_tested = sum(keep)))
}
