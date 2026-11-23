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
  paths <- c(results = .ba_write_tsv(ft$table, file.path(outdir, "results.tsv")))
  truth <- .ba_entry_truth(entry, design$dir)
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

# feature_counts: files$counts = feature x sample count matrix (first column the feature id).
.ba_method_feature_counts <- function(entry, design, outdir, alpha = 0.05, args = list()) {
  m <- .ba_read_matrix(file.path(design$dir, entry$files$counts))
  lc <- .ba_logcpm(m, args$min_cpm %||% 1)
  ft <- .ba_feature_test(lc$y, design$meta, design$parsed, alpha)
  if (lc$dropped) ft$warnings <- c(ft$warnings, sprintf("%d low-count features not tested (CPM < %g in most samples)", lc$dropped, args$min_cpm %||% 1))
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
  lc <- .ba_logcpm(m, args$min_cpm %||% 1)
  ft <- .ba_feature_test(lc$y, design$meta, design$parsed, alpha)
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
  lc <- .ba_logcpm(m, args$min_cpm %||% 1)
  ft <- .ba_feature_test(lc$y, design$meta, design$parsed, alpha)
  feats <- data.frame(feature_id = rownames(m), chr = sub(":.*", "", rownames(m)),
                      start = as.integer(sub(".*:", "", rownames(m))), stringsAsFactors = FALSE)
  feats$end <- feats$start + bin - 1L
  .ba_finish_features(ft, entry, design, outdir, alpha, list(scale = sprintf("log2 CPM of reads in %d bp windows", bin), window_bp = bin),
                      features = feats, large = 1)
}

# bisulfite_counts: sites (feature_id, chr, pos[, context]) with methylated and coverage matrices.
# Sites are pooled into windows; each window's methylated share per sample on the logit scale (M-value).
.ba_method_bisulfite_counts <- function(entry, design, outdir, alpha = 0.05, args = list()) {
  win <- as.integer(args$window_bp %||% 1000); min_cov <- as.numeric(args$min_coverage %||% 10)
  sites <- utils::read.delim(file.path(design$dir, entry$files$sites), stringsAsFactors = FALSE)
  me <- .ba_read_matrix(file.path(design$dir, entry$files$methylated))
  cv <- .ba_read_matrix(file.path(design$dir, entry$files$coverage))
  sites <- sites[match(rownames(cv), sites$feature_id), , drop = FALSE]
  ctx <- if (!is.null(sites$context) && !isTRUE(args$pool_contexts)) sites$context else "all"
  w <- paste0(sites$chr, ":", (sites$pos %/% win) * win + 1L, if (length(ctx) > 1) paste0(":", ctx) else "")
  M <- rowsum(me[, colnames(cv), drop = FALSE], w); C <- rowsum(cv, w)
  Y <- log2((M + 0.5) / (C - M + 0.5))
  Y[C < min_cov] <- NA
  ft <- .ba_feature_test(Y, design$meta, design$parsed, alpha)
  n_low <- sum(!stats::complete.cases(Y))
  if (n_low) ft$warnings <- c(ft$warnings, sprintf("%d windows below %g reads in some sample not tested", n_low, min_cov))
  # report the change also as a difference in methylated share at the mean level
  prop <- C; prop[] <- M / pmax(C, 1)
  feats <- data.frame(feature_id = rownames(Y), chr = sub(":.*", "", rownames(Y)),
                      start = as.integer(sub("^[^:]*:([0-9]+).*$", "\\1", rownames(Y))), stringsAsFactors = FALSE)
  feats$end <- feats$start + win - 1L
  ft$table$mean_methylation <- round(rowMeans(prop)[ft$table$feature_id], 4)
  .ba_finish_features(ft, entry, design, outdir, alpha,
                      list(scale = sprintf("M-value (log2 methylated/unmethylated reads) per %d bp window%s", win, if (length(ctx) > 1) " and context" else ""),
                           window_bp = win, min_coverage = min_cov),
                      features = feats)
}

# methylation_beta: array beta values (probe x sample) tested as M-values.
.ba_method_methylation_beta <- function(entry, design, outdir, alpha = 0.05, args = list()) {
  b <- .ba_read_matrix(file.path(design$dir, entry$files$beta))
  b <- pmin(pmax(b, 0.001), 0.999)
  ft <- .ba_feature_test(log2(b / (1 - b)), design$meta, design$parsed, alpha)
  .ba_finish_features(ft, entry, design, outdir, alpha, list(scale = "M-value (log2 beta/(1-beta))"))
}

# ip_input_counts (MeRIP-seq, MeDIP, RIP): IP and input counts per window. Each library is scaled to
# its size; the tested value is log2(IP / input) per window, so expression changes (both libraries)
# cancel and only enrichment changes remain.
.ba_method_ip_input <- function(entry, design, outdir, alpha = 0.05, args = list()) {
  ip <- .ba_read_matrix(file.path(design$dir, entry$files$ip))
  inp <- .ba_read_matrix(file.path(design$dir, entry$files$input))
  ids <- intersect(colnames(ip), colnames(inp))
  ip <- ip[, ids, drop = FALSE]; inp <- inp[rownames(ip), ids, drop = FALSE]
  keep <- rowMeans(inp) >= (args$min_input %||% 5)
  # IP efficiency differs between units: scale IP by its median ratio to input over all windows kept
  sf_in <- colSums(inp) / mean(colSums(inp))
  r <- (ip[keep, , drop = FALSE] + 1) / sweep(inp[keep, , drop = FALSE] + 1, 2, sf_in, "/")
  sf_ip <- apply(r, 2, stats::median)
  Y <- log2(sweep(ip[keep, , drop = FALSE] + 1, 2, sf_ip, "/") / sweep(inp[keep, , drop = FALSE] + 1, 2, sf_in, "/"))
  ft <- .ba_feature_test(Y, design$meta, design$parsed, alpha)
  ft$table$mean_log2_enrichment <- round(rowMeans(Y)[ft$table$feature_id], 3)
  .ba_finish_features(ft, entry, design, outdir, alpha,
                      list(scale = "log2(IP / input) per window, IP scaled by its median ratio to input", windows_tested = sum(keep)))
}
