# ChIP / ATAC peak analysis on read positions (moved from Simitall): windowed peak calling, differential peaks, scoring against truth peaks, results figure.

#' Call peaks with a simple window-based Poisson test
#'
#' Counts reads in fixed windows, tests each window against a background rate
#' (scaled control reads when a control assay is given, otherwise a trimmed
#' genome-wide mean), applies Benjamini-Hochberg correction, and merges
#' adjacent significant windows. This is a transparent baseline caller, not a
#' replacement for MACS.
#'
#' @param reads Read-position table or path (as written by the simulators).
#' @param genome_fa Reference FASTA (for chromosome lengths).
#' @param signal_assay Assay label of signal reads (`"chip"` or `"atac"`).
#' @param control_assay Optional control assay label (e.g. `"input"`).
#' @param window Window width in bp.
#' @param fdr BH-adjusted threshold.
#' @param min_count Minimum reads in a window.
#' @return A data frame of called peaks.
#' @export
call_peaks_windowed <- function(reads, genome_fa, signal_assay = "chip", control_assay = NULL,
                                window = 200L, fdr = 0.05, min_count = 5L) {
  reads <- .bioana_peak_reads(reads)
  if ("is_duplicate" %in% names(reads)) reads <- reads[!as.logical(reads$is_duplicate), , drop = FALSE]
  genome_lengths <- nchar(.bioana_chip_read_fasta(genome_fa))
  sig <- reads[reads$assay == signal_assay, , drop = FALSE]
  ctl <- if (!is.null(control_assay)) reads[reads$assay == control_assay, , drop = FALSE] else NULL
  if (!nrow(sig)) stop("No reads for signal assay '", signal_assay, "'")
  bins <- lapply(names(genome_lengths), function(chr) {
    nb <- ceiling(genome_lengths[[chr]] / window)
    s <- tabulate(floor((sig$read_position[sig$seqname == chr] - 1) / window) + 1L, nb)
    cc <- if (!is.null(ctl)) tabulate(floor((ctl$read_position[ctl$seqname == chr] - 1) / window) + 1L, nb) else NULL
    data.frame(seqname = chr, bin = seq_len(nb), signal = s, control = if (is.null(cc)) NA_real_ else cc)
  })
  bins <- do.call(rbind, bins)
  # Background: trimmed mean of signal windows (peaks excluded by the trim).
  background <- mean(bins$signal[bins$signal <= stats::quantile(bins$signal, 0.9)])
  lambda <- rep(max(background, 0.5), nrow(bins))
  if (!is.null(ctl) && nrow(ctl)) {
    scale <- nrow(sig) / nrow(ctl)
    local <- stats::filter(bins$control, rep(1 / 11, 11), sides = 2)
    local[is.na(local)] <- bins$control[is.na(local)]
    lambda <- pmax(lambda, scale * pmax(bins$control, local))
  }
  p <- stats::ppois(bins$signal - 1, lambda, lower.tail = FALSE)
  bins$q <- stats::p.adjust(p, "BH")
  hit <- bins$q <= fdr & bins$signal >= min_count
  if (!any(hit)) {
    return(data.frame(peak_id = character(), seqname = character(), start = integer(), end = integer(),
                      max_count = integer(), min_q = numeric()))
  }
  b <- bins[hit, , drop = FALSE]
  run <- cumsum(c(TRUE, diff(b$bin) != 1L | b$seqname[-1L] != b$seqname[-nrow(b)]))
  out <- do.call(rbind, lapply(split(b, run), function(x) {
    data.frame(seqname = x$seqname[1L], start = as.integer((min(x$bin) - 1L) * window + 1L),
               end = as.integer(max(x$bin) * window), max_count = max(x$signal), min_q = min(x$q),
               stringsAsFactors = FALSE)
  }))
  out$peak_id <- sprintf("called_%05d", seq_len(nrow(out)))
  rownames(out) <- NULL
  out[, c("peak_id", "seqname", "start", "end", "max_count", "min_q")]
}

#' Score called peaks against true peaks
#'
#' @param called Called peaks (seqname, start, end).
#' @param truth True peaks (seqname, start, end), or a path.
#' @return One-row data frame with recall, precision, and F1.
#' @export
score_peak_calls <- function(called, truth) {
  if (is.character(truth)) truth <- utils::read.delim(truth, stringsAsFactors = FALSE)
  recovered <- if (nrow(called)) sum(.bioana_peak_overlaps(truth, called)) else 0L
  true_calls <- if (nrow(called)) sum(.bioana_peak_overlaps(called, truth)) else 0L
  recall <- recovered / nrow(truth)
  precision <- if (nrow(called)) true_calls / nrow(called) else NA_real_
  data.frame(analysis = "peak calling", true_peaks = nrow(truth), called_peaks = nrow(called),
             recovered_peaks = recovered, false_peaks = nrow(called) - true_calls,
             recall = recall, precision = precision,
             f1 = if (is.na(precision) || recall + precision == 0) 0 else 2 * recall * precision / (recall + precision))
}

#' Test peaks for differential binding or accessibility between two conditions
#'
#' Counts reads per peak and library, then tests the second condition against
#' the first with a per-peak t-test on log-CPM, edgeR (exact test), or DESeq2.
#'
#' @param reads Read-position table or path.
#' @param peaks Peaks to test (seqname, start, end, peak_id).
#' @param assay Assay label of the reads to count.
#' @param method `"ttest"`, `"edger"`, or `"deseq2"`.
#' @param fdr BH-adjusted threshold for calling a differential peak.
#' @return Data frame with peak_id, log2_fold_change, p_value, q_value, significant.
#' @export
test_differential_peaks <- function(reads, peaks, assay = "chip", method = c("ttest", "edger", "deseq2"), fdr = 0.05) {
  method <- match.arg(method)
  reads <- .bioana_peak_reads(reads)
  reads <- reads[reads$assay == assay, , drop = FALSE]
  if (!nrow(peaks)) stop("No peaks to test")
  libraries <- unique(reads$library_id)
  condition <- reads$condition[match(libraries, reads$library_id)]
  levels <- unique(condition)
  if (length(levels) < 2L) stop("Differential testing needs two conditions")
  counts <- matrix(0L, nrow(peaks), length(libraries), dimnames = list(peaks$peak_id, libraries))
  for (chr in unique(peaks$seqname)) {
    pk <- which(peaks$seqname == chr)
    pk <- pk[order(peaks$start[pk])]
    r <- reads[reads$seqname == chr, , drop = FALSE]
    idx <- findInterval(r$read_position, peaks$start[pk])
    ok <- idx > 0L
    ok[ok] <- r$read_position[ok] <= peaks$end[pk][idx[ok]]
    if (any(ok)) {
      tab <- table(factor(peaks$peak_id[pk][idx[ok]], levels = peaks$peak_id), factor(r$library_id[ok], levels = libraries))
      counts <- counts + matrix(as.integer(tab), nrow(peaks))
    }
  }
  group <- factor(condition, levels = levels[1:2])
  keep <- !is.na(group)
  counts <- counts[, keep, drop = FALSE]; group <- droplevels(group[keep])
  if (min(table(group)) < 2L) stop("Each condition needs at least two replicates")
  result <- switch(method,
    ttest = {
      lib <- colSums(counts)
      logcpm <- log2(t(t(counts + 0.5) / (lib + 1)) * 1e6)
      a <- group == levels(group)[1L]
      stats <- t(apply(logcpm, 1L, function(x) {
        if (stats::sd(x[a]) + stats::sd(x[!a]) == 0) return(c(mean(x[!a]) - mean(x[a]), 1))
        c(mean(x[!a]) - mean(x[a]), stats::t.test(x[!a], x[a])$p.value)
      }))
      data.frame(log2_fold_change = stats[, 1L], p_value = stats[, 2L])
    },
    edger = {
      if (!requireNamespace("edgeR", quietly = TRUE)) stop("method 'edger' needs the edgeR package")
      y <- edgeR::DGEList(counts, group = group)
      y <- edgeR::calcNormFactors(y)
      y <- edgeR::estimateDisp(y, stats::model.matrix(~group))
      et <- edgeR::exactTest(y)$table
      data.frame(log2_fold_change = et$logFC, p_value = et$PValue)
    },
    deseq2 = {
      if (!requireNamespace("DESeq2", quietly = TRUE)) stop("method 'deseq2' needs the DESeq2 package")
      dds <- DESeq2::DESeqDataSetFromMatrix(counts, data.frame(group = group), ~group)
      dds <- suppressMessages(DESeq2::DESeq(dds, quiet = TRUE))
      res <- as.data.frame(DESeq2::results(dds))
      data.frame(log2_fold_change = res$log2FoldChange, p_value = ifelse(is.na(res$pvalue), 1, res$pvalue))
    })
  result <- cbind(data.frame(peak_id = peaks$peak_id, stringsAsFactors = FALSE), result)
  result$q_value <- stats::p.adjust(result$p_value, "BH")
  result$significant <- result$q_value <= fdr
  result$method <- method
  result
}

#' Score differential peaks against the simulated truth
#'
#' A called peak is truly differential when it overlaps a true peak with a
#' differential effect. Recall is over true differential peaks; precision is
#' over called peaks declared significant.
#'
#' @param differential Output of test_differential_peaks().
#' @param called Called peaks used in the test.
#' @param truth_peaks True peaks (table or path).
#' @param differential_truth True differential effects (table or path).
#' @return One-row data frame with recall, precision, and F1.
#' @export
score_differential_peaks <- function(differential, called, truth_peaks, differential_truth) {
  if (is.character(truth_peaks)) truth_peaks <- utils::read.delim(truth_peaks, stringsAsFactors = FALSE)
  if (is.character(differential_truth)) differential_truth <- utils::read.delim(differential_truth, stringsAsFactors = FALSE)
  true_diff <- truth_peaks[truth_peaks$peak_id %in% differential_truth$peak_id, , drop = FALSE]
  significant <- called[called$peak_id %in% differential$peak_id[differential$significant], , drop = FALSE]
  recovered <- if (nrow(significant) && nrow(true_diff)) sum(.bioana_peak_overlaps(true_diff, significant)) else 0L
  correct <- if (nrow(significant) && nrow(true_diff)) sum(.bioana_peak_overlaps(significant, true_diff)) else 0L
  recall <- if (nrow(true_diff)) recovered / nrow(true_diff) else NA_real_
  precision <- if (nrow(significant)) correct / nrow(significant) else NA_real_
  data.frame(analysis = "differential peaks", method = differential$method[1L], true_differential = nrow(true_diff),
             called_differential = nrow(significant), recovered = recovered, false_differential = nrow(significant) - correct,
             recall = recall, precision = precision,
             f1 = if (is.na(recall) || is.na(precision) || recall + precision == 0) 0 else 2 * recall * precision / (recall + precision))
}

#' Plot ChIP-seq simulation results
#'
#' Generate a six-panel publication figure containing genomic peak locations,
#' aggregate peak-centered read enrichment, FRiP, library complexity, peak
#' signal by annotation class, and differential-binding recovery. Panel source
#' tables are exported alongside PDF and PNG files.
#'
#' @param chip_prefix Prefix passed to `simitall::simulate_chipseq()`.
#' @param out_prefix Prefix for figures, summary, and source-data files.
#' @param profile_window Distance on either side of peak centers in the
#'   aggregate enrichment profile.
#' @param profile_bin Width of peak-centered profile bins.
#' @param max_profile_reads Maximum reads sampled for the profile panel.
#' @param width Figure width in inches.
#' @param height Figure height in inches.
#' @param dpi PNG resolution.
#' @param seed Seed used for read subsampling.
#' @param individual_panels_dir Optional directory for standalone PNG/PDF
#'   panels. Defaults to a `figures` folder beside `out_prefix`.
#'
#' @return Invisibly returns the combined figure, panels, source data, and
#'   output paths.
#' @examples
#' \dontrun{
#' plot_chipseq_results(
#'   chip_prefix = "results/chipseq/ctcf",
#'   out_prefix = "results/figures/ctcf_chipseq"
#' )
#' }
#' @export
plot_chipseq_results <- function(
    chip_prefix,
    out_prefix,
    profile_window = 1000L,
    profile_bin = 50L,
    max_profile_reads = 200000L,
    width = 15,
    height = 10,
    dpi = 300,
    seed = 1L,
    individual_panels_dir = file.path(dirname(out_prefix), "figures")) {
  if (!requireNamespace("ggplot2", quietly = TRUE) ||
      !requireNamespace("patchwork", quietly = TRUE)) {
    stop("ggplot2 and patchwork are required for ChIP-seq figures")
  }
  if (profile_window < 1L || profile_bin < 1L || max_profile_reads < 100L) {
    stop("Profile plotting limits must be positive")
  }
  set.seed(seed)
  read_table <- function(suffix, label) {
    path <- paste0(chip_prefix, suffix)
    if (!file.exists(path)) stop(label, " does not exist: ", path)
    utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
  }
  peaks <- read_table(".truth_peaks.tsv", "Peak truth")
  reads <- read_table(".read_positions.tsv", "Read positions")
  metadata <- read_table(".sample_metadata.tsv", "Sample metadata")
  qc <- read_table(".qc.tsv", "QC table")
  peak_counts <- read_table(".peak_counts.tsv", "Peak counts")
  differential <- read_table(
    ".differential_binding_truth.tsv", "Differential-binding truth"
  )
  # Native simulation tables retain condition and read_position; BED-like
  # exports retain start/end instead. Accept either form so the lightweight
  # agent recipe can produce complete QC figures without duplicating reads.
  if (!all(c("start", "end") %in% names(reads)) &&
      "read_position" %in% names(reads)) {
    reads$start <- pmax(0L, as.integer(reads$read_position) - 1L)
    reads$end <- reads$start + 1L
  }
  if (!all(c("start", "end", "seqname", "library_id", "assay") %in%
           names(reads))) {
    stop("Read-position table does not have the expected BED-like columns")
  }
  reads$position <- round((reads$start + reads$end) / 2)
  reads$condition <- metadata$condition[
    match(reads$library_id, metadata$library_id)
  ]

  palette <- c(
    navy = "#153243", blue = "#2D6A8A", cyan = "#4B9DA9",
    green = "#6A994E", sand = "#E5B25D", orange = "#D9783D",
    coral = "#C8553D", gray = "#A8B1B7"
  )
  theme <- ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        face = "bold", color = palette[["navy"]], size = 11
      ),
      plot.subtitle = ggplot2::element_text(color = "#52616B", size = 8.5),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "bottom"
    )

  panel_a_data <- peaks
  panel_a <- ggplot2::ggplot(
    panel_a_data,
    ggplot2::aes(
      x = start / 1e6, xend = end / 1e6,
      y = feature_type, yend = feature_type, color = feature_type
    )
  ) +
    ggplot2::geom_segment(linewidth = 3.5, alpha = 0.75) +
    ggplot2::facet_wrap(~seqname, scales = "free_x") +
    ggplot2::labs(
      title = "A. Programmed binding landscape",
      subtitle = paste(nrow(peaks), "truth peaks"),
      x = "Genomic position (Mb)", y = NULL, color = "Target"
    ) + theme + ggplot2::guides(color = "none")

  profile_reads <- reads
  if (nrow(profile_reads) > max_profile_reads) {
    profile_reads <- profile_reads[
      sort(sample(seq_len(nrow(profile_reads)), max_profile_reads)), ,
      drop = FALSE
    ]
  }
  profile_reads$distance <- .bioana_chip_nearest_distance(
    profile_reads$seqname, profile_reads$position, peaks
  )
  profile_reads <- profile_reads[
    is.finite(profile_reads$distance) &
      abs(profile_reads$distance) <= profile_window, , drop = FALSE
  ]
  breaks <- seq(-profile_window, profile_window + profile_bin, by = profile_bin)
  profile_reads$distance_bin <- breaks[pmax(
    1L, pmin(length(breaks), findInterval(profile_reads$distance, breaks))
  )]
  profile <- stats::aggregate(
    list(reads = rep(1L, nrow(profile_reads))),
    profile_reads[c("distance_bin", "assay")], sum
  )
  assay_totals <- table(reads$assay)
  profile$reads_per_million_per_peak <- profile$reads /
    as.numeric(assay_totals[profile$assay]) * 1e6 / nrow(peaks)
  panel_b_data <- profile
  panel_b <- ggplot2::ggplot(
    panel_b_data,
    ggplot2::aes(
      x = distance_bin, y = reads_per_million_per_peak, color = assay
    )
  ) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::scale_color_manual(values = c(chip = palette[["coral"]],
                                           input = palette[["gray"]])) +
    ggplot2::labs(
      title = "B. Peak-centered enrichment",
      subtitle = "Aggregate read profile around truth peak centers",
      x = "Distance from peak center (bp)",
      y = "Reads per million per peak", color = NULL
    ) + theme

  panel_c_data <- qc[qc$assay == "chip", , drop = FALSE]
  panel_c <- ggplot2::ggplot(
    panel_c_data,
    ggplot2::aes(x = condition, y = frip, color = condition)
  ) +
    ggplot2::geom_boxplot(width = 0.45, outlier.shape = NA) +
    ggplot2::geom_jitter(width = 0.08, size = 1.6) +
    ggplot2::scale_y_continuous(labels = function(x) paste0(round(100 * x), "%")) +
    ggplot2::labs(
      title = "C. Fraction of reads in peaks",
      subtitle = "ChIP libraries only", x = NULL, y = "FRiP", color = NULL
    ) + theme + ggplot2::guides(color = "none")

  panel_d_data <- qc
  panel_d_data$unique_fraction <- panel_d_data$unique_reads /
    pmax(1, panel_d_data$total_reads)
  panel_d <- ggplot2::ggplot(
    panel_d_data,
    ggplot2::aes(x = assay, y = unique_fraction, color = condition)
  ) +
    ggplot2::geom_boxplot(width = 0.45, outlier.shape = NA) +
    ggplot2::geom_jitter(width = 0.08, size = 1.4) +
    ggplot2::scale_y_continuous(labels = function(x) paste0(round(100 * x), "%")) +
    ggplot2::labs(
      title = "D. Library complexity",
      subtitle = "Unique genomic read positions", x = NULL,
      y = "Unique read fraction", color = "Condition"
    ) + theme

  chip_counts <- peak_counts[peak_counts$assay == "chip", , drop = FALSE]
  chip_counts$feature_type <- peaks$feature_type[
    match(chip_counts$peak_id, peaks$peak_id)
  ]
  chip_qc <- qc[qc$assay == "chip", , drop = FALSE]
  library_sizes <- chip_qc$total_reads[
    match(chip_counts$library_id, chip_qc$library_id)
  ]
  chip_counts$counts_per_million <- chip_counts$count /
    pmax(1, library_sizes) * 1e6
  panel_e_data <- chip_counts
  panel_e <- ggplot2::ggplot(
    panel_e_data,
    ggplot2::aes(x = feature_type, y = log2(counts_per_million + 1),
                 fill = feature_type)
  ) +
    ggplot2::geom_violin(trim = FALSE, alpha = 0.65, color = NA) +
    ggplot2::geom_boxplot(width = 0.14, outlier.shape = NA, fill = "white") +
    ggplot2::labs(
      title = "E. Peak signal by target class",
      subtitle = "Read counts within programmed intervals",
      x = NULL, y = "log2 counts per million + 1", fill = NULL
    ) + theme + ggplot2::guides(fill = "none") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 25, hjust = 1))

  if (nrow(differential) && length(unique(chip_counts$condition)) > 1L) {
    reference <- unique(metadata$condition)[1L]
    condition_means <- stats::aggregate(
      counts_per_million ~ peak_id + condition, chip_counts, mean
    )
    observed <- do.call(rbind, lapply(
      unique(differential$condition), function(condition) {
        ref <- condition_means[condition_means$condition == reference, ]
        alt <- condition_means[condition_means$condition == condition, ]
        merged <- merge(ref, alt, by = "peak_id", suffixes = c("_ref", "_alt"))
        merged$condition <- condition
        merged$observed_effect_log2 <- log2(
          (merged$counts_per_million_alt + 0.5) /
            (merged$counts_per_million_ref + 0.5)
        )
        merged
      }
    ))
    panel_f_data <- merge(
      differential, observed[, c("peak_id", "condition", "observed_effect_log2")],
      by = c("peak_id", "condition"), all.x = TRUE
    )
    correlation <- if (nrow(panel_f_data) >= 3L) {
      stats::cor(
        panel_f_data$effect_log2, panel_f_data$observed_effect_log2,
        use = "complete.obs"
      )
    } else NA_real_
    panel_f <- ggplot2::ggplot(
      panel_f_data,
      ggplot2::aes(x = effect_log2, y = observed_effect_log2,
                   color = direction)
    ) +
      ggplot2::geom_abline(slope = 1, intercept = 0, color = palette[["gray"]]) +
      ggplot2::geom_point(size = 1.7, alpha = 0.75) +
      ggplot2::scale_color_manual(values = c(
        gain = palette[["green"]], loss = palette[["coral"]]
      )) +
      ggplot2::labs(
        title = "F. Differential-binding recovery",
        subtitle = sprintf("Truth versus simulated peak counts; r = %.2f", correlation),
        x = "True log2 effect", y = "Observed log2 fold change", color = NULL
      ) + theme
  } else {
    panel_f_data <- peaks
    panel_f <- ggplot2::ggplot(
      panel_f_data,
      ggplot2::aes(x = baseline_enrichment, fill = feature_type)
    ) +
      ggplot2::geom_histogram(bins = 25, alpha = 0.75) +
      ggplot2::labs(
        title = "F. Peak enrichment truth",
        subtitle = "Add multiple conditions for differential binding",
        x = "Baseline enrichment", y = "Peaks", fill = "Target"
      ) + theme
  }

  combined <- (panel_a | panel_b | panel_c) /
    (panel_d | panel_e | panel_f) +
    patchwork::plot_annotation(
      title = "Annotation-aware ChIP-seq simulation",
      subtitle = paste0(
        unique(peaks$assay_type)[1L], " assay | ", nrow(peaks), " peaks | ",
        length(unique(metadata$biological_id)), " biological replicates | ",
        sum(qc$total_reads), " reads"
      ),
      theme = ggplot2::theme(
        plot.title = ggplot2::element_text(
          face = "bold", size = 16, color = palette[["navy"]]
        ),
        plot.subtitle = ggplot2::element_text(size = 10, color = "#52616B")
      )
    )
  source_data <- list(
    panel_a_peak_landscape = panel_a_data,
    panel_b_peak_centered_profile = panel_b_data,
    panel_c_frip = panel_c_data,
    panel_d_library_complexity = panel_d_data,
    panel_e_peak_signal = panel_e_data,
    panel_f_differential_binding = panel_f_data
  )
  source_directory <- paste0(out_prefix, "_source_data")
  dir.create(source_directory, recursive = TRUE, showWarnings = FALSE)
  for (name in names(source_data)) {
    utils::write.csv(
      source_data[[name]], file.path(source_directory, paste0(name, ".csv")),
      row.names = FALSE
    )
  }
  paths <- list(
    pdf = paste0(out_prefix, ".pdf"),
    png = paste0(out_prefix, ".png"),
    summary = paste0(out_prefix, ".summary.tsv"),
    source_data = source_directory
  )
  dir.create(dirname(out_prefix), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(paths$pdf, combined, width = width, height = height)
  ggplot2::ggsave(paths$png, combined, width = width, height = height, dpi = dpi)
  summary <- data.frame(
    peaks = nrow(peaks), libraries = nrow(metadata),
    biological_replicates = length(unique(metadata$biological_id)),
    total_reads = sum(qc$total_reads), mean_chip_frip = mean(panel_c_data$frip),
    mean_unique_fraction = mean(panel_d_data$unique_fraction),
    stringsAsFactors = FALSE
  )
  utils::write.table(
    summary, paths$summary, sep = "\t", row.names = FALSE, quote = FALSE
  )
  individual_panels <- if (is.null(individual_panels_dir)) {
    NULL
  } else {
    export_bioana_plot_panels(
      list(panels = list(panel_a, panel_b, panel_c, panel_d, panel_e, panel_f)),
      individual_panels_dir,
      c("chip_peak_landscape", "chip_peak_profile", "chip_frip",
        "chip_library_complexity", "chip_peak_signal", "chip_differential_binding"),
      width = 7, height = 5, dpi = dpi
    )
  }
  paths$individual_panels <- individual_panels
  invisible(list(
    figure = combined,
    panels = list(panel_a, panel_b, panel_c, panel_d, panel_e, panel_f),
    source_data = source_data,
    paths = paths
  ))
}

.bioana_chip_read_fasta <- function(path) {
  if (!file.exists(path)) stop("Genome FASTA does not exist: ", path)
  lines <- readLines(path, warn = FALSE)
  headers <- which(startsWith(lines, ">"))
  if (!length(headers)) stop("Genome FASTA has no records")
  ends <- c(headers[-1L] - 1L, length(lines))
  ids <- sub("^>([^[:space:]]+).*$", "\\1", lines[headers])
  sequences <- vapply(seq_along(headers), function(i) {
    paste(lines[(headers[i] + 1L):ends[i]], collapse = "")
  }, character(1L))
  sequences <- toupper(gsub("[[:space:]]", "", sequences))
  if (any(!nzchar(sequences))) stop("Genome FASTA contains an empty record")
  stats::setNames(sequences, make.unique(ids))
}

.bioana_peak_reads <- function(reads) {
  if (is.character(reads)) reads <- utils::read.delim(reads, check.names = FALSE, stringsAsFactors = FALSE)
  need <- c("seqname", "read_position", "library_id", "condition", "assay")
  if (any(!need %in% names(reads))) stop("Read table needs columns: ", paste(need, collapse = ", "))
  reads
}

.bioana_peak_overlaps <- function(a, b) {
  vapply(seq_len(nrow(a)), function(i) {
    any(b$seqname == a$seqname[i] & b$start <= a$end[i] & b$end >= a$start[i])
  }, logical(1))
}

.bioana_chip_nearest_distance <- function(seqname, position, peaks) {
  distance <- rep(NA_real_, length(position))
  for (chr in intersect(unique(seqname), unique(peaks$seqname))) {
    read_index <- which(seqname == chr)
    centers <- sort(peaks$center[peaks$seqname == chr])
    insertion <- findInterval(position[read_index], centers)
    left <- pmax(1L, insertion)
    right <- pmin(length(centers), insertion + 1L)
    left_distance <- position[read_index] - centers[left]
    right_distance <- position[read_index] - centers[right]
    use_right <- abs(right_distance) < abs(left_distance)
    distance[read_index] <- ifelse(use_right, right_distance, left_distance)
  }
  distance
}
