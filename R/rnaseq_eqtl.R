# RNA-seq and eQTL results (moved from Simitall): eQTL benchmark, results summary and figure.

#' Benchmark eQTL recovery
#'
#' Test genotype-expression associations and compare discoveries with a known
#' eQTL truth set. Candidate pairs can be restricted to neighborhoods around
#' simulated causal genes, all cis pairs, or all variant-gene pairs. Batch and
#' condition covariates are included when available, and a two-condition model
#' can also test genotype-by-condition interactions.
#'
#' @param genotype_file Genotype TSV or VCF/VCF.GZ.
#' @param expression_file Expression TSV or matrix with genes in rows and
#'   samples in columns. The `.expression.tsv` output from
#'   `simitall::simulate_rnaseq_from_gwas()` is ready to use.
#' @param out_prefix Prefix for association results and benchmark metrics.
#' @param sample_metadata Optional data frame or TSV with sample covariates.
#' @param truth_file Optional eQTL truth TSV written by
#'   `simitall::simulate_rnaseq_from_gwas()`.
#' @param gene_metadata Optional data frame or TSV with `gene_id`, `seqname`,
#'   and `tss`, required for `pair_mode = "cis"`.
#' @param pair_mode Candidate-pair strategy: `"truth_neighborhood"` tests cis
#'   neighborhoods for truth genes plus causal trans pairs, `"cis"` tests all
#'   local pairs, and `"all"` tests every pair.
#' @param cis_window_bp Cis testing window around each TSS.
#' @param covariates Metadata columns to include in each regression.
#' @param test_gxe Test genotype-by-condition interaction when exactly two
#'   conditions are present.
#' @param expression_scale `"log2"`, `"counts"`, or `"auto"`. Count data are
#'   converted to log2 counts per million.
#' @param fdr_threshold Adjusted-p-value threshold used for benchmark metrics.
#' @param max_tests Safety limit on the number of variant-gene tests.
#'
#' @return Invisibly returns a list with result and metric paths and tables.
#' @examples
#' \dontrun{
#' benchmark_eqtl(
#'   genotype_file = "results/gwas/demo.geno.tsv",
#'   expression_file = "results/rnaseq/demo.expression.tsv",
#'   sample_metadata = "results/rnaseq/demo.sample_metadata.tsv",
#'   truth_file = "results/rnaseq/demo.eqtl_truth.tsv",
#'   out_prefix = "results/eqtl/demo"
#' )
#' }
#' @export
benchmark_eqtl <- function(
    genotype_file,
    expression_file,
    out_prefix,
    sample_metadata = NULL,
    truth_file = NULL,
    gene_metadata = NULL,
    pair_mode = c("truth_neighborhood", "cis", "all"),
    cis_window_bp = 1e6,
    covariates = c("condition", "batch"),
    test_gxe = TRUE,
    expression_scale = c("auto", "log2", "counts"),
    fdr_threshold = 0.05,
    max_tests = 1e6) {
  pair_mode <- match.arg(pair_mode)
  expression_scale <- match.arg(expression_scale)
  genotype_data <- .bioana_read_genotypes(genotype_file)
  variants <- genotype_data$variants
  genotype <- genotype_data$genotype
  expression <- .bioana_read_expression_matrix(expression_file)
  sample_ids <- intersect(colnames(genotype), colnames(expression))
  if (length(sample_ids) < 4L) {
    stop("At least four shared genotype and expression samples are required")
  }
  genotype <- genotype[, sample_ids, drop = FALSE]
  expression <- expression[, sample_ids, drop = FALSE]
  metadata <- .bioana_read_sample_metadata(sample_metadata, sample_ids)

  if (expression_scale == "auto") {
    finite_values <- expression[is.finite(expression)]
    looks_like_counts <- length(finite_values) &&
      all(finite_values >= 0) && all(abs(finite_values - round(finite_values)) < 1e-8) &&
      max(finite_values) > 50
    expression_scale <- if (looks_like_counts) "counts" else "log2"
  }
  if (expression_scale == "counts") {
    library_sizes <- pmax(1, colSums(expression, na.rm = TRUE))
    expression <- log2(sweep(expression, 2L, library_sizes, "/") * 1e6 + 1)
  }

  read_table <- function(value, label) {
    if (is.null(value)) return(NULL)
    if (is.character(value) && length(value) == 1L) {
      if (!file.exists(value)) stop(label, " does not exist: ", value)
      value <- utils::read.delim(
        value, header = TRUE, sep = "\t", check.names = FALSE,
        stringsAsFactors = FALSE
      )
    }
    as.data.frame(value, stringsAsFactors = FALSE)
  }
  truth <- read_table(truth_file, "Truth file")
  gene_table <- read_table(gene_metadata, "Gene metadata")

  if (pair_mode == "truth_neighborhood" && is.null(truth)) {
    stop("truth_file is required for pair_mode = 'truth_neighborhood'")
  }
  if (pair_mode == "cis" && is.null(gene_table)) {
    stop("gene_metadata is required for pair_mode = 'cis'")
  }

  make_cis_pairs <- function(gene_rows) {
    rows <- vector("list", nrow(gene_rows))
    for (i in seq_len(nrow(gene_rows))) {
      variant_index <- which(
        variants$seqname == gene_rows$seqname[i] &
          abs(variants$pos - gene_rows$tss[i]) <= cis_window_bp
      )
      rows[[i]] <- if (length(variant_index)) {
        data.frame(
          gene_id = gene_rows$gene_id[i],
          variant_id = variants$id[variant_index],
          stringsAsFactors = FALSE
        )
      } else NULL
    }
    rows <- rows[!vapply(rows, is.null, logical(1L))]
    if (length(rows)) do.call(rbind, rows) else data.frame(
      gene_id = character(), variant_id = character(), stringsAsFactors = FALSE
    )
  }

  if (pair_mode == "truth_neighborhood") {
    required_truth <- c("gene_id", "variant_id", "seqname", "gene_tss")
    if (!all(required_truth %in% names(truth))) {
      stop("truth_file is missing required columns: ", paste(required_truth, collapse = ", "))
    }
    truth_genes <- unique(truth[, c("gene_id", "seqname", "gene_tss")])
    names(truth_genes)[names(truth_genes) == "gene_tss"] <- "tss"
    pairs <- make_cis_pairs(truth_genes)
    pairs <- unique(rbind(pairs, truth[, c("gene_id", "variant_id")]))
  } else if (pair_mode == "cis") {
    required_genes <- c("gene_id", "seqname", "tss")
    if (!all(required_genes %in% names(gene_table))) {
      stop("gene_metadata is missing required columns: ", paste(required_genes, collapse = ", "))
    }
    gene_table <- gene_table[gene_table$gene_id %in% rownames(expression), , drop = FALSE]
    pairs <- make_cis_pairs(gene_table)
  } else {
    pairs <- expand.grid(
      gene_id = rownames(expression),
      variant_id = variants$id,
      KEEP.OUT.ATTRS = FALSE,
      stringsAsFactors = FALSE
    )
  }
  pairs <- pairs[
    pairs$gene_id %in% rownames(expression) &
      pairs$variant_id %in% variants$id,
    , drop = FALSE
  ]
  pairs <- unique(pairs)
  if (!nrow(pairs)) {
    stop("No testable variant-gene pairs were found")
  }
  if (nrow(pairs) > max_tests) {
    stop(
      "Candidate set contains ", nrow(pairs), " tests, exceeding max_tests = ",
      max_tests, ". Reduce the cis window or raise max_tests explicitly."
    )
  }

  usable_covariates <- intersect(covariates, names(metadata))
  if (length(usable_covariates)) {
    covariate_data <- metadata[, usable_covariates, drop = FALSE]
    for (name in names(covariate_data)) {
      if (is.character(covariate_data[[name]])) {
        preferred_code <- paste0(name, "_code")
        if (preferred_code %in% names(metadata)) {
          level_order <- unique(as.character(metadata[[name]])[
            order(metadata[[preferred_code]])
          ])
          covariate_data[[name]] <- factor(
            covariate_data[[name]], levels = level_order
          )
        } else {
          covariate_data[[name]] <- factor(covariate_data[[name]])
        }
      }
    }
    covariate_matrix <- stats::model.matrix(~ ., data = covariate_data)
  } else {
    covariate_matrix <- matrix(1, nrow = length(sample_ids), ncol = 1L)
    colnames(covariate_matrix) <- "(Intercept)"
  }
  interaction <- NULL
  if (isTRUE(test_gxe) && "condition" %in% names(metadata) &&
      length(unique(metadata$condition)) == 2L) {
    if ("condition_code" %in% names(metadata)) {
      interaction <- as.numeric(metadata$condition_code > min(metadata$condition_code))
    } else {
      condition_values <- unique(as.character(metadata$condition))
      interaction <- as.numeric(metadata$condition == condition_values[2L])
    }
  }

  gene_index <- match(pairs$gene_id, rownames(expression))
  variant_index <- match(pairs$variant_id, variants$id)
  fits <- matrix(NA_real_, nrow = nrow(pairs), ncol = 5L)
  colnames(fits) <- c("beta", "se", "p_value", "gxe_beta", "gxe_p_value")
  for (i in seq_len(nrow(pairs))) {
    fits[i, ] <- .bioana_fit_eqtl(
      y = expression[gene_index[i], ],
      genotype = genotype[variant_index[i], ],
      covariate_matrix = covariate_matrix,
      interaction = interaction
    )
  }
  results <- cbind(pairs, as.data.frame(fits, stringsAsFactors = FALSE))
  results$q_value <- stats::p.adjust(results$p_value, method = "BH")
  results$gxe_q_value <- stats::p.adjust(results$gxe_p_value, method = "BH")
  results$variant_seqname <- variants$seqname[variant_index]
  results$variant_pos <- variants$pos[variant_index]

  truth_keys <- character()
  if (!is.null(truth)) {
    truth_keys <- unique(paste(truth$gene_id, truth$variant_id, sep = "::"))
  }
  result_keys <- paste(results$gene_id, results$variant_id, sep = "::")
  results$is_causal <- result_keys %in% truth_keys
  discovered <- !is.na(results$q_value) & results$q_value <= fdr_threshold
  true_positive <- sum(discovered & results$is_causal)
  false_positive <- sum(discovered & !results$is_causal)
  truth_in_tests <- sum(results$is_causal)
  metrics <- data.frame(
    tests = nrow(results),
    truth_pairs_in_tests = truth_in_tests,
    discoveries = sum(discovered),
    true_positives = true_positive,
    false_positives = false_positive,
    false_negatives = truth_in_tests - true_positive,
    precision = if (sum(discovered)) true_positive / sum(discovered) else NA_real_,
    recall = if (truth_in_tests) true_positive / truth_in_tests else NA_real_,
    empirical_fdr = if (sum(discovered)) false_positive / sum(discovered) else NA_real_,
    fdr_threshold = fdr_threshold,
    stringsAsFactors = FALSE
  )

  dir.create(dirname(out_prefix), recursive = TRUE, showWarnings = FALSE)
  result_path <- paste0(out_prefix, ".eqtl_results.tsv")
  metric_path <- paste0(out_prefix, ".eqtl_metrics.tsv")
  utils::write.table(
    results, result_path, sep = "\t", row.names = FALSE, quote = FALSE
  )
  utils::write.table(
    metrics, metric_path, sep = "\t", row.names = FALSE, quote = FALSE
  )
  message(
    "eQTL benchmark complete: ", nrow(results), " tests and ",
    sum(discovered), " discoveries at FDR <= ", fdr_threshold, "."
  )
  invisible(list(
    results_path = result_path,
    metrics_path = metric_path,
    results = results,
    metrics = metrics
  ))
}

#' Summarize RNA-seq and eQTL simulation results
#'
#' Collect key dimensions, count quality metrics, causal-feature counts, and
#' eQTL recovery statistics from the files written by
#' `simitall::simulate_rnaseq_from_gwas()` and [benchmark_eqtl()]. The resulting
#' one-row table is suitable for result reports and simulation manifests.
#'
#' @param rnaseq_prefix Prefix passed to `simitall::simulate_rnaseq_from_gwas()`.
#' @param eqtl_prefix Optional prefix passed to [benchmark_eqtl()].
#' @param out_tsv Optional output TSV path.
#'
#' @return A one-row data frame, invisibly when `out_tsv` is supplied.
#' @examples
#' \dontrun{
#' summarize_rnaseq_results(
#'   rnaseq_prefix = "results/rnaseq/demo",
#'   eqtl_prefix = "results/eqtl/demo",
#'   out_tsv = "results/figures/rnaseq_summary.tsv"
#' )
#' }
#' @export
summarize_rnaseq_results <- function(
    rnaseq_prefix,
    eqtl_prefix = NULL,
    out_tsv = NULL) {
  counts <- .bioana_read_expression_matrix(
    paste0(rnaseq_prefix, ".counts.tsv")
  )
  metadata <- .bioana_results_table(
    paste0(rnaseq_prefix, ".sample_metadata.tsv"),
    "RNA-seq sample metadata"
  )
  truth <- .bioana_results_table(
    paste0(rnaseq_prefix, ".eqtl_truth.tsv"),
    "eQTL truth table"
  )
  ase <- .bioana_results_table(
    paste0(rnaseq_prefix, ".ase_counts.tsv"),
    "ASE count table"
  )

  observed_library <- if ("observed_library_size" %in% names(metadata)) {
    metadata$observed_library_size
  } else {
    colSums(counts)
  }
  result <- data.frame(
    samples = ncol(counts),
    genes = nrow(counts),
    median_library_size = stats::median(observed_library),
    zero_fraction = mean(counts == 0),
    cis_eqtls = sum(truth$type == "cis"),
    trans_eqtls = sum(truth$type == "trans"),
    gxe_eqtls = sum(truth$gxe_effect_log2 != 0),
    ase_pairs = length(unique(paste(ase$gene_id, ase$variant_id))),
    ase_observations = nrow(ase),
    eqtl_tests = NA_integer_,
    discoveries = NA_integer_,
    precision = NA_real_,
    recall = NA_real_,
    empirical_fdr = NA_real_,
    stringsAsFactors = FALSE
  )

  if (!is.null(eqtl_prefix)) {
    metrics <- .bioana_results_table(
      paste0(eqtl_prefix, ".eqtl_metrics.tsv"),
      "eQTL benchmark metrics"
    )
    if (nrow(metrics)) {
      result$eqtl_tests <- metrics$tests[1L]
      result$discoveries <- metrics$discoveries[1L]
      result$precision <- metrics$precision[1L]
      result$recall <- metrics$recall[1L]
      result$empirical_fdr <- metrics$empirical_fdr[1L]
    }
  }

  if (!is.null(out_tsv)) {
    dir.create(dirname(out_tsv), recursive = TRUE, showWarnings = FALSE)
    utils::write.table(
      result, out_tsv, sep = "\t", row.names = FALSE, quote = FALSE
    )
    return(invisible(result))
  }
  result
}

#' Plot RNA-seq and eQTL simulation results
#'
#' Build a six-panel publication figure directly from simulated RNA-seq,
#' causal-truth, ASE, and eQTL benchmark files. Panels show library sizes,
#' expression PCA, association evidence, causal-effect recovery, ASE recovery,
#' and benchmark performance. Source data for every panel are exported beside
#' the PDF and PNG figures.
#'
#' @param rnaseq_prefix Prefix passed to `simitall::simulate_rnaseq_from_gwas()`.
#' @param eqtl_prefix Prefix passed to [benchmark_eqtl()].
#' @param out_prefix Prefix for PDF, PNG, summary TSV, and source-data outputs.
#' @param top_variable_genes Number of variable genes used for expression PCA.
#' @param fdr_threshold FDR threshold used to label discoveries.
#' @param width Figure width in inches.
#' @param height Figure height in inches.
#' @param dpi PNG resolution.
#' @param individual_panels_dir Optional directory for standalone PNG/PDF
#'   panels. Defaults to a `figures` folder beside `out_prefix`.
#'
#' @return Invisibly returns a list containing the combined figure, individual
#'   panels, source data, and output paths.
#' @examples
#' \dontrun{
#' plot_rnaseq_eqtl_results(
#'   rnaseq_prefix = "results/rnaseq/demo",
#'   eqtl_prefix = "results/eqtl/demo",
#'   out_prefix = "results/figures/rnaseq_eqtl"
#' )
#' }
#' @export
plot_rnaseq_eqtl_results <- function(
    rnaseq_prefix,
    eqtl_prefix,
    out_prefix,
    top_variable_genes = 500L,
    fdr_threshold = 0.05,
    width = 15,
    height = 10,
    dpi = 300,
    individual_panels_dir = file.path(dirname(out_prefix), "figures")) {
  if (!requireNamespace("ggplot2", quietly = TRUE) ||
      !requireNamespace("patchwork", quietly = TRUE)) {
    stop(
      "ggplot2 and patchwork are required for result figures. Install them ",
      "with install.packages(c('ggplot2', 'patchwork'))."
    )
  }
  if (top_variable_genes < 2L) {
    stop("top_variable_genes must be at least two")
  }
  if (fdr_threshold <= 0 || fdr_threshold > 1) {
    stop("fdr_threshold must be greater than zero and at most one")
  }

  counts <- .bioana_read_expression_matrix(
    paste0(rnaseq_prefix, ".counts.tsv")
  )
  expression <- .bioana_read_expression_matrix(
    paste0(rnaseq_prefix, ".expression.tsv")
  )
  metadata <- .bioana_results_table(
    paste0(rnaseq_prefix, ".sample_metadata.tsv"),
    "RNA-seq sample metadata"
  )
  truth <- .bioana_results_table(
    paste0(rnaseq_prefix, ".eqtl_truth.tsv"),
    "eQTL truth table"
  )
  ase_truth <- .bioana_results_table(
    paste0(rnaseq_prefix, ".ase_truth.tsv"),
    "ASE truth table"
  )
  ase_counts <- .bioana_results_table(
    paste0(rnaseq_prefix, ".ase_counts.tsv"),
    "ASE count table"
  )
  associations <- .bioana_results_table(
    paste0(eqtl_prefix, ".eqtl_results.tsv"),
    "eQTL association results"
  )
  metrics <- .bioana_results_table(
    paste0(eqtl_prefix, ".eqtl_metrics.tsv"),
    "eQTL benchmark metrics"
  )

  shared_samples <- Reduce(
    intersect,
    list(colnames(counts), colnames(expression), metadata$sample)
  )
  if (length(shared_samples) < 3L) {
    stop("At least three samples must be shared across result files")
  }
  counts <- counts[, shared_samples, drop = FALSE]
  expression <- expression[, shared_samples, drop = FALSE]
  metadata <- metadata[match(shared_samples, metadata$sample), , drop = FALSE]
  metadata$condition <- if ("condition" %in% names(metadata)) {
    as.character(metadata$condition)
  } else {
    "all samples"
  }
  metadata$batch <- if ("batch" %in% names(metadata)) {
    as.character(metadata$batch)
  } else {
    "batch 1"
  }
  metadata$library_size <- if ("observed_library_size" %in% names(metadata)) {
    metadata$observed_library_size
  } else {
    colSums(counts)
  }

  palette <- c(
    navy = "#17324D",
    blue = "#2F6690",
    teal = "#3A7D78",
    sand = "#E7C27D",
    coral = "#D95D4F",
    gray = "#A7B0B7"
  )
  result_theme <- ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        face = "bold", color = palette[["navy"]], size = 11
      ),
      plot.subtitle = ggplot2::element_text(color = "#52616B", size = 8.5),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "bottom",
      legend.title = ggplot2::element_text(face = "bold")
    )

  panel_a_data <- metadata
  panel_a <- ggplot2::ggplot(
    panel_a_data,
    ggplot2::aes(x = condition, y = library_size, fill = condition)
  ) +
    ggplot2::geom_violin(trim = FALSE, alpha = 0.75, color = NA) +
    ggplot2::geom_boxplot(width = 0.16, outlier.shape = NA, fill = "white") +
    ggplot2::geom_jitter(width = 0.08, size = 0.8, alpha = 0.35) +
    ggplot2::scale_y_log10(labels = function(value) format(value, scientific = TRUE)) +
    ggplot2::scale_fill_manual(values = rep(
      c(palette[["blue"]], palette[["coral"]], palette[["teal"]], palette[["sand"]]),
      length.out = length(unique(panel_a_data$condition))
    )) +
    ggplot2::labs(
      title = "A. RNA-seq library sizes",
      x = NULL,
      y = "Observed counts (log scale)"
    ) +
    result_theme +
    ggplot2::theme(legend.position = "none")

  gene_variance <- apply(expression, 1L, stats::var, na.rm = TRUE)
  variable_genes <- which(is.finite(gene_variance) & gene_variance > 0)
  if (length(variable_genes) < 2L) {
    stop("Expression matrix needs at least two variable genes for PCA")
  }
  variable_genes <- variable_genes[
    order(gene_variance[variable_genes], decreasing = TRUE)
  ]
  selected_genes <- variable_genes[
    seq_len(min(top_variable_genes, length(variable_genes)))
  ]
  pca <- stats::prcomp(t(expression[selected_genes, , drop = FALSE]), scale. = TRUE)
  explained <- 100 * pca$sdev^2 / sum(pca$sdev^2)
  panel_b_data <- data.frame(
    sample = rownames(pca$x),
    PC1 = pca$x[, 1L],
    PC2 = pca$x[, 2L],
    condition = metadata$condition[match(rownames(pca$x), metadata$sample)],
    batch = metadata$batch[match(rownames(pca$x), metadata$sample)],
    stringsAsFactors = FALSE
  )
  panel_b <- ggplot2::ggplot(
    panel_b_data,
    ggplot2::aes(x = PC1, y = PC2, color = condition, shape = batch)
  ) +
    ggplot2::geom_point(size = 2, alpha = 0.8) +
    ggplot2::scale_color_manual(values = rep(
      c(palette[["blue"]], palette[["coral"]], palette[["teal"]], palette[["sand"]]),
      length.out = length(unique(panel_b_data$condition))
    )) +
    ggplot2::labs(
      title = "B. Expression structure",
      subtitle = paste(length(selected_genes), "most-variable genes"),
      x = sprintf("PC1 (%.1f%%)", explained[1L]),
      y = sprintf("PC2 (%.1f%%)", explained[2L]),
      color = "Condition",
      shape = "Batch"
    ) +
    result_theme

  associations$is_causal <- as.logical(associations$is_causal)
  associations$neglog10_p <- -log10(pmax(associations$p_value, .Machine$double.xmin))
  associations$discovery <- !is.na(associations$q_value) &
    associations$q_value <= fdr_threshold
  associations$causal_label <- ifelse(associations$is_causal, "Causal", "Background")
  associations$variant_mb <- associations$variant_pos / 1e6
  panel_c_data <- associations
  panel_c <- ggplot2::ggplot(
    panel_c_data,
    ggplot2::aes(
      x = variant_mb,
      y = neglog10_p,
      color = causal_label,
      shape = discovery
    )
  ) +
    ggplot2::geom_point(alpha = 0.65, size = 1.25, na.rm = TRUE) +
    ggplot2::facet_wrap(~variant_seqname, scales = "free_x") +
    ggplot2::scale_color_manual(values = c(
      Background = palette[["gray"]], Causal = palette[["coral"]]
    )) +
    ggplot2::scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 17)) +
    ggplot2::labs(
      title = "C. eQTL association landscape",
      subtitle = paste0("Triangles: FDR <= ", fdr_threshold),
      x = "Variant position (Mb)",
      y = expression(-log[10](p)),
      color = "Truth",
      shape = "Discovery"
    ) +
    result_theme

  effect_data <- merge(
    truth[, c("gene_id", "variant_id", "type", "effect_log2")],
    associations[, c("gene_id", "variant_id", "beta", "q_value")],
    by = c("gene_id", "variant_id"),
    all.x = TRUE
  )
  effect_complete <- is.finite(effect_data$effect_log2) & is.finite(effect_data$beta)
  effect_correlation <- if (sum(effect_complete) >= 3L) {
    stats::cor(effect_data$effect_log2[effect_complete], effect_data$beta[effect_complete])
  } else {
    NA_real_
  }
  effect_rmse <- if (any(effect_complete)) {
    sqrt(mean((effect_data$effect_log2[effect_complete] - effect_data$beta[effect_complete])^2))
  } else {
    NA_real_
  }
  panel_d_data <- effect_data
  if (any(effect_complete)) {
    panel_d <- ggplot2::ggplot(
      panel_d_data[effect_complete, , drop = FALSE],
      ggplot2::aes(x = effect_log2, y = beta, color = type)
    ) +
      ggplot2::geom_abline(
        slope = 1, intercept = 0, linetype = "dashed", color = "#52616B"
      ) +
      ggplot2::geom_point(size = 2, alpha = 0.8) +
      ggplot2::scale_color_manual(values = c(
        cis = palette[["blue"]], trans = palette[["coral"]]
      )) +
      ggplot2::coord_equal() +
      ggplot2::labs(
        title = "D. Causal effect recovery",
        subtitle = sprintf("Correlation = %.2f; RMSE = %.2f", effect_correlation, effect_rmse),
        x = "True effect (log2)",
        y = "Estimated effect",
        color = "eQTL type"
      ) +
      result_theme
  } else {
    panel_d <- .bioana_empty_panel(
      "D. Causal effect recovery", "No finite causal estimates"
    )
  }

  ase_summary <- data.frame()
  if (nrow(ase_counts) && nrow(ase_truth)) {
    ase_aggregate <- stats::aggregate(
      cbind(alt_count, total_count) ~ gene_id + variant_id,
      data = ase_counts,
      FUN = sum,
      na.rm = TRUE
    )
    ase_aggregate$observed_alt_fraction <- with(
      ase_aggregate,
      ifelse(total_count > 0, alt_count / total_count, NA_real_)
    )
    ase_summary <- merge(
      ase_truth,
      ase_aggregate,
      by = c("gene_id", "variant_id"),
      all.x = TRUE
    )
  }
  ase_complete <- nrow(ase_summary) && any(
    is.finite(ase_summary$alt_probability) &
      is.finite(ase_summary$observed_alt_fraction)
  )
  panel_e_data <- ase_summary
  if (ase_complete) {
    panel_e <- ggplot2::ggplot(
      ase_summary,
      ggplot2::aes(
        x = alt_probability,
        y = observed_alt_fraction,
        size = total_count
      )
    ) +
      ggplot2::geom_abline(
        slope = 1, intercept = 0, linetype = "dashed", color = "#52616B"
      ) +
      ggplot2::geom_point(
        color = palette[["teal"]], alpha = 0.75, na.rm = TRUE
      ) +
      ggplot2::coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
      ggplot2::labs(
        title = "E. Allele-specific expression",
        x = "Expected alternate fraction",
        y = "Observed alternate fraction",
        size = "Allelic depth"
      ) +
      result_theme
  } else {
    panel_e <- .bioana_empty_panel(
      "E. Allele-specific expression", "No heterozygous ASE observations"
    )
  }

  performance_data <- data.frame(
    metric = c("Precision", "Recall", "1 - empirical FDR"),
    value = c(
      metrics$precision[1L],
      metrics$recall[1L],
      1 - metrics$empirical_fdr[1L]
    ),
    stringsAsFactors = FALSE
  )
  performance_data$metric <- factor(
    performance_data$metric,
    levels = performance_data$metric
  )
  panel_f_data <- performance_data
  panel_f <- ggplot2::ggplot(
    performance_data,
    ggplot2::aes(x = metric, y = value, fill = metric)
  ) +
    ggplot2::geom_col(width = 0.68, color = "white", na.rm = TRUE) +
    ggplot2::geom_text(
      ggplot2::aes(label = ifelse(is.finite(value), sprintf("%.2f", value), "NA")),
      vjust = -0.35,
      fontface = "bold"
    ) +
    ggplot2::scale_fill_manual(values = c(
      palette[["blue"]], palette[["teal"]], palette[["sand"]]
    )) +
    ggplot2::coord_cartesian(ylim = c(0, 1.12), clip = "off") +
    ggplot2::labs(
      title = "F. eQTL benchmark performance",
      subtitle = paste(
        metrics$discoveries[1L], "discoveries from", metrics$tests[1L], "tests"
      ),
      x = NULL,
      y = "Score"
    ) +
    result_theme +
    ggplot2::theme(
      legend.position = "none",
      axis.text.x = ggplot2::element_text(angle = 18, hjust = 1)
    )

  panels <- list(panel_a, panel_b, panel_c, panel_d, panel_e, panel_f)
  figure <- patchwork::wrap_plots(panels, ncol = 3, guides = "collect") +
    patchwork::plot_annotation(
      title = "GWAS-linked RNA-seq and eQTL simulation",
      subtitle = paste(
        ncol(counts), "individuals,", nrow(counts), "genes,",
        nrow(truth), "causal eQTLs"
      ),
      theme = ggplot2::theme(
        plot.title = ggplot2::element_text(
          face = "bold", size = 16, color = palette[["navy"]]
        ),
        plot.subtitle = ggplot2::element_text(size = 10, color = "#52616B")
      )
    )

  output_paths <- list(
    pdf = paste0(out_prefix, ".pdf"),
    png = paste0(out_prefix, ".png"),
    summary = paste0(out_prefix, ".summary.tsv"),
    source_data = paste0(out_prefix, "_source_data")
  )
  dir.create(dirname(out_prefix), recursive = TRUE, showWarnings = FALSE)
  dir.create(output_paths$source_data, recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(
    output_paths$pdf, figure, width = width, height = height, units = "in"
  )
  ggplot2::ggsave(
    output_paths$png, figure, width = width, height = height,
    units = "in", dpi = dpi
  )
  summarize_rnaseq_results(rnaseq_prefix, eqtl_prefix, output_paths$summary)

  source_data <- list(
    library_sizes = panel_a_data,
    expression_pca = panel_b_data,
    eqtl_associations = panel_c_data,
    effect_recovery = panel_d_data,
    ase_recovery = panel_e_data,
    benchmark_performance = panel_f_data
  )
  for (name in names(source_data)) {
    utils::write.csv(
      source_data[[name]],
      file.path(output_paths$source_data, paste0(name, ".csv")),
      row.names = FALSE
    )
  }
  individual_panels <- if (is.null(individual_panels_dir)) {
    NULL
  } else {
    export_bioana_plot_panels(
      list(panels = panels), individual_panels_dir,
      c("rnaseq_library_size", "rnaseq_pca", "eqtl_associations",
        "eqtl_effect_recovery", "ase_recovery", "eqtl_benchmark"),
      width = 7, height = 5, dpi = dpi
    )
  }

  message("RNA-seq/eQTL figure written to: ", output_paths$pdf)
  invisible(list(
    figure = figure,
    panels = panels,
    source_data = source_data,
    paths = c(output_paths, list(individual_panels = individual_panels))
  ))
}

.bioana_read_genotypes <- function(path) {
  if (!file.exists(path)) {
    stop("Genotype file does not exist: ", path)
  }

  is_vcf <- grepl("\\.vcf(?:\\.gz)?$", path, ignore.case = TRUE)
  if (is_vcf) {
    con <- if (grepl("\\.gz$", path, ignore.case = TRUE)) {
      gzfile(path, open = "rt")
    } else {
      file(path, open = "rt")
    }
    on.exit(close(con), add = TRUE)
    lines <- readLines(con, warn = FALSE)
    header_index <- grep("^#CHROM\\t", lines)
    if (length(header_index) != 1L) {
      stop("VCF must contain exactly one #CHROM header line")
    }

    header <- strsplit(sub("^#", "", lines[header_index]), "\\t")[[1L]]
    sample_ids <- header[-seq_len(9L)]
    body <- if (header_index < length(lines)) {
      lines[seq.int(header_index + 1L, length(lines))]
    } else {
      character()
    }
    body <- body[nzchar(body) & !startsWith(body, "#")]
    if (!length(body)) {
      stop("VCF contains no variant records")
    }

    fields <- strsplit(body, "\\t", fixed = FALSE)
    if (any(lengths(fields) < 9L + length(sample_ids))) {
      stop("VCF record has fewer columns than the #CHROM header")
    }

    variants <- data.frame(
      seqname = vapply(fields, `[[`, character(1L), 1L),
      pos = as.integer(vapply(fields, `[[`, character(1L), 2L)),
      id = vapply(fields, `[[`, character(1L), 3L),
      ref = vapply(fields, `[[`, character(1L), 4L),
      alt = vapply(fields, `[[`, character(1L), 5L),
      stringsAsFactors = FALSE
    )
    missing_id <- is.na(variants$id) | variants$id %in% c("", ".")
    variants$id[missing_id] <- paste0(
      variants$seqname[missing_id], ":", variants$pos[missing_id]
    )

    parse_gt <- function(value) {
      gt <- strsplit(value, ":", fixed = TRUE)[[1L]][1L]
      if (!nzchar(gt) || grepl("\\.", gt)) {
        return(NA_real_)
      }
      alleles <- strsplit(gt, "[/|]")[[1L]]
      alleles <- suppressWarnings(as.numeric(alleles))
      if (anyNA(alleles)) NA_real_ else sum(alleles > 0)
    }

    genotype <- matrix(
      NA_real_,
      nrow = length(fields),
      ncol = length(sample_ids),
      dimnames = list(variants$id, sample_ids)
    )
    for (i in seq_along(fields)) {
      genotype[i, ] <- vapply(
        fields[[i]][seq.int(10L, 9L + length(sample_ids))],
        parse_gt,
        numeric(1L)
      )
    }
  } else {
    table <- utils::read.delim(
      path,
      header = TRUE,
      sep = "\t",
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
    if (!nrow(table)) {
      stop("Genotype table contains no variants")
    }

    lower_names <- tolower(names(table))
    find_column <- function(candidates, required = FALSE) {
      hit <- match(candidates, lower_names, nomatch = 0L)
      hit <- hit[hit > 0L]
      if (!length(hit)) {
        if (required) {
          stop("Genotype table is missing column: ", candidates[1L])
        }
        return(NULL)
      }
      names(table)[hit[1L]]
    }

    id_col <- find_column(c("id", "variant_id", "marker"), required = TRUE)
    pos_col <- find_column(c("pos", "position", "position_bp"), required = TRUE)
    seq_col <- find_column(c("seqname", "chrom", "chr", "#chrom"))
    ref_col <- find_column("ref")
    alt_col <- find_column("alt")
    metadata_names <- unique(stats::na.omit(c(
      id_col, pos_col, seq_col, ref_col, alt_col,
      names(table)[lower_names %in% c("type", "len", "length", "qual", "filter", "info")]
    )))
    sample_ids <- setdiff(names(table), metadata_names)
    if (!length(sample_ids)) {
      stop("Genotype table contains no sample columns")
    }

    genotype <- as.matrix(table[, sample_ids, drop = FALSE])
    suppressWarnings(storage.mode(genotype) <- "numeric")
    if (all(is.na(genotype))) {
      stop("Sample columns in the genotype table are not numeric")
    }
    variants <- data.frame(
      seqname = if (is.null(seq_col)) "chr1" else as.character(table[[seq_col]]),
      pos = as.integer(table[[pos_col]]),
      id = as.character(table[[id_col]]),
      ref = if (is.null(ref_col)) "N" else as.character(table[[ref_col]]),
      alt = if (is.null(alt_col)) "N" else as.character(table[[alt_col]]),
      stringsAsFactors = FALSE
    )
    rownames(genotype) <- variants$id
    colnames(genotype) <- sample_ids
  }

  if (anyDuplicated(variants$id)) {
    variants$id <- make.unique(variants$id, sep = "_")
    rownames(genotype) <- variants$id
  }
  if (anyDuplicated(colnames(genotype))) {
    stop("Sample IDs must be unique")
  }
  if (anyNA(variants$pos)) {
    stop("Variant positions must be numeric and non-missing")
  }

  list(variants = variants, genotype = genotype)
}

.bioana_empty_panel <- function(title, message) {
  ggplot2::ggplot() +
    ggplot2::annotate(
      "text", x = 0.5, y = 0.5, label = message,
      color = "#52616B", size = 4
    ) +
    ggplot2::labs(title = title) +
    ggplot2::xlim(0, 1) +
    ggplot2::ylim(0, 1) +
    ggplot2::theme_void() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", color = "#17324D")
    )
}

.bioana_read_expression_matrix <- function(expression) {
  if (is.character(expression) && length(expression) == 1L) {
    if (!file.exists(expression)) {
      stop("Expression file does not exist: ", expression)
    }
    expression <- utils::read.delim(
      expression,
      header = TRUE,
      sep = "\t",
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
  }

  if (is.data.frame(expression)) {
    id_candidates <- which(tolower(names(expression)) %in% c(
      "gene_id", "transcript_id", "feature_id", "id"
    ))
    if (length(id_candidates)) {
      ids <- as.character(expression[[id_candidates[1L]]])
      expression <- expression[, -id_candidates[1L], drop = FALSE]
    } else if (!is.null(rownames(expression))) {
      ids <- rownames(expression)
    } else {
      stop("Expression data frame needs a gene_id column or row names")
    }
    expression <- as.matrix(expression)
    rownames(expression) <- ids
  } else {
    expression <- as.matrix(expression)
  }

  if (is.null(rownames(expression)) || is.null(colnames(expression))) {
    stop("Expression matrix must have gene and sample names")
  }
  suppressWarnings(storage.mode(expression) <- "numeric")
  if (all(is.na(expression))) {
    stop("Expression matrix is not numeric")
  }
  expression
}

.bioana_results_table <- function(path, label) {
  if (!file.exists(path)) {
    stop(label, " does not exist: ", path)
  }
  result <- utils::read.delim(
    path,
    header = TRUE,
    sep = "\t",
    check.names = FALSE,
    stringsAsFactors = FALSE,
    colClasses = "character"
  )
  result[] <- lapply(
    result,
    utils::type.convert,
    as.is = TRUE,
    tryLogical = FALSE
  )
  result
}

.bioana_fit_eqtl <- function(y, genotype, covariate_matrix, interaction = NULL) {
  x <- cbind(covariate_matrix, genotype = genotype)
  if (!is.null(interaction)) {
    x <- cbind(x, genotype_by_condition = genotype * interaction)
  }
  complete <- is.finite(y) & is.finite(genotype) & apply(is.finite(x), 1L, all)
  x <- x[complete, , drop = FALSE]
  y <- y[complete]
  if (length(y) <= ncol(x) || stats::var(genotype[complete]) == 0) {
    return(c(beta = NA, se = NA, p_value = NA, gxe_beta = NA, gxe_p_value = NA))
  }

  cross_inverse <- tryCatch(
    solve(crossprod(x)),
    error = function(error) NULL
  )
  if (is.null(cross_inverse)) {
    return(c(beta = NA, se = NA, p_value = NA, gxe_beta = NA, gxe_p_value = NA))
  }
  coefficients <- as.numeric(cross_inverse %*% crossprod(x, y))
  names(coefficients) <- colnames(x)
  residuals <- y - as.numeric(x %*% coefficients)
  degrees_freedom <- length(y) - ncol(x)
  sigma_squared <- sum(residuals^2) / degrees_freedom
  standard_errors <- sqrt(diag(cross_inverse) * sigma_squared)
  names(standard_errors) <- colnames(x)
  statistics <- coefficients / standard_errors
  p_values <- 2 * stats::pt(abs(statistics), df = degrees_freedom, lower.tail = FALSE)

  c(
    beta = unname(coefficients["genotype"]),
    se = unname(standard_errors["genotype"]),
    p_value = unname(p_values["genotype"]),
    gxe_beta = if ("genotype_by_condition" %in% names(coefficients)) {
      unname(coefficients["genotype_by_condition"])
    } else NA_real_,
    gxe_p_value = if ("genotype_by_condition" %in% names(p_values)) {
      unname(p_values["genotype_by_condition"])
    } else NA_real_
  )
}

.bioana_read_sample_metadata <- function(sample_metadata, sample_ids) {
  if (is.null(sample_metadata)) {
    return(data.frame(sample = sample_ids, stringsAsFactors = FALSE))
  }
  if (is.character(sample_metadata) && length(sample_metadata) == 1L) {
    if (!file.exists(sample_metadata)) {
      stop("Sample metadata file does not exist: ", sample_metadata)
    }
    sample_metadata <- utils::read.delim(
      sample_metadata,
      header = TRUE,
      sep = "\t",
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
  }
  sample_metadata <- as.data.frame(sample_metadata, stringsAsFactors = FALSE)
  id_column <- names(sample_metadata)[
    tolower(names(sample_metadata)) %in% c("sample", "sample_id", "individual", "id")
  ][1L]
  if (!length(id_column) || is.na(id_column)) {
    stop("Sample metadata needs a sample or sample_id column")
  }
  names(sample_metadata)[names(sample_metadata) == id_column] <- "sample"
  if (anyDuplicated(sample_metadata$sample)) {
    stop("Sample metadata contains duplicate sample IDs")
  }
  missing_samples <- setdiff(sample_ids, sample_metadata$sample)
  if (length(missing_samples)) {
    stop(
      "Sample metadata is missing genotype samples: ",
      paste(utils::head(missing_samples, 5L), collapse = ", ")
    )
  }
  sample_metadata[match(sample_ids, sample_metadata$sample), , drop = FALSE]
}
