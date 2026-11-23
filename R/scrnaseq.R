# Single-cell analysis (moved from Simitall): clustering, cluster scoring, condition tests, cell-type eQTL benchmarks, results figure.

#' Cluster single cells and report the clusters
#'
#' Log-normalizes counts, keeps the most variable genes, runs PCA, and
#' clusters cells with k-means on the leading PCs. The number of clusters is
#' chosen by average silhouette width over a range of k.
#'
#' @param object single-cell object (from simitall::simulate_scrnaseq_experiment() or the same format) or its prefix.
#' @param n_pcs Number of principal components used for clustering.
#' @param n_hvg Number of highly variable genes; `NULL` keeps 20% of genes
#'   (between 50 and 2,000).
#' @param k_range Candidate cluster numbers.
#' @param seed Random seed.
#' @return A list with `clusters` (data frame of cell_id and cluster), `pcs`,
#'   and `chosen_k`.
#' @export
cluster_scrnaseq <- function(object, n_pcs = 20L, n_hvg = NULL, k_range = 2:10, seed = 1L) {
  x <- .bioana_scrna_counts(object)
  logn <- .bioana_scrna_lognorm(x$counts)
  # Keep a minority of genes, as standard pipelines do (~2,000 of 20,000):
  # this enriches for cell-type markers over genome-wide variation.
  if (is.null(n_hvg)) n_hvg <- min(2000L, max(50L, round(0.2 * nrow(logn))))
  vars <- apply(logn, 1L, stats::var)
  keep <- order(vars, decreasing = TRUE)[seq_len(min(n_hvg, sum(vars > 0)))]
  scaled <- scale(t(logn[keep, , drop = FALSE]))
  scaled[!is.finite(scaled)] <- 0
  n_pcs <- max(2L, min(as.integer(n_pcs), ncol(scaled) - 1L, nrow(scaled) - 1L))
  pcs <- stats::prcomp(scaled, center = FALSE, rank. = n_pcs)$x
  set.seed(seed)
  sub <- if (nrow(pcs) > 1000L) sample(nrow(pcs), 1000L) else seq_len(nrow(pcs))
  d <- stats::dist(pcs[sub, , drop = FALSE])
  fits <- lapply(k_range[k_range < nrow(pcs)], function(k) stats::kmeans(pcs, centers = k, nstart = 10, iter.max = 50))
  dm <- as.matrix(d)
  widths <- vapply(fits, function(fit) .bioana_mean_silhouette(fit$cluster[sub], dm), numeric(1))
  best <- if (all(is.na(widths))) 1L else which.max(widths)
  list(clusters = data.frame(cell_id = rownames(pcs) %||% x$cells$cell_id, cluster = paste0("cluster_", fits[[best]]$cluster),
                             stringsAsFactors = FALSE),
       pcs = pcs, chosen_k = k_range[best], silhouette = stats::setNames(widths, k_range[seq_along(widths)]))
}

#' Score clusters against the true cell types
#'
#' @param clusters Output of cluster_scrnaseq() (or its `clusters` table).
#' @param cell_metadata Cell metadata with `cell_id` and `cell_type`.
#' @return One-row data frame with the adjusted Rand index and NMI.
#' @export
score_scrnaseq_clusters <- function(clusters, cell_metadata) {
  if (is.list(clusters) && !is.data.frame(clusters)) clusters <- clusters$clusters
  truth <- cell_metadata$cell_type[match(clusters$cell_id, cell_metadata$cell_id)]
  data.frame(analysis = "clustering", true_cell_types = length(unique(truth)),
             clusters_found = length(unique(clusters$cluster)),
             adjusted_rand_index = .bioana_adjusted_rand(clusters$cluster, truth),
             nmi = .bioana_nmi(clusters$cluster, truth))
}

#' Test genes for condition differences in single-cell data
#'
#' Pseudobulk methods sum counts per donor and test donors as replicates
#' (edgeR, DESeq2, or a t-test on log-CPM). The per-cell Wilcoxon test treats
#' every cell as a replicate, which ignores donor variation and inflates false
#' positives; it is included as a contrast.
#'
#' @param object single-cell object (from simitall::simulate_scrnaseq_experiment() or the same format) or its prefix.
#' @param method `"pseudobulk_deseq2"`, `"pseudobulk_edger"`,
#'   `"pseudobulk_ttest"`, or `"cell_wilcoxon"`.
#' @param fdr BH-adjusted threshold.
#' @return Data frame with gene_id, log2_fold_change, p_value, q_value, significant.
#' @export
test_scrnaseq_conditions <- function(object, method = c("pseudobulk_deseq2", "pseudobulk_edger", "pseudobulk_ttest", "cell_wilcoxon"),
                                     fdr = 0.05) {
  method <- match.arg(method)
  x <- .bioana_scrna_counts(object)
  cells <- x$cells[match(colnames(x$counts), x$cells$cell_id), , drop = FALSE]
  levels <- unique(cells$condition)
  if (length(levels) < 2L) stop("Condition testing needs two conditions")
  levels <- levels[1:2]
  if (identical(method, "cell_wilcoxon")) {
    logn <- .bioana_scrna_lognorm(x$counts)
    a <- cells$condition == levels[1L]; b <- cells$condition == levels[2L]
    stat <- t(apply(logn, 1L, function(v) {
      p <- if (stats::var(v[a | b]) == 0) 1 else suppressWarnings(stats::wilcox.test(v[b], v[a])$p.value)
      c(log2(mean(expm1(v[b])) + 1) - log2(mean(expm1(v[a])) + 1), p)
    }))
    result <- data.frame(gene_id = rownames(logn), log2_fold_change = stat[, 1L], p_value = stat[, 2L], stringsAsFactors = FALSE)
  } else {
    donors <- unique(cells$donor)
    pb <- vapply(donors, function(d) Matrix::rowSums(x$counts[, cells$donor == d, drop = FALSE]), numeric(nrow(x$counts)))
    rownames(pb) <- rownames(x$counts)
    group <- factor(cells$condition[match(donors, cells$donor)], levels = levels)
    keep <- !is.na(group); pb <- pb[, keep, drop = FALSE]; group <- droplevels(group[keep])
    if (min(table(group)) < 2L) stop("Pseudobulk testing needs at least two donors per condition")
    pb <- round(pb)
    result <- switch(method,
      pseudobulk_ttest = {
        logcpm <- log2(t(t(pb + 0.5) / (colSums(pb) + 1)) * 1e6)
        a <- group == levels(group)[1L]
        s <- t(apply(logcpm, 1L, function(v) {
          if (stats::sd(v[a]) + stats::sd(v[!a]) == 0) return(c(mean(v[!a]) - mean(v[a]), 1))
          c(mean(v[!a]) - mean(v[a]), stats::t.test(v[!a], v[a])$p.value)
        }))
        data.frame(log2_fold_change = s[, 1L], p_value = s[, 2L])
      },
      pseudobulk_edger = {
        if (!requireNamespace("edgeR", quietly = TRUE)) stop("method 'pseudobulk_edger' needs the edgeR package")
        y <- edgeR::calcNormFactors(edgeR::DGEList(pb, group = group))
        y <- edgeR::estimateDisp(y, stats::model.matrix(~group))
        et <- edgeR::exactTest(y)$table
        data.frame(log2_fold_change = et$logFC, p_value = et$PValue)
      },
      pseudobulk_deseq2 = {
        if (!requireNamespace("DESeq2", quietly = TRUE)) stop("method 'pseudobulk_deseq2' needs the DESeq2 package")
        dds <- DESeq2::DESeqDataSetFromMatrix(pb, data.frame(group = group), ~group)
        res <- as.data.frame(DESeq2::results(suppressMessages(DESeq2::DESeq(dds, quiet = TRUE))))
        data.frame(log2_fold_change = res$log2FoldChange, p_value = ifelse(is.na(res$pvalue), 1, res$pvalue))
      })
    result <- cbind(data.frame(gene_id = rownames(pb), stringsAsFactors = FALSE), result)
  }
  result$q_value <- stats::p.adjust(result$p_value, "BH")
  result$significant <- result$q_value <= fdr
  result$method <- method
  result
}

#' Score condition DE calls against the simulated condition truth
#'
#' @param de Output of test_scrnaseq_conditions().
#' @param condition_truth Condition truth table or path (gene_id, effect_log2).
#' @param min_effect Smallest absolute true log2 effect counted as a real
#'   difference (tiny simulated effects are practically undetectable).
#' @return One-row data frame with recall, precision, F1, and false positives.
#' @export
score_scrnaseq_conditions <- function(de, condition_truth, min_effect = 0.25) {
  if (is.character(condition_truth)) condition_truth <- utils::read.delim(condition_truth, stringsAsFactors = FALSE)
  real <- unique(condition_truth$gene_id[abs(condition_truth$effect_log2) >= min_effect])
  tiny <- setdiff(unique(condition_truth$gene_id), real)
  called <- setdiff(de$gene_id[de$significant], tiny)
  tp <- length(intersect(called, real)); fp <- length(setdiff(called, real))
  recall <- if (length(real)) tp / length(real) else NA_real_
  precision <- if (length(called)) tp / length(called) else NA_real_
  data.frame(analysis = "condition DE", method = de$method[1L], true_de_genes = length(real),
             called = length(called), recovered = tp, false_positives = fp, recall = recall, precision = precision,
             f1 = if (is.na(recall) || is.na(precision) || recall + precision == 0) 0 else 2 * recall * precision / (recall + precision),
             note = paste0("genes with |true effect| < ", min_effect, " are excluded from both counts"))
}

#' Plot donor-aware single-cell RNA-seq simulation results
#'
#' Build a six-panel publication figure showing cells per donor, single-cell
#' PCA, marker-program recovery, pseudobulk structure, cell-type eQTL effect
#' recovery, and benchmark precision/recall. Every panel's source data are
#' exported beside PDF and PNG figures.
#'
#' @param scrna_prefix Prefix passed to `simitall::simulate_scrnaseq_from_gwas()`.
#' @param eqtl_prefix Optional prefix passed to [benchmark_celltype_eqtls()].
#'   When omitted, the two eQTL panels explain that no benchmark was supplied.
#' @param out_prefix Prefix for figures, summary TSV, and source-data CSVs.
#' @param top_variable_genes Number of variable genes used for PCA.
#' @param markers_per_cell_type Number of marker genes shown per cell type.
#' @param max_plot_cells Maximum cells displayed in the single-cell PCA.
#' @param width Figure width in inches.
#' @param height Figure height in inches.
#' @param dpi PNG resolution.
#' @param seed Seed used for plot-cell subsampling.
#' @param individual_panels_dir Optional directory for standalone PNG/PDF
#'   panels. Defaults to a `figures` folder beside `out_prefix`.
#'
#' @return Invisibly returns the combined figure, panel objects, source data,
#'   and output paths.
#' @examples
#' \dontrun{
#' plot_scrnaseq_results(
#'   scrna_prefix = "results/scrna/demo",
#'   eqtl_prefix = "results/scrna/demo_celltype_eqtl",
#'   out_prefix = "results/figures/scrna_demo"
#' )
#' }
#' @export
plot_scrnaseq_results <- function(
    scrna_prefix,
    eqtl_prefix = NULL,
    out_prefix,
    top_variable_genes = 500L,
    markers_per_cell_type = 4L,
    max_plot_cells = 4000L,
    width = 15,
    height = 10,
    dpi = 300,
    seed = 1L,
    individual_panels_dir = file.path(dirname(out_prefix), "figures")) {
  if (!requireNamespace("ggplot2", quietly = TRUE) ||
      !requireNamespace("patchwork", quietly = TRUE)) {
    stop(
      "ggplot2 and patchwork are required. Install them with ",
      "install.packages(c('ggplot2', 'patchwork'))."
    )
  }
  if (top_variable_genes < 2L || markers_per_cell_type < 1L ||
      max_plot_cells < 10L) {
    stop("PCA, marker, and cell plotting limits are too small")
  }
  set.seed(seed)
  object <- .bioana_read_scrna_object(scrna_prefix)
  counts <- object$counts
  cells <- object$cell_metadata
  genes <- object$gene_metadata
  marker_truth <- object$marker_truth
  if (!all(colnames(counts) == cells$cell_id)) {
    cells <- cells[match(colnames(counts), cells$cell_id), , drop = FALSE]
  }

  palette <- c(
    navy = "#153243",
    blue = "#2D6A8A",
    cyan = "#4B9DA9",
    green = "#6A994E",
    sand = "#E5B25D",
    orange = "#D9783D",
    coral = "#C8553D",
    gray = "#A8B1B7"
  )
  cell_type_colors <- stats::setNames(
    rep(
      c(
        palette[["blue"]], palette[["coral"]], palette[["green"]],
        palette[["sand"]], palette[["cyan"]], palette[["orange"]]
      ),
      length.out = length(unique(cells$cell_type))
    ),
    unique(cells$cell_type)
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

  cell_table <- as.data.frame(
    table(cells$donor, cells$cell_type),
    stringsAsFactors = FALSE
  )
  names(cell_table) <- c("donor", "cell_type", "n_cells")
  panel_a_data <- cell_table
  panel_a <- ggplot2::ggplot(
    panel_a_data,
    ggplot2::aes(x = cell_type, y = n_cells, fill = cell_type)
  )
  if (length(unique(cells$donor)) >= 2L) {
    panel_a <- panel_a +
      ggplot2::geom_violin(trim = FALSE, alpha = 0.72, color = NA) +
      ggplot2::geom_boxplot(
        width = 0.16, outlier.shape = NA, fill = "white"
      ) +
      ggplot2::geom_jitter(width = 0.08, size = 0.7, alpha = 0.35)
  } else {
    panel_a <- panel_a +
      ggplot2::geom_col(width = 0.7, alpha = 0.8) +
      ggplot2::geom_text(
        ggplot2::aes(label = n_cells), vjust = -0.35, size = 3
      )
  }
  donor_label <- if (identical(object$experiment_type, "standalone")) {
    "biological samples"
  } else {
    "GWAS donors"
  }
  panel_a <- panel_a +
    ggplot2::scale_fill_manual(values = cell_type_colors) +
    ggplot2::labs(
      title = "A. Cells per donor",
      subtitle = paste(length(unique(cells$donor)), donor_label),
      x = NULL,
      y = "Cells"
    ) +
    result_theme +
    ggplot2::theme(
      legend.position = "none",
      axis.text.x = ggplot2::element_text(angle = 25, hjust = 1)
    )

  plot_cell_index <- seq_len(ncol(counts))
  if (length(plot_cell_index) > max_plot_cells) {
    plot_cell_index <- sort(sample(plot_cell_index, max_plot_cells))
  }
  gene_mean <- Matrix::rowMeans(counts)
  gene_mean_square <- Matrix::rowMeans(counts^2)
  gene_variance <- pmax(0, gene_mean_square - gene_mean^2)
  variable_gene_index <- which(is.finite(gene_variance) & gene_variance > 0)
  if (length(variable_gene_index) < 2L) {
    stop("At least two variable genes are required for single-cell PCA")
  }
  variable_gene_index <- variable_gene_index[
    order(gene_variance[variable_gene_index], decreasing = TRUE)
  ]
  variable_gene_index <- variable_gene_index[
    seq_len(min(top_variable_genes, length(variable_gene_index)))
  ]
  selected_counts <- as.matrix(
    counts[variable_gene_index, plot_cell_index, drop = FALSE]
  )
  selected_libraries <- pmax(1, Matrix::colSums(counts[, plot_cell_index, drop = FALSE]))
  selected_logcounts <- log2(
    sweep(selected_counts, 2L, selected_libraries, "/") * 1e4 + 1
  )
  pca <- stats::prcomp(t(selected_logcounts), center = TRUE, scale. = FALSE)
  explained <- 100 * pca$sdev^2 / sum(pca$sdev^2)
  panel_b_data <- data.frame(
    cell_id = rownames(pca$x),
    PC1 = pca$x[, 1L],
    PC2 = pca$x[, 2L],
    cell_type = cells$cell_type[
      match(rownames(pca$x), cells$cell_id)
    ],
    donor = cells$donor[match(rownames(pca$x), cells$cell_id)],
    condition = cells$condition[match(rownames(pca$x), cells$cell_id)],
    stringsAsFactors = FALSE
  )
  panel_b <- ggplot2::ggplot(
    panel_b_data,
    ggplot2::aes(x = PC1, y = PC2, color = cell_type)
  ) +
    ggplot2::geom_point(size = 0.65, alpha = 0.55) +
    ggplot2::scale_color_manual(values = cell_type_colors) +
    ggplot2::labs(
      title = "B. Single-cell expression structure",
      subtitle = paste(nrow(panel_b_data), "cells;", length(variable_gene_index),
                       "variable genes"),
      x = sprintf("PC1 (%.1f%%)", explained[1L]),
      y = sprintf("PC2 (%.1f%%)", explained[2L]),
      color = "Cell type"
    ) +
    result_theme

  marker_truth <- marker_truth[order(
    marker_truth$cell_type,
    -marker_truth$effect_log2
  ), , drop = FALSE]
  selected_marker_rows <- unlist(lapply(
    split(seq_len(nrow(marker_truth)), marker_truth$cell_type),
    function(index) utils::head(index, markers_per_cell_type)
  ))
  marker_subset <- marker_truth[selected_marker_rows, , drop = FALSE]
  marker_index <- match(marker_subset$gene_id, rownames(counts))
  marker_counts <- as.matrix(counts[marker_index, , drop = FALSE])
  marker_logcounts <- log2(
    sweep(marker_counts, 2L, pmax(1, Matrix::colSums(counts)), "/") * 1e4 + 1
  )
  type_levels <- unique(cells$cell_type)
  marker_means <- sapply(type_levels, function(cell_type) {
    rowMeans(marker_logcounts[, cells$cell_type == cell_type, drop = FALSE])
  })
  if (is.null(dim(marker_means))) {
    marker_means <- matrix(marker_means, ncol = 1L)
    colnames(marker_means) <- type_levels
  }
  rownames(marker_means) <- marker_subset$gene_id
  panel_c_data <- as.data.frame(
    as.table(marker_means),
    stringsAsFactors = FALSE,
    responseName = "mean_log_expression"
  )
  names(panel_c_data)[1:2] <- c("gene_id", "observed_cell_type")
  panel_c_data$marker_cell_type <- marker_subset$cell_type[
    match(panel_c_data$gene_id, marker_subset$gene_id)
  ]
  panel_c_data$gene_label <- paste0(
    panel_c_data$gene_id,
    " [",
    panel_c_data$marker_cell_type,
    "]"
  )
  panel_c <- ggplot2::ggplot(
    panel_c_data,
    ggplot2::aes(
      x = observed_cell_type,
      y = gene_label,
      fill = mean_log_expression
    )
  ) +
    ggplot2::geom_tile(color = "white", linewidth = 0.25) +
    ggplot2::scale_fill_gradientn(colors = c("#F5F2E8", palette[["sand"]], palette[["coral"]])) +
    ggplot2::labs(
      title = "C. Cell-type marker recovery",
      subtitle = "Labels show the programmed marker cell type",
      x = "Observed cell type",
      y = NULL,
      fill = "Mean log2\nCP10K"
    ) +
    result_theme +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 25, hjust = 1),
      axis.text.y = ggplot2::element_text(size = 6.5)
    )

  pseudobulk_expression <- .bioana_read_expression_matrix(
    paste0(scrna_prefix, ".pseudobulk.expression.tsv")
  )
  pseudobulk_metadata <- .bioana_results_table(
    paste0(scrna_prefix, ".pseudobulk.metadata.tsv"),
    "Pseudobulk metadata"
  )
  pseudobulk_metadata <- pseudobulk_metadata[
    match(colnames(pseudobulk_expression), pseudobulk_metadata$pseudobulk_id),
    , drop = FALSE
  ]
  pseudobulk_variance <- apply(pseudobulk_expression, 1L, stats::var)
  pseudobulk_genes <- order(pseudobulk_variance, decreasing = TRUE)
  pseudobulk_genes <- pseudobulk_genes[
    seq_len(min(top_variable_genes, length(pseudobulk_genes)))
  ]
  pseudobulk_pca <- stats::prcomp(
    t(pseudobulk_expression[pseudobulk_genes, , drop = FALSE]),
    center = TRUE,
    scale. = FALSE
  )
  pseudobulk_explained <- 100 * pseudobulk_pca$sdev^2 /
    sum(pseudobulk_pca$sdev^2)
  panel_d_data <- data.frame(
    pseudobulk_id = rownames(pseudobulk_pca$x),
    PC1 = pseudobulk_pca$x[, 1L],
    PC2 = pseudobulk_pca$x[, 2L],
    donor = pseudobulk_metadata$donor,
    cell_type = pseudobulk_metadata$cell_type,
    condition = pseudobulk_metadata$condition,
    stringsAsFactors = FALSE
  )
  panel_d <- ggplot2::ggplot(
    panel_d_data,
    ggplot2::aes(
      x = PC1,
      y = PC2,
      color = cell_type,
      shape = condition
    )
  ) +
    ggplot2::geom_point(size = 1.8, alpha = 0.75) +
    ggplot2::scale_color_manual(values = cell_type_colors) +
    ggplot2::labs(
      title = "D. Donor-level pseudobulk structure",
      subtitle = "Donors remain the biological replication unit",
      x = sprintf("PC1 (%.1f%%)", pseudobulk_explained[1L]),
      y = sprintf("PC2 (%.1f%%)", pseudobulk_explained[2L]),
      color = "Cell type",
      shape = "Condition"
    ) +
    result_theme +
    ggplot2::guides(color = "none")

  if (!is.null(eqtl_prefix)) {
    associations <- .bioana_results_table(
      paste0(eqtl_prefix, ".celltype_eqtl_results.tsv"),
      "Cell-type eQTL results"
    )
    metrics <- .bioana_results_table(
      paste0(eqtl_prefix, ".celltype_eqtl_metrics.tsv"),
      "Cell-type eQTL metrics"
    )
    truth <- object$eqtl_truth
    truth_rows <- do.call(rbind, lapply(unique(cells$cell_type), function(cell_type) {
      active <- truth$scope == "shared" |
        (truth$scope == "cell_type_specific" & truth$cell_type == cell_type)
      data.frame(
        truth[active, c("gene_id", "variant_id", "effect_log2", "scope")],
        cell_type = cell_type,
        stringsAsFactors = FALSE
      )
    }))
    effect_data <- merge(
      truth_rows,
      associations[, c("gene_id", "variant_id", "cell_type", "beta", "q_value")],
      by = c("gene_id", "variant_id", "cell_type"),
      all.x = TRUE
    )
    complete <- is.finite(effect_data$effect_log2) & is.finite(effect_data$beta)
    effect_correlation <- if (sum(complete) >= 3L) {
      stats::cor(effect_data$effect_log2[complete], effect_data$beta[complete])
    } else NA_real_
    panel_e_data <- effect_data
    panel_e <- ggplot2::ggplot(
      panel_e_data,
      ggplot2::aes(
        x = effect_log2,
        y = beta,
        color = cell_type,
        shape = scope
      )
    ) +
      ggplot2::geom_abline(slope = 1, intercept = 0, color = palette[["gray"]]) +
      ggplot2::geom_point(size = 1.5, alpha = 0.7, na.rm = TRUE) +
      ggplot2::scale_color_manual(values = cell_type_colors) +
      ggplot2::labs(
        title = "E. Cell-type eQTL effect recovery",
        subtitle = sprintf("Truth versus estimate; r = %.2f", effect_correlation),
        x = "True effect (log2)",
        y = "Estimated effect",
        color = "Cell type",
        shape = "Scope"
      ) +
      result_theme +
      ggplot2::guides(color = "none")

    metric_names <- c("precision", "recall", "empirical_fdr")
    metric_labels <- c(
      precision = "Precision",
      recall = "Recall",
      empirical_fdr = "Empirical FDR"
    )
    panel_f_data <- do.call(rbind, lapply(metric_names, function(metric) {
      data.frame(
        cell_type = metrics$cell_type,
        metric = metric_labels[[metric]],
        value = metrics[[metric]],
        stringsAsFactors = FALSE
      )
    }))
    panel_f <- ggplot2::ggplot(
      panel_f_data,
      ggplot2::aes(x = cell_type, y = value, fill = metric)
    ) +
      ggplot2::geom_col(position = "dodge", width = 0.72, na.rm = TRUE) +
      ggplot2::scale_fill_manual(values = c(
        Precision = palette[["blue"]],
        Recall = palette[["green"]],
        `Empirical FDR` = palette[["coral"]]
      )) +
      ggplot2::scale_y_continuous(limits = c(0, 1)) +
      ggplot2::labs(
        title = "F. Donor-level eQTL benchmark",
        x = NULL,
        y = "Metric value",
        fill = NULL
      ) +
      result_theme +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 25, hjust = 1))
  } else {
    panel_e_data <- data.frame(message = "Run benchmark_celltype_eqtls()")
    panel_f_data <- data.frame(message = "Run benchmark_celltype_eqtls()")
    panel_e <- .bioana_empty_panel(
      "E. Cell-type eQTL effect recovery",
      "Run benchmark_celltype_eqtls() to add this panel"
    )
    panel_f <- .bioana_empty_panel(
      "F. Donor-level eQTL benchmark",
      "Run benchmark_celltype_eqtls() to add this panel"
    )
  }

  figure_title <- if (identical(object$experiment_type, "standalone")) {
    "Standalone single-cell RNA-seq experiment simulation"
  } else {
    "Donor-aware single-cell RNA-seq and cell-type eQTL simulation"
  }
  combined <- (panel_a | panel_b | panel_c) /
    (panel_d | panel_e | panel_f) +
    patchwork::plot_annotation(
      title = figure_title,
      subtitle = paste0(
        length(unique(cells$donor)), " donors | ", nrow(cells), " cells | ",
        nrow(genes), " genes | ", object$backend, " backend"
      ),
      theme = ggplot2::theme(
        plot.title = ggplot2::element_text(
          face = "bold", size = 16, color = palette[["navy"]]
        ),
        plot.subtitle = ggplot2::element_text(size = 10, color = "#52616B")
      )
    )
  source_data <- list(
    panel_a_cells_per_donor = panel_a_data,
    panel_b_single_cell_pca = panel_b_data,
    panel_c_marker_expression = panel_c_data,
    panel_d_pseudobulk_pca = panel_d_data,
    panel_e_eqtl_effects = panel_e_data,
    panel_f_eqtl_metrics = panel_f_data
  )
  source_directory <- paste0(out_prefix, "_source_data")
  dir.create(source_directory, recursive = TRUE, showWarnings = FALSE)
  for (name in names(source_data)) {
    utils::write.csv(
      source_data[[name]],
      file.path(source_directory, paste0(name, ".csv")),
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
  ggplot2::ggsave(
    paths$png, combined, width = width, height = height, dpi = dpi
  )
  summary <- data.frame(
    donors = length(unique(cells$donor)),
    cells = nrow(cells),
    genes = nrow(genes),
    cell_types = length(unique(cells$cell_type)),
    median_umis = stats::median(cells$observed_library_size),
    median_detected_genes = stats::median(cells$detected_genes),
    doublet_fraction = mean(cells$is_doublet),
    zero_fraction = 1 - Matrix::nnzero(counts) / length(counts),
    backend = object$backend,
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
      c("scrnaseq_cells_per_donor", "scrnaseq_cell_pca",
        "scrnaseq_marker_expression", "scrnaseq_pseudobulk_pca",
        "scrnaseq_eqtl_effects", "scrnaseq_eqtl_metrics"),
      width = 7, height = 5, dpi = dpi
    )
  }
  paths$individual_panels <- individual_panels
  invisible(list(
    figure = combined,
    panels = list(
      panel_a, panel_b, panel_c, panel_d, panel_e, panel_f
    ),
    source_data = source_data,
    paths = paths
  ))
}

#' Benchmark cell-type-specific eQTL recovery
#'
#' Test eQTLs separately in each cell type using donor-level pseudobulk
#' expression. Shared eQTLs are causal in every cell type; cell-type-specific
#' eQTLs are causal only in their assigned cell type. This avoids treating
#' individual cells as independent biological replicates.
#'
#' @param genotype_file Genotype TSV or VCF/VCF.GZ.
#' @param pseudobulk_expression Gene-by-pseudobulk log2-expression matrix or
#'   TSV from `simitall::aggregate_scrnaseq_pseudobulk()`.
#' @param pseudobulk_metadata Metadata data frame or TSV with `pseudobulk_id`,
#'   `donor`, and `cell_type` columns.
#' @param truth_file Cell-type eQTL truth table from
#'   `simitall::simulate_scrnaseq_from_gwas()`.
#' @param out_prefix Output prefix.
#' @param cis_window_bp Cis candidate window.
#' @param covariates Donor metadata columns included in regression models.
#' @param test_gxe Test genotype-by-condition interactions.
#' @param fdr_threshold BH-adjusted significance threshold.
#' @param min_donors Minimum donors required per cell type.
#' @param max_tests Maximum association tests per cell type.
#'
#' @return Invisibly returns combined result and metric paths and tables.
#' @examples
#' \dontrun{
#' benchmark_celltype_eqtls(
#'   genotype_file = "results/gwas/demo.geno.tsv",
#'   pseudobulk_expression = sc$pseudobulk_expression,
#'   pseudobulk_metadata = sc$pseudobulk_metadata,
#'   truth_file = sc$eqtl_truth,
#'   out_prefix = "results/scrna/demo_celltype_eqtl"
#' )
#' }
#' @export
benchmark_celltype_eqtls <- function(
    genotype_file,
    pseudobulk_expression,
    pseudobulk_metadata,
    truth_file,
    out_prefix,
    cis_window_bp = 1e6,
    covariates = c("condition", "batch", "pop"),
    test_gxe = TRUE,
    fdr_threshold = 0.05,
    min_donors = 8L,
    max_tests = 1e6) {
  expression <- .bioana_read_expression_matrix(pseudobulk_expression)
  metadata <- .bioana_read_table_or_data(
    pseudobulk_metadata, "Pseudobulk metadata"
  )
  truth <- .bioana_read_table_or_data(truth_file, "eQTL truth")
  required_metadata <- c("pseudobulk_id", "donor", "cell_type")
  if (!all(required_metadata %in% names(metadata))) {
    stop(
      "Pseudobulk metadata is missing: ",
      paste(setdiff(required_metadata, names(metadata)), collapse = ", ")
    )
  }
  required_truth <- c(
    "gene_id", "variant_id", "seqname", "gene_tss", "scope", "cell_type"
  )
  if (!all(required_truth %in% names(truth))) {
    stop(
      "eQTL truth is missing: ",
      paste(setdiff(required_truth, names(truth)), collapse = ", ")
    )
  }
  if (fdr_threshold <= 0 || fdr_threshold > 1) {
    stop("fdr_threshold must be greater than zero and at most one")
  }
  cell_types <- unique(metadata$cell_type)
  results_list <- list()
  metrics_list <- list()
  result_number <- 0L

  for (cell_type in cell_types) {
    cell_metadata <- metadata[metadata$cell_type == cell_type, , drop = FALSE]
    shared_libraries <- intersect(cell_metadata$pseudobulk_id, colnames(expression))
    cell_metadata <- cell_metadata[
      match(shared_libraries, cell_metadata$pseudobulk_id), , drop = FALSE
    ]
    if (length(unique(cell_metadata$donor)) < min_donors) {
      warning(
        "Skipping ", cell_type, ": fewer than ", min_donors, " donors",
        call. = FALSE
      )
      next
    }
    truth_for_type <- truth[
      truth$scope == "shared" |
        (truth$scope == "cell_type_specific" & truth$cell_type == cell_type),
      , drop = FALSE
    ]
    if (!nrow(truth_for_type)) next
    cell_expression <- expression[, shared_libraries, drop = FALSE]
    colnames(cell_expression) <- cell_metadata$donor
    analysis_metadata <- cell_metadata
    analysis_metadata$sample <- analysis_metadata$donor
    analysis_metadata <- analysis_metadata[
      !duplicated(analysis_metadata$sample), , drop = FALSE
    ]
    cell_expression <- cell_expression[
      , analysis_metadata$sample, drop = FALSE
    ]
    usable_covariates <- intersect(covariates, names(analysis_metadata))
    usable_covariates <- usable_covariates[vapply(
      analysis_metadata[, usable_covariates, drop = FALSE],
      function(value) length(unique(value[!is.na(value)])) > 1L,
      logical(1L)
    )]
    temporary_prefix <- tempfile(paste0("bioana-", cell_type, "-"))
    fit <- benchmark_eqtl(
      genotype_file = genotype_file,
      expression_file = cell_expression,
      out_prefix = temporary_prefix,
      sample_metadata = analysis_metadata,
      truth_file = truth_for_type,
      pair_mode = "truth_neighborhood",
      cis_window_bp = cis_window_bp,
      covariates = usable_covariates,
      test_gxe = test_gxe,
      expression_scale = "log2",
      fdr_threshold = fdr_threshold,
      max_tests = max_tests
    )
    unlink(paste0(temporary_prefix, c(".eqtl_results.tsv", ".eqtl_metrics.tsv")))
    result_number <- result_number + 1L
    fit$results$cell_type <- cell_type
    fit$metrics$cell_type <- cell_type
    results_list[[result_number]] <- fit$results
    metrics_list[[result_number]] <- fit$metrics
  }
  if (!length(results_list)) {
    stop("No cell type had enough donors and causal truth for benchmarking")
  }
  results <- do.call(rbind, results_list)
  metrics <- do.call(rbind, metrics_list)
  overall <- data.frame(
    tests = sum(metrics$tests),
    truth_pairs_in_tests = sum(metrics$truth_pairs_in_tests),
    discoveries = sum(metrics$discoveries),
    true_positives = sum(metrics$true_positives),
    false_positives = sum(metrics$false_positives),
    false_negatives = sum(metrics$false_negatives),
    precision = if (sum(metrics$discoveries)) {
      sum(metrics$true_positives) / sum(metrics$discoveries)
    } else NA_real_,
    recall = if (sum(metrics$truth_pairs_in_tests)) {
      sum(metrics$true_positives) / sum(metrics$truth_pairs_in_tests)
    } else NA_real_,
    empirical_fdr = if (sum(metrics$discoveries)) {
      sum(metrics$false_positives) / sum(metrics$discoveries)
    } else NA_real_,
    fdr_threshold = fdr_threshold,
    cell_type = "Overall",
    stringsAsFactors = FALSE
  )
  metrics <- rbind(metrics, overall)

  dir.create(dirname(out_prefix), recursive = TRUE, showWarnings = FALSE)
  paths <- list(
    results = paste0(out_prefix, ".celltype_eqtl_results.tsv"),
    metrics = paste0(out_prefix, ".celltype_eqtl_metrics.tsv")
  )
  utils::write.table(
    results, paths$results, sep = "\t", row.names = FALSE, quote = FALSE
  )
  utils::write.table(
    metrics, paths$metrics, sep = "\t", row.names = FALSE, quote = FALSE
  )
  message(
    "Cell-type eQTL benchmark complete: ", nrow(results), " tests across ",
    length(results_list), " cell types."
  )
  invisible(list(
    results_path = paths$results,
    metrics_path = paths$metrics,
    results = results,
    metrics = metrics
  ))
}

# Mean silhouette width from a full distance matrix (no extra packages).
.bioana_mean_silhouette <- function(labels, dm) {
  groups <- unique(labels)
  if (length(groups) < 2L) return(NA_real_)
  widths <- vapply(seq_along(labels), function(i) {
    same <- labels == labels[i]; same[i] <- FALSE
    if (!any(same)) return(0)
    a <- mean(dm[i, same])
    b <- min(vapply(setdiff(groups, labels[i]), function(g) mean(dm[i, labels == g]), numeric(1)))
    (b - a) / max(a, b)
  }, numeric(1))
  mean(widths)
}

.bioana_scrna_counts <- function(object) {
  if (is.character(object)) object <- .bioana_read_scrna_object(object)
  counts <- object$counts
  if (is.null(rownames(counts))) rownames(counts) <- object$gene_metadata$gene_id
  if (is.null(colnames(counts))) colnames(counts) <- object$cell_metadata$cell_id
  list(counts = counts, cells = object$cell_metadata)
}

.bioana_scrna_lognorm <- function(counts) {
  lib <- Matrix::colSums(counts)
  log1p(t(t(as.matrix(counts)) / pmax(lib, 1)) * 1e4)
}

.bioana_adjusted_rand <- function(a, b) {
  tab <- table(a, b)
  comb <- function(n) n * (n - 1) / 2
  sum_ij <- sum(comb(tab)); sum_a <- sum(comb(rowSums(tab))); sum_b <- sum(comb(colSums(tab)))
  expected <- sum_a * sum_b / comb(sum(tab))
  max_index <- (sum_a + sum_b) / 2
  if (max_index == expected) return(1)
  (sum_ij - expected) / (max_index - expected)
}

.bioana_nmi <- function(a, b) {
  tab <- table(a, b) / length(a)
  pa <- rowSums(tab); pb <- colSums(tab)
  mi <- sum(tab[tab > 0] * log(tab[tab > 0] / outer(pa, pb)[tab > 0]))
  h <- function(p) -sum(p[p > 0] * log(p[p > 0]))
  if (h(pa) + h(pb) == 0) return(1)
  2 * mi / (h(pa) + h(pb))
}

.bioana_read_scrna_object <- function(prefix) {
  path <- if (grepl("\\.scrnaseq\\.rds$", prefix)) {
    prefix
  } else {
    paste0(prefix, ".scrnaseq.rds")
  }
  if (!file.exists(path)) {
    stop("Single-cell result object does not exist: ", path)
  }
  object <- readRDS(path)
  required <- c("counts", "cell_metadata", "gene_metadata")
  if (!all(required %in% names(object))) {
    stop("Single-cell result object is missing required components")
  }
  object
}

.bioana_read_table_or_data <- function(value, label) {
  if (is.character(value) && length(value) == 1L) {
    if (!file.exists(value)) stop(label, " does not exist: ", value)
    value <- utils::read.delim(
      value,
      header = TRUE,
      sep = "\t",
      check.names = FALSE,
      stringsAsFactors = FALSE,
      colClasses = "character"
    )
    value[] <- lapply(
      value,
      utils::type.convert,
      as.is = TRUE,
      tryLogical = FALSE
    )
  }
  as.data.frame(value, stringsAsFactors = FALSE)
}
