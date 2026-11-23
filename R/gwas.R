# GWAS on line genotypes (moved from Simitall): mixed-model association, binary traits, truth benchmarks.

#' Run a mixed-model GWAS with rrBLUP
#'
#' Test marker-trait associations with the mixed model implemented by
#' [rrBLUP::GWAS()]. Genotypes can be supplied as a VCF or as the marker-by-
#' sample TSV written by `simitall::simulate_gwas_cohort()`. Samples are aligned by ID,
#' population structure can be represented by principal components or fixed
#' effects, and a realized additive relationship matrix controls relatedness.
#'
#' @param genotype_file VCF, VCF.GZ, or genotype TSV path.
#' @param phenotype Data frame or phenotype TSV containing one row per sample.
#' @param out_prefix Prefix for `.gwas.tsv` and `.gwas_summary.json` outputs.
#' @param trait Numeric phenotype column to analyze.
#' @param sample_id Sample-ID column in `phenotype`.
#' @param fixed_effects Optional phenotype columns treated as categorical fixed
#'   effects by `rrBLUP`.
#' @param kinship Optional sample-by-sample kinship matrix. When `NULL`, an
#'   additive relationship matrix is estimated with [rrBLUP::A.mat()].
#' @param n_pcs Number of marker-derived principal components included as fixed
#'   effects.
#' @param min_maf Minimum minor-allele frequency passed to [rrBLUP::GWAS()].
#' @param n_cores Number of CPU cores used by `rrBLUP` where supported.
#' @param p3d Estimate variance components once before scanning markers. This is
#'   faster than marker-specific REML but can be less significant.
#' @param fdr_threshold BH-adjusted p-value threshold used to label discoveries.
#' @param backend `"auto"` (rrBLUP when installed, else native), `"rrBLUP"`,
#'   `"GAPIT"` (compressed MLM, kinship + PCs; returns allele effects), or
#'   `"native"`: the same model (kinship random effect, PCs and `fixed_effects`
#'   as fixed covariates, variance components estimated once by REML, as
#'   EMMAX / rrBLUP with `P3D = TRUE`), then GLS per marker. The native scan
#'   also returns each marker's effect (`beta`, `se`, per alternative allele);
#'   rrBLUP returns only the score. `p3d = FALSE` needs rrBLUP.
#'
#' @return Invisibly returns a list containing the GWAS table, kinship matrix,
#'   summary, and output paths.
#' @examples
#' \dontrun{
#' result <- analyze_gwas(
#'   genotype_file = "results/gwas/demo.geno.tsv",
#'   phenotype = "results/gwas/demo.pheno.tsv",
#'   out_prefix = "results/gwas/demo_analysis",
#'   trait = "trait",
#'   fixed_effects = "pop",
#'   n_pcs = 2
#' )
#' }
#' @export
analyze_gwas <- function(
    genotype_file,
    phenotype,
    out_prefix,
    trait = "trait",
    sample_id = "sample",
    fixed_effects = NULL,
    kinship = NULL,
    n_pcs = 0L,
    min_maf = 0.05,
    n_cores = 1L,
    p3d = TRUE,
    fdr_threshold = 0.05,
    backend = c("auto", "rrBLUP", "GAPIT", "native")) {
  backend <- match.arg(backend)
  if (backend == "auto") backend <- if (requireNamespace("rrBLUP", quietly = TRUE)) "rrBLUP" else "native"
  if (backend == "rrBLUP") .bioana_require_package("rrBLUP", "GWAS analysis")
  if (backend == "GAPIT") .bioana_require_package("GAPIT", "GWAS analysis (GAPIT)")
  if (backend == "native" && !isTRUE(p3d)) stop("p3d = FALSE (marker-specific REML) needs backend = 'rrBLUP'")
  genotype_data <- .bioana_read_genotypes(genotype_file)
  sample_ids <- colnames(genotype_data$genotype)
  phenotype <- .bioana_prepare_phenotypes(
    phenotype,
    sample_ids,
    sample_id,
    trait,
    required_columns = fixed_effects
  )
  observed <- is.finite(phenotype[[trait]])
  if (sum(observed) < 10L) {
    stop("At least ten samples with finite phenotypes are required for GWAS")
  }
  if (!is.numeric(min_maf) || length(min_maf) != 1L ||
      !is.finite(min_maf) || min_maf < 0 || min_maf >= 0.5) {
    stop("min_maf must be one number in [0, 0.5)")
  }
  if (!is.numeric(fdr_threshold) || length(fdr_threshold) != 1L ||
      !is.finite(fdr_threshold) || fdr_threshold <= 0 ||
      fdr_threshold >= 1) {
    stop("fdr_threshold must be one number between zero and one")
  }

  observed_ids <- sample_ids[observed]
  dosage <- genotype_data$genotype[, observed, drop = FALSE]
  rr_dosage <- dosage - 1
  variants <- genotype_data$variants
  geno <- data.frame(
    marker = variants$id,
    chrom = variants$seqname,
    pos = variants$pos,
    rr_dosage,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  names(geno)[seq.int(4L, ncol(geno))] <- observed_ids

  pheno_columns <- unique(c(sample_id, trait, fixed_effects))
  pheno <- phenotype[observed, pheno_columns, drop = FALSE]
  names(pheno)[names(pheno) == sample_id] <- "line"
  if (is.null(kinship) && backend != "GAPIT") {
    kinship <- if (backend == "rrBLUP") rrBLUP::A.mat(t(rr_dosage)) else .bioana_kinship(rr_dosage)
  } else {
    kinship <- as.matrix(kinship)
    if (is.null(rownames(kinship)) || is.null(colnames(kinship))) {
      stop("kinship must have sample IDs as row and column names")
    }
    if (any(!observed_ids %in% rownames(kinship)) ||
        any(!observed_ids %in% colnames(kinship))) {
      stop("kinship does not contain every phenotyped sample")
    }
    kinship <- kinship[observed_ids, observed_ids, drop = FALSE]
  }

  if (backend == "rrBLUP") {
    scores <- rrBLUP::GWAS(
      pheno = pheno,
      geno = geno,
      fixed = fixed_effects,
      K = kinship,
      n.PC = as.integer(n_pcs),
      min.MAF = min_maf,
      n.core = as.integer(n_cores),
      P3D = isTRUE(p3d),
      plot = FALSE
    )
    score_column <- setdiff(names(scores), c("marker", "chrom", "pos"))
    if (length(score_column) != 1L) {
      stop("rrBLUP returned an unexpected number of trait columns")
    }
    score <- as.numeric(scores[[score_column]])
    beta <- se <- rep(NA_real_, length(score))
    h2_null <- NA_real_
  } else if (backend == "GAPIT") {
    gp <- .bioana_gwas_gapit(rr_dosage + 1, variants, pheno, trait, as.integer(n_pcs))
    scores <- gp[, c("marker", "chrom", "pos")]
    score <- -log10(gp$p)
    beta <- gp$beta
    se <- rep(NA_real_, nrow(gp))
    h2_null <- NA_real_
  } else {
    covariates <- if (length(fixed_effects)) {
      stats::model.matrix(stats::reformulate(fixed_effects), data = pheno)[, -1L, drop = FALSE]
    } else {
      NULL
    }
    native <- .bioana_gwas_native(rr_dosage, as.numeric(pheno[[trait]]), covariates,
                                  kinship, as.integer(n_pcs), min_maf)
    scores <- data.frame(marker = variants$id[native$index], chrom = variants$seqname[native$index],
                         pos = variants$pos[native$index], stringsAsFactors = FALSE)
    score <- -log10(native$p_value)
    beta <- native$beta
    se <- native$se
    h2_null <- native$h2
  }
  p_value <- pmin(1, 10^(-score))
  q_value <- stats::p.adjust(p_value, method = "BH")
  allele_frequency <- rowMeans(dosage, na.rm = TRUE) / 2
  maf <- pmin(allele_frequency, 1 - allele_frequency)
  result <- data.frame(
    marker_id = as.character(scores$marker),
    seqname = as.character(scores$chrom),
    pos = as.numeric(scores$pos),
    score = score,
    beta = beta,
    se = se,
    p_value = p_value,
    q_value = q_value,
    maf = maf[match(as.character(scores$marker), variants$id)],
    significant = !is.na(q_value) & q_value <= fdr_threshold,
    stringsAsFactors = FALSE
  )
  result <- result[order(result$p_value, na.last = TRUE), , drop = FALSE]
  rownames(result) <- NULL

  result_path <- paste0(out_prefix, ".gwas.tsv")
  summary_path <- paste0(out_prefix, ".gwas_summary.json")
  dir.create(dirname(result_path), recursive = TRUE, showWarnings = FALSE)
  utils::write.table(
    result,
    result_path,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  summary <- list(
    backend = switch(backend, rrBLUP = "rrBLUP::GWAS", GAPIT = "GAPIT::GAPIT (MLM)", native = "native EMMAX (P3D)"),
    trait = trait,
    samples = sum(observed),
    markers_tested = nrow(result),
    fixed_effects = fixed_effects %||% character(),
    principal_components = as.integer(n_pcs),
    minimum_maf = min_maf,
    p3d = isTRUE(p3d),
    fdr_threshold = fdr_threshold,
    discoveries = sum(result$significant, na.rm = TRUE),
    null_model_h2 = h2_null
  )
  .bioana_write_json(summary, summary_path)
  invisible(list(
    results = result,
    kinship = kinship,
    summary = summary,
    paths = list(results = result_path, summary = summary_path)
  ))
}

#' Run a fixed-effect logistic GWAS for a binary trait
#'
#' Test binary case/control, presence/absence, or other Bernoulli phenotypes
#' marker by marker with logistic regression. Optional categorical covariates
#' and marker-derived principal components can be included to reduce measured
#' structure confounding. The function writes the same standard association
#' columns as [analyze_gwas()], so [plot_gwas_results()] can create Manhattan
#' and QQ plots from its output.
#'
#' This is a transparent fixed-effect baseline, not a replacement for a
#' scalable binary mixed-model GWAS. For strongly related cohorts, large human
#' biobanks, severe case-control imbalance, or separation, use an appropriate
#' externally validated mixed-model method such as SAIGE, GMMAT, or GENESIS and
#' retain the model specification in the result manifest.
#'
#' @param genotype_file VCF, VCF.GZ, or marker-by-sample genotype TSV path.
#' @param phenotype Data frame or phenotype TSV containing one row per sample.
#' @param out_prefix Prefix for `.gwas.tsv` and `.gwas_summary.json` outputs.
#' @param trait Binary phenotype column encoded as `0/1`, `FALSE/TRUE`, or a
#'   two-level factor/character vector.
#' @param sample_id Sample-ID column in `phenotype`.
#' @param fixed_effects Optional phenotype columns included as fixed effects.
#' @param n_pcs Number of marker-derived principal components to include.
#' @param min_maf Minimum minor allele frequency for a marker to be tested.
#' @param fdr_threshold BH-adjusted p-value cutoff used to label discoveries.
#'
#' @return Invisibly returns the association results, summary, and output paths.
#' @examples
#' \dontrun{
#' binary_gwas <- analyze_binary_gwas(
#'   "results/nam.vcf", "results/case_control.tsv", "results/binary_gwas",
#'   trait = "case", fixed_effects = "family", n_pcs = 2
#' )
#' plot_gwas_results(binary_gwas$paths$results, "results/figures/binary_gwas")
#' }
#' @export
analyze_binary_gwas <- function(
    genotype_file,
    phenotype,
    out_prefix,
    trait = "trait",
    sample_id = "sample",
    fixed_effects = NULL,
    n_pcs = 0L,
    min_maf = 0.05,
    fdr_threshold = 0.05) {
  if (!is.numeric(min_maf) || length(min_maf) != 1L ||
      !is.finite(min_maf) || min_maf < 0 || min_maf >= 0.5) {
    stop("min_maf must be one number in [0, 0.5)")
  }
  if (!is.numeric(fdr_threshold) || length(fdr_threshold) != 1L ||
      !is.finite(fdr_threshold) || fdr_threshold <= 0 || fdr_threshold >= 1) {
    stop("fdr_threshold must be one number between zero and one")
  }
  genotype_data <- .bioana_read_genotypes(genotype_file)
  sample_ids <- colnames(genotype_data$genotype)
  phenotype <- .bioana_read_analysis_table(phenotype, "Phenotype")
  required_columns <- unique(c(sample_id, trait, fixed_effects))
  missing_columns <- setdiff(required_columns, names(phenotype))
  if (length(missing_columns)) {
    stop("Phenotype data are missing columns: ", paste(missing_columns, collapse = ", "))
  }
  ids <- as.character(phenotype[[sample_id]])
  if (anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids)) {
    stop("Phenotype sample IDs must be present and unique")
  }
  phenotype[[sample_id]] <- ids
  phenotype <- phenotype[match(sample_ids, ids), , drop = FALSE]
  phenotype[[sample_id]] <- sample_ids
  raw_trait <- phenotype[[trait]]
  if (is.logical(raw_trait)) {
    outcome <- as.integer(raw_trait)
  } else if (is.numeric(raw_trait)) {
    outcome <- as.integer(raw_trait)
  } else {
    levels <- unique(stats::na.omit(as.character(raw_trait)))
    if (length(levels) != 2L) {
      stop("A non-numeric binary trait must have exactly two observed levels")
    }
    outcome <- match(as.character(raw_trait), levels) - 1L
  }
  if (any(!is.na(outcome) & !outcome %in% c(0L, 1L))) {
    stop("Binary trait values must be 0/1, FALSE/TRUE, or exactly two levels")
  }
  covariates <- phenotype[, fixed_effects %||% character(), drop = FALSE]
  observed <- !is.na(outcome)
  if (ncol(covariates)) {
    observed <- observed & stats::complete.cases(covariates)
  }
  if (sum(observed) < 10L || length(unique(outcome[observed])) != 2L) {
    stop("Binary GWAS requires at least ten observed samples and both outcome classes")
  }

  dosage <- genotype_data$genotype[, observed, drop = FALSE]
  samples <- sample_ids[observed]
  outcome <- outcome[observed]
  covariates <- covariates[observed, , drop = FALSE]
  if (ncol(covariates)) {
    covariates[] <- lapply(covariates, function(x) {
      if (is.character(x)) factor(x) else x
    })
  }
  maf <- rowMeans(dosage, na.rm = TRUE) / 2
  maf <- pmin(maf, 1 - maf)
  keep <- is.finite(maf) & maf >= min_maf
  if (!any(keep)) stop("No markers meet min_maf after phenotype alignment")

  n_pcs <- as.integer(n_pcs)
  pcs <- NULL
  if (n_pcs > 0L) {
    centered <- t(dosage[keep, , drop = FALSE])
    center_values <- colMeans(centered, na.rm = TRUE)
    centered[!is.finite(centered)] <- center_values[col(centered)[!is.finite(centered)]]
    pca <- stats::prcomp(centered, center = TRUE, scale. = TRUE)
    n_pcs <- min(n_pcs, ncol(pca$x))
    pcs <- as.data.frame(pca$x[, seq_len(n_pcs), drop = FALSE])
    names(pcs) <- paste0("PC", seq_len(n_pcs))
  }

  tested <- which(keep)
  scan <- lapply(tested, function(i) {
    marker <- dosage[i, ]
    marker[!is.finite(marker)] <- mean(marker, na.rm = TRUE)
    if (!is.finite(stats::var(marker)) || stats::var(marker) == 0) {
      return(c(beta = NA_real_, se = NA_real_, p_value = NA_real_))
    }
    frame <- data.frame(outcome = outcome, dosage = marker, covariates,
                        check.names = FALSE)
    if (!is.null(pcs)) frame <- cbind(frame, pcs)
    fit <- suppressWarnings(tryCatch(
      stats::glm(outcome ~ ., data = frame, family = stats::binomial()),
      error = function(e) NULL
    ))
    if (is.null(fit) || !isTRUE(fit$converged)) {
      return(c(beta = NA_real_, se = NA_real_, p_value = NA_real_))
    }
    coefficients <- summary(fit)$coefficients
    if (!"dosage" %in% rownames(coefficients)) {
      return(c(beta = NA_real_, se = NA_real_, p_value = NA_real_))
    }
    c(beta = coefficients["dosage", "Estimate"],
      se = coefficients["dosage", "Std. Error"],
      p_value = coefficients["dosage", "Pr(>|z|)"])
  })
  scan <- as.data.frame(do.call(rbind, scan))
  variants <- genotype_data$variants[tested, , drop = FALSE]
  result <- data.frame(
    marker_id = as.character(variants$id),
    seqname = as.character(variants$seqname),
    pos = as.numeric(variants$pos),
    beta = scan$beta,
    se = scan$se,
    odds_ratio = exp(scan$beta),
    p_value = scan$p_value,
    q_value = stats::p.adjust(scan$p_value, method = "BH"),
    maf = maf[tested],
    samples = length(samples),
    stringsAsFactors = FALSE
  )
  result$significant <- !is.na(result$q_value) & result$q_value <= fdr_threshold
  result <- result[order(result$p_value, na.last = TRUE), , drop = FALSE]
  result_path <- paste0(out_prefix, ".gwas.tsv")
  summary_path <- paste0(out_prefix, ".gwas_summary.json")
  dir.create(dirname(result_path), recursive = TRUE, showWarnings = FALSE)
  utils::write.table(result, result_path, sep = "\t", quote = FALSE, row.names = FALSE)
  summary <- list(
    backend = "stats::glm(binomial)", trait = trait, samples = length(samples),
    cases = sum(outcome), controls = sum(outcome == 0L), markers_tested = nrow(result),
    fixed_effects = fixed_effects %||% character(), principal_components = n_pcs,
    minimum_maf = min_maf, fdr_threshold = fdr_threshold,
    discoveries = sum(result$significant, na.rm = TRUE),
    limitation = "Fixed-effect logistic baseline; no kinship random effect."
  )
  .bioana_write_json(summary, summary_path)
  invisible(list(results = result, summary = summary,
                 paths = list(results = result_path, summary = summary_path)))
}

#' Benchmark GWAS discoveries against simulated causal variants
#'
#' Compare an association table from [analyze_gwas()] with the causal-variant
#' truth table written by `simitall::simulate_gwas_cohort()`. Recovery can be evaluated
#' by exact marker identity or by a chromosome-specific distance window.
#'
#' @param gwas_results GWAS result data frame or `.gwas.tsv` path.
#' @param truth Causal-variant data frame or `.causal.tsv` path.
#' @param out_prefix Prefix for annotated results and benchmark metrics.
#' @param fdr_threshold Discovery threshold applied to `q_value`.
#' @param window_bp Non-negative distance around each causal variant. A value of
#'   zero requires an exact marker or coordinate match.
#' @param ld_r2 Optional r2 threshold. With `genotypes`, a marker on the same
#'   chromosome whose r2 with a causal marker is at least `ld_r2` also counts as
#'   hitting it (for crosses, where linked markers carry the same signal).
#'   `NULL` (default) uses `window_bp` only.
#' @param genotypes Genotype file (VCF or TSV) or a marker-by-sample dosage
#'   matrix with marker IDs as row names; needed with `ld_r2`.
#'
#' @return Invisibly returns annotated associations, metrics, and output paths.
#' @examples
#' \dontrun{
#' benchmark_gwas(
#'   "results/gwas/demo_analysis.gwas.tsv",
#'   "results/gwas/demo.causal.tsv",
#'   "results/gwas/demo_benchmark"
#' )
#' }
#' @export
benchmark_gwas <- function(
    gwas_results,
    truth,
    out_prefix,
    fdr_threshold = 0.05,
    window_bp = 0L,
    ld_r2 = NULL,
    genotypes = NULL) {
  results <- .bioana_read_analysis_table(gwas_results, "GWAS result")
  truth <- .bioana_read_analysis_table(truth, "GWAS truth")
  required_results <- c("marker_id", "seqname", "pos", "q_value")
  required_truth <- c("marker_id", "seqname", "pos")
  if (any(!required_results %in% names(results))) {
    stop("GWAS results need marker_id, seqname, pos, and q_value columns")
  }
  if (any(!required_truth %in% names(truth))) {
    stop("GWAS truth needs marker_id, seqname, and pos columns")
  }
  if (!nrow(truth)) stop("GWAS truth contains no causal variants")
  window_bp <- as.integer(window_bp)
  if (length(window_bp) != 1L || is.na(window_bp) || window_bp < 0L) {
    stop("window_bp must be one non-negative integer")
  }

  exact <- results$marker_id %in% truth$marker_id
  # linked[i, j]: result marker i hits causal variant j (same chromosome, within window_bp or r2 >= ld_r2)
  linked <- outer(as.character(results$seqname), as.character(truth$seqname), "==") &
    abs(outer(as.numeric(results$pos), as.numeric(truth$pos), "-")) <= window_bp
  if (!is.null(ld_r2)) {
    if (is.null(genotypes)) stop("ld_r2 needs genotypes")
    dosage <- if (is.matrix(genotypes)) genotypes else .bioana_read_genotypes(genotypes)$genotype
    r2 <- .bioana_ld_r2(dosage, as.character(results$marker_id), as.character(truth$marker_id))
    linked <- linked | (r2 >= ld_r2 & outer(as.character(results$seqname), as.character(truth$seqname), "=="))
  }
  near <- rowSums(linked) > 0
  results$is_causal <- exact
  results$within_causal_window <- near
  results$discovery <- !is.na(results$q_value) &
    results$q_value <= fdr_threshold

  recovered_truth <- colSums(linked & results$discovery) > 0
  true_positive <- sum(results$discovery & near)
  false_positive <- sum(results$discovery & !near)
  false_negative <- sum(!recovered_truth)
  precision <- if (true_positive + false_positive == 0L) {
    NA_real_
  } else {
    true_positive / (true_positive + false_positive)
  }
  recall <- mean(recovered_truth)
  f1 <- if (!is.finite(precision) || precision + recall == 0) {
    NA_real_
  } else {
    2 * precision * recall / (precision + recall)
  }
  metrics <- data.frame(
    fdr_threshold = fdr_threshold,
    window_bp = window_bp,
    ld_r2 = ld_r2 %||% NA_real_,
    causal_variants = nrow(truth),
    discoveries = sum(results$discovery),
    true_positive_markers = true_positive,
    false_positive_markers = false_positive,
    unrecovered_causal_variants = false_negative,
    precision = precision,
    recall = recall,
    f1 = f1,
    stringsAsFactors = FALSE
  )

  result_path <- paste0(out_prefix, ".benchmark.tsv")
  metric_path <- paste0(out_prefix, ".benchmark_metrics.tsv")
  summary_path <- paste0(out_prefix, ".benchmark_summary.json")
  dir.create(dirname(result_path), recursive = TRUE, showWarnings = FALSE)
  utils::write.table(
    results, result_path, sep = "\t", quote = FALSE, row.names = FALSE
  )
  utils::write.table(
    metrics, metric_path, sep = "\t", quote = FALSE, row.names = FALSE
  )
  .bioana_write_json(as.list(metrics[1L, ]), summary_path)
  invisible(list(
    results = results,
    truth = transform(truth, recovered = recovered_truth),
    metrics = metrics,
    paths = list(
      results = result_path,
      metrics = metric_path,
      summary = summary_path
    )
  ))
}

#' Plot GWAS Manhattan and quantile-quantile panels
#'
#' Create publication-ready Manhattan and QQ plots from [analyze_gwas()] or
#' [benchmark_gwas()] output. Benchmark annotations are used to highlight
#' causal regions when present.
#'
#' @param gwas_results GWAS or benchmark result data frame or TSV path.
#' @param out_prefix Prefix for PDF and PNG figures.
#' @param fdr_threshold Adjusted-p-value threshold shown as a horizontal line.
#' @param width,height Figure dimensions in inches.
#' @param dpi PNG resolution.
#' @param individual_panels_dir Optional directory for standalone Manhattan and
#'   QQ PNG/PDF panels. Defaults to a `figures` folder beside `out_prefix`.
#'
#' @return Invisibly returns the combined plot, source data, and output paths.
#' @examples
#' \dontrun{
#' plot_gwas_results(
#'   "results/gwas/demo_analysis.gwas.tsv",
#'   "results/figures/demo_gwas"
#' )
#' }
#' @export
plot_gwas_results <- function(
    gwas_results,
    out_prefix,
    fdr_threshold = 0.05,
    width = 13,
    height = 6,
    dpi = 300,
    individual_panels_dir = file.path(dirname(out_prefix), "figures")) {
  .bioana_require_package("ggplot2", "GWAS plotting")
  .bioana_require_package("patchwork", "GWAS plotting")
  results <- .bioana_read_analysis_table(gwas_results, "GWAS result")
  required <- c("marker_id", "seqname", "pos", "p_value", "q_value")
  if (any(!required %in% names(results))) {
    stop("GWAS results need marker_id, seqname, pos, p_value, and q_value")
  }
  results <- results[is.finite(results$pos) & is.finite(results$p_value), ]
  if (!nrow(results)) stop("GWAS results contain no finite association tests")
  chromosome_levels <- unique(as.character(results$seqname))
  results$seqname <- factor(results$seqname, levels = chromosome_levels)
  chromosome_max <- tapply(results$pos, results$seqname, max, na.rm = TRUE)
  offsets <- c(0, head(cumsum(as.numeric(chromosome_max)), -1L))
  names(offsets) <- names(chromosome_max)
  results$cumulative_pos <- results$pos + offsets[as.character(results$seqname)]
  axis_positions <- tapply(
    results$cumulative_pos,
    results$seqname,
    function(value) mean(range(value, na.rm = TRUE))
  )
  results$neglog10_p <- -log10(pmax(results$p_value, .Machine$double.xmin))
  results$highlight <- if ("within_causal_window" %in% names(results)) {
    ifelse(results$within_causal_window, "Causal region", "Background")
  } else {
    "Background"
  }
  threshold_candidates <- results$p_value[
    !is.na(results$q_value) & results$q_value <= fdr_threshold
  ]
  threshold <- if (length(threshold_candidates)) {
    -log10(max(threshold_candidates))
  } else {
    -log10(0.05 / nrow(results))
  }

  manhattan <- ggplot2::ggplot(
    results,
    ggplot2::aes(x = cumulative_pos, y = neglog10_p)
  ) +
    ggplot2::geom_point(
      ggplot2::aes(color = highlight),
      alpha = 0.75,
      size = 1.2
    ) +
    ggplot2::geom_hline(
      yintercept = threshold,
      linetype = 2,
      color = "#B23A48"
    ) +
    ggplot2::scale_x_continuous(
      breaks = as.numeric(axis_positions),
      labels = names(axis_positions)
    ) +
    ggplot2::scale_color_manual(
      values = c("Background" = "#315C6B", "Causal region" = "#D97706")
    ) +
    ggplot2::labs(
      title = "Genome-wide association landscape",
      x = "Chromosome",
      y = expression(-log[10](p)),
      color = NULL
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(legend.position = "bottom")

  observed <- sort(results$p_value)
  expected <- stats::ppoints(length(observed))
  qq_data <- data.frame(
    expected = -log10(expected),
    observed = -log10(pmax(observed, .Machine$double.xmin))
  )
  qq <- ggplot2::ggplot(
    qq_data,
    ggplot2::aes(x = expected, y = observed)
  ) +
    ggplot2::geom_abline(slope = 1, intercept = 0, color = "#B23A48") +
    ggplot2::geom_point(color = "#315C6B", alpha = 0.75, size = 1.3) +
    ggplot2::coord_equal() +
    ggplot2::labs(
      title = "Quantile-quantile plot",
      x = expression(Expected~-log[10](p)),
      y = expression(Observed~-log[10](p))
    ) +
    ggplot2::theme_minimal(base_size = 11)

  combined <- manhattan + qq + patchwork::plot_layout(widths = c(2, 1))
  pdf_path <- paste0(out_prefix, ".pdf")
  png_path <- paste0(out_prefix, ".png")
  source_path <- paste0(out_prefix, ".source_data.tsv")
  dir.create(dirname(pdf_path), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(pdf_path, combined, width = width, height = height)
  ggplot2::ggsave(
    png_path,
    combined,
    width = width,
    height = height,
    dpi = dpi
  )
  utils::write.table(
    results, source_path, sep = "\t", quote = FALSE, row.names = FALSE
  )
  individual_panels <- if (is.null(individual_panels_dir)) {
    NULL
  } else {
    export_bioana_plot_panels(
      list(manhattan = manhattan, qq = qq), individual_panels_dir,
      c("gwas_manhattan", "gwas_qq"), width = 7, height = 5, dpi = dpi
    )
  }
  invisible(list(
    plot = combined,
    manhattan = manhattan,
    qq = qq,
    source_data = list(manhattan = results, qq = qq_data),
    paths = list(
      pdf = pdf_path, png = png_path, source_data = source_path,
      individual_panels = individual_panels
    )
  ))
}

.bioana_prepare_phenotypes <- function(
    phenotype,
    genotype_samples,
    sample_id,
    trait,
    required_columns = character()) {
  phenotype <- .bioana_read_analysis_table(phenotype, "Phenotype")
  columns <- unique(c(sample_id, trait, required_columns))
  missing_columns <- setdiff(columns, names(phenotype))
  if (length(missing_columns)) {
    stop(
      "Phenotype data are missing columns: ",
      paste(missing_columns, collapse = ", ")
    )
  }
  ids <- as.character(phenotype[[sample_id]])
  if (anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids)) {
    stop("Phenotype sample IDs must be present and unique")
  }
  phenotype[[sample_id]] <- ids
  phenotype <- phenotype[match(genotype_samples, ids), , drop = FALSE]
  phenotype[[sample_id]] <- genotype_samples
  suppressWarnings(phenotype[[trait]] <- as.numeric(phenotype[[trait]]))
  phenotype
}

.bioana_require_package <- function(package, feature) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop(
      "Package '", package, "' is required for ", feature, ". ",
      "Install the population dependency profile or run install.packages('",
      package, "').",
      call. = FALSE
    )
  }
}

.bioana_write_json <- function(value, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  jsonlite::write_json(
    value,
    path,
    auto_unbox = TRUE,
    pretty = TRUE,
    na = "null"
  )
  path
}

.bioana_read_analysis_table <- function(data, label) {
  if (is.character(data) && length(data) == 1L) {
    if (!file.exists(data)) stop(label, " file does not exist: ", data)
    data <- utils::read.delim(
      data,
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
  }
  if (!is.data.frame(data)) {
    stop(label, " must be a data frame or a tab-separated file path")
  }
  as.data.frame(data, check.names = FALSE, stringsAsFactors = FALSE)
}

# Additive relationship matrix from -1/0/1 dosages (markers x samples), as rrBLUP::A.mat
# without its marker filters: W = allele count - 2p, A = W W' / (2 sum p(1 - p)).
.bioana_kinship <- function(rr_dosage) {
  counts <- t(rr_dosage) + 1
  freq <- colMeans(counts, na.rm = TRUE) / 2
  w <- sweep(counts, 2L, 2 * freq)
  w[is.na(w)] <- 0
  keep <- freq > 0 & freq < 1
  tcrossprod(w[, keep, drop = FALSE]) / (2 * sum(freq[keep] * (1 - freq[keep])))
}

# EMMAX with P3D: y = X b + g + e, var(y) = s2 (K + delta I). delta by REML on the
# eigenbasis of K, then GLS of y on X and each marker. Markers are -1/0/1 dosages
# (markers x samples, missing imputed with the marker mean). Returns the tested
# marker indices with beta (per alternative allele), se and p.
.bioana_gwas_native <- function(rr_dosage, y, covariates, kinship, n_pcs = 0L, min_maf = 0.05, chunk = 5000L) {
  n <- length(y)
  eig <- eigen(kinship, symmetric = TRUE)
  X <- matrix(1, n, 1L)
  if (!is.null(covariates) && ncol(covariates)) X <- cbind(X, covariates)
  if (n_pcs > 0L) X <- cbind(X, eig$vectors[, seq_len(min(n_pcs, n - ncol(X) - 2L)), drop = FALSE])
  qx <- qr(X)
  X <- X[, qx$pivot[seq_len(qx$rank)], drop = FALSE]
  p <- ncol(X)
  if (n - p - 1L < 2L) stop("Too few lines for the covariates in the GWAS model")
  d <- pmax(eig$values, 0)
  U <- eig$vectors
  ys <- drop(crossprod(U, y))
  Xs <- crossprod(U, X)
  reml <- function(log_delta) {
    w <- 1 / (d + exp(log_delta))
    xtwx <- crossprod(Xs, Xs * w)
    b <- solve(xtwx, crossprod(Xs, w * ys))
    r <- ys - drop(Xs %*% b)
    s2 <- sum(w * r^2) / (n - p)
    -0.5 * ((n - p) * log(s2) + sum(log(d + exp(log_delta))) + as.numeric(determinant(xtwx)$modulus))
  }
  delta <- exp(stats::optimize(reml, c(-10, 10), maximum = TRUE)$maximum)
  w <- 1 / (d + delta)
  xtwx_inv <- solve(crossprod(Xs, Xs * w))
  resid_w <- function(m) m - Xs %*% (xtwx_inv %*% crossprod(Xs, w * m))
  yr <- drop(resid_w(ys))
  ywy <- sum(w * yr^2)
  counts <- rr_dosage + 1
  freq <- rowMeans(counts, na.rm = TRUE) / 2
  tested <- which(is.finite(freq) & pmin(freq, 1 - freq) >= min_maf)
  if (!length(tested)) stop("No markers meet min_maf among the phenotyped samples")
  df <- n - p - 1L
  beta <- se <- p_value <- rep(NA_real_, length(tested))
  for (start in seq(1L, length(tested), by = chunk)) {
    k <- seq.int(start, min(start + chunk - 1L, length(tested)))
    z <- rr_dosage[tested[k], , drop = FALSE]
    miss <- is.na(z)
    if (any(miss)) z[miss] <- (2 * freq[tested[k]] - 1)[row(z)[miss]]
    zr <- resid_w(crossprod(U, t(z)))
    zwz <- colSums(w * zr^2)
    zwy <- colSums(w * zr * yr)
    ok <- zwz > 1e-8 * max(zwz)
    b <- zwy / zwz
    s2 <- (ywy - b * zwy) / df
    s <- sqrt(s2 / zwz)
    beta[k] <- ifelse(ok, b, NA_real_)
    se[k] <- ifelse(ok, s, NA_real_)
    p_value[k] <- ifelse(ok, 2 * stats::pt(-abs(b / s), df), NA_real_)
  }
  list(index = tested, beta = beta, se = se, p_value = p_value,
       h2 = mean(d) / (mean(d) + delta))
}

# r2 between result markers (rows) and truth markers (columns); NA-safe, 0 when a marker is missing or monomorphic.
.bioana_ld_r2 <- function(dosage, result_ids, truth_ids) {
  out <- matrix(0, length(result_ids), length(truth_ids))
  tr <- match(truth_ids, rownames(dosage))
  rs <- match(result_ids, rownames(dosage))
  ok_t <- which(!is.na(tr)); ok_r <- which(!is.na(rs))
  if (!length(ok_t) || !length(ok_r)) return(out)
  r <- suppressWarnings(stats::cor(t(dosage[rs[ok_r], , drop = FALSE]), t(dosage[tr[ok_t], , drop = FALSE]),
                                   use = "pairwise.complete.obs"))
  r[is.na(r)] <- 0
  out[ok_r, ok_t] <- r^2
  out
}

# GAPIT MLM on 0/1/2 dosages (markers x samples). GAPIT's output column names have changed
# between releases, so the p value and effect columns are found by name pattern.
.bioana_gwas_gapit <- function(dosage, variants, pheno, trait, n_pcs) {
  ids <- colnames(dosage)
  Y <- data.frame(Taxa = pheno$line, trait = pheno[[trait]], stringsAsFactors = FALSE)
  GD <- data.frame(Taxa = ids, t(dosage), check.names = FALSE, stringsAsFactors = FALSE)
  GM <- data.frame(SNP = variants$id, Chromosome = sub("^chr", "", variants$seqname), Position = variants$pos,
                   stringsAsFactors = FALSE)
  out <- tempfile("bioana_gapit_"); dir.create(out); old <- setwd(out); on.exit(setwd(old), add = TRUE)
  res <- suppressWarnings(utils::capture.output(g <- getExportedValue("GAPIT", "GAPIT")(Y = Y, GD = GD, GM = GM, model = "MLM",
                                                                PCA.total = n_pcs, file.output = FALSE)))
  tab <- g$GWAS
  if (is.null(tab) || !nrow(tab)) stop("GAPIT returned no GWAS table")
  pcol <- grep("^p[._]?value$", names(tab), ignore.case = TRUE, value = TRUE)[1]
  ecol <- grep("effect", names(tab), ignore.case = TRUE, value = TRUE)[1]
  if (is.na(pcol)) stop("GAPIT GWAS table has no P.value column: ", paste(names(tab), collapse = ", "))
  m <- match(as.character(tab$SNP), variants$id)
  data.frame(marker = variants$id[m], chrom = variants$seqname[m], pos = variants$pos[m],
             p = as.numeric(tab[[pcol]]), beta = if (is.na(ecol)) NA_real_ else as.numeric(tab[[ecol]]),
             stringsAsFactors = FALSE)
}
