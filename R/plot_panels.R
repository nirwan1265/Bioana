# Saving multi-panel figures (moved from Simitall).

#' Export individual panels from a bioana result figure
#'
#' Save each component of a figure returned by [plot_gwas_results()],
#' [plot_rnaseq_eqtl_results()], [plot_chipseq_results()], or
#' [plot_scrnaseq_results()] as its own PNG and PDF. This lets a workflow keep
#' a compact multi-panel overview while also providing individual plots for a
#' report, Shiny download, or agent-generated result folder.
#'
#' @param figure_result Result object returned by a bioana plotting function.
#' @param out_dir Directory for individual figure files, normally
#'   `results/<run>/figures`.
#' @param panel_names Optional file stems. Supply one name per panel.
#' @param width,height Dimensions in inches for each panel.
#' @param dpi PNG resolution.
#'
#' @return A named list of PNG and PDF paths, invisibly.
#' @examples
#' \dontrun{
#' gwas_plot <- plot_gwas_results("results/gwas.benchmark.tsv", "results/gwas")
#' export_bioana_plot_panels(
#'   gwas_plot, "results/figures", c("gwas_manhattan", "gwas_qq")
#' )
#' }
#' @export
export_bioana_plot_panels <- function(
    figure_result,
    out_dir,
    panel_names = NULL,
    width = 7,
    height = 5,
    dpi = 300) {
  .bioana_require_figure_packages("individual figure export")
  if (!is.list(figure_result)) {
    stop("figure_result must be the list returned by a bioana plotting function")
  }
  panels <- figure_result$panels
  if (is.null(panels)) {
    panels <- Filter(Negate(is.null), list(
      manhattan = figure_result$manhattan,
      qq = figure_result$qq
    ))
  }
  if (!length(panels)) stop("No individual panels were found in figure_result")
  if (is.null(panel_names)) {
    panel_names <- names(panels)
    if (is.null(panel_names) || any(!nzchar(panel_names))) {
      panel_names <- paste0("panel_", seq_along(panels))
    }
  }
  if (length(panel_names) != length(panels) || any(!nzchar(panel_names))) {
    stop("panel_names must contain one non-empty name per panel")
  }
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  paths <- lapply(seq_along(panels), function(i) {
    prefix <- file.path(out_dir, panel_names[i])
    png <- paste0(prefix, ".png")
    pdf <- paste0(prefix, ".pdf")
    ggplot2::ggsave(png, panels[[i]], width = width, height = height, dpi = dpi)
    ggplot2::ggsave(pdf, panels[[i]], width = width, height = height)
    c(png = png, pdf = pdf)
  })
  names(paths) <- panel_names
  invisible(paths)
}

.bioana_require_figure_packages <- function(feature) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for ", feature, ".", call. = FALSE)
  }
}
