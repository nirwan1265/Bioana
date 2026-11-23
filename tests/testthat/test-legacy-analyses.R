# Moved from Simitall: analyses of simulated RNA-seq, single-cell, ChIP-seq data and GO annotation.
test_that("eQTL benchmarking writes association results and metrics", {
  skip_if_not_installed("simitall")
  out_dir <- tempfile("simitall-eqtl-")
  dir.create(out_dir)
  genome <- system.file(
    "extdata", "examples", "demo_genome.fa",
    package = "simitall"
  )
  gwas_prefix <- file.path(out_dir, "cohort")
  rna_prefix <- file.path(out_dir, "rna")

  simulate_gwas_cohort(
    genome_fa = genome,
    out_prefix = gwas_prefix,
    n_samples = 16,
    snp_rate = 0.01,
    indel_rate = 0,
    n_causal = 2,
    seed = 31
  )
  paths <- simulate_rnaseq_from_gwas(
    genotype_file = paste0(gwas_prefix, ".geno.tsv"),
    out_prefix = rna_prefix,
    n_genes = 16,
    n_cis_eqtl = 4,
    n_trans_eqtl = 2,
    cis_window_bp = 400,
    n_latent = 1,
    library_size_mean = 1e5,
    seed = 32
  )
  benchmark <- benchmark_eqtl(
    genotype_file = paste0(gwas_prefix, ".geno.tsv"),
    expression_file = paths$expression,
    out_prefix = file.path(out_dir, "benchmark"),
    sample_metadata = paths$sample_metadata,
    truth_file = paths$eqtl_truth,
    cis_window_bp = 400
  )

  expect_true(file.exists(benchmark$results_path))
  expect_true(file.exists(benchmark$metrics_path))
  expect_gt(nrow(benchmark$results), 0)
  expect_true(all(c("p_value", "q_value", "is_causal") %in%
                    names(benchmark$results)))
  expect_equal(benchmark$metrics$truth_pairs_in_tests, 6)

  summary <- summarize_rnaseq_results(rna_prefix, file.path(out_dir, "benchmark"))
  expect_equal(summary$samples, 16)
  expect_equal(summary$genes, 16)
  expect_equal(summary$cis_eqtls, 4)

  skip_if_not_installed("ggplot2")
  skip_if_not_installed("patchwork")
  figure <- plot_rnaseq_eqtl_results(
    rnaseq_prefix = rna_prefix,
    eqtl_prefix = file.path(out_dir, "benchmark"),
    out_prefix = file.path(out_dir, "rnaseq_eqtl_figure"),
    top_variable_genes = 12,
    width = 9,
    height = 6,
    dpi = 72
  )
  expect_true(file.exists(figure$paths$pdf))
  expect_true(file.exists(figure$paths$png))
  expect_true(file.exists(figure$paths$summary))
  expect_length(list.files(figure$paths$source_data, pattern = "\\.csv$"), 6)
})

test_that("cell-type eQTL benchmarking and figures use donor pseudobulk", {
  skip_if_not_installed("simitall")
  out_dir <- tempfile("simitall-scrna-benchmark-")
  dir.create(out_dir)
  genome <- system.file(
    "extdata", "examples", "demo_genome.fa",
    package = "simitall"
  )
  gwas_prefix <- file.path(out_dir, "cohort")
  scrna_prefix <- file.path(out_dir, "scrna")
  eqtl_prefix <- file.path(out_dir, "celltype_eqtl")

  simulate_gwas_cohort(
    genome_fa = genome,
    out_prefix = gwas_prefix,
    n_samples = 12,
    snp_rate = 0.02,
    indel_rate = 0,
    n_causal = 2,
    seed = 81
  )
  paths <- simulate_scrnaseq_from_gwas(
    genotype_file = paste0(gwas_prefix, ".geno.tsv"),
    out_prefix = scrna_prefix,
    sample_metadata = paste0(gwas_prefix, ".pheno.tsv"),
    n_genes = 24,
    cells_per_donor = 15,
    cell_type_proportions = c(T_cell = 0.4, B_cell = 0.3, Monocyte = 0.3),
    n_cis_eqtl = 5,
    n_trans_eqtl = 2,
    cis_window_bp = 700,
    markers_per_cell_type = 3,
    n_latent = 1,
    library_size_mean = 1000,
    backend = "native",
    write_sce = FALSE,
    min_cells_pseudobulk = 3,
    seed = 82
  )
  benchmark <- benchmark_celltype_eqtls(
    genotype_file = paste0(gwas_prefix, ".geno.tsv"),
    pseudobulk_expression = paths$pseudobulk_expression,
    pseudobulk_metadata = paths$pseudobulk_metadata,
    truth_file = paths$eqtl_truth,
    out_prefix = eqtl_prefix,
    cis_window_bp = 700,
    min_donors = 8
  )

  expect_true(file.exists(benchmark$results_path))
  expect_true(file.exists(benchmark$metrics_path))
  expect_setequal(
    unique(benchmark$metrics$cell_type),
    c("T_cell", "B_cell", "Monocyte", "Overall")
  )
  expect_true(all(c("cell_type", "is_causal", "q_value") %in%
                    names(benchmark$results)))

  skip_if_not_installed("ggplot2")
  skip_if_not_installed("patchwork")
  figure <- plot_scrnaseq_results(
    scrna_prefix = scrna_prefix,
    eqtl_prefix = eqtl_prefix,
    out_prefix = file.path(out_dir, "scrna_figure"),
    top_variable_genes = 20,
    markers_per_cell_type = 2,
    max_plot_cells = 200,
    width = 10,
    height = 7,
    dpi = 72
  )
  expect_true(file.exists(figure$paths$pdf))
  expect_true(file.exists(figure$paths$png))
  expect_true(file.exists(figure$paths$summary))
  expect_length(
    list.files(figure$paths$source_data, pattern = "\\.csv$"),
    6
  )
})

test_that("ChIP-seq results figure from simulated reads", {
  skip_if_not_installed("simitall")
  out_dir <- tempfile("simitall-chipseq-")
  dir.create(out_dir)
  genome <- system.file(
    "extdata", "examples", "demo_genome.fa", package = "simitall"
  )
  seqname <- sub("^>([^ ]+).*", "\\1", readLines(genome, n = 1L))
  gff <- file.path(out_dir, "targets.gff3")
  writeLines(c(
    "##gff-version 3",
    paste(seqname, "test", "promoter", 500, 650, ".", "+", ".",
          "ID=promoter1", sep = "\t"),
    paste(seqname, "test", "enhancer", 1500, 1700, ".", "+", ".",
          "ID=enhancer1", sep = "\t"),
    paste(seqname, "test", "promoter", 3000, 3150, ".", "-", ".",
          "ID=promoter2", sep = "\t"),
    paste(seqname, "test", "enhancer", 5000, 5200, ".", "+", ".",
          "ID=enhancer2", sep = "\t")
  ), gff)
  prefix <- file.path(out_dir, "tf")

  paths <- simulate_chipseq(
    genome_fa = genome,
    annotation_gff3 = gff,
    out_prefix = prefix,
    assay_type = "TF",
    n_peaks = 6,
    peak_width = 120,
    peak_width_sd = 10,
    conditions = c("control", "treated"),
    biological_replicates = 2,
    technical_replicates = 1,
    differential_binding_fraction = 0.5,
    n_reads = 300,
    input_reads = 300,
    read_length = 40,
    fragment_mean = 100,
    fragment_sd = 15,
    fragment_min = 60,
    fragment_max = 150,
    backend = "native",
    write_fastq = TRUE,
    seed = 111
  )

  expect_identical(attr(paths, "backend"), "native")
  expect_true(all(file.exists(unlist(paths))))
  peaks <- read.delim(paths$peaks, check.names = FALSE)
  metadata <- read.delim(paths$sample_metadata, check.names = FALSE)
  qc <- read.delim(paths$qc, check.names = FALSE)
  manifest <- read.delim(paths$read_manifest, check.names = FALSE)
  expect_equal(nrow(peaks), 6)
  expect_equal(nrow(metadata), 4)
  expect_equal(nrow(qc), 8)
  expect_equal(nrow(manifest), 8)
  expect_true(all(file.exists(manifest$fastq)))
  expect_true(all(qc$frip >= 0 & qc$frip <= 1))
  expect_gt(mean(qc$frip[qc$assay == "chip"]),
            mean(qc$frip[qc$assay == "input"]))

  skip_if_not_installed("ggplot2")
  skip_if_not_installed("patchwork")
  figure <- plot_chipseq_results(
    chip_prefix = prefix,
    out_prefix = file.path(out_dir, "figure"),
    profile_window = 400,
    profile_bin = 25,
    max_profile_reads = 2400,
    width = 10,
    height = 7,
    dpi = 72
  )
  expect_true(file.exists(figure$paths$pdf))
  expect_true(file.exists(figure$paths$png))
  expect_true(file.exists(figure$paths$summary))
  expect_length(list.files(figure$paths$source_data, pattern = "\\.csv$"), 6)
})

test_that("GAF annotations join to compact gene coordinates", {
  gaf <- system.file("extdata", "go", "arabidopsis_tair10_go_demo.gaf.gz", package = "bioana")
  expect_true(file.exists(gaf))
  annotations <- read_go_gaf(gaf)
  expect_gt(nrow(annotations), 0L)
  expect_true(all(grepl("^GO:", annotations$go_id)))

  coordinates <- read_gene_coordinates("Arabidopsis_thaliana")
  joined <- join_gene_go(coordinates, annotations)
  expect_gt(nrow(joined), 0L)
  expect_true(all(joined$match_method == "stable_id"))
  summary <- summarize_go_terms(joined, by = "go_id_aspect")
  expect_gt(nrow(summary), 0L)
  expect_true(all(summary$gene_count <= summary$annotation_count))
})

test_that("compact gene-coordinate catalogues retain complete gene intervals", {
  manifest_path <- system.file("extdata", "gene_coordinates", "manifest.tsv", package = "bioana")
  manifest <- read.delim(manifest_path, check.names = FALSE, stringsAsFactors = FALSE)
  expect_equal(nrow(manifest), 6L)

  for (i in seq_len(nrow(manifest))) {
    path <- system.file("extdata", "gene_coordinates", manifest$coordinate_file[[i]], package = "bioana")
    expect_true(file.exists(path))
    genes <- read.delim(gzfile(path), check.names = FALSE, stringsAsFactors = FALSE)
    expect_equal(nrow(genes), manifest$gene_records[[i]])
    expect_true(all(nzchar(genes$seqname)))
    expect_true(all(nzchar(genes$gene_id)))
    expect_true(all(genes$start_bp <= genes$end_bp))
  }
})
