# Real data: the airway RNA-seq data set (Himes et al. 2014), the DESeq2 vignette's example.
# Bioana's run must equal the vignette code run directly, and find the known glucocorticoid-response genes.

test_that("airway: Bioana = the DESeq2 vignette code, and known dexamethasone genes go up", {
  skip_if_not_installed("DESeq2")
  d <- bioana_example("airway")
  man <- bioana_manifest(d)
  expect_equal(man$data[[1]]$data_type, "feature_counts")
  expect_equal(man$data[[1]]$source$format, "matrix")
  utils::capture.output(r <- bioana_analyze(d, secondary = FALSE))
  a <- r$summary$analyses$counts
  expect_equal(a$engine$package, "DESeq2")
  expect_equal(a$contrasts, "treatmenttrt")                 # untrt is the reference (spec level order)
  res <- r$results$counts
  cnt <- as.matrix(utils::read.delim(file.path(d, man$data[[1]]$files$counts), row.names = 1, check.names = FALSE))
  s <- utils::read.delim(file.path(d, "sample_metadata_template.tsv"), stringsAsFactors = FALSE)
  coldata <- data.frame(cell = factor(s$cell), dex = factor(s$dex), row.names = s$sample_id)
  dds <- DESeq2::DESeqDataSetFromMatrix(cnt[res$feature_id, rownames(coldata)], coldata, design = ~ cell + dex)
  dds$dex <- stats::relevel(dds$dex, ref = "untrt")
  dds <- suppressMessages(DESeq2::DESeq(dds, quiet = TRUE))
  v <- DESeq2::results(dds, name = "dex_trt_vs_untrt", alpha = 0.05)
  m <- res[match(rownames(v), res$feature_id), ]
  expect_equal(m$estimate, v$log2FoldChange, tolerance = 1e-6)
  expect_equal(m$p, v$pvalue, tolerance = 1e-4)
  expect_identical(is.na(m$padj), is.na(v$padj))
  expect_equal(sum(m$padj < 0.05, na.rm = TRUE), sum(v$padj < 0.05, na.rm = TRUE), tolerance = 5)
  known <- c(CRISPLD2 = "ENSG00000103196", DUSP1 = "ENSG00000120129", FKBP5 = "ENSG00000096060", TSC22D3 = "ENSG00000157514",
             KLF15 = "ENSG00000163884", PER1 = "ENSG00000179094", ZBTB16 = "ENSG00000109906")
  k <- res[match(known, res$feature_id), ]
  expect_true(all(k$estimate > 1))
  expect_true(all(k$padj < 1e-10))
})
