# Readers: real file formats -> Bioana data sets. Synthetic files in each format.

write_lines_gz <- function(lines, path) { con <- if (grepl("\\.gz$", path)) gzfile(path, "w") else file(path, "w"); writeLines(lines, con); close(con); path }

test_that("formats are detected from names, headers and columns", {
  td <- tempfile("fmt_"); dir.create(file.path(td, "S01"), recursive = TRUE)
  fc <- write_lines_gz(c("# Program:featureCounts v2.0.6; Command:...", "Geneid\tChr\tStart\tEnd\tStrand\tLength\t/x/S01.bam", "g1\tchr1\t1\t100\t+\t100\t5"), file.path(td, "counts.txt"))
  sf <- write_lines_gz(c("Name\tLength\tEffectiveLength\tTPM\tNumReads", "tx1\t1000\t800\t10\t50"), file.path(td, "S01", "quant.sf"))
  st <- write_lines_gz(c("N_unmapped\t10\t10\t10", "N_multimapping\t5\t5\t5", "N_noFeature\t1\t1\t1", "N_ambiguous\t0\t0\t0", "g1\t50\t2\t48"), file.path(td, "S01ReadsPerGene.out.tab"))
  bc <- write_lines_gz("chr1\t100\t100\t50\t5\t5", file.path(td, "S01_bismark_bt2.bismark.cov.gz"))
  vc <- write_lines_gz(c("##fileformat=VCFv4.2", "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tL1"), file.path(td, "geno.vcf"))
  np <- write_lines_gz("chr1\t99\t300\tp1\t10\t.\t5\t6\t7\t50", file.path(td, "S01_peaks.narrowPeak"))
  mx <- write_lines_gz(c("feature_id\tS01", "g1\t4"), file.path(td, "table.tsv"))
  d <- bioana_detect_format(c(fc, sf, st, bc, vc, np, mx))
  expect_equal(d$format, c("featurecounts", "salmon", "star", "bismark_cov", "vcf", "narrowpeak", "matrix"))
  expect_equal(bioana:::.ba_sample_name(c("/a/b/S01.sorted.bam", "S02Aligned.sortedByCoord.out.bam", "S03_bismark_bt2_pe.deduplicated.bismark.cov.gz")),
               c("S01", "S02", "S03"))
})

test_that("featureCounts with BAM-path columns -> feature_counts, analysed by DESeq2, source recorded", {
  skip_if_not_installed("DESeq2")
  d <- rcbd_folder(); m <- read_meta(d); trt <- m$treatment == "treatment_2"
  set.seed(2); G <- 200; mu <- 2^stats::rnorm(G, 7, 1); lfc <- c(rep(2, 20), rep(0, G - 20))
  cnt <- sapply(seq_len(nrow(m)), function(j) stats::rnbinom(G, mu = mu * 2^(lfc * trt[j]), size = 10))
  fc <- data.frame(Geneid = sprintf("g%03d", 1:G), Chr = "chr1", Start = 1, End = 100, Strand = "+", Length = 100, cnt, check.names = FALSE)
  names(fc)[-(1:6)] <- sprintf("/scratch/run1/bam/%s.sorted.bam", m$sample_id)
  f <- file.path(tempdir(), "featurecounts_out.txt")
  writeLines("# Program:featureCounts v2.0.6", f)
  suppressWarnings(utils::write.table(fc, f, sep = "\t", row.names = FALSE, quote = FALSE, append = TRUE))
  out <- utils::capture.output(e <- bioana_manifest_from_files(d, f, assay = "rnaseq"))
  expect_equal(e$counts$data_type, "feature_counts")
  expect_equal(sort(e$counts$source$samples), sort(m$sample_id))
  expect_match(out[1], "featurecounts .* feature_counts, 12 samples")
  utils::capture.output(r <- bioana_analyze(d, only = "counts", secondary = FALSE))
  a <- r$summary$analyses$counts
  expect_equal(a$engine$package, "DESeq2")
  expect_equal(a$source$format, "featurecounts")
  expect_gte(a$significant$treatmenttreatment_2, 15)
})

test_that("a count table becomes peaks for ChIP and IP/input for MeRIP; beta tables become methylation_beta", {
  d <- rcbd_folder(); m <- read_meta(d); s <- m$sample_id
  set.seed(1)
  cnt <- matrix(stats::rpois(50 * 12, 30), 50, dimnames = list(paste0("peak", 1:50), s))
  f1 <- file.path(tempdir(), "chip_counts.tsv"); write_mat(cnt, f1)
  utils::capture.output(e <- bioana_manifest_from_files(d, f1, assay = "chipseq", append = FALSE))
  expect_equal(e$peaks$data_type, "peak_counts_long")
  ipin <- cbind(cnt, cnt); colnames(ipin) <- c(paste0(s, "_IP"), paste0(s, "_input"))
  f2 <- file.path(tempdir(), "merip.tsv"); write_mat(ipin, f2)
  utils::capture.output(e <- bioana_manifest_from_files(d, f2, assay = "merip_seq"))
  expect_equal(e$ip_input$data_type, "ip_input_counts")
  expect_equal(sort(e$ip_input$source$samples), sort(s))
  b <- matrix(stats::runif(40 * 12), 40, dimnames = list(sprintf("cg%05d", 1:40), s))
  f3 <- file.path(tempdir(), "beta.tsv"); write_mat(b, f3)
  utils::capture.output(e <- bioana_manifest_from_files(d, f3))
  expect_equal(e$methylation_array$data_type, "methylation_beta")
  man <- bioana_manifest(d)
  expect_setequal(vapply(man$data, `[[`, "", "id"), c("peaks", "ip_input", "methylation_array"))
  expect_error(bioana_manifest_from_files(d, { f <- file.path(tempdir(), "tpm.tsv"); x <- cnt / 7.3; write_mat(x, f); f }), "neither counts")
})

test_that("Bismark coverage files (one per sample) -> bisulfite_counts that DSS analyses", {
  skip_if_not_installed("DSS")
  d <- rcbd_folder(); m <- read_meta(d); trt <- m$treatment == "treatment_2"
  set.seed(5); pos <- rep(seq(1000, 20000, by = 1000), each = 5) + rep(0:4 * 20, 20)
  td <- tempfile("bis_"); dir.create(td)
  files <- vapply(seq_len(nrow(m)), function(j) {
    cv <- stats::rpois(length(pos), 25) + 1
    p <- ifelse(pos < 6000 & trt[j], 0.8, 0.3)
    me <- stats::rbinom(length(pos), cv, p)
    write_lines_gz(sprintf("chr1\t%d\t%d\t%.2f\t%d\t%d", pos, pos, 100 * me / cv, me, cv - me),
                   file.path(td, paste0(m$sample_id[j], "_bismark_bt2_pe.deduplicated.bismark.cov.gz")))
  }, "")
  utils::capture.output(e <- bioana_manifest_from_files(d, files, append = FALSE))
  expect_equal(e$methylation$data_type, "bisulfite_counts")
  utils::capture.output(r <- bioana_analyze(d, only = "methylation", secondary = FALSE))
  a <- r$summary$analyses$methylation
  expect_equal(a$engine$package, "DSS")
  res <- r$results$methylation
  expect_true(all(res$padj[res$feature_id %in% c("chr1:1001", "chr1:2001", "chr1:3001", "chr1:4001", "chr1:5001")] < 0.05))
})

test_that("Salmon folders go through tximport; STAR strandedness is detected; peak calls give a consensus SAF", {
  skip_if_not_installed("tximport")
  d <- rcbd_folder(); m <- read_meta(d); td <- tempfile("sal_")
  set.seed(3)
  files <- vapply(m$sample_id, function(s) {
    dir.create(file.path(td, s), recursive = TRUE)
    n <- stats::rpois(6, 100)
    utils::write.table(data.frame(Name = paste0("tx", 1:6), Length = 1000, EffectiveLength = 800, TPM = n / sum(n) * 1e6, NumReads = n),
                       file.path(td, s, "quant.sf"), sep = "\t", row.names = FALSE, quote = FALSE)
    file.path(td, s, "quant.sf")
  }, "")
  t2g <- data.frame(tx = paste0("tx", 1:6), gene = rep(c("gA", "gB", "gC"), each = 2))
  utils::capture.output(e <- bioana_manifest_from_files(d, files, assay = "rnaseq", tx2gene = t2g, append = FALSE))
  expect_match(e$counts$source$reader, "^tximport .*lengthScaledTPM")
  cnt <- utils::read.delim(file.path(d, e$counts$files$counts))
  expect_equal(cnt$feature_id, c("gA", "gB", "gC"))
  st <- vapply(m$sample_id[1:3], function(s) write_lines_gz(c("N_unmapped\t1\t1\t1", "N_multimapping\t1\t1\t1", "N_noFeature\t1\t1\t1", "N_ambiguous\t1\t1\t1",
                                                              "g1\t100\t3\t97", "g2\t50\t2\t48"), file.path(td, paste0(s, "ReadsPerGene.out.tab"))), "")
  x <- bioana:::.ba_read_star(st)
  expect_match(x$reader, "reverse")
  expect_equal(unname(x$counts["g1", ]), c(97, 97, 97))
  pk <- c(write_lines_gz(c("chr1\t100\t300\tp\t1\t.\t1\t1\t1\t1", "chr1\t5000\t5200\tp\t1\t.\t1\t1\t1\t1"), file.path(td, "A_peaks.narrowPeak")),
          write_lines_gz(c("chr1\t250\t400\tp\t1\t.\t1\t1\t1\t1"), file.path(td, "B_peaks.narrowPeak")))
  saf <- utils::read.delim(bioana_consensus_peaks(pk, file.path(td, "consensus.saf")))
  expect_equal(nrow(saf), 1)
  expect_equal(c(saf$Start, saf$End), c(101, 400))
  expect_error(bioana_manifest_from_files(d, pk), "no read counts")
})
