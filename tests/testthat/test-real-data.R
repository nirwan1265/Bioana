# Real data from published packages: each run through Bioana and checked against the package called
# directly or against its own published results. Data: github.com/nirwan1265/bioana-testdata (README there).

test_that("airway BAMs (airway package): featureCounts counts equal the package's gene counts", {
  skip_if_not_installed("Rsubread")
  td <- skip_without_testdata("rnaseq_airway")
  d <- bioana_example("airway")
  utils::capture.output(e <- bioana_manifest_from_files(d, list.files(td, "_subset\\.bam$", full.names = TRUE), assay = "rnaseq",
                                                        annotation = file.path(td, "Homo_sapiens.GRCh37.75_subset.gtf")))
  expect_equal(e$counts2$data_type, "feature_counts")
  expect_match(e$counts2$source$reader, "Rsubread::featureCounts .*isPairedEnd = TRUE")
  man <- bioana_manifest(d)
  bc <- as.matrix(utils::read.delim(file.path(d, e$counts2$files$counts), row.names = 1, check.names = FALSE))
  pc <- as.matrix(utils::read.delim(file.path(d, man$data[[1]]$files$counts), row.names = 1, check.names = FALSE))
  g <- intersect(rownames(bc)[rowSums(bc) > 0], rownames(pc)); s <- colnames(pc)
  expect_gte(length(g), 10)
  expect_equal(stats::median(bc[g, s] / pmax(pc[g, s], 1)), 1)
  expect_gt(stats::cor(as.vector(log1p(bc[g, s])), as.vector(log1p(pc[g, s]))), 0.99)
})

test_that("DiffBind tamoxifen ChIP (chr18 1-20 Mb): Bioana's calls are DiffBind's calls, same direction", {
  skip_if_not_installed("Rsubread"); skip_if_not_installed("DESeq2")
  td <- skip_without_testdata("chip_diffbind_chr18_20mb")
  ss <- utils::read.csv(text = gsub("\r", "\n", readChar(file.path(td, "tamoxifen.csv"), 1e5)), stringsAsFactors = FALSE)
  env <- new.env(); load(file.path(td, "dba_norm_comp.RData"), envir = env); o <- env$DBAobject
  s <- data.frame(o$merged); names(s) <- c("chr", "start", "end"); s$chr <- o$chrmap[s$chr]; s <- s[s$end <= 2e7, ]
  meta <- data.frame(sample_id = ss$SampleID, Condition = ss$Condition, Tissue = ss$Tissue)
  d <- real_design(meta, "Condition", c("Responsive", "Resistant"), "chipseq", "y ~ treatment",
                   list(treatment = "Condition", layout_columns = list("Condition")))
  saf <- file.path(d, "diffbind_sites.saf")
  utils::write.table(data.frame(GeneID = sprintf("site_%d_%d", s$start, s$end), Chr = s$chr, Start = s$start, End = s$end, Strand = "+"),
                     saf, sep = "\t", row.names = FALSE, quote = FALSE)
  ctl <- unique(ss$bamControl)
  libs <- rbind(data.frame(file = ss$bamReads, sample_id = ss$SampleID, library = "chip"),
                data.frame(file = ctl, sample_id = basename(ctl), library = "control"))
  utils::capture.output(e <- bioana_manifest_from_files(d, file.path(td, c(ss$bamReads, ctl)), assay = "chipseq", annotation = saf, libraries = libs))
  expect_equal(e$peaks$data_type, "peak_counts_long")
  expect_match(e$peaks$source$reader, "reads extended to 125 bp fragments")
  utils::capture.output(r <- bioana_analyze(d, secondary = TRUE))
  expect_equal(r$summary$analyses$peaks$engine$package, "edgeR")
  for (res in list(r$results$peaks, utils::read.delim(file.path(d, "bioana", "peaks", "results_deseq2.tsv")))) {
    db <- data.frame(o$peaks[[3]]); key <- sprintf("site_%d_%d", db$Start, db$End)     # DiffBind DESeq2, RLE on reads in peaks
    call <- res$feature_id[!is.na(res$padj) & res$padj < 0.05]
    expect_gt(length(call), 30)
    expect_gt(mean(call %in% key), 0.9)                      # precision against DiffBind's calls
    sc <- db$Score[match(call, key)]; ok <- !is.na(sc)
    expect_equal(mean(sign(sc[ok]) == sign(res$estimate[match(call, res$feature_id)][ok])), 1)
  }
})

test_that("maize 282 panel (TASSEL): line means over 2 locations, Bioana's GWAS = rrBLUP called directly", {
  skip_if_not_installed("rrBLUP")
  td <- skip_without_testdata("maize_tassel")
  ph <- utils::read.delim(file.path(td, "mdp_phenotype.txt"), skip = 2, stringsAsFactors = FALSE, na.strings = c("NaN", "-999"))
  ph$sample_id <- paste(ph$Taxa, ph$location, sep = "_")
  meta <- data.frame(ph[, c("sample_id", "Taxa", "location")], primary_outcome = ph$dpoll)
  d <- real_design(meta, "Taxa", character(), "gwas", "y ~ treatment + block",
                   list(treatment = "Taxa", block = "location", layout_columns = list("Taxa", "location")),
                   list(primary_outcome = list(name = "dpoll", type = "continuous"), genetics = list(unit = "line")))
  utils::capture.output(e <- bioana_manifest_from_files(d, file.path(td, "mdp_genotype.hmp.txt")))
  expect_match(e$genotypes$source$reader, "3093 of 3093 biallelic SNPs")
  m <- jsonlite::read_json(file.path(d, "data_manifest.json")); m$data[[1]]$genotype_term <- "treatment"
  jsonlite::write_json(m, file.path(d, "data_manifest.json"), auto_unbox = TRUE)
  utils::capture.output(r <- bioana_analyze(d, secondary = FALSE))
  a <- r$summary$analyses$genotypes
  expect_equal(a$engine$package, "rrBLUP")
  lm <- utils::read.delim(file.path(d, "bioana", "genotypes", "line_means.tsv"), colClasses = c(line = "character"))
  tr <- utils::read.delim(file.path(td, "mdp_traits.txt"), skip = 1, header = FALSE, na.strings = "-999",
                          colClasses = c("character", rep("numeric", 3)))
  expect_gt(stats::cor(lm$mean, tr$V3[match(lm$line, tr$V1)], use = "complete.obs"), 0.95)   # TASSEL's own dpoll values
  g <- bioana:::.bioana_read_genotypes(file.path(d, m$data[[1]]$files$vcf))
  lm <- lm[lm$line %in% colnames(g$genotype) & is.finite(lm$mean), ]
  X <- t(g$genotype[, lm$line] - 1)
  geno <- data.frame(marker = rownames(g$genotype), chrom = g$variants$seqname, pos = g$variants$pos, t(X), check.names = FALSE)
  utils::capture.output(sc <- rrBLUP::GWAS(pheno = data.frame(line = lm$line, mean = lm$mean), geno = geno, K = rrBLUP::A.mat(X),
                                           n.PC = a$n_pcs, min.MAF = 0.05, P3D = TRUE, plot = FALSE))
  res <- r$results$genotypes; k <- match(res$marker_id, sc$marker)
  expect_lt(max(abs(-log10(res$p) - sc$mean[k]), na.rm = TRUE), 1e-6)
})

test_that("DSS example counts (DSS package): Bioana at site level = DSS called directly", {
  skip_if_not_installed("DSS")
  f <- system.file("extdata", c("cond1_1.txt", "cond1_2.txt", "cond2_1.txt", "cond2_2.txt"), package = "DSS")
  meta <- data.frame(sample_id = c("cond1_1", "cond1_2", "cond2_1", "cond2_2"), group = c("c1", "c1", "c2", "c2"))
  d <- real_design(meta, "group", c("c1", "c2"), "wgbs", "y ~ treatment", list(treatment = "group", layout_columns = list("group")))
  utils::capture.output(e <- bioana_manifest_from_files(d, f))
  expect_equal(e$methylation$source$format, "dss_txt")
  utils::capture.output(r <- bioana_analyze(d, params = list(methylation = list(window_bp = 1, min_coverage = 1)), secondary = FALSE))
  res <- r$results$methylation
  dat <- lapply(f, utils::read.table, header = TRUE)
  bs <- DSS::makeBSseqData(dat, meta$sample_id)
  keep <- rowSums(bsseq::getCoverage(bs, type = "Cov") >= 1) == 4
  utils::capture.output(fit <- DSS::DMLfit.multiFactor(bs[keep], design = data.frame(group = factor(meta$group)), formula = ~ group))
  t2 <- DSS::DMLtest.multiFactor(fit, coef = "groupc2")
  k <- match(paste0(t2$chr, ":", t2$pos), res$feature_id)
  expect_false(anyNA(k))
  expect_equal(res$statistic[k], t2$stat, tolerance = 1e-8)
  expect_equal(sum(res$padj < 0.05, na.rm = TRUE), sum(t2$fdrs < 0.05, na.rm = TRUE))
})

test_that("bsseq's real Bismark output files are recognised and read", {
  skip_if_not_installed("bsseq")
  bf <- system.file("extdata", c("test_data.fastq_bismark.bismark.cov.gz", "test_data.cytosineReport.gz"), package = "bsseq")
  expect_equal(bioana_detect_format(bf)$format, c("bismark_cov", "bismark_report"))
  x <- bioana:::.ba_read_bismark(bf[1])
  expect_gt(nrow(x$sites), 1000)
  y <- bioana:::.ba_read_bismark(bf[2])
  expect_true("context" %in% names(y$sites))
})

test_that("exomePeak2 MeRIP BAMs: exon windows, IP paired with input, the IP x treatment model runs", {
  skip_if_not_installed("Rsubread"); skip_if_not_installed("edgeR")
  td <- skip_without_testdata("merip_exomepeak2")
  libs <- data.frame(file = c(paste0("IP", 1:4, ".bam"), paste0("Input", 1:3, ".bam"), "treated_IP1.bam", "treated_Input1.bam"),
                     sample_id = c("c1", "c2", "c3", "c4", "c1", "c2", "c3", "t1", "t1"),
                     library = c(rep("ip", 4), rep("input", 3), "ip", "input"))
  meta <- data.frame(sample_id = c("c1", "c2", "c3", "t1"), group = c("control", "control", "control", "treated"))
  d <- real_design(meta, "group", c("control", "treated"), "merip_seq", "y ~ treatment", list(treatment = "group", layout_columns = list("group")))
  utils::capture.output(e <- bioana_manifest_from_files(d, file.path(td, libs$file), assay = "merip_seq",
                                                        annotation = file.path(td, "example.gtf"), libraries = libs))
  expect_equal(e$ip_input$data_type, "ip_input_counts")
  expect_match(e$ip_input$source$warnings, "c4")                       # IP4 has no input: left out
  utils::capture.output(r <- bioana_analyze(d, params = list(ip_input = list(min_input = 1)), secondary = FALSE))
  expect_equal(r$summary$analyses$ip_input$engine$engine, "edger_ip_input")
  expect_gt(max(r$results$ip_input$mean_log2_enrichment), 3)            # m6A peaks: IP well above input
})
