# A Biodes-style GWAS folder: n_lines inbred lines x n_blocks blocks (RCBD), lines as the treatment factor,
# a random 0/2 VCF with no LD, and additive QTN with known effects (-1/0/1 coding, as Simitall).
gwas_folder <- function(n_lines = 80, n_blocks = 2, n_markers = 400, qtn_effect = c(1.2, -1, 0.8, -0.6),
                        sim_truth = TRUE, seed = 11) {
  set.seed(seed)
  d <- file.path(tempfile("bioana_gwas_"), "design"); dir.create(d, recursive = TRUE)
  lines <- sprintf("L%03d", seq_len(n_lines))
  spec <- list(handoff_version = "0.4.0", system = "field", organism = "Sorghum bicolor",
               factors = list(list(name = "genotype", n_levels = n_lines)),
               design = list(chosen = "rcbd", reps = n_blocks), assays = list("gwas"),
               primary_outcome = list(name = "height", type = "continuous"),
               analysis = list(model = "y ~ treatment + block",
                               terms = list(treatment = "genotype", block = "block", layout_columns = list("genotype", "block"))))
  jsonlite::write_json(spec, file.path(d, "study_spec.json"), auto_unbox = TRUE, pretty = TRUE)
  lay <- do.call(rbind, lapply(seq_len(n_blocks), function(b) data.frame(plot = b * 1000 + seq_len(n_lines), block = b, genotype = sample(lines))))
  utils::write.csv(lay, file.path(d, "field_layout.csv"), row.names = FALSE)
  freq <- stats::runif(n_markers, 0.2, 0.8)
  dos <- 2L * (matrix(stats::runif(n_markers * n_lines), n_markers) < freq)
  chr <- sort(sample(1:4, n_markers, TRUE)); pos <- unlist(lapply(1:4, function(k) sort(sample.int(1e6, sum(chr == k)))))
  ids <- paste0("chr", chr, "_", pos)
  dimnames(dos) <- list(ids, lines)
  gt <- matrix(c("0/0", "0/1", "1/1")[dos + 1L], nrow(dos))
  writeLines(c("##fileformat=VCFv4.2", paste(c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", lines), collapse = "\t"),
               paste(paste0("chr", chr), pos, ids, "A", "G", ".", "PASS", ".", "GT", apply(gt, 1, paste, collapse = "\t"), sep = "\t")),
             file.path(d, "genotypes.vcf"))
  qi <- sample(n_markers, length(qtn_effect))
  g <- drop(crossprod(dos[qi, , drop = FALSE] - 1, qtn_effect)); names(g) <- lines
  meta <- data.frame(lay, sample_id = sprintf("S%04d", seq_len(nrow(lay))))
  meta$primary_outcome <- 20 + g[meta$genotype] + 0.5 * meta$block + stats::rnorm(nrow(meta), 0, 1)
  utils::write.table(meta, file.path(d, "sample_metadata_template.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
  dir.create(file.path(d, "genetic_truth"))
  utils::write.table(data.frame(marker_id = ids[qi], seqname = paste0("chr", chr[qi]), pos = pos[qi], type = "additive", effect = qtn_effect),
                     file.path(d, "genetic_truth", "qtn_truth.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
  utils::write.table(data.frame(line = lines, genetic_value = g), file.path(d, "genetic_truth", "genetic_values.tsv"),
                     sep = "\t", row.names = FALSE, quote = FALSE)
  if (sim_truth) {
    jsonlite::write_json(list(genetics = list(genotype_term = "treatment", genetic_values = "genetic_truth/genetic_values.tsv")),
                         file.path(d, "sim_truth.json"), auto_unbox = TRUE)
  }
  bioana_write_manifest(d, list(id = "genetics", data_type = "genotypes", assay = "genotypes", files = list(vcf = "genotypes.vcf"),
                                truth = list(kind = "qtn", file = "genetic_truth/qtn_truth.tsv")))
  list(dir = d, g = g, dosage = dos, qtn = ids[qi])
}

test_that("genotypes: line means from the design's model, native GWAS, QTN scored", {
  f <- gwas_folder()
  out <- utils::capture.output(r <- bioana_analyze(f$dir, params = list(genotypes = list(engine = "bioana_emmax"))))
  expect_null(r$summary$not_analysed$genetics)
  a <- r$summary$analyses$genetics
  expect_equal(a$genotype_term, "treatment")
  expect_equal(a$n_lines, 80)
  expect_match(a$model_fitted, "treatment \\+ block")
  expect_equal(a$n_pcs, 3L)
  lm <- utils::read.delim(file.path(f$dir, "bioana", "genetics", "line_means.tsv"), colClasses = c(line = "character"))
  expect_gt(stats::cor(lm$mean, f$g[lm$line]), 0.9)
  expect_gt(a$line_mean_cor, 0.9)
  rec <- a$recovery
  expect_equal(rec$qtn, 4)
  expect_equal(rec$testable, 4)
  expect_gte(rec$true_positives, 2)          # the two large QTN; 0.8 and -0.6 are weak with 80 lines x 2 blocks
  expect_lt(rec$fdr, 0.35)
  expect_equal(rec$sign_agreement, 1)
  expect_gt(rec$estimate_cor, 0.8)
  q <- utils::read.delim(file.path(f$dir, "bioana", "genetics", "qtn_recovery.tsv"))
  expect_true(all(c("marker_id", "type", "effect", "tested", "recovered", "estimate", "p", "padj") %in% names(q)))
  expect_true(all(q$recovered[abs(q$effect) >= 1]))
  expect_true(any(grepl("\\[genetics\\].* 80 lines x .* QTN found", out)))
})

test_that("the genotype term is found from the rules when the folder does not name it", {
  f <- gwas_folder(n_lines = 30, n_markers = 60, sim_truth = FALSE)
  design <- bioana:::.ba_read_biodes(f$dir)
  expect_equal(bioana:::.ba_genotype_term(list(), design), "treatment")
  expect_equal(bioana:::.ba_genotype_term(list(genotype_term = "treatment"), design), "treatment")
  expect_error(bioana:::.ba_genotype_term(list(genotype_term = "env"), design), "not in analysis.model")
})

test_that("line means average over blocks and match the cell means in a balanced design", {
  f <- gwas_folder(n_lines = 20, n_markers = 40)
  design <- bioana:::.ba_read_biodes(f$dir)
  lm <- bioana:::.ba_line_means(design, "treatment", list(), list())$means
  meta <- design$meta
  raw <- tapply(as.numeric(meta$primary_outcome), meta$genotype, mean)
  expect_equal(unname(lm$mean), as.vector(raw[lm$line]), tolerance = 1e-8)
  expect_true(all(lm$n_units == 2))
})

test_that("benchmark_gwas: a hit in LD with a QTN counts with ld_r2, not without", {
  set.seed(5)
  n <- 100
  hap <- stats::rbinom(n, 1, 0.5) * 2
  dos <- rbind(qtn = hap, linked = ifelse(stats::runif(n) < 0.95, hap, 2 - hap), far = stats::rbinom(n, 1, 0.5) * 2)
  res <- data.frame(marker_id = c("qtn", "linked", "far"), seqname = "chr1", pos = c(1000, 50000, 900000),
                    q_value = c(0.5, 0.001, 0.5))
  truth <- data.frame(marker_id = "qtn", seqname = "chr1", pos = 1000)
  out <- tempfile()
  plain <- benchmark_gwas(res, truth, out)
  expect_equal(plain$metrics$recall, 0)
  expect_equal(plain$metrics$false_positive_markers, 1)
  ld <- benchmark_gwas(res, truth, out, ld_r2 = 0.5, genotypes = dos)
  expect_equal(ld$metrics$recall, 1)
  expect_equal(ld$metrics$false_positive_markers, 0)
})

test_that("native and rrBLUP scans agree", {
  skip_if_not_installed("rrBLUP")
  f <- gwas_folder(n_lines = 60, n_markers = 200)
  ph <- data.frame(line = names(f$g), y = f$g + stats::rnorm(60, 0, 0.7))
  a <- analyze_gwas(file.path(f$dir, "genotypes.vcf"), ph, tempfile(), trait = "y", sample_id = "line", n_pcs = 2, backend = "native")$results
  b <- analyze_gwas(file.path(f$dir, "genotypes.vcf"), ph, tempfile(), trait = "y", sample_id = "line", n_pcs = 2, backend = "rrBLUP")$results
  m <- match(a$marker_id, b$marker_id)
  expect_gt(stats::cor(a$score, b$score[m], use = "complete.obs"), 0.95)
})

test_that("a Simitall GWAS design runs end to end", {
  skip_if_not_installed("simitall")
  f <- gwas_folder(n_lines = 60, n_markers = 10)
  unlink(file.path(f$dir, c("genotypes.vcf", "sim_truth.json", "data_manifest.json")))
  unlink(file.path(f$dir, "genetic_truth"), recursive = TRUE)
  meta <- read_meta(f$dir); meta$primary_outcome <- ""; write_meta(meta, f$dir)
  utils::capture.output(simitall::simitall_from_biodes(f$dir, seed = 3,
    params = list(suite_genetics = list(heritability = 0.8), suite_genotype_panel = list(n_markers = 600))))
  man <- bioana_manifest(f$dir)
  expect_true("genotypes" %in% vapply(man$data, `[[`, "", "data_type"))
  utils::capture.output(r <- bioana_analyze(f$dir, only = "genotypes", params = list(genotypes = list(engine = "bioana_emmax"))))
  a <- r$summary$analyses$genetics
  expect_equal(a$n_lines, 60)
  expect_gt(a$line_mean_cor, 0.7)
  expect_gt(a$recovery$tpr, 0)
  expect_gte(a$recovery$sign_agreement, 0.5)
})
