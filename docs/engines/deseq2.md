# DESeq2 (`deseq2`)

**What it fits.** A negative-binomial GLM per gene, using median-of-ratios
size factors. Gene-wise dispersions are shrunk toward a fitted mean-dispersion
trend, and each coefficient gets a Wald test. Adjusted p values are
Benjamini-Hochberg, with independent filtering on the mean count (the default
of `results()`). Genes removed by that filter get `padj = NA`.

**When Bioana uses it.** It is rank 1 for `feature_counts` when
`analysis.model` has no random terms, and rank 2 for peak and window counts.
It is never used with random terms (DESeq2 has none); dream or limma-voom take
those designs.

**Input.** Raw integer counts, one column per experimental unit (`sample_id`),
and one library per unit. Do not give it CPM, TPM or normalized values. Bioana
first drops features below `min_cpm` (1) in most samples, the same filter for
every engine.

**Assumptions.** The counts follow NB mean-variance; units are independent;
only fixed effects.

**Estimate.** log2 fold change (`log2FoldChange`), standard error (`lfcSE`).
This is the unshrunken MLE, so it can be compared with the truth.

## Tutorial: DESeq2 directly

```r
library(DESeq2)
# counts: genes x samples integer matrix; coldata: one row per sample, with factors
coldata <- data.frame(treatment = factor(rep(c("treatment_1", "treatment_2"), 6)),
                      block = factor(rep(1:6, each = 2)))
dds <- DESeqDataSetFromMatrix(countData = counts, colData = coldata, design = ~ treatment + block)
dds <- DESeq(dds)
resultsNames(dds)        # "Intercept" "treatment_treatment_2_vs_treatment_1" "block_2_vs_1" ...
res <- results(dds, name = "treatment_treatment_2_vs_treatment_1", alpha = 0.05)
summary(res)
head(res[order(res$padj), ])
```

## The same through Bioana

```r
bioana_write_manifest(dir, list(id = "rna", data_type = "feature_counts", assay = "rnaseq",
                                files = list(counts = "counts.tsv")))
r <- bioana_analyze(dir)
r$summary$analyses$rna$engine     # package, version, exact call, model, test, ...
head(r$results$rna)               # feature_id, contrast, estimate, se, statistic, p, padj
```

Bioana passes DESeq2 the design's model matrix (`DESeq(dds, full = X)`) and
does not rebuild a formula. The coefficient names (`treatmenttreatment_2`) are
then the ones in Simitall's truth files. The parametrisation is the same as
`design = ~ treatment + block`.

**Gold check.** `test-engines.R`, "DESeq2 adapter = DESeq2 called directly".
It runs the tutorial code above on the same counts and requires equal
`log2FoldChange`, `lfcSE`, `pvalue` and `padj`.

**Reference.** Love MI, Huber W, Anders S (2014). Moderated estimation of fold
change and dispersion for RNA-seq data with DESeq2. *Genome Biology* 15:550.
Vignette: "Analyzing RNA-seq data with DESeq2" (Bioconductor).
