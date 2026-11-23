# edgeR (`edger`, `edger_ip_input`)

## `edger`: counts, quasi-likelihood

**What it fits.** TMM normalisation (`normLibSizes`, or `calcNormFactors` in
older edgeR), NB dispersions (`estimateDisp`), then a quasi-likelihood NB GLM
(`glmQLFit`). Each coefficient gets a QL F test (`glmQLFTest`), with
Benjamini-Hochberg adjustment (`topTags`).

**When Bioana uses it.**
- Rank 1 for `peak_counts_long` (ChIP, CUT&Tag) and `read_positions` (ATAC
  windows) in fixed designs. csaw and DiffBind fit edgeR or DESeq2 on peak or
  window counts, so this is their statistical core; adapters for those two
  packages are planned.
- Rank 2 for `feature_counts`, so it is the default secondary to DESeq2.

**Input.** Raw integer counts per feature and unit.

**Assumptions.** NB counts with QL dispersion; independent units; fixed effects.

**Estimate.** `logFC` (log2). edgeR reports no standard error for QL tests, so
`se` is `NA`. The `statistic` column holds the F value.

### Tutorial: edgeR directly

```r
library(edgeR)
y <- DGEList(counts)
y <- normLibSizes(y)                       # calcNormFactors(y) in edgeR < 4
X <- model.matrix(~ treatment + block, data = coldata)
y <- estimateDisp(y, X)
fit <- glmQLFit(y, X)
qlf <- glmQLFTest(fit, coef = "treatmenttreatment_2")
topTags(qlf)                               # logFC, logCPM, F, PValue, FDR
```

### The same through Bioana

```r
bioana_analyze(dir, params = list(feature_counts = list(engine = "edger")))
bioana_analyze(dir, secondary = TRUE)      # DESeq2 primary, edgeR secondary, agreement reported
```

**Gold check.** "edgeR adapter = edgeR QL pipeline called directly": equal
`logFC`, `PValue` and `FDR`.

## `edger_ip_input`: m6A MeRIP-seq, MeDIP, RIP

**What it fits.** One NB GLM on the IP and input libraries together. Design:
one effect per unit, an `ip` indicator, and `ip x` the design's fixed terms.
The test is the **IP × treatment interaction**: the change in enrichment.

- A change in RNA expression raises IP and input together, so it is absorbed
  and is not called as methylation.
- Library offsets are as in the native method: input libraries by library
  size, and IP libraries by their median IP/input ratio, which absorbs IP
  efficiency per unit.

This is the interaction model used by exomePeak2 and RADAR for differential
modification. exomePeak2 itself needs BAM files; its adapter comes with the
raw-data readers.

**Input.** IP and input counts per window for the same units. Windows with
mean input below `min_input` are dropped.

**Estimate.** log2 change in IP/input enrichment.

### Tutorial: the model with edgeR directly

```r
units <- colnames(ip)
U <- diag(length(units))
Xf <- model.matrix(~ treatment + block, coldata)[, -1]
X2 <- rbind(cbind(U, ip = 1, Xf), cbind(U, ip = 0, Xf * 0))
colnames(X2) <- c(paste0("unit_", units), "ip", paste0("ip:", colnames(Xf)))
y <- DGEList(cbind(ip, input))
fit <- glmQLFit(estimateDisp(y, X2), X2)
topTags(glmQLFTest(fit, coef = "ip:treatmenttreatment_2"))
```

**Check.** "edgeR IP x input: enrichment changes are found, expression changes
are not": windows with an IP-only change are found (>80%). Windows where
expression rises in both libraries are not (<15%).

**References.** Chen Y, Lun ATL, Smyth GK (2016). From reads to genes to
pathways: differential expression analysis of RNA-Seq experiments using
Rsubread and the edgeR quasi-likelihood pipeline. *F1000Research* 5:1438.
edgeR User's Guide (Bioconductor).
