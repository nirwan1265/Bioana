# DSS (`dss`)

**What it fits.** A beta-binomial regression with arcsine link for a general
design (`DMLfit.multiFactor`), on methylated and total reads, followed by a
Wald test per coefficient (`DMLtest.multiFactor`). DSS's `fdrs` are
Benjamini-Hochberg.

**When Bioana uses it.** It is rank 1 for `bisulfite_counts` (WGBS, RRBS,
EM-seq) in fixed designs. With random terms limma on M-values takes over,
because DSS has no random effects.

**Input.** Bioana pools sites into windows (`window_bp`, 1000, per context)
and gives DSS the pooled methylated and total reads per window and unit. Each
window counts as one locus. That keeps the tested unit (the window) the same as
in Simitall's truth and the other engines. To fit sites instead, set
`window_bp = 1`. Windows below `min_coverage` reads in any unit are dropped.

**Assumptions.** Read counts are beta-binomial, with overdispersion shared
through DSS's prior.

**Estimate.** The coefficient on DSS's arcsine scale, with its standard error.
Its sign is the direction of the methylation change. For a share-of-reads
difference, see the `mean_methylation` column.

## Tutorial: DSS directly

```r
library(DSS); library(bsseq)
bs <- BSseq(chr = chr, pos = pos, M = methylated, Cov = coverage, sampleNames = colnames(coverage))
design <- data.frame(treatment = factor(...), block = factor(...))
fit <- DMLfit.multiFactor(bs, design = design, formula = ~ treatment + block)
res <- DMLtest.multiFactor(fit, coef = "treatmenttreatment_2")
head(res[order(res$pvals), ])     # chr, pos, stat, pvals, fdrs
```

## The same through Bioana

```r
bioana_write_manifest(dir, list(id = "meth", data_type = "bisulfite_counts",
  files = list(sites = "sites.tsv", methylated = "methylated.tsv", coverage = "coverage.tsv")))
bioana_analyze(dir)
```

**Gold check.** "DSS adapter = DMLfit.multiFactor / DMLtest.multiFactor called
directly": equal `stat` and `fdrs`, same feature order, planted changes found.

**Reference.** Park Y, Wu H (2016). Differential methylation analysis for BS-seq
data under general experimental design. *Bioinformatics* 32:1446. DSS
vignette, section on general experimental design.
