# dream (`dream`, package variancePartition)

**What it fits.** A precision-weighted **linear mixed model** per feature. The
weights come from `voomWithDreamWeights`, the model is fitted with lme4 inside
`dream`, and the random terms are exactly those of `analysis.model`. Each
coefficient gets a moderated t with Satterthwaite degrees of freedom, after
`eBayes`, with Benjamini-Hochberg adjustment (`topTable`).

**When Bioana uses it.** It is rank 1 for counts, peaks and windows whenever
`analysis.model` has random terms. Examples: the same plant or animal sampled
over time (`(1 | subject)`), split plots (`(1 | block:main_plot)`), or units
nested in blocks. Use it instead of DESeq2 when the units are not independent.
DESeq2's design formula has no random terms, and treating repeated samples as
replicates inflates false positives.

**Input.** Raw integer counts; the unit metadata with every model column.

**Assumptions.** Normal random effects. log2-CPM is roughly normal given the
voom/dream weights.

## Tutorial: dream directly

```r
library(variancePartition); library(edgeR)
y <- normLibSizes(DGEList(counts))
form <- ~ treatment * time + (1 | subject)     # same individuals at every time
vobj <- voomWithDreamWeights(y, form, metadata)
fit <- dream(vobj, form, metadata)
fit <- eBayes(fit)
topTable(fit, coef = "treatmenttreatment_2:timet2", number = Inf)
```

## The same through Bioana

Biodes writes `analysis.model = "y ~ treatment * time + (1 | subject)"`, and
`bioana_analyze(dir)` picks dream on its own. Without variancePartition,
limma-voom with `duplicateCorrelation` takes over (one random term). Without
limma, the native fallback runs; it drops the random terms and says so.

**Gold check.** "dream adapter = variancePartition::dream called directly":
equal t statistics. This test is skipped when variancePartition is not
installed.

**Reference.** Hoffman GE, Roussos P (2021). Dream: powerful differential
expression analysis for repeated measures designs. *Bioinformatics* 37:192.
Vignette: "dream: differential expression testing with linear mixed models for
repeated measures" (variancePartition).
