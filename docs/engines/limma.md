# limma (`limma_voom`, `limma`)

## `limma_voom`: counts

**What it fits.** TMM-normalised counts are turned into log2-CPM. `voom`
estimates the mean-variance trend and gives every observation a precision
weight. Then a weighted linear model per feature (`lmFit`), empirical-Bayes
moderated variances (`eBayes`), and a moderated t per coefficient, with
Benjamini-Hochberg adjustment (`topTable`).

**One random grouping term** (for example `(1 | subject)` or `(1 | block)`) is
modelled as a correlation within groups (`duplicateCorrelation`), with two
passes of voom as the limma User's Guide recommends. With more than one random
term, use dream.

**When Bioana uses it.** Rank 3 for counts in fixed designs. Rank 2 for counts
with random terms, when dream (variancePartition) is not installed or as the
secondary.

### Tutorial: limma-voom directly

```r
library(limma); library(edgeR)
y <- normLibSizes(DGEList(counts))
X <- model.matrix(~ treatment, coldata)
v <- voom(y, X)
# same subjects measured repeatedly: correlation within subject
cor <- duplicateCorrelation(v, X, block = coldata$subject)$consensus.correlation
v <- voom(y, X, block = coldata$subject, correlation = cor)
cor <- duplicateCorrelation(v, X, block = coldata$subject)$consensus.correlation
fit <- eBayes(lmFit(v, X, block = coldata$subject, correlation = cor))
topTable(fit, coef = "treatmenttreatment_2")
```

## `limma`: values on a linear scale

**What it fits.** `lmFit` + `eBayes` on values that are already roughly normal:
methylation-array M-values, bisulfite window M-values (when the design has a
random term, which DSS cannot model), and log2(IP/input) ratios.
`duplicateCorrelation` handles one random grouping term.

**When Bioana uses it.**
- Rank 1 for `methylation_beta`, the minfi/limma convention.
- Rank 1 for bisulfite counts and IP/input data with a random term.
- Rank 2 for bisulfite counts and IP/input data otherwise.

### Tutorial: limma on M-values

```r
M <- log2(beta / (1 - beta))
X <- model.matrix(~ treatment + block, coldata)
fit <- eBayes(lmFit(M, X))
topTable(fit, coef = "treatmenttreatment_2", number = Inf)
```

### The same through Bioana

```r
bioana_analyze(dir)   # the rules pick limma_voom / limma where they rank first
bioana_analyze(dir, params = list(feature_counts = list(engine = "limma_voom")))
```

**Gold checks.** All three require equal t statistics and adjusted p values:
- "limma-voom adapter = voom/lmFit/eBayes called directly"
- "limma with one random term = duplicateCorrelation called directly"
- "limma on values = lmFit/eBayes"

**References.** Law CW, Chen Y, Shi W, Smyth GK (2014). voom: precision weights
unlock linear model analysis tools for RNA-seq read counts. *Genome Biology*
15:R29. Ritchie ME et al. (2015). limma powers differential expression
analyses for RNA-sequencing and microarray studies. *Nucleic Acids Research*
43:e47. Smyth GK (2004). *Stat Appl Genet Mol Biol* 3:3.
