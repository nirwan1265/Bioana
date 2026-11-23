# GWAS engines (`rrblup`, `gapit_mlm`, `bioana_emmax`)

Every engine scans the same phenotype: the **line means (BLUEs)** from the
design's `analysis.model`. The genotype term is fixed, and each line's
prediction is averaged over the other fixed factors, as emmeans does. Two-stage
GWAS on BLUEs is the plant-breeding standard; BLUPs would shrink the
phenotypes and count the kinship twice.

| rank | engine | model |
|---|---|---|
| 1 | `rrBLUP::GWAS(P3D = TRUE)` | mixed model, kinship (`A.mat`) random effect, PCs fixed, variances fitted once |
| 2 | `GAPIT::GAPIT(model = "MLM")` | compressed MLM, kinship + PCs; reports allele effects |
| 9 | Bioana EMMAX (fallback) | the same model as rrBLUP P3D, written in Bioana; also gives effects |

`n_pcs = auto`: 0 for lines from a cross (`population.mating` in
`genetics.json` `cross_mating`), 3 otherwise. Markers below `min_maf` (0.05)
are not tested. Adjusted p values are Benjamini-Hochberg.

**Truth (simulated).** Additive and dominance QTN count. A hit is the QTN
marker, or a marker with r2 >= `ld_r2` (0.5) with it on the same chromosome.
Scores: tpr, fdr, sign agreement, estimate r, and line means vs true genetic
values.

## Tutorial: rrBLUP directly

```r
library(rrBLUP)
# geno: data.frame(marker, chrom, pos, <one -1/0/1 column per line>); pheno: data.frame(line, mean)
K <- A.mat(t(as.matrix(geno[, -(1:3)])))
scores <- GWAS(pheno = pheno, geno = geno, K = K, n.PC = 3, min.MAF = 0.05, P3D = TRUE, plot = FALSE)
head(scores[order(-scores$mean), ])       # -log10 p
```

## Tutorial: GAPIT directly

```r
library(GAPIT)
myGAPIT <- GAPIT(Y = data.frame(Taxa = lines, trait = means),
                 GD = data.frame(Taxa = lines, t(dosage012)),
                 GM = data.frame(SNP = ids, Chromosome = chr, Position = pos),
                 model = "MLM", PCA.total = 3, file.output = FALSE)
head(myGAPIT$GWAS)
```

## Through Bioana

```r
bioana_analyze(dir)                                              # rrBLUP if installed
bioana_analyze(dir, params = list(genotypes = list(engine = "gapit_mlm")))
analyze_gwas("genotypes.vcf", line_means, "out/gwas", trait = "mean", sample_id = "line", backend = "rrBLUP")
```

**Checks.** "native and rrBLUP scans agree": the -log10 p scores correlate
> 0.95. GAPIT's adapter finds its output columns by name, because they have
changed between releases. Run `bioana_analyze(dir, params = list(genotypes =
list(engine = "gapit_mlm")))` on a Simitall GWAS folder after installing GAPIT.

**References.** Endelman JB (2011). Ridge regression and other kernels for
genomic selection with R package rrBLUP. *Plant Genome* 4:250. Wang J, Zhang Z
(2021). GAPIT version 3. *Genomics Proteomics Bioinformatics* 19:629. Kang HM
et al. (2010). Variance component model to account for sample structure in
genome-wide association studies. *Nature Genetics* 42:348.
