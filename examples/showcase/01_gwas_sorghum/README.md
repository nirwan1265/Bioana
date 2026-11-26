# 01 GWAS of plant height in a sorghum panel

**Question** (to Biodes): "GWAS of plant height in 300 sorghum accessions grown in an alpha lattice with 2
replicates at one location in one season."

Design: [Biodes](https://github.com/nirwan1265/BioDesignGPT/tree/main/examples/showcase/01_gwas_sorghum) ·
data and truth: [Simitall](https://github.com/nirwan1265/simitall/tree/main/examples/showcase/01_gwas_sorghum)
(300 inbred lines x 1,905 markers with population structure, 10 additive QTNs, heritability 0.6) ·
figures: Biopaint (planned).

## Model

1. Line means: `y ~ treatment + rep + (1 | rep:block)` (lmerTest), one BLUE per accession.
2. Scan: `rrBLUP::GWAS` on the BLUEs, kinship (`A.mat`) as a random effect, 3 principal components as fixed
   covariates, markers with MAF >= 0.05, Benjamini-Hochberg FDR 0.05.

## Result

6 of 1,905 markers are significant at FDR 0.05 (`results/results.tsv`).

## Against the truth

| | |
|---|---|
| QTNs found (FDR 0.05) | 6 of 10 (recall 0.60) |
| False positives | 0 (precision 1.00) |
| Line BLUEs vs true genetic values | r = 0.858 |

The four QTNs missed are four of the five smallest (planted effects 5.79, 5.21, 4.69 and 3.80; p between
0.001 and 0.12). The fifth small one, effect 4.22 on chr4, was found. Every QTN found has its planted sign.
Per QTN (`results/qtn_recovery.tsv`):

| Marker | Planted effect | MAF | Estimate | adj. p | Found |
|---|---|---|---|---|---|
| chr3_48370 | 9.81 | 0.18 | 9.56 | 3.2e-7 | yes |
| chr2_724691 | 8.83 | 0.33 | 8.43 | 5.0e-8 | yes |
| chr3_46778 | 7.94 | 0.25 | 6.65 | 0.0010 | yes |
| chr4_552758 | 7.15 | 0.48 | 5.18 | 0.0069 | yes |
| chr2_238608 | 6.43 | 0.41 | 5.23 | 0.0098 | yes |
| chr1_482382 | 5.79 | 0.21 | 4.66 | 0.50 | no |
| chr4_433883 | 5.21 | 0.35 | 4.48 | 0.25 | no |
| chr3_930637 | 4.69 | 0.30 | 3.63 | 0.69 | no |
| chr4_423997 | 4.22 | 0.48 | 5.92 | 0.0012 | yes |
| chr3_975805 | 3.80 | 0.34 | 2.33 | 1.00 | no |

The outcome entry in `report.txt` is the line-mean model: its coefficients cover the truth (298 of 299 in their
CI). The study goal (association) is scored on the scan.

## results/

`results.tsv` (marker scan: p, adjusted p, MAF), `qtn_recovery.tsv`, `gwas.benchmark_metrics.tsv`
(precision, recall, F1), `line_means.tsv` (BLUEs), `report.txt`, `bioana_summary.json`, `truth_recovery.json`.
For Biopaint: a Manhattan and QQ plot from `results.tsv` with the QTNs from `qtn_recovery.tsv` marked.
