# Showcase: Bioana analyses three simulated experiments and scores itself against the truth

Each folder follows one question through the suite. [Biodes](https://github.com/nirwan1265/BioDesignGPT/tree/main/examples/showcase)
turns the question into a design, [Simitall](https://github.com/nirwan1265/simitall/tree/main/examples/showcase)
simulates data for that layout and records what it planted, and Bioana analyses the folder and compares its
results with what was planted. Biopaint (planned) will draw the figures from the tables in each `results/` folder.
The folder names are the same in every repository.

| Example | Bioana's model (default engine) | Headline | Against the truth |
|---|---|---|---|
| [01 GWAS, sorghum panel](01_gwas_sorghum) | line BLUEs from the alpha lattice, then a kinship + 3 PC mixed-model scan (rrBLUP) | 6 markers at FDR 0.05 | 6 of 10 QTNs found, 0 false positives; line means vs true genetic values r = 0.86 |
| [02 RNA-seq, genotype x drought](02_rnaseq_drought) | `~ treatment * treatment2 + block`, Wald test of the interaction (DESeq2) | 273 interaction genes at FDR 0.05 | 215 of 500 goal genes found; 51 calls (19%) have no planted effect, against 1 of 91 calls with edgeR QL |
| [03 Seed priming, time to germination](03_germination_survival) | cause-specific Cox with tray and cell frailties (coxme), competing risks (cmprsk) | weighted log-rank FH(0,1) p = 3.6e-5 | HR after day 6: line 0.63 [0.45, 0.87], priming 0.67 [0.49, 0.91], true 0.6; 24 of 24 group measures in their CI |

Each folder's README gives the question, the model, the headline result, the recovery against the truth, and
what is in `results/`.

Rebuild all three from the Bioana repository root (seed 2026 data from Simitall):

```sh
Rscript examples/showcase/make_showcase.R                      # reads ../simitall/examples/showcase, else GitHub
SIMITALL_SHOWCASE=/path/to/simitall/examples/showcase Rscript examples/showcase/make_showcase.R
```

The script copies each Simitall folder to a temporary directory, runs `bioana_analyze()` with the default
engines (edgeR is also run for 02 as the second engine), and writes `results/`: the main tables, `report.txt`
(the printed report), `bioana_summary.json` and `truth_recovery.json`.
