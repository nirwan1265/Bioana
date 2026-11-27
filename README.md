# Bioana: analyses biological data and checks its own answers

**Type a research question in plain words; get a designed experiment, data, and a statistically sound answer
that has been checked against the known truth.** Bioana is the analysis step of a three-part suite for biology
experiments:

| Step | Tool | What it does |
|---|---|---|
| 1. Design | [**Biodes**](https://github.com/nirwan1265/BioDesignGPT) | turns the question into a randomized experimental design: layout, sample sheet, model, sample size |
| 2. Simulate | [**Simitall**](https://github.com/nirwan1265/simitall) | creates realistic data for that exact layout, with effects whose true size is known |
| 3. Analyse | **Bioana** (this repository) | picks the standard R package for the data type (rrBLUP for GWAS, DESeq2 for RNA-seq, survival for time-to-event data, and many more), fits the design's model, and grades the answer against the truth |

Because the truth is known, every analysis is graded: did it find the real effects, did it report effects
that are not there, and how close are the estimates? On real data, Bioana runs the same analysis without the
grading.

## Example: does seed priming delay germination?

A researcher asks:

> *"Seeds of 2 sorghum lines primed or not in 12 trays in a growth chamber, 10 seeds per treatment per tray,
> days to germination checked every 2 days, followed for 21 days; seed rot can occur before germination; the
> priming effect is delayed; hazard ratio 0.6; germination by day 7 and 14."*

**Biodes designs it:** 12 trays, each with the 4 line x priming combinations once, 480 seeds. It also chooses a
test suited to an effect that starts late and works out how many germinations are needed (121).

![Biodes tray layout](examples/showcase/03_germination_survival/figures/biodes_design.png)

**Simitall simulates it:** a germination day for each seed. After day 6, priming and line 2 each lower the
germination rate to 60% (the truth). Some seeds rot first.

![Simulated germination](examples/showcase/03_germination_survival/figures/simitall_data.png)

**Bioana analyses it in one line of R:**

```r
library(bioana)
res <- bioana_analyze("03_germination_survival")
```

```
      Kaplan-Meier [survival 3.5.8] median: genotype_1 / not primed 14 [10, 16] (75/120 events); genotype_2 / not primed 14 [10, 18] (67/120 events); genotype_1 / primed 16 [14, 18] (71/120 events); genotype_2 / primed NA [18, NA] (47/120 events)
      weighted log-rank chi-squared (fh_0_1, prespecified)(3) = 23.2, p = 3.6e-05
      Cox HR treatment2primed 0.83 [0.57, 1.21], p = 0.331
      ...
      goal truth (survival, terms treatment, treatment2): detected (p = 3.6e-05), estimate 0.63 vs true 0.6, in CI; piecewise HR after t = 5.99: treatment 0.63 [0.454, 0.873] vs 0.6, treatment2 0.667 [0.488, 0.911] vs 0.6; ...; per group 23/24 measures in CI
```

**The answer, graded:** yes. After day 6, primed seeds germinate at 0.67 times the rate of unprimed seeds
(95% CI 0.49 to 0.91), against a truth of 0.6. A standard Cox model, which averages early and late, would have
reported 0.83 with p = 0.33 and missed it. 23 of 24 group summaries (medians, share germinated by day 7 and 14,
mean time to germination) match the truth within their confidence intervals, about what 95% intervals should
give.

## Two more examples

| Question | Bioana's answer | Graded against the truth |
|---|---|---|
| Which DNA markers control plant height in 300 sorghum lines? (GWAS, rrBLUP) | 6 markers | 6 of the 10 real genes found, no false hits |
| Which genes respond to drought differently in a mutant? (RNA-seq, DESeq2 and edgeR) | 273 genes (DESeq2), 91 (edgeR) | DESeq2: 215 of 500 real genes, but 1 in 5 calls false; edgeR: 89 real, 1 false |

Each example (question, design figure, simulated data, R code, R output and graded answer) is in
**[examples/showcase](examples/showcase)**.

---

# Technical reference

Bioana takes a **Biodes design folder filled with data** (simulated by Simitall or measured by you) and runs,
for every data set in it, the analysis that data type calls for, with the model the design calls for.

```
question --Biodes--> design folder --Simitall or the lab--> data + data_manifest.json --Bioana--> results (+ truth scores)
```

Bioana never reads the question: Biodes already turned it into a spec. It never simulates: it treats simulated
and measured data the same way.

## Quick start

Simulated data (Simitall writes `data_manifest.json` for you):

```r
simitall::simitall_from_biodes("~/biodes_runs/my_design")
bioana::bioana_analyze("~/biodes_runs/my_design")
#> bioana: design 'rcbd' (y ~ treatment + block)
#>   [outcome] continuous outcome, lm ; truth: 1/1 coefficients in their CI
#>   [rnaseq] 500 features tested, treatmenttreatment_2 10 significant; truth: treatmenttreatment_2 found 9/50 (FDR 0.1)
#>   [methylation] 4809 features tested, ...; truth: found 80/119 (FDR 0.059)
```

Your own data: put the files in the design folder (matrix columns = the
`sample_id` values of `sample_metadata_template.tsv`), then say what they are:

```r
library(bioana)
bioana_write_manifest("~/biodes_runs/my_design",
  list(id = "outcome", data_type = "outcome", files = list(table = "sample_metadata_template.tsv"),
       column = "primary_outcome", outcome_type = "continuous"),
  list(id = "rnaseq", data_type = "feature_counts", assay = "rnaseq", files = list(counts = "counts.tsv")))
bioana_analyze("~/biodes_runs/my_design")
```

## Install

```r
# install.packages("remotes")
remotes::install_github("nirwan1265/bioana")
```

With every engine (from a clone):

```bash
git clone https://github.com/nirwan1265/bioana && cd bioana
Rscript install_bioana.R            # bioana + engines, readers, genetics (the "analysis" profile)
Rscript install_bioana.R full       # + plots and test tools;  "suite" adds Simitall
```

Or from R: `bioana_install_dependencies("analysis")`, and `bioana_check_dependencies()` to see
what is missing. Every package Bioana uses, with its source (CRAN, Bioconductor, GitHub) and what it is
for, is in `inst/rules/requirements.tsv`.

## Your own data

For real data the flow is **Biodes -> Bioana** (Simitall is not needed).
Biodes describes the experiment you ran and writes the design folder (the
`biodes_from_data()` contract: [docs/contracts/biodes_from_data.md](docs/contracts/biodes_from_data.md)).
Then:

```r
bioana_manifest_from_files(dir, c("featurecounts.txt", "genotypes.vcf"), assay = "rnaseq")
bioana_analyze(dir)
```

A worked real example, the DESeq2 vignette's airway data (built the way
Biodes' describe-my-data mode should write it):

```r
d <- bioana_example("airway")      # 4 cell lines x dexamethasone, Himes et al. 2014
bioana_analyze(d, secondary = TRUE) # DESeq2 (= the vignette, checked in test-real-airway.R), edgeR as secondary
```

`bioana_manifest_from_files()` detects each file's format, converts it to
Bioana's files under `<dir>/data/<id>/`, matches sample names to
`sample_id` (BAM paths and tool suffixes are stripped), and writes
`data_manifest.json` with the source files and tool recorded.

| format | detected by | becomes |
|---|---|---|
| BAM (one per library) | file name | reads counted with **Rsubread featureCounts** in genes (`annotation` = GTF), peaks (SAF or peak files) or exon windows (IP/input assays); `libraries` says which file is which unit and library (`chip`, `control`, `ip`, `input`); single-end ChIP reads extended to 125 bp fragments as DiffBind does |
| featureCounts | `# Program:featureCounts` or its columns | counts: `feature_counts`, or peaks for ChIP/ATAC/CUT&Tag assays |
| Salmon `quant.sf`, kallisto `abundance.tsv` (one folder per sample) | columns | gene counts via **tximport** (`tx2gene`, `lengthScaledTPM`) |
| STAR `ReadsPerGene.out.tab` | `N_unmapped` rows | counts; strandedness detected or given |
| Bismark `.cov` / cytosine report (one file per sample) | line pattern | `bisulfite_counts` (context kept from reports) |
| VCF | `##fileformat=VCF` | `genotypes` |
| HapMap (`.hmp.txt`, TASSEL / GAPIT) | `rs#  alleles  chrom  pos` header | `genotypes` (converted to VCF) |
| DSS / BSmooth text (`chr pos N X`, one per sample) | columns | `bisulfite_counts` |
| plain feature x sample table | file type | integers: counts; 0-1: `methylation_beta`; `<id>_IP` / `<id>_input` columns: `ip_input_counts` |
| MACS narrowPeak / broadPeak | file name | peak calls have no counts: `bioana_consensus_peaks()` writes a SAF to count with `featureCounts -F SAF` |

Formats live in `inst/rules/file_formats.tsv` (pattern, header, columns,
reader) and `file_formats.json` (sample-name rules, assay -> data type). A new
format is a row and a reader. Normalised values (TPM, CPM) are refused: the
engines need raw counts.

### Checked on real data

`tests/testthat/test-real-data.R` runs published data sets through Bioana. The
data live in [bioana-testdata](https://github.com/nirwan1265/bioana-testdata);
set `BIOANA_TESTDATA` or clone it next to `bioana/`.

| data | check | result |
|---|---|---|
| airway counts (DESeq2 vignette) | Bioana = vignette code | log2FC to 5e-9; known dexamethasone genes up |
| airway BAMs + GTF | featureCounts = the package's counts | median ratio 1, r = 0.998 |
| DiffBind tamoxifen ChIP, chr18 1-20 Mb | Bioana's DB calls vs DiffBind's own (DESeq2) | 91-96% of Bioana's calls are DiffBind's, same direction at every one; DiffBind calls more sites (it subtracts input reads, which raises its folds about 1.5x) |
| maize 282 panel, 2 locations (TASSEL) | line means; GWAS = rrBLUP called directly | r = 0.97 with TASSEL's line values; -log10 p to 1e-11 |
| DSS example counts | Bioana at site level = DSS called directly | identical statistics, same 247 DML |
| bsseq Bismark files | formats read | coverage and cytosine report |
| exomePeak2 MeRIP BAMs | windows, IP/input pairing, model runs | runs; too small for a differential check |

## Analysis goals: what is tested

Biodes writes the goal of the study into `study_spec.json` (`analysis.goals`,
`analysis.primary_goal`, `analysis.goal_terms`; contract:
[docs/contracts/analysis_goals.md](docs/contracts/analysis_goals.md); the vocabulary lives in Biodes'
`inst/registry/analysis_goals.json`). `engine_selection.tsv` is keyed by goal (`goals` column) and each
engine says which tests it can run (`engines.tsv` `tests`: `wald`, `lrt`).

| goal | test | how |
|---|---|---|
| `treatment_effect` (default) | Wald, each treatment coefficient | as before |
| `interaction` | Wald if one coefficient, joint test if several | DESeq2 LRT (reduced: without the interaction), edgeR QL F, limma/voom/dream moderated F, outcome `anova(reduced, full)` |
| `time_course` | joint test: the model without every term containing time | same engines |
| `trend` | Wald on the slope of a numeric term | the term is fitted as numeric |
| `association` | per marker | GWAS engines |

Engines without a joint test (DSS, the IP x input model) are skipped for joint-test goals. With no goal
in `study_spec`, Bioana tests treatment effects and says so when the model has time or interaction
terms. `bioana_goals()` lists what is available per data type; `params = list(<id> = list(goal = ...))`
picks one. Each analysis records its goal, test and the coefficients tested in `summary$...$goal`;
outcome joint tests go to `goal_test.tsv`.

## Safeguards: what Bioana refuses to do

Before fitting a data set Bioana checks the design (`bioana_check_design(dir)` runs the same checks on
the whole sample sheet; thresholds in `inst/rules/design_checks.json`). It **stops** when:

- a treatment coefficient is confounded with another term (e.g. every control in one batch);
- there are no residual degrees of freedom;
- a subject (plant, animal, donor, ...) is measured more than once but `analysis.model` has no random term
  for it: the samples would be treated as independent when they are not;
- several libraries come from one unit under the same treatment (technical replicates);
- a random term has one observation per level, or one level;
- a count matrix holds normalised values (CPM, TPM), or the outcome values do not match its type
  (binary not 0/1, negative counts, survival events not 0/1).

It **warns** (in the report and `summary$diagnostics`) about few residual df, treatment levels with one
unit, random terms with few levels, and every warning a package raised while fitting (singular fits,
DESeq2 genes that did not converge, glmmTMB Hessian problems).

It **never swaps in a method that does not fit the design**. Every engine row in
`engine_selection.tsv` says how it fits the design (`design_fit`): `full`, `approx` (e.g. limma
duplicateCorrelation for one random intercept instead of dream's mixed model; marked "approximates the
design" in the report) or `drops_random` (only with `allow_fallback = TRUE`; marked "RANDOM TERMS
DROPPED"). An engine that cannot hold the design's random terms (`max_random` in `engines.tsv`) is
skipped. If an engine fails while fitting, the next one runs only if it fits the design fully;
otherwise the analysis stops and says which engine you could choose yourself.

## Which packages it runs (gold standards, not home-made statistics)

Every analysis calls a published package, chosen by two rule tables:
`inst/rules/engines.tsv` (package, function, model, test statistic, multiple
testing, assumptions, reference) and `inst/rules/engine_selection.tsv` (engine
order per data type and per design: `fixed` = no random terms in
`analysis.model`, `random` = repeated or nested units).

| data | fixed design | random terms (same individuals over time, split plots, ...) |
|---|---|---|
| gene / miRNA / taxon counts | DESeq2 → edgeR → limma-voom | dream → limma-voom + duplicateCorrelation |
| ChIP peaks, ATAC windows | edgeR (as csaw) → DESeq2 | dream → limma-voom |
| bisulfite counts | DSS (beta-binomial) → limma on M-values | limma on M-values |
| methylation arrays | limma on M-values | limma + duplicateCorrelation |
| m6A / MeDIP / RIP | edgeR IP x treatment interaction → limma on log2(IP/input) | limma |
| primary outcome | lm, glm, glm.nb, clm, coxph → glmmTMB | lmerTest, glmer, glmer.nb, clmm, coxme (survival; frailty Cox as approximation) → glmmTMB |
| genotypes (GWAS on line means) | rrBLUP → GAPIT MLM | (kinship random effect) |

The first installed engine is the **primary**; the next is the **secondary**.
Only published packages run by default (`allow_fallback = FALSE`): if none
is installed for a data set, Bioana stops and names the packages to install,
so every result can be cited.
`bioana_analyze(dir, secondary = TRUE)` runs both and reports how well they
agree (in an interactive session Bioana asks). Bioana's own code runs only
with `allow_fallback = TRUE` (or when named as the engine), and is marked
`FALLBACK` in the report. Every data set's
`engine` block in `bioana_summary.json` records package, version, exact call,
formula, model, input, assumptions, test statistic, multiple testing,
experimental unit and whether a fallback was used.

```r
bioana_engines()                                  # the registry, with what is installed
bioana_engine_plan("feature_counts", "random")    # what would run for a repeated-measures RNA-seq design
bioana_analyze(dir, params = list(feature_counts = list(engine = "edger")))   # choose an engine
```

Tutorials (each package's own example code, then the same through Bioana) and
the gold-standard checks (Bioana's adapter must give the same numbers as the
package called directly) are in [docs/engines/](docs/engines/README.md).

## What it analyses

`bioana_analyses()` prints `inst/rules/analyses.tsv`:

| data_type | files | assays | analysis |
|---|---|---|---|
| `outcome` | table | primary outcome | the design's `analysis.model`: lm / lmer (continuous), glm / glmer binomial (binary), negative binomial (count), cumulative logit (ordinal), Cox (survival, needs an `event` column) |
| `feature_counts` | counts | RNA-seq, 3' tag, small RNA, microbiome, single-cell pseudobulk | log2-CPM, per-feature linear model with the design's fixed terms, empirical-Bayes moderated t (as limma) |
| `peak_counts_long` | counts, libraries | ChIP-seq, DAP-seq, CUT&Tag | IP reads per peak and library, mapped to units; then as above |
| `read_positions` | reads, libraries | ATAC-seq, DNase-seq, FAIRE-seq | reads counted in windows (no peak calling); then as above |
| `bisulfite_counts` | sites, methylated, coverage | WGBS, RRBS, EM-seq | sites pooled in windows per context, M-values, moderated t |
| `methylation_beta` | beta | 450k / EPIC arrays | M-values, moderated t per probe |
| `ip_input_counts` | ip, input | m6A MeRIP-seq, MeDIP, RIP | log2(IP/input) per window, so expression changes cancel |
| `genotypes` | vcf | GWAS, GBS, WGS, QTL mapping | line means (BLUEs) of the primary outcome per line from the design's `analysis.model`, then a kinship mixed-model GWAS (rrBLUP, or native EMMAX without it); scored against `qtn_truth.tsv` by exact marker or LD |
| `single_cell_counts` | counts, cells | scRNA-seq (cell level) | planned in `bioana_analyze()`; `cluster_scrnaseq()` works directly |

Settings per data type: `inst/rules/method_parameters.tsv`, overridden with
`bioana_analyze(dir, params = list(bisulfite_counts = list(window_bp = 500)))`.

## Other analyses (moved from Simitall)

These work on files directly (Simitall's simulator outputs or the same formats
from real data) and are not yet wired into `bioana_analyze()`:

| area | functions |
|---|---|
| GWAS | `analyze_gwas()` (mixed model with kinship and PCs, rrBLUP), `analyze_binary_gwas()`, `benchmark_gwas()` (recovery of causal variants), `plot_gwas_results()` (Manhattan, QQ) |
| Genomic prediction and selection | `calculate_genomic_kinship()`, `fit_genomic_model()` (GBLUP, RR-BLUP, BGLR Bayesian models, random forest), `predict_genomic_values()`, `cross_validate_genomic_prediction()`, `select_parents()`, `design_parent_crosses()`, `run_genomic_selection()` |
| eQTL and RNA-seq | `benchmark_eqtl()`, `benchmark_celltype_eqtls()` (donor pseudobulk per cell type), `summarize_rnaseq_results()`, `plot_rnaseq_eqtl_results()` |
| Single cell | `cluster_scrnaseq()`, `score_scrnaseq_clusters()`, `test_scrnaseq_conditions()`, `score_scrnaseq_conditions()`, `plot_scrnaseq_results()` |
| ChIP / ATAC peaks | `call_peaks_windowed()`, `score_peak_calls()`, `test_differential_peaks()`, `score_differential_peaks()`, `plot_chipseq_results()` |
| Genes and GO | `read_gene_coordinates()` (bundled catalogues: Arabidopsis, barley, rice, wheat, maize, fly), `read_go_gaf()`, `join_gene_go()`, `summarize_go_terms()` |
| Genome assembly | `run_unicycler_grid()`, `run_quast_grid()`, `summarize_quast()` |
| Figures | `export_bioana_plot_panels()` saves each panel of a results figure |

```r
g <- bioana::analyze_gwas("cohort.vcf", phenotypes, "gwas_out", trait = "yield", n_pcs = 3)
bioana::benchmark_gwas(g$paths$results, "cohort.causal.tsv")
fit <- bioana::fit_genomic_model("cohort.geno.tsv", phenotypes, trait = "yield", model = "gblup")
```

## Output

`<folder>/bioana/`:

- `<id>/results.tsv` — outcome: one row per coefficient (`term, estimate, se, df, statistic, p, role, padj`);
  features: one row per feature and treatment contrast (`feature_id, contrast, estimate, se, statistic, p, padj`;
  `statistic` is the engine's test statistic: Wald z for DESeq2, QL F for edgeR, moderated t for limma/dream).
- `<id>/results_<engine>.tsv` — the secondary engine's table, when `secondary = TRUE`.
- `<id>/recovery.tsv` — outcome only, simulated data: true vs estimated coefficients.
- genotypes: `<id>/results.tsv` one row per tested marker (`marker_id, seqname, pos, estimate, se, p, padj, maf, significant`;
  `estimate` = effect per alternative allele, on the outcome scale); `<id>/line_means.tsv` (`line, mean, se, n_units`);
  `<id>/qtn_recovery.tsv`, simulated data: one row per QTN (`marker_id, seqname, pos, type, effect, maf, tested, recovered,
  estimate, se, p, padj`); `<id>/gwas.gwas.tsv` is `analyze_gwas()`'s own table, for `plot_gwas_results()`.
- `bioana_summary.json` — per data set: method, model fitted, counts tested and significant,
  warnings, the truth scores (`tpr`, `fdr`, `sign_agreement`, `estimate_cor`, CI coverage), the `engine`
  provenance block, `secondary_engine`, and with `secondary = TRUE` the secondary's provenance,
  `concordance` (estimate r, sign agreement, significant in both / one only, Jaccard) and truth scores.

## The link to Simitall: `data_manifest.json`

```json
{
  "manifest_version": "0.1",
  "written_by": "simitall 0.2.0",
  "data": [
    {"id": "outcome", "data_type": "outcome", "files": {"table": "sample_metadata_template.tsv"},
     "column": "primary_outcome", "outcome_type": "continuous",
     "truth": {"kind": "outcome", "file": "sim_truth.json"}},
    {"id": "rnaseq", "data_type": "feature_counts", "assay": "rnaseq", "files": {"counts": "counts.tsv"},
     "truth": {"kind": "features", "file": "truth/rnaseq.tsv", "term": "treatment", "large": 1}},
    {"id": "methylation", "data_type": "bisulfite_counts", "assay": "methylation",
     "files": {"sites": "methylation/methylation.sites.tsv.gz", "methylated": "...", "coverage": "..."},
     "truth": {"kind": "regions", "file": "truth/methylation.tsv", "term": "treatment"}}
  ]
}
```

Each simulated data type in Simitall maps to one `data_type` here; a new
simulator in Simitall adds a manifest entry, and a new analysis here adds a
row to `analyses.tsv`. Truth files: `features` = `feature_id, condition,
effect`; `regions` = `region_id, chr, start, end, condition, effect`; the
tested contrast is `term` + `condition` (e.g. `treatmenttreatment_2`).

## Future analyses to check

Each needs an analysis here, a row in `inst/rules/analyses.tsv`, a test, and
(for simulated data) the matching Simitall manifest entry. "Questions" = how
many of Biodes' 4,695 paper questions need it.

**Wire existing functions into `bioana_analyze()`**

- [x] `genotypes`: GWAS / QTL mapping with `analyze_gwas()` on the design's
      line means, scored with `benchmark_gwas()` against `genetic_truth/qtn_truth.tsv`.
      Next for it: GWAS per environment under GxE, ordinal and survival outcomes,
      QTL mapping with cross-aware genotype probabilities.
- [x] Real data readers: featureCounts, Salmon/kallisto (tximport), STAR, Bismark, VCF, plain
      tables, MACS peak calls (consensus SAF); `bioana_manifest_from_files()`. Still to do: BAM
      input (counting, exomePeak2) and a run on Simitall's real templates.
- [ ] eQTL designs: `benchmark_eqtl()` on `counts.tsv` + `genotypes.vcf`.
- [ ] Single cell at cell level: `cluster_scrnaseq()`, cell-type proportions,
      `test_scrnaseq_conditions()`; cell-type eQTL.
- [ ] Peaks: `call_peaks_windowed()` + `test_differential_peaks()` as the
      peak-calling alternative to fixed windows for ChIP / ATAC.
- [ ] Genomic prediction for breeding designs (`fit_genomic_model()`,
      `cross_validate_genomic_prediction()`), scored against the true genetic values.

**Analyses for data types Simitall will simulate**

- [ ] Metabolomics / lipidomics (54 questions): log intensities, missing
      values, moderated t; pathway / class summaries.
- [ ] Genetics with no lines (36): BSA / QTL-seq (allele-frequency difference
      between bulks, G' statistic), pool-seq (allele frequencies over
      generations), population scans (Fst, selection).
- [ ] Ribo-seq (42): translational efficiency, the RNA x footprint interaction.
- [ ] scATAC-seq (41): peaks x cells, cluster-level differential accessibility.
- [ ] Hi-C (38): differential contacts, loops and TAD boundaries.
- [ ] MeDIP (31) and RIP/CLIP (14): already `ip_input_counts`; check on their
      simulated data.
- [ ] Degradome (22): cleavage sites and miRNA targets.
- [ ] Spatial (19): spatially variable genes, domains.
- [ ] MNase-seq (16): nucleosome positions and occupancy changes.
- [ ] CRISPR screens: guide-level counts to gene scores (MAGeCK-like).
- [ ] Amplicon genotyping: editing efficiency, allele frequencies per sample.
- [ ] TCR/BCR repertoires: clonotype abundance and diversity.
- [ ] Proteomics: log intensities with missing values (limma-like).
- [ ] STARR-seq / MPRA: RNA over DNA barcode counts per element.
- [ ] Compositional outcomes; survival with random effects (frailty / coxme).

**Methods and checks**

- [x] Package adapters with the native methods as fallback: DESeq2, edgeR,
      limma-voom, limma, dream, DSS, glmmTMB, rrBLUP, GAPIT (engine registry,
      selection rules, provenance, tutorials and gold checks).
- [ ] Next adapters: qtl2, muscat, csaw / DiffBind, ANCOM-BC2 / ALDEx2,
      exomePeak2 (BAM), methylKit, MSstats / proDA, clusterProfiler / fgsea,
      mixOmics; then Python engines (Scanpy/scVI-tools, PyDESeq2, statsmodels,
      scikit-learn, COBRApy).
- [x] Random terms in feature tests: dream (any lme4 terms) or limma
      duplicateCorrelation (one grouping term).
- [ ] Contrasts beyond "each level vs the reference": pairwise and interaction
      contrasts for multi-level factors.
- [ ] Power for every data type: Simitall's power sweep uses Bioana's scores
      for the outcome and counts; add methylation (DMRs), m6A, peaks, GWAS.
- [ ] Real-data readers: featureCounts / Salmon tables, Bismark coverage,
      MACS2 peaks, VCF; QC and normalisation reports.
- [ ] Check low power where the simulation says so: small RNA and microbiome
      templates vary a lot between real samples (0 found with 5 per group);
      m6A at 5 per group (2/40).

## Status

Version 0.5: `bioana_analyze()` for the data types in the first table, and
all of Simitall's former analysis functions (second table), tested on
synthetic folders and on Simitall output (`devtools::test()`; the moved tests
run when Simitall is installed). See `AGENTS.md` for how the code is organised.
