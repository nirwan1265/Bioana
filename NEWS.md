# bioana 0.2.0

- All of Simitall's analysis functions moved here: GWAS (`analyze_gwas()`,
  `analyze_binary_gwas()`, `benchmark_gwas()`, `plot_gwas_results()`), genomic
  prediction and selection, peak calling and differential peaks, single-cell
  clustering and condition tests, eQTL benchmarks, RNA-seq results summary,
  GO annotation with the gene-coordinate catalogues, and assembly evaluation.
- `bioana_analyze()` returns the truth tables (`recovery`) and fits, and takes
  `quiet`.

# bioana 0.1.0

- `bioana_analyze()`: every data set in a Biodes design folder, by data type
  from `data_manifest.json`, with the design's model; scored against the truth
  when simulated.
