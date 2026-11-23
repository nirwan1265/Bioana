# Native fallback (`bioana_native`)

Bioana's own code, **not citable**. It never runs by default. It runs only
with `bioana_analyze(dir, allow_fallback = TRUE)` when no gold-standard
engine for the data type is installed, or when you name it:
`params = list(<data_type> = list(engine = "bioana_native"))`. When it runs:
the report line says `FALLBACK`, `engine$fallback_used` is `true`, and
`engine$not_used` lists the engines skipped and why.

**What it fits.** A linear model per feature with the design's **fixed** terms
(random terms are dropped, and a warning says so). It then applies
empirical-Bayes variance moderation: a re-implementation of limma's
`squeezeVar`/`fitFDist` (a scaled-F prior fitted by the moments of log s2),
then moderated t and Benjamini-Hochberg. The inputs are log2-CPM (counts),
M-values (methylation) or log2(IP/input) (m6A).

**Gold check.** In "limma on values = lmFit/eBayes; the native fallback
reproduces limma", the native estimates equal limma's exactly. The moderated
t statistics correlate > 0.999 and differ by less than 5% of the largest
statistic.

Install the engines in [README.md](README.md) so the fallback is not needed.
