# Primary outcome engines

The primary outcome is fitted with the formula Biodes wrote
(`analysis.model`), with the family the outcome type calls for:

| type | fixed design (engine) | random terms (engine) | secondary |
|---|---|---|---|
| continuous | `stats::lm` (`stats_lm`) | `lmerTest::lmer`, Satterthwaite df (`lmer`) | glmmTMB (gaussian) |
| binary | `stats::glm(binomial)` (`glm_binomial`) | `lme4::glmer(binomial)` (`glmer_binomial`) | glmmTMB (binomial) |
| count | `MASS::glm.nb` (`glm_nb`) | `lme4::glmer.nb` (`glmer_nb`) | glmmTMB (`nbinom2`) |
| ordinal | `ordinal::clm` (`clm`) | `ordinal::clmm` | none |
| survival | `survival::coxph` (`coxph`) | random terms dropped, with a warning | none |

Treatment coefficients are adjusted with Benjamini-Hochberg. Subsamples of one
experimental unit are averaged (`subsamples = "aggregate"`, default) or nested
as a random effect (`"nest"`, continuous only). They are never treated as
replicates.

## Tutorial: the same model directly

```r
library(lmerTest)
fit <- lmer(yield ~ treatment + (1 | block), data = d)    # random blocks
summary(fit); anova(fit)
library(glmmTMB)
fit2 <- glmmTMB(count ~ treatment + (1 | block), family = nbinom2, data = d)
summary(fit2)$coefficients$cond
```

## Through Bioana

```r
r <- bioana_analyze(dir, secondary = TRUE)
r$summary$analyses$outcome$engine                    # lm / lmer / ... with version and call
r$summary$analyses$outcome$secondary$concordance     # agreement with glmmTMB
```

**Gold check.** "glmmTMB outcome engine matches glmmTMB called directly, and
agrees with lm": the glmmTMB estimates equal a direct `glmmTMB()` call, and
correlate > 0.999 with lm's.

**References.** Bates D, Mächler M, Bolker B, Walker S (2015). *J Stat Softw*
67(1). Kuznetsova A, Brockhoff PB, Christensen RHB (2017). *J Stat Softw*
82(13). Brooks ME et al. (2017). glmmTMB. *The R Journal* 9:378. Christensen
RHB, ordinal (R package). Therneau TM, Grambsch PM (2000). *Modeling Survival
Data*.
