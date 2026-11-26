# Survival layer (R/survival_tests.R) on the survival package's data sets: equality with survival, cmprsk,
# survMisc, survminer and KMsurv, and the numbers published in the survival tutorials (datanovia).

strata <- survival::strata   # the strata() special in survdiff formulas below
surv_run <- function(d, args = NULL) { utils::capture.output(r <- bioana_analyze(d, secondary = FALSE, params = if (!is.null(args)) list(outcome = args)))
  r$summary$analyses$outcome$stat_test }
trow <- function(st, pat) Filter(function(r) grepl(pat, r$statistic_name), st$table)
cluster <- survival::cluster   # the cluster() special in coxph formulas below
Surv <- survival::Surv
lung_df <- function() { e <- new.env(); utils::data("cancer", package = "survival", envir = e); x <- e$lung
  x$sex <- factor(c("male", "female")[x$sex], c("male", "female")); x }

test_that("lung: KM medians, log-rank, Cox HR, landmarks = survival and the tutorial; weighted and stratified log-rank on request", {
  st <- surv_run(bioana_example("lung"), list(survival = list(logrank_weight = "fh_0_1", strata = list("ph.ecog"))))
  expect_equal(st$test, "survival")
  km <- st$km; names(km) <- vapply(km, `[[`, "", "group")
  expect_equal(c(km$male$median, km$male$median_lcl, km$male$median_ucl), c(270, 212, 310))       # tutorial: 270 (212-310)
  expect_equal(c(km$female$median, km$female$median_lcl, km$female$median_ucl), c(426, 348, 550))  # 426 (348-550)
  x <- lung_df()
  lr <- survival::survdiff(survival::Surv(time, status) ~ sex, data = x)
  expect_equal(trow(st, "^log-rank")[[1]]$statistic, lr$chisq); expect_equal(round(lr$chisq, 1), 10.3)
  st0 <- surv_run(bioana_example("lung"))                                                            # no strata: the tutorial's Cox
  cx <- survival::coxph(survival::Surv(time, status) ~ sex, data = x)
  expect_equal(trow(st0, "Cox")[[1]]$estimate, unname(exp(stats::coef(cx))), tolerance = 1e-8)
  expect_equal(round(trow(st0, "Cox")[[1]]$estimate, 2), 0.59)                                       # tutorial: HR 0.59
  z <- survival::cox.zph(cx); expect_equal(st0$ph_check[[1]]$p, unname(z$table[1, "p"]), tolerance = 1e-8)
  cs <- survival::coxph(survival::Surv(time, status) ~ sex + strata(ph.ecog), data = x)               # strata asked: Cox stratified too
  expect_equal(trow(st, "Cox")[[1]]$estimate, unname(exp(stats::coef(cs))), tolerance = 1e-8)
  ls <- survival::survdiff(survival::Surv(time, status) ~ sex + strata(ph.ecog), data = x)
  expect_equal(trow(st, "stratified")[[1]]$statistic, ls$chisq); expect_equal(round(ls$chisq, 1), 10.8)  # tutorial: 10.8
  expect_equal(round(trow(st, "weighted")[[1]]$statistic, 2), 3.45)                                       # tutorial: FH(0,1) 3.45
  st2 <- surv_run(bioana_example("lung"), list(survival = list(logrank_weight = "peto")))
  expect_equal(trow(st2, "weighted")[[1]]$statistic, survival::survdiff(survival::Surv(time, status) ~ sex, data = x, rho = 1)$chisq)  # 12.7
  d <- bioana_example("lung"); utils::capture.output(bioana_analyze(d, secondary = FALSE))
  lm_ <- utils::read.delim(file.path(d, "bioana", "outcome", "landmarks.tsv"))
  ref <- summary(survival::survfit(survival::Surv(time, status) ~ sex, data = x), times = c(180, 360))
  expect_equal(lm_$surv, ref$surv, tolerance = 1e-10)
  kc <- utils::read.delim(file.path(d, "bioana", "outcome", "km_curve.tsv"))
  expect_equal(nrow(kc), length(summary(survival::survfit(survival::Surv(time, status) ~ sex, data = x), censored = TRUE)$time))
})

test_that("RMST (colon, Obs vs Lev+5FU, tau 1825) = the tutorial's survRM2 numbers; AFT time ratios = survreg", {
  d <- bioana_example("colon")
  m <- utils::read.delim(file.path(d, "sample_metadata_template.tsv"), check.names = FALSE); m <- m[m$rx %in% c("Obs", "Lev+5FU"), ]
  utils::write.table(m, file.path(d, "sample_metadata_template.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
  s <- jsonlite::read_json(file.path(d, "study_spec.json")); s$factors[[1]]$levels <- list("Obs", "Lev+5FU"); s$factors[[1]]$n_levels <- 2L
  s$analysis$posthoc <- NULL; s$analysis$survival <- list(rmst_tau = 1825, aft = TRUE)
  jsonlite::write_json(s, file.path(d, "study_spec.json"), auto_unbox = TRUE, pretty = TRUE)
  st <- surv_run(d)
  r <- st$rmst
  expect_equal(r$tau_source, "analysis.survival.rmst_tau")
  expect_equal(round(r$difference, 1), 111.3); expect_equal(round(r$ci, 1), c(19.3, 203.4)); expect_equal(round(r$p, 3), 0.018)   # tutorial
  expect_equal(round(r$ratio, 3), 1.083); expect_equal(round(r$ratio_ci, 3), c(1.014, 1.157))
  x <- m; x$rx <- factor(x$rx, c("Obs", "Lev+5FU"))
  best <- st$aft$chosen
  sr <- survival::survreg(survival::Surv(time, event) ~ rx, data = x, dist = best)
  expect_equal(trow(st, "time ratio")[[1]]$estimate, unname(exp(stats::coef(sr)[2])), tolerance = 1e-6)
  expect_equal(unname(unlist(st$aft$aic)[best]), stats::AIC(sr), tolerance = 1e-6)
})

test_that("colon: pairwise log-rank (BH) = survminer::pairwise_survdiff and the tutorial", {
  skip_if_not_installed("survminer")
  d <- bioana_example("colon"); st <- surv_run(d)
  ph <- utils::read.delim(file.path(d, "bioana", "outcome", "posthoc.tsv"))
  e <- new.env(); utils::data("cancer", package = "survival", envir = e); x <- e$colon[e$colon$etype == 2, ]
  pw <- survminer::pairwise_survdiff(survival::Surv(time, status) ~ rx, data = x, p.adjust.method = "BH")$p.value
  get <- function(a, b) ph$p_adj[ph$contrast %in% c(sprintf("%s - %s", a, b), sprintf("%s - %s", b, a))]
  expect_equal(get("Lev", "Obs"), pw["Lev", "Obs"], tolerance = 1e-8)
  expect_equal(get("Lev+5FU", "Obs"), pw["Lev+5FU", "Obs"], tolerance = 1e-8)
  expect_equal(get("Lev+5FU", "Lev"), pw["Lev+5FU", "Lev"], tolerance = 1e-8)
  expect_equal(round(get("Lev+5FU", "Obs"), 3), 0.005); expect_equal(round(get("Lev", "Obs"), 2), 0.81)      # tutorial
})

test_that("mgus2 competing risks: Gray's test and Fine-Gray = cmprsk; Fine-Gray ~ survival::finegray; CIF 10% vs 1 - KM 21% at 20 years", {
  d <- bioana_example("mgus2"); utils::capture.output(r <- bioana_analyze(d, secondary = FALSE))
  a <- r$summary$analyses$outcome; st <- a$stat_test
  expect_equal(st$test, "competing_risks"); expect_match(paste(a$warnings, collapse = " "), "cause-specific")
  m <- utils::read.delim(file.path(d, "sample_metadata_template.tsv")); m$sex <- factor(m$sex, c("female", "male"))
  ci <- cmprsk::cuminc(m$time, m$event, m$sex, cencode = 0)
  expect_equal(trow(st, "Gray.*progression")[[1]]$p, unname(ci$Tests["1", "pv"]))
  f <- cmprsk::crr(m$time, m$event, cov1 = cbind(male = as.numeric(m$sex == "male")), failcode = 1, cencode = 0)
  expect_equal(trow(st, "Fine-Gray")[[1]]$estimate, unname(exp(f$coef)), tolerance = 1e-8)
  m$ev <- factor(m$event, 0:2, c("censor", "progression", "death"))
  fg <- survival::finegray(survival::Surv(time, ev) ~ ., data = m[, c("time", "ev", "sex")], etype = "progression")
  fc <- survival::coxph(survival::Surv(fgstart, fgstop, fgstatus) ~ sex, data = fg, weights = fgwt)
  expect_equal(trow(st, "Fine-Gray")[[1]]$estimate, unname(exp(stats::coef(fc))), tolerance = 0.01)
  all240 <- Filter(function(z) z$group == "all" && z$time == 240, st$cif_vs_naive)[[1]]
  expect_equal(round(all240$cif, 2), 0.10); expect_equal(round(all240$one_minus_km, 2), 0.21)                      # tutorial
  # the engine's Cox model is cause-specific: progression, deaths censored
  cs <- survival::coxph(survival::Surv(time, event == 1) ~ sex, data = m)
  expect_equal(r$results$outcome$estimate, unname(stats::coef(cs)), tolerance = 1e-8)
})

test_that("rats: litter frailty fit, log-rank stratified by litter = survdiff(strata); life table = KMsurv::lifetab", {
  skip_if_not_installed("KMsurv")
  d <- bioana_example("rats"); st <- surv_run(d, list(survival = list(life_table = 20)))
  e <- new.env(); utils::data("cancer", package = "survival", envir = e); x <- e$rats
  ls <- survival::survdiff(survival::Surv(time, status) ~ rx + strata(litter), data = x)
  expect_equal(trow(st, "stratified")[[1]]$statistic, ls$chisq)
  lt <- utils::read.delim(file.path(d, "bioana", "outcome", "life_table.tsv"), check.names = FALSE)
  g <- lt[lt$group == "control", ]
  s <- x[x$rx == 0, ]; br <- c(g$start, g$end[nrow(g)])
  bin <- cut(s$time, br, right = FALSE, labels = FALSE)
  ref <- KMsurv::lifetab(br, nrow(s), tabulate(bin[s$status == 0], length(br) - 1), tabulate(bin[s$status == 1], length(br) - 1))
  expect_equal(g$surv, ref$surv); expect_equal(g$nrisk, ref$nrisk)
})

test_that("CDISC ADTTE: event = 1 - CNSR, one PARAMCD, arm as treatment; same log-rank as the source data", {
  e <- new.env(); utils::data("cancer", package = "survival", envir = e); x <- e$colon[e$colon$etype == 2 & e$colon$rx %in% c("Obs", "Lev+5FU"), ]
  ad <- data.frame(STUDYID = "COLON", USUBJID = sprintf("C-%04d", x$id), PARAMCD = "OS", PARAM = "Overall Survival", AVAL = x$time, AVALU = "DAYS",
                   CNSR = 1L - x$status, TRT01P = as.character(x$rx))
  ad <- rbind(ad, transform(ad, PARAMCD = "PFS", PARAM = "Progression-free survival"))
  f <- tempfile(fileext = ".csv"); utils::write.csv(ad, f, row.names = FALSE)
  expect_error(bioana_adtte_design(f, file.path(tempfile(), "d")), "several endpoints")
  d <- bioana_adtte_design(f, file.path(tempfile(), "d"), paramcd = "OS", reference = "Obs")
  st <- surv_run(d)
  lr <- survival::survdiff(survival::Surv(time, status) ~ rx, data = x)
  expect_equal(trow(st, "^log-rank")[[1]]$statistic, lr$chisq)
  expect_equal(sum(utils::read.delim(file.path(d, "sample_metadata_template.tsv"))$event), sum(ad$CNSR[ad$PARAMCD == "OS"] == 0))
})

test_that("goal truth as a time ratio (Simitall lognormal AFT) is scored with survreg(lognormal), not the Cox log hazard ratio", {
  d <- file.path(tempfile("aft_"), "design"); dir.create(d, recursive = TRUE)
  set.seed(11); n <- 80; g <- rep(c("ctrl", "trt"), each = n / 2)
  t <- 20 * exp(0.5 * stats::rnorm(n) + log(1.6) * (g == "trt")); fu <- stats::quantile(t, 0.85)
  m <- data.frame(unit = sprintf("U%02d", 1:n), sample_id = sprintf("U%02d", 1:n), trt = g, time = pmin(t, fu), event = as.integer(t <= fu))
  spec <- list(handoff_version = "0.6.5", origin = "planned", organism = "test", system = "lab", factors = list(list(name = "trt", n_levels = 2L, levels = list("ctrl", "trt"))),
               design = list(chosen = "crd", reps = n / 2), primary_outcome = list(name = "days", type = "survival"),
               analysis = list(model = "y ~ treatment", terms = list(treatment = "trt", layout_columns = list("trt")),
                               goals = list("survival"), primary_goal = "survival", goal_terms = list(survival = list("treatment"))))
  jsonlite::write_json(spec, file.path(d, "study_spec.json"), auto_unbox = TRUE, pretty = TRUE)
  utils::write.csv(m[, c("unit", "trt")], file.path(d, "field_layout.csv"), row.names = FALSE)
  utils::write.table(m, file.path(d, "sample_metadata_template.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
  jsonlite::write_json(list(outcome_model = list(family = "lognormal_aft"), goal_truth = list(goal = "survival", model_terms = list("treatment"), status = "planted",
                       true_effect = 1.6, effect_scale = "time ratio (accelerated failure time), last level vs reference")), file.path(d, "sim_truth.json"), auto_unbox = TRUE)
  bioana_write_manifest(d, list(id = "outcome", data_type = "outcome", files = list(table = "sample_metadata_template.tsv"), column = "time", event_column = "event",
                                truth = list(kind = "outcome", file = "sim_truth.json", goal = list(goal = "survival", field = "goal_truth", terms = list("treatment"), test = "wald"))))
  utils::capture.output(r <- bioana_analyze(d, secondary = FALSE))
  gr <- r$summary$analyses$outcome$goal_recovery
  m$trt <- factor(m$trt, c("ctrl", "trt"))
  sr <- survival::survreg(survival::Surv(time, event) ~ trt, data = m, dist = "lognormal")
  expect_equal(gr$estimate, unname(exp(stats::coef(sr)[2])), tolerance = 1e-6)
  expect_true(gr$sign_agreement); expect_match(gr$scale, "survreg\\(dist = 'lognormal'\\)")
})

seed_folder <- function(trays = 24, seeds = 10, hr = 0.6, theta = 0.5, power = NULL, truth = TRUE, seed = 5) {
  d <- file.path(tempfile("seeds_"), "design"); dir.create(d, recursive = TRUE)
  set.seed(seed)
  tr <- data.frame(unit = sprintf("T%02d", seq_len(trays)), lot = rep(c("untreated", "primed"), length.out = trays), stringsAsFactors = FALSE)
  fr <- stats::rgamma(trays, 1 / theta, 1 / theta)                                     # tray frailty, mean 1
  m <- tr[rep(seq_len(trays), each = seeds), ]; m$unit_id <- m$unit; m$sample_id <- sprintf("%s_s%02d", m$unit, seq_len(seeds))
  h <- 0.05 * fr[match(m$unit, tr$unit)] * ifelse(m$lot == "primed", hr, 1)
  t <- (-log(stats::runif(nrow(m))) / h)^(1 / 1.5)                                     # Weibull PH, shape 1.5
  m$time <- round(pmin(t, 21), 2); m$event <- as.integer(t <= 21)
  spec <- list(handoff_version = "0.6.6", origin = "planned", organism = "test", system = "germination chamber", experimental_unit = "tray",
               subsamples_per_unit = seeds, factors = list(list(name = "lot", n_levels = 2L, levels = list("untreated", "primed"))),
               design = list(chosen = "crd", reps = trays / 2), primary_outcome = list(name = "days to germination", type = "survival", unit = "days", event_column = "event"),
               analysis = list(model = "y ~ treatment", terms = list(treatment = "lot", layout_columns = list("lot")),
                               goals = list("survival"), primary_goal = "survival", goal_terms = list(survival = list("treatment"))))
  if (!is.null(power)) spec$power <- power
  jsonlite::write_json(spec, file.path(d, "study_spec.json"), auto_unbox = TRUE, pretty = TRUE)
  utils::write.csv(tr, file.path(d, "field_layout.csv"), row.names = FALSE)
  utils::write.table(m, file.path(d, "sample_metadata_template.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
  e <- list(id = "outcome", data_type = "outcome", files = list(table = "sample_metadata_template.tsv"), column = "time", event_column = "event")
  if (truth) {
    jsonlite::write_json(list(outcome_model = list(family = "weibull_ph"), goal_truth = list(goal = "survival", model_terms = list("treatment"), status = "planted",
      true_effect = hr, true_log_hazard_ratio = log(hr), effect_scale = "hazard ratio, last level vs reference")), file.path(d, "sim_truth.json"), auto_unbox = TRUE, digits = NA)
    e$truth <- list(kind = "outcome", file = "sim_truth.json", goal = list(goal = "survival", field = "goal_truth", terms = list("treatment"), test = "wald"))
  }
  bioana_write_manifest(d, e)
  list(dir = d, data = m)
}

test_that("seeds on trays: subsamples kept (no averaging), tray = shared frailty, robust (Lin-Wei) score test = coxph(cluster); HR truth on one scale", {
  f <- seed_folder(power = list(target_hazard_ratio = 0.6, events_required = 10000))
  utils::capture.output(r <- bioana_analyze(f$dir, secondary = FALSE))
  a <- r$summary$analyses$outcome; st <- a$stat_test
  expect_equal(a$n_units_fitted, nrow(f$data))
  expect_match(paste(a$warnings, collapse = " "), "not averaged")
  expect_true(st$clustered); expect_equal(st$primary, "robust score (clustered)")
  m <- f$data; m$lot <- factor(m$lot, c("untreated", "primed"))
  rc <- summary(survival::coxph(survival::Surv(time, event) ~ lot + cluster(unit), data = m))$robscore
  expect_equal(trow(st, "robust score")[[1]]$statistic, unname(rc["test"]), tolerance = 1e-8)
  expect_match(st$note, "events observed, 10000 required")
  gr <- a$goal_recovery
  expect_equal(gr$true_effect, 0.6); expect_match(gr$scale, "hazard ratio")
  expect_equal(gr$estimate, exp(r$results$outcome$estimate[r$results$outcome$term == "treatmentprimed"]))
  expect_true(gr$sign_agreement)
})

block_surv_folder <- function(nblocks, per_block_line, hr = 0.7, seed = 1) {
  d <- file.path(tempfile("blk_"), "design"); dir.create(d, recursive = TRUE); set.seed(seed)
  g <- expand.grid(rep = seq_len(per_block_line), line = c("B73", "Mo17"), block = sprintf("B%02d", seq_len(nblocks)), stringsAsFactors = FALSE)
  g$unit <- sprintf("P%03d", seq_len(nrow(g))); g$sample_id <- g$unit
  be <- stats::rnorm(nblocks, 0, 0.3)[match(g$block, sprintf("B%02d", seq_len(nblocks)))]
  t <- (-log(stats::runif(nrow(g))) / (0.01 * exp(be) * ifelse(g$line == "Mo17", hr, 1)))^(1 / 1.5)
  g$time <- pmin(round(t, 1), 90); g$event <- as.integer(t <= 90)
  spec <- list(handoff_version = "0.6.6", origin = "planned", organism = "Zea mays", system = "field", factors = list(list(name = "line", n_levels = 2L, levels = list("B73", "Mo17"))),
               design = list(chosen = "rcbd"), primary_outcome = list(name = "days to flowering", type = "survival", event_column = "event"),
               analysis = list(model = "y ~ treatment + block", terms = list(treatment = "line", block = "block", layout_columns = list("line", "block"))))
  jsonlite::write_json(spec, file.path(d, "study_spec.json"), auto_unbox = TRUE, pretty = TRUE)
  utils::write.csv(g[c("unit", "line", "block")], file.path(d, "field_layout.csv"), row.names = FALSE)
  utils::write.table(g, file.path(d, "sample_metadata_template.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
  bioana_write_manifest(d, list(id = "outcome", data_type = "outcome", files = list(table = "sample_metadata_template.tsv"), column = "time", event_column = "event"))
  list(dir = d, data = g)
}

test_that("Cox with fixed blocks: few blocks -> strata(block) = coxph(strata); many small blocks -> shared frailty; stratified log-rank kept", {
  f <- block_surv_folder(4, 15)
  utils::capture.output(r <- bioana_analyze(f$dir, secondary = FALSE)); a <- r$summary$analyses$outcome
  expect_match(paste(a$warnings, collapse = " "), "fitted as strata\\(block\\)")
  m <- f$data; m$line <- factor(m$line, c("B73", "Mo17"))
  ref <- survival::coxph(survival::Surv(time, event) ~ line + strata(block), data = m)
  expect_equal(r$results$outcome$estimate[r$results$outcome$term == "treatmentMo17"], unname(stats::coef(ref)), tolerance = 1e-8)
  f2 <- block_surv_folder(30, 1)
  utils::capture.output(r2 <- bioana_analyze(f2$dir, secondary = FALSE)); a2 <- r2$summary$analyses$outcome
  expect_match(paste(a2$warnings, collapse = " "), "shared frailty \\(1 \\| block\\)")
  expect_true(a2$engine$engine %in% c("coxme", "coxph_frailty"))
  expect_true(any(grepl("stratified log-rank", vapply(a2$stat_test$table, `[[`, "", "statistic_name"))))
})

test_that("delayed effect (non-PH truth): goal detected by the prespecified weighted log-rank; estimate = piecewise Cox HR after change_time", {
  d <- file.path(tempfile("delay_"), "design"); dir.create(d, recursive = TRUE); set.seed(3)
  n <- 240; g <- rep(c("unprimed", "primed"), each = n / 2); t0 <- 6; h0 <- 0.08; hr <- 0.6
  e <- -log(stats::runif(n)); hz <- ifelse(g == "primed", hr, 1)
  t <- ifelse(e <= h0 * t0, e / h0, t0 + (e - h0 * t0) / (h0 * hz))
  m <- data.frame(unit = sprintf("S%03d", 1:n), sample_id = sprintf("S%03d", 1:n), priming = g, time = round(pmin(t, 21), 2), event = as.integer(t <= 21))
  spec <- list(handoff_version = "0.6.6", origin = "planned", organism = "test", system = "lab", factors = list(list(name = "priming", n_levels = 2L, levels = list("unprimed", "primed"))),
               design = list(chosen = "crd"), primary_outcome = list(name = "days", type = "survival", event_column = "event"),
               analysis = list(model = "y ~ treatment", terms = list(treatment = "priming", layout_columns = list("priming")), goals = list("survival"), primary_goal = "survival",
                               goal_terms = list(survival = list("treatment")), goal_effect = list(survival = list(size = 0.6, scale = "hazard_ratio", proportional_hazards = FALSE, test = "logrank_fh_0_1")),
                               survival = list(logrank_weight = "fh_0_1")))
  jsonlite::write_json(spec, file.path(d, "study_spec.json"), auto_unbox = TRUE, pretty = TRUE)
  utils::write.csv(m[c("unit", "priming")], file.path(d, "field_layout.csv"), row.names = FALSE)
  utils::write.table(m, file.path(d, "sample_metadata_template.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
  jsonlite::write_json(list(outcome_model = list(family = "delayed_effect", proportional_hazards = FALSE, change_time = t0),
                            goal_truth = list(goal = "survival", model_terms = list("treatment"), status = "planted", true_effect = hr, true_log_hazard_ratio = log(hr),
                                              proportional_hazards = FALSE, effect_scale = "hazard ratio, last level vs reference (after t = 6; 1 before)")),
                       file.path(d, "sim_truth.json"), auto_unbox = TRUE, digits = NA)
  bioana_write_manifest(d, list(id = "outcome", data_type = "outcome", files = list(table = "sample_metadata_template.tsv"), column = "time", event_column = "event",
                                truth = list(kind = "outcome", file = "sim_truth.json", goal = list(goal = "survival", field = "goal_truth", terms = list("treatment"), test = "wald"))))
  utils::capture.output(r <- bioana_analyze(d, secondary = FALSE)); gr <- r$summary$analyses$outcome$goal_recovery
  expect_match(gr$test, "weighted log-rank fh_0_1"); expect_false(gr$proportional_hazards)
  m$grp_ <- factor(m$priming, c("unprimed", "primed")); m$y_ <- m$time; m$event_ <- m$event
  expect_equal(gr$p, bioana:::.ba_weighted_logrank(m, "fh_0_1")$p)
  sp <- survival::survSplit(Surv(time, event) ~ priming, data = transform(m[c("time", "event", "priming")], priming = factor(priming, c("unprimed", "primed"))), cut = t0, episode = "ph")
  sp$after <- (sp$priming == "primed") * (sp$ph == 2); sp$before <- (sp$priming == "primed") * (sp$ph == 1)
  pc <- survival::coxph(survival::Surv(tstart, time, event) ~ before + after, data = sp)
  expect_equal(gr$estimate, unname(exp(stats::coef(pc)["after"])), tolerance = 1e-8)
  expect_equal(gr$true_effect, 0.6); expect_equal(gr$piecewise$affected_phase, "after")
})

test_that("a Cox fit that does not converge (separation) is flagged and its goal truth not scored", {
  d <- file.path(tempfile("sep_"), "design"); dir.create(d, recursive = TRUE)
  m <- data.frame(unit = sprintf("U%02d", 1:30), sample_id = sprintf("U%02d", 1:30), g = rep(c("a", "b"), each = 15),
                  time = c(seq(2, 30, 2), rep(40, 15)), event = c(rep(1L, 15), rep(0L, 15)))
  spec <- list(handoff_version = "0.6.6", origin = "planned", organism = "test", system = "lab", factors = list(list(name = "g", n_levels = 2L, levels = list("a", "b"))),
               design = list(chosen = "crd"), primary_outcome = list(name = "days", type = "survival", event_column = "event"),
               analysis = list(model = "y ~ treatment", terms = list(treatment = "g", layout_columns = list("g")), goals = list("survival"), primary_goal = "survival",
                               goal_terms = list(survival = list("treatment"))))
  jsonlite::write_json(spec, file.path(d, "study_spec.json"), auto_unbox = TRUE, pretty = TRUE)
  utils::write.csv(m[c("unit", "g")], file.path(d, "field_layout.csv"), row.names = FALSE)
  utils::write.table(m, file.path(d, "sample_metadata_template.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
  jsonlite::write_json(list(outcome_model = list(family = "weibull_ph"), goal_truth = list(goal = "survival", model_terms = list("treatment"), status = "planted",
                       true_effect = 0.5, true_log_hazard_ratio = log(0.5), effect_scale = "hazard ratio")), file.path(d, "sim_truth.json"), auto_unbox = TRUE)
  bioana_write_manifest(d, list(id = "outcome", data_type = "outcome", files = list(table = "sample_metadata_template.tsv"), column = "time", event_column = "event",
                                truth = list(kind = "outcome", file = "sim_truth.json", goal = list(goal = "survival", field = "goal_truth", terms = list("treatment"), test = "wald"))))
  utils::capture.output(r <- bioana_analyze(d, secondary = FALSE)); a <- r$summary$analyses$outcome
  expect_match(paste(a$warnings, collapse = " "), "did not converge")
  expect_false(a$goal_recovery$scored)
})

test_that("factorial_crd survival with blocks only in analysis.survival.strata: Cox and log-rank stratified, no block coefficient (Biodes 01d41ef)", {
  d <- file.path(tempfile("fcrd_"), "design"); dir.create(d, recursive = TRUE); set.seed(5)
  g <- expand.grid(rep = 1:2, line = c("B73", "Mo17"), n = c("low", "high"), block = sprintf("R%02d", 1:12), stringsAsFactors = FALSE)
  g$unit <- sprintf("P%03d", seq_len(nrow(g))); g$sample_id <- g$unit
  be <- stats::rnorm(12, 0, 0.4)[match(g$block, sprintf("R%02d", 1:12))]
  t <- -log(stats::runif(nrow(g))) / (0.03 * exp(be) * ifelse(g$line == "Mo17", 0.7, 1) * ifelse(g$n == "high", 1.3, 1))
  g$time <- pmin(round(t, 1), 90); g$event <- as.integer(t <= 90)
  spec <- list(handoff_version = "0.6.6", origin = "planned", organism = "Zea mays", system = "field",
               factors = list(list(name = "line", n_levels = 2L, levels = list("B73", "Mo17")), list(name = "n", n_levels = 2L, levels = list("low", "high"))),
               design = list(chosen = "factorial_crd"), primary_outcome = list(name = "days to flowering", type = "survival", event_column = "event"),
               analysis = list(model = "y ~ treatment * treatment2", terms = list(treatment = "line", treatment2 = "n", layout_columns = list("line", "n", "block")),
                               survival = list(strata = list("block"))))
  jsonlite::write_json(spec, file.path(d, "study_spec.json"), auto_unbox = TRUE, pretty = TRUE)
  utils::write.csv(g[c("unit", "line", "n", "block")], file.path(d, "field_layout.csv"), row.names = FALSE)
  utils::write.table(g, file.path(d, "sample_metadata_template.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
  bioana_write_manifest(d, list(id = "outcome", data_type = "outcome", files = list(table = "sample_metadata_template.tsv"), column = "time", event_column = "event"))
  utils::capture.output(r <- bioana_analyze(d, secondary = FALSE)); a <- r$summary$analyses$outcome
  expect_match(a$model_fitted, "strata\\(block\\)")
  expect_match(paste(a$warnings, collapse = " "), "stratified on block")
  expect_false(any(grepl("^block", r$results$outcome$term)))
  m <- g; m$line <- factor(m$line, c("B73", "Mo17")); m$n <- factor(m$n, c("low", "high"))
  ref <- survival::coxph(survival::Surv(time, event) ~ line * n + strata(block), data = m)
  expect_equal(unname(r$results$outcome$estimate), unname(stats::coef(ref)), tolerance = 1e-8)
  expect_true(any(grepl("stratified log-rank chi-squared \\(strata: block\\)", vapply(a$stat_test$table, `[[`, "", "statistic_name"))))
})
