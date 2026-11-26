# Survival outcomes: what Bioana reads, and what Biodes and Simitall should write (2026-11-26)

Bioana 0.10.0 runs a full survival analysis for `primary_outcome.type = "survival"`
(`docs/engines/survival.md`). It runs without anything new; these fields make it match the study.

## Biodes

1. Endpoint (in `primary_outcome`): `unit` (days, weeks), `event` (what counts as the event, e.g.
   "radicle >= 2 mm", "first lesion"), `origin` (sowing, inoculation, transplant), `censoring`
   (end of the experiment = administrative; lost or damaged = loss to follow-up), `event_column`
   (default "event": 1 event, 0 censored).
2. Competing events: `event_codes` {"0": "censored", "1": <event of interest>, "2": <competing>}
   when something else can end the observation (seed rot before germination, death before
   flowering). Bioana then runs cumulative incidence, Gray's test and Fine-Gray, and the Cox model
   becomes cause-specific.
3. `analysis.survival` (optional, only what the question says):
   - `landmarks`: times named in the question ("germination by day 7 and 14").
   - `rmst_tau`: the horizon when the question has one ("within the first 21 days").
   - `logrank_weight`: only for a stated early or delayed effect (early -> `fh_1_0`, delayed ->
     `fh_0_1`); prespecified, one weight.
   - `aft`: true when the question asks how much longer or shorter (time ratios), or for
     extrapolation.
   - `strata`: randomization strata that are not blocks in the model.
   - `life_table`: interval width when the event is scored at intervals (checked every 2 days).
4. Clusters (trays, litters, plots): a random intercept `(1 | block)` gives a shared frailty;
   with few clusters (< ~10) put them in `analysis.survival.strata` instead.
5. Power: size on events, not units (Schoenfeld: d = (z_a + z_b)^2 / (p (1 - p) (log HR)^2), 247
   events for HR 0.7, 80% power, alpha 0.05, 1:1), then N = d / P(event), inflated for loss. Please
   write `power.target_hazard_ratio`, `power.events_required`, `power.event_probability`,
   `power.loss_fraction`; Bioana will report observed events against `events_required`.
6. Goal: the registry's `survival` goal. Please say in `goal_effect.survival` which estimand the
   size is on: `scale = "hazard_ratio"` or `"time_ratio"` (or `"rmst_difference"` with `tau`).

## Read since (Bioana 0.10.1)

- Biodes 0.6.6: all fields above. Subsamples in a unit (seeds on a tray: `subsamples_per_unit` > 1,
  `unit_id` column): Bioana keeps every seed, makes the unit a shared frailty and uses the robust
  (Lin-Wei) score test as the primary test; Biodes need not write a frailty term for it.
- Simitall 3f82312: hazard-ratio truth is compared with the Cox coefficient on the log scale
  (`true_log_hazard_ratio`), time-ratio truth with survreg of `outcome_model.family`.

## Read since (Bioana 0.11.2)

- Biodes 825068a `goal_effect.survival.proportional_hazards: false`, `test: "logrank_fh_0_1"`; Simitall
  76d66c0 `outcome_model.change_time` / family: non-PH goals are scored with the weighted log-rank and a
  piecewise Cox HR (see NEWS 0.11.2).
- RCBD blocks as fixed terms in Cox models: Bioana fits >= 10 blocks as a shared frailty, fewer as strata.
  Biodes may keep writing `y ~ treatment + block`.

- Biodes 01d41ef: survival specs without blocks use `crd` / `factorial_crd`; blocks listed in
  `analysis.survival.strata` are not in the model. Bioana stratifies both the log-rank and the Cox model on them
  (Bioana 0.11.3).

## Simitall

1. Survival is now log-normal AFT only (hazards not proportional). Bioana scores that goal truth
   with `survreg(dist = "lognormal")` (time ratio) and the Cox HR is reported but not scored.
   Please add `suite_outcome survival_model`: `weibull_ph` (shape from the spec or 1.5; both PH and
   AFT: truth gives the HR and the time ratio), `exponential`, `lognormal_aft` (now), and non-PH
   cases for the weighted log-rank and RMST: `delayed_effect` (HR 1 until t0, then HR) and
   `crossing`.
2. Censoring: administrative at the end of follow-up (now: a quantile) plus random loss
   (`loss_rate`), recorded separately in `sim_truth.json`.
3. Competing events when the spec has `event_codes`: a cause-specific hazard per code; truth the
   true CIF per group at the landmarks and the cause-specific HRs.
4. Clusters: a gamma or log-normal frailty on the block term, truth its variance.
5. Truth for scoring: per group the true median, S(t) at `analysis.survival.landmarks`, RMST at
   tau, the true HR (PH models) or time ratio (AFT), and `proportional_hazards: true/false`.
   `outcome_model.family` names the model (Bioana matches the AFT distribution to it).
6. Interval scoring (`analysis.survival.life_table`): event times rounded up to the scoring grid.
