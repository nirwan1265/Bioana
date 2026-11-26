# Analysis goals (docs/contracts/analysis_goals.md). Biodes owns the vocabulary
# (studyforge inst/registry/analysis_goals.json) and writes analysis.goals /
# primary_goal / goal_terms into study_spec.json. Bioana keys engine_selection.tsv
# (`goals` column) to those ids and turns a goal into a test plan:
#   wald: the coefficients tested one at a time (treatment_effect, trend, association,
#         survival, an interaction with one coefficient)
#   lrt:  a joint test of a set of coefficients, full vs reduced model
#         (time_course; an interaction with several coefficients)

#' Goals Bioana can test, per data type
#'
#' @return A data frame: `data_type`, `goal`, and whether the Biodes registry
#'   (`studyforge`, when installed) lists the goal.
#' @export
bioana_goals <- function() {
  s <- .ba_tsv("engine_selection")
  g <- do.call(rbind, lapply(seq_len(nrow(s)), function(i) data.frame(data_type = s$data_type[i], goal = strsplit(s$goals[i], ",")[[1]])))
  g <- unique(g); rownames(g) <- NULL
  reg <- .ba_goal_registry()
  g$in_biodes_registry <- if (is.null(reg)) NA else g$goal %in% reg$id
  g
}

# The Biodes goal registry, when studyforge is installed (NULL otherwise).
.ba_goal_registry <- function() {
  f <- if (requireNamespace("studyforge", quietly = TRUE)) system.file("registry", "analysis_goals.json", package = "studyforge") else ""
  if (!nzchar(f) || !file.exists(f)) return(NULL)
  x <- jsonlite::read_json(f, simplifyVector = FALSE)
  items <- if (!is.null(x$goals)) x$goals else x
  ids <- vapply(items, function(i) i$id %||% "", "")
  list(id = ids[nzchar(ids)], items = items)
}

# Which goal a data set is analysed for. First goal of analysis.goals (primary_goal first)
# that engine_selection.tsv has rows for; otherwise the default for the data type.
.ba_resolve_goal <- function(key, design, args = list()) {
  a <- design$spec$analysis %||% list()
  listed <- unique(c(unlist(a$primary_goal), unlist(a$goals)))
  s <- .ba_tsv("engine_selection")
  avail <- unique(unlist(strsplit(s$goals[s$data_type == key], ",", fixed = TRUE)))
  default <- if (key == "genotypes") "association" else "treatment_effect"
  asm <- design$spec$assumptions %||% list()
  # Biodes writes an object keyed by field path: assumptions["analysis.goals"] = "why"
  assumed <- any(grepl("^analysis\\.(primary_goal|goals)$", names(asm))) ||
    any(vapply(asm, function(x) is.list(x) && grepl("^analysis\\.(primary_goal|goals)$", x$field %||% ""), logical(1)))
  why <- asm[["analysis.goals"]] %||% asm[["analysis.primary_goal"]]
  note <- NULL
  if (!is.null(args$goal)) {
    if (!args$goal %in% avail) stop(sprintf("Goal '%s' is not available for %s. Available: %s.", args$goal, key, paste(avail, collapse = ", ")), call. = FALSE)
    id <- args$goal; src <- "params"
  } else if (length(listed)) {
    hit <- listed[listed %in% avail]
    if (length(hit)) { id <- hit[1]; src <- "study_spec" }
    else {
      id <- default; src <- "default"
      note <- sprintf("the study's goals (%s) do not apply to this data set in Bioana yet; tested %s", paste(listed, collapse = ", "), default)
    }
  } else {
    id <- default; src <- "default"
    vars <- design$parsed$vars
    if (key != "genotypes" && ("time" %in% vars || any(grepl(":", design$parsed$fixed_labels)))) {
      note <- sprintf("no analysis.goals in study_spec: tested %s, although the model has %s; Biodes should state the goal (time_course, interaction, ...)",
                      default, paste(c(if ("time" %in% vars) "a time term", if (any(grepl(":", design$parsed$fixed_labels))) "interaction terms"), collapse = " and "))
    }
  }
  reg <- .ba_goal_registry()
  unknown_to_biodes <- !is.null(reg) && !id %in% reg$id
  eff <- a$goal_effect[[id]]
  # Biodes 0.6.7: one-group studies (single_group, y ~ 1) write analysis.goals: [] -- no effect goal; the test comes
  # from expected_ratio / correlation / agreement. The id stays for engine selection only.
  pv <- design$parsed
  no_terms <- !length(pv$fixed_labels) && !length(pv$random)
  none <- is.null(args$goal) && (("goals" %in% names(a) && !length(listed)) || (no_terms && !length(listed)))
  list(id = id, source = src, assumed = assumed && src == "study_spec", assumed_why = if (assumed && src == "study_spec") why,
       terms = unlist(a$goal_terms[[id]]), effect = eff, listed = listed, note = if (length(note)) paste(note, collapse = "; "),
       unknown_to_biodes = unknown_to_biodes, none = none, adjust = .ba_adjust_vars(pv, design$spec))
}

# time_course within one group (goal_effect.time_course group / group_factor): the model term and level,
# or a note saying why time is tested over all groups.
.ba_within_group <- function(parsed, eff, meta) {
  gf <- eff$group_factor; g <- eff$group
  if (is.null(g) || is.null(gf)) return("goal_effect scope 'within_group' names no group (group, group_factor): time tested over all groups")
  hit_v <- vapply(parsed$vars, function(z) identical(z, gf) || identical(as.character(parsed$col[[z]]), gf), logical(1))
  if (!any(hit_v)) return(sprintf("within_group: group_factor '%s' is not a term or column of analysis.model: time tested over all groups", gf))
  v <- parsed$vars[hit_v][1]
  lev <- unique(c(parsed$levels[[v]], if (!is.null(meta)) as.character(meta[[parsed$col[[v]]]])))
  lev <- lev[!is.na(lev)]
  hit <- lev[tolower(lev) == tolower(g)]
  if (!length(hit)) hit <- lev[grepl(tolower(g), tolower(lev), fixed = TRUE)]
  if (length(hit) != 1L) return(sprintf("within_group: group '%s' is not one level of %s (%s): time tested over all groups", g, v, paste(lev, collapse = ", ")))
  list(term = v, column = as.character(parsed$col[[v]]), level = hit, group_factor = gf)
}

# Design copy for the goal: trend fits its term as numeric; time_course within a group makes that
# group the reference level.
.ba_goal_design <- function(design, goal) {
  if (identical(goal$id, "trend")) {
    tv <- goal$terms[1]
    if (is.null(tv) || is.na(tv)) stop("Goal 'trend' needs analysis.goal_terms$trend: the numeric term.", call. = FALSE)
    if (!tv %in% design$parsed$vars) stop(sprintf("Goal 'trend': term '%s' is not in analysis.model (%s).", tv, design$model), call. = FALSE)
    design$parsed$numeric <- unique(c(design$parsed$numeric, tv))
    vv <- design$parsed$values[[tv]]
    unit <- if (!is.null(vv) && nzchar(vv$unit)) paste0(" ", vv$unit) else ""
    if (!is.null(vv) && identical(vv$scale, "log")) {
      gr <- .ba_rules("goal_rules")
      base <- as.numeric(vv$log_base %||% gr$trend_log_base_default %||% 10)
      off <- vv$log_offset; src <- "factors[].log_offset"
      if (is.null(off)) { off <- if (any(vv$values <= 0)) min(vv$values[vv$values > 0]) else 0; src <- "the smallest positive value (goal_rules.json)" }
      off <- as.numeric(off)
      if (any(vv$values + off <= 0)) stop(sprintf("Goal 'trend': log scale with log_offset %g leaves values <= 0 for '%s'.", off, tv), call. = FALSE)
      vv$transform <- function(x) log(x + off, base)
      design$parsed$values[[tv]] <- vv
      lg <- if (base == 10) "log10" else if (base == 2) "log2" else if (abs(base - exp(1)) < 1e-9) "ln" else sprintf("log%g", base)
      goal$predictor <- list(scale = "log", log_base = base, log_offset = off)
      goal$scale_note <- if (off > 0) sprintf("slope per %s unit of %s + %g%s (offset: %s)", lg, tv, off, unit, src)
                         else sprintf("slope per %s unit of %s%s", lg, tv, unit)
    } else {
      goal$predictor <- list(scale = "linear")
      goal$scale_note <- sprintf("slope per%s of %s (linear scale)", if (nzchar(unit)) unit else " unit", tv)
    }
    pu <- goal$effect$per_unit
    if (!is.null(pu)) {
      goal$scale_note <- sprintf("%s; goal_effect size per %s", goal$scale_note, pu)
      log_pu <- grepl("^log", pu)
      if (log_pu != identical(goal$predictor$scale, "log"))
        goal$note <- paste(c(goal$note, sprintf("goal_effect.trend.per_unit '%s' does not match factor scale '%s': Bioana follows the factor scale",
                                              pu, goal$predictor$scale)), collapse = "; ")
    }
  }
  if (identical(goal$id, "time_course") && identical(goal$effect$scope, "within_group")) {
    w <- .ba_within_group(design$parsed, goal$effect, design$meta)
    if (is.character(w)) goal$note <- paste(c(goal$note, w), collapse = "; ") else {
      design$parsed$levels[[w$term]] <- c(w$level, setdiff(design$parsed$levels[[w$term]], w$level))
      goal$within <- w
    }
  }
  design$goal <- goal
  design
}

# Model-matrix column -> the model term it belongs to.
.ba_column_terms <- function(rhs, frame) {
  X <- stats::model.matrix(rhs, data = frame)
  labs <- c("(Intercept)", attr(stats::terms(rhs), "term.labels"))
  stats::setNames(labs[attr(X, "assign") + 1L], colnames(X))
}

# The test plan for a goal, given the model matrix columns (after estimability) and their terms.
.ba_goal_plan <- function(goal, cols, col_term, p) {
  id <- goal$id %||% "treatment_effect"
  vars_of <- function(t) if (t == "(Intercept)") character() else all.vars(str2lang(t))
  ct <- col_term[cols]
  treat <- cols[vapply(cols, function(t) .ba_coef_role(t, p), "") == "treatment"]
  if (id %in% c("treatment_effect", "association", "survival")) {
    # Biodes goal_terms (e.g. ["treatment"]): exactly those terms; adjustment factors (analysis.adjust_for,
    # factors[].role adjustment) are never part of the goal (Biodes 412c0e6, Simitall a1503f4)
    adj_term <- vapply(ct[treat], function(t) any(all.vars(str2lang(t)) %in% goal$adjust), logical(1))
    k <- treat[!adj_term]
    if (length(goal$terms)) { kk <- k[ct[k] %in% goal$terms]; if (length(kk)) k <- kk }
    return(list(test = "wald", coefs = k, drop = character(), drop_terms = character(), label = NULL, tested_terms = unique(ct[k])))
  }
  if (id == "interaction") {
    tt <- goal$terms
    if (length(tt) != 2) stop("Goal 'interaction' needs analysis.goal_terms$interaction: two model terms, e.g. [\"treatment\", \"time\"].", call. = FALSE)
    terms_hit <- unique(ct[vapply(ct, function(t) setequal(vars_of(t), tt) && grepl(":", t), logical(1))])
    if (!length(terms_hit)) stop(sprintf("Goal 'interaction': analysis.model (%s) has no %s:%s term.", deparse(p$formula), tt[1], tt[2]), call. = FALSE)
    k <- cols[ct %in% terms_hit]
    lab <- paste(tt, collapse = ":")
    if (length(k) == 1L) return(list(test = "wald", coefs = k, drop = character(), drop_terms = character(), label = NULL, tested_terms = terms_hit))
    return(list(test = "lrt", coefs = character(), drop = k, drop_terms = terms_hit, label = sprintf("interaction %s (%d coefficients)", lab, length(k)), tested_terms = terms_hit))
  }
  if (id == "time_course") {
    tv <- goal$terms[1] %||% "time"
    if (is.na(tv)) tv <- "time"
    terms_hit <- unique(ct[vapply(ct, function(t) tv %in% vars_of(t), logical(1))])
    if (!length(terms_hit)) stop(sprintf("Goal 'time_course': analysis.model has no term with '%s'.", tv), call. = FALSE)
    w <- goal$within; pnote <- NULL
    if (!is.null(w)) {
      inter <- unique(ct[vapply(ct, function(t) all(c(tv, w$term) %in% vars_of(t)) && grepl(":", t), logical(1))])
      k <- cols[ct == tv]
      # group = reference level: the time main-effect coefficients are the time effects in that group
      if (length(inter) && length(k)) return(list(test = "lrt", coefs = character(), drop = k, drop_terms = tv, tested_terms = inter, within = w,
        label = sprintf("time_course within %s = %s: %s (%d coefficients)", w$term, w$level, tv, length(k))))
      pnote <- sprintf("within_group: analysis.model has no %s:%s interaction, so time has the same effect in every group: tested over all groups", w$term, tv)
    }
    k <- cols[ct %in% terms_hit]
    return(list(test = "lrt", coefs = character(), drop = k, drop_terms = terms_hit, tested_terms = terms_hit, note = pnote,
                label = sprintf("time_course: %s (%d coefficients)", paste(terms_hit, collapse = ", "), length(k))))
  }
  if (id == "regression") {
    # slopes of the covariates (role covariate): the goal's terms, else every covariate main effect
    cv <- p$vars[vapply(p$vars, function(v) identical(p$role[[v]], "covariate"), logical(1))]
    tv <- intersect(goal$terms %||% cv, cv); if (!length(tv)) tv <- cv
    if (!length(tv)) stop("Goal 'regression': analysis.model has no covariate term (covariate, covariate2, ...).", call. = FALSE)
    k <- cols[ct %in% tv]
    if (!length(k)) stop(sprintf("Goal 'regression': no main-effect coefficient for %s.", paste(tv, collapse = ", ")), call. = FALSE)
    return(list(test = "wald", coefs = k, drop = character(), drop_terms = character(), label = NULL, tested_terms = unique(ct[k])))
  }
  if (id == "trend") {
    tv <- goal$terms[1]
    k <- cols[ct == tv]
    if (!length(k)) stop(sprintf("Goal 'trend': analysis.model has no main effect of '%s'.", tv), call. = FALSE)
    return(list(test = "wald", coefs = k, drop = character(), drop_terms = character(), label = NULL, tested_terms = tv))
  }
  stop(sprintf("Goal '%s' has no test plan in Bioana yet.", id), call. = FALSE)
}

# Summary block for the goal of one data set.
.ba_goal_summary <- function(goal, plan) {
  list(id = goal$id, source = goal$source, assumed_by_biodes = isTRUE(goal$assumed), assumed_why = goal$assumed_why, test = plan$test,
       effect = goal$effect, scale_note = goal$scale_note, predictor = goal$predictor,
       coefficients = if (plan$test == "wald") as.list(plan$coefs) else as.list(plan$drop),
       dropped_terms = if (plan$test == "lrt") as.list(plan$drop_terms), tested_terms = as.list(plan$tested_terms),
       within = plan$within, label = plan$label,
       note = if (length(c(goal$note, plan$note))) paste(c(goal$note, plan$note), collapse = "; "),
       unknown_to_biodes_registry = if (isTRUE(goal$unknown_to_biodes)) TRUE)
}

# ---- goal recovery (simulated data): Bioana's goal test against Simitall's goal_truth ----
# Keyed to the data set's manifest truth.goal and sim_truth.json goal_truth (model_terms, status).

.ba_goal_recovery_base <- function(gt, g) {
  vs <- function(x) lapply(unlist(x), function(z) sort(all.vars(str2lang(z))))
  truth_terms <- unlist(gt$model_terms); tested <- unlist(g$tested_terms)
  match <- length(truth_terms) > 0 && length(tested) > 0 &&
    all(vapply(vs(truth_terms), function(a) any(vapply(vs(tested), identical, logical(1), a)), logical(1)))
  out <- list(goal = gt$goal %||% NA_character_, bioana_goal = g$id, status = gt$status %||% "planted",
              model_terms = as.list(truth_terms), tested_terms = as.list(tested), terms_match = match, test = g$test)
  if (!identical(out$status, "planted")) { out$scored <- FALSE; out$note <- sprintf("goal truth %s: not scored", out$status) }
  else if (!identical(gt$goal, g$id)) { out$scored <- FALSE; out$note <- sprintf("Simitall planted %s, Bioana tested %s: not scored", gt$goal, g$id) }
  else out$scored <- TRUE
  if (!is.null(g$within)) {
    out$within <- g$within$level
    if (length(gt$groups_changing)) out$group_match <- g$within$level %in% unlist(gt$groups_changing)
  }
  if (identical(g$id, "trend") && length(gt$values) && length(gt$x)) {
    pr <- g$predictor %||% list(scale = "linear")
    raw <- as.numeric(unlist(gt$values))
    mine <- if (identical(pr$scale, "log")) log(raw + pr$log_offset, pr$log_base) else raw
    out$predictor_match <- isTRUE(all.equal(unname(mine), unname(as.numeric(unlist(gt$x))), tolerance = 1e-3))
    if (!out$predictor_match) out$note <- paste(c(out$note, sprintf("Simitall planted the slope on another predictor (x = %s for values %s); Bioana fitted %s: slopes not comparable",
      paste(signif(as.numeric(unlist(gt$x)), 4), collapse = ", "), paste(raw, collapse = ", "), g$scale_note)), collapse = "; ")
  }
  out
}

# Simitall hidden factors for one assay: scores (units x factors) and loadings (features x factors).
.ba_hidden_factors <- function(dir, id) {
  st <- file.path(dir, "sim_truth.json")
  if (!file.exists(st)) return(NULL)
  t <- tryCatch(jsonlite::read_json(st, simplifyVector = FALSE), error = function(e) NULL)
  hf <- t$assays[[id]]$counts$hidden_factors %||% t$assays[[paste0("assay_", id)]]$counts$hidden_factors %||% t$counts$hidden_factors
  if (is.null(hf$scores) || is.null(hf$loadings)) return(NULL)
  rd <- function(f) { f <- file.path(dir, f); if (!file.exists(f)) return(NULL)
    x <- utils::read.delim(f, check.names = FALSE, stringsAsFactors = FALSE); m <- as.matrix(x[, -1, drop = FALSE]); rownames(m) <- x[[1]]; m }
  S <- rd(hf$scores); L <- rd(hf$loadings)
  if (is.null(S) || is.null(L)) return(NULL)
  k <- intersect(colnames(S), colnames(L))
  list(S = S[, k, drop = FALSE], L = L[, k, drop = FALSE])
}

# The design model fitted to each feature's hidden part (L %*% t(S)): coefficients x features (log2).
.ba_hidden_goal_effect <- function(hf, gd, cols) {
  ids <- gd$ids
  if (!all(ids %in% rownames(hf$S))) return(NULL)
  C <- hf$L %*% t(hf$S[ids, , drop = FALSE])
  X <- gd$X
  B <- qr.coef(qr(X), t(C))
  B[is.na(B)] <- 0
  if (!length(intersect(cols, rownames(B)))) return(NULL)
  B
}

# Outcome: the goal test detects the planted effect (detected), and a one-coefficient goal's estimate
# against true_effect (continuous outcomes: bias and CI coverage; other types: sign).
.ba_goal_recovery_outcome <- function(gt, g, res, alpha, type) {
  out <- .ba_goal_recovery_base(gt, g)
  if (!isTRUE(out$scored)) return(out)
  coefs <- unlist(g$coefficients)
  p <- if (identical(g$test, "lrt")) g$result$p %||% NA_real_ else {
    r <- res[res$term %in% coefs, , drop = FALSE]; if (nrow(r)) min(r$padj, na.rm = TRUE) else NA_real_
  }
  out$p <- p; out$detected <- isTRUE(p < alpha)
  te <- suppressWarnings(as.numeric(unlist(gt$true_effect)))
  if (identical(type, "survival") && identical(g$test, "wald") && length(coefs) == 1L) {
    # Cox coefficient = log hazard ratio: compared with Simitall's true_log_hazard_ratio (or log of a hazard-ratio truth)
    tl <- suppressWarnings(as.numeric(unlist(gt$true_log_hazard_ratio)))
    if (!length(tl) && grepl("hazard ratio", gt$effect_scale %||% "") && length(te) == 1L) tl <- log(te)
    if (length(tl) == 1L) {
      r <- res[res$term == coefs, , drop = FALSE]
      out$estimate <- exp(r$estimate); out$true_effect <- exp(tl); out$scale <- "hazard ratio (compared on the log scale, Cox coefficient)"
      out$bias <- r$estimate - tl; out$covered <- abs(r$estimate - tl) <= stats::qnorm(1 - alpha / 2) * r$se
      out$sign_agreement <- sign(r$estimate) == sign(tl)
      return(out)
    }
  }
  if (identical(type, "aft") && identical(g$test, "wald") && length(coefs) == 1L && length(te) == 1L) {
    # time ratio: compared on the log scale (survreg coefficient and SE)
    r <- res[res$term == coefs, , drop = FALSE]
    out$estimate <- exp(r$estimate); out$true_effect <- te
    out$bias <- r$estimate - log(te); out$covered <- abs(r$estimate - log(te)) <= stats::qnorm(1 - alpha / 2) * r$se
    out$sign_agreement <- sign(r$estimate) == sign(log(te))
    return(out)
  }
  if (identical(g$test, "wald") && length(coefs) == 1L && length(te) == 1L && !isFALSE(out$predictor_match)) {
    r <- res[res$term == coefs, , drop = FALSE]
    out$estimate <- r$estimate; out$true_effect <- te
    if (type == "continuous" && !is.na(r$se)) {
      q <- if (!is.null(r$df) && is.finite(r$df)) stats::qt(1 - alpha / 2, r$df) else stats::qnorm(1 - alpha / 2)
      out$bias <- r$estimate - te; out$covered <- abs(r$estimate - te) <= q * r$se
    } else out$sign_agreement <- sign(r$estimate) == sign(te)
  }
  out
}

# Features: share of Simitall's goal features (truth.goal file: feature_id, term, effect, change) that
# the goal test finds (tpr); one-coefficient goals also compare estimates with the planted effects.
.ba_goal_recovery_features <- function(entry, dir, g, table, alpha, gd = NULL) {
  tg <- entry$truth$goal
  f <- file.path(dir, tg$file %||% "")
  if (is.null(tg$file) || !file.exists(f)) return(list(goal = tg$goal, scored = FALSE, note = "truth.goal file missing: not scored"))
  tf <- utils::read.delim(f, stringsAsFactors = FALSE, check.names = FALSE)
  st <- file.path(dir, "sim_truth.json")
  gt <- if (file.exists(st)) tryCatch(jsonlite::read_json(st, simplifyVector = FALSE)$goal_truth, error = function(e) NULL) else NULL
  if (is.null(gt)) gt <- list(goal = tg$goal, model_terms = list(tg$term), status = "planted")
  out <- .ba_goal_recovery_base(gt, g)
  if (!isTRUE(out$scored)) return(out)
  rows <- if (identical(g$test, "lrt")) table[table$contrast == g$label, , drop = FALSE] else table[table$contrast %in% unlist(g$coefficients), , drop = FALSE]
  pad <- rows$padj; pad[is.na(pad)] <- 1
  sig <- tapply(pad, rows$feature_id, min) < alpha
  goal_ids <- unique(tf$feature_id); tested <- intersect(goal_ids, names(sig))
  out$n_goal_features <- length(goal_ids); out$n_tested <- length(tested)
  out$found <- sum(sig[tested]); out$tpr <- if (length(tested)) mean(sig[tested]) else NA_real_
  out$n_called <- sum(sig); out$called_outside_goal <- sum(sig[setdiff(names(sig), goal_ids)])
  # split the rest: features Simitall planted with another effect (truth.file) vs features with no planted effect
  of <- file.path(dir, entry$truth$file %||% "")
  if (!is.null(entry$truth$file) && file.exists(of)) {
    other <- setdiff(unique(utils::read.delim(of, stringsAsFactors = FALSE)$feature_id), goal_ids)
    out$called_other_planted <- sum(sig[intersect(names(sig), other)])
    out$called_not_planted <- sum(sig[setdiff(names(sig), c(goal_ids, other))])
    out$called_hidden_factor <- 0L
    hf <- .ba_hidden_factors(dir, entry$id)
    np <- setdiff(names(sig)[sig], c(goal_ids, other))
    if (!is.null(hf) && !is.null(gd) && length(np)) {
      gr <- .ba_rules("goal_rules")
      cols <- unlist(g$coefficients)
      B <- .ba_hidden_goal_effect(hf, gd, cols)
      if (!is.null(B)) {
        np <- intersect(np, colnames(B))
        if (identical(g$test, "wald") && length(cols) == 1L) {
          est <- rows$estimate[match(np, rows$feature_id)]; hb <- B[cols, np]
          hit <- np[!is.na(est) & sign(hb) == sign(est) & abs(hb) >= (gr$hidden_factor_share %||% 0.5) * abs(est)]
          out$hidden_factor_rule <- sprintf("hidden-factor effect on %s has the estimate's sign and >= %g of its size", cols, gr$hidden_factor_share %||% 0.5)
        } else {
          big <- apply(abs(B[intersect(cols, rownames(B)), np, drop = FALSE]), 2, max)
          hit <- np[big >= (gr$hidden_factor_effect_large %||% 0.25)]
          out$hidden_factor_rule <- sprintf("largest hidden-factor effect on the tested coefficients >= %g log2", gr$hidden_factor_effect_large %||% 0.25)
        }
        out$called_hidden_factor <- length(hit)
        out$called_not_planted <- out$called_not_planted - length(hit)
      }
    }
    out$false_call_rate <- if (out$n_called) out$called_not_planted / out$n_called else 0
  }
  out$effect_scale <- tg$effect_scale
  if (identical(g$test, "wald") && length(unlist(g$coefficients)) == 1L && !isFALSE(out$predictor_match) && "effect" %in% names(tf)) {
    m <- rows[match(tested, rows$feature_id), , drop = FALSE]; e <- tf$effect[match(tested, tf$feature_id)]
    ok <- !is.na(m$estimate) & !is.na(e)
    if (sum(ok) >= 3) {
      if (stats::sd(e[ok]) > 0 && stats::sd(m$estimate[ok]) > 0) out$estimate_cor <- stats::cor(m$estimate[ok], e[ok])
      out$sign_agreement <- mean(sign(m$estimate[ok]) == sign(e[ok]))
    }
  }
  out
}

# LRT between two outcome fits: statistic, df, p (anova methods differ in column names).
.ba_lrt <- function(fit0, fit) {
  a <- if (inherits(fit, "glm") && !inherits(fit, "negbin")) stats::anova(fit0, fit, test = "LRT") else stats::anova(fit0, fit)
  a <- as.data.frame(a)
  pcol <- grep("^Pr", names(a), value = TRUE)[1]
  scol <- grep("^(F|Chisq|Deviance|LR|LR\\.stat|LR stat\\.)$", names(a), value = TRUE)[1]
  dcol <- grep("^(Df|Chi Df|df|Resid\\. Df)$", names(a), value = TRUE)
  dcol <- if ("Df" %in% dcol) "Df" else if ("Chi Df" %in% dcol) "Chi Df" else dcol[1]
  last <- nrow(a)
  list(statistic = as.numeric(a[last, scol]), df = abs(as.numeric(a[last, dcol])), p = as.numeric(a[last, pcol]),
       method = sprintf("anova(reduced, full)%s: %s", if (inherits(fit, "glm") && !inherits(fit, "negbin")) ", test = 'LRT'" else "", scol))
}

# Class probabilities per treatment group (nominal, and binary on the logit scale): Simitall's
# outcome_model.group_probability {group: {category: p}} (binary: {group: P(y = 1)}), groups = the treatment cells
# "a:b" in the model's treatment order. Estimate = the observed share in the group, with a Clopper-Pearson interval.
# Simitall's probabilities are at an average unit; with block or unit effects the observed (marginal) share is
# slightly closer to the overall mean.
.ba_group_probability_recovery <- function(gp, frame, p, alpha = 0.05) {
  if (!length(gp)) return(NULL)
  tv <- intersect(p$vars[vapply(p$vars, function(v) identical(p$role[[v]] %||% "treatment", "treatment"), logical(1))], names(frame))
  tv <- setdiff(tv, p$numeric)
  grp <- if (length(tv)) do.call(paste, c(lapply(frame[tv], as.character), sep = ":")) else rep("all", nrow(frame))
  y <- if (is.factor(frame$y_)) as.character(frame$y_) else as.character(frame$y_)
  ok <- !is.na(y)
  rows <- list()
  for (g in names(gp)) {
    sel <- ok & grp == g; n <- sum(sel)
    if (!n) next
    tp <- gp[[g]]
    if (!is.list(tp)) tp <- list(`1` = tp)          # binary: P(y = 1)
    for (k in names(tp)) {
      x <- sum(y[sel] == k); ci <- stats::binom.test(x, n, conf.level = 1 - alpha)$conf.int
      rows[[length(rows) + 1L]] <- data.frame(group = g, category = k, n = n, true_probability = as.numeric(tp[[k]]), estimate = x / n,
                                              ci_low = ci[1], ci_high = ci[2], stringsAsFactors = FALSE)
    }
  }
  if (!length(rows)) return(list(scored = FALSE, note = sprintf("group_probability groups (%s) not found in the data (treatment cells %s)",
                                                                  paste(names(gp), collapse = ", "), paste(utils::head(unique(grp), 4), collapse = ", "))))
  tb <- do.call(rbind, rows)
  tb$covered <- tb$true_probability >= tb$ci_low & tb$true_probability <= tb$ci_high
  list(scored = TRUE, table = tb, max_abs_diff = max(abs(tb$estimate - tb$true_probability)), covered = sum(tb$covered), cells = nrow(tb),
       note = "observed share per group vs Simitall's class probability (at an average unit); Clopper-Pearson interval")
}

# Segregation ratio truth (Simitall outcome_model.true_ratio, expected_ratio, distorted) against the goodness-of-fit
# test: detected = test rejects the expected ratio (p < alpha), correct when that matches `distorted`; each class's
# observed share vs its true share (Clopper-Pearson).
.ba_ratio_recovery <- function(om, st, frame, alpha = 0.05) {
  g <- if (identical(st$test, "goodness_of_fit")) st else Filter(function(x) identical(x$test, "goodness_of_fit"), st$also %||% list())[1][[1]]
  tr <- unlist(om$true_ratio); cats <- names(tr)
  y <- as.character(frame$y_); y <- y[!is.na(y)]; n <- length(y)
  tb <- do.call(rbind, lapply(cats, function(k) { x <- sum(y == k); ci <- stats::binom.test(x, n, conf.level = 1 - alpha)$conf.int
    data.frame(category = k, n = n, true_share = tr[[k]], expected_share = as.numeric(unlist(om$expected_ratio)[k]) / sum(as.numeric(unlist(om$expected_ratio))),
               estimate = x / n, ci_low = ci[1], ci_high = ci[2], stringsAsFactors = FALSE) }))
  tb$covered <- tb$true_share >= tb$ci_low & tb$true_share <= tb$ci_high
  pv <- if (!is.null(g)) { r <- Filter(function(x) identical(x$contrast, "all"), g$table); if (length(r)) r[[1]]$p else g$p } else NA_real_
  dist <- isTRUE(om$distorted)
  list(goal = "segregation_ratio", model_terms = list(), scored = !is.null(g), test = "goodness_of_fit",
       p = pv, detected = isTRUE(pv < alpha), distorted = dist, correct = isTRUE(pv < alpha) == dist,
       covered = all(tb$covered), largest_difference = max(abs(tb$estimate - tb$true_share)), table = tb,
       note = sprintf("%s ratio %s; %s", if (dist) "distorted: true" else "true", paste(round(tr, 3), collapse = ":"),
                      if (is.null(g)) "goodness-of-fit not run" else if (isTRUE(pv < alpha) == dist) "test call correct" else if (dist) "distortion missed" else "false rejection"))
}

# Simitall's binary logit truth: conditional_odds_ratios {term: {level: OR vs the reference}}, conditional on the
# unit / block effects (random_effect_sd_logit). Compared with the model's log odds ratio (glmer with a random block
# is conditional on the block; glm without unit effects is the same OR) and its Wald interval.
.ba_odds_ratio_recovery <- function(cor, res, alpha = 0.05, conditional = TRUE) {
  rows <- list()
  for (v in names(cor)) { lv <- cor[[v]]
    for (l in names(lv)[-1]) { r <- res[res$term == paste0(v, l), , drop = FALSE]
      if (!nrow(r)) next
      q <- stats::qnorm(1 - alpha / 2); tv <- as.numeric(lv[[l]])
      rows[[length(rows) + 1L]] <- data.frame(term = r$term, true_or = tv, estimate = exp(r$estimate), ci_low = exp(r$estimate - q * r$se), ci_high = exp(r$estimate + q * r$se),
                                              covered = abs(r$estimate - log(tv)) <= q * r$se, stringsAsFactors = FALSE) } }
  if (!length(rows)) return(NULL)
  list(table = do.call(rbind, rows), scale = sprintf("odds ratio (Simitall conditional_odds_ratios; the model's OR is %s)",
                                                       if (conditional) "conditional on the random block / unit effect" else "marginal: no unit effects in the model"))
}

# Regression goal truth (Simitall goal_truth.true_slopes {term: slope per unit} for continuous outcomes, or
# true_log_odds_ratios {term: log OR per unit} for binary): each covariate's coefficient, its interval, covered.
.ba_slope_recovery <- function(gt, res, alpha, type, p) {
  tv <- unlist(if (type == "binary") gt$true_log_odds_ratios %||% gt$true_slopes else gt$true_slopes)
  if (!length(tv)) return(NULL)
  rows <- lapply(names(tv), function(t) {
    r <- res[res$term == t, , drop = FALSE]; if (!nrow(r)) return(NULL)
    q <- if (type == "continuous" && is.finite(r$df[1])) stats::qt(1 - alpha / 2, r$df[1]) else stats::qnorm(1 - alpha / 2)
    lo <- r$estimate - q * r$se; hi <- r$estimate + q * r$se; tr <- as.numeric(tv[[t]])
    f <- if (type == "binary") exp else identity
    data.frame(term = t, predictor = as.character(p$values[[t]]$name %||% t), true = f(tr), estimate = f(r$estimate), ci_low = f(lo), ci_high = f(hi),
               covered = tr >= lo & tr <= hi, p = r$p, stringsAsFactors = FALSE)
  })
  do.call(rbind, Filter(Negate(is.null), rows))
}
