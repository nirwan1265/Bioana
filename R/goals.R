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
  list(id = id, source = src, assumed = assumed && src == "study_spec", assumed_why = if (assumed && src == "study_spec") why,
       terms = unlist(a$goal_terms[[id]]), effect = eff, listed = listed, note = if (length(note)) paste(note, collapse = "; "),
       unknown_to_biodes = unknown_to_biodes)
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
    return(list(test = "wald", coefs = treat, drop = character(), drop_terms = character(), label = NULL, tested_terms = unique(ct[treat])))
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
.ba_goal_recovery_features <- function(entry, dir, g, table, alpha) {
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
  out$effect_scale <- tg$effect_scale
  if (identical(g$test, "wald") && length(unlist(g$coefficients)) == 1L && !isFALSE(out$predictor_match) && "effect" %in% names(tf)) {
    m <- rows[match(tested, rows$feature_id), , drop = FALSE]; e <- tf$effect[match(tested, tf$feature_id)]
    ok <- !is.na(m$estimate) & !is.na(e)
    if (sum(ok) >= 3) { out$estimate_cor <- stats::cor(m$estimate[ok], e[ok]); out$sign_agreement <- mean(sign(m$estimate[ok]) == sign(e[ok])) }
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
