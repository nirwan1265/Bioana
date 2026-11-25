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
  assumed <- any(vapply(design$spec$assumptions %||% list(), function(x) {
    f <- if (is.list(x)) x$field %||% "" else as.character(x)
    grepl("analysis\\.(primary_goal|goals)", f)
  }, logical(1)))
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
  list(id = id, source = src, assumed = assumed && src == "study_spec", terms = unlist(a$goal_terms[[id]]),
       listed = listed, note = note, unknown_to_biodes = unknown_to_biodes)
}

# Design copy for the goal: trend fits its term as numeric.
.ba_goal_design <- function(design, goal) {
  if (identical(goal$id, "trend")) {
    tv <- goal$terms[1]
    if (is.null(tv) || is.na(tv)) stop("Goal 'trend' needs analysis.goal_terms$trend: the numeric term.", call. = FALSE)
    if (!tv %in% design$parsed$vars) stop(sprintf("Goal 'trend': term '%s' is not in analysis.model (%s).", tv, design$model), call. = FALSE)
    design$parsed$numeric <- unique(c(design$parsed$numeric, tv))
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
    return(list(test = "wald", coefs = treat, drop = character(), drop_terms = character(), label = NULL))
  }
  if (id == "interaction") {
    tt <- goal$terms
    if (length(tt) != 2) stop("Goal 'interaction' needs analysis.goal_terms$interaction: two model terms, e.g. [\"treatment\", \"time\"].", call. = FALSE)
    terms_hit <- unique(ct[vapply(ct, function(t) setequal(vars_of(t), tt) && grepl(":", t), logical(1))])
    if (!length(terms_hit)) stop(sprintf("Goal 'interaction': analysis.model (%s) has no %s:%s term.", deparse(p$formula), tt[1], tt[2]), call. = FALSE)
    k <- cols[ct %in% terms_hit]
    lab <- paste(tt, collapse = ":")
    if (length(k) == 1L) return(list(test = "wald", coefs = k, drop = character(), drop_terms = character(), label = NULL))
    return(list(test = "lrt", coefs = character(), drop = k, drop_terms = terms_hit, label = sprintf("interaction %s (%d coefficients)", lab, length(k))))
  }
  if (id == "time_course") {
    tv <- goal$terms[1] %||% "time"
    if (is.na(tv)) tv <- "time"
    terms_hit <- unique(ct[vapply(ct, function(t) tv %in% vars_of(t), logical(1))])
    if (!length(terms_hit)) stop(sprintf("Goal 'time_course': analysis.model has no term with '%s'.", tv), call. = FALSE)
    k <- cols[ct %in% terms_hit]
    return(list(test = "lrt", coefs = character(), drop = k, drop_terms = terms_hit,
                label = sprintf("time_course: %s (%d coefficients)", paste(terms_hit, collapse = ", "), length(k))))
  }
  if (id == "trend") {
    tv <- goal$terms[1]
    k <- cols[ct == tv]
    if (!length(k)) stop(sprintf("Goal 'trend': analysis.model has no main effect of '%s'.", tv), call. = FALSE)
    return(list(test = "wald", coefs = k, drop = character(), drop_terms = character(), label = NULL))
  }
  stop(sprintf("Goal '%s' has no test plan in Bioana yet.", id), call. = FALSE)
}

# Summary block for the goal of one data set.
.ba_goal_summary <- function(goal, plan) {
  list(id = goal$id, source = goal$source, assumed_by_biodes = isTRUE(goal$assumed), test = plan$test,
       coefficients = if (plan$test == "wald") as.list(plan$coefs) else as.list(plan$drop),
       dropped_terms = if (plan$test == "lrt") as.list(plan$drop_terms), label = plan$label,
       note = goal$note, unknown_to_biodes_registry = if (isTRUE(goal$unknown_to_biodes)) TRUE)
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
