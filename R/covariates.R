# covariates.R
#
# The mean of the model: an intercept, optionally a linear trend in time,
# and covariates that vary by area, by time step, or by both. The user
# gives a one-sided formula naming columns of covariate tables keyed by
# area id and/or time step; build_design() turns them into a compact
# design spec, and design_eval() gives the design rows of any cells. The
# spec is the only source of the fixed-effect design: fit_setup() uses it
# for the observed cells, predict() for the full table (design_rows) and
# crossval() for held-out cells (design_cells).
#
# Covariates are stored once per area, per time step or per area-time cell,
# never expanded to the full table, so memory stays at the size the user
# gave.
#
# ------------------------------------------------------------------

# names a design column may not take: the model's own parameters, and the
# names of the latent cells in the draws
reserved_names <- c("intercept", "trend", "sigma", "rho_s", "rho_t", "sigma_n")

# Build the design spec. formula: one-sided formula; covariates: NULL, a
# data frame or a list of data frames; scale: standardise numeric
# covariates (mean 0, sd 1 over the areas, time steps or cells they are
# given for).
build_design <- function(formula, covariates, frame, time, scale = TRUE) {
  if (!inherits(formula, "formula") || length(formula) != 2)
    stop("`formula` must be a one-sided formula, e.g. ~ trend() + rain", call. = FALSE)
  tt <- stats::terms(formula, specials = "trend")
  if (attr(tt, "intercept") != 1) stop("the formula must keep its intercept (no `- 1` or `+ 0`)", call. = FALSE)
  if (!is.null(attr(tt, "offset")))
    stop("offset() is not supported in the formula; give the offset as a covariate column, or leave it out",
         call. = FALSE)
  labels <- attr(tt, "term.labels")
  trend <- "trend()" %in% labels
  vars <- setdiff(labels, "trend()")
  if ("trend" %in% vars) stop("the time trend is written trend(), with parentheses", call. = FALSE)
  odd <- vars[make.names(vars) != vars]
  if (length(odd))
    stop(sprintf(paste("terms of the formula must be trend() or column names of `covariates`, not: %s.",
                       "Compute transforms and interactions as columns first."), paste(odd, collapse = ", ")),
         call. = FALSE)
  m <- frame$m; T <- time$n_steps
  if (trend && T == 1) {
    message("one time step: trend() dropped from the formula")
    trend <- FALSE
  }
  tabs <- covariate_tables(covariates, frame, time)
  cols <- list()
  add <- function(nm, kind, values, centre, scale, source) {
    if (nm %in% c(reserved_names, names(cols)) || grepl("^(logit|prev)\\[", nm))
      stop(sprintf("design column `%s` (from covariate `%s`) clashes with %s; rename the covariate", nm, source,
                   if (nm %in% names(cols)) "another covariate" else "a model parameter"), call. = FALSE)
    cols[[nm]] <<- list(name = nm, kind = kind, values = shape_values(values, kind, m, T),
                        centre = centre, scale = scale, source = source)
  }
  for (v in vars) {
    home <- which(vapply(tabs, function(tb) v %in% names(tb$values), logical(1)))
    if (!length(home)) stop(sprintf("covariate `%s` is not a column of `covariates`", v), call. = FALSE)
    if (length(home) > 1) stop(sprintf("covariate `%s` is in more than one covariate table", v), call. = FALSE)
    tb <- tabs[[home]]
    x <- tb$values[[v]]
    if (is.numeric(x)) {
      if (any(!is.finite(x))) stop(sprintf("covariate `%s` has missing or non-finite values", v), call. = FALSE)
      ctr <- mean(x); scl <- stats::sd(x)
      if (!is.finite(scl) || scl <= sqrt(.Machine$double.eps) * max(1, abs(ctr)))
        stop(sprintf("covariate `%s` is constant", v), call. = FALSE)
      if (!scale) { ctr <- 0; scl <- 1 }
      add(v, tb$kind, (x - ctr) / scl, ctr, scl, v)
    } else {
      if (anyNA(x)) stop(sprintf("covariate `%s` has missing values", v), call. = FALSE)
      f <- factor(x)
      if (nlevels(f) < 2) stop(sprintf("covariate `%s` has only one level", v), call. = FALSE)
      for (lv in levels(f)[-1])                           # treatment contrasts: first level is the baseline
        add(paste0(v, lv), tb$kind, as.numeric(f == lv), 0, 1, v)
    }
  }
  environment(formula) <- emptyenv()                      # do not carry the caller's frame (and data) along
  structure(list(formula = formula, trend = trend, cols = unname(cols), m = m, T = T,
                 names = c("intercept", if (trend) "trend", unname(vapply(cols, `[[`, "", "name"))),
                 scaled = scale),
            class = "stareal_design")
}

# the covariate tables as list(kind, values: data frame of covariate
# columns in canonical order) with kind "area" (frame order), "time" (steps
# 1..T) or "cell" (area-major: (a - 1) * T + t)
covariate_tables <- function(covariates, frame, time) {
  if (is.null(covariates)) return(list())
  if (is.data.frame(covariates)) covariates <- list(covariates)
  if (!is.list(covariates) || !all(vapply(covariates, is.data.frame, logical(1))))
    stop("`covariates` must be a data frame or a list of data frames", call. = FALSE)
  m <- frame$m; T <- time$n_steps
  lapply(seq_along(covariates), function(i) {
    tb <- as.data.frame(covariates[[i]])
    what <- if (length(covariates) > 1) sprintf("covariate table %d", i) else "`covariates`"
    has_id <- "id" %in% names(tb); has_t <- any(c("t", "time") %in% names(tb))
    if (!has_id && !has_t) stop(sprintf("%s needs an `id` column, a `t` or `time` column, or both", what),
                                call. = FALSE)
    a <- if (has_id) match(id_chr(tb$id), frame$ids) else NULL
    if (has_id && anyNA(a)) stop(sprintf("%s has ids not in the frame, e.g. %s", what,
                                         tb$id[which(is.na(a))[1]]), call. = FALSE)
    t <- NULL
    if (has_t) t <- if ("t" %in% names(tb)) time_steps(tb$t, time, what, "t") else time_steps(tb$time, time, what, "time")
    kind <- if (has_id && has_t) "cell" else if (has_id) "area" else "time"
    key <- switch(kind, area = a, time = t, cell = (a - 1L) * T + t)
    n_key <- switch(kind, area = m, time = T, cell = m * T)
    if (anyDuplicated(key)) stop(sprintf("%s has more than one row for some %s", what,
                                         switch(kind, area = "areas", time = "time steps", cell = "area-time cells")),
                                 call. = FALSE)
    if (length(key) < n_key)
      stop(sprintf("%s covers %d of the %d %s; covariates must cover every one", what, length(key), n_key,
                   switch(kind, area = "areas", time = "time steps", cell = "area-time cells")), call. = FALSE)
    vals <- tb[order(key), setdiff(names(tb), c("id", "t", "time")), drop = FALSE]
    rownames(vals) <- NULL
    list(kind = kind, values = vals)
  })
}

# Time steps from a key column. key = "time": the time frame's labels
# (numbers are matched as text, so a year 2005 matches the label "2005").
# key = "t": step numbers, whole, 1..T (text that matches a label is taken
# as that label). na_ok: return NA for values that do not match, instead of
# stopping with an error naming `what`.
time_steps <- function(v, time, what, key = c("t", "time"), na_ok = FALSE) {
  key <- match.arg(key)
  T <- time$n_steps
  if (key == "time") {
    t <- match(trimws(as.character(v)), time$labels)
  } else {
    num <- as_num(v)
    t <- ifelse(!is.na(num) & num == round(num) & num >= 1 & num <= T, num, NA)
    if (is.character(v) || is.factor(v)) t <- ifelse(is.na(t), match(as.character(v), time$labels), t)
  }
  if (!na_ok && anyNA(t)) {
    bad <- as.character(v[which(is.na(t))[1]])
    stop(sprintf("%s has time %s not in the time frame, e.g. %s (%s)", what,
                 if (key == "time") "labels" else "steps", bad,
                 if (key == "time") sprintf("labels run %s to %s", time$labels[1], time$labels[T]) else
                   sprintf("`t` takes step numbers 1..%d; give labels such as \"%s\" in a `time` column", T,
                           time$labels[1])), call. = FALSE)
  }
  as.integer(t)
}

# canonical vector -> what design_eval() indexes: area (length m), time
# (length T), cell (m x T matrix)
shape_values <- function(x, kind, m, T) if (kind == "cell") matrix(x, m, T, byrow = TRUE) else x

# the design rows of cells (area a, time step t): n x p
design_eval <- function(des, a, t) {
  X <- matrix(1, length(a), length(des$names), dimnames = list(NULL, des$names))
  j <- 1L
  if (des$trend) { j <- j + 1L; X[, j] <- trend_value(t, des$T) }
  for (cl in des$cols) {
    j <- j + 1L
    X[, j] <- switch(cl$kind, area = cl$values[a], time = cl$values[t], cell = cl$values[cbind(a, t)])
  }
  X
}

# The largest prior variance the fixed effects give any cell,
# max_cells x' B x, and the column that contributes most. Above ~1e8 the
# field (variance ~1) is lost in rounding once X B X' is added to the
# observed-cell covariance.
design_prior_var <- function(des, beta_sd) {
  colmax <- c(1, if (des$trend) 1, vapply(des$cols, function(cl) max(abs(cl$values)), 1))
  contrib <- (colmax * beta_sd)^2
  list(var = sum(contrib), worst = des$names[which.max(contrib)], sd = sqrt(max(contrib)))
}

#' A linear trend in time
#'
#' For use in the `formula` of [fit_areal()] and [simulate_areal()] only: a
#' linear effect of time, running from -1 at the first time step to 1 at
#' the last, so its coefficient is half the change in logit prevalence
#' over the window. With one time step it is dropped.
#' @return Not meant to be called.
#' @export
trend <- function() stop("trend() is only meaningful inside the formula of fit_areal()", call. = FALSE)

#' @export
print.stareal_design <- function(x, ...) {
  cat(sprintf("<stareal_design> %s: %d fixed effects (%s)\n",
              paste(deparse(x$formula, width.cutoff = 500L), collapse = " "), length(x$names),
              paste(x$names, collapse = ", ")))
  for (cl in x$cols) if (cl$scale != 1 || cl$centre != 0)
    cat(sprintf("  %s (by %s): standardised, centre %.4g, scale %.4g\n", cl$name, cl$kind, cl$centre, cl$scale))
  invisible(x)
}
