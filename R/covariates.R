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

# Build the design spec. mean: one-sided formula; covariates: NULL, a data
# frame or a list of data frames; scale: standardise numeric covariates
# (mean 0, sd 1 over the areas, time steps or cells they are given for).
build_design <- function(mean, covariates, frame, time, scale = TRUE) {
  if (!inherits(mean, "formula") || length(mean) != 2)
    stop("`mean` must be a one-sided formula, e.g. ~ trend() + rain", call. = FALSE)
  tt <- stats::terms(mean, specials = "trend")
  if (attr(tt, "intercept") != 1) stop("the mean must keep its intercept (no `- 1` or `+ 0`)", call. = FALSE)
  labels <- attr(tt, "term.labels")
  trend <- "trend()" %in% labels
  vars <- setdiff(labels, "trend()")
  odd <- vars[make.names(vars) != vars]
  if (length(odd))
    stop(sprintf(paste("terms of `mean` must be trend() or column names of `covariates`, not: %s.",
                       "Compute transforms and interactions as columns first."), paste(odd, collapse = ", ")),
         call. = FALSE)
  m <- frame$m; T <- time$T
  tabs <- covariate_tables(covariates, frame, time)
  cols <- list()
  for (v in vars) {
    home <- which(vapply(tabs, function(tb) v %in% names(tb$values), logical(1)))
    if (!length(home)) stop(sprintf("covariate `%s` is not a column of `covariates`", v), call. = FALSE)
    if (length(home) > 1) stop(sprintf("covariate `%s` is in more than one covariate table", v), call. = FALSE)
    tb <- tabs[[home]]
    x <- tb$values[[v]]
    if (anyNA(x)) stop(sprintf("covariate `%s` has missing values", v), call. = FALSE)
    if (is.numeric(x)) {
      ctr <- 0; scl <- 1
      if (scale) {
        ctr <- mean(x); scl <- stats::sd(x)
        if (!is.finite(scl) || scl == 0) stop(sprintf("covariate `%s` is constant", v), call. = FALSE)
      }
      cols[[v]] <- list(name = v, kind = tb$kind, values = shape_values((x - ctr) / scl, tb$kind, m, T),
                        centre = ctr, scale = scl, source = v)
    } else {
      f <- factor(x)
      if (nlevels(f) < 2) stop(sprintf("covariate `%s` has only one level", v), call. = FALSE)
      for (lv in levels(f)[-1]) {                         # treatment contrasts: first level is the baseline
        nm <- paste0(v, lv)
        cols[[nm]] <- list(name = nm, kind = tb$kind, values = shape_values(as.numeric(f == lv), tb$kind, m, T),
                           centre = 0, scale = 1, source = v)
      }
    }
  }
  structure(list(formula = mean, trend = trend, cols = unname(cols), m = m, T = T,
                 names = c("intercept", if (trend) "trend", unname(vapply(cols, `[[`, "", "name"))), scaled = scale),
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
  m <- frame$m; T <- time$T
  lapply(seq_along(covariates), function(i) {
    tb <- as.data.frame(covariates[[i]])
    what <- if (length(covariates) > 1) sprintf("covariate table %d", i) else "`covariates`"
    has_id <- "id" %in% names(tb); has_t <- any(c("t", "time") %in% names(tb))
    if (!has_id && !has_t) stop(sprintf("%s needs an `id` column, a `t` or `time` column, or both", what),
                                call. = FALSE)
    a <- if (has_id) match(as.character(tb$id), frame$ids) else NULL
    if (has_id && anyNA(a)) stop(sprintf("%s has ids not in the frame, e.g. %s", what,
                                         tb$id[which(is.na(a))[1]]), call. = FALSE)
    t <- NULL
    if (has_t) {
      t <- if ("t" %in% names(tb)) suppressWarnings(as.integer(tb$t)) else match(as.character(tb$time), time$labels)
      if (anyNA(t) || any(t < 1 | t > T))
        stop(sprintf("%s has time steps outside the time frame (give `t` as 1..%d or `time` as a label)", what, T),
             call. = FALSE)
    }
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

# canonical vector -> what design_eval() indexes: area (length m), time
# (length T), cell (m x T matrix)
shape_values <- function(x, kind, m, T) if (kind == "cell") matrix(x, m, T, byrow = TRUE) else x

# the design rows of cells (area a, time step t): n x p
design_eval <- function(des, a, t) {
  T <- des$T
  X <- matrix(1, length(a), length(des$names), dimnames = list(NULL, des$names))
  j <- 1L
  if (des$trend) { j <- j + 1L; X[, j] <- (t - (T + 1) / 2) / (T / 2) }   # = design_X
  for (cl in des$cols) {
    j <- j + 1L
    X[, j] <- switch(cl$kind, area = cl$values[a], time = cl$values[t], cell = cl$values[cbind(a, t)])
  }
  X
}

#' A linear trend in time
#'
#' For use in the `mean` formula of [fit_areal()] and [simulate_areal()]
#' only: a linear effect of time, scaled to run from -1 at the first time
#' step to 1 at the last.
#' @return Not meant to be called.
#' @export
trend <- function() stop("trend() is only meaningful inside the `mean` formula of fit_areal()", call. = FALSE)

#' @export
print.stareal_design <- function(x, ...) {
  cat(sprintf("<stareal_design> mean %s: %d fixed effects (%s)\n", deparse(x$formula), length(x$names),
              paste(x$names, collapse = ", ")))
  for (cl in x$cols) if (cl$scale != 1 || cl$centre != 0)
    cat(sprintf("  %s (by %s): standardised, centre %.4g, scale %.4g\n", cl$name, cl$kind, cl$centre, cl$scale))
  invisible(x)
}
