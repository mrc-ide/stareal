# loo.R
#
# Approximate leave-one-out cross-validation from one fit, and a table that
# compares fits.
#
# Plain PSIS-LOO on the pointwise likelihood p(y_i | x_i) does not work
# here: every observed cell has its own latent value, so leaving a cell out
# moves the posterior of x_i a long way, and the importance weights have
# infinite variance (crossval.R says so too). The fix (Vehtari, Mononen,
# Tolvanen, Sivula and Winther 2016, "Bayesian leave-one-out
# cross-validation approximations for Gaussian latent variable models",
# JMLR 17(103), section 3.6, "integrated importance sampling") is to
# integrate x_i out first. For each stored draw (theta, x_o),
#   x_i | x_{-i}, theta ~ N(x_i - (P r)_i / P_ii, 1 / P_ii),
#   P = Sigma_oo^{-1}, r = x_o - mu_o,
# (Sigma_oo with the fixed effects and the nugget, as the sampler uses it),
# and p(y_i | x_{-i}, theta) is the binomial likelihood integrated over that
# Gaussian (log_binom_gauss(), crossval.R). These are the importance ratios'
# denominators: PSIS-LOO on the matrix log p(y_i | x_{-i}, theta) is
# integrated importance-sampling LOO, smoothed by the loo package.
#
# The unit left out is an observed cell (all the surveys of one area-time
# cell), as in crossval(). It is the leave-one-out of a cell given every
# other cell, so it measures interpolation, not extrapolation to unsurveyed
# areas: use crossval() with folds by area or region for that.
#
# ------------------------------------------------------------------

#' Approximate leave-one-out cross-validation of a fit
#'
#' Estimates how well the model predicts each observed cell when that cell
#' is left out, from the one fit: no refits. Each observed cell's own
#' latent value is integrated out exactly before Pareto-smoothed importance
#' sampling (PSIS; Vehtari et al. 2016, 2017), which makes the
#' approximation stable where plain PSIS-LOO on this model is not. The
#' result is a `loo` object, so `loo::loo_compare()` and the `loo` package's
#' other tools work on it; [compare_models()] puts several fits side by
#' side.
#'
#' The unit left out is an observed cell (all the surveys of one area and
#' time step), predicted as a new survey of the cell would be, nugget
#' included, as in [crossval()]. Leaving one cell out keeps its neighbours
#' in, so this measures how well the model fills gaps between data; it says
#' little about extrapolation to areas with no data. For that, use
#' [crossval()] with folds by area or region.
#'
#' Cells with a high Pareto \eqn{\hat k} (above 0.7, or less with fewer
#' than about 2,200 draws: the `loo` package's threshold
#' \eqn{\min(1 - 1/\log_{10} S, 0.7)} for \eqn{S} draws) have unreliable
#' estimates; the print-out counts them. A few are normal for cells far from
#' the rest of the data, or with very large samples. If there are many,
#' check against [crossval()], or use more draws.
#'
#' @param x A [fit_areal()] object.
#' @param n_draws Number of stored draws to use, taken evenly from every
#'   chain (default: all).
#' @param verbose Report progress.
#' @param ... Not used.
#' @return A `psis_loo` object from the `loo` package (see [loo::loo()]),
#'   with class `stareal_loo` added and an extra table `cells`: one row per
#'   observed cell, with `cell_id`, `id`, `t`, `time`, `y`, `k`, `elpd_loo`
#'   (its log predictive density when left out) and `pareto_k`.
#' @references Vehtari, A., Mononen, T., Tolvanen, V., Sivula, T. and
#'   Winther, O. (2016). Bayesian leave-one-out cross-validation
#'   approximations for Gaussian latent variable models. *Journal of Machine
#'   Learning Research* 17(103), 1-38.
#'
#'   Vehtari, A., Gelman, A. and Gabry, J. (2017). Practical Bayesian model
#'   evaluation using leave-one-out cross-validation and WAIC. *Statistics
#'   and Computing* 27, 1413-1432.
#' @examples
#' sim <- simulate_areal(areal_lattice(4, 4), time_frame(n_steps = 5), n_cells = 40, seed = 1)
#' fit <- fit_areal(sim, control = list(n_iter = 200, n_burn = 100, n_store = 50, n_chains = 2, n_cores = 1,
#'                                      verbose = FALSE))
#' loo(fit)
#' @exportS3Method loo::loo
loo.stareal_fit <- function(x, n_draws = NULL, verbose = FALSE, ...) {
  check_dots(...)
  fit <- x
  ll <- loo_loglik(fit, n_draws, verbose)
  # (r_eff is the same for the likelihood scaled per cell, which keeps exp() in range)
  r_eff <- suppressWarnings(loo::relative_eff(exp(ll$ll - matrix(apply(ll$ll, 2, max), nrow(ll$ll), ncol(ll$ll),
                                                                 byrow = TRUE)), chain_id = ll$chain))
  out <- suppressWarnings(loo::loo(ll$ll, r_eff = r_eff))
  cells <- fit$data$cells
  out$cells <- data.frame(cells[c("cell_id", "id", "t", "time", "y", "k")],
                          elpd_loo = out$pointwise[, "elpd_loo"], pareto_k = out$diagnostics$pareto_k,
                          row.names = NULL)
  out$version <- version_info()
  class(out) <- c("stareal_loo", class(out))
  out
}

#' @export
print.stareal_loo <- function(x, ...) {
  e <- x$estimates
  thr <- k_threshold(x)
  n_high <- sum(x$diagnostics$pareto_k > thr)
  cat(sprintf("<stareal_loo> integrated PSIS-LOO over %d observed cells, %d draws\n", nrow(x$pointwise),
              attr(x, "dims")[1]))
  cat(sprintf("  elpd_loo %.1f (se %.1f); p_loo %.1f; %d cell%s with Pareto k > %.2f\n", e["elpd_loo", 1],
              e["elpd_loo", 2], e["p_loo", 1], n_high, if (n_high == 1) "" else "s", thr))
  invisible(x)
}

# The Pareto k above which a cell's estimate is unreliable: 0.7, or lower
# with fewer draws, min(1 - 1 / log10(S), 0.7) (Vehtari et al. 2024, the
# loo package's rule)
k_threshold <- function(l) min(1 - 1 / log10(attr(l, "dims")[1]), 0.7)

# The matrix log p(y_i | x_{-i}, theta): stored draws (rows, chain by
# chain) x observed cells, and the chain of each row.
loo_loglik <- function(fit, n_draws = NULL, verbose = FALSE) {
  st <- fit_setup_of(fit)
  obs <- st$obs; eig <- st$eig; cells <- fit$data$cells
  S <- dim(fit$draws$x)[1]; C <- dim(fit$draws$x)[2]
  keep <- if (is.null(n_draws)) seq_len(S) else {
    n_draws <- check_whole(n_draws, "n_draws", 2)
    if (n_draws > S * C) stop(sprintf("`n_draws` must be at most the %d stored draws", S * C), call. = FALSE)
    unique(round(seq(1, S, length.out = ceiling(n_draws / C))))      # the same draws from every chain
  }
  rows <- as.vector(outer(keep, (seq_len(C) - 1) * S, `+`))           # draw fastest, chain by chain
  pool <- function(a) matrix(a, S * C, dim(a)[3])[rows, , drop = FALSE]
  th <- pool(fit$draws$theta); xo <- pool(fit$draws$x)
  D <- nrow(xo); n <- ncol(xo)
  M <- V <- matrix(NA_real_, D, n)
  every <- max(1L, D %/% 10L); t0 <- proc.time()[["elapsed"]]
  for (d in seq_len(D)) {
    mo <- loo_moments(th[d, ], xo[d, ], obs, eig, st$B, st$mu_o)
    M[d, ] <- mo$mean; V[d, ] <- mo$var
    if (verbose && d %% every == 0)
      say(TRUE, "loo: %d of %d draws (%s left)", d, D,
          format_seconds((proc.time()[["elapsed"]] - t0) / d * (D - d)))
  }
  list(ll = log_binom_gauss(cells$y, cells$k, M, V), chain = rep(seq_len(C), each = length(keep)))
}

# For one draw: the mean and variance of each observed cell's latent value
# given all the others, x_i | x_{-i}, theta (fixed effects integrated out,
# nugget included)
loo_moments <- function(theta, x_o, obs, eig, B, mu_o, jitter = 1e-8) {
  sig <- build_sigma_oo(theta, obs, eig, B, jitter = jitter)
  P <- chol2inv(t(sig$L))                                   # Sigma_oo^{-1}
  pd <- P[seq.int(1L, length(P), by = nrow(P) + 1L)]
  list(mean = x_o - as.vector(P %*% (x_o - mu_o)) / pd, var = 1 / pd)
}

#' Compare models fitted to the same data
#'
#' One row per fit: its model, whether it passes [diagnose()], its
#' approximate leave-one-out score from [loo()][loo.stareal_fit], and, if
#' given, its [crossval()] scores. All fits must be made to the same
#' observed cells (the same data, whatever the formula, spatial prior or
#' nugget).
#'
#' Scores are expected log predictive densities: higher is better. `elpd_diff`
#' is each model's difference from the best (0 for the best, negative for
#' the others), summed over cells, and `se_diff` its standard error from the
#' cell-by-cell differences, as in `loo::loo_compare()`. A difference
#' smaller than about two standard errors is not a clear preference. The
#' cross-validation columns are per cell (`cv_log_score` is the mean log
#' score of [compare_cv()]), as there.
#'
#' Why not AIC, BIC or DIC: AIC needs a number of parameters, which a model
#' with a latent value per cell does not have (the effective number,
#' `p_loo` here, depends on the hyperparameters); BIC approximates a
#' marginal likelihood that, in a hierarchical model, is driven by the
#' priors; and DIC is unreliable for latent Gaussian models. Predictive
#' scores avoid all three problems.
#'
#' @param ... Two or more [fit_areal()] objects, named (e.g.
#'   `compare_models(leroux = fit1, bym2 = fit2)`); unnamed ones are named by
#'   their expressions.
#' @param cv Optional: [crossval()] results, a list named as the fits, with
#'   the same folds and seed for every model (see [compare_cv()]).
#' @param loo Optional: [loo()][loo.stareal_fit] results already computed, a
#'   list named as the fits (any fit left out is computed here).
#' @return A data frame, best first (by `elpd_loo`): `model`, `spatial` (the
#'   spatial prior and any fixed hyperparameters), `nugget`, `n_hyper` (free
#'   hyperparameters), `n_beta` (fixed effects), `pass` ([diagnose()]),
#'   `elpd_loo`, `se_elpd_loo`, `elpd_diff`, `se_diff`, `p_loo`, `n_k_high`
#'   (cells with a high Pareto k, see [loo()][loo.stareal_fit]); and with `cv`: `cv_log_score`,
#'   `cv_diff`, `cv_se_diff`, `cv_pct_outside_95`, `cv_rmse`. The `loo`
#'   results are attached as attribute `loo`.
#' @examples
#' \donttest{
#' sim <- simulate_areal(areal_lattice(4, 4), time_frame(n_steps = 5), n_cells = 40, seed = 1)
#' ctl <- list(n_iter = 200, n_burn = 100, n_store = 50, n_chains = 2, n_cores = 1, verbose = FALSE)
#' compare_models(leroux = fit_areal(sim, control = ctl),
#'                bym2 = fit_areal(sim, spatial = "bym2", control = ctl),
#'                no_nugget = fit_areal(sim, nugget = FALSE, control = ctl))
#' }
#' @export
compare_models <- function(..., cv = NULL, loo = NULL) {
  fits <- list(...)
  nm <- names(fits)
  exprs <- vapply(as.list(substitute(list(...)))[-1], function(e) paste(deparse(e), collapse = ""), "")
  if (is.null(nm)) nm <- exprs
  nm[nm == ""] <- exprs[nm == ""]
  if (length(fits) < 2) stop("give at least two fits", call. = FALSE)
  if (anyDuplicated(nm)) stop("name every fit differently", call. = FALSE)
  if (!all(vapply(fits, inherits, logical(1), "stareal_fit"))) stop("every fit must be made by fit_areal()", call. = FALSE)
  names(fits) <- nm
  key <- function(f) f$data$cells[c("cell_id", "y", "k")]
  if (!all(vapply(fits, function(f) identical(key(f), key(fits[[1]])), logical(1))))
    stop("the fits were not made to the same observed cells (the same data)", call. = FALSE)
  for (a in c("cv", "loo")) {
    v <- get(a)
    if (!is.null(v) && (!is.list(v) || is.null(names(v)) || length(setdiff(names(v), nm))))
      stop(sprintf("`%s` must be a list named as the fits (%s)", a, paste(nm, collapse = ", ")), call. = FALSE)
  }
  if (!is.null(cv) && !setequal(names(cv), nm)) stop("`cv` needs a crossval() result for every fit", call. = FALSE)
  loos <- lapply(nm, function(n) if (!is.null(loo[[n]])) loo[[n]] else loo.stareal_fit(fits[[n]]))
  names(loos) <- nm
  pw <- vapply(loos, function(l) l$pointwise[, "elpd_loo"], numeric(nrow(fits[[1]]$data$cells)))
  if (is.null(dim(pw))) pw <- matrix(pw, nrow = 1)
  best <- which.max(colSums(pw))
  dd <- pw - pw[, best]
  out <- data.frame(model = nm,
                    spatial = vapply(fits, function(f) model_label(fit_model(f)), ""),
                    nugget = vapply(fits, function(f) f$nugget, logical(1)),
                    n_hyper = vapply(fits, function(f) length(fit_theta_names(f)), 1L),
                    n_beta = vapply(fits, function(f) length(f$design$names), 1L),
                    pass = vapply(fits, function(f) isTRUE(suppressWarnings(diagnose(f))$pass), logical(1)),
                    elpd_loo = colSums(pw),
                    se_elpd_loo = vapply(loos, function(l) l$estimates["elpd_loo", 2], 1),
                    elpd_diff = colSums(dd), se_diff = sqrt(nrow(dd)) * apply(dd, 2, stats::sd),
                    p_loo = vapply(loos, function(l) l$estimates["p_loo", 1], 1),
                    n_k_high = vapply(loos, function(l) sum(l$diagnostics$pareto_k > k_threshold(l)), 1L),
                    stringsAsFactors = FALSE, row.names = NULL)
  if (!is.null(cv)) {
    cc <- do.call(compare_cv, cv[nm])
    i <- match(nm, cc$model)
    out$cv_log_score <- cc$mean_log_score[i]; out$cv_diff <- cc$diff[i]; out$cv_se_diff <- cc$se_diff[i]
    out$cv_pct_outside_95 <- cc$pct_outside_95[i]; out$cv_rmse <- cc$rmse[i]
  }
  out <- out[order(-out$elpd_loo), ]
  rownames(out) <- NULL
  attr(out, "loo") <- loos
  out
}
