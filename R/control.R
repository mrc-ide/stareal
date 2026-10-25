# control.R
#
# The user's settings for a fit: the priors (areal_priors) and the sampler
# (mcmc_control). Both return plain, validated lists. fit_areal() also
# accepts a partial list, which it completes with these defaults
# (complete_settings), so priors = list(sigma_sd = 1) works.
#
# ------------------------------------------------------------------

#' Priors for an areal model
#'
#' Priors are on the natural scale of each parameter. The field has
#' precision \eqn{\tau = 1/\sigma^2}; the fixed effects are integrated out
#' analytically, so their prior must be Gaussian.
#'
#' `fit_areal()` also accepts a plain list of any of these arguments, e.g.
#' `priors = list(sigma_sd = 1)`; the others take their defaults.
#'
#' @param sigma_sd Scale of the half-normal prior on `sigma`, the scale of
#'   the space-time field: \eqn{\sigma \sim N^+(0, \code{sigma\_sd}^2)}.
#' @param rho_s Shape parameters of the Beta prior on `rho_s`, the spatial
#'   dependence of the Leroux prior (0: independent areas, 1: intrinsic
#'   CAR). The default Beta(2, 1) favours strong spatial dependence.
#' @param rho_t Shape parameters of the Beta prior on `rho_t`, the lag-one
#'   correlation of the AR1 in time. Default uniform. With one time step
#'   the data say nothing about `rho_t`, and its posterior is this prior.
#' @param sigma_n_sd Scale of the half-normal prior on `sigma_n`, the
#'   standard deviation of the cell-level nugget.
#' @param beta_mean,beta_sd Means and standard deviations of the
#'   independent normal priors on the fixed effects (logit scale). Give one
#'   value (used for every fixed effect), one value per fixed effect in the
#'   order of the design, or a vector named by design column, e.g.
#'   `c(intercept = 0, rain = 0)`, where unnamed columns take the default.
#'   The design columns are the intercept, `trend` (if `trend()` is in the
#'   formula), then the covariates in the order of the formula, a factor
#'   giving one column per level after the first (named by the covariate
#'   and the level, e.g. `zonehigh`). They are listed in
#'   `fit$design$names`.
#' @return A list of the settings.
#' @examples
#' areal_priors()
#' areal_priors(sigma_sd = 1, rho_t = c(4, 1))
#' @export
areal_priors <- function(sigma_sd = 2, rho_s = c(2, 1), rho_t = c(1, 1), sigma_n_sd = 1,
                         beta_mean = 0, beta_sd = 3) {
  pos <- function(x, nm, len = 1) {
    if (!is.numeric(x) || length(x) != len || any(!is.finite(x)) || any(x <= 0))
      stop(sprintf("`%s` must be %s", nm, if (len == 1) "one positive, finite number" else "two positive, finite numbers"),
           call. = FALSE)
  }
  pos(sigma_sd, "sigma_sd"); pos(sigma_n_sd, "sigma_n_sd")
  pos(rho_s, "rho_s", 2); pos(rho_t, "rho_t", 2)
  if (!is.numeric(beta_mean) || !length(beta_mean) || any(!is.finite(beta_mean)))
    stop("`beta_mean` must be finite numbers", call. = FALSE)
  if (!is.numeric(beta_sd) || !length(beta_sd) || any(!is.finite(beta_sd)) || any(beta_sd <= 0))
    stop("`beta_sd` must be positive, finite numbers (for a vague prior use a large value, e.g. 10)", call. = FALSE)
  for (v in c("beta_mean", "beta_sd")) {
    nm <- names(get(v))
    if (!is.null(nm) && (anyNA(nm) || any(nm == "") || anyDuplicated(nm)))
      stop(sprintf("`%s`: name every value, once (or name none)", v), call. = FALSE)
  }
  list(sigma_sd = sigma_sd, rho_s = rho_s, rho_t = rho_t, sigma_n_sd = sigma_n_sd,
       beta_mean = beta_mean, beta_sd = beta_sd)
}

#' Settings for the MCMC sampler
#'
#' Each iteration makes `ess_sweeps` updates of the latent values of the
#' observed cells (elliptical slice sampling) and `theta_moves` pairs of
#' Metropolis moves on the hyperparameters. During burn-in, the proposals
#' adapt and the Gaussian approximation that steers the latent updates is
#' refitted to the chain (twice, at half and three-quarters of burn-in);
#' both are fixed afterwards, so the draws after burn-in come from an exact
#' MCMC. Chains start from dispersed hyperparameters and latent values, so
#' that R-hat in [diagnose()] can reveal chains that have not mixed.
#'
#' The defaults are a pilot run. A run for publication is typically much
#' longer (e.g. 40,000 iterations with 8,000 burn-in and `theta_thin = 10`).
#' Use [diagnose()] to judge whether a run is long enough.
#'
#' `fit_areal()` also accepts a plain list of any of these arguments, e.g.
#' `control = list(n_iter = 500)`; the others take their defaults.
#'
#' @param n_iter Iterations per chain, burn-in included.
#' @param n_burn Burn-in iterations (adaptation; discarded).
#' @param n_store Posterior draws kept per chain, spread evenly over the
#'   iterations after burn-in (the last iteration is always kept). These
#'   are the draws used for prediction.
#' @param n_chains Number of chains.
#' @param ess_sweeps Latent updates per iteration.
#' @param theta_moves Hyperparameter move pairs per iteration.
#' @param theta_thin Keep every `theta_thin`-th iteration of the
#'   hyperparameter trace (burn-in included). At least 4 traced iterations
#'   must fall after burn-in.
#' @param seed Random seed, any whole number >= 0. Each chain gets its own
#'   seed, drawn from this one, and R's default generator is used whatever
#'   the session's `RNGkind()`. The caller's random-number stream is left
#'   untouched.
#' @param n_cores Number of worker processes; chains run one per worker.
#'   Default: one per chain, at most the number of cores less one. Use 1 to
#'   run the chains one after another in this R session.
#' @param blas_threads BLAS threads per worker. One is safest; with an
#'   optimised BLAS (MKL, OpenBLAS, Accelerate) and large fits, more may be
#'   faster, as long as `n_cores * blas_threads` does not exceed the cores.
#' @param checkpoint_every Iterations between checkpoints, when
#'   `fit_areal()` is given a `checkpoint_dir`.
#' @param verbose Print progress: a line per chain after 10 iterations,
#'   then every tenth of the run (at most every 500 iterations), with the
#'   time left.
#' @return A list of the settings.
#' @examples
#' mcmc_control(n_iter = 500, n_burn = 200, n_store = 100, n_chains = 2)
#' @export
mcmc_control <- function(n_iter = 3000, n_burn = 1000, n_store = 250, n_chains = 4, ess_sweeps = 10,
                         theta_moves = 2, theta_thin = 1, seed = 1, n_cores = NULL, blas_threads = 1,
                         checkpoint_every = 500, verbose = TRUE) {
  n_iter <- check_whole(n_iter, "n_iter"); n_burn <- check_whole(n_burn, "n_burn", 0)
  n_store <- check_whole(n_store, "n_store"); n_chains <- check_whole(n_chains, "n_chains")
  ess_sweeps <- check_whole(ess_sweeps, "ess_sweeps"); theta_moves <- check_whole(theta_moves, "theta_moves")
  theta_thin <- check_whole(theta_thin, "theta_thin"); seed <- check_whole(seed, "seed", 0)
  blas_threads <- check_whole(blas_threads, "blas_threads")
  checkpoint_every <- check_whole(checkpoint_every, "checkpoint_every")
  n_cores <- check_whole(n_cores, "n_cores", null_ok = TRUE)
  if (n_burn >= n_iter) stop("`n_burn` must be less than `n_iter`", call. = FALSE)
  if (n_store > n_iter - n_burn)
    stop(sprintf("`n_store` (%d) cannot exceed the %d iterations after burn-in", n_store, n_iter - n_burn),
         call. = FALSE)
  n_traced <- n_iter %/% theta_thin - n_burn %/% theta_thin
  if (n_traced < 4)
    stop(sprintf("`theta_thin` = %d leaves %d traced iterations after burn-in; at least 4 are needed (use theta_thin <= %d)",
                 theta_thin, n_traced, max(1L, (n_iter - n_burn) %/% 4L)), call. = FALSE)
  if (is.null(n_cores)) n_cores <- default_cores(n_chains)
  if (!is.logical(verbose) || length(verbose) != 1 || is.na(verbose)) stop("`verbose` must be TRUE or FALSE", call. = FALSE)
  list(n_iter = n_iter, n_burn = n_burn, n_store = n_store, n_chains = n_chains, ess_sweeps = ess_sweeps,
       theta_moves = theta_moves, theta_thin = theta_thin, seed = seed, n_cores = n_cores,
       blas_threads = blas_threads, checkpoint_every = checkpoint_every, verbose = verbose)
}

# A settings list completed by its constructor: NULL gives the defaults; a
# partial list is passed to the constructor, so unknown names are errors and
# every value is validated.
complete_settings <- function(x, constructor, what) {
  if (is.null(x)) return(constructor())
  if (!is.list(x)) stop(sprintf("`%s` must be a list (see ?%s)", what, deparse(substitute(constructor))), call. = FALSE)
  bad <- setdiff(names(x), names(formals(constructor)))
  if (is.null(names(x)) || any(names(x) == "") || length(bad))
    stop(sprintf("`%s`: unknown setting(s) %s; valid: %s", what,
                 if (length(bad)) paste(bad, collapse = ", ") else "(unnamed)",
                 paste(names(formals(constructor)), collapse = ", ")), call. = FALSE)
  do.call(constructor, x)
}

# the iterations at which draws are stored: evenly spread after burn-in,
# the last iteration always included
store_iterations <- function(n_iter, n_burn, n_store) {
  n_burn + round(seq_len(n_store) * (n_iter - n_burn) / n_store)
}
