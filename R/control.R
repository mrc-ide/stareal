# control.R
#
# The user's settings for a fit: the priors (areal_priors) and the sampler
# (mcmc_control). Defaults are those of Africa_areal: the priors of
# config/settings.yml (model:) and the pilot profile of config/profiles.yml.
#
# ------------------------------------------------------------------

#' Priors for an areal model
#'
#' Priors are on the natural scale of each parameter. The field has
#' precision \eqn{\tau = 1/\sigma^2}; the fixed effects are integrated out
#' analytically, so their prior must be Gaussian.
#'
#' @param sigma_sd Scale of the half-normal prior on `sigma`, the scale of
#'   the space-time field: \eqn{\sigma \sim N^+(0, \code{sigma\_sd}^2)}.
#' @param rho_s Shape parameters of the Beta prior on `rho_s`, the spatial
#'   dependence of the Leroux prior (0: independent areas, 1: intrinsic
#'   CAR). The default Beta(2, 1) favours strong spatial dependence.
#' @param rho_t Shape parameters of the Beta prior on `rho_t`, the lag-one
#'   correlation of the AR1 in time. Default uniform.
#' @param sigma_n_sd Scale of the half-normal prior on `sigma_n`, the
#'   standard deviation of the cell-level nugget.
#' @param beta_mean,beta_sd Means and standard deviations of the
#'   independent normal priors on the fixed effects (logit scale), in the
#'   order of the design: intercept, then trend (scaled to run from -1 to 1
#'   over the time window). A single value is used for every fixed effect.
#' @return An object of class `stareal_priors`.
#' @examples
#' areal_priors()
#' areal_priors(sigma_sd = 1, rho_t = c(4, 1))
#' @export
areal_priors <- function(sigma_sd = 2, rho_s = c(2, 1), rho_t = c(1, 1), sigma_n_sd = 1,
                         beta_mean = 0, beta_sd = 3) {
  pos <- function(x, nm, len = 1) {
    if (!is.numeric(x) || length(x) != len || anyNA(x) || any(x <= 0))
      stop(sprintf("`%s` must be %s", nm, if (len == 1) "one positive number" else "two positive numbers"),
           call. = FALSE)
  }
  pos(sigma_sd, "sigma_sd"); pos(sigma_n_sd, "sigma_n_sd")
  pos(rho_s, "rho_s", 2); pos(rho_t, "rho_t", 2)
  if (!is.numeric(beta_mean) || !length(beta_mean) || anyNA(beta_mean))
    stop("`beta_mean` must be numbers", call. = FALSE)
  if (!is.numeric(beta_sd) || !length(beta_sd) || anyNA(beta_sd) || any(beta_sd <= 0))
    stop("`beta_sd` must be positive numbers", call. = FALSE)
  structure(list(sigma_sd = sigma_sd, rho_s = rho_s, rho_t = rho_t, sigma_n_sd = sigma_n_sd,
                 beta_mean = beta_mean, beta_sd = beta_sd),
            class = "stareal_priors")
}

#' @export
print.stareal_priors <- function(x, ...) {
  cat("<stareal_priors>\n")
  cat(sprintf("  sigma   ~ half-N(0, %g^2)\n", x$sigma_sd))
  cat(sprintf("  rho_s   ~ Beta(%g, %g)\n", x$rho_s[1], x$rho_s[2]))
  cat(sprintf("  rho_t   ~ Beta(%g, %g)\n", x$rho_t[1], x$rho_t[2]))
  cat(sprintf("  sigma_n ~ half-N(0, %g^2)   (if nugget = TRUE)\n", x$sigma_n_sd))
  cat(sprintf("  beta    ~ N(%s, %s^2)\n", paste(format(x$beta_mean), collapse = "/"),
              paste(format(x$beta_sd), collapse = "/")))
  invisible(x)
}

#' Settings for the MCMC sampler
#'
#' Each iteration makes `ess_sweeps` updates of the latent values of the
#' observed cells (elliptical slice sampling) and `theta_moves` pairs of
#' Metropolis moves on the hyperparameters. The proposals adapt during
#' burn-in only. The defaults are a pilot run; a run for publication is
#' typically much longer (Africa_areal used 40,000 iterations with 8,000
#' burn-in and `theta_thin = 10`). Use [diagnose()] to judge whether a run
#' is long enough.
#'
#' @param n_iter Iterations per chain, burn-in included.
#' @param n_burn Burn-in iterations (adaptation; discarded).
#' @param n_store Posterior draws kept per chain, evenly spaced after
#'   burn-in. These are the draws used for prediction.
#' @param n_chains Number of chains.
#' @param ess_sweeps Latent updates per iteration.
#' @param theta_moves Hyperparameter move pairs per iteration.
#' @param theta_thin Keep every `theta_thin`-th iteration of the
#'   hyperparameter trace (burn-in included).
#' @param seed Chain `c` starts from seed `seed * 1000 + c`. The caller's
#'   random-number stream is left untouched.
#' @param parallel Run chains at the same time, one per worker process.
#' @param n_cores Maximum number of worker processes (default: one per
#'   chain, at most the number of cores less one).
#' @param blas_threads BLAS threads per worker. One is safest; with an
#'   optimised BLAS (MKL, OpenBLAS, Accelerate) and large fits, more may be
#'   faster, as long as `n_cores * blas_threads` does not exceed the cores.
#' @param checkpoint_dir A folder to save each chain's state in. If a fit is
#'   interrupted, running it again with the same folder resumes every chain
#'   where it stopped, and gives exactly the draws of an uninterrupted run.
#'   A finished fit run again is read back from the folder. `NULL`: no
#'   checkpoints.
#' @param checkpoint_every Iterations between checkpoints.
#' @param progress Print progress messages.
#' @return An object of class `stareal_control`.
#' @examples
#' mcmc_control(n_iter = 500, n_burn = 200, n_store = 100, n_chains = 2)
#' @export
mcmc_control <- function(n_iter = 3000, n_burn = 1000, n_store = 250, n_chains = 4, ess_sweeps = 10,
                         theta_moves = 2, theta_thin = 1, seed = 1, parallel = TRUE, n_cores = NULL,
                         blas_threads = 1, checkpoint_dir = NULL, checkpoint_every = 500, progress = TRUE) {
  whole <- function(x, nm, min = 1) {
    if (!is.numeric(x) || length(x) != 1 || is.na(x) || x != round(x) || x < min)
      stop(sprintf("`%s` must be a whole number >= %d", nm, min), call. = FALSE)
    as.integer(x)
  }
  n_iter <- whole(n_iter, "n_iter"); n_burn <- whole(n_burn, "n_burn", 0)
  n_store <- whole(n_store, "n_store"); n_chains <- whole(n_chains, "n_chains")
  ess_sweeps <- whole(ess_sweeps, "ess_sweeps"); theta_moves <- whole(theta_moves, "theta_moves")
  theta_thin <- whole(theta_thin, "theta_thin"); seed <- whole(seed, "seed", 0)
  blas_threads <- whole(blas_threads, "blas_threads"); checkpoint_every <- whole(checkpoint_every, "checkpoint_every")
  if (n_burn >= n_iter) stop("`n_burn` must be less than `n_iter`", call. = FALSE)
  if (n_store > n_iter - n_burn)
    stop(sprintf("`n_store` (%d) cannot exceed the %d iterations after burn-in", n_store, n_iter - n_burn),
         call. = FALSE)
  if (theta_thin > n_iter) stop("`theta_thin` cannot exceed `n_iter`", call. = FALSE)
  if (is.null(n_cores)) n_cores <- max(1L, min(n_chains, parallel::detectCores() - 1L, na.rm = TRUE))
  n_cores <- whole(n_cores, "n_cores")
  if (!is.null(checkpoint_dir) && (!is.character(checkpoint_dir) || length(checkpoint_dir) != 1))
    stop("`checkpoint_dir` must be one folder path", call. = FALSE)
  structure(list(n_iter = n_iter, n_burn = n_burn, n_store = n_store, n_chains = n_chains,
                 ess_sweeps = ess_sweeps, theta_moves = theta_moves, theta_thin = theta_thin, seed = seed,
                 parallel = isTRUE(parallel), n_cores = n_cores, blas_threads = blas_threads,
                 checkpoint_dir = checkpoint_dir, checkpoint_every = checkpoint_every,
                 progress = isTRUE(progress)),
            class = "stareal_control")
}

#' @export
print.stareal_control <- function(x, ...) {
  cat("<stareal_control>\n")
  cat(sprintf("  %d chain%s x %d iterations (%d burn-in), %d draws kept per chain\n", x$n_chains,
              if (x$n_chains > 1) "s" else "", x$n_iter, x$n_burn, x$n_store))
  cat(sprintf("  per iteration: %d latent sweeps, %d theta move pairs; theta trace every %d\n",
              x$ess_sweeps, x$theta_moves, x$theta_thin))
  cat(sprintf("  seed %d; %s; checkpoints %s\n", x$seed,
              if (x$parallel && x$n_chains > 1) sprintf("parallel on up to %d workers (%d BLAS thread%s each)",
                                                       x$n_cores, x$blas_threads, if (x$blas_threads > 1) "s" else "")
              else "chains run one after another",
              if (is.null(x$checkpoint_dir)) "off" else sprintf("in %s every %d iterations", x$checkpoint_dir,
                                                               x$checkpoint_every)))
  invisible(x)
}
