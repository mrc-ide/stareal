# crossval.R
#
# K-fold cross-validation over observed cells: out-of-sample checks and
# model comparison. Cells are split at random into K folds; each fold is
# held out in turn, the model is refitted to the rest with the same priors
# and sampler settings, and the held-out cells are predicted by kriging
# from each stored draw. A held-out cell is an observation, so its
# prediction carries the nugget.
#
# Ported from Africa_areal (pipeline/09_crossval.R; predict_heldout and
# score_heldout in R/predict.R), 2026-10-12. Changes: all folds (or chosen
# ones) in one call, returned rather than written; the held-out design
# comes from design_cells(), so covariates carry through; the fold assignment
# seed and the per-fold seed (k, as in Africa_areal) restore the caller's
# random-number stream. The pointwise log-likelihood of fitted cells is not
# used for comparison: with one latent value per cell, PSIS-LOO on it is
# unreliable.
#
# ------------------------------------------------------------------

#' Cross-validate a fit
#'
#' Splits the observed cells at random into `folds` groups. Each group is
#' held out in turn: the model is refitted to the other cells, with the
#' fit's mean, covariates, priors, nugget and sampler settings, and then predicts the
#' held-out cells. Use it to check calibration out of sample and to compare
#' models fitted to the same cells. The same `seed` gives the same folds,
#' so models are compared on the same splits. The mean log score is higher
#' for the better model; the share of observations outside their 95%
#' posterior-predictive interval should be near 5%.
#'
#' Each fold is a full refit, so this costs `folds` (or `length(fold)`)
#' times the fit.
#'
#' @param fit A [fit_areal()] object.
#' @param folds Number of folds.
#' @param fold Which folds to run (default all), e.g. to run them in
#'   separate sessions.
#' @param control Sampler settings for the refits (default: the fit's). If
#'   it has a `checkpoint_dir`, each fold checkpoints in its own subfolder.
#' @param seed Seed of the random split into folds.
#' @return An object of class `stareal_cv`:
#'   \describe{
#'     \item{cells}{one row per held-out cell: `fold`, `cell_id`, `id`,
#'       `t`, `time`, `y`, `k`, `obs_prev` (`y / k`), `pred_mean`, `pred_lo`,
#'       `pred_hi` (95% interval of the cell's prevalence), `rep_lo`,
#'       `rep_hi` (95% posterior-predictive interval of `y / k`), `outside`
#'       (the observation is outside it), and `log_score`, the log
#'       predictive density of `y`}
#'     \item{folds}{one row per fold, and `overall`: the number held out,
#'       RMSE of `pred_mean` against `obs_prev`, percent outside the 95%
#'       interval, mean log score and correlation}
#'   }
#' @examples
#' \donttest{
#' fr <- areal_lattice(4, 4)
#' tf <- time_frame(T = 5)
#' sim <- simulate_areal(fr, tf, n_cells = 40, seed = 1)
#' d <- areal_data(sim$data, fr, tf, y = "positive", k = "tested", area = "id", t = "t",
#'                 verbose = FALSE)
#' ctl <- mcmc_control(n_iter = 100, n_burn = 50, n_store = 25, n_chains = 1,
#'                     ess_sweeps = 2, progress = FALSE)
#' fit <- fit_areal(d, control = ctl)
#' cv <- crossval(fit, folds = 4)
#' cv
#' }
#' @export
crossval <- function(fit, folds = 10, fold = NULL, control = fit$control, seed = 20261012) {
  if (!inherits(fit, "stareal_fit")) stop("`fit` must be made by fit_areal()", call. = FALSE)
  if (!inherits(control, "stareal_control")) stop("`control` must be made by mcmc_control()", call. = FALSE)
  cells <- fit$data$cells; n <- nrow(cells)
  folds <- as.integer(folds)
  if (is.na(folds) || folds < 2 || folds > n) stop(sprintf("`folds` must be between 2 and the %d cells", n),
                                                   call. = FALSE)
  if (is.null(fold)) fold <- seq_len(folds)
  fold <- as.integer(fold)
  if (any(!fold %in% seq_len(folds))) stop(sprintf("`fold` must be among 1..%d", folds), call. = FALSE)
  fold_of <- with_seed(seed, sample(rep_len(seq_len(folds), n)))   # the same folds for every model
  base_dir <- control$checkpoint_dir
  scored <- list()
  for (k in fold) {
    held <- cells[fold_of == k, ]
    train <- fit$data
    train$cells <- cells[fold_of != k, ]
    ctl <- control
    if (!is.null(base_dir)) ctl$checkpoint_dir <- file.path(base_dir, sprintf("fold_%d_of_%d", k, folds))
    if (control$progress) message(sprintf("fold %d of %d: fitting %d cells, holding out %d", k, folds,
                                          nrow(train$cells), nrow(held)))
    fk <- fit_with_design(train, fit$design, fit$nugget, fit$priors, ctl)
    sc <- with_seed(k, {
      xh <- predict_heldout(fk, held)
      score_heldout(xh, held)
    })
    scored[[length(scored) + 1]] <- cbind(fold = k, sc)
  }
  out <- do.call(rbind, scored)
  rownames(out) <- NULL
  summ <- function(s) data.frame(n_held = nrow(s), rmse = sqrt(mean((s$obs_prev - s$pred_mean)^2)),
                                 pct_outside_95 = 100 * mean(s$outside), mean_log_score = mean(s$log_score),
                                 correlation = if (nrow(s) > 2) stats::cor(s$obs_prev, s$pred_mean) else NA_real_)
  per_fold <- do.call(rbind, lapply(split(out, out$fold), function(s) cbind(fold = s$fold[1], summ(s))))
  structure(list(cells = out, folds = per_fold, overall = summ(out), n_folds = folds, seed = seed),
            class = "stareal_cv")
}

# For each stored draw of `fit`, a draw of the held-out cells' logit values
# given (x_o, theta), by the kriging equations, with the cell nugget.
# Returns a D x h matrix (draws pooled chain by chain).
predict_heldout <- function(fit, held) {
  obs <- fit$setup$obs; eig <- fit$data$frame$eig
  pool <- function(a) matrix(a, prod(dim(a)[1:2]), dim(a)[3])
  th <- pool(fit$draws$theta); xo <- pool(fit$draws$x)
  Xh <- design_cells(fit, held$area, held$t)
  t(vapply(seq_len(nrow(xo)), function(d) {
    mo <- heldout_moments(th[d, ], xo[d, ], obs, eig, fit$setup$B, fit$setup$b0, Xh, held$area, held$t)
    mo$mean + as.vector(t(chol(mo$cov + diag(1e-8, nrow(mo$cov)))) %*% rnorm(length(mo$mean)))
  }, numeric(nrow(held))))
}

# the conditional mean and covariance of held-out cells (areas a_h, steps
# t_h, design Xh), nugget included, given the observed cells' x_o and theta
heldout_moments <- function(theta, x_o, obs, eig, B, b0, Xh, a_h, t_h, jitter = 1e-8) {
  sig <- build_sigma_oo(theta, obs, eig, B, jitter = jitter)
  mu_o <- as.vector(obs$X %*% b0); mu_h <- as.vector(Xh %*% b0)
  all_a <- sort(unique(c(obs$ua, a_h)))
  ia <- match(a_h, all_a); io <- match(obs$a_o, all_a)
  Ss <- sigma_s_block(eig$vectors, sig$g, all_a)
  S_ho <- Ss[ia, io, drop = FALSE] * sig$St[t_h, obs$t_o, drop = FALSE] + Xh %*% B %*% t(obs$X)
  S_hh <- Ss[ia, ia, drop = FALSE] * sig$St[t_h, t_h, drop = FALSE] + Xh %*% B %*% t(Xh) +
    diag(sig$par$sigma_n^2, length(a_h))
  A <- t(backsolve(t(sig$L), forwardsolve(sig$L, t(S_ho))))       # S_ho Sigma_oo^{-1}
  list(mean = mu_h + as.vector(A %*% (x_o - mu_o)), cov = S_hh - A %*% t(S_ho))
}

# Scores for held-out cells from a D x h matrix of logit draws: predictive
# mean and 95% interval of prevalence, 95% posterior-predictive interval of
# y/k, whether the observation falls outside it, and the log predictive
# density log p(y | data), estimated by averaging the binomial likelihood
# over draws.
score_heldout <- function(xh, held) {
  ph <- stats::plogis(xh)
  yrep <- matrix(rbinom(length(ph), size = rep(held$k, each = nrow(ph)), prob = ph), nrow(ph))
  ld <- matrix(dbinom(rep(held$y, each = nrow(ph)), rep(held$k, each = nrow(ph)), ph, log = TRUE), nrow(ph))
  mx <- apply(ld, 2, max)
  out <- data.frame(held[c("cell_id", "id", "t", "time", "y", "k")], obs_prev = held$y / held$k,
                    pred_mean = colMeans(ph),
                    pred_lo = apply(ph, 2, quantile, 0.025), pred_hi = apply(ph, 2, quantile, 0.975),
                    rep_lo = apply(yrep, 2, quantile, 0.025) / held$k,
                    rep_hi = apply(yrep, 2, quantile, 0.975) / held$k,
                    log_score = mx + log(colMeans(exp(sweep(ld, 2, mx)))))
  out$outside <- out$obs_prev < out$rep_lo | out$obs_prev > out$rep_hi
  out
}

#' @export
print.stareal_cv <- function(x, ...) {
  cat(sprintf("<stareal_cv> %d-fold cross-validation, %d fold%s run, %d cells held out\n", x$n_folds,
              nrow(x$folds), if (nrow(x$folds) > 1) "s" else "", nrow(x$cells)))
  o <- x$overall
  cat(sprintf("  RMSE %.4f; %.1f%% outside the 95%% predictive interval; mean log score %.3f; correlation %.2f\n",
              o$rmse, o$pct_outside_95, o$mean_log_score, o$correlation))
  invisible(x)
}
