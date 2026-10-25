# crossval.R
#
# K-fold cross-validation over observed cells: out-of-sample checks and
# model comparison. Cells are split at random into K folds; each fold is
# held out in turn, the model is refitted to the rest with the same model
# and sampler settings, and the held-out cells are predicted by kriging
# from each stored draw. A held-out cell is an observation, so its
# prediction carries the nugget.
#
# Ported from Africa_areal (pipeline/09_crossval.R; predict_heldout and
# score_heldout in R/predict.R), 2026-10-12. Changes: all folds (or chosen
# ones) in one call, returned rather than written; the held-out design
# comes from design_cells(), so covariates carry through; the log score
# integrates the held-out cell's Gaussian exactly (Gauss-Hermite) instead of
# averaging the likelihood over one random draw per posterior draw, which
# was biased low when k is large; whether an observation falls outside its
# predictive interval is judged by a randomised PIT, so that a calibrated
# model gives 5% whatever k (with discrete counts, the plain interval test
# gives less). The pointwise log-likelihood of fitted cells is not used for
# comparison: with one latent value per cell, PSIS-LOO on it is unreliable.
#
# ------------------------------------------------------------------

#' Cross-validate a fit
#'
#' Splits the observed cells at random into `folds` groups. Each group is
#' held out in turn: the model is refitted to the other cells, with the
#' fit's formula, covariates, priors, nugget and sampler settings, and then
#' predicts the held-out cells. Use it to check calibration out of sample
#' and to compare models fitted to the same cells. The same `seed` gives
#' the same folds, so models are compared on the same splits. The mean log
#' score is higher for the better model; about 5% of observations should
#' fall outside their 95% predictive interval.
#'
#' Held-out cells are predicted as observations, with the cell nugget, so
#' `pred_*` describe the prevalence a new survey of the cell would see, not
#' the smoother surface of [predict()][predict.stareal_fit].
#'
#' Each fold is a full refit, so this costs `folds` (or `length(fold)`)
#' times the fit.
#'
#' @param fit A [fit_areal()] object.
#' @param folds The number of folds, assigned at random; or a fold label
#'   for each observed cell, in the order of `fit$data$cells`, for folds
#'   of your own: by site, by country (spatial cross-validation, which tests
#'   extrapolation to unsurveyed places), or by period.
#' @param fold Which folds to run (default all), as numbers or labels, e.g.
#'   to run them in separate sessions.
#' @param control Sampler settings for the refits: [mcmc_control()] or a
#'   list of some of its arguments (default: the fit's).
#' @param checkpoint_dir A folder for the refits' checkpoints: fold `k` of
#'   `K` with split seed `s` checkpoints in
#'   `<checkpoint_dir>/cv_seed<s>_fold_<k>_of_<K>`. `NULL`: none.
#' @param seed Seed of the random split into folds.
#' @return An object of class `stareal_cv`:
#'   \describe{
#'     \item{cells}{one row per held-out cell: `fold`, `cell_id`, `id`,
#'       `t`, `time`, `y`, `k`, `obs_prev` (`y / k`), `pred_mean`, `pred_lo`,
#'       `pred_hi` (mean and 95% interval of the cell's prevalence, with the
#'       nugget), `rep_lo`, `rep_hi` (95% posterior-predictive interval of
#'       `y / k`), `pit` (randomised probability integral transform of
#'       `y`), `outside` (`pit` outside 0.025-0.975), and `log_score`, the
#'       log predictive density of `y`}
#'     \item{folds}{one row per fold, and `overall`: the number held out,
#'       RMSE of `pred_mean` against `obs_prev`, percent outside the 95%
#'       interval, mean log score and correlation}
#'     \item{fold_of}{the fold of every observed cell}
#'   }
#'   Compare models with [compare_cv()].
#' @examples
#' \donttest{
#' fr <- areal_lattice(4, 4)
#' tf <- time_frame(n_steps = 5)
#' sim <- simulate_areal(fr, tf, n_cells = 40, seed = 1)
#' fit <- fit_areal(sim, control = list(n_iter = 100, n_burn = 50, n_store = 25, n_chains = 1,
#'                                      ess_sweeps = 2, verbose = FALSE))
#' cv <- crossval(fit, folds = 4)
#' cv
#' }
#' @export
crossval <- function(fit, folds = 10, fold = NULL, control = fit$control, checkpoint_dir = NULL, seed = 1) {
  if (!inherits(fit, "stareal_fit")) stop("`fit` must be made by fit_areal()", call. = FALSE)
  control <- complete_settings(control, mcmc_control, "control")
  cells <- fit$data$cells; n <- nrow(cells)
  seed <- check_whole(seed, "seed", 0)
  if (length(folds) == 1) {
    folds <- check_whole(folds, "folds", 2)
    if (folds > n) stop(sprintf("`folds` must be between 2 and the %d cells", n), call. = FALSE)
    fold_of <- with_seed(seed, sample(rep_len(seq_len(folds), n)))   # the same folds for every model
  } else {
    if (length(folds) != n || anyNA(folds))
      stop(sprintf("`folds` must be a number of folds, or a fold label for each of the %d observed cells (in the order of fit$data$cells)", n),
           call. = FALSE)
    fold_of <- folds
  }
  labels <- sort(unique(fold_of))
  if (length(labels) < 2) stop("`folds` must give at least two folds", call. = FALSE)
  if (is.null(fold)) fold <- labels
  if (!all(fold %in% labels))
    stop(sprintf("`fold` must be among the fold labels (%s)", paste(utils::head(labels, 10), collapse = ", ")), call. = FALSE)
  fold <- labels[labels %in% fold]
  scored <- list()
  for (k in fold) {
    i <- match(k, labels)
    held <- cells[fold_of == k, ]
    train <- fit$data
    train$cells <- cells[fold_of != k, ]
    ck <- if (is.null(checkpoint_dir)) NULL else
      file.path(checkpoint_dir, sprintf("cv_seed%d_fold_%s_of_%d", seed, gsub("[^A-Za-z0-9_.-]", "_", k), length(labels)))
    if (control$verbose) message(sprintf("fold %s (%d of %d): fitting %d cells, holding out %d", k, i, length(labels),
                                         nrow(train$cells), nrow(held)))
    fk <- fit_with_design(train, fit$design, fit$nugget, fit$priors, control, ck)
    sc <- with_seed(seed + i, score_heldout(predict_heldout(fk, held), held))
    scored[[length(scored) + 1]] <- cbind(fold = rep(k, nrow(sc)), sc)
  }
  out <- do.call(rbind, scored)
  out <- out[order(match(out$fold, labels), out$cell_id), ]
  rownames(out) <- NULL
  summ <- function(s) data.frame(n_held = nrow(s), rmse = sqrt(mean((s$obs_prev - s$pred_mean)^2)),
                                 pct_outside_95 = 100 * mean(s$outside), mean_log_score = mean(s$log_score),
                                 correlation = if (nrow(s) > 2) stats::cor(s$obs_prev, s$pred_mean) else NA_real_)
  per_fold <- do.call(rbind, lapply(fold, function(k) cbind(fold = k, summ(out[out$fold == k, ]))))
  rownames(per_fold) <- NULL
  structure(list(cells = out, folds = per_fold, overall = summ(out), n_folds = length(labels), seed = seed,
                 fold_of = fold_of),
            class = "stareal_cv")
}

#' Compare cross-validated models
#'
#' Compares models cross-validated on the same cells with the same folds
#' (the same `folds` and `seed` in [crossval()]). The log scores are compared
#' cell by cell, which removes most of the variation between cells: the
#' difference from the best model is given with its standard error, the
#' standard deviation of the per-cell differences over the square root of
#' the number of cells. A difference of less than about two standard errors
#' is not a clear preference.
#'
#' @param ... Two or more `stareal_cv` objects, named (e.g.
#'   `compare_cv(default = cv1, no_nugget = cv2)`); unnamed ones are named
#'   by their expressions.
#' @return A data frame, one row per model, best first: `model`, `n_cells`,
#'   `mean_log_score`, `diff` (from the best model), `se_diff`,
#'   `pct_outside_95`, `rmse`, `correlation`.
#' @examples
#' \donttest{
#' sim <- simulate_areal(areal_lattice(4, 4), time_frame(n_steps = 5), n_cells = 40, seed = 1)
#' ctl <- list(n_iter = 100, n_burn = 50, n_store = 25, n_chains = 1, ess_sweeps = 2, verbose = FALSE)
#' fit1 <- fit_areal(sim, control = ctl)
#' fit2 <- fit_areal(sim, nugget = FALSE, control = ctl)
#' compare_cv(nugget = crossval(fit1, folds = 4), no_nugget = crossval(fit2, folds = 4))
#' }
#' @export
compare_cv <- function(...) {
  cvs <- list(...)
  nm <- names(cvs)
  exprs <- vapply(as.list(substitute(list(...)))[-1], function(e) paste(deparse(e), collapse = ""), "")
  if (is.null(nm)) nm <- exprs
  nm[nm == ""] <- exprs[nm == ""]
  if (length(cvs) < 2) stop("give at least two `stareal_cv` objects", call. = FALSE)
  if (!all(vapply(cvs, inherits, logical(1), "stareal_cv"))) stop("every argument must be made by crossval()", call. = FALSE)
  key <- function(cv) paste(cv$cells$fold, cv$cells$cell_id)
  ref <- sort(key(cvs[[1]]))
  if (!all(vapply(cvs, function(cv) identical(sort(key(cv)), ref), logical(1))))
    stop("the models were not cross-validated on the same cells with the same folds; use the same `folds` and `seed`",
         call. = FALSE)
  ls <- vapply(cvs, function(cv) cv$cells$log_score[match(ref, key(cv))], numeric(length(ref)))
  if (is.null(dim(ls))) ls <- matrix(ls, nrow = 1)
  best <- which.max(colMeans(ls))
  d <- ls - ls[, best]
  out <- data.frame(model = nm, n_cells = length(ref), mean_log_score = colMeans(ls), diff = colMeans(d),
                    se_diff = apply(d, 2, stats::sd) / sqrt(length(ref)),
                    pct_outside_95 = vapply(cvs, function(cv) cv$overall$pct_outside_95, 1),
                    rmse = vapply(cvs, function(cv) cv$overall$rmse, 1),
                    correlation = vapply(cvs, function(cv) cv$overall$correlation, 1), stringsAsFactors = FALSE)
  out <- out[order(-out$mean_log_score), ]
  rownames(out) <- NULL
  out
}

# For each stored draw of `fit`, the conditional mean and variance of the
# held-out cells' logits given (x_o, theta), and one draw from that
# Gaussian, by the kriging equations, with the cell nugget. Returns D x h
# matrices mean, var and draw (draws pooled chain by chain).
predict_heldout <- function(fit, held) {
  obs <- make_obs(fit$data$cells$area, fit$data$cells$t, fit$data$time$n_steps,
                  X = design_eval(fit$design, fit$data$cells$area, fit$data$cells$t))
  eig <- fit$data$frame$eig
  pool <- function(a) matrix(a, prod(dim(a)[1:2]), dim(a)[3])
  th <- pool(fit$draws$theta); xo <- pool(fit$draws$x)
  Xh <- design_cells(fit, held$area, held$t)
  B <- diag(fit$priors$beta_sd^2, length(fit$priors$beta_sd))
  D <- nrow(xo); h <- nrow(held)
  out <- list(mean = matrix(NA_real_, D, h), var = matrix(NA_real_, D, h), draw = matrix(NA_real_, D, h))
  for (d in seq_len(D)) {
    mo <- heldout_moments(th[d, ], xo[d, ], obs, eig, B, fit$priors$beta_mean, Xh, held$area, held$t)
    out$mean[d, ] <- mo$mean
    out$var[d, ] <- pmax(diag(mo$cov), 0)
    out$draw[d, ] <- mo$mean + as.vector(t(chol(mo$cov + diag(1e-8, h))) %*% rnorm(h))
  }
  out
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

# Gauss-Hermite nodes and weights for integrals of f(z) exp(-z^2)
# (Golub-Welsch)
gauss_hermite <- function(n) {
  J <- matrix(0, n, n)
  off <- sqrt(seq_len(n - 1) / 2)
  J[cbind(1:(n - 1), 2:n)] <- off; J[cbind(2:n, 1:(n - 1))] <- off
  e <- eigen(J, symmetric = TRUE)
  list(nodes = e$values, weights = sqrt(pi) * e$vectors[1, ]^2)
}

# log of the integral of Binom(y | k, plogis(x)) N(x; mu, s2) dx, for
# matrices mu, s2 (draws x cells) and cell vectors y, k: adaptive
# Gauss-Hermite, centred on the integrand's mode and scaled by its
# curvature, so a narrow likelihood (large k) is resolved.
#
# The mode is found by Newton's method with backtracking: a step is halved
# until the (strictly concave) log integrand rises. Plain Newton overshoots
# when the observed prevalence is far from mu (e.g. y/k = 0.16, mu = 1.8,
# s2 = 4): it oscillated, stopped far from the mode, and the nodes missed
# the mass, giving log scores of -10^4 where the truth is about -6 (found
# 2026-10-25 in the cross-validation tutorial).
log_binom_gauss <- function(y, k, mu, s2, n_nodes = 20) {
  D <- nrow(mu)
  Y <- matrix(y, D, length(y), byrow = TRUE); K <- matrix(k, D, length(k), byrow = TRUE)
  s2 <- pmax(s2, 1e-12)
  lc <- lchoose(K, Y)
  h <- function(x) lc + Y * x - K * softplus(x) - 0.5 * log(2 * pi * s2) - (x - mu)^2 / (2 * s2)
  x <- mu
  hx <- h(x)
  for (it in 1:100) {                                         # damped Newton steps to the mode
    p <- stats::plogis(x)
    g <- Y - K * p - (x - mu) / s2
    H <- -K * p * (1 - p) - 1 / s2
    step <- g / H
    t <- array(1, dim(x))
    for (half in 1:60) {
      x_new <- x - t * step
      h_new <- h(x_new)
      worse <- !(h_new >= hx - 1e-12 * abs(hx))              # NaN counts as worse
      if (!any(worse)) break
      t[worse] <- t[worse] / 2
    }
    moved <- max(abs(t * step))
    x <- x_new; hx <- h_new
    if (moved < 1e-10) break
  }
  p <- stats::plogis(x)
  sd <- 1 / sqrt(K * p * (1 - p) + 1 / s2)
  gh <- gauss_hermite(n_nodes)
  terms <- vapply(seq_len(n_nodes), function(i) {
    z <- gh$nodes[i]
    h(x + sqrt(2) * sd * z) + z^2 + log(gh$weights[i])
  }, matrix(0, D, length(y)))
  terms <- array(terms, c(D, length(y), n_nodes))
  mx <- apply(terms, c(1, 2), max)
  log(sqrt(2) * sd) + mx + log(apply(exp(terms - as.vector(mx)), c(1, 2), sum))
}

# Scores for held-out cells from predict_heldout(): predictive mean and 95%
# interval of prevalence (from the draws), 95% posterior-predictive
# interval of y/k, the randomised PIT of y and whether it is outside
# 0.025-0.975, and the log predictive density log p(y | data), averaged
# over posterior draws with each draw's Gaussian integrated exactly.
score_heldout <- function(ph, held) {
  pd <- stats::plogis(ph$draw)
  D <- nrow(pd)
  yrep <- matrix(rbinom(length(pd), size = rep(held$k, each = D), prob = pd), D)
  lp <- log_binom_gauss(held$y, held$k, ph$mean, ph$var)          # D x h
  mx <- apply(lp, 2, max)
  Fy <- colMeans(matrix(stats::pbinom(rep(held$y, each = D), rep(held$k, each = D), pd), D))
  Fy1 <- colMeans(matrix(stats::pbinom(rep(held$y - 1, each = D), rep(held$k, each = D), pd), D))
  pit <- Fy1 + stats::runif(nrow(held)) * (Fy - Fy1)
  q <- function(M, p) apply(M, 2, stats::quantile, p, names = FALSE)
  data.frame(held[c("cell_id", "id", "t", "time", "y", "k")], obs_prev = held$y / held$k,
             pred_mean = colMeans(pd), pred_lo = q(pd, 0.025), pred_hi = q(pd, 0.975),
             rep_lo = q(yrep, 0.025) / held$k, rep_hi = q(yrep, 0.975) / held$k,
             pit = pit, outside = pit < 0.025 | pit > 0.975,
             log_score = mx + log(colMeans(exp(sweep(lp, 2, mx)))), row.names = NULL)
}

#' @export
print.stareal_cv <- function(x, ...) {
  cat(sprintf("<stareal_cv> %d-fold cross-validation (seed %d), %d fold%s run, %d cells held out\n", x$n_folds,
              x$seed, nrow(x$folds), if (nrow(x$folds) > 1) "s" else "", nrow(x$cells)))
  o <- x$overall
  cat(sprintf("  RMSE %.4f; %.1f%% outside the 95%% predictive interval; mean log score %.3f; correlation %.2f\n",
              o$rmse, o$pct_outside_95, o$mean_log_score, o$correlation))
  invisible(x)
}
