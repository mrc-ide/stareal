# likelihood.R
#
# The observation model as one object. The sampler touches the data only
# through these pieces:
#   loglik(x)       log-likelihood of the observed cells given latent values x
#   pseudo          Gaussian pseudo-data g* on the latent scale, and
#   pseudo_var      their noise variances S_aux,
#   pseudo_at(x)    pseudo-data refitted at latent values x (a Gaussian
#                   approximation of the likelihood around x)
# g* and S_aux set the posterior-scaled ellipse of the latent ESS and the
# surrogate data of the theta move. Any FIXED g* and S_aux leave both moves
# exact; they only set how well the moves fit the posterior. The chain
# starts from the empirical pseudo-data (pseudo, pseudo_var) and refits
# them with pseudo_at() during burn-in only (chains.R), so the chain after
# burn-in is exact. A new likelihood (Poisson, Gaussian) is a new
# constructor here; nothing in sampler.R changes.
#
# Why refit: the empirical logit of a cell with y = 0 (or y = k) is a
# two-sided Gaussian claim, at log(0.5 / (k + 0.5)) with variance about 2
# whatever k, about a likelihood that is one-sided and nearly flat beyond
# the data. With many such cells (rare markers) the ellipse sits away from
# the posterior and the latent level mixes slowly. A Gaussian approximation
# at the chain's own mean of x is close to the posterior for every cell.
#
# Ported from Africa_areal/R/model.R (loglik_binom, pseudo_logits,
# surrogate_noise, make_loglik), 2026-10-12. Changes: a stable softplus in
# the log-likelihood (log1p(exp(x)) overflowed for x > 709); pseudo_at().
#
# ------------------------------------------------------------------

# log(1 + exp(x)), without overflow
softplus <- function(x) pmax(x, 0) + log1p(exp(-abs(x)))

# binomial log-likelihood of the pooled cell counts given logit-prevalence x
loglik_binom <- function(x, y, k) {
  sum(y * x - k * softplus(x))
}

# empirical logits, and the inverse binomial information at a shrunk
# prevalence estimate
pseudo_logits <- function(y, k) log((y + 0.5) / (k - y + 0.5))
surrogate_noise <- function(y, k) {
  p <- (y + 0.5) / (k + 1)
  1 / (k * p * (1 - p))
}

# The binomial log-likelihood expanded to second order about x0:
# curvature w = k p (1 - p) at p = plogis(x0), centre x0 + (y - k p) / w.
# A cell whose likelihood is nearly flat there gets a large variance
# (capped), so it barely shapes the ellipse.
pseudo_binom_at <- function(x0, y, k, var_max = 1e4) {
  p <- stats::plogis(x0)
  w <- k * p * (1 - p)
  w <- pmax(w, 1 / var_max)
  list(pseudo = x0 + (y - k * p) / w, pseudo_var = pmin(1 / w, var_max))
}

# the Gaussian pseudo-likelihood log N(g*; x, S_aux), up to a constant
log_gauss_pseudo <- function(x, g_star, S_aux) -0.5 * sum((g_star - x)^2 / S_aux)

# y positives out of k tested, one entry per observed cell (pooled)
lik_binomial <- function(y, k) {
  stopifnot(length(y) == length(k), all(k >= 1), all(y >= 0), all(y <= k))
  force(y); force(k)
  structure(list(family = "binomial", n = length(y), y = y, k = k,
                 # closures over y and k only, so they are cheap to send to workers
                 loglik = function(x) loglik_binom(x, y, k),
                 pseudo_at = function(x) pseudo_binom_at(x, y, k),
                 pseudo = pseudo_logits(y, k), pseudo_var = surrogate_noise(y, k)),
            class = "stareal_lik")
}
