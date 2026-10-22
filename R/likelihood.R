# likelihood.R
#
# The observation model as one object. The sampler touches the data only
# through these three pieces:
#   loglik(x)   log-likelihood of the observed cells given latent values x
#   pseudo      fixed Gaussian pseudo-data g* on the latent scale, and
#   pseudo_var  their noise variances S_aux,
# which set the posterior-scaled ellipse of the latent ESS and the surrogate
# data of the theta move. Any fixed g* and S_aux leave both moves exact; they
# only set how well the moves fit the posterior. A new likelihood (Poisson,
# Gaussian) is a new constructor here; nothing in sampler.R changes.
#
# Ported from Africa_areal/R/model.R (loglik_binom, pseudo_logits,
# surrogate_noise, make_loglik), 2026-10-12.
#
# ------------------------------------------------------------------

# binomial log-likelihood of the pooled cell counts given logit-prevalence x
loglik_binom <- function(x, y, k) {
  sum(y * x - k * log1p(exp(x)))
}

# empirical logits, and the inverse binomial information at a shrunk
# prevalence estimate
pseudo_logits <- function(y, k) log((y + 0.5) / (k - y + 0.5))
surrogate_noise <- function(y, k) {
  p <- (y + 0.5) / (k + 1)
  1 / (k * p * (1 - p))
}

# the Gaussian pseudo-likelihood log N(g*; x, S_aux), up to a constant
log_gauss_pseudo <- function(x, g_star, S_aux) -0.5 * sum((g_star - x)^2 / S_aux)

# y positives out of k tested, one entry per observed cell (pooled)
lik_binomial <- function(y, k) {
  stopifnot(length(y) == length(k), all(k >= 1), all(y >= 0), all(y <= k))
  force(y); force(k)
  structure(list(family = "binomial", n = length(y), y = y, k = k,
                 # a closure over y and k only, so it is cheap to send to workers
                 loglik = function(x) loglik_binom(x, y, k),
                 pseudo = pseudo_logits(y, k), pseudo_var = surrogate_noise(y, k)),
            class = "stareal_lik")
}
