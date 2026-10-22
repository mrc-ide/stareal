# model.R
#
# The model: covariance pieces, the observed-cell covariance, priors and the
# fixed-effect posterior. Internal; no I/O.
#
# Ported from Africa_areal/R/model.R (2026-10-12). Changes: the design
# matrix X of the observed cells is an input (any fixed effects, not only
# intercept + trend); the prior closure takes the package's prior names; the
# likelihood moved to likelihood.R.
#
# Notation:
#   m areas, T time steps, N = m*T cells; n observed cells with area index a_o
#   and time index t_o (both 1-based). Sigma = Sigma_t (x) Sigma_s, so
#   Sigma_oo = Sigma_t[t_o, t_o] * Sigma_s[a_o, a_o] elementwise.
#   Spatial: Leroux, Q_s = tau (rho_s R + (1 - rho_s) I), R = D - W = V L V'.
#   Temporal: AR1, Sigma_t[t, t'] = rho_t^|t - t'|, variance fixed at 1 (only
#   the product of the two variances is identifiable; tau carries the scale).
#   Mean: x = X beta + field, beta ~ N(b0, B), marginalised into the
#   covariance as X B X'. theta = (log sigma, logit rho_s, logit rho_t
#   [, log sigma_n]) with sigma^2 = 1 / tau; sigma_n is a cell-level nugget
#   on observed cells only (the predicted surface is the field without it).
#
# The model is fitted on the observed cells only (n latent values); the full
# area x time table is reconstructed afterwards by conditional simulation
# (predict.R). This is exact, not an approximation (tests/testthat/test-model.R).
#
# ------------------------------------------------------------------

# ---- covariance pieces ------------------------------------------------------

# eigenvalues of Sigma_s = Q_s^{-1}, given the Laplacian eigenvalues
leroux_g <- function(lambda, tau, rho_s) {
  1 / (tau * (rho_s * lambda + 1 - rho_s))
}

# Sigma_s[rows, rows] = V[rows,] diag(g) V[rows,]'  (only the areas asked for)
sigma_s_block <- function(V, g, rows) {
  Vr <- V[rows, , drop = FALSE]
  (Vr * rep(g, each = length(rows))) %*% t(Vr)
}

# Sigma_t for a stationary AR1 with unit variance
ar1_cov <- function(T, rho_t) {
  toeplitz(rho_t^(0:(T - 1)))      # = rho_t^|t - t'|, without T^2 calls to pow
}

# lower Cholesky factor of ar1_cov(T, rho_t) in closed form:
# L[1,1] = 1, L[t,1] = rho^(t-1), and for j >= 2, L[t,j] = rho^(t-j) sqrt(1-rho^2)
ar1_chol <- function(T, rho_t) {
  L <- rho_t^pmax(outer(seq_len(T), seq_len(T), "-"), 0)
  L[lower.tri(L, diag = TRUE) & col(L) > 1] <-
    L[lower.tri(L, diag = TRUE) & col(L) > 1] * sqrt(1 - rho_t^2)
  L[upper.tri(L)] <- 0
  L
}

# the default fixed effects: intercept + centred, scaled time index (-1 to 1
# over the window). Observed cells and the full table use the same function
# so the two can never disagree.
design_X <- function(t, T) {
  cbind(1, (t - (T + 1) / 2) / (T / 2))
}

# ---- parameters -------------------------------------------------------------

# theta on the unconstrained scale -> natural parameters
theta_to_par <- function(theta) {
  list(sigma = exp(theta[1]), tau = exp(-2 * theta[1]),
       rho_s = plogis(theta[2]), rho_t = plogis(theta[3]),
       sigma_n = if (length(theta) >= 4) exp(theta[4]) else 0)
}

theta_names <- function(nugget) {
  c("log_sigma", "logit_rho_s", "logit_rho_t", if (nugget) "log_sigma_n")
}

# the same, as columns of natural-scale values (for traces and summaries)
theta_natural <- function(th) {
  th <- as.matrix(th)
  out <- cbind(sigma = exp(th[, 1]), rho_s = plogis(th[, 2]), rho_t = plogis(th[, 3]))
  if (ncol(th) >= 4) out <- cbind(out, sigma_n = exp(th[, 4]))
  out
}

# ---- the observed-cell covariance -------------------------------------------

# the observed-cell bookkeeping: unique areas, inverse index, design matrix
# (X: n x p fixed-effect design of the observed cells; intercept + trend if
# not given)
make_obs <- function(a_o, t_o, T, X = design_X(t_o, T)) {
  stopifnot(length(a_o) == length(t_o), nrow(X) == length(a_o))
  ua <- sort(unique(a_o))
  list(a_o = a_o, t_o = t_o, T = T, n = length(a_o),
       ua = ua, inv = match(a_o, ua), X = X)
}

# Everything the sampler needs for one value of theta: the marginal
# covariance of x_o (field + nugget + fixed effects), its lower Cholesky
# factor and log-determinant.
build_sigma_oo <- function(theta, obs, eig, B, jitter = 1e-8) {
  p <- theta_to_par(theta)
  g <- leroux_g(eig$values, p$tau, p$rho_s)
  Ss_u <- sigma_s_block(eig$vectors, g, obs$ua)          # observed areas only
  St <- ar1_cov(obs$T, p$rho_t)
  S <- Ss_u[obs$inv, obs$inv] * St[obs$t_o, obs$t_o] + obs$X %*% B %*% t(obs$X)
  dg <- seq.int(1L, length(S), by = nrow(S) + 1L)        # diagonal, written in place
  S[dg] <- S[dg] + jitter + p$sigma_n^2                  # (diag<- would copy the n x n matrix)
  L <- t(chol(S))
  list(S = S, L = L, logdet = 2 * sum(log(diag(L))), g = g, St = St, par = p)
}

# log N(x | mu, S) from a lower Cholesky factor
log_dmvnorm_chol <- function(x, mu, L, logdet) {
  z <- forwardsolve(L, x - mu)
  -0.5 * (length(x) * log(2 * pi) + logdet + sum(z^2))
}

# ---- prior ------------------------------------------------------------------

# log prior density of theta on the unconstrained scale, Jacobians included,
# as a closure over the prior settings: sigma ~ half-N(0, sigma_sd),
# rho_s ~ Beta(rho_s[1], rho_s[2]), rho_t ~ Beta(rho_t[1], rho_t[2]),
# sigma_n ~ half-N(0, sigma_n_sd)
make_log_prior <- function(priors) {
  sd_s <- priors$sigma_sd; rs <- priors$rho_s; rt <- priors$rho_t
  sd_n <- priors$sigma_n_sd
  function(theta) {
    p <- theta_to_par(theta)
    lp <- dnorm(p$sigma, 0, sd_s, log = TRUE) + log(2) + theta[1] +          # half-normal, + log sigma
      dbeta(p$rho_s, rs[1], rs[2], log = TRUE) + log(p$rho_s * (1 - p$rho_s)) +
      dbeta(p$rho_t, rt[1], rt[2], log = TRUE) + log(p$rho_t * (1 - p$rho_t))
    if (length(theta) >= 4) lp <- lp + dnorm(p$sigma_n, 0, sd_n, log = TRUE) + log(2) + theta[4]
    lp
  }
}

# ---- fixed effects, recovered after the fact --------------------------------

# beta | x_o, theta ~ N(m, V) with prior N(b0, B) and x_o = X beta + e,
# e ~ N(0, S_field). Takes the field-only covariance (no X B X' term).
beta_posterior <- function(x_o, X, S_field, b0, B) {
  Si_X <- solve(S_field, X)
  V <- solve(solve(B) + t(X) %*% Si_X)
  m <- V %*% (solve(B, b0) + t(Si_X) %*% x_o)
  list(mean = as.vector(m), cov = V)
}
