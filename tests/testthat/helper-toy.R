# helper-toy.R
#
# A small frame shared by the tests: a 5 x 5 lattice of areas, 7 time steps,
# 30 observed cells, and the dense full-precision form of the same model,
# against which the observed-cell machinery is checked. Ported from
# Africa_areal/tests/helper_toy.R (2026-10-12), plus a general design matrix
# (intercept, trend and a covariate varying by area and time) for the
# covariate paths.
#
# ------------------------------------------------------------------

make_toy <- function(seed = 1) {
  set.seed(seed)
  k_side <- 5; m <- k_side^2; T <- 7
  W <- matrix(0, m, m)
  for (i in 1:k_side) for (j in 1:k_side) {
    a <- (i - 1) * k_side + j
    if (i < k_side) W[a, a + k_side] <- W[a + k_side, a] <- 1
    if (j < k_side) W[a, a + 1] <- W[a + 1, a] <- 1
  }
  R_lap <- diag(rowSums(W)) - W
  eig <- eigen(R_lap, symmetric = TRUE)
  theta <- c(log(0.8), qlogis(0.9), qlogis(0.7))
  p <- theta_to_par(theta)
  b0 <- c(-1, 0.5); B <- diag(c(1, 0.5))

  # dense full precision, cell index c = (t - 1) * m + a  (time-major)
  Q_s <- p$tau * (p$rho_s * R_lap + (1 - p$rho_s) * diag(m))
  Q_t <- solve(ar1_cov(T, p$rho_t))
  Q <- kronecker(Q_t, Q_s)
  Sigma_full <- solve(Q)
  X_full <- design_X(rep(seq_len(T), each = m), T)

  n <- 30
  obs_cells <- sort(sample(m * T, n))
  a_o <- (obs_cells - 1) %% m + 1; t_o <- (obs_cells - 1) %/% m + 1
  obs <- make_obs(a_o, t_o, T)
  sig <- build_sigma_oo(theta, obs, eig, B, jitter = 0)
  mu_o <- as.vector(obs$X %*% b0)

  # the general case: a third fixed effect varying by area and time
  Xg_full <- cbind(X_full, rnorm(m * T))
  b0g <- c(b0, 0.3); Bg <- diag(c(1, 0.5, 0.8))
  design_g <- function(rows) Xg_full[as.vector(outer(rows, (seq_len(T) - 1) * m, "+")), , drop = FALSE]
  obs_g <- make_obs(a_o, t_o, T, X = Xg_full[obs_cells, , drop = FALSE])

  list(k_side = k_side, m = m, T = T, W = W, R_lap = R_lap, eig = eig, theta = theta, p = p,
       b0 = b0, B = B, Q = Q, Sigma_full = Sigma_full, X_full = X_full,
       Sigma_full_fx = Sigma_full + X_full %*% B %*% t(X_full),
       n = n, obs_cells = obs_cells, a_o = a_o, t_o = t_o, obs = obs, sig = sig, mu_o = mu_o,
       Xg_full = Xg_full, b0g = b0g, Bg = Bg, design_g = design_g, obs_g = obs_g,
       Sigma_full_g = Sigma_full + Xg_full %*% Bg %*% t(Xg_full),
       priors = list(sigma_sd = 2, rho_s = c(2, 1), rho_t = c(1, 1), sigma_n_sd = 1))
}

# zero noise in prediction: Matheron's rule then returns the conditional mean
zero_noise <- function(seed, m, T, p = 0, n = 0) {
  list(Z = matrix(0, m, T), bz = rep(0, p), ez = rep(0, n))
}
