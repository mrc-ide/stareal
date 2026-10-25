# Checks of model.R and likelihood.R against the dense full-precision model
# (helper-toy.R). Ported from Africa_areal/tests/test_model.R, plus the
# general design matrix and the likelihood object.

tz <- make_toy()

test_that("Sigma_oo from Kronecker indexing equals the inverse Schur complement of the full precision", {
  with(tz, {
    Q_oo <- Q[obs_cells, obs_cells]; Q_ou <- Q[obs_cells, -obs_cells]; Q_uu <- Q[-obs_cells, -obs_cells]
    schur_inv <- solve(Q_oo - Q_ou %*% solve(Q_uu, t(Q_ou)))
    err <- max(abs(schur_inv + obs$X %*% B %*% t(obs$X) - sig$S)) / max(abs(sig$S))
    expect_lt(err, 1e-10)
  })
})

test_that("Sigma_oo is exact with a general design matrix (covariate varying by area and time)", {
  with(tz, {
    sig_g <- build_sigma_oo(theta, obs_g, eig, Bg, jitter = 0)
    err <- max(abs(Sigma_full_g[obs_cells, obs_cells] - sig_g$S)) / max(abs(sig_g$S))
    expect_lt(err, 1e-10)
  })
})

test_that("the nugget adds exactly sigma_n^2 I", {
  with(tz, {
    sn <- 0.7
    sig4 <- build_sigma_oo(c(theta, log(sn)), obs, eig, B, jitter = 0)
    expect_lt(max(abs(sig4$S - (sig$S + diag(sn^2, n)))) / max(abs(sig$S)), 1e-12)
  })
})

test_that("ar1_chol is the Cholesky factor of ar1_cov", {
  Lt <- ar1_chol(tz$T, tz$p$rho_t)
  expect_lt(max(abs(Lt %*% t(Lt) - ar1_cov(tz$T, tz$p$rho_t))), 1e-12)
})

test_that("the log prior matches the stated densities, Jacobians included", {
  lp <- make_log_prior(tz$priors)
  th4 <- c(log(1.3), qlogis(0.8), qlogis(0.6), log(0.4))
  direct <- log(2 * dnorm(1.3, 0, 2)) + log(1.3) +
    log(dbeta(0.8, 2, 1)) + log(0.8 * 0.2) +
    log(dbeta(0.6, 1, 1)) + log(0.6 * 0.4) +
    log(2 * dnorm(0.4, 0, 1)) + log(0.4)
  expect_lt(abs(lp(th4) - direct), 1e-12)
})

test_that("beta_posterior matches the dense conjugate posterior", {
  with(tz, {
    x_o <- mu_o + as.vector(sig$L %*% rnorm(n))
    S_field <- sig$S - obs$X %*% B %*% t(obs$X)
    bp <- beta_posterior(x_o, obs$X, S_field, b0, B)
    # joint Gaussian of (beta, x_o), conditioned on x_o
    C_bx <- B %*% t(obs$X)
    m_direct <- b0 + C_bx %*% solve(sig$S, x_o - mu_o)
    V_direct <- B - C_bx %*% solve(sig$S, t(C_bx))
    expect_lt(max(abs(bp$mean - m_direct)), 1e-8)
    expect_lt(max(abs(bp$cov - V_direct)), 1e-8)
  })
})

test_that("the binomial likelihood object", {
  y <- c(0, 3, 10); k <- c(5, 10, 10); x <- c(-2, 0.1, 3)
  lk <- lik_binomial(y, k)
  expect_s3_class(lk, "stareal_lik")
  # log-likelihood up to the binomial coefficients
  expect_equal(lk$loglik(x), sum(dbinom(y, k, plogis(x), log = TRUE) - lchoose(k, y)))
  expect_true(all(is.finite(lk$pseudo)) && all(lk$pseudo_var > 0))
  expect_error(lik_binomial(c(1, 11), c(10, 10)))
  expect_error(lik_binomial(1, 0))
  # no overflow for large |x|
  expect_true(is.finite(lik_binomial(5, 10)$loglik(800)) && is.finite(lik_binomial(5, 10)$loglik(-800)))
})

test_that("the refitted pseudo-data are the likelihood's second-order expansion", {
  y <- c(0, 3, 10, 40); k <- c(50, 10, 10, 60); x0 <- c(-4, 0.2, 2.5, 0.4)
  lk <- lik_binomial(y, k)
  ps <- lk$pseudo_at(x0)
  # numerical first and second derivatives of each cell's log-likelihood at x0
  f <- function(x) y * x - k * softplus(x)
  h <- 1e-4
  d1 <- (f(x0 + h) - f(x0 - h)) / (2 * h); d2 <- (f(x0 + h) - 2 * f(x0) + f(x0 - h)) / h^2
  # the Gaussian -(x - g)^2 / (2 S) has slope (g - x0) / S and curvature -1 / S at x0
  expect_equal(ps$pseudo_var, -1 / d2, tolerance = 1e-5)
  expect_equal((ps$pseudo - x0) / ps$pseudo_var, d1, tolerance = 1e-5)
  # where the likelihood is flat the variance is large (capped), so the cell barely shapes the ellipse
  expect_equal(lk$pseudo_at(c(-30, 0, 0, 0))$pseudo_var[1], 1e4)
})
