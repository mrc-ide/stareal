# Checks of predict_core.R against the dense model (helper-toy.R). Ported
# from Africa_areal/tests/test_predict.R checks 1-3 (check 4, the two-phase
# prediction and aggregation, arrives with predict()), plus the general
# design matrix and the random-number state.
#   1. ar1_filter(Z) == Z Lt' (the AR1 recursion replaces the T x T product).
#   2. Matheron's rule with all random inputs zero returns the exact
#      conditional (kriging) mean, with and without a nugget, and with a
#      covariate varying by area and time; conditional variances match by
#      Monte Carlo; without a nugget the observed cells are reproduced.
#   3. Rows computed for a subset of areas equal those rows of the whole
#      surface (same draw).

tz <- make_toy()
ctx <- with(tz, list(eig = eig, obs = obs, b0 = b0, B = B, T = T, jitter = 0))
ctx_g <- with(tz, list(eig = eig, obs = obs_g, b0 = b0g, B = Bg, T = T, jitter = 0, design = design_g))
set.seed(21)
x_o <- with(tz, mu_o + as.vector(sig$L %*% rnorm(n)))
theta4 <- c(tz$theta, log(0.7))

test_that("the AR1 recursion equals Z Lt'", {
  Z <- matrix(rnorm(tz$m * tz$T), tz$m, tz$T)
  expect_lt(max(abs(ar1_filter(Z, tz$p$rho_t) - Z %*% t(ar1_chol(tz$T, tz$p$rho_t)))), 1e-12)
})

test_that("Matheron's rule with zero noise gives the exact kriging mean", {
  local_mocked_bindings(draw_noise = zero_noise)
  with(tz, {
    mu_full <- as.vector(X_full %*% b0)
    S_uo <- Sigma_full_fx[, obs_cells]
    m_exact <- mu_full + S_uo %*% solve(Sigma_full_fx[obs_cells, obs_cells], x_o - mu_full[obs_cells])
    expect_lt(max(abs(as.vector(predict_surface(x_o, theta, 1, ctx)) - m_exact)), 1e-8)

    S_oo4 <- Sigma_full_fx[obs_cells, obs_cells] + diag(0.7^2, n)
    m_exact4 <- mu_full + S_uo %*% solve(S_oo4, x_o - mu_full[obs_cells])
    m_pred4 <- as.vector(predict_surface(x_o, theta4, 1, ctx))
    expect_lt(max(abs(m_pred4 - m_exact4)), 1e-8)
    expect_gt(max(abs(m_pred4[obs_cells] - x_o)), 1e-3)   # nugget: data not reproduced

    # general design: intercept, trend and a covariate by area and time
    mu_g <- as.vector(Xg_full %*% b0g)
    S_oo_g <- Sigma_full_g[obs_cells, obs_cells] + diag(0.7^2, n)
    m_exact_g <- mu_g + Sigma_full_g[, obs_cells] %*% solve(S_oo_g, x_o - mu_g[obs_cells])
    expect_lt(max(abs(as.vector(predict_surface(x_o, theta4, 1, ctx_g)) - m_exact_g)), 1e-8)
  })
})

test_that("conditional variances match exact kriging (Monte Carlo); no nugget reproduces the data", {
  skip_on_cran()
  with(tz, {
    S <- 4000
    draws <- vapply(seq_len(S), function(s) as.vector(predict_surface(x_o, theta, s, ctx)), numeric(m * T))
    S_uo <- Sigma_full_fx[, obs_cells]
    C_exact <- Sigma_full_fx - S_uo %*% solve(Sigma_full_fx[obs_cells, obs_cells], t(S_uo))
    expect_lt(max(abs(apply(draws, 1, var) - diag(C_exact))), 0.1 * max(diag(C_exact)))
    expect_lt(max(abs(draws[obs_cells, ] - x_o)), 1e-8)

    draws_g <- vapply(seq_len(S), function(s) as.vector(predict_surface(x_o, theta, s, ctx_g)), numeric(m * T))
    C_g <- Sigma_full_g - Sigma_full_g[, obs_cells] %*%
      solve(Sigma_full_g[obs_cells, obs_cells], t(Sigma_full_g[, obs_cells]))
    expect_lt(max(abs(apply(draws_g, 1, var) - diag(C_g))), 0.1 * max(diag(C_g)))

    # with a nugget: the observed cells are noisy observations of the field,
    # so the conditional variance uses Sigma_oo + sigma_n^2 I (the nugget
    # draw in Matheron's rule; dropping it would shrink these variances)
    draws4 <- vapply(seq_len(S), function(s) as.vector(predict_surface(x_o, theta4, s, ctx)), numeric(m * T))
    C4 <- Sigma_full_fx - S_uo %*% solve(Sigma_full_fx[obs_cells, obs_cells] + diag(0.7^2, n), t(S_uo))
    expect_lt(max(abs(apply(draws4, 1, var) - diag(C4))), 0.1 * max(diag(C4)))
    expect_gt(min(apply(draws4[obs_cells, ], 1, var)), 0.01)        # the data are not reproduced
  })
})

test_that("rows predicted for a subset equal those rows of the whole surface", {
  rows <- c(3, 7, 8, 20)
  for (cx in list(ctx, ctx_g)) {
    prep <- prepare_draw(x_o, theta4, 42, cx)
    whole <- predict_rows(seq_len(tz$m), prep, cx)
    expect_lt(max(abs(predict_rows(rows, prep, cx) - whole[rows, ])), 1e-10)
  }
})

test_that("prediction leaves the caller's random-number stream alone", {
  set.seed(99); a <- runif(3)
  set.seed(99); invisible(predict_surface(x_o, theta4, 5, ctx)); b <- runif(3)
  expect_identical(a, b)
  # and a given seed always gives the same surface
  expect_identical(predict_surface(x_o, theta4, 5, ctx), predict_surface(x_o, theta4, 5, ctx))
})
