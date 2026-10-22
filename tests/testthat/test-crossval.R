# Checks of crossval():
#   1. The held-out kriging moments equal the dense model's conditional
#      mean and covariance, with a nugget on the held-out cells too.
#   2. crossval() runs, scores every held-out cell once, gives the same
#      folds for the same seed, and leaves the user's stream untouched.

test_that("held-out moments equal the dense conditional distribution", {
  tz <- make_toy()
  with(tz, {
    sn <- 0.6; th4 <- c(theta, log(sn))
    set.seed(41)
    x_o <- mu_o + as.vector(build_sigma_oo(th4, obs, eig, B, jitter = 0)$L %*% rnorm(n))
    held_cells <- sort(sample(setdiff(seq_len(m * T), obs_cells), 6))
    a_h <- (held_cells - 1) %% m + 1; t_h <- (held_cells - 1) %/% m + 1
    Xh <- design_X(t_h, T)
    mo <- heldout_moments(th4, x_o, obs, eig, B, b0, Xh, a_h, t_h, jitter = 0)
    S <- Sigma_full_fx
    S_oo <- S[obs_cells, obs_cells] + diag(sn^2, n)
    S_ho <- S[held_cells, obs_cells]
    mu_full <- as.vector(X_full %*% b0)
    m_exact <- mu_full[held_cells] + S_ho %*% solve(S_oo, x_o - mu_full[obs_cells])
    C_exact <- S[held_cells, held_cells] + diag(sn^2, 6) - S_ho %*% solve(S_oo, t(S_ho))
    expect_lt(max(abs(mo$mean - m_exact)), 1e-8)
    expect_lt(max(abs(mo$cov - C_exact)), 1e-8)
  })
})

test_that("crossval() scores every cell once, reproducibly", {
  fr <- areal_lattice(5, 5); tf <- time_frame(T = 6)
  sim <- simulate_areal(fr, tf, n_cells = 40, k = 30, seed = 2)
  d <- areal_data(sim$data, fr, tf, y = "positive", k = "tested", area = "id", t = "t", verbose = FALSE)
  ctl <- mcmc_control(n_iter = 40, n_burn = 20, n_store = 10, n_chains = 1, ess_sweeps = 2, theta_moves = 1,
                      progress = FALSE)
  fit <- fit_areal(d, control = ctl)
  set.seed(5); before <- runif(3)
  set.seed(5); cv <- crossval(fit, folds = 4); after <- runif(3)
  expect_identical(after, before)
  expect_setequal(cv$cells$cell_id, d$cells$cell_id)
  expect_equal(nrow(cv$cells), nrow(d$cells))
  expect_equal(nrow(cv$folds), 4)
  expect_true(all(cv$cells$pred_lo <= cv$cells$pred_mean & cv$cells$pred_mean <= cv$cells$pred_hi))
  expect_true(all(is.finite(cv$cells$log_score) & cv$cells$log_score <= 0))
  expect_type(cv$cells$outside, "logical")
  cv2 <- crossval(fit, folds = 4, fold = 2)
  expect_identical(cv2$cells, cv$cells[cv$cells$fold == 2, ], ignore_attr = TRUE)
  expect_output(print(cv), "stareal_cv")
  expect_error(crossval(fit, folds = 1), "folds")
})
