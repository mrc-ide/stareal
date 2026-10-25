# Checks of crossval():
#   1. The held-out kriging moments equal the dense model's conditional
#      mean and covariance, with a nugget on the held-out cells too.
#   2. The quadrature log predictive density equals a direct numerical
#      integral, also when k is large and the likelihood is narrow, and
#      when the observation is far from the predictive mean.
#   3. crossval() runs, scores every held-out cell once (one-cell folds
#      included), gives the same folds for the same seed, and leaves the
#      user's stream untouched; the randomised PIT is uniform under the
#      model.

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

test_that("the quadrature log predictive density equals a direct integral", {
  y <- c(0, 3, 400, 7); k <- c(20, 10, 1000, 7)
  mu <- matrix(c(-2, 0.3, -0.4, 1.5, -1, 0, -0.5, 3), 2, byrow = TRUE)
  s2 <- matrix(c(0.5, 1, 0.3, 2, 0.2, 0.05, 0.6, 0.4), 2, byrow = TRUE)
  lp <- log_binom_gauss(y, k, mu, s2)
  direct <- matrix(NA, 2, 4)
  for (d in 1:2) for (j in 1:4)
    direct[d, j] <- log(stats::integrate(function(x) dbinom(y[j], k[j], plogis(x)) * dnorm(x, mu[d, j], sqrt(s2[d, j])),
                                         -Inf, Inf, rel.tol = 1e-10)$value)
  expect_lt(max(abs(lp - direct)), 1e-6)
  # observations far from the predictive mean, where undamped Newton steps
  # to the mode overshot (log scores of -10^4 for a truth near -6); the
  # reference is a fine grid on the log scale
  grid_lp <- function(y, k, mu, s2) {
    x <- seq(mu - 40, mu + 40, by = 2e-4)
    l <- dbinom(y, k, plogis(x), log = TRUE) + dnorm(x, mu, sqrt(s2), log = TRUE)
    max(l) + log(sum(exp(l - max(l))) * 2e-4)
  }
  far <- data.frame(y = c(31, 31, 31, 485, 26, 0, 1729), k = c(195, 195, 195, 599, 85, 1729, 1729),
                    mu = c(1.8, 2.5, 4, -1, 2.5, 3, 2), s2 = c(4, 25, 0.5, 3, 4, 0.1, 0.1))
  lp_far <- log_binom_gauss(far$y, far$k, matrix(far$mu, 1), matrix(far$s2, 1))[1, ]
  expect_lt(max(abs(lp_far - mapply(grid_lp, far$y, far$k, far$mu, far$s2))), 1e-6)
})

test_that("crossval() scores every cell once, reproducibly", {
  sim <- simulate_areal(areal_lattice(5, 5), time_frame(n_steps = 6), n_cells = 40, k = 30, seed = 2)
  ctl <- list(n_iter = 40, n_burn = 20, n_store = 10, n_chains = 1, ess_sweeps = 2, theta_moves = 1, verbose = FALSE)
  fit <- fit_areal(sim, control = ctl)
  set.seed(5); before <- runif(3)
  set.seed(5); cv <- crossval(fit, folds = 4); after <- runif(3)
  expect_identical(after, before)
  expect_setequal(cv$cells$cell_id, sim$cells$cell_id)
  expect_equal(nrow(cv$cells), nrow(sim$cells))
  expect_equal(nrow(cv$folds), 4)
  expect_true(all(cv$cells$pred_lo <= cv$cells$pred_mean & cv$cells$pred_mean <= cv$cells$pred_hi))
  expect_true(all(is.finite(cv$cells$log_score) & cv$cells$log_score <= 0))
  expect_true(all(cv$cells$pit >= 0 & cv$cells$pit <= 1))
  expect_type(cv$cells$outside, "logical")
  cv2 <- crossval(fit, folds = 4, fold = 2)
  expect_identical(cv2$cells, cv$cells[cv$cells$fold == 2, ], ignore_attr = TRUE)
  expect_output(print(cv), "stareal_cv")
  expect_error(crossval(fit, folds = 1), "folds")
  # one-cell folds (leave-one-out over a few cells)
  loo <- crossval(fit, folds = nrow(sim$cells), fold = 1:3)
  expect_equal(nrow(loo$cells), 3)
  expect_true(all(loo$cells$pred_lo < loo$cells$pred_hi))
  # checkpoints go to a subfolder named by the split seed
  dir <- tempfile("cvck")
  crossval(fit, folds = 4, fold = 1, checkpoint_dir = dir, seed = 9)
  expect_true(dir.exists(file.path(dir, "cv_seed9_fold_1_of_4")))
  unlink(dir, recursive = TRUE)
})

test_that("the randomised PIT is uniform when the predictive distribution is right", {
  # y drawn from the predictive distribution itself: PIT ~ U(0, 1), even for small k
  set.seed(3)
  h <- 4000; D <- 50
  mu <- matrix(rnorm(D * h, -1, 0.2), D, h); v <- matrix(0.1, D, h)
  ph <- list(mean = mu, var = v, draw = mu + matrix(rnorm(D * h, 0, sqrt(0.1)), D, h))
  k <- sample(1:5, h, replace = TRUE)
  y <- rbinom(h, k, plogis(ph$draw[cbind(sample.int(D, h, replace = TRUE), seq_len(h))]))
  held <- data.frame(cell_id = seq_len(h), id = "a", t = 1L, time = "1", y = y, k = k)
  sc <- score_heldout(ph, held)
  expect_lt(abs(mean(sc$outside) - 0.05), 0.015)
  expect_lt(abs(mean(sc$pit) - 0.5), 0.02)
})

test_that("folds can be given as labels, and compare_cv compares models on the same folds", {
  fr <- areal_lattice(5, 5, blocks = c(5, 1))               # five column bands, b1_1 to b1_5
  sim <- simulate_areal(fr, time_frame(n_steps = 6), n_cells = 40, k = 30, seed = 2)
  ctl <- list(n_iter = 40, n_burn = 20, n_store = 10, n_chains = 1, ess_sweeps = 2, theta_moves = 1, verbose = FALSE)
  fit <- fit_areal(sim, control = ctl)
  band <- fr$groups$block[sim$cells$area]                    # spatial folds: hold out a band of areas
  cv_b <- crossval(fit, folds = band, fold = c("b1_2", "b1_4"))
  expect_setequal(unique(cv_b$cells$fold), c("b1_2", "b1_4"))
  expect_equal(nrow(cv_b$cells), sum(band %in% c("b1_2", "b1_4")))
  expect_equal(cv_b$n_folds, length(unique(band)))
  expect_error(crossval(fit, folds = band[-1]), "fold label for each")
  expect_error(crossval(fit, folds = band, fold = "zz"), "among the fold labels")
  # model comparison on the same folds
  fit_nn <- fit_areal(sim, nugget = FALSE, control = ctl)
  cv1 <- crossval(fit, folds = 4, seed = 3); cv2 <- crossval(fit_nn, folds = 4, seed = 3)
  cmp <- compare_cv(nugget = cv1, no_nugget = cv2)
  expect_identical(sort(cmp$model), c("no_nugget", "nugget"))
  expect_equal(cmp$diff[1], 0)
  expect_true(all(cmp$se_diff >= 0))
  expect_equal(cmp$mean_log_score, sort(c(cv1$overall$mean_log_score, cv2$overall$mean_log_score), decreasing = TRUE))
  expect_error(compare_cv(cv1, crossval(fit_nn, folds = 4, seed = 4)), "same folds")
})
