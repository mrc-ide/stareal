# Checks of the spatial priors and of fixed hyperparameters:
#   1. Each family's basis gives the covariance it is defined by, against a
#      dense construction made independently: BYM2 from the generalised
#      inverse (R + J/m)^{-1} - J/m, the proper CAR by inverting
#      tau (D - rho_s W), the ICAR as R^+ / tau, iid as I / tau.
#   2. For each family, Sigma_oo equals the dense model's observed block,
#      and Matheron's rule with zero noise gives the exact kriging mean (the
#      proper CAR's basis is not orthonormal; nothing may assume it is).
#   3. Fixed hyperparameters: theta_full() puts them back in place; the
#      log prior covers the free ones only; the default model's prior and
#      basis are untouched.
#   4. model_spec() refuses what it should; fits with fixed
#      hyperparameters leave them out of the draws, the trace and summary();
#      the checkpoint fingerprint changes with the model, but not for the
#      default one.

tz <- make_toy()
toy_frame <- function() areal_frame(adjacency = tz$W, ids = sprintf("u%02d", seq_len(tz$m)), verbose = FALSE)
fr_toy <- toy_frame()
dense_sigma_s <- function(eig, p) eig$vectors %*% diag(spatial_g(eig, p)) %*% t(eig$vectors)

test_that("each spatial family's basis gives its defining covariance", {
  m <- tz$m; R <- tz$R_lap; D <- diag(rowSums(tz$W)); J <- matrix(1 / m, m, m)
  R_plus <- solve(R + J) - J                                # the generalised inverse, connected graph
  p <- list(sigma = 0.8, tau = 1 / 0.64, rho_s = 0.7)
  # BYM2: scaled ICAR + iid
  b <- spatial_basis(fr_toy, "bym2")
  s <- exp(mean(log(diag(R_plus))))
  expect_lt(max(abs(dense_sigma_s(b, p) - 0.64 * (0.7 * R_plus / s + 0.3 * diag(m)))), 1e-10)
  # proper CAR
  cr <- spatial_basis(fr_toy, "car")
  expect_lt(max(abs(dense_sigma_s(cr, p) - solve(p$tau * (D - 0.7 * tz$W)))), 1e-10)
  # ICAR (Leroux at rho_s = 1) and iid (rho_s = 0)
  lr <- spatial_basis(fr_toy, "leroux")
  expect_lt(max(abs(dense_sigma_s(lr, modifyList(p, list(rho_s = 1))) - R_plus / p$tau)), 1e-10)
  expect_lt(max(abs(dense_sigma_s(lr, modifyList(p, list(rho_s = 0))) - diag(m) / p$tau)), 1e-10)
  # the default model's basis is the frame's eigendecomposition, untouched
  expect_identical(spatial_basis(fr_toy), fr_toy$eig)
})

test_that("Sigma_oo and Matheron's rule are exact for every family", {
  local_mocked_bindings(draw_noise = zero_noise)
  th <- c(log(0.8), qlogis(0.7), qlogis(0.6), log(0.5))
  pp <- theta_to_par(th)
  # (names other than tz's: inside with(tz, ...) its eig and p would win)
  for (fam in c("leroux", "bym2", "car")) {
    bas <- spatial_basis(fr_toy, fam)
    with(tz, {
      S_full <- kronecker(ar1_cov(T, pp$rho_t), dense_sigma_s(bas, pp)) + X_full %*% B %*% t(X_full)
      S_oo <- S_full[obs_cells, obs_cells] + diag(pp$sigma_n^2, n)
      sig <- build_sigma_oo(th, obs, bas, B, jitter = 0)
      expect_lt(max(abs(sig$S - S_oo)), 1e-10)
      set.seed(5)
      x_o <- mu_o + as.vector(sig$L %*% rnorm(n))
      mu_full <- as.vector(X_full %*% b0)
      m_exact <- mu_full + S_full[, obs_cells] %*% solve(S_oo, x_o - mu_full[obs_cells])
      ctx <- list(eig = bas, obs = obs, b0 = b0, B = B, T = T, jitter = 0)
      expect_lt(max(abs(as.vector(predict_surface(x_o, th, 1, ctx)) - m_exact)), 1e-8)
    })
  }
})

test_that("fixed hyperparameters are put back in place, and only free ones have a prior", {
  fx <- c(logit_rho_s = Inf, logit_rho_t = qlogis(0.3))
  expect_equal(theta_full(c(0.1, -0.5), fx), c(0.1, Inf, qlogis(0.3), -0.5))
  expect_identical(theta_full(c(1, 2, 3)), c(1, 2, 3))
  p <- theta_to_par(c(log(1.3), log(0.4)), fx)
  expect_equal(c(p$sigma, p$rho_s, p$rho_t, p$sigma_n), c(1.3, 1, 0.3, 0.4))
  expect_identical(theta_names(TRUE, fx), c("log_sigma", "log_sigma_n"))
  lp <- make_log_prior(tz$priors, fx)
  direct <- log(2 * dnorm(1.3, 0, 2)) + log(1.3) + log(2 * dnorm(0.4, 0, 1)) + log(0.4)
  expect_lt(abs(lp(c(log(1.3), log(0.4))) - direct), 1e-12)
  # natural-scale columns by name
  nat <- theta_natural(matrix(c(log(2), log(0.5)), 1, dimnames = list(NULL, c("log_sigma", "log_sigma_n"))))
  expect_equal(nat, matrix(c(2, 0.5), 1, dimnames = list(NULL, c("sigma", "sigma_n"))))
})

test_that("model_spec() checks the spatial prior and fixed values", {
  expect_identical(model_spec()$fixed, NULL)
  expect_equal(model_spec("icar")$fixed, c(rho_s = 1))
  expect_equal(model_spec("iid", list(rho_t = 0))$fixed, c(rho_s = 0, rho_t = 0))
  expect_equal(model_spec(fixed = c(rho_t = 0.5, sigma = 2))$fixed, c(sigma = 2, rho_t = 0.5))   # ordered
  expect_error(model_spec("besag"), "must be one of")
  expect_error(model_spec("icar", list(rho_s = 0.5)), "fixes rho_s already")
  expect_error(model_spec(fixed = list(phi = 1)), "can hold")
  expect_error(model_spec(fixed = list(rho_t = 1)), "below 1")
  expect_error(model_spec("car", list(rho_s = 1)), "below 1 for the proper CAR")
  expect_error(model_spec(fixed = list(sigma_n = 0.2), nugget = FALSE), "needs nugget = TRUE")
  expect_error(model_spec(fixed = list(sigma = 1, rho_s = 0.5, rho_t = 0.5), nugget = FALSE), "leave at least one free")
  expect_error(model_spec(fixed = list(0.5)), "name each")
  # graphs the priors cannot use
  W2 <- matrix(0, 4, 4); W2[1, 2] <- W2[2, 1] <- W2[3, 4] <- W2[4, 3] <- 1
  fr2 <- suppressWarnings(areal_frame(adjacency = W2, ids = letters[1:4], verbose = FALSE))
  expect_error(model_spec("bym2", frame = fr2), "connected graph")
  expect_error(model_spec("icar", frame = fr2), "connected graph")
  expect_silent(model_spec("car", frame = fr2))
  W3 <- W2; W3[3, 4] <- W3[4, 3] <- 0
  fr3 <- suppressWarnings(areal_frame(adjacency = W3, ids = letters[1:4], verbose = FALSE))
  expect_error(model_spec("car", frame = fr3), "every area to have a neighbour")
})

test_that("fits with other priors and fixed hyperparameters run and report only the free ones", {
  sim <- simulate_areal(fr_toy, time_frame(n_steps = 5), n_cells = 30, seed = 2, spatial = "bym2")
  expect_identical(sim$truth$spatial, "bym2")
  ctl <- list(n_iter = 60, n_burn = 20, n_store = 10, n_chains = 2, n_cores = 1, verbose = FALSE)
  f <- fit_areal(sim, spatial = "car", fixed = list(rho_t = 0.5), control = ctl)
  expect_identical(dimnames(f$draws$theta)[[3]], c("log_sigma", "logit_rho_s", "log_sigma_n"))
  expect_false("logit_rho_t" %in% names(f$trace))
  s <- summary(f)
  expect_identical(s$table$variable[s$table$type == "theta"], c("sigma", "rho_s", "sigma_n"))
  expect_true(all(c("sigma", "rho_s", "sigma_n") %in% posterior::variables(posterior::as_draws(f, latent = FALSE))))
  expect_output(print(f), "spatial prior: car, rho_t = 0.5")
  expect_identical(f$model$family, "car")
  expect_identical(f$version$stareal, as.character(utils::packageVersion("stareal")))
  fi <- fit_areal(sim, spatial = "icar", control = ctl)
  expect_identical(dimnames(fi$draws$theta)[[3]], c("log_sigma", "logit_rho_t", "log_sigma_n"))
  pr <- predict(fi, n_draws = 5, verbose = FALSE)
  expect_true(all(is.finite(pr$area$mean)))
})

test_that("the checkpoint fingerprint covers the model, and not for the default one", {
  sim <- simulate_areal(fr_toy, time_frame(n_steps = 5), n_cells = 30, seed = 2)
  cells <- sim$cells; X <- design_X(cells$t, 5)
  pri <- complete_settings(areal_priors(), areal_priors, "priors")
  pri$beta_mean <- c(0, 0); pri$beta_sd <- c(3, 3)
  ctl <- mcmc_control(verbose = FALSE)
  st0 <- fit_setup(cells, fr_toy$eig, 5, pri, TRUE, X)
  expect_null(st0$eig$family)                                 # the default: no model key
  st1 <- fit_setup(cells, spatial_basis(fr_toy, "bym2"), 5, pri, TRUE, X)
  d <- tempfile("ck"); dir.create(d); on.exit(unlink(d, recursive = TRUE))
  fp0 <- checkpoint_guard(file.path(d, "a"), st0, ctl)
  fp1 <- checkpoint_guard(file.path(d, "b"), st1, ctl)
  expect_false(identical(fp0, fp1))
})
