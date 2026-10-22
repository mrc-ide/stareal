# Checks of the chain runner, fit_areal() and diagnose():
#   - settings are validated;
#   - run_chain() returns draws of the right shape, and diagnose() covers
#     every parameter (ported from Africa_areal/tests/test_sampler.R,
#     check 4);
#   - a chain stopped at a checkpoint and resumed is identical to one run
#     straight through, stopped during burn-in (adaptation running) or after
#     it (check 5); a finished fit is read back from its checkpoints, and a
#     checkpoint folder from a different fit is refused;
#   - each stored fixed-effect draw is a draw from beta | x_o, theta at that
#     iteration (whitened, the draws are standard normal);
#   - chains in parallel give the draws of chains run one after another;
#   - the user's random-number stream is untouched;
#   - recovery of known truth: data simulated from the model, fitted, and
#     the 95% intervals checked against the truth.

toy_data <- function(n_cells = 40, seed = 2) {
  fr <- areal_lattice(5, 5); tf <- time_frame(T = 6)
  sim <- simulate_areal(fr, tf, n_cells = n_cells, k = 30, seed = seed)
  areal_data(sim$data, fr, tf, y = "positive", k = "tested", area = "id", t = "t", verbose = FALSE)
}
quick <- function(...) mcmc_control(n_iter = 60, n_burn = 20, n_store = 10, n_chains = 2, ess_sweeps = 2,
                                    theta_moves = 1, theta_thin = 3, parallel = FALSE, progress = FALSE, ...)

test_that("priors and control are validated", {
  expect_error(areal_priors(sigma_sd = -1), "sigma_sd")
  expect_error(areal_priors(rho_s = 2), "rho_s")
  expect_error(mcmc_control(n_iter = 100, n_burn = 100), "n_burn")
  expect_error(mcmc_control(n_iter = 100, n_burn = 50, n_store = 60), "n_store")
  expect_error(mcmc_control(n_chains = 1.5), "n_chains")
  expect_s3_class(mcmc_control(), "stareal_control")
  d <- toy_data()
  expect_error(fit_areal(d, priors = areal_priors(beta_sd = c(1, 2, 3)), control = quick()), "beta_sd")
  expect_error(fit_areal(d$cells), "areal_data")
})

test_that("run_chain returns draws of the right shape and diagnose covers every parameter", {
  d <- toy_data()
  fit <- fit_areal(d, control = quick())
  n <- nrow(d$cells)
  expect_equal(dim(fit$draws$x), c(10, 2, n))
  expect_equal(dim(fit$draws$theta), c(10, 2, 4))
  expect_equal(dim(fit$draws$beta), c(10, 2, 2))
  expect_equal(nrow(fit$trace), 2 * 20)
  expect_false(anyNA(fit$draws$x))
  dg <- suppressWarnings(diagnose(fit))                 # posterior warns on such short chains
  expect_equal(nrow(dg$table), 4 + 2 + n)
  expect_equal(dg$table$cell_id[dg$table$type == "latent"], d$cells$cell_id)
  expect_false(dg$pass)                                  # 60 iterations: far too short
  da <- posterior::as_draws_array(fit)
  expect_equal(posterior::nvariables(da), 4 + 2 + n)
  expect_equal(posterior::ndraws(da), 20)
  expect_s3_class(suppressWarnings(summary(fit)), "summary.stareal_fit")
  expect_output(suppressWarnings(print(fit)), "stareal_fit")
  expect_output(print(dg), "FAIL")
  # no nugget: three hyperparameters
  fit3 <- fit_areal(d, nugget = FALSE, control = quick())
  expect_equal(dim(fit3$draws$theta)[3], 3)
  expect_equal(nrow(suppressWarnings(diagnose(fit3))$table), 3 + 2 + n)
})

test_that("a chain resumed from a checkpoint is identical to one run straight through", {
  d <- toy_data()
  ctl <- quick(checkpoint_every = 7)
  setup <- fit_setup(d$cells, d$frame$eig, d$time$T, areal_priors(), TRUE)
  straight <- run_chain(1, setup, ctl)
  ck <- tempfile(fileext = ".rds")
  for (stop_at in c(13, 41)) {                           # during burn-in, and after it
    unlink(ck)
    expect_null(run_chain(1, setup, ctl, checkpoint = ck, stop_after = stop_at))
    set.seed(stop_at); runif(5)                          # move the stream: the resume must restore it
    resumed <- run_chain(1, setup, ctl, checkpoint = ck)
    expect_identical(resumed$x, straight$x)
    expect_identical(resumed$draw_theta, straight$draw_theta)
    expect_identical(resumed$beta, straight$beta)
    expect_identical(resumed$theta, straight$theta)
  }
  unlink(ck)
})

test_that("fit_areal resumes from, and is read back from, its checkpoint folder", {
  d <- toy_data()
  dir <- tempfile("ck")
  plain <- fit_areal(d, control = quick())
  first <- fit_areal(d, control = quick(checkpoint_dir = dir, checkpoint_every = 7))
  expect_identical(first$draws, plain$draws)
  expect_true(all(file.exists(file.path(dir, c("chain_1.rds", "chain_2.rds", "fingerprint.txt")))))
  again <- suppressMessages(fit_areal(d, control = quick(checkpoint_dir = dir, checkpoint_every = 7)))
  expect_identical(again$draws, plain$draws)
  expect_identical(again$trace, plain$trace)
  # a different fit (other priors) must not resume these checkpoints
  expect_error(fit_areal(d, priors = areal_priors(sigma_sd = 1),
                         control = quick(checkpoint_dir = dir, checkpoint_every = 7)), "different fit")
  unlink(dir, recursive = TRUE)
})

test_that("the stored fixed effects are draws from beta | x_o, theta", {
  d <- toy_data()
  fit <- fit_areal(d, control = mcmc_control(n_iter = 200, n_burn = 50, n_store = 100, n_chains = 2,
                                             ess_sweeps = 2, theta_moves = 1, parallel = FALSE, progress = FALSE))
  st <- fit$setup
  XBX <- st$obs$X %*% st$B %*% t(st$obs$X)
  z <- NULL
  for (c in 1:2) for (i in 1:100) {
    sig <- build_sigma_oo(fit$draws$theta[i, c, ], st$obs, d$frame$eig, st$B)
    bp <- beta_posterior(fit$draws$x[i, c, ], st$obs$X, sig$S - XBX, st$b0, st$B)
    z <- rbind(z, forwardsolve(t(chol(bp$cov)), fit$draws$beta[i, c, ] - bp$mean))
  }
  # 200 whitened draws: column means have sd 0.07, variances about 0.1
  expect_lt(max(abs(colMeans(z))), 0.3)
  expect_true(all(abs(apply(z, 2, var) - 1) < 0.4))
})

test_that("a fit leaves the user's random-number stream untouched", {
  d <- toy_data()
  set.seed(5); before <- runif(3)
  set.seed(5); fit_areal(d, control = quick()); after <- runif(3)
  expect_identical(after, before)
})

test_that("chains in parallel give the draws of chains run one after another", {
  skip_on_cran()
  d <- toy_data()
  seq_fit <- fit_areal(d, control = quick())
  par_fit <- fit_areal(d, control = mcmc_control(n_iter = 60, n_burn = 20, n_store = 10, n_chains = 2,
                                                 ess_sweeps = 2, theta_moves = 1, theta_thin = 3,
                                                 parallel = TRUE, n_cores = 2, progress = FALSE))
  # the same chains; workers run BLAS on one thread, so equal up to rounding
  expect_equal(par_fit$draws, seq_fit$draws, tolerance = 1e-8)
  expect_equal(par_fit$trace, seq_fit$trace, tolerance = 1e-8)
})

test_that("a fit recovers the parameters and latent values it was simulated from", {
  skip_on_cran()
  # 64 units x 10 steps, 200 observed cells; about 15 s on two workers.
  # Seed 1 was the first tried. Over seeds 1-6 the hyperparameter
  # intervals covered the truth 23 times in 24 and the latent coverage
  # was 0.92-0.99; the fixed-effect intervals missed in seeds 3-5, because
  # the posteriors of rho_s and rho_t sit below the truth with T = 10, which
  # narrows them (CLAUDE.md, "Verification").
  fr <- areal_lattice(8, 8); tf <- time_frame(T = 10)
  sim <- simulate_areal(fr, tf, sigma = 1, rho_s = 0.9, rho_t = 0.8, sigma_n = 0.3, beta = c(-1, 1),
                        n_cells = 200, k = 100, seed = 1)
  d <- areal_data(sim$data, fr, tf, y = "positive", k = "tested", area = "id", t = "t", verbose = FALSE)
  fit <- fit_areal(d, control = mcmc_control(n_iter = 2000, n_burn = 700, n_store = 200, n_chains = 2,
                                             ess_sweeps = 5, n_cores = 2, progress = FALSE))
  tab <- diagnose(fit)$table
  par <- tab[tab$type != "latent", ]
  truth <- c(sim$truth$theta, sim$truth$beta)[par$variable]
  expect_true(all(par$q2.5 <= truth & truth <= par$q97.5),
              label = paste("95% intervals cover the truth:", paste(par$variable, collapse = ", ")))
  expect_lt(max(tab$rhat), 1.1)
  # latent cells, logit scale, nugget included
  xt <- sim$truth$x[match(paste(d$cells$id, d$cells$t), paste(sim$data$id, sim$data$t))]
  q <- apply(matrix(fit$draws$x, ncol = nrow(d$cells)), 2, quantile, c(0.025, 0.975))
  coverage <- mean(q[1, ] <= xt & xt <= q[2, ])
  expect_gt(coverage, 0.85)
  expect_lte(coverage, 1)
})
