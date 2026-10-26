# Checks of the chain runner, fit_areal(), summary() and diagnose():
#   - settings are validated, and plain lists are completed;
#   - run_chain() returns draws of the right shape, and diagnose() covers
#     every parameter;
#   - chains start from dispersed points;
#   - a chain stopped at a checkpoint and resumed is identical to one run
#     straight through, stopped during burn-in (adaptation and pseudo-data
#     refits running) or after it, and resumed with a different
#     checkpoint_every; a finished fit is read back from its checkpoints;
#     a checkpoint folder from a different fit (other priors, other graph)
#     is refused;
#   - a failing chain says which, with the resume hint;
#   - each stored fixed-effect draw is a draw from beta | x_o, theta at that
#     iteration (whitened, the draws are standard normal);
#   - chains in parallel give the draws of chains run one after another;
#   - the user's random-number stream is untouched, and the draws do not
#     depend on the session's RNGkind;
#   - recovery of known truth: data simulated from the model, fitted, and
#     every posterior mean within 3 posterior sd of the truth.

toy_data <- function(n_cells = 40, seed = 2) {
  simulate_areal(areal_lattice(5, 5), time_frame(n_steps = 6), n_cells = n_cells, k = 30, seed = seed)
}
quick <- function(...) utils::modifyList(list(n_iter = 60, n_burn = 20, n_store = 10, n_chains = 2, ess_sweeps = 2,
                                              theta_moves = 1, theta_thin = 3, n_cores = 1, verbose = FALSE),
                                         list(...))

test_that("priors and control are validated, and plain lists are completed", {
  expect_error(areal_priors(sigma_sd = -1), "sigma_sd")
  expect_error(areal_priors(rho_s = 2), "rho_s")
  expect_error(areal_priors(beta_sd = Inf), "beta_sd")
  expect_error(mcmc_control(n_iter = 100, n_burn = 100), "n_burn")
  expect_error(mcmc_control(n_iter = 100, n_burn = 50, n_store = 60), "n_store")
  expect_error(mcmc_control(n_chains = 1.5), "n_chains")
  expect_error(mcmc_control(n_iter = 100, n_burn = 70, n_store = 5, theta_thin = 60), "theta_thin")
  expect_type(mcmc_control(), "list")
  expect_equal(mcmc_control(seed = 20261012)$seed, 20261012)        # date-style seeds are fine
  d <- toy_data()
  expect_error(fit_areal(d, priors = list(sigma = 1), control = quick()), "unknown setting")
  expect_error(fit_areal(d, control = list(n_iters = 10)), "unknown setting")
  expect_error(fit_areal(d, priors = list(beta_sd = c(1, 2, 3)), control = quick()), "beta_sd")
  expect_error(fit_areal(d, priors = list(beta_sd = c(slope = 1)), control = quick()), "unknown design column")
  expect_error(fit_areal(d$cells), "areal_data")
  # named priors are matched by design column
  f <- fit_areal(d, priors = list(beta_sd = c(trend = 0.5)), control = quick(n_chains = 1))
  expect_equal(f$priors$beta_sd, c(intercept = 3, trend = 0.5))
})

test_that("a fit has draws of the right shape, and diagnose covers every parameter", {
  d <- toy_data()
  fit <- fit_areal(d, control = quick())
  n <- nrow(d$cells)
  expect_equal(dim(fit$draws$x), c(10, 2, n))
  expect_equal(dim(fit$draws$theta), c(10, 2, 4))
  expect_equal(dim(fit$draws$beta), c(10, 2, 2))
  expect_equal(nrow(fit$trace), 2 * 20)
  expect_false(anyNA(fit$draws$x))
  expect_equal(fit$stored_iter[10], 60)                    # the last iteration is stored
  dg <- suppressWarnings(diagnose(fit))                    # posterior warns on such short chains
  expect_s3_class(dg, "summary.stareal_fit")
  expect_equal(nrow(dg$table), 4 + 2 + n)
  expect_equal(dg$table$cell_id[dg$table$type == "latent"], d$cells$cell_id)
  expect_false(dg$pass)                                    # 60 iterations: far too short
  expect_true("beta bulk ESS (min)" %in% dg$checks$check)
  s <- suppressWarnings(summary(fit))
  expect_equal(nrow(s$table), 4 + 2)
  expect_identical(s$table, dg$table[dg$table$type != "latent", ], ignore_attr = TRUE)
  expect_null(s$pass)
  expect_equal(posterior::nvariables(posterior::as_draws_array(fit)), 4 + 2 + n)
  expect_equal(posterior::ndraws(posterior::as_draws_df(fit, latent = FALSE)), 20)
  expect_s3_class(posterior::as_draws_rvars(fit, latent = FALSE), "draws_rvars")
  expect_output(suppressWarnings(print(fit)), "stareal_fit")
  expect_output(print(dg), "FAIL")
  expect_output(print(s), "use diagnose")
  expect_true(is.call(fit$call))
  expect_null(fit$setup)
  # no nugget: three hyperparameters
  fit3 <- fit_areal(d, nugget = FALSE, control = quick())
  expect_equal(dim(fit3$draws$theta)[3], 3)
  expect_equal(nrow(suppressWarnings(diagnose(fit3))$table), 3 + 2 + n)
  pdf(NULL); on.exit(grDevices::dev.off())
  expect_silent(plot(fit3, col = c("red", "blue"), main = "x"))
  # plot() has a panel for every parameter of summary(), fixed effects included
  panels <- trace_panels(fit3)
  expect_identical(names(panels), suppressWarnings(summary(fit3))$table$variable)
  expect_equal(panels$intercept$y[, 1], fit3$draws$beta[, 1, "intercept"])
  expect_length(panels$intercept$x, dim(fit3$draws$beta)[1])
})

test_that("chains start from dispersed points", {
  fit <- fit_areal(toy_data(), control = quick(n_chains = 4, theta_thin = 1))
  first <- fit$trace[fit$trace$iteration == 1, ]
  expect_gt(stats::sd(first$logit_rho_s), 0.2)
  expect_gt(stats::sd(first$log_sigma), 0.1)
})

test_that("a chain resumed from a checkpoint is identical to one run straight through", {
  d <- toy_data()
  ctl <- mcmc_control(n_iter = 60, n_burn = 20, n_store = 10, n_chains = 1, ess_sweeps = 2, theta_moves = 1,
                      theta_thin = 3, n_cores = 1, verbose = FALSE, checkpoint_every = 7)
  setup <- fit_setup_of(fit_areal(d, control = quick(n_chains = 1)))
  straight <- run_chain(1, setup, ctl, seed = 17)
  ck <- tempfile(fileext = ".rds")
  for (stop_at in c(8, 13, 41)) {                          # before, between and after the pseudo-data refits
    unlink(ck)
    expect_null(run_chain(1, setup, ctl, seed = 17, checkpoint = ck, stop_after = stop_at))
    set.seed(stop_at); runif(5)                            # move the stream: the resume must restore it
    ctl2 <- ctl; ctl2$checkpoint_every <- 5                # does not change the draws
    resumed <- run_chain(1, setup, ctl2, seed = 17, checkpoint = ck)
    expect_identical(resumed$x, straight$x)
    expect_identical(resumed$draw_theta, straight$draw_theta)
    expect_identical(resumed$beta, straight$beta)
    expect_identical(resumed$theta, straight$theta)
  }
  unlink(ck)
})

test_that("fit_areal resumes from, and is read back from, its checkpoint folder; other fits are refused", {
  d <- toy_data()
  dir <- tempfile("ck")
  plain <- fit_areal(d, control = quick())
  first <- fit_areal(d, control = quick(checkpoint_every = 7), checkpoint_dir = dir)
  expect_identical(first$draws, plain$draws)
  expect_true(all(file.exists(file.path(dir, c("chain_1.rds", "chain_2.rds", "fingerprint.txt")))))
  again <- fit_areal(d, control = quick(checkpoint_every = 5), checkpoint_dir = dir)
  expect_identical(again$draws, plain$draws)
  expect_identical(again$trace, plain$trace)
  # a different fit (other priors) must not resume these checkpoints
  expect_error(fit_areal(d, priors = list(sigma_sd = 1), control = quick(), checkpoint_dir = dir), "different fit")
  # nor a fit with the same cells on another graph
  dq <- simulate_areal(areal_lattice(5, 5, queen = FALSE), time_frame(n_steps = 6), n_cells = 40, k = 30, seed = 2)
  dq$cells <- d$cells
  expect_error(fit_areal(dq, control = quick(), checkpoint_dir = dir), "different fit")
  unlink(dir, recursive = TRUE)
})

test_that("a failing chain says which, with the resume hint", {
  setup <- fit_setup_of(fit_areal(toy_data(), control = quick(n_chains = 1)))
  setup$lik$loglik <- function(x) stop("boom")
  ctl <- mcmc_control(n_iter = 60, n_burn = 20, n_store = 10, n_chains = 2, n_cores = 1, verbose = FALSE)
  expect_error(run_chains(setup, ctl), "chain 1 failed: boom")
  expect_error(run_chains(setup, ctl, checkpoint_dir = tempfile("ck")), "resume from the last checkpoints")
  skip_on_cran()
  ctl$n_cores <- 2
  expect_error(run_chains(setup, ctl), "chain [12] failed: boom")
})

test_that("the stored fixed effects are draws from beta | x_o, theta", {
  d <- toy_data()
  fit <- fit_areal(d, control = list(n_iter = 200, n_burn = 50, n_store = 100, n_chains = 2,
                                     ess_sweeps = 2, theta_moves = 1, n_cores = 1, verbose = FALSE))
  st <- fit_setup_of(fit)
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

test_that("a fit leaves the user's stream untouched and does not depend on RNGkind", {
  d <- toy_data()
  set.seed(5); before <- runif(3)
  set.seed(5); f1 <- fit_areal(d, control = quick()); after <- runif(3)
  expect_identical(after, before)
  old <- RNGkind("L'Ecuyer-CMRG"); on.exit(do.call(RNGkind, as.list(old)))
  f2 <- fit_areal(d, control = quick())
  expect_identical(f2$draws, f1$draws)
  expect_identical(RNGkind()[1], "L'Ecuyer-CMRG")
})

test_that("chains in parallel give the draws of chains run one after another", {
  skip_on_cran()
  d <- toy_data()
  seq_fit <- fit_areal(d, control = quick())
  par_fit <- fit_areal(d, control = quick(n_cores = 2))
  # the same chains; workers run BLAS on one thread, so equal up to rounding
  expect_equal(par_fit$draws, seq_fit$draws, tolerance = 1e-8)
  expect_equal(par_fit$trace, seq_fit$trace, tolerance = 1e-8)
})

test_that("a fit recovers the parameters and latent values it was simulated from", {
  skip_on_cran()
  # 64 units x 10 steps, 200 observed cells; about 15 s on two workers. A
  # bias check rather than a coverage check, so that it does not depend on
  # one seed: a sampler fault moves posterior means by more than 3 sd.
  sim <- simulate_areal(areal_lattice(8, 8), time_frame(n_steps = 10), sigma = 1, rho_s = 0.9, rho_t = 0.8,
                        sigma_n = 0.3, beta = c(-1, 1), n_cells = 200, k = 100, seed = 1)
  fit <- fit_areal(sim, control = list(n_iter = 2000, n_burn = 700, n_store = 200, n_chains = 2,
                                       ess_sweeps = 5, n_cores = 2, verbose = FALSE))
  tab <- diagnose(fit)$table
  par <- tab[tab$type != "latent", ]
  truth <- c(sim$truth$theta, sim$truth$beta)[par$variable]
  expect_true(all(abs(par$mean - truth) < 3 * par$sd),
              label = paste("posterior means within 3 sd of the truth:", paste(par$variable, collapse = ", ")))
  expect_lt(max(tab$rhat), 1.1)
  # latent cells, logit scale, nugget included (truth$x is in cell order)
  q <- apply(matrix(fit$draws$x, ncol = nrow(sim$cells)), 2, quantile, c(0.025, 0.975))
  coverage <- mean(q[1, ] <= sim$truth$x & sim$truth$x <= q[2, ])
  expect_gt(coverage, 0.85)
})
