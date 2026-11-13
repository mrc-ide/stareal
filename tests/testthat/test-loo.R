# Checks of loo() and compare_models():
#   1. The leave-one-out moments of each observed cell, x_i | x_{-i}, theta,
#      equal the held-out kriging moments of crossval() for that cell with
#      the others observed (an independent route: the cell is removed from
#      the observed set and predicted), for every spatial family, with and
#      without the nugget and with a covariate.
#   2. loo() gives a loo object over the observed cells, with the cells
#      table, and its pointwise values are the PSIS estimates from the
#      integrated log-likelihood matrix.
#   3. compare_models() lines up fits of the same data, refuses others,
#      takes cross-validation results and precomputed loo results.

tz <- make_toy()
fr_toy <- areal_frame(adjacency = tz$W, ids = sprintf("u%02d", seq_len(tz$m)), verbose = FALSE)

test_that("leave-one-out moments equal the held-out kriging moments", {
  th <- c(log(0.8), qlogis(0.7), qlogis(0.6), log(0.5))
  for (fam in c("leroux", "bym2", "car")) for (nug in c(TRUE, FALSE)) for (gen in c(FALSE, TRUE)) {
    bas <- spatial_basis(fr_toy, fam)          # (not `eig`: inside with(tz, ...) tz$eig would win)
    th_i <- if (nug) th else th[1:3]
    with(tz, {
      ob <- if (gen) obs_g else obs; bb <- if (gen) b0g else b0; BB <- if (gen) Bg else B
      mu <- as.vector(ob$X %*% bb)
      set.seed(9)
      x_o <- mu + as.vector(build_sigma_oo(th_i, ob, bas, BB, jitter = 0)$L %*% rnorm(n))
      lm <- loo_moments(th_i, x_o, ob, bas, BB, mu, jitter = 0)
      for (i in c(1, 7, n)) {
        rest <- make_obs(a_o[-i], t_o[-i], T, X = ob$X[-i, , drop = FALSE])
        ho <- heldout_moments(th_i, x_o[-i], rest, bas, BB, bb, ob$X[i, , drop = FALSE], a_o[i], t_o[i], jitter = 0)
        expect_lt(abs(lm$mean[i] - ho$mean), 1e-8)
        expect_lt(abs(lm$var[i] - ho$cov[1, 1]), 1e-8)
      }
    })
  }
})

ctl <- list(n_iter = 120, n_burn = 60, n_store = 20, n_chains = 2, n_cores = 1, verbose = FALSE)
sim <- simulate_areal(fr_toy, time_frame(n_steps = 5), n_cells = 30, seed = 3)
f1 <- fit_areal(sim, control = ctl)
f2 <- fit_areal(sim, spatial = "bym2", control = ctl)

test_that("loo() is PSIS-LOO on the integrated log-likelihood", {
  l <- loo(f1)
  expect_s3_class(l, c("stareal_loo", "psis_loo", "loo"))
  expect_identical(nrow(l$cells), nrow(f1$data$cells))
  expect_identical(l$cells$cell_id, f1$data$cells$cell_id)
  ll <- loo_loglik(f1)
  expect_identical(dim(ll$ll), c(40L, nrow(f1$data$cells)))
  ref <- suppressWarnings(loo::loo(ll$ll, r_eff = suppressWarnings(
    loo::relative_eff(exp(ll$ll - matrix(apply(ll$ll, 2, max), nrow(ll$ll), ncol(ll$ll), byrow = TRUE)),
                      chain_id = ll$chain))))
  expect_equal(l$estimates, ref$estimates)
  expect_output(print(l), "integrated PSIS-LOO over 30 observed cells, 40 draws")
  expect_identical(attr(loo(f1, n_draws = 20), "dims"), c(20L, 30L))
})

test_that("compare_models() lines up fits of the same data", {
  cm <- compare_models(leroux = f1, bym2 = f2)
  expect_identical(sort(cm$model), c("bym2", "leroux"))
  expect_equal(cm$elpd_diff[1], 0)
  expect_true(cm$elpd_diff[2] <= 0)
  expect_equal(cm$elpd_loo[2] - cm$elpd_loo[1], cm$elpd_diff[2])
  expect_identical(names(attr(cm, "loo")), c("leroux", "bym2"))
  # precomputed loo results are used as given
  l1 <- loo(f1)
  cm2 <- compare_models(leroux = f1, bym2 = f2, loo = list(leroux = l1))
  expect_equal(cm2$elpd_loo, cm$elpd_loo)
  # with cross-validation
  cvc <- list(n_iter = 60, n_burn = 30, n_store = 10, n_chains = 1, verbose = FALSE)
  cv <- list(leroux = crossval(f1, folds = 2, control = cvc), bym2 = crossval(f2, folds = 2, control = cvc))
  cm3 <- compare_models(leroux = f1, bym2 = f2, cv = cv)
  expect_true(all(c("cv_log_score", "cv_diff", "cv_se_diff", "cv_rmse") %in% names(cm3)))
  # other data are refused
  sim2 <- simulate_areal(fr_toy, time_frame(n_steps = 5), n_cells = 30, seed = 4)
  f3 <- fit_areal(sim2, control = ctl)
  expect_error(compare_models(f1, f3), "same observed cells")
  expect_error(compare_models(f1), "at least two")
  expect_error(compare_models(a = f1, b = f2, cv = list(a = cv$leroux)), "for every fit")
})
