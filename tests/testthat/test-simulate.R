# Checks of simulate.R: it returns data ready to fit, with the truth in
# cell order; reproducibility; the caller's stream untouched; validation;
# and the simulated field's covariance against the dense Kronecker
# covariance of the model.

fr <- areal_lattice(3, 3, queen = FALSE)
ti <- time_frame(n_steps = 4)

test_that("a simulation is data ready to fit, with the truth in cell order", {
  sim <- simulate_areal(fr, ti, n_cells = 20, k = 30, seed = 1)
  expect_s3_class(sim, "stareal_data")
  expect_equal(nrow(sim$cells), 20)
  expect_equal(nrow(sim$surveys), 20)
  expect_true(all(sim$cells$y <= sim$cells$k))
  expect_equal(dim(sim$truth$logit), c(9, 4))
  expect_length(sim$truth$x, 20)
  # truth$x lines up with the cells: the cell's logit plus its nugget
  nug <- sim$truth$x - sim$truth$logit[cbind(sim$cells$area, sim$cells$t)]
  expect_lt(abs(mean(nug)), 0.5)
  expect_identical(names(sim$truth$beta), c("intercept", "trend"))
  expect_output(print(sim), "20 of 20 surveys used.*\\(simulated\\)")
  expect_output(print(sim), "truth: sigma 1")
  sim_d <- simulate_areal(fr, time_frame(2001, 2004, "year"), n_cells = 5, seed = 1)
  expect_s3_class(sim_d$surveys$date, "Date")
})

test_that("inputs are validated", {
  expect_error(simulate_areal(fr, ti, n_cells = 37), "exceeds")
  expect_error(simulate_areal(fr, ti, n_cells = 5, k = 0), "whole numbers")
  expect_error(simulate_areal(fr, ti, sigma = c(1, 2)), "`sigma` must be one number")
  expect_error(simulate_areal(fr, ti, rho_t = 1), "`rho_t` must be one number")
  expect_error(simulate_areal(fr, ti, n_cells = 5, beta = 1), "needs 2 values")
  # a named beta is matched by name, not position
  a <- simulate_areal(fr, ti, n_cells = 5, beta = c(trend = 1, intercept = -2), seed = 4)
  b <- simulate_areal(fr, ti, n_cells = 5, beta = c(-2, 1), seed = 4)
  expect_identical(a$cells, b$cells)
  expect_error(simulate_areal(fr, ti, n_cells = 5, beta = c(intercept = 1, slope = 2)), "named `beta`")
})

test_that("a seed reproduces the simulation and leaves the caller's stream alone", {
  set.seed(5); a <- runif(2)
  set.seed(5); s1 <- simulate_areal(fr, ti, n_cells = 10, seed = 3); b <- runif(2)
  expect_identical(a, b)
  expect_identical(s1$cells, simulate_areal(fr, ti, n_cells = 10, seed = 3)$cells)
  # and does not depend on the session's generator
  old <- RNGkind("L'Ecuyer-CMRG"); on.exit(do.call(RNGkind, as.list(old)))
  expect_identical(s1$cells, simulate_areal(fr, ti, n_cells = 10, seed = 3)$cells)
  expect_identical(RNGkind()[1], "L'Ecuyer-CMRG")
})

test_that("the simulated field has the model's covariance", {
  skip_on_cran()
  sigma <- 0.8; rho_s <- 0.9; rho_t <- 0.7
  draws <- vapply(1:4000, function(i)
    as.vector(simulate_areal(fr, ti, sigma = sigma, rho_s = rho_s, rho_t = rho_t, beta = c(0, 0),
                             n_cells = 1, seed = i)$truth$logit), numeric(36))
  R <- as.matrix(Matrix::Diagonal(x = Matrix::rowSums(fr$W)) - fr$W)
  Ss <- solve((rho_s * R + (1 - rho_s) * diag(9)) / sigma^2)
  C <- kronecker(ar1_cov(4, rho_t), Ss)            # cell order (t - 1) * m + a, as as.vector(logit)
  expect_lt(max(abs(stats::cov(t(draws)) - C)), 0.12 * max(diag(C)))
  expect_lt(max(abs(rowMeans(draws))), 0.1 * sqrt(max(diag(C))))
})

test_that("a simulation can copy a survey design", {
  d0 <- simulate_areal(fr, ti, n_cells = 12, k = 1:12, seed = 5)
  s <- simulate_areal(fr, ti, cells = d0$cells, seed = 6)
  expect_identical(s$cells$cell_id, d0$cells$cell_id)
  expect_equal(s$cells$k, d0$cells$k)
  lab <- data.frame(id = c("r1c1", "r2c2"), time = c("2", "4"), k = c(10, 20))
  s2 <- simulate_areal(fr, ti, cells = lab, seed = 1)
  expect_equal(s2$cells$t, c(2L, 4L)); expect_equal(s2$cells$k, c(10, 20))
  expect_error(simulate_areal(fr, ti, cells = rbind(lab, lab)), "more than one row")
  expect_error(simulate_areal(fr, ti, cells = data.frame(id = "zz", t = 1)), "not in the frame")
})
