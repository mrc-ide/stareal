# Checks of simulate.R: output shape, reproducibility, the caller's stream
# untouched, and the simulated field's covariance against the dense
# Kronecker covariance of the model.

fr <- areal_lattice(3, 3, queen = FALSE)
ti <- time_frame(T = 4)

test_that("simulated data have the right shape and feed areal_data()", {
  sim <- simulate_areal(fr, ti, n_cells = 20, k = 30, seed = 1)
  expect_s3_class(sim, "stareal_sim")
  expect_equal(nrow(sim$data), 20)
  expect_true(all(sim$data$positive <= sim$data$tested))
  expect_equal(dim(sim$truth$logit), c(9, 4))
  d <- areal_data(sim$data, fr, ti, y = "positive", k = "tested", area = "id", t = "t", verbose = FALSE)
  expect_equal(nrow(d$cells), 20)
  expect_output(print(sim), "20 surveys")
  sim_d <- simulate_areal(fr, time_frame(2001, 2004, "year"), n_cells = 5, seed = 1)
  expect_s3_class(sim_d$data$date, "Date")
  expect_error(simulate_areal(fr, ti, n_cells = 37), "exceeds")
  expect_error(simulate_areal(fr, ti, n_cells = 5, k = 0), "whole numbers")
})

test_that("a seed reproduces the simulation and leaves the caller's stream alone", {
  set.seed(5); a <- runif(2)
  set.seed(5); s1 <- simulate_areal(fr, ti, n_cells = 10, seed = 3); b <- runif(2)
  expect_identical(a, b)
  expect_identical(s1$data, simulate_areal(fr, ti, n_cells = 10, seed = 3)$data)
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
