# Checks of covariates.R and covariates through fit, predict and crossval:
#   1. The design rows equal a design built by hand from the raw tables
#      (area, time and area-time covariates, a factor, standardisation),
#      and the full-table and per-cell forms agree.
#   2. The default mean reproduces design_X exactly (the Africa_areal
#      design), so fits without covariates are unchanged.
#   3. Bad inputs fail with clear messages: missing coverage, unknown
#      columns, duplicates, transforms, no intercept.
#   4. predict() with covariates equals brute force (Matheron's rule per
#      draw with the hand-built design), and crossval() runs with them.
#   5. Recovery: data simulated with a covariate effect; the fit's 95%
#      interval covers it.

cov_inputs <- function(fr, T) {
  set.seed(51)
  list(area = data.frame(id = rev(fr$ids), alt = rnorm(fr$m, 500, 100),
                         zone = rep_len(c("lo", "mid", "hi"), fr$m)),
       time = data.frame(t = T:1, rain = runif(T, 0, 200)),
       cell = data.frame(id = rep(fr$ids, each = T), t = rep(seq_len(T), fr$m), ndvi = rnorm(fr$m * T)))
}

test_that("design rows equal a design built by hand", {
  fr <- areal_lattice(4, 4); tf <- time_frame(T = 5); m <- fr$m; T <- tf$T
  cv <- cov_inputs(fr, T)
  des <- build_design(~ trend() + alt + zone + rain + ndvi, cv, fr, tf)
  expect_identical(des$names, c("intercept", "trend", "alt", "zonelo", "zonemid", "rain", "ndvi"))
  # by hand, cells in the order (t - 1) * m + a
  a <- rep(seq_len(m), T); t <- rep(seq_len(T), each = m)
  alt <- cv$area$alt[match(fr$ids, cv$area$id)]; zone <- cv$area$zone[match(fr$ids, cv$area$id)]
  rain <- cv$time$rain[match(seq_len(T), cv$time$t)]
  ndvi <- cv$cell$ndvi[match(paste(fr$ids[a], t), paste(cv$cell$id, cv$cell$t))]
  st <- function(x) (x - mean(x)) / sd(x)
  X <- cbind(1, (t - (T + 1) / 2) / (T / 2), st(alt)[a], (zone == "lo")[a], (zone == "mid")[a],
             st(rain)[t], st(ndvi))
  expect_equal(unname(design_eval(des, a, t)), unname(X), tolerance = 1e-14)
  # the full-table form (predict) and the per-cell form (fit, crossval) agree
  fit_like <- list(design = des)
  rows <- c(2, 7, 11)
  Xr <- design_rows(fit_like)(rows)
  expect_equal(Xr, design_cells(fit_like, rep(rows, T), rep(seq_len(T), each = 3)))
  # unscaled
  des0 <- build_design(~ alt, cv, fr, tf, scale = FALSE)
  expect_equal(design_eval(des0, 1:m, rep(1, m))[, "alt"], alt)
})

test_that("the default mean is the Africa_areal design", {
  fr <- areal_lattice(3, 3); tf <- time_frame(T = 7)
  des <- build_design(~ trend(), NULL, fr, tf)
  a <- rep(1:9, 7); t <- rep(1:7, each = 9)
  expect_identical(unname(design_eval(des, a, t)), unname(design_X(t, 7)))
  expect_identical(build_design(~ 1, NULL, fr, tf)$names, "intercept")
})

test_that("bad covariates fail clearly", {
  fr <- areal_lattice(3, 3); tf <- time_frame(T = 4)
  cv <- cov_inputs(fr, 4)
  expect_error(build_design(~ alt, cv$area[-1, ], fr, tf), "covers 8 of the 9 areas")
  expect_error(build_design(~ ndvi, cv$cell[-1, ], fr, tf), "covers 35 of the 36")
  expect_error(build_design(~ elev, cv, fr, tf), "`elev` is not a column")
  expect_error(build_design(~ log(alt), cv, fr, tf), "Compute transforms")
  expect_error(build_design(~ alt - 1, cv, fr, tf), "intercept")
  expect_error(build_design(y ~ alt, cv, fr, tf), "one-sided")
  expect_error(build_design(~ alt, rbind(cv$area, cv$area[1, ]), fr, tf), "more than one row")
  expect_error(build_design(~ alt, list(cv$area, cv$area), fr, tf), "more than one covariate table")
  expect_error(build_design(~ rain, data.frame(t = 1:4, rain = 1), fr, tf), "constant")
  na <- cv$area; na$alt[2] <- NA
  expect_error(build_design(~ alt, na, fr, tf), "missing values")
  expect_error(build_design(~ alt, data.frame(alt = 1:9), fr, tf), "needs an `id` column")
  expect_error(build_design(~ rain, data.frame(time = c("1", "2", "3", "9"), rain = 1:4), fr, tf), "outside")
})

test_that("predict() with covariates equals brute force, and crossval() runs with them", {
  fr <- areal_lattice(4, 4); tf <- time_frame(T = 5); m <- fr$m; T <- tf$T
  cv <- cov_inputs(fr, T)
  mean <- ~ trend() + alt + zone + rain + ndvi
  sim <- simulate_areal(fr, tf, n_cells = 40, k = 30, seed = 3, mean = mean, covariates = cv,
                        beta = c(-1, 0.5, 0.3, 0.2, -0.2, 0.4, 0.3))
  d <- areal_data(sim$data, fr, tf, y = "positive", k = "tested", area = "id", t = "t", verbose = FALSE)
  ctl <- mcmc_control(n_iter = 40, n_burn = 20, n_store = 10, n_chains = 1, ess_sweeps = 2, theta_moves = 1,
                      progress = FALSE)
  fit <- fit_areal(d, mean = mean, covariates = cv, control = ctl)
  expect_identical(dimnames(fit$draws$beta)[[3]], sim$design$names)
  pr <- predict(fit, n_draws = 6, seed = 7, parallel = FALSE, progress = FALSE)
  # brute force: whole surfaces with the hand-checked design
  Xg <- design_eval(sim$design, rep(seq_len(m), T), rep(seq_len(T), each = m))
  ctx <- list(eig = fr$eig, obs = fit$setup$obs, b0 = fit$setup$b0, B = fit$setup$B, T = T,
              design = function(rows) Xg[as.vector(outer(rows, (seq_len(T) - 1) * m, "+")), , drop = FALSE])
  pool <- function(a) matrix(a, prod(dim(a)[1:2]), dim(a)[3])
  th <- pool(fit$draws$theta); xo <- pool(fit$draws$x); idx <- pr$settings$draws_used
  P <- sapply(seq_along(idx), function(d) as.vector(t(plogis(predict_surface(xo[idx[d], ], th[idx[d], ], 7 + d, ctx)))))
  expect_lt(max(abs(pr$area$mean - rowMeans(P))), 1e-10)
  # the fit's observed design is the same design
  expect_equal(unname(fit$setup$obs$X), unname(Xg[(d$cells$t - 1) * m + d$cells$area, ]))
  cvres <- crossval(fit, folds = 3)
  expect_equal(nrow(cvres$cells), nrow(d$cells))
  expect_output(print(summary(fit)), "standardised covariates")
})

test_that("a fit recovers a covariate effect", {
  skip_on_cran()
  # 64 units x 10 steps, 250 cells; an area covariate with effect 0.8 per
  # sd and a time covariate with effect -0.5 per sd; about 15 s. Seed 1
  # was the first tried.
  fr <- areal_lattice(8, 8); tf <- time_frame(T = 10)
  set.seed(61)
  cv <- list(data.frame(id = fr$ids, x_area = rnorm(fr$m)), data.frame(t = 1:10, x_time = rnorm(10)))
  mean <- ~ x_area + x_time
  sim <- simulate_areal(fr, tf, sigma = 0.7, rho_s = 0.8, rho_t = 0.8, sigma_n = 0.3, beta = c(-0.5, 0.8, -0.5),
                        n_cells = 250, k = 100, seed = 1, mean = mean, covariates = cv)
  d <- areal_data(sim$data, fr, tf, y = "positive", k = "tested", area = "id", t = "t", verbose = FALSE)
  fit <- fit_areal(d, mean = mean, covariates = cv,
                   control = mcmc_control(n_iter = 2000, n_burn = 700, n_store = 200, n_chains = 2,
                                          ess_sweeps = 5, n_cores = 2, progress = FALSE))
  s <- summary(fit)$table
  b <- s[s$type == "beta", ]
  truth <- sim$truth$beta[b$variable]
  expect_true(all(b$q2.5 <= truth & truth <= b$q97.5))
  # the area effect is well identified: its interval excludes 0
  expect_gt(b$q2.5[b$variable == "x_area"], 0)
})
