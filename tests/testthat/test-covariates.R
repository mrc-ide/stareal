# Checks of covariates.R and covariates through fit, predict and crossval:
#   1. The design rows equal a design built by hand from the raw tables
#      (area, time and area-time covariates, a factor, standardisation),
#      in both forms the code uses (full table for predict, cells for
#      fit and crossval).
#   2. The trend runs from -1 at the first step to 1 at the last, and is
#      dropped (with a message) when there is one step.
#   3. Bad inputs fail with clear messages: missing coverage, unknown
#      columns, duplicates, transforms, offsets, no intercept, name
#      clashes, non-finite or constant values.
#   4. Fixed effects the observed cells cannot inform are flagged; a prior
#      so wide it swamps the field is refused.
#   5. predict() with covariates equals brute force (Matheron's rule per
#      draw with the hand-built design), and crossval() runs with them.
#   6. Recovery: data simulated with covariate effects; posterior means
#      within 3 sd of the truth.

cov_inputs <- function(fr, T) {
  set.seed(51)
  list(area = data.frame(id = rev(fr$ids), alt = rnorm(fr$m, 500, 100),
                         zone = rep_len(c("lo", "mid", "hi"), fr$m)),
       time = data.frame(t = T:1, rain = runif(T, 0, 200)),
       cell = data.frame(id = rep(fr$ids, each = T), t = rep(seq_len(T), fr$m), ndvi = rnorm(fr$m * T)))
}

test_that("design rows equal a design built by hand", {
  fr <- areal_lattice(4, 4); tf <- time_frame(n_steps = 5); m <- fr$m; T <- tf$n_steps
  cv <- cov_inputs(fr, T)
  des <- build_design(~ trend() + alt + zone + rain + ndvi, cv, fr, tf)
  expect_identical(des$names, c("intercept", "trend", "alt", "zonelo", "zonemid", "rain", "ndvi"))
  # by hand, cells in the order (t - 1) * m + a
  a <- rep(seq_len(m), T); t <- rep(seq_len(T), each = m)
  alt <- cv$area$alt[match(fr$ids, cv$area$id)]; zone <- cv$area$zone[match(fr$ids, cv$area$id)]
  rain <- cv$time$rain[match(seq_len(T), cv$time$t)]
  ndvi <- cv$cell$ndvi[match(paste(fr$ids[a], t), paste(cv$cell$id, cv$cell$t))]
  st <- function(x) (x - mean(x)) / sd(x)
  X <- cbind(1, (2 * t - T - 1) / (T - 1), st(alt)[a], (zone == "lo")[a], (zone == "mid")[a],
             st(rain)[t], st(ndvi))
  expect_equal(unname(design_eval(des, a, t)), unname(X), tolerance = 1e-14)
  # the full-table form (predict) has rows in the order (t - 1) * nr + i
  rows <- c(2, 7, 11)
  Xr <- design_rows(des)(rows)
  expect_equal(unname(Xr), unname(X[as.vector(outer(rows, (seq_len(T) - 1) * m, "+")), ]), tolerance = 1e-14)
  # unscaled
  des0 <- build_design(~ alt, cv, fr, tf, scale = FALSE)
  expect_equal(design_eval(des0, 1:m, rep(1, m))[, "alt"], alt)
  # the stored formula does not carry the caller's environment
  expect_identical(environment(des$formula), emptyenv())
})

test_that("the trend runs from -1 to 1, and is dropped with one time step", {
  fr <- areal_lattice(3, 3); tf <- time_frame(n_steps = 7)
  des <- build_design(~ trend(), NULL, fr, tf)
  tr <- design_eval(des, rep(1, 7), 1:7)[, "trend"]
  expect_equal(tr[c(1, 7)], c(-1, 1)); expect_equal(diff(tr), rep(1 / 3, 6))
  expect_identical(build_design(~ 1, NULL, fr, tf)$names, "intercept")
  expect_message(d1 <- build_design(~ trend(), NULL, fr, time_frame(n_steps = 1)), "trend\\(\\) dropped")
  expect_identical(d1$names, "intercept")
})

test_that("bad covariates fail clearly", {
  fr <- areal_lattice(3, 3); tf <- time_frame(n_steps = 4)
  cv <- cov_inputs(fr, 4)
  expect_error(build_design(~ alt, cv$area[-1, ], fr, tf), "covers 8 of the 9 areas")
  expect_error(build_design(~ ndvi, cv$cell[-1, ], fr, tf), "covers 35 of the 36")
  expect_error(build_design(~ elev, cv, fr, tf), "`elev` is not a column")
  expect_error(build_design(~ log(alt), cv, fr, tf), "Compute transforms")
  expect_error(build_design(~ trend, cv, fr, tf), "parentheses")
  expect_error(build_design(~ alt + offset(rain), cv, fr, tf), "offset")
  expect_error(build_design(~ alt - 1, cv, fr, tf), "intercept")
  expect_error(build_design(y ~ alt, cv, fr, tf), "one-sided")
  expect_error(build_design(~ alt, rbind(cv$area, cv$area[1, ]), fr, tf), "more than one row")
  expect_error(build_design(~ alt, list(cv$area, cv$area), fr, tf), "more than one covariate table")
  expect_error(build_design(~ rain, data.frame(t = 1:4, rain = 1), fr, tf), "constant")
  expect_error(build_design(~ rain, data.frame(t = 1:4, rain = 1 + c(0, 1e-15, 0, 0)), fr, tf), "constant")
  na <- cv$area; na$alt[2] <- NA
  expect_error(build_design(~ alt, na, fr, tf), "missing or non-finite")
  na$alt[2] <- Inf
  expect_error(build_design(~ alt, na, fr, tf), "missing or non-finite")
  expect_error(build_design(~ alt, data.frame(alt = 1:9), fr, tf), "needs an `id` column")
  expect_error(build_design(~ rain, data.frame(time = c("1", "2", "3", "9"), rain = 1:4), fr, tf), "not in the time frame")
  expect_error(build_design(~ rain, data.frame(t = c(1, 2, 3, 3.5), rain = 1:4), fr, tf), "not in the time frame")
  # names that clash with the model's parameters or with each other
  expect_error(build_design(~ sigma, data.frame(id = fr$ids, sigma = rnorm(9)), fr, tf), "model parameter")
  clash <- data.frame(id = fr$ids, zone = rep(c("a", "b", "c"), 3), zoneb = rnorm(9))
  expect_error(build_design(~ zone + zoneb, clash, fr, tf), "another covariate")
  # a factor `t` key is read by its labels
  ft <- data.frame(t = factor(c("4", "1", "2", "3")), rain = c(4, 1, 2, 3))
  expect_equal(design_eval(build_design(~ rain, ft, fr, tf, scale = FALSE), rep(1, 4), 1:4)[, "rain"], c(1, 2, 3, 4))
})

test_that("uninformed fixed effects are flagged, and a swamping prior is refused", {
  fr <- areal_lattice(4, 4); tf <- time_frame(n_steps = 5)
  zone <- data.frame(id = fr$ids, zone = rep(c("a", "b"), each = 8))
  sim <- simulate_areal(fr, tf, n_cells = 30, seed = 1)
  sim$cells <- sim$cells[zone$zone[sim$cells$area] == "a", ]           # no data in zone b
  ctl <- list(n_iter = 40, n_burn = 20, n_store = 10, n_chains = 1, ess_sweeps = 2, verbose = FALSE)
  expect_warning(fit_areal(sim, formula = ~ trend() + zone, covariates = zone, control = ctl), "zoneb")
  big <- data.frame(id = fr$ids, pop = seq(1e5, 1e6, length.out = 16))
  expect_error(fit_areal(sim, formula = ~ pop, covariates = big, scale_covariates = FALSE, control = ctl), "swamps")
})

test_that("predict() with covariates equals brute force, and crossval() runs with them", {
  fr <- areal_lattice(4, 4); tf <- time_frame(n_steps = 5); m <- fr$m; T <- tf$n_steps
  cv <- cov_inputs(fr, T)
  f <- ~ trend() + alt + zone + rain + ndvi
  sim <- simulate_areal(fr, tf, formula = f, covariates = cv, beta = c(-1, 0.5, 0.3, 0.2, -0.2, 0.4, 0.3),
                        n_cells = 40, k = 30, seed = 3)
  ctl <- list(n_iter = 40, n_burn = 20, n_store = 10, n_chains = 1, ess_sweeps = 2, theta_moves = 1, verbose = FALSE)
  fit <- fit_areal(sim, formula = f, covariates = cv, control = ctl)
  expect_identical(dimnames(fit$draws$beta)[[3]], sim$truth$design$names)
  pr <- predict(fit, n_draws = 6, seed = 7, verbose = FALSE)
  # brute force: whole surfaces with the hand-checked design
  Xg <- design_eval(sim$truth$design, rep(seq_len(m), T), rep(seq_len(T), each = m))
  st <- fit_setup_of(fit)
  ctx <- list(eig = fr$eig, obs = st$obs, b0 = st$b0, B = st$B, T = T,
              design = function(rows) Xg[as.vector(outer(rows, (seq_len(T) - 1) * m, "+")), , drop = FALSE])
  pool <- function(a) matrix(a, prod(dim(a)[1:2]), dim(a)[3])
  th <- pool(fit$draws$theta); xo <- pool(fit$draws$x); idx <- pr$settings$draws_used
  seeds <- with_seed(7, sample.int(.Machine$integer.max - 1L, length(idx)))
  P <- sapply(seq_along(idx), function(d) as.vector(t(plogis(predict_surface(xo[idx[d], ], th[idx[d], ], seeds[d], ctx)))))
  expect_lt(max(abs(pr$area$mean - rowMeans(P))), 1e-10)
  # the fit's observed design is the same design
  expect_equal(unname(st$obs$X), unname(Xg[(sim$cells$t - 1) * m + sim$cells$area, ]))
  cvres <- crossval(fit, folds = 3)
  expect_equal(nrow(cvres$cells), nrow(sim$cells))
  expect_output(print(summary(fit)), "standardised covariates")
})

test_that("a fit recovers covariate effects", {
  skip_on_cran()
  # 64 units x 10 steps, 250 cells; an area covariate with effect 0.8 per
  # sd and a time covariate with effect -0.5 per sd; about 15 s.
  fr <- areal_lattice(8, 8); tf <- time_frame(n_steps = 10)
  set.seed(61)
  cv <- list(data.frame(id = fr$ids, x_area = rnorm(fr$m)), data.frame(t = 1:10, x_time = rnorm(10)))
  f <- ~ x_area + x_time
  sim <- simulate_areal(fr, tf, formula = f, covariates = cv, beta = c(-0.5, 0.8, -0.5), sigma = 0.7,
                        rho_s = 0.8, rho_t = 0.8, sigma_n = 0.3, n_cells = 250, k = 100, seed = 1)
  fit <- fit_areal(sim, formula = f, covariates = cv,
                   control = list(n_iter = 2000, n_burn = 700, n_store = 200, n_chains = 2, ess_sweeps = 5,
                                  n_cores = 2, verbose = FALSE))
  s <- summary(fit)$table
  b <- s[s$type == "beta", ]
  truth <- sim$truth$beta[b$variable]
  expect_true(all(abs(b$mean - truth) < 3 * b$sd))
  expect_gt(b$q2.5[b$variable == "x_area"], 0)          # well identified: its interval excludes 0
})
