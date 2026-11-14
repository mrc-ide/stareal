# Checks of the additive components (main effects in space and time, the
# interaction types of Knorr-Held 2000, random walks in time):
#   1. The random-walk bases: the generalised inverse against an SVD
#      pseudo-inverse, the scaling, and the constraints (an RW1 sums to zero,
#      an RW2 also has zero slope), for the basis and for draws.
#   2. For each kind of model, against a dense full-covariance construction
#      made independently (Kronecker products of dense factors):
#      - Sigma_oo equals the dense model's observed block;
#      - Matheron's rule gives the exact kriging mean (zero noise) and the
#        exact conditional covariance (the surface is linear in its normals,
#        so unit-noise draws give it with no Monte Carlo error);
#      - crossval()'s held-out moments and loo()'s leave-one-out moments
#        equal dense conditioning.
#   3. model_spec() refuses what it should; theta is matched by name; the
#      prior has one term per free hyperparameter.
#   4. Fits, predictions and simulations of additive models: the right
#      hyperparameters, recovery of known values, and a checkpoint
#      fingerprint that tells models apart.

tz <- make_toy()
fr_add <- areal_frame(adjacency = tz$W, ids = sprintf("u%02d", seq_len(tz$m)), verbose = FALSE)

# ---- dense constructions, independent of the package's bases -----------------

pinv <- function(R) with(svd(R), v %*% diag(ifelse(d > 1e-9 * max(d), 1 / d, 0)) %*% t(u))
dense_rw <- function(T, k) {
  Rp <- pinv(crossprod(diff(diag(T), differences = k)))
  Rp / exp(mean(log(diag(Rp))))
}
dense_ar1 <- function(T, rho) outer(seq_len(T), seq_len(T), function(a, b) rho^abs(a - b))
dense_spatial <- function(fam, sigma, rho) {
  m <- tz$m; R <- tz$R_lap; D <- diag(rowSums(tz$W))
  switch(fam,
         leroux = solve((rho * R + (1 - rho) * diag(m)) / sigma^2),
         bym2 = { Rp <- pinv(R); sigma^2 * (rho * Rp / exp(mean(log(diag(Rp)))) + (1 - rho) * diag(m)) },
         car = solve((D - rho * tz$W) / sigma^2))
}
nat <- c(sigma = 0.8, rho_s = 0.7, rho_t = 0.6, sigma_n = 0.5, sigma_space = 0.6, rho_space = 0.8,
         sigma_time = 0.7, rho_time = 0.5)
# the full prior covariance of the latent field + fixed effects (cells time-major)
dense_full <- function(md) {
  T <- tz$T; m <- tz$m; int <- md$interaction
  S <- tz$X_full %*% tz$B %*% t(tz$X_full)
  if (int != "none") {
    St <- if (int %in% c("II", "IV")) dense_ar1(T, nat[["rho_t"]]) else diag(T)
    Ss <- if (int %in% c("III", "IV")) dense_spatial(md$family, nat[["sigma"]], nat[["rho_s"]]) else nat[["sigma"]]^2 * diag(m)
    S <- S + kronecker(St, Ss)
  }
  if (md$main %in% c("space", "both"))
    S <- S + kronecker(matrix(1, T, T), dense_spatial(md$family, nat[["sigma_space"]], nat[["rho_space"]]))
  if (md$main %in% c("time", "both")) {
    ST <- nat[["sigma_time"]]^2 * switch(md$temporal, ar1 = dense_ar1(T, nat[["rho_time"]]),
                                         rw1 = dense_rw(T, 1), rw2 = dense_rw(T, 2))
    S <- S + kronecker(ST, matrix(1, m, m))
  }
  S
}
# theta of a model at the values in `nat`
theta_at <- function(all) {
  v <- nat[unname(theta_natural_names[all])]
  ifelse(startsWith(all, "log_"), log(v), qlogis(v))
}

models <- list(
  default = list(),
  time_rw2 = list(main = "time", temporal = "rw2"),
  time_rw1_III_bym2 = list(main = "time", temporal = "rw1", interaction = "III", spatial = "bym2"),
  time_ar1_II_car = list(main = "time", interaction = "II", spatial = "car"),
  space_I = list(main = "space", interaction = "I", nugget = FALSE),
  both_rw2_car = list(main = "both", temporal = "rw2", spatial = "car"),
  space_none_bym2 = list(main = "space", interaction = "none", spatial = "bym2")
)
make_case <- function(a) {
  nug <- if (is.null(a$nugget)) TRUE else a$nugget
  md <- do.call(model_spec, c(a[setdiff(names(a), "nugget")], list(nugget = nug, frame = fr_add, n_steps = tz$T)))
  bas <- model_basis(fr_add, tz$T, md, nug)
  all <- if (is.null(bas$all)) theta_all[seq_len(3 + nug)] else bas$all
  list(md = md, bas = bas, nug = nug, th = theta_at(all), S_full = dense_full(md))
}

# ---- 1. random walks ------------------------------------------------------------

test_that("random-walk bases are the scaled generalised inverse, with their constraints", {
  for (T in c(7, 30)) for (k in 1:2) {
    tb <- time_basis(T, paste0("rw", k))
    expect_lt(max(abs(tb$cov - dense_rw(T, k))), 1e-8)
    expect_equal(exp(mean(log(diag(tb$cov)))), 1, tolerance = 1e-10)
    expect_lt(max(abs(tb$cov %*% rep(1, T))), 1e-8)                               # sums to zero
    if (k == 2) expect_lt(max(abs(tb$cov %*% seq_len(T))), 1e-7)                  # zero slope
    set.seed(T + k)
    f <- time_draw(tb, list(sigma_time = 1.3), rnorm(T))
    expect_lt(abs(sum(f)), 1e-8)
    if (k == 2) expect_lt(abs(sum(f * seq_len(T))), 1e-7)
  }
  expect_error(time_basis(2, "rw2"), "at least 3")
  # an AR1 time effect's draws have the AR1 covariance (the filter is exact)
  tb <- time_basis(5, "ar1")
  L <- vapply(1:5, function(i) time_draw(tb, list(sigma_time = 2, rho_time = 0.4), diag(5)[, i]), numeric(5))
  expect_lt(max(abs(L %*% t(L) - 4 * dense_ar1(5, 0.4))), 1e-12)
})

# ---- 2. exactness against the dense model ---------------------------------------

test_that("Sigma_oo, the kriging mean and Matheron's covariance are exact for every kind of model", {
  local_mocked_bindings(draw_noise = unit_noise)
  for (nm in names(models)) {
    cs <- make_case(models[[nm]])
    with(tz, {
      sn2 <- if (cs$nug) nat[["sigma_n"]]^2 else 0
      S_oo <- cs$S_full[obs_cells, obs_cells] + diag(sn2, n)
      sig <- build_sigma_oo(cs$th, obs, cs$bas, B, jitter = 0)
      expect_lt(max(abs(sig$S - S_oo)), 1e-9, label = paste(nm, "Sigma_oo"))
      set.seed(5)
      x_o <- mu_o + as.vector(sig$L %*% rnorm(n))
      mu_full <- as.vector(X_full %*% b0)
      K <- cs$S_full[, obs_cells] %*% solve(S_oo)
      m_exact <- mu_full + K %*% (x_o - mu_full[obs_cells])
      C_exact <- cs$S_full - K %*% t(cs$S_full[, obs_cells])
      ctx <- list(eig = cs$bas, obs = obs, b0 = b0, B = B, T = T, jitter = 0)
      s0 <- as.vector(predict_surface(x_o, cs$th, 0, ctx))
      expect_lt(max(abs(s0 - m_exact)), 1e-8, label = paste(nm, "kriging mean"))
      n_noise <- m * T + length(b0) + n + if (is.null(cs$bas$all)) 0 else m + T
      D <- vapply(seq_len(n_noise), function(s) as.vector(predict_surface(x_o, cs$th, s, ctx)) - s0, numeric(m * T))
      expect_lt(max(abs(D %*% t(D) - C_exact)), 1e-8, label = paste(nm, "Matheron covariance"))
    })
  }
})

test_that("held-out and leave-one-out moments equal dense conditioning for every kind of model", {
  for (nm in names(models)) {
    cs <- make_case(models[[nm]])
    with(tz, {
      sn2 <- if (cs$nug) nat[["sigma_n"]]^2 else 0
      S_all <- cs$S_full + diag(sn2, m * T)                       # every cell with its own nugget
      S_oo <- S_all[obs_cells, obs_cells]
      set.seed(7)
      x_o <- mu_o + as.vector(t(chol(S_oo)) %*% rnorm(n))
      mu_full <- as.vector(X_full %*% b0)
      held <- setdiff(seq_len(m * T), obs_cells)[c(1, 9, 40, 77)]
      a_h <- (held - 1) %% m + 1; t_h <- (held - 1) %/% m + 1
      mo <- heldout_moments(cs$th, x_o, obs, cs$bas, B, b0, X_full[held, , drop = FALSE], a_h, t_h, jitter = 0)
      A <- S_all[held, obs_cells] %*% solve(S_oo)
      expect_lt(max(abs(mo$mean - (mu_full[held] + A %*% (x_o - mu_full[obs_cells])))), 1e-9, label = paste(nm, "held mean"))
      expect_lt(max(abs(mo$cov - (S_all[held, held] - A %*% S_all[obs_cells, held]))), 1e-9, label = paste(nm, "held cov"))
      lm <- loo_moments(cs$th, x_o, obs, cs$bas, B, mu_o, jitter = 0)
      for (i in c(1, 13, n)) {
        a_i <- S_oo[i, -i] %*% solve(S_oo[-i, -i])
        expect_lt(abs(lm$mean[i] - (mu_o[i] + a_i %*% (x_o[-i] - mu_o[-i]))), 1e-9, label = paste(nm, "loo mean"))
        expect_lt(abs(lm$var[i] - (S_oo[i, i] - a_i %*% S_oo[-i, i])), 1e-9, label = paste(nm, "loo var"))
      }
    })
  }
})

# ---- 3. the model specification and theta ---------------------------------------

test_that("model_spec() sets up the components and refuses what it should", {
  md <- model_spec(main = "both", temporal = "rw2", interaction = "III")
  expect_identical(model_theta_all(md, TRUE),
                   c("log_sigma", "logit_rho_s", "log_sigma_n", "log_sigma_space", "logit_rho_space", "log_sigma_time"))
  expect_identical(model_theta_all(model_spec(main = "time", interaction = "none"), FALSE),
                   c("log_sigma_time", "logit_rho_time"))
  expect_equal(model_spec("icar", main = "space")$fixed, c(rho_s = 1, rho_space = 1))
  expect_equal(model_spec("iid", main = "space", interaction = "II")$fixed, c(rho_space = 0))
  expect_false(is_additive(model_spec()))
  expect_identical(model_spec()[c("main", "interaction", "temporal")], list(main = "none", interaction = "IV", temporal = "ar1"))
  expect_error(model_spec(interaction = "I"), "use nugget = FALSE")
  expect_silent(model_spec(interaction = "I", nugget = FALSE))
  expect_error(model_spec(interaction = "none"), "needs a main effect")
  expect_error(model_spec(temporal = "rw1"), "time main effect")
  expect_error(model_spec(main = "sideways"), "must be one of")
  expect_error(model_spec(interaction = "V"), "must be one of")
  expect_error(model_spec(interaction = "III", fixed = list(rho_t = 0.5)), "has no rho_t")
  expect_error(model_spec(main = "time", temporal = "rw1", fixed = list(rho_time = 0.5)), "has no rho_time")
  expect_error(model_spec(fixed = list(sigma_space = 1)), "has no sigma_space")
  expect_error(model_spec(main = "time", n_steps = 1), "more than one time step")
  expect_error(model_spec(main = "space", n_steps = 1), "duplicates")
  expect_silent(model_spec(main = "space", interaction = "II", n_steps = 1))
  expect_error(model_spec(main = "time", temporal = "rw2", n_steps = 2), "at least three")
  expect_error(model_spec(main = "space", fixed = list(rho_space = 1.5)), "between 0 and 1")
  expect_error(model_spec(main = "time", fixed = list(sigma_time = -1)), "must be positive")
  expect_error(model_spec(main = "time", interaction = "none", nugget = FALSE, fixed = list(sigma_time = 1, rho_time = 0.5)),
               "leave at least one free")
  expect_identical(model_label(model_spec("bym2", main = "both", temporal = "rw2", interaction = "III")),
                   "bym2; main space + time (rw2); interaction III")
  expect_identical(model_label(model_spec("car", fixed = list(rho_t = 0.5))), "car, rho_t = 0.5")
})

test_that("theta is matched by name in additive models, and the prior has one term per free value", {
  all <- c("log_sigma", "log_sigma_n", "log_sigma_space", "logit_rho_space", "log_sigma_time")
  fx <- c(logit_rho_space = qlogis(0.3))
  full <- theta_full(c(0.1, 0.2, 0.3, 0.4), fx, all)
  expect_equal(unname(full), c(0.1, 0.2, 0.3, qlogis(0.3), 0.4))
  p <- theta_to_par(c(log(1.2), log(0.4), log(0.6), log(0.9)), fx, all)
  expect_equal(c(p$sigma, p$sigma_n, p$sigma_space, p$rho_space, p$sigma_time), c(1.2, 0.4, 0.6, 0.3, 0.9))
  expect_identical(c(p$rho_s, p$rho_t), c(0, 0))                     # absent: iid in space and time
  expect_identical(theta_names(TRUE, fx, all), setdiff(all, names(fx)))
  pri <- areal_priors(sigma_space_sd = 0.7, sigma_time_sd = 1.5)
  lp <- make_log_prior(pri, fx, all)
  hn <- function(s, sd) log(2 * dnorm(s, 0, sd)) + log(s)
  direct <- hn(1.2, 2) + hn(0.4, 1) + hn(0.6, 0.7) + hn(0.9, 1.5)
  expect_lt(abs(lp(c(log(1.2), log(0.4), log(0.6), log(0.9))) - direct), 1e-12)
  lp2 <- make_log_prior(pri, NULL, c("logit_rho_time", "log_sigma_time"))
  expect_lt(abs(lp2(c(qlogis(0.4), log(0.9))) - (dbeta(0.4, 1, 1, log = TRUE) + log(0.4 * 0.6) + hn(0.9, 1.5))), 1e-12)
  expect_error(areal_priors(rho_time = 1), "two positive")
})

# ---- 4. fits, predictions, simulations -------------------------------------------

test_that("simulations with components keep their truth, and default simulations are unchanged", {
  tf <- time_frame(n_steps = 9)
  s <- simulate_areal(fr_add, tf, n_cells = 40, seed = 3, main = "both", temporal = "rw2", interaction = "II",
                      sigma_time = 0.8)
  expect_identical(names(s$truth$theta), c("sigma", "rho_t", "sigma_n", "sigma_space", "rho_space", "sigma_time"))
  expect_length(s$truth$time, 9); expect_length(s$truth$space, tz$m)
  expect_lt(abs(sum(s$truth$time)), 1e-8)
  expect_lt(abs(sum(s$truth$time * 1:9)), 1e-7)
  s0 <- simulate_areal(fr_add, tf, n_cells = 40, seed = 3)
  expect_null(s0$truth$main)
  expect_identical(names(s0$truth$theta), c("sigma", "rho_s", "rho_t", "sigma_n"))
  expect_error(simulate_areal(fr_add, tf, temporal = "rw1"), "time main effect")
  expect_error(simulate_areal(fr_add, tf, interaction = "none"), "needs a main effect")
})

test_that("additive fits report their hyperparameters, predict, and have their own checkpoint key", {
  sim <- simulate_areal(fr_add, time_frame(n_steps = 6), n_cells = 40, seed = 4, main = "time", temporal = "rw1")
  ctl <- list(n_iter = 60, n_burn = 20, n_store = 10, n_chains = 2, n_cores = 1, verbose = FALSE)
  f <- fit_areal(sim, main = "time", temporal = "rw1", interaction = "III", control = ctl)
  expect_identical(dimnames(f$draws$theta)[[3]], c("log_sigma", "logit_rho_s", "log_sigma_n", "log_sigma_time"))
  s <- summary(f)
  expect_identical(s$table$variable[s$table$type == "theta"], c("sigma", "rho_s", "sigma_n", "sigma_time"))
  expect_output(print(f), "model: leroux; main time \\(rw1\\); interaction III")
  pr <- predict(f, n_draws = 5, verbose = FALSE)
  expect_true(all(is.finite(pr$area$mean)))
  cv <- crossval(f, folds = 2)
  expect_true(all(is.finite(cv$cells$log_score)))
  # the checkpoint key: default, additive, and another additive model differ
  cells <- sim$cells; X <- design_X(cells$t, 6)
  pri <- areal_priors(); pri$beta_mean <- c(0, 0); pri$beta_sd <- c(3, 3)
  mk <- function(...) fit_setup(cells, model_basis(fr_add, 6, model_spec(..., n_steps = 6), TRUE), 6, pri, TRUE, X)
  d <- tempfile("ck"); dir.create(d); on.exit(unlink(d, recursive = TRUE))
  cg <- function(st, sub) checkpoint_guard(file.path(d, sub), st, mcmc_control(verbose = FALSE))
  fps <- c(cg(mk(), "a"), cg(mk(main = "time"), "b"), cg(mk(main = "time", temporal = "rw2"), "c"),
           cg(mk(main = "time", interaction = "III"), "e"))
  expect_identical(anyDuplicated(fps), 0L)
  # the default model's key leaves out the main effects' priors, so it is what it was before they existed
  st0 <- mk()
  pri_old <- pri[setdiff(names(pri), additive_prior_names)]
  st_old <- fit_setup(cells, fr_add$eig, 6, pri_old, TRUE, X)
  expect_identical(cg(st0, "f"), cg(st_old, "g"))
})

test_that("an additive fit recovers known hyperparameters and fixed effects", {
  skip_on_cran()
  fr <- areal_lattice(6, 6)
  sim <- simulate_areal(fr, time_frame(n_steps = 12), n_cells = 260, k = 100, seed = 11, beta = c(-0.5, 0.8),
                        main = "both", temporal = "rw2", interaction = "II", sigma = 0.5, rho_t = 0.7,
                        sigma_n = 0.2, sigma_space = 0.8, rho_space = 0.8, sigma_time = 0.6)
  f <- fit_areal(sim, main = "both", temporal = "rw2", interaction = "II",
                 control = list(n_iter = 2500, n_burn = 1000, n_store = 300, n_chains = 2, n_cores = 1,
                                verbose = FALSE))
  ta <- trace_array(f)
  post <- apply(ta, 3, function(v) c(mean(v), stats::sd(v)))
  truth <- sim$truth$theta[colnames(post)]
  z <- (post[1, ] - truth) / post[2, ]
  expect_true(all(abs(z) < 3), info = paste(names(z), round(z, 2), collapse = ", "))
  b <- f$draws$beta
  zb <- (apply(b, 3, mean) - sim$truth$beta) / apply(b, 3, stats::sd)
  expect_true(all(abs(zb) < 3), info = paste(names(zb), round(zb, 2), collapse = ", "))
})
