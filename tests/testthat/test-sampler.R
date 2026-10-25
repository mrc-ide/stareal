# Checks of sampler.R (ported from Africa_areal/tests/test_sampler.R,
# checks 1-3; the chain runner and checkpoint tests arrive with fit_areal()):
#   1. Plain ESS with a flat likelihood leaves its Gaussian invariant.
#   2. ESS against the posterior-scaled Gaussian targets the right
#      distribution: with the Gaussian pseudo-likelihood at half strength as
#      the likelihood, the exact target is Gaussian and differs from both
#      the prior and the surrogate Gaussian, so the moments check the algebra
#      of m* and Sigma_post and the sign of the correction.
#   3. The centred and surrogate-data theta moves leave the prior invariant
#      (flat likelihood, with ESS refreshing x_o; compares deciles). On 8
#      cells: the centred move is pinned by x_o, which is the very property
#      the surrogate move exists to get round.
# These are long Monte Carlo runs (about a minute in all), so not on CRAN.

tz <- make_toy()

test_that("plain ESS with a flat likelihood leaves its Gaussian invariant", {
  skip_on_cran()
  set.seed(11)
  with(tz, {
    x <- mu_o; ll <- 0
    xs <- matrix(NA, 20000, n)
    for (i in seq_len(nrow(xs))) {
      st <- ess_step(x, ll, mu_o, sig$L, function(x) 0)
      x <- st$x; xs[i, ] <- x
    }
    expect_lt(max(abs(colMeans(xs) - mu_o)) / max(sqrt(diag(sig$S))), 0.1)
    expect_lt(max(abs(apply(xs, 2, var) / diag(sig$S) - 1)), 0.25)
  })
})

test_that("ESS against the posterior-scaled Gaussian targets the exact posterior", {
  skip_on_cran()
  set.seed(12)
  with(tz, {
    g_star <- rnorm(n, 0, 2); S_aux <- runif(n, 0.2, 2)
    sig_e <- ess_surrogate_prep(sig, mu_o, g_star, S_aux)
    ll_half <- function(x) 0.5 * log_gauss_pseudo(x, g_star, S_aux)
    Q2 <- solve(sig$S) + diag(0.5 / S_aux); S2 <- solve(Q2)
    m2 <- as.vector(S2 %*% (solve(sig$S, mu_o) + 0.5 * g_star / S_aux))
    x <- sig_e$m_star; ll <- ll_half(x)
    xs <- matrix(NA, 30000, n)
    for (i in seq_len(nrow(xs))) {
      st <- ess_step_surrogate(x, ll, sig_e, g_star, S_aux, ll_half)
      x <- st$x; ll <- st$ll; xs[i, ] <- x
    }
    xs <- xs[-(1:3000), ]
    expect_lt(max(abs(colMeans(xs) - m2)) / max(sqrt(diag(S2))), 0.15)
    expect_lt(max(abs(apply(xs, 2, var) / diag(S2) - 1)), 0.25)
  })
})

test_that("the centred and surrogate-data theta moves leave the prior invariant", {
  skip_on_cran()
  set.seed(13)
  with(tz, {
    obs5 <- make_obs(a_o[1:8], t_o[1:8], T); n5 <- obs5$n
    mu5 <- as.vector(obs5$X %*% b0)
    S5 <- rep(0.5, n5); flat <- function(x) 0; prior <- make_log_prior(priors)
    prop <- diag(c(0.4, 0.8, 0.8))
    run_move <- function(move, n_iter = 50000) {
      th0 <- c(0, qlogis(0.5), qlogis(0.5))
      cur <- list(theta = th0, sig = build_sigma_oo(th0, obs5, eig, B))
      x <- mu5 + as.vector(cur$sig$L %*% rnorm(n5)); ll <- 0
      th <- matrix(NA, n_iter, 3)
      for (i in seq_len(n_iter)) {
        st <- ess_step(x, ll, mu5, cur$sig$L, flat); x <- st$x; ll <- st$ll
        if (move == "centred") {
          cur$lp <- prior(cur$theta) + log_dmvnorm_chol(x, mu5, cur$sig$L, cur$sig$logdet)
          cur <- mh_step_theta(cur, x, mu5, obs5, eig, B, prop, prior)
        } else {
          st <- mh_step_theta_surrogate(cur, x, ll, mu5, obs5, eig, B, S5, prop, prior, flat)
          cur <- st$cur; x <- st$x; ll <- st$ll
        }
        th[i, ] <- cur$theta
      }
      th[-(1:4000), ]
    }
    prior_draws <- cbind(log(abs(rnorm(2e5, 0, 2))), qlogis(rbeta(2e5, 2, 1)), qlogis(rbeta(2e5, 1, 1)))
    qp <- apply(prior_draws, 2, quantile, c(.1, .5, .9))
    for (mv in c("centred", "surrogate")) {
      err <- max(abs(apply(run_move(mv), 2, quantile, c(.1, .5, .9)) - qp))
      expect_lt(err, 0.35, label = paste(mv, "move: max |decile error| vs prior"))
    }
  })
})

test_that("adaptation keeps the proposal a valid Cholesky factor", {
  am <- am_init(c(0, 1, 2), c(0.1, 0.2, 0.3))
  set.seed(14)
  for (i in 1:200) am <- am_update(am, rnorm(3), runif(1))
  L <- am_prop_chol(am)
  expect_true(all(L[upper.tri(L)] == 0) && all(diag(L) > 0))
})

test_that("a proposal with a non-finite acceptance ratio is rejected, not an error", {
  with(tz, {
    nan_prior <- function(theta) if (identical(theta, tz$theta)) 0 else NaN
    cur <- list(theta = theta, sig = sig, lp = 0)
    st <- mh_step_theta(cur, mu_o, mu_o, obs, eig, B, diag(0.1, 3), nan_prior)
    expect_false(st$accepted)
    expect_identical(st$theta, theta)
    expect_equal(st$alpha, 0)
  })
})
