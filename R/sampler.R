# sampler.R
#
# The MCMC moves on the observed cells. Per iteration (the loop is in
# chains.R):
#   1. `ess_sweeps` elliptical slice sampling updates of the latent values
#      x_o, against the posterior-scaled Gaussian (not the prior);
#   2. `theta_moves` pairs of adaptive random-walk Metropolis moves on theta:
#      (a) centred: theta | x_o;
#      (b) surrogate-data (Murray & Adams 2010): x_o moves with theta.
# Internal; no I/O.
#
# Ported verbatim from Africa_areal/R/sampler.R (2026-10-12).
#
# Why this sampler (Africa_areal plan_ADMIN1-crt76T_2026-10-10.md, section 10a):
# - the centred move alone left chains stuck where they started (theta is
#   pinned by n latent values); a whitened move connected the regions but
#   mixed at ~0.05% ESS per iteration (pinned by the well-sampled cells);
#   the surrogate-data move gives 1-3% and chains that agree;
# - plain ESS proposes on the prior-scaled ellipse, which the heavily
#   sampled cells pin: latent median ESS 14 of 3,200 draws in a 40k run
#   whose theta looked converged (2026-10-11). Proposing on the
#   posterior-scaled ellipse instead gave ~13x more latent mixing per
#   iteration at ADMIN2.
#
# ------------------------------------------------------------------

# ---- latent update: elliptical slice sampling -------------------------------

# One ESS update of x given a Gaussian "prior" N(mu, L L') and log-likelihood
# ll_fun (Murray, Adams & MacKay 2010). Returns the new x and its ll.
ess_step <- function(x, ll, mu, L, ll_fun) {
  nu <- as.vector(L %*% rnorm(length(x)))
  log_u <- ll + log(runif(1))
  th <- runif(1, 0, 2 * pi)
  lo <- th - 2 * pi; hi <- th
  repeat {
    xp <- mu + (x - mu) * cos(th) + nu * sin(th)
    llp <- ll_fun(xp)
    if (llp > log_u) return(list(x = xp, ll = llp))
    if (th < 0) lo <- th else hi <- th
    th <- runif(1, lo, hi)
  }
}

# ESS against the posterior-scaled Gaussian. With FIXED pseudo-data g* (the
# empirical logits) and noise variances S_aux,
#   N(x; mu, Sigma) = N(x; m*, Sigma_post) N(g*; mu, A) / N(g*; x, S_aux),
# so the target is N(x; m*, Sigma_post) x [L(x) / N(g*; x, S_aux)] up to a
# constant: ESS with the posterior-scaled Gaussian as the "prior" and the
# bracket as the "likelihood" (Nishihara, Murray & Adams 2014, section 3).
# Exact for any fixed g* (tests/testthat/test-sampler.R); g* only sets how
# well the ellipse fits. m* and the lower factor of Sigma_post depend on
# theta only and are cached on the sig object.
ess_surrogate_prep <- function(sig, mu, g_star, S_aux) {
  if (is.null(sig$sur)) sig$sur <- surrogate_parts(sig, S_aux)
  if (is.null(sig$m_star)) {
    sig$m_star <- surrogate_eval(sig$sur, g_star, mu, S_aux)$m
    sig$L_post <- t(sig$sur$U_post)
  }
  sig
}

# one update; `ll` is log L(x) in and out, the Gaussian pseudo-likelihood is
# handled inside
ess_step_surrogate <- function(x, ll, sig, g_star, S_aux, ll_fun) {
  llp_fun <- function(xp) ll_fun(xp) - log_gauss_pseudo(xp, g_star, S_aux)
  st <- ess_step(x, ll - log_gauss_pseudo(x, g_star, S_aux), sig$m_star, sig$L_post, llp_fun)
  list(x = st$x, ll = st$ll + log_gauss_pseudo(st$x, g_star, S_aux))
}

# ---- hyperparameter move (a): centred ---------------------------------------

# One joint MH step on theta given x_o. `cur` carries the current theta and
# its build_sigma_oo() output so nothing is recomputed on rejection.
# `prop_chol` is the lower Cholesky factor of the proposal covariance.
# Returns the acceptance probability too, for the adaptation.
mh_step_theta <- function(cur, x_o, mu_o, obs, eig, B, prop_chol, prior_fun) {
  # (the rejection branch must rebuild the list: c(cur, list(...)) once
  # appended a second `accepted` behind the stale first one, which `$` kept
  # returning — the "acceptance 1.00" of the first pilots)
  reject <- function(alpha) list(theta = cur$theta, sig = cur$sig, lp = cur$lp,
                                 accepted = FALSE, alpha = alpha)
  theta_p <- cur$theta + as.vector(prop_chol %*% rnorm(length(cur$theta)))
  sig_p <- tryCatch(build_sigma_oo(theta_p, obs, eig, B), error = function(e) NULL)
  if (is.null(sig_p)) return(reject(0))
  lp_p <- prior_fun(theta_p) + log_dmvnorm_chol(x_o, mu_o, sig_p$L, sig_p$logdet)
  alpha <- min(1, exp(lp_p - cur$lp))
  if (runif(1) < alpha) {
    list(theta = theta_p, sig = sig_p, lp = lp_p, accepted = TRUE, alpha = alpha)
  } else {
    reject(alpha)
  }
}

# ---- hyperparameter move (b): surrogate data ---------------------------------
# Murray & Adams (2010), section 3, with a Metropolis proposal instead of
# slice sampling. Surrogate Gaussian data g ~ N(x_o, S_aux), with S_aux set
# from the likelihood curvature of each cell (data only, so the move is
# exact), carry roughly the information the likelihood carries. theta then
# moves with the whitened residual of x_o GIVEN g held fixed:
#   Sigma_post = (Sigma^{-1} + S_aux^{-1})^{-1},  m = mu + Sigma (Sigma+S)^{-1} (g - mu)
#   eta = R^{-1} (x_o - m) with R R' = Sigma_post;  x_o' = m' + R' eta
# accepted with prob min(1, prior ratio x lik ratio x N(g; mu, Sigma'+S) / N(g; mu, Sigma+S)).

# The theta-dependent pieces, cached with the sig object once computed.
# With A = Sigma + D, D = diag(S_aux): Sigma_post = D - D A^{-1} D (equal to
# Sigma - Sigma A^{-1} Sigma), which needs only A^{-1} — chol2inv gets it in a
# third of the flops of the n x n solves the Sigma form needs. Factors are
# kept upper-triangular to avoid transposes.
surrogate_parts <- function(sig, S_aux, jitter = 1e-8) {
  dg <- seq.int(1L, length(sig$S), by = nrow(sig$S) + 1L)
  A <- sig$S; A[dg] <- A[dg] + S_aux
  U_A <- chol(A)
  Sp <- chol2inv(U_A)
  Sp <- Sp * tcrossprod(-S_aux, S_aux)
  Sp[dg] <- Sp[dg] + S_aux + jitter
  list(U_A = U_A, logdet_A = 2 * sum(log(diag(U_A))), U_post = chol(Sp))
}

# log N(g; mu, A) and the posterior mean m = mu + Sigma A^{-1} (g - mu)
# = g - D A^{-1} (g - mu), sharing one triangular solve
surrogate_eval <- function(sur, g, mu, S_aux) {
  z <- forwardsolve(sur$U_A, g - mu, upper.tri = TRUE, transpose = TRUE)
  list(m = g - S_aux * backsolve(sur$U_A, z),
       lg = -0.5 * (length(g) * log(2 * pi) + sur$logdet_A + sum(z^2)))
}

mh_step_theta_surrogate <- function(cur, x_o, ll, mu_o, obs, eig, B, S_aux, prop_chol, prior_fun, ll_fun) {
  if (is.null(cur$sig$sur)) cur$sig$sur <- surrogate_parts(cur$sig, S_aux)
  g <- x_o + rnorm(length(x_o), 0, sqrt(S_aux))
  e_cur <- surrogate_eval(cur$sig$sur, g, mu_o, S_aux)
  eta <- forwardsolve(cur$sig$sur$U_post, x_o - e_cur$m, upper.tri = TRUE, transpose = TRUE)
  stay <- list(cur = cur, x = x_o, ll = ll, accepted = FALSE, alpha = 0)

  theta_p <- cur$theta + as.vector(prop_chol %*% rnorm(length(cur$theta)))
  sig_p <- tryCatch(build_sigma_oo(theta_p, obs, eig, B), error = function(e) NULL)
  if (is.null(sig_p)) return(stay)
  sig_p$sur <- tryCatch(surrogate_parts(sig_p, S_aux), error = function(e) NULL)
  if (is.null(sig_p$sur)) return(stay)
  e_p <- surrogate_eval(sig_p$sur, g, mu_o, S_aux)
  x_p <- e_p$m + as.vector(crossprod(sig_p$sur$U_post, eta))
  ll_p <- ll_fun(x_p)
  alpha <- min(1, exp(prior_fun(theta_p) - prior_fun(cur$theta) + ll_p - ll + e_p$lg - e_cur$lg))
  if (runif(1) < alpha) {
    lp_p <- prior_fun(theta_p) + log_dmvnorm_chol(x_p, mu_o, sig_p$L, sig_p$logdet)
    list(cur = list(theta = theta_p, sig = sig_p, lp = lp_p, accepted = TRUE, alpha = alpha),
         x = x_p, ll = ll_p, accepted = TRUE, alpha = alpha)
  } else {
    stay$alpha <- alpha
    stay
  }
}

# ---- proposal adaptation ----------------------------------------------------
# Adaptive Metropolis with a global scale (Andrieu & Thoms 2008, algorithm 4):
# running mean and covariance of theta give the proposal's shape, and a
# log-scale is nudged towards acceptance 0.234 using each proposal's
# acceptance probability. Used during burn-in only, so the chain after
# burn-in is a plain Metropolis chain.
am_init <- function(theta0, sd0) {
  list(mu = theta0, cov = diag(sd0^2, length(theta0)), log_lambda = 0, t = 0)
}
am_update <- function(am, theta, alpha, target = 0.234) {
  am$t <- am$t + 1
  gamma <- (am$t + 20)^-0.7
  d <- theta - am$mu
  am$mu <- am$mu + gamma * d
  am$cov <- am$cov + gamma * (tcrossprod(d) - am$cov)
  am$log_lambda <- am$log_lambda + gamma * (alpha - target)
  am
}
am_prop_chol <- function(am, eps = 1e-10) {
  d <- length(am$mu)
  t(chol(exp(am$log_lambda) * (2.38^2 / d) * am$cov + diag(eps, d)))
}
