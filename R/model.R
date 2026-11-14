# model.R
#
# The model: covariance pieces, the observed-cell covariance, priors and the
# fixed-effect posterior. Internal; no I/O.
#
# Ported from Africa_areal/R/model.R (2026-10-12). Changes: the design
# matrix X of the observed cells is an input (any fixed effects, not only
# intercept + trend); the trend runs from -1 at the first time step to 1 at
# the last (Africa_areal: +/-(T-1)/T); the prior closure takes the
# package's prior names; the likelihood moved to likelihood.R.
#
# Notation:
#   m areas, T time steps, N = m*T cells; n observed cells with area index a_o
#   and time index t_o (both 1-based). Sigma = Sigma_t (x) Sigma_s, so
#   Sigma_oo = Sigma_t[t_o, t_o] * Sigma_s[a_o, a_o] elementwise.
#   Spatial: Leroux, Q_s = tau (rho_s R + (1 - rho_s) I), R = D - W = V L V'
#   (the default; BYM2 and proper CAR below, "spatial families").
#   Temporal: AR1, Sigma_t[t, t'] = rho_t^|t - t'|, variance fixed at 1 (only
#   the product of the two variances is identifiable; tau carries the scale).
#   Mean: x = X beta + field, beta ~ N(b0, B), marginalised into the
#   covariance as X B X'. theta = (log sigma, logit rho_s, logit rho_t
#   [, log sigma_n]) with sigma^2 = 1 / tau; sigma_n is a cell-level nugget
#   on observed cells only (the predicted surface is the field without it).
#
# The model is fitted on the observed cells only (n latent values); the full
# area x time table is reconstructed afterwards by conditional simulation
# (predict.R). This is exact, not an approximation (tests/testthat/test-model.R).
#
# ------------------------------------------------------------------

# ---- covariance pieces ------------------------------------------------------

# eigenvalues of Sigma_s = Q_s^{-1}, given the Laplacian eigenvalues. At
# rho_s = 1 (ICAR) the constant direction has zero precision; it gets zero
# variance instead (the field sums to zero over the areas, and the intercept
# carries the level), so g is 0 there rather than Inf.
leroux_g <- function(lambda, tau, rho_s) {
  g <- 1 / (tau * (rho_s * lambda + 1 - rho_s))
  g[is.infinite(g)] <- 0
  g
}

# ---- spatial families -------------------------------------------------------
# Every spatial prior is Sigma_s = V diag(g) V' for a basis (V, values) made
# once from the graph; only g depends on theta, so everything downstream
# (Sigma_oo, Matheron's rule, the held-out moments) works for any family.
#   leroux  V, values: eigenvectors and eigenvalues lambda of R = D - W;
#           g = 1 / (tau (rho_s lambda + 1 - rho_s)).
#   bym2    (Riebler et al. 2016) the same basis; Sigma_s = sigma^2 (rho_s R*^+
#           + (1 - rho_s) I), R*^+ = R^+ / s scaled so that the geometric mean
#           of its diagonal (the ICAR part's marginal variances) is 1;
#           g = sigma^2 (rho_s lambda^+ / s + 1 - rho_s), lambda^+ = 1 / lambda
#           (0 for lambda = 0). sigma is then the marginal sd, roughly, and
#           rho_s the share of the variance that is spatially structured.
#   car     proper CAR, Q_s = tau (D - rho_s W) = tau D^{1/2} (I - rho_s Wn)
#           D^{1/2}, Wn = D^{-1/2} W D^{-1/2} = U diag(mu) U' (mu <= 1); so
#           V = D^{-1/2} U (not orthonormal; nothing needs it to be),
#           values = mu, g = 1 / (tau (1 - rho_s mu)).
# In all three rho_s = 0 is independent areas (with variance 1 / tau for
# leroux and bym2, 1 / (tau d_a) for car) and rho_s -> 1 the ICAR.
#
# The basis travels as `eig` (list(values, vectors), plus `family`, `scale`
# for bym2, and `fixed`, the hyperparameters held fixed, on the theta
# scale). A plain list(values, vectors) is the Leroux model with every
# hyperparameter free, so the default model runs exactly as it always has.
spatial_g <- function(eig, p) {
  switch(if (is.null(eig$family)) "leroux" else eig$family,
         leroux = leroux_g(eig$values, p$tau, p$rho_s),
         bym2 = p$sigma^2 * (p$rho_s * ifelse(eig$values > 0, 1 / eig$values, 0) / eig$scale + 1 - p$rho_s),
         car = 1 / (p$tau * (1 - p$rho_s * eig$values)))
}

# the basis of a family on a frame (W, its Laplacian eigendecomposition eig
# and n_components), with the fixed hyperparameters on the theta scale
spatial_basis <- function(frame, family = "leroux", fixed = NULL) {
  if (family == "leroux" && !length(fixed)) return(frame$eig)          # the default, untouched
  out <- switch(family,
                leroux = frame$eig,
                bym2 = {
                  e <- frame$eig
                  lp <- ifelse(e$values > 0, 1 / e$values, 0)
                  c(e, list(scale = exp(mean(log(rowSums(e$vectors^2 * rep(lp, each = nrow(e$vectors))))))))
                },
                car = {
                  d <- Matrix::rowSums(frame$W)
                  Wn <- as.matrix(frame$W) / sqrt(outer(d, d))
                  e <- eigen(Wn, symmetric = TRUE)
                  list(values = pmin(e$values, 1), vectors = e$vectors / sqrt(d))
                })
  out$family <- family
  out$fixed <- if (length(fixed)) fixed else NULL
  out
}

# Sigma_s[rows, rows] = V[rows,] diag(g) V[rows,]'  (only the areas asked for)
sigma_s_block <- function(V, g, rows) {
  Vr <- V[rows, , drop = FALSE]
  (Vr * rep(g, each = length(rows))) %*% t(Vr)
}

# Sigma_t for a stationary AR1 with unit variance
ar1_cov <- function(T, rho_t) {
  toeplitz(rho_t^(0:(T - 1)))      # = rho_t^|t - t'|, without T^2 calls to pow
}

# lower Cholesky factor of ar1_cov(T, rho_t) in closed form:
# L[1,1] = 1, L[t,1] = rho^(t-1), and for j >= 2, L[t,j] = rho^(t-j) sqrt(1-rho^2)
ar1_chol <- function(T, rho_t) {
  L <- rho_t^pmax(outer(seq_len(T), seq_len(T), "-"), 0)
  L[lower.tri(L, diag = TRUE) & col(L) > 1] <-
    L[lower.tri(L, diag = TRUE) & col(L) > 1] * sqrt(1 - rho_t^2)
  L[upper.tri(L)] <- 0
  L
}

# the time trend: -1 at the first time step, 1 at the last, linear in
# between (0 throughout when T = 1)
trend_value <- function(t, T) if (T > 1) (2 * t - T - 1) / (T - 1) else 0 * t

# the default fixed effects: intercept + trend. Observed cells and the full
# table use the same function (via design_eval) so the two never disagree.
design_X <- function(t, T) cbind(1, trend_value(t, T))

# ---- parameters -------------------------------------------------------------

# The sampler moves only the free hyperparameters. `fixed` holds the others,
# named, on the unconstrained scale (logit_rho_s = Inf for rho_s = 1);
# theta_full() puts the two back in the order log_sigma, logit_rho_s,
# logit_rho_t[, log_sigma_n].
theta_all <- c("log_sigma", "logit_rho_s", "logit_rho_t", "log_sigma_n")
# Additive models (main effects, or an interaction other than type IV;
# model_basis() below) name every hyperparameter they have in `all`, from
# this list in this order; theta is then matched to it by name. For the
# default model `all` is NULL and everything is positional, as it always was.
theta_all_ext <- c(theta_all, "log_sigma_space", "logit_rho_space", "log_sigma_time", "logit_rho_time")
theta_full <- function(theta, fixed = NULL, all = NULL) {
  if (!is.null(all)) {
    out <- stats::setNames(numeric(length(all)), all)
    out[names(fixed)] <- fixed
    out[setdiff(all, names(fixed))] <- theta
    return(out)
  }
  if (is.null(fixed)) return(theta)
  nm <- theta_all[seq_len(length(theta) + length(fixed))]
  out <- numeric(length(nm))
  is_fixed <- nm %in% names(fixed)
  out[is_fixed] <- fixed[nm[is_fixed]]
  out[!is_fixed] <- theta
  out
}

# theta on the unconstrained scale -> natural parameters. In an additive
# model a component that is absent has zero scale (an interaction of type
# "none"), and an absent correlation is 0 (iid time in types I and III,
# iid space in types I and II).
theta_to_par <- function(theta, fixed = NULL, all = NULL) {
  if (!is.null(all)) {
    th <- theta_full(theta, fixed, all)
    v <- function(nm, absent) if (nm %in% all) th[[nm]] else absent
    ls <- v("log_sigma", -Inf)
    return(list(sigma = exp(ls), tau = exp(-2 * ls), rho_s = plogis(v("logit_rho_s", -Inf)),
                rho_t = plogis(v("logit_rho_t", -Inf)), sigma_n = exp(v("log_sigma_n", -Inf)),
                sigma_space = exp(v("log_sigma_space", -Inf)), rho_space = plogis(v("logit_rho_space", -Inf)),
                sigma_time = exp(v("log_sigma_time", -Inf)), rho_time = plogis(v("logit_rho_time", -Inf))))
  }
  theta <- theta_full(theta, fixed)
  list(sigma = exp(theta[1]), tau = exp(-2 * theta[1]),
       rho_s = plogis(theta[2]), rho_t = plogis(theta[3]),
       sigma_n = if (length(theta) >= 4) exp(theta[4]) else 0)
}

# the free hyperparameters' names
theta_names <- function(nugget, fixed = NULL, all = NULL) {
  if (!is.null(all)) return(setdiff(all, names(fixed)))
  setdiff(c("log_sigma", "logit_rho_s", "logit_rho_t", if (nugget) "log_sigma_n"), names(fixed))
}

# natural-scale names and transforms of the theta names
theta_natural_names <- c(log_sigma = "sigma", logit_rho_s = "rho_s", logit_rho_t = "rho_t", log_sigma_n = "sigma_n",
                         log_sigma_space = "sigma_space", logit_rho_space = "rho_space",
                         log_sigma_time = "sigma_time", logit_rho_time = "rho_time")

# the same, as columns of natural-scale values (for traces and summaries).
# Columns are matched by name when named (a fit with fixed hyperparameters
# has only the free ones), by position otherwise.
theta_natural <- function(th) {
  th <- as.matrix(th)
  nm <- colnames(th)
  if (is.null(nm)) nm <- theta_all[seq_len(ncol(th))]
  out <- vapply(nm, function(v) if (startsWith(v, "log_")) exp(th[, v == nm]) else plogis(th[, v == nm]),
                numeric(nrow(th)))
  out <- matrix(out, nrow(th), length(nm), dimnames = list(NULL, unname(theta_natural_names[nm])))
  out
}

# ---- the observed-cell covariance -------------------------------------------

# the observed-cell bookkeeping: unique areas, inverse index, design matrix
# (X: n x p fixed-effect design of the observed cells; intercept + trend if
# not given)
make_obs <- function(a_o, t_o, T, X = design_X(t_o, T)) {
  stopifnot(length(a_o) == length(t_o), nrow(X) == length(a_o))
  ua <- sort(unique(a_o))
  list(a_o = a_o, t_o = t_o, T = T, n = length(a_o),
       ua = ua, inv = match(a_o, ua), X = X)
}

# Everything the sampler needs for one value of theta: the marginal
# covariance of x_o (field + nugget + fixed effects), its lower Cholesky
# factor and log-determinant.
build_sigma_oo <- function(theta, obs, eig, B, jitter = 1e-8) {
  p <- theta_to_par(theta, eig$fixed, eig$all)
  g <- spatial_g(eig, p)
  Ss_u <- sigma_s_block(eig$vectors, g, obs$ua)          # observed areas only
  St <- ar1_cov(obs$T, p$rho_t)
  S <- Ss_u[obs$inv, obs$inv] * St[obs$t_o, obs$t_o] + obs$X %*% B %*% t(obs$X)
  # main effects (additive models only; the default adds nothing here)
  main <- main_parts(eig, p)
  if (!is.null(main$gS)) S <- S + sigma_s_block(eig$space$vectors, main$gS, obs$ua)[obs$inv, obs$inv]
  if (!is.null(main$ST)) S <- S + main$ST[obs$t_o, obs$t_o]
  dg <- seq.int(1L, length(S), by = nrow(S) + 1L)        # diagonal, written in place
  S[dg] <- S[dg] + jitter + p$sigma_n^2                  # (diag<- would copy the n x n matrix)
  L <- t(chol(S))
  out <- list(S = S, L = L, logdet = 2 * sum(log(diag(L))), g = g, St = St, par = p)
  if (!is.null(eig$all)) out$main <- main
  out
}

# log N(x | mu, S) from a lower Cholesky factor
log_dmvnorm_chol <- function(x, mu, L, logdet) {
  z <- forwardsolve(L, x - mu)
  -0.5 * (length(x) * log(2 * pi) + logdet + sum(z^2))
}

# ---- prior ------------------------------------------------------------------

# log prior density of theta on the unconstrained scale, Jacobians included,
# as a closure over the prior settings: sigma ~ half-N(0, sigma_sd),
# rho_s ~ Beta(rho_s[1], rho_s[2]), rho_t ~ Beta(rho_t[1], rho_t[2]),
# sigma_n ~ half-N(0, sigma_n_sd). With `fixed` (theta scale, named), only
# the free hyperparameters have a prior; theta is then the free ones.
make_log_prior <- function(priors, fixed = NULL, all = NULL) {
  sd_s <- priors$sigma_sd; rs <- priors$rho_s; rt <- priors$rho_t
  sd_n <- priors$sigma_n_sd
  if (!is.null(all)) {
    # additive models: one term per free hyperparameter, by name. A scale
    # has a half-normal prior (log sigma + log 2 is the Jacobian and the
    # fold), a correlation a Beta prior (log rho (1 - rho) the Jacobian).
    half_n <- c(log_sigma = "sigma_sd", log_sigma_n = "sigma_n_sd", log_sigma_space = "sigma_space_sd",
                log_sigma_time = "sigma_time_sd")
    beta_p <- c(logit_rho_s = "rho_s", logit_rho_t = "rho_t", logit_rho_space = "rho_space",
                logit_rho_time = "rho_time")
    free <- setdiff(all, names(fixed))
    return(function(theta) {
      names(theta) <- free
      lp <- 0
      for (v in free) {
        th <- theta[[v]]
        lp <- lp + if (v %in% names(half_n)) {
          stats::dnorm(exp(th), 0, priors[[half_n[[v]]]], log = TRUE) + log(2) + th
        } else {
          r <- stats::plogis(th); ab <- priors[[beta_p[[v]]]]
          stats::dbeta(r, ab[1], ab[2], log = TRUE) + log(r * (1 - r))
        }
      }
      lp
    })
  }
  if (!is.null(fixed)) return(function(theta) {
    full <- theta_full(theta, fixed); p <- theta_to_par(full)
    free <- !theta_all[seq_along(full)] %in% names(fixed)
    terms <- c(dnorm(p$sigma, 0, sd_s, log = TRUE) + log(2) + full[1],
               dbeta(p$rho_s, rs[1], rs[2], log = TRUE) + log(p$rho_s * (1 - p$rho_s)),
               dbeta(p$rho_t, rt[1], rt[2], log = TRUE) + log(p$rho_t * (1 - p$rho_t)),
               if (length(full) >= 4) dnorm(p$sigma_n, 0, sd_n, log = TRUE) + log(2) + full[4])
    sum(terms[free])
  })
  function(theta) {
    p <- theta_to_par(theta)
    lp <- dnorm(p$sigma, 0, sd_s, log = TRUE) + log(2) + theta[1] +          # half-normal, + log sigma
      dbeta(p$rho_s, rs[1], rs[2], log = TRUE) + log(p$rho_s * (1 - p$rho_s)) +
      dbeta(p$rho_t, rt[1], rt[2], log = TRUE) + log(p$rho_t * (1 - p$rho_t))
    if (length(theta) >= 4) lp <- lp + dnorm(p$sigma_n, 0, sd_n, log = TRUE) + log(2) + theta[4]
    lp
  }
}

# ---- fixed effects, recovered after the fact --------------------------------

# beta | x_o, theta ~ N(m, V) with prior N(b0, B) and x_o = X beta + e,
# e ~ N(0, S_field). Takes the field-only covariance (no X B X' term).
beta_posterior <- function(x_o, X, S_field, b0, B) {
  Si_X <- solve(S_field, X)
  V <- solve(solve(B) + t(X) %*% Si_X)
  m <- V %*% (solve(B, b0) + t(Si_X) %*% x_o)
  list(mean = as.vector(m), cov = V)
}

# ---- additive components (main effects and interaction types) ---------------
# The latent field is a sum of independent components (Knorr-Held 2000):
#   f_ST  the space-time interaction, Sigma_t^int (x) Sigma_s^int, as above;
#         its time factor is AR1 (types II, IV) or iid (I, III), its space
#         factor the spatial prior (III, IV) or iid (I, II); absent for "none"
#   f_S   a space main effect, 1_T 1_T' (x) Sigma_S, Sigma_S from the spatial
#         prior with its own scale and correlation (sigma_space, rho_space)
#   f_T   a time main effect, Sigma_T (x) 1_m 1_m': AR1 (sigma_time,
#         rho_time) or a scaled first- or second-order random walk
# so on the observed cells
#   Sigma_oo = Sigma_t^int[t_o, t_o] * Sigma_s^int[a_o, a_o]
#              + Sigma_S[a_o, a_o] + Sigma_T[t_o, t_o] + X B X' + sigma_n^2 I,
# still one dense n x n matrix.
#
# Constraints go through the basis, as for the ICAR: an improper direction
# gets zero variance, never infinite precision. A random walk of order k has
# precision R_k = D_k' D_k (D_k the k-th differences), whose null space is
# the constant (RW1) or the constant and the linear (RW2); its covariance is
# the generalised inverse R_k^+, computed as (R_k + N N')^{-1} - N N' with N
# an orthonormal basis of the null space. So an RW1 sums to zero over time,
# and an RW2 also has zero slope: the intercept and trend() carry the level
# and the straight line. R_k^+ is scaled so that the geometric mean of its
# diagonal is 1 (Sorbye & Rue 2014, as for BYM2), so sigma_time is roughly
# the marginal sd whatever the number of time steps.
#
# The default model (main "none", interaction "IV") is not additive: its
# basis is the plain spatial basis, and none of this is used.

main_choices <- c("none", "space", "time", "both")
interaction_choices <- c("IV", "I", "II", "III", "none")
temporal_choices <- c("ar1", "rw1", "rw2")

# a model spec with the component fields filled in (fits made before they
# existed have the default components)
complete_model <- function(md) {
  if (is.null(md$main)) md$main <- "none"
  if (is.null(md$interaction)) md$interaction <- "IV"
  if (is.null(md$temporal)) md$temporal <- "ar1"
  md
}
is_additive <- function(md) {
  md <- complete_model(md)
  md$main != "none" || md$interaction != "IV"
}

# every hyperparameter of an additive model (free or fixed), in the order of
# theta_all_ext
model_theta_all <- function(md, nugget) {
  int <- md$interaction
  space <- md$main %in% c("space", "both"); time <- md$main %in% c("time", "both")
  c(if (int != "none") "log_sigma", if (int %in% c("III", "IV")) "logit_rho_s",
    if (int %in% c("II", "IV")) "logit_rho_t", if (isTRUE(nugget)) "log_sigma_n",
    if (space) c("log_sigma_space", "logit_rho_space"),
    if (time) c("log_sigma_time", if (md$temporal == "ar1") "logit_rho_time"))
}

# The basis of a model on a frame with T time steps: for the default model,
# spatial_basis() exactly as before. For an additive model, the same list
# (the interaction's space factor: the spatial prior for types III and IV,
# the Leroux basis, used at rho_s = 0, for I, II and "none") plus `space`
# (the space main effect's basis), `time` (time_basis()), the component
# types, and `all` (model_theta_all()).
model_basis <- function(frame, T, model, nugget = TRUE) {
  md <- complete_model(model)
  fx <- fixed_theta(md$fixed)
  if (!is_additive(md)) return(spatial_basis(frame, md$family, fx))
  structured <- md$interaction %in% c("III", "IV")
  fam_int <- if (structured) md$family else "leroux"
  out <- spatial_basis(frame, fam_int, NULL)
  out$family <- fam_int
  if (md$main %in% c("space", "both")) {
    sb <- spatial_basis(frame, md$family, NULL)
    out$space <- list(vectors = sb$vectors, values = sb$values, family = md$family, scale = sb$scale)
  }
  if (md$main %in% c("time", "both")) out$time <- time_basis(T, md$temporal)
  out$interaction <- md$interaction; out$main <- md$main; out$temporal <- md$temporal
  out$all <- model_theta_all(md, nugget)
  out$fixed <- fx
  out
}

# the time main effect's basis: AR1 (nothing to precompute), or the scaled
# generalised inverse of a random walk's precision and its eigenvectors
time_basis <- function(T, type) {
  if (type == "ar1") return(list(type = "ar1", T = T))
  k <- if (type == "rw1") 1L else 2L
  if (T <= k) stop(sprintf("a %s needs at least %d time steps", toupper(type), k + 1L), call. = FALSE)
  D <- diff(diag(T), differences = k)
  R <- crossprod(D)
  N <- if (k == 1L) matrix(1 / sqrt(T), T, 1) else qr.Q(qr(cbind(1, seq_len(T))))
  NN <- tcrossprod(N)
  Rp <- solve(R + NN) - NN
  Rp <- (Rp + t(Rp)) / 2
  scale <- exp(mean(log(diag(Rp))))
  Rp <- Rp / scale
  e <- eigen(Rp, symmetric = TRUE)
  vals <- e$values
  vals[vals < 1e-10 * max(vals)] <- 0                    # the constrained directions: exactly zero
  Rp <- (e$vectors * rep(vals, each = T)) %*% t(e$vectors)
  list(type = type, T = T, cov = Rp, vectors = e$vectors, values = vals, scale = scale)
}

# the main effects' covariance pieces at the natural parameters p: gS (the
# space main effect's eigenvalue weights) and ST (T x T); NULL for the
# default model
main_parts <- function(eig, p) {
  if (is.null(eig$all)) return(NULL)
  out <- list()
  if (!is.null(eig$space))
    out$gS <- spatial_g(eig$space, list(sigma = p$sigma_space, tau = p$sigma_space^-2, rho_s = p$rho_space))
  if (!is.null(eig$time)) out$ST <- time_cov(eig$time, p)
  out
}

time_cov <- function(tb, p) {
  if (tb$type == "ar1") p$sigma_time^2 * ar1_cov(tb$T, p$rho_time) else p$sigma_time^2 * tb$cov
}

# unconditional draws of the main effects from standard normals z
space_draw <- function(sb, gS, z) as.vector(sb$vectors %*% (sqrt(gS) * z))
time_draw <- function(tb, p, z) {
  if (tb$type == "ar1") p$sigma_time * as.vector(ar1_filter(matrix(z, 1), p$rho_time)) else
    p$sigma_time * as.vector(tb$vectors %*% (sqrt(tb$values) * z))
}
