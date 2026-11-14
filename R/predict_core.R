# predict_core.R
#
# From one posterior draw of (theta, x_o) on the observed cells to the full
# area x time surface, by conditional simulation (Matheron's rule). Internal;
# no I/O. The user-facing predict() (draws in parallel, summaries, draw-wise
# aggregation to groups) is built on these.
#
# Ported from Africa_areal/R/predict.R (2026-10-12). Changes:
# - the fixed effects of the full table come from ctx$design(rows), a
#   function returning the (rows x T) x p design matrix, so the mean may
#   vary by area as well as by time (covariates); the Africa code had
#   intercept + trend by time only;
# - draw_noise() restores the caller's random-number state, so prediction
#   does not move the user's stream. The numbers drawn, and their order,
#   are unchanged.
#
# Matheron's rule, for one posterior draw (x_o, theta):
#   1. unconditional draw: beta* ~ N(b0, B); field* = V diag(sqrt g) E, where
#      E = Z Lt' (Z: m x T standard normals, Lt the AR1 Cholesky factor) is
#      computed as an AR1 recursion along each row of Z, O(mT) not O(mT^2);
#      at the observed cells add a nugget draw e* ~ N(0, sigma_n^2)
#   2. residual at the observed cells r = x_o - x*_o - e*; w = Sigma_oo^{-1} r
#   3. surface = x* + Sigma_s[, a_o] diag(w) Sigma_t[t_o, ] + X B X_o' w
# The surface is the FIELD plus fixed effects, without the nugget.
#
# Memory is bounded at any size by splitting the work in two: step 2 needs
# only the observed areas, and is done once per draw (prepare_draw), keeping
# w, beta* and the seed; step 3 is done for any subset of rows
# (predict_rows), regenerating the same Z from the seed.
#
# ctx: list(eig = spatial_basis() (the frame's list(values, vectors) for
#           the default model), obs = make_obs() output, b0, B, T,
#           design = function(rows) -> (length(rows) * T) x p matrix, rows in
#           the order (t - 1) * length(rows) + i [default: design_time(T)],
#           jitter [default 1e-8; the exactness tests set 0])
#
# ------------------------------------------------------------------

# ---- design of the full table -----------------------------------------------

# the default: intercept + trend, the same for every area (design_X)
design_time <- function(T) {
  force(T)
  function(rows) design_X(rep(seq_len(T), each = length(rows)), T)
}

ctx_design <- function(ctx) if (is.null(ctx$design)) design_time(ctx$T) else ctx$design

# ---- the unconditional field ------------------------------------------------

# E = Z Lt' for the unit-variance AR1 (= ar1_chol), row by row:
# e_1 = z_1, e_t = rho e_{t-1} + sqrt(1 - rho^2) z_t
ar1_filter <- function(Z, rho) {
  s <- sqrt(1 - rho^2)
  if (ncol(Z) > 1) for (t in 2:ncol(Z)) Z[, t] <- rho * Z[, t - 1] + s * Z[, t]
  Z
}

# the standard normals of draw `seed`, in a fixed order: Z first (shared by
# both phases), then beta and nugget noise (phase A only). The caller's
# random-number state is put back afterwards.
draw_noise <- function(seed, m, T, p = 0, n = 0, k = 0) {
  if (k == 0) return(with_seed(seed, list(Z = matrix(rnorm(m * T), m, T), bz = rnorm(p), ez = rnorm(n))))
  # additive models: k more, for the main effects, after all the others
  with_seed(seed, list(Z = matrix(rnorm(m * T), m, T), bz = rnorm(p), ez = rnorm(n), mz = rnorm(k)))
}

# ---- phase A: the conditioning step for one draw ----------------------------

# Returns what phase B needs for this draw.
prepare_draw <- function(x_o, theta, seed, ctx) {
  obs <- ctx$obs; V <- ctx$eig$vectors; m <- nrow(V)
  sig <- build_sigma_oo(theta, obs, ctx$eig, ctx$B, jitter = if (is.null(ctx$jitter)) 1e-8 else ctx$jitter)
  nz <- draw_noise(seed, m, ctx$T, length(ctx$b0), obs$n, k = if (is.null(ctx$eig$all)) 0 else m + ctx$T)
  E <- ar1_filter(nz$Z, sig$par$rho_t)
  beta_star <- ctx$b0 + as.vector(t(chol(ctx$B)) %*% nz$bz)
  field_ua <- (V[obs$ua, , drop = FALSE] * rep(sqrt(sig$g), each = length(obs$ua))) %*% E
  x_star_o <- field_ua[cbind(obs$inv, obs$t_o)] + as.vector(obs$X %*% beta_star)
  r <- x_o - x_star_o - nz$ez * sig$par$sigma_n
  if (!is.null(ctx$eig$all)) {
    # main effects: their unconditional draws, kept for phase B (m + T numbers)
    um <- main_draws(ctx$eig, sig, nz$mz, m, ctx$T)
    r <- r - um$fS[obs$a_o] - um$fT[obs$t_o]
  }
  w <- backsolve(t(sig$L), forwardsolve(sig$L, r))
  out <- list(seed = seed, theta = theta, w = as.vector(w), beta_star = beta_star)
  if (!is.null(ctx$eig$all)) out$main <- um
  out
}

# the main effects' unconditional draws for one posterior draw, from the
# standard normals mz (m for space, then T for time): fS (m, or 0) and fT
# (T, or 0)
main_draws <- function(eig, sig, mz, m, T) {
  fS <- if (!is.null(eig$space)) space_draw(eig$space, sig$main$gS, mz[seq_len(m)]) else numeric(m)
  fT <- if (!is.null(eig$time)) time_draw(eig$time, sig$par, mz[m + seq_len(T)]) else numeric(T)
  list(fS = fS, fT = fT)
}

# ---- phase B: rows of the surface for one draw ------------------------------

# the parts of predict_rows() that depend only on the rows, not the draw
rows_parts <- function(rows, ctx) {
  list(Vr = ctx$eig$vectors[rows, , drop = FALSE],
       Vt_ua = t(ctx$eig$vectors[ctx$obs$ua, , drop = FALSE]),
       Xr = ctx_design(ctx)(rows))
}

# logit surface for the areas `rows` (all time steps) given prepare_draw()
# output; a rows x T matrix. `parts` (rows_parts) can be passed in when the
# same rows are predicted for many draws.
predict_rows <- function(rows, prep, ctx, parts = rows_parts(rows, ctx)) {
  obs <- ctx$obs; m <- nrow(ctx$eig$vectors); T <- ctx$T; nr <- length(rows)
  p <- theta_to_par(prep$theta, ctx$eig$fixed, ctx$eig$all)
  g <- spatial_g(ctx$eig, p)
  E <- ar1_filter(draw_noise(prep$seed, m, T)$Z, p$rho_t)
  Vr <- parts$Vr
  field <- (Vr * rep(sqrt(g), each = nr)) %*% E                      # rows x T
  Ss_r <- (Vr * rep(g, each = nr)) %*% parts$Vt_ua                   # rows x |ua|
  St_o <- ar1_cov(T, p$rho_t)[obs$t_o, , drop = FALSE]               # n x T
  corr_field <- (Ss_r[, obs$inv, drop = FALSE] * rep(prep$w, each = nr)) %*% St_o
  # x* mean + fixed-effect correction, cell (i, t) at (t - 1) * nr + i
  mean_rt <- as.vector(parts$Xr %*% (prep$beta_star + ctx$B %*% (t(obs$X) %*% prep$w)))
  out <- field + corr_field + mean_rt
  if (!is.null(ctx$eig$all)) out <- out + main_rows(rows, prep, p, ctx)
  out
}

# the main effects' part of the surface for the areas `rows` (rows x T):
# unconditional draw + Sigma_S[rows, a_o] w (the same for every step) +
# Sigma_T[, t_o] w (the same for every area)
main_rows <- function(rows, prep, p, ctx) {
  obs <- ctx$obs; nr <- length(rows); T <- ctx$T; eig <- ctx$eig
  out <- matrix(0, nr, T)
  if (!is.null(eig$space)) {
    gS <- main_parts(eig, p)$gS
    w_area <- as.vector(rowsum(prep$w, obs$inv, reorder = TRUE))           # w summed by observed area
    Vr <- eig$space$vectors[rows, , drop = FALSE]
    Ss_r <- (Vr * rep(gS, each = nr)) %*% t(eig$space$vectors[obs$ua, , drop = FALSE])
    out <- out + (prep$main$fS[rows] + as.vector(Ss_r %*% w_area))
  }
  if (!is.null(eig$time)) {
    ST <- time_cov(eig$time, p)
    out <- out + rep(prep$main$fT + as.vector(ST[, obs$t_o, drop = FALSE] %*% prep$w), each = nr)
  }
  out
}

# the whole surface for one posterior draw (m x T, logit); used by the tests
# and for small frames
predict_surface <- function(x_o, theta, seed, ctx) {
  predict_rows(seq_len(nrow(ctx$eig$vectors)), prepare_draw(x_o, theta, seed, ctx), ctx)
}
