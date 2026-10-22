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
# ctx: list(eig = list(values, vectors), obs = make_obs() output, b0, B, T,
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
draw_noise <- function(seed, m, T, p = 0, n = 0) {
  old <- get0(".Random.seed", envir = globalenv(), inherits = FALSE)
  on.exit(if (is.null(old)) rm(".Random.seed", envir = globalenv())
          else assign(".Random.seed", old, envir = globalenv()))
  set.seed(seed)
  list(Z = matrix(rnorm(m * T), m, T), bz = rnorm(p), ez = rnorm(n))
}

# ---- phase A: the conditioning step for one draw ----------------------------

# Returns what phase B needs for this draw.
prepare_draw <- function(x_o, theta, seed, ctx) {
  obs <- ctx$obs; V <- ctx$eig$vectors; m <- nrow(V)
  sig <- build_sigma_oo(theta, obs, ctx$eig, ctx$B, jitter = if (is.null(ctx$jitter)) 1e-8 else ctx$jitter)
  nz <- draw_noise(seed, m, ctx$T, length(ctx$b0), obs$n)
  E <- ar1_filter(nz$Z, sig$par$rho_t)
  beta_star <- ctx$b0 + as.vector(t(chol(ctx$B)) %*% nz$bz)
  field_ua <- (V[obs$ua, , drop = FALSE] * rep(sqrt(sig$g), each = length(obs$ua))) %*% E
  x_star_o <- field_ua[cbind(obs$inv, obs$t_o)] + as.vector(obs$X %*% beta_star)
  r <- x_o - x_star_o - nz$ez * sig$par$sigma_n
  w <- backsolve(t(sig$L), forwardsolve(sig$L, r))
  list(seed = seed, theta = theta, w = as.vector(w), beta_star = beta_star)
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
  p <- theta_to_par(prep$theta)
  g <- leroux_g(ctx$eig$values, p$tau, p$rho_s)
  E <- ar1_filter(draw_noise(prep$seed, m, T)$Z, p$rho_t)
  Vr <- parts$Vr
  field <- (Vr * rep(sqrt(g), each = nr)) %*% E                      # rows x T
  Ss_r <- (Vr * rep(g, each = nr)) %*% parts$Vt_ua                   # rows x |ua|
  St_o <- ar1_cov(T, p$rho_t)[obs$t_o, , drop = FALSE]               # n x T
  corr_field <- (Ss_r[, obs$inv, drop = FALSE] * rep(prep$w, each = nr)) %*% St_o
  # x* mean + fixed-effect correction, cell (i, t) at (t - 1) * nr + i
  mean_rt <- as.vector(parts$Xr %*% (prep$beta_star + ctx$B %*% (t(obs$X) %*% prep$w)))
  field + corr_field + mean_rt
}

# the whole surface for one posterior draw (m x T, logit); used by the tests
# and for small frames
predict_surface <- function(x_o, theta, seed, ctx) {
  predict_rows(seq_len(nrow(ctx$eig$vectors)), prepare_draw(x_o, theta, seed, ctx), ctx)
}
