# simulate.R
#
# Synthetic data from the model itself, with the truth kept: for examples,
# the tutorials, and checks that a fit recovers known parameters. The field
# is drawn as prediction draws it (V diag(sqrt g) Z Lt', the AR1 factor
# applied as a recursion). The surveys then go through areal_data(), so a
# simulation is the same object a real data set is, and fit_areal(sim)
# works.
#
# ------------------------------------------------------------------

#' Simulate survey data from the model
#'
#' Draws one realisation of the latent surface on a frame and time frame,
#' then binomial surveys at randomly chosen area-time cells (one survey per
#' cell). The logit prevalence is the fixed effects set by `formula` (by
#' default an intercept and a linear trend running from -1 at the first
#' time step to 1 at the last) + a Leroux x AR1 field, with a cell-level
#' nugget added at each surveyed cell.
#'
#' @param frame A [areal_frame()].
#' @param time A [time_frame()].
#' @param formula,covariates,scale_covariates The fixed effects, as in
#'   [fit_areal()]. Simulating and fitting with the same settings gives
#'   effects on the same scale.
#' @param beta Fixed effects on the logit scale: one per design column, in
#'   order (intercept, `trend` if in the formula, then the covariates; see
#'   [areal_priors()]), or named by design column.
#' @param sigma Scale of the field: the standard deviation of the Leroux
#'   prior is `sigma` times a factor set by the graph (\eqn{\tau = 1/\sigma^2}).
#' @param rho_s Spatial dependence of the Leroux prior, in (0, 1).
#' @param rho_t Lag-one temporal correlation of the AR1, in (0, 1).
#' @param sigma_n Standard deviation of the cell-level nugget (0 for none).
#' @param n_cells Number of distinct area-time cells surveyed, chosen at
#'   random.
#' @param k Number tested per survey: one number, or one per cell.
#' @param cells Instead of `n_cells` random cells: the cells to survey, as a
#'   data frame with an area `id` column, a `t` (step number) or `time`
#'   (label) column, and optionally `k` (number tested). The `$cells` of an
#'   [areal_data()] object qualifies, so a simulation can copy a real survey
#'   design: its areas, time steps and sample sizes.
#' @param seed Random seed. If given, the caller's random-number stream is
#'   left untouched.
#' @return A [areal_data()] object (`stareal_data`), ready for
#'   [fit_areal()], with `$truth`: `theta` (natural scale), `beta` (named
#'   by design column), `logit` and `prevalence` (areas x time steps
#'   matrices of the surface without the nugget; prevalence in percent),
#'   `x` (the logit at each surveyed cell with its nugget, in the order of
#'   `$cells`: what a fit
#'   estimates for the observed cells) and `design`. `$surveys` has columns
#'   `id`, `t`, `time`, `date` (calendar frames), `y`, `k`.
#' @examples
#' fr <- areal_lattice(6, 6)
#' tf <- time_frame(2001, 2020, "year")
#' sim <- simulate_areal(fr, tf, n_cells = 60, seed = 1)
#' sim
#' @export
simulate_areal <- function(frame, time, formula = ~ trend(), covariates = NULL, scale_covariates = TRUE,
                           beta = c(-1, 1), sigma = 1, rho_s = 0.9, rho_t = 0.9, sigma_n = 0.3,
                           n_cells = 100, k = 50, seed = NULL, cells = NULL) {
  if (!inherits(frame, "stareal_frame")) stop("`frame` must be made by areal_frame()", call. = FALSE)
  if (!inherits(time, "stareal_time")) stop("`time` must be made by time_frame()", call. = FALSE)
  scalar <- function(x, nm, lo, hi, lo_ok = FALSE) {
    if (!is.numeric(x) || length(x) != 1 || !is.finite(x) || x < lo || (!lo_ok && x == lo) || x >= hi)
      stop(sprintf("`%s` must be one number in %s%g, %g)", nm, if (lo_ok) "[" else "(", lo, hi), call. = FALSE)
  }
  scalar(sigma, "sigma", 0, Inf); scalar(rho_s, "rho_s", 0, 1); scalar(rho_t, "rho_t", 0, 1)
  scalar(sigma_n, "sigma_n", 0, Inf, lo_ok = TRUE)
  m <- frame$m; T <- time$n_steps
  given <- NULL
  if (!is.null(cells)) {
    if (!is.data.frame(cells) || !"id" %in% names(cells) || !any(c("t", "time") %in% names(cells)))
      stop("`cells` must be a data frame with an `id` column and a `t` or `time` column", call. = FALSE)
    a_g <- match(id_chr(cells$id), frame$ids)
    if (anyNA(a_g)) stop(sprintf("`cells` has ids not in the frame, e.g. %s", cells$id[which(is.na(a_g))[1]]), call. = FALSE)
    t_g <- if ("t" %in% names(cells)) time_steps(cells$t, time, "`cells`", "t") else time_steps(cells$time, time, "`cells`", "time")
    given <- (t_g - 1L) * m + a_g                     # cell index (t - 1) * m + a
    if (anyDuplicated(given)) stop("`cells` has more than one row for some area-time cells", call. = FALSE)
    if ("k" %in% names(cells)) k <- cells$k
    n_cells <- length(given)
  }
  n_cells <- check_whole(n_cells, "n_cells")
  if (n_cells > m * T) stop(sprintf("n_cells = %d exceeds the %d area-time cells", n_cells, m * T), call. = FALSE)
  if (!is.numeric(k) || !length(k) %in% c(1, n_cells) || any(!is.finite(k)) || any(k < 1) || any(k != round(k)))
    stop("`k` must be whole numbers >= 1, one or one per cell", call. = FALSE)
  des <- build_design(formula, covariates, frame, time, scale = scale_covariates)
  if (!is.numeric(beta) || any(!is.finite(beta))) stop("`beta` must be finite numbers", call. = FALSE)
  if (!is.null(names(beta))) {
    if (!setequal(names(beta), des$names) || anyDuplicated(names(beta)))
      stop(sprintf("named `beta` must name each design column once: %s", paste(des$names, collapse = ", ")),
           call. = FALSE)
    beta <- beta[des$names]
  }
  if (length(beta) != length(des$names))
    stop(sprintf("`beta` needs %d values (%s)", length(des$names), paste(des$names, collapse = ", ")), call. = FALSE)
  with_seed(seed, {
    # the cells first, so the field noise is never the first draw after a
    # seed (as a prediction draw's is) and the two cannot coincide
    if (is.null(given)) {
      cells <- sort(sample.int(m * T, n_cells))        # cell index (t - 1) * m + a
    } else {
      stats::runif(1)                                  # (as above: the field noise is not the first draw)
      o_g <- order(given)
      cells <- given[o_g]
      if (length(k) == n_cells) k <- k[o_g]
    }
    g <- leroux_g(frame$eig$values, 1 / sigma^2, rho_s)
    E <- ar1_filter(matrix(rnorm(m * T), m, T), rho_t)
    field <- (frame$eig$vectors * rep(sqrt(g), each = m)) %*% E
    logit <- field + matrix(design_eval(des, rep(seq_len(m), T), rep(seq_len(T), each = m)) %*% beta, m, T)
    a <- (cells - 1L) %% m + 1L; t <- (cells - 1L) %/% m + 1L
    kk <- rep_len(as.numeric(k), n_cells)
    x <- logit[cbind(a, t)] + rnorm(n_cells, 0, sigma_n)
    pos <- rbinom(n_cells, kk, stats::plogis(x))
  })
  o <- order(a, t)                                     # cell_id order: area, then time step
  surveys <- data.frame(id = frame$ids[a], t = t, time = time$labels[t], stringsAsFactors = FALSE)
  if (!is.null(time$starts)) surveys$date <- time$starts[t]
  surveys$y <- pos; surveys$k <- kk
  surveys <- surveys[o, ]
  rownames(surveys) <- NULL
  d <- areal_data(surveys, frame, time, y = "y", k = "k", area = "id", t = "t", verbose = FALSE)
  dimnames(logit) <- list(frame$ids, time$labels)
  d$truth <- list(theta = c(sigma = sigma, rho_s = rho_s, rho_t = rho_t, sigma_n = sigma_n),
                  beta = stats::setNames(as.numeric(beta), des$names),
                  logit = logit, prevalence = 100 * stats::plogis(logit), x = x[o], design = des)
  d
}
