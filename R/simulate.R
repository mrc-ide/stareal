# simulate.R
#
# Synthetic data from the model itself, with the truth kept: for examples,
# the vignette, and checks that a fit recovers known parameters. The field
# is drawn exactly as prediction draws it (V diag(sqrt g) Z Lt', the AR1
# factor applied as a recursion), so simulation and prediction share one
# code path.
#
# ------------------------------------------------------------------

#' Simulate survey data from the model
#'
#' Draws one realisation of the latent surface on a frame and time frame,
#' then binomial surveys at randomly chosen area-time cells. The logit
#' prevalence is intercept + linear trend (scaled to run from -1 to 1 over
#' the window) + a Leroux x AR1 field, with a cell-level nugget added at
#' each surveyed cell.
#'
#' @param frame A [areal_frame()].
#' @param time A [time_frame()].
#' @param sigma Scale of the field: the standard deviation of the Leroux
#'   prior is `sigma` times a factor set by the graph (\eqn{\tau = 1/\sigma^2}).
#' @param rho_s Spatial dependence of the Leroux prior, in (0, 1).
#' @param rho_t Lag-one temporal correlation of the AR1, in (0, 1).
#' @param sigma_n Standard deviation of the cell-level nugget (0 for none).
#' @param beta Intercept and trend on the logit scale.
#' @param n_cells Number of distinct area-time cells surveyed.
#' @param k Number tested per survey: one number, or one per cell.
#' @param seed Random seed. If given, the caller's random-number stream is
#'   left untouched.
#' @return An object of class `stareal_sim`:
#'   \describe{
#'     \item{data}{one row per survey: `id`, `t`, `time` (label), `date`
#'       (start of the step; calendar frames only), `tested`, `positive`}
#'     \item{truth}{`theta` (natural scale), `beta`, `logit` and
#'       `prevalence` (m x T matrices of the surface without the nugget),
#'       and `x`, the logit at each surveyed cell with its nugget, in the
#'       row order of `data` (what a fit estimates for the observed cells)}
#'     \item{frame, time}{as given}
#'   }
#'   Pass `sim$data` to [areal_data()] with `area = "id"` and `t = "t"`.
#' @examples
#' fr <- areal_lattice(6, 6)
#' tf <- time_frame(2001, 2020, "year")
#' sim <- simulate_areal(fr, tf, n_cells = 60, seed = 1)
#' d <- areal_data(sim$data, fr, tf, y = "positive", k = "tested", area = "id", t = "t")
#' @export
simulate_areal <- function(frame, time, sigma = 1, rho_s = 0.9, rho_t = 0.9, sigma_n = 0.3,
                           beta = c(-1, 1), n_cells = 100, k = 50, seed = NULL) {
  if (!inherits(frame, "stareal_frame")) stop("`frame` must be made by areal_frame()", call. = FALSE)
  if (!inherits(time, "stareal_time")) stop("`time` must be made by time_frame()", call. = FALSE)
  stopifnot(sigma > 0, rho_s > 0, rho_s < 1, rho_t > 0, rho_t < 1, sigma_n >= 0, length(beta) == 2)
  m <- frame$m; T <- time$T
  if (n_cells > m * T) stop(sprintf("n_cells = %d exceeds the %d area-time cells", n_cells, m * T), call. = FALSE)
  if (!length(k) %in% c(1, n_cells) || any(k < 1) || any(k != round(k)))
    stop("`k` must be whole numbers >= 1, one or one per cell", call. = FALSE)
  with_seed(seed, {
    g <- leroux_g(frame$eig$values, 1 / sigma^2, rho_s)
    E <- ar1_filter(matrix(rnorm(m * T), m, T), rho_t)
    field <- (frame$eig$vectors * rep(sqrt(g), each = m)) %*% E
    logit <- field + rep(as.vector(design_X(seq_len(T), T) %*% beta), each = m)
    cells <- sort(sample.int(m * T, n_cells))          # cell index (t - 1) * m + a
    a <- (cells - 1L) %% m + 1L; t <- (cells - 1L) %/% m + 1L
    kk <- rep_len(as.integer(k), n_cells)
    x <- logit[cbind(a, t)] + rnorm(n_cells, 0, sigma_n)
    pos <- rbinom(n_cells, kk, plogis(x))
  })
  data <- data.frame(id = frame$ids[a], t = t, time = time$labels[t], stringsAsFactors = FALSE)
  if (!is.null(time$starts)) data$date <- time$starts[t]
  data$tested <- kk; data$positive <- as.integer(pos)
  data$.x <- x                                       # the truth at each surveyed cell, nugget included
  data <- data[order(data$t, data$id), ]
  rownames(data) <- NULL
  x <- data$.x; data$.x <- NULL
  dimnames(logit) <- list(frame$ids, time$labels)
  structure(list(data = data,
                 truth = list(theta = c(sigma = sigma, rho_s = rho_s, rho_t = rho_t, sigma_n = sigma_n),
                              beta = c(intercept = beta[1], trend = beta[2]),
                              logit = logit, prevalence = plogis(logit), x = x),
                 frame = frame, time = time),
            class = "stareal_sim")
}

#' @export
print.stareal_sim <- function(x, ...) {
  th <- x$truth$theta
  cat(sprintf("<stareal_sim> %d surveys on %d units x %d time steps\n", nrow(x$data), x$frame$m, x$time$T))
  cat(sprintf("  truth: sigma %.3g, rho_s %.3g, rho_t %.3g, sigma_n %.3g; intercept %.3g, trend %.3g\n",
              th[["sigma"]], th[["rho_s"]], th[["rho_t"]], th[["sigma_n"]], x$truth$beta[[1]], x$truth$beta[[2]]))
  cat(sprintf("  observed prevalence: median %.3f, range %.3f-%.3f\n", stats::median(x$data$positive / x$data$tested),
              min(x$data$positive / x$data$tested), max(x$data$positive / x$data$tested)))
  invisible(x)
}
