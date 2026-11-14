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
#' time step to 1 at the last) + a spatial x AR1 field (the spatial prior
#' set by `spatial`, as in [fit_areal()]), with a cell-level nugget added at
#' each surveyed cell.
#'
#' @param frame A [areal_frame()].
#' @param time A [time_frame()].
#' @param formula,covariates,scale_covariates The fixed effects, as in
#'   [fit_areal()]. Simulating and fitting with the same settings gives
#'   effects on the same scale.
#' @param beta Fixed effects on the logit scale: one per design column, in
#'   order (intercept, `trend` if in the formula, then the covariates; see
#'   [areal_priors()]), or named by design column.
#' @param sigma Scale of the field. For `"leroux"` and `"car"`, the
#'   standard deviation is `sigma` times a factor set by the graph
#'   (\eqn{\tau = 1/\sigma^2}); for `"bym2"`, roughly the marginal standard
#'   deviation.
#' @param rho_s Spatial dependence, between 0 and 1 (below 1 for `"car"`;
#'   at 0 the areas are independent, at 1 the Leroux prior is the ICAR).
#' @param spatial The spatial prior: `"leroux"`, `"bym2"` or `"car"` (see
#'   "Spatial priors" in [fit_areal()]).
#' @param rho_t Lag-one temporal correlation of the AR1, at least 0 and below 1.
#' @param sigma_n Standard deviation of the cell-level nugget (0 for none).
#' @param main,interaction,temporal The components of the field, as in
#'   [fit_areal()]: main effects in space and/or time (`"none"`, `"space"`,
#'   `"time"`, `"both"`), the type of space-time interaction (`"IV"`,
#'   `"I"`, `"II"`, `"III"`, `"none"`), and the prior of the time main
#'   effect (`"ar1"`, `"rw1"`, `"rw2"`). The defaults give the separable
#'   field described above. `sigma`, `rho_s` and `rho_t` belong to the
#'   interaction (each only where its type has it).
#' @param sigma_space,rho_space Scale and spatial dependence of the space
#'   main effect (under the prior `spatial`).
#' @param sigma_time,rho_time Scale of the time main effect (for a random
#'   walk, roughly its marginal standard deviation) and, for `"ar1"`, its
#'   lag-one correlation.
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
#'   [fit_areal()], with `$truth`: `spatial`, `theta` (natural scale; for a
#'   model with main effects or another interaction type, only the
#'   hyperparameters it has), `beta` (named
#'   by design column), `logit` and `prevalence` (areas x time steps
#'   matrices of the surface without the nugget; prevalence in percent),
#'   `x` (the logit at each surveyed cell with its nugget, in the order of
#'   `$cells`: what a fit
#'   estimates for the observed cells) and `design`. `$surveys` has columns
#'   `id`, `t`, `time`, `date` (calendar frames), `y`, `k`. A simulation
#'   with main effects or another interaction type also has, in `$truth`,
#'   `main`, `interaction`, `temporal`, `space` (the space main effect, one
#'   value per area, or `NULL`) and `time` (the time main effect, one value
#'   per step, or `NULL`).
#' @examples
#' fr <- areal_lattice(6, 6)
#' tf <- time_frame(2001, 2020, "year")
#' sim <- simulate_areal(fr, tf, n_cells = 60, seed = 1)
#' sim
#' @export
simulate_areal <- function(frame, time, formula = ~ trend(), covariates = NULL, scale_covariates = TRUE,
                           beta = c(-1, 1), sigma = 1, rho_s = 0.9, rho_t = 0.9, sigma_n = 0.3,
                           n_cells = 100, k = 50, seed = NULL, cells = NULL, spatial = "leroux",
                           main = "none", interaction = "IV", temporal = "ar1",
                           sigma_space = 0.5, rho_space = 0.9, sigma_time = 0.5, rho_time = 0.9) {
  if (!inherits(frame, "stareal_frame")) stop("`frame` must be made by areal_frame()", call. = FALSE)
  if (!inherits(time, "stareal_time")) stop("`time` must be made by time_frame()", call. = FALSE)
  scalar <- function(x, nm, lo, hi, lo_ok = FALSE) {
    if (!is.numeric(x) || length(x) != 1 || !is.finite(x) || x < lo || (!lo_ok && x == lo) || x >= hi)
      stop(sprintf("`%s` must be one number in %s%g, %g)", nm, if (lo_ok) "[" else "(", lo, hi), call. = FALSE)
  }
  if (!is.character(spatial) || length(spatial) != 1 || !spatial %in% c("leroux", "bym2", "car"))
    stop("`spatial` must be \"leroux\", \"bym2\" or \"car\" (for iid or ICAR, \"leroux\" with rho_s 0 or 1)",
         call. = FALSE)
  scalar(sigma, "sigma", 0, Inf); scalar(rho_t, "rho_t", 0, 1, lo_ok = TRUE)
  if (!is.numeric(rho_s) || length(rho_s) != 1 || !is.finite(rho_s) || rho_s < 0 || rho_s > 1 ||
      (spatial == "car" && rho_s == 1))
    stop(sprintf("`rho_s` must be one number between 0 and 1%s", if (spatial == "car") ", below 1 for \"car\"" else ""),
         call. = FALSE)
  basis <- spatial_basis(frame, spatial)
  one_of <- function(x, what, choices)
    if (!is.character(x) || length(x) != 1 || !x %in% choices)
      stop(sprintf("`%s` must be one of %s", what, paste(sprintf('"%s"', choices), collapse = ", ")), call. = FALSE)
  one_of(main, "main", main_choices); one_of(interaction, "interaction", interaction_choices)
  one_of(temporal, "temporal", temporal_choices)
  additive <- is_additive(list(main = main, interaction = interaction))
  has_space <- main %in% c("space", "both"); has_time <- main %in% c("time", "both")
  if (temporal != "ar1" && !has_time)
    stop("`temporal` is the prior of the time main effect: use it with main = \"time\" or \"both\"", call. = FALSE)
  if (interaction == "none" && main == "none")
    stop("interaction = \"none\" needs a main effect (main = \"space\", \"time\" or \"both\")", call. = FALSE)
  if (has_time && time$n_steps < 2) stop("a time main effect needs more than one time step", call. = FALSE)
  if (has_space) {
    scalar(sigma_space, "sigma_space", 0, Inf)
    if (!is.numeric(rho_space) || length(rho_space) != 1 || !is.finite(rho_space) || rho_space < 0 ||
        rho_space > 1 || (spatial == "car" && rho_space == 1))
      stop(sprintf("`rho_space` must be one number between 0 and 1%s",
                   if (spatial == "car") ", below 1 for \"car\"" else ""), call. = FALSE)
  }
  if (has_time) {
    scalar(sigma_time, "sigma_time", 0, Inf)
    if (temporal == "ar1") scalar(rho_time, "rho_time", 0, 1, lo_ok = TRUE)
  }
  uses_spatial <- !additive || interaction %in% c("III", "IV") || has_space
  if (uses_spatial && (spatial == "bym2" || (rho_s == 1 && (!additive || interaction %in% c("III", "IV"))) ||
                       (has_space && rho_space == 1)) && frame$n_components > 1)
    stop("this spatial prior needs a connected graph; link the frame's components (see areal_frame())", call. = FALSE)
  if (uses_spatial && spatial == "car" && any(Matrix::rowSums(frame$W) == 0))
    stop("the proper CAR needs every area to have a neighbour", call. = FALSE)
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
    if (!additive) {
      g <- spatial_g(basis, list(sigma = sigma, tau = 1 / sigma^2, rho_s = rho_s))
      E <- ar1_filter(matrix(rnorm(m * T), m, T), rho_t)
      field <- (basis$vectors * rep(sqrt(g), each = m)) %*% E
    } else {
      # the interaction (its space factor iid for types I and II, its time
      # factor iid for I and III, absent for "none"), then the main effects,
      # whose noise comes after the interaction's
      struct_s <- interaction %in% c("III", "IV"); struct_t <- interaction %in% c("II", "IV")
      b_int <- if (struct_s) basis else frame$eig
      g <- if (interaction == "none") numeric(m) else if (struct_s)
        spatial_g(basis, list(sigma = sigma, tau = 1 / sigma^2, rho_s = rho_s)) else
          leroux_g(frame$eig$values, 1 / sigma^2, 0)
      E <- ar1_filter(matrix(rnorm(m * T), m, T), if (struct_t) rho_t else 0)
      field <- (b_int$vectors * rep(sqrt(g), each = m)) %*% E
      mz <- rnorm(m + T)
      f_space <- f_time <- NULL
      if (has_space) {
        gS <- spatial_g(basis, list(sigma = sigma_space, tau = 1 / sigma_space^2, rho_s = rho_space))
        f_space <- space_draw(basis, gS, mz[seq_len(m)])
        field <- field + f_space
      }
      if (has_time) {
        f_time <- time_draw(time_basis(T, temporal), list(sigma_time = sigma_time, rho_time = rho_time),
                            mz[m + seq_len(T)])
        field <- field + rep(f_time, each = m)
      }
    }
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
  d$truth <- list(spatial = spatial, theta = c(sigma = sigma, rho_s = rho_s, rho_t = rho_t, sigma_n = sigma_n),
                  beta = stats::setNames(as.numeric(beta), des$names),
                  logit = logit, prevalence = 100 * stats::plogis(logit), x = x[o], design = des)
  if (additive) {
    all <- model_theta_all(list(main = main, interaction = interaction, temporal = temporal), sigma_n > 0)
    th <- c(sigma = sigma, rho_s = rho_s, rho_t = rho_t, sigma_n = sigma_n, sigma_space = sigma_space,
            rho_space = rho_space, sigma_time = sigma_time, rho_time = rho_time)
    d$truth$theta <- th[unname(theta_natural_names[all])]
    d$truth <- c(d$truth, list(main = main, interaction = interaction, temporal = temporal,
                               space = if (!is.null(f_space)) stats::setNames(f_space, frame$ids),
                               time = if (!is.null(f_time)) stats::setNames(f_time, time$labels)))
  }
  d
}
