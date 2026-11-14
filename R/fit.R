# fit.R
#
# fit_areal(): the user-facing fit, a thin layer over fit_setup() and
# run_chains() (chains.R); the fit object's print and plot methods and its
# posterior draws (posterior::as_draws*). summary() and diagnose() are in
# diagnose.R.
#
# ------------------------------------------------------------------

#' Fit the space-time areal model
#'
#' Fits the model to the observed area-time cells of `data` by MCMC. The
#' logit prevalence of cell (area \eqn{a}, time \eqn{t}) is
#' \deqn{x_{at} = z_{at}^\top \beta + f_{at} + \epsilon_{at},}
#' where \eqn{z_{at}} holds the fixed effects set by `formula` (by default
#' an intercept and a linear trend in time), \eqn{f} is a Gaussian field
#' with separable covariance (a conditional autoregressive prior on the
#' frame's graph in space, with scale `sigma` and dependence `rho_s`; a
#' unit-variance AR1 with correlation `rho_t` in time), and \eqn{\epsilon}
#' is a cell-level nugget with standard deviation `sigma_n`. Positives are
#' binomial given \eqn{x}.
#'
#' @section Spatial priors: `spatial` chooses the prior on the graph. In
#'   each, `rho_s = 0` makes the areas independent and `rho_s` near 1 makes
#'   neighbours alike (the intrinsic CAR, ICAR, at 1). \eqn{D} is the
#'   diagonal matrix of each area's number of neighbours and \eqn{W} the
#'   adjacency, so \eqn{R = D - W} is the graph Laplacian.
#'   \describe{
#'     \item{`"leroux"`}{(the default; Leroux et al. 2000) precision
#'       \eqn{\tau(\rho_s R + (1 - \rho_s) I)}, \eqn{\sigma^2 = 1/\tau}.}
#'     \item{`"bym2"`}{(Riebler et al. 2016) covariance
#'       \eqn{\sigma^2(\rho_s R_*^+ + (1 - \rho_s) I)}: a mix of a scaled
#'       ICAR (\eqn{R_*^+}, the generalised inverse of \eqn{R} scaled so
#'       that the geometric mean of its variances is 1) and independent
#'       noise. `sigma` is then roughly the marginal sd and `rho_s` the
#'       share of the variance that is spatially structured. Needs a
#'       connected graph.}
#'     \item{`"car"`}{the proper CAR, precision \eqn{\tau(D - \rho_s W)}.
#'       Its marginal variances differ between areas (larger where there
#'       are fewer neighbours). `rho_s` must stay below 1, and every area
#'       needs a neighbour.}
#'     \item{`"iid"`}{independent areas: the Leroux prior with `rho_s`
#'       fixed at 0.}
#'     \item{`"icar"`}{the intrinsic CAR: the Leroux prior with `rho_s`
#'       fixed at 1. The field then sums to zero over the areas at every
#'       time step, so the regional level and its changes over time come
#'       from the fixed effects alone. Needs a connected graph.}
#'   }
#'   Any hyperparameter can also be held at a value with `fixed`, e.g.
#'   `fixed = list(rho_t = 0)` for a field with no memory from one time step
#'   to the next. Fixed hyperparameters are not sampled and do not appear in
#'   [summary()][summary.stareal_fit] or [diagnose()].
#'
#' @section Main effects and interaction types: By default the field
#'   \eqn{f} is one separable space-time component. `main` and
#'   `interaction` make it a sum of independent components (Knorr-Held
#'   2000, *Statistics in Medicine* 19, 2555-2567):
#'   \deqn{f_{at} = f^{ST}_{at} + f^S_a + f^T_t.}
#'   \describe{
#'     \item{`main`}{`"none"` (the default), `"space"` (a spatial effect
#'       \eqn{f^S}, the same at every time step, under the prior `spatial`,
#'       with its own scale `sigma_space` and dependence `rho_space`),
#'       `"time"` (a temporal effect \eqn{f^T}, the same in every area,
#'       with scale `sigma_time`) or `"both"`.}
#'     \item{`temporal`}{the prior of the time main effect: `"ar1"` (the
#'       default; lag-one correlation `rho_time`), `"rw1"` or `"rw2"` (first-
#'       or second-order random walks). The random walks are scaled so that
#'       `sigma_time` is roughly the marginal standard deviation whatever the
#'       number of time steps (Sorbye and Rue 2014). An RW1 sums to zero over
#'       time and an RW2 also has zero slope, so the intercept and `trend()`
#'       carry the level and the straight line, and the RW2 the curve about
#'       it: the natural pair. An RW1 or AR1 time effect overlaps `trend()`
#'       in part, which is not an error but weakens what the trend
#'       coefficient means. Random walks are for the time main effect only.}
#'     \item{`interaction`}{the space-time interaction \eqn{f^{ST}}, with scale
#'       `sigma`, by whether its time and space factors are structured:
#'       `"IV"` (the default; AR1 in time with `rho_t` x the spatial prior
#'       with `rho_s`), `"III"` (independent in time x the spatial prior),
#'       `"II"` (AR1 in time x independent areas), `"I"` (independent cells)
#'       or `"none"`. Type I needs `nugget = FALSE`: on the observed cells
#'       an independent cell effect is the nugget again. Type III with a
#'       nugget comes close to that when `rho_s` is small, and then mixes
#'       slowly: check [diagnose()] and run longer.}
#'   }
#'   With the defaults the model is exactly the separable one above. With
#'   the ICAR, BYM2 or a random walk a component sums to zero, so the
#'   regional level always comes from the fixed effects. More components
#'   mean more hyperparameters, each less well identified than in a simpler
#'   model: expect to need longer runs, and compare models with
#'   [loo()][loo.stareal_fit] or [crossval()].
#'
#' The unobserved cells and the fixed effects are integrated out exactly, so
#' the sampler works on the observed cells only; its cost per iteration
#' grows as the cube of the number of observed cells. The full area-time
#' surface is reconstructed afterwards by [predict()][predict.stareal_fit].
#'
#' @param data A [areal_data()] object (or a simulation from
#'   [simulate_areal()]).
#' @param formula A one-sided formula for the fixed effects. It always has
#'   an intercept; `trend()` adds a linear trend in time (from -1 at the
#'   first time step to 1 at the last); other terms are names of columns of
#'   `covariates`. Factor or character covariates get one dummy per level
#'   after the first. Transforms and interactions must be made as columns.
#' @param covariates A data frame, or a list of data frames, of covariates.
#'   Each is keyed by an `id` column (areas, for covariates that do not
#'   change over time), a `t` (step number) or `time` (label) column (for
#'   covariates the same in every area), or both (covariates by area and
#'   time). Each must cover every area, time step or area-time cell it is
#'   keyed by, since predictions are made everywhere.
#' @param scale_covariates Standardise numeric covariates to mean 0 and
#'   sd 1, over the areas, time steps or cells they are given for, so that
#'   the default prior on their effects (sd 3, see [areal_priors()]) is
#'   weakly informative whatever their units. Effects are then per sd of
#'   the covariate; the centres and scales are in `fit$design` and in
#'   `summary()`. Unscaled covariates in large units (e.g. population
#'   counts) can make the prior variance of the mean so large that the
#'   field is lost in rounding; the fit then stops and says so.
#' @param nugget Include the cell-level nugget.
#' @param spatial The spatial prior: `"leroux"`, `"bym2"`, `"car"`, `"iid"`
#'   or `"icar"` (see "Spatial priors"), for every component that is
#'   structured in space. `"iid"` and `"icar"` fix `rho_space` too.
#' @param main Main effects: `"none"`, `"space"`, `"time"` or `"both"` (see
#'   "Main effects and interaction types").
#' @param interaction The space-time interaction: `"IV"`, `"III"`, `"II"`,
#'   `"I"` or `"none"`.
#' @param temporal The prior of the time main effect: `"ar1"`, `"rw1"` or
#'   `"rw2"`.
#' @param fixed Hyperparameters held at given values instead of sampled: a
#'   named list (or vector) with any of the model's hyperparameters:
#'   `sigma`, `rho_s` (between 0 and 1; below 1 for `"car"`), `rho_t` (at
#'   least 0 and below 1), `sigma_n` (with `nugget = TRUE`), and with main
#'   effects `sigma_space`, `rho_space`, `sigma_time` and `rho_time`. At
#'   least one hyperparameter must stay free.
#' @param priors [areal_priors()], or a list of some of its arguments.
#' @param control [mcmc_control()], or a list of some of its arguments.
#' @param checkpoint_dir A folder to save each chain's state in. If a fit is
#'   interrupted, running it again with the same folder resumes every chain
#'   where it stopped, and gives exactly the draws of an uninterrupted run;
#'   a finished fit run again is read back from the folder. A folder left
#'   by a different fit (other data, frame, model or sampler settings) is
#'   refused. Use a folder outside synced storage (Dropbox, OneDrive), which
#'   can lock files. `NULL`: no checkpoints.
#' @return An object of class `stareal_fit`:
#'   \describe{
#'     \item{trace}{the hyperparameter trace, burn-in included: chain,
#'       iteration, the free hyperparameters on the unconstrained scale
#'       (`log_sigma`, `logit_rho_s`, `logit_rho_t`, `log_sigma_n`, and
#'       with main effects `log_sigma_space`, `logit_rho_space`,
#'       `log_sigma_time`, `logit_rho_time`), `lp`
#'       (the log prior of the hyperparameters plus the log density of the
#'       latent values given them; `lp + ll` is the unnormalised log
#'       joint), the log-likelihood `ll`, and the acceptance rates of the
#'       two hyperparameter moves}
#'     \item{draws}{the stored posterior draws, as arrays (draw x chain x
#'       variable): `theta` (unconstrained), `beta`, and `x` (logit
#'       prevalence of each observed cell, nugget included, in the order of
#'       `data$cells`)}
#'     \item{stored_iter, seconds}{the iterations the draws were stored at;
#'       run time per chain}
#'     \item{design}{the fixed-effect design: formula, names, and the
#'       centre and scale of each standardised covariate}
#'     \item{model}{the model: `spatial` (as given), `family` (the prior it
#'       is: `"leroux"`, `"bym2"` or `"car"`), `fixed` (the fixed
#'       hyperparameters, natural scale, named), `main`, `interaction` and
#'       `temporal`}
#'     \item{data, priors, control, nugget, checkpoint_dir, call}{as given
#'       (priors with one prior per fixed effect)}
#'     \item{version}{what made the fit: the *stareal* version, its commit
#'       (`sha`, when installed from GitHub), the sampler version, R and
#'       the BLAS library}
#'   }
#'   Use [summary()][summary.stareal_fit], [diagnose()],
#'   [plot()][plot.stareal_fit], [predict()][predict.stareal_fit],
#'   [loo()][loo.stareal_fit], [compare_models()] and
#'   `posterior::as_draws_df(fit)` (or `as_draws_array()`, ...).
#' @examples
#' fr <- areal_lattice(4, 4)
#' tf <- time_frame(n_steps = 5)
#' sim <- simulate_areal(fr, tf, n_cells = 30, seed = 1)
#' fit <- fit_areal(sim, control = list(n_iter = 100, n_burn = 50, n_store = 25, n_chains = 1,
#'                                      ess_sweeps = 2, verbose = FALSE))
#' fit
#'
#' # the same data under a BYM2 prior, and with no memory in time
#' fit_bym2 <- fit_areal(sim, spatial = "bym2", control = list(n_iter = 100, n_burn = 50, n_store = 25,
#'                                                            n_chains = 1, ess_sweeps = 2, verbose = FALSE))
#' fit_iid_t <- fit_areal(sim, fixed = list(rho_t = 0), control = list(n_iter = 100, n_burn = 50, n_store = 25,
#'                                                                    n_chains = 1, ess_sweeps = 2, verbose = FALSE))
#'
#' # a regional curve over time (an RW2 about the trend) plus the space-time field
#' fit_rw2 <- fit_areal(sim, main = "time", temporal = "rw2",
#'                      control = list(n_iter = 100, n_burn = 50, n_store = 25, n_chains = 1, ess_sweeps = 2,
#'                                     verbose = FALSE))
#' @export
fit_areal <- function(data, formula = ~ trend(), covariates = NULL, scale_covariates = TRUE, nugget = TRUE,
                      spatial = "leroux", main = "none", interaction = "IV", temporal = "ar1",
                      fixed = NULL, priors = areal_priors(), control = mcmc_control(),
                      checkpoint_dir = NULL) {
  if (!inherits(data, "stareal_data")) stop("`data` must be made by areal_data() or simulate_areal()", call. = FALSE)
  model <- model_spec(spatial, fixed, nugget, data$frame, main = main, interaction = interaction,
                      temporal = temporal, n_steps = data$time$n_steps)
  priors <- complete_settings(priors, areal_priors, "priors")
  control <- complete_settings(control, mcmc_control, "control")
  if (!is.null(checkpoint_dir) && (!is.character(checkpoint_dir) || length(checkpoint_dir) != 1))
    stop("`checkpoint_dir` must be one folder path", call. = FALSE)
  design <- build_design(formula, covariates, data$frame, data$time, scale = scale_covariates)
  fit <- fit_with_design(data, design, nugget, priors, control, checkpoint_dir, model = model)
  fit$call <- match.call()
  fit
}

# ---- the spatial prior and fixed hyperparameters -----------------------------

spatial_choices <- c("leroux", "bym2", "car", "iid", "icar")

# The model as fit_areal() was asked for it, checked: spatial (as given),
# family (leroux, bym2 or car) and fixed (natural scale, named, in the order
# sigma, rho_s, rho_t, sigma_n).
model_spec <- function(spatial = "leroux", fixed = NULL, nugget = TRUE, frame = NULL,
                       main = "none", interaction = "IV", temporal = "ar1", n_steps = NULL) {
  one_of <- function(x, what, choices)
    if (!is.character(x) || length(x) != 1 || !x %in% choices)
      stop(sprintf("`%s` must be one of %s", what, paste(sprintf('"%s"', choices), collapse = ", ")), call. = FALSE)
  one_of(spatial, "spatial", spatial_choices)
  one_of(main, "main", main_choices)
  one_of(interaction, "interaction", interaction_choices)
  one_of(temporal, "temporal", temporal_choices)
  md0 <- list(main = main, interaction = interaction, temporal = temporal)
  additive <- is_additive(md0)
  has_space <- main %in% c("space", "both"); has_time <- main %in% c("time", "both")
  if (temporal != "ar1" && !has_time)
    stop("`temporal` is the prior of the time main effect: use it with main = \"time\" or \"both\"", call. = FALSE)
  if (interaction == "none" && main == "none")
    stop("interaction = \"none\" needs a main effect (main = \"space\", \"time\" or \"both\")", call. = FALSE)
  if (interaction == "I" && isTRUE(nugget))
    stop(paste("interaction = \"I\" (independent cells) is the nugget again on the observed cells, so the two",
               "cannot be told apart: use nugget = FALSE with it"), call. = FALSE)
  if (!is.null(n_steps) && n_steps == 1) {
    if (has_time) stop("a time main effect needs more than one time step", call. = FALSE)
    if (has_space && interaction %in% c("III", "IV"))
      stop(paste("with one time step a space main effect duplicates a type III or IV interaction; use",
                 "interaction = \"I\", \"II\" or \"none\""), call. = FALSE)
  }
  if (!is.null(n_steps) && temporal == "rw2" && n_steps < 3)
    stop("an RW2 needs at least three time steps", call. = FALSE)
  if (is.null(fixed)) fixed <- list()
  if (!is.list(fixed) && !is.numeric(fixed)) stop("`fixed` must be a named list, e.g. list(rho_t = 0)", call. = FALSE)
  fixed <- as.list(fixed)
  nm <- names(fixed)
  if (length(fixed) && (is.null(nm) || any(nm == "") || anyDuplicated(nm)))
    stop("`fixed` must name each hyperparameter once, e.g. list(rho_t = 0)", call. = FALSE)
  # the hyperparameters this model has (natural-scale names)
  known <- if (additive) unname(theta_natural_names[model_theta_all(c(md0, list(family = "leroux")), nugget)]) else
    c("sigma", "rho_s", "rho_t", if (isTRUE(nugget)) "sigma_n")
  all_known <- unname(theta_natural_names)
  if (length(setdiff(nm, all_known)))
    stop(sprintf("`fixed` can hold %s, not %s", paste(known, collapse = ", "),
                 paste(setdiff(nm, all_known), collapse = ", ")), call. = FALSE)
  if ("sigma_n" %in% nm && !isTRUE(nugget)) stop("`fixed$sigma_n` needs nugget = TRUE", call. = FALSE)
  if (length(setdiff(nm, known)))
    stop(sprintf("this model has no %s to fix (its hyperparameters: %s)", paste(setdiff(nm, known), collapse = ", "),
                 paste(known, collapse = ", ")), call. = FALSE)
  for (v in nm) if (!is.numeric(fixed[[v]]) || length(fixed[[v]]) != 1 || !is.finite(fixed[[v]]))
    stop(sprintf("`fixed$%s` must be one number", v), call. = FALSE)
  family <- switch(spatial, iid = "leroux", icar = "leroux", spatial)
  if (spatial %in% c("iid", "icar")) {
    # the shortcut fixes the correlation of every component that uses the spatial prior
    for (v in intersect(c("rho_s", "rho_space"), known)) {
      if (v %in% nm) stop(sprintf("spatial = \"%s\" fixes %s already; leave it out of `fixed`", spatial, v),
                          call. = FALSE)
      fixed[[v]] <- if (spatial == "iid") 0 else 1
    }
  }
  fx <- unlist(fixed)
  if (is.null(fx)) fx <- stats::setNames(numeric(0), character(0))
  for (v in intersect(c("sigma", "sigma_space", "sigma_time"), names(fx)))
    if (fx[[v]] <= 0) stop(sprintf("`fixed$%s` must be positive", v), call. = FALSE)
  if ("sigma_n" %in% names(fx) && fx[["sigma_n"]] <= 0)
    stop("`fixed$sigma_n` must be positive (for no nugget, use nugget = FALSE)", call. = FALSE)
  for (v in intersect(c("rho_t", "rho_time"), names(fx)))
    if (fx[[v]] < 0 || fx[[v]] >= 1) stop(sprintf("`fixed$%s` must be at least 0 and below 1", v), call. = FALSE)
  for (v in intersect(c("rho_s", "rho_space"), names(fx))) {
    hi_ok <- family != "car"
    if (fx[[v]] < 0 || fx[[v]] > 1 || (!hi_ok && fx[[v]] == 1))
      stop(sprintf("`fixed$%s` must be between 0 and 1%s", v, if (hi_ok) "" else " and below 1 for the proper CAR"),
           call. = FALSE)
  }
  if (length(fx) >= length(known)) stop("`fixed` fixes every hyperparameter; leave at least one free", call. = FALSE)
  fx <- fx[intersect(all_known, names(fx))]
  if (!length(fx)) fx <- NULL
  uses_spatial <- !additive || interaction %in% c("III", "IV") || has_space
  if (!is.null(frame) && uses_spatial) {
    icar_like <- any(c("rho_s", "rho_space") %in% names(fx) & fx[c("rho_s", "rho_space")] %in% 1)
    if ((family == "bym2" || icar_like) && frame$n_components > 1)
      stop(sprintf("%s needs a connected graph, and this frame has %d components; link them (see areal_frame())",
                   if (family == "bym2") "the BYM2 prior" else "the ICAR (rho_s = 1)", frame$n_components),
           call. = FALSE)
    if (family == "car" && any(Matrix::rowSums(frame$W) == 0))
      stop("the proper CAR needs every area to have a neighbour; link the frame's islands (see areal_frame())",
           call. = FALSE)
  }
  list(spatial = spatial, family = family, fixed = fx, main = main, interaction = interaction, temporal = temporal)
}

# the fixed hyperparameters on the theta scale, named as theta is
fixed_theta <- function(fixed) {
  if (!length(fixed)) return(NULL)
  tr <- stats::setNames(names(theta_natural_names), theta_natural_names)
  stats::setNames(vapply(names(fixed), function(v) if (startsWith(v, "sigma")) log(fixed[[v]]) else
    stats::qlogis(fixed[[v]]), numeric(1)), tr[names(fixed)])
}

# a fit's model (fits made before `spatial` existed are Leroux, all free;
# fits made before the components existed have the default ones), and its
# basis
fit_model <- function(fit) complete_model(if (is.null(fit$model)) model_spec() else fit$model)
fit_basis <- function(fit) model_basis(fit$data$frame, fit$data$time$n_steps, fit_model(fit), fit$nugget)

# a one-line description of a model, e.g. "bym2, rho_t = 0", or for an
# additive model "leroux; main space + time (rw2); interaction III"
model_label <- function(md) {
  md <- complete_model(md)
  fx <- md$fixed
  if (md$spatial %in% c("iid", "icar")) fx <- fx[!names(fx) %in% c("rho_s", "rho_space")]
  if (md$interaction == "III") fx <- fx[names(fx) != "rho_t"]
  head <- md$spatial
  if (is_additive(md)) {
    tm <- if (md$temporal != "ar1" || md$main %in% c("time", "both")) sprintf("time (%s)", md$temporal)
    parts <- switch(md$main, none = NULL, space = "space", time = tm, both = paste("space +", tm))
    head <- paste(c(md$spatial, if (length(parts)) paste("main", parts), paste("interaction", md$interaction)),
                  collapse = "; ")
  }
  paste(c(head, if (length(fx)) sprintf("%s = %g", names(fx), fx)), collapse = ", ")
}

# what made a fit, prediction or cross-validation
version_info <- function() {
  d <- suppressWarnings(utils::packageDescription("stareal"))
  sha <- if (is.list(d) && !is.null(d$RemoteSha)) d$RemoteSha else NA_character_
  blas <- tryCatch(extSoftVersion()[["BLAS"]], error = function(e) "")
  if (!nzchar(blas)) blas <- tryCatch(La_library(), error = function(e) "")
  if (!nzchar(blas) && .Platform$OS.type == "windows")        # R for Windows reports neither; MKL replaces Rblas.dll
    blas <- if (length(list.files(R.home("bin"), pattern = "^mkl_rt"))) "mkl_rt (Rblas.dll)" else "Rblas.dll"
  list(stareal = as.character(utils::packageVersion("stareal")), sha = sha, sampler_version = SAMPLER_VERSION,
       r = R.version.string, blas = if (nzchar(blas)) basename(blas) else NA_character_)
}

# fit_areal() given a built design and completed settings (crossval()
# refits with the fit's own)
fit_with_design <- function(data, design, nugget, priors, control, checkpoint_dir = NULL, model = model_spec()) {
  cells <- data$cells; frame <- data$frame; T <- data$time$n_steps
  n <- nrow(cells); p <- length(design$names)
  if (n < 2) stop("at least two observed cells are needed", call. = FALSE)
  priors$beta_mean <- expand_beta_prior(priors$beta_mean, design$names, 0, "beta_mean")
  priors$beta_sd <- expand_beta_prior(priors$beta_sd, design$names, 3, "beta_sd")
  pv <- design_prior_var(design, priors$beta_sd)
  if (pv$var > 1e8)
    stop(sprintf(paste("the fixed effects' prior gives the logit prevalence a prior sd of up to %.3g (mostly from",
                       "`%s`), which swamps the field in rounding. Use scale_covariates = TRUE, or a smaller",
                       "beta_sd for `%s`."), sqrt(pv$var), pv$worst, pv$worst), call. = FALSE)
  X <- design_eval(design, cells$area, cells$t)
  q <- qr(X)
  if (q$rank < p)
    warning(sprintf(paste("the observed cells do not inform the fixed effect(s) %s (collinear with the others",
                          "over the observed cells, e.g. a factor level or time step with no data): their posterior",
                          "is their prior, and predictions where they apply carry that uncertainty"),
                    paste(colnames(X)[q$pivot[(q$rank + 1):p]], collapse = ", ")), call. = FALSE, immediate. = TRUE)
  setup <- fit_setup(cells, model_basis(frame, T, model, nugget), T, priors, nugget, X)

  # scale guard-rails: the eigendecomposition is m x m, each iteration O(n^3)
  if (n > 4000) warning(sprintf(paste("n = %d observed cells: each iteration costs O(n^3) and will be slow;",
                                      "use an optimised BLAS (MKL, OpenBLAS or Accelerate)"), n),
                        call. = FALSE, immediate. = TRUE)
  n_workers <- min(control$n_chains, control$n_cores)
  if (control$verbose)
    message(sprintf("fitting %d observed cells (%d units x %d time steps): %d chain%s x %d iterations, %s",
                    n, frame$m, T, control$n_chains, if (control$n_chains > 1) "s" else "", control$n_iter,
                    if (n_workers > 1) sprintf("in parallel on %d workers", n_workers) else "in this session"))

  res <- with_seed(control$seed, run_chains(setup, control, checkpoint_dir, frame$W))

  to_array <- function(l) {
    a <- array(NA_real_, c(nrow(l[[1]]), length(l), ncol(l[[1]])),
               dimnames = list(NULL, NULL, colnames(l[[1]])))
    for (i in seq_along(l)) a[, i, ] <- l[[i]]
    a
  }
  xs <- lapply(res, function(r) { x <- r$x; colnames(x) <- cells$cell_id; x })
  trace <- do.call(rbind, lapply(res, function(r) cbind(chain = r$chain, r$theta)))
  rownames(trace) <- NULL
  structure(list(trace = trace,
                 draws = list(theta = to_array(lapply(res, `[[`, "draw_theta")),
                              beta = to_array(lapply(res, `[[`, "beta")), x = to_array(xs)),
                 stored_iter = res[[1]]$stored_iter,
                 seconds = vapply(res, `[[`, numeric(1), "seconds"),
                 design = design, model = model, data = data, priors = priors, control = control,
                 nugget = isTRUE(nugget), checkpoint_dir = checkpoint_dir, call = NULL, version = version_info()),
            class = "stareal_fit")
}

# a prior given as one value, one per design column, or named by design
# column (others take `default`) -> one per design column, named
expand_beta_prior <- function(v, names, default, what) {
  p <- length(names)
  if (!is.null(names(v))) {
    bad <- setdiff(names(v), names)
    if (length(bad)) stop(sprintf("`%s` names unknown design column(s) %s; the columns are %s", what,
                                  paste(bad, collapse = ", "), paste(names, collapse = ", ")), call. = FALSE)
    out <- rep(default, p); out[match(names(v), names)] <- v
  } else if (length(v) == 1) {
    out <- rep(v, p)
  } else if (length(v) == p) {
    out <- v
  } else {
    stop(sprintf("`%s` needs 1 or %d values, or values named by design column (%s)", what, p,
                 paste(names, collapse = ", ")), call. = FALSE)
  }
  stats::setNames(as.numeric(out), names)
}

# the sampler's setup for a fit, rebuilt (deterministically) from the fit
fit_setup_of <- function(fit) {
  cells <- fit$data$cells
  fit_setup(cells, fit_basis(fit), fit$data$time$n_steps, fit$priors, fit$nugget,
            design_eval(fit$design, cells$area, cells$t))
}

# the names of a fit's free hyperparameters (theta scale)
fit_theta_names <- function(fit) dimnames(fit$draws$theta)[[3]]

# ---- the fixed-effect design --------------------------------------------------

# The design in the two forms the code needs: rows of the full table (for
# predict_rows: (length(rows) * T) x p, rows in the order
# (t - 1) * length(rows) + i), and arbitrary cells (area a, step t; for
# held-out cells). The closure holds only the design, so it is small to send
# to workers.
design_rows <- function(des) {
  force(des)
  function(rows) design_eval(des, rep(rows, des$T), rep(seq_len(des$T), each = length(rows)))
}
design_cells <- function(fit, a, t) design_eval(fit$design, a, t)

# ---- draws ------------------------------------------------------------------

# the hyperparameter trace after burn-in, natural scale, as an
# iterations x chains x variables array (the full trace, not only the stored
# draws: more draws for theta's R-hat and ESS)
trace_array <- function(fit) {
  tn <- fit_theta_names(fit)
  post <- fit$trace[fit$trace$iteration > fit$control$n_burn, ]
  ch <- sort(unique(post$chain))
  l <- lapply(ch, function(c) theta_natural(as.matrix(post[post$chain == c, tn, drop = FALSE])))
  a <- array(NA_real_, c(nrow(l[[1]]), length(l), ncol(l[[1]])), dimnames = list(NULL, NULL, colnames(l[[1]])))
  for (i in seq_along(l)) a[, i, ] <- l[[i]]
  a
}

# the stored draws: hyperparameters (natural scale), fixed effects, and
# optionally the latent logits of the observed cells
stored_array <- function(fit, latent = TRUE) {
  d <- fit$draws
  dn <- dim(d$theta)
  nat <- theta_natural(matrix(d$theta, prod(dn[1:2]), dn[3], dimnames = list(NULL, dimnames(d$theta)[[3]])))
                                                                   # (draw, chain) rows, draw fastest
  th <- array(nat, c(dn[1:2], ncol(nat)), dimnames = list(NULL, NULL, colnames(nat)))
  parts <- list(th, d$beta)
  if (latent) {
    x <- d$x
    dimnames(x)[[3]] <- sprintf("logit[%s]", dimnames(x)[[3]])
    parts <- c(parts, list(x))
  }
  out <- array(NA_real_, c(dim(th)[1:2], sum(vapply(parts, function(a) dim(a)[3], 1))),
               dimnames = list(NULL, NULL, unlist(lapply(parts, function(a) dimnames(a)[[3]]))))
  k <- 0
  for (a in parts) { out[, , k + seq_len(dim(a)[3])] <- a; k <- k + dim(a)[3] }
  out
}

#' Posterior draws of a fit
#'
#' The stored draws as a `posterior` draws object: the free
#' hyperparameters (`sigma`, `rho_s`, `rho_t`, `sigma_n`), the fixed effects, and, if
#' `latent`, the logit prevalence of each observed cell, nugget included,
#' named `logit[<cell_id>]`. All the `posterior` conversions work:
#' `as_draws_df()`, `as_draws_array()`, `as_draws_matrix()`,
#' `as_draws_list()`, `as_draws_rvars()`.
#'
#' @param x A `stareal_fit`.
#' @param latent Include the observed cells' latent values.
#' @param ... Not used.
#' @return A `posterior` draws object.
#' @name as_draws.stareal_fit
#' @exportS3Method posterior::as_draws
as_draws.stareal_fit <- function(x, latent = TRUE, ...) as_draws_array.stareal_fit(x, latent = latent)

#' @rdname as_draws.stareal_fit
#' @exportS3Method posterior::as_draws_array
as_draws_array.stareal_fit <- function(x, latent = TRUE, ...) {
  a <- stored_array(x, latent)
  dimnames(a) <- list(iteration = NULL, chain = NULL, variable = dimnames(a)[[3]])
  posterior::as_draws_array(a)
}

#' @rdname as_draws.stareal_fit
#' @exportS3Method posterior::as_draws_df
as_draws_df.stareal_fit <- function(x, latent = TRUE, ...) posterior::as_draws_df(as_draws_array.stareal_fit(x, latent))

#' @rdname as_draws.stareal_fit
#' @exportS3Method posterior::as_draws_matrix
as_draws_matrix.stareal_fit <- function(x, latent = TRUE, ...) posterior::as_draws_matrix(as_draws_array.stareal_fit(x, latent))

#' @rdname as_draws.stareal_fit
#' @exportS3Method posterior::as_draws_list
as_draws_list.stareal_fit <- function(x, latent = TRUE, ...) posterior::as_draws_list(as_draws_array.stareal_fit(x, latent))

#' @rdname as_draws.stareal_fit
#' @exportS3Method posterior::as_draws_rvars
as_draws_rvars.stareal_fit <- function(x, latent = TRUE, ...) posterior::as_draws_rvars(as_draws_array.stareal_fit(x, latent))

# ---- print and plot ---------------------------------------------------------

#' @export
print.stareal_fit <- function(x, ...) {
  ct <- x$control
  cat(sprintf("<stareal_fit> %d observed cells on %d units x %d time steps%s\n", nrow(x$data$cells),
              x$data$frame$m, x$data$time$n_steps, if (x$nugget) ", with nugget" else ""))
  cat(sprintf("  %s: %s\n", if (is_additive(fit_model(x))) "model" else "spatial prior", model_label(fit_model(x))))
  cat(sprintf("  formula %s: %s\n", paste(deparse(x$design$formula, width.cutoff = 500L), collapse = " "),
              paste(x$design$names, collapse = ", ")))
  cat(sprintf("  %d chain%s x %d iterations (%d burn-in); %d stored draws in all; %s per chain\n",
              ct$n_chains, if (ct$n_chains > 1) "s" else "", ct$n_iter, ct$n_burn, ct$n_store * ct$n_chains,
              format_seconds(mean(x$seconds))))
  ta <- trace_array(x)
  med <- apply(ta, 3, stats::median)
  cat("  posterior medians: ", paste(sprintf("%s %.3g", names(med), med), collapse = ", "), "\n", sep = "")
  cat("  see summary(), diagnose(), plot() and predict()\n")
  invisible(x)
}

#' Trace plots of a fit
#'
#' One panel per parameter in [summary()][summary.stareal_fit], one colour
#' per chain. The hyperparameters are drawn from the full trace on the
#' natural scale, burn-in included (left of the dashed line). The fixed
#' effects are drawn from the stored draws, which start after burn-in.
#'
#' @param x A `stareal_fit`.
#' @param ... Passed to [graphics::matplot()] (e.g. `col`, `lty`, `lwd`).
#' @return `x`, invisibly.
#' @export
plot.stareal_fit <- function(x, ...) {
  panels <- trace_panels(x)
  n_col <- if (length(panels) <= 6) 2 else 3
  op <- graphics::par(mfrow = c(ceiling(length(panels) / n_col), n_col), mar = c(4, 4, 1.5, 1))
  on.exit(graphics::par(op))
  cols <- grDevices::hcl.colors(max(2, ncol(panels[[1]]$y)), "Dark 3")[seq_len(ncol(panels[[1]]$y))]
  for (v in names(panels)) {
    args <- utils::modifyList(list(type = "l", lty = 1, col = cols, xlab = "iteration", ylab = v, main = v),
                              list(...))
    do.call(graphics::matplot, c(list(panels[[v]]$x, panels[[v]]$y), args))
    if (panels[[v]]$burn_in) graphics::abline(v = x$control$n_burn, lty = 2, col = "grey40")
  }
  invisible(x)
}

# What plot() draws, one panel per parameter of summary(): x (iterations),
# y (iterations x chains), and whether the burn-in is included. The
# hyperparameters come from the full trace, the fixed effects from the
# stored draws.
trace_panels <- function(x) {
  tn <- fit_theta_names(x)
  ch <- sort(unique(x$trace$chain))
  nat <- theta_natural(as.matrix(x$trace[tn]))
  it <- x$trace$iteration[x$trace$chain == ch[1]]
  out <- list()
  for (v in colnames(nat))
    out[[v]] <- list(x = it, y = matrix(vapply(ch, function(c) nat[x$trace$chain == c, v], numeric(length(it))),
                                        nrow = length(it)), burn_in = TRUE)
  beta <- x$draws$beta                                     # draw x chain x fixed effect
  for (j in seq_len(dim(beta)[3]))
    out[[dimnames(beta)[[3]][j]]] <- list(x = x$stored_iter, y = matrix(beta[, , j], nrow = dim(beta)[1]),
                                          burn_in = FALSE)
  out
}
