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
#' with separable covariance (a Leroux prior on the frame's graph in space,
#' with scale `sigma` and dependence `rho_s`; a unit-variance AR1 with
#' correlation `rho_t` in time), and \eqn{\epsilon} is a cell-level nugget
#' with standard deviation `sigma_n`. Positives are binomial given \eqn{x}.
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
#'       iteration, the hyperparameters on the unconstrained scale
#'       (`log_sigma`, `logit_rho_s`, `logit_rho_t`, `log_sigma_n`), `lp`
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
#'     \item{data, priors, control, nugget, checkpoint_dir, call}{as given
#'       (priors with one prior per fixed effect)}
#'   }
#'   Use [summary()][summary.stareal_fit], [diagnose()],
#'   [plot()][plot.stareal_fit], [predict()][predict.stareal_fit] and
#'   `posterior::as_draws_df(fit)` (or `as_draws_array()`, ...).
#' @examples
#' fr <- areal_lattice(4, 4)
#' tf <- time_frame(n_steps = 5)
#' sim <- simulate_areal(fr, tf, n_cells = 30, seed = 1)
#' fit <- fit_areal(sim, control = list(n_iter = 100, n_burn = 50, n_store = 25, n_chains = 1,
#'                                      ess_sweeps = 2, verbose = FALSE))
#' fit
#' @export
fit_areal <- function(data, formula = ~ trend(), covariates = NULL, scale_covariates = TRUE, nugget = TRUE,
                      priors = areal_priors(), control = mcmc_control(), checkpoint_dir = NULL) {
  if (!inherits(data, "stareal_data")) stop("`data` must be made by areal_data() or simulate_areal()", call. = FALSE)
  priors <- complete_settings(priors, areal_priors, "priors")
  control <- complete_settings(control, mcmc_control, "control")
  if (!is.null(checkpoint_dir) && (!is.character(checkpoint_dir) || length(checkpoint_dir) != 1))
    stop("`checkpoint_dir` must be one folder path", call. = FALSE)
  design <- build_design(formula, covariates, data$frame, data$time, scale = scale_covariates)
  fit <- fit_with_design(data, design, nugget, priors, control, checkpoint_dir)
  fit$call <- match.call()
  fit
}

# fit_areal() given a built design and completed settings (crossval()
# refits with the fit's own)
fit_with_design <- function(data, design, nugget, priors, control, checkpoint_dir = NULL) {
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
  setup <- fit_setup(cells, frame$eig, T, priors, nugget, X)

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
                 design = design, data = data, priors = priors, control = control, nugget = isTRUE(nugget),
                 checkpoint_dir = checkpoint_dir, call = NULL),
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
  fit_setup(cells, fit$data$frame$eig, fit$data$time$n_steps, fit$priors, fit$nugget,
            design_eval(fit$design, cells$area, cells$t))
}

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
  tn <- theta_names(fit$nugget)
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
  nat <- theta_natural(matrix(d$theta, prod(dn[1:2]), dn[3]))   # (draw, chain) rows, draw fastest
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
#' The stored draws as a `posterior` draws object: the hyperparameters
#' (`sigma`, `rho_s`, `rho_t`, `sigma_n`), the fixed effects, and, if
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
  tn <- theta_names(x$nugget)
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
