# fit.R
#
# fit_areal(): the user-facing fit, a thin layer over fit_setup() and
# run_chains() (chains.R); and the fit object's print, summary, plot and
# posterior::as_draws_array methods.
#
# ------------------------------------------------------------------

#' Fit the space-time areal model
#'
#' Fits the model to the observed area-time cells of `data` by MCMC. The
#' logit prevalence of cell (area \eqn{a}, time \eqn{t}) is
#' \deqn{x_{at} = \beta_0 + \beta_1 s_t + f_{at} + \epsilon_{at},}
#' where \eqn{s_t} runs from -1 to 1 over the time window, \eqn{f} is a
#' Gaussian field with separable covariance (a Leroux prior on the frame's
#' graph in space, with scale `sigma` and dependence `rho_s`; a
#' unit-variance AR1 with correlation `rho_t` in time), and
#' \eqn{\epsilon} is a cell-level nugget with standard deviation `sigma_n`.
#' Positives are binomial given \eqn{x}.
#'
#' The unobserved cells and the fixed effects are integrated out exactly, so
#' the sampler works on the observed cells only; its cost per iteration
#' grows as the cube of the number of observed cells. The full area-time
#' surface is reconstructed afterwards by `predict()`.
#'
#' @param data A [areal_data()] object.
#' @param nugget Include the cell-level nugget.
#' @param priors An [areal_priors()] object.
#' @param control An [mcmc_control()] object.
#' @return An object of class `stareal_fit`. Its elements include
#'   \describe{
#'     \item{trace}{the hyperparameter trace, burn-in included: chain,
#'       iteration, the hyperparameters on the unconstrained scale
#'       (`log_sigma`, `logit_rho_s`, `logit_rho_t`, `log_sigma_n`), the log
#'       posterior `lp`, the log-likelihood `ll`, and the acceptance rates
#'       of the two hyperparameter moves}
#'     \item{draws}{the stored posterior draws, as arrays (draw x chain x
#'       variable): `theta` (unconstrained), `beta`, and `x` (logit
#'       prevalence of each observed cell, nugget included, in the order of
#'       `data$cells`)}
#'     \item{data, priors, control, nugget}{as given}
#'   }
#'   Use [summary()][summary.stareal_fit], [diagnose()],
#'   [plot()][plot.stareal_fit], and `posterior::as_draws_array(fit)`.
#' @examples
#' fr <- areal_lattice(4, 4)
#' tf <- time_frame(T = 5)
#' sim <- simulate_areal(fr, tf, n_cells = 30, seed = 1)
#' d <- areal_data(sim$data, fr, tf, y = "positive", k = "tested", area = "id", t = "t",
#'                 verbose = FALSE)
#' fit <- fit_areal(d, control = mcmc_control(n_iter = 100, n_burn = 50, n_store = 25,
#'                                            n_chains = 1, ess_sweeps = 2, progress = FALSE))
#' fit
#' @export
fit_areal <- function(data, nugget = TRUE, priors = areal_priors(), control = mcmc_control()) {
  if (!inherits(data, "stareal_data")) stop("`data` must be made by areal_data()", call. = FALSE)
  if (!inherits(priors, "stareal_priors")) stop("`priors` must be made by areal_priors()", call. = FALSE)
  if (!inherits(control, "stareal_control")) stop("`control` must be made by mcmc_control()", call. = FALSE)
  cells <- data$cells; frame <- data$frame; T <- data$time$T
  n <- nrow(cells)
  if (n < 2) stop("at least two observed cells are needed", call. = FALSE)
  setup <- fit_setup(cells, frame$eig, T, priors, nugget)
  p <- ncol(setup$obs$X)
  if (!length(priors$beta_mean) %in% c(1, p) || !length(priors$beta_sd) %in% c(1, p))
    stop(sprintf("`beta_mean` and `beta_sd` in areal_priors() need 1 or %d values (%s)", p,
                 paste(setup$beta_names, collapse = ", ")), call. = FALSE)

  # scale guard-rails (plan section 6): the eigendecomposition is m x m, and
  # each iteration costs O(n^3)
  if (frame$m > 8000) warning(sprintf("m = %d units: the dense eigendecomposition and prediction will be slow and large",
                                      frame$m), call. = FALSE)
  if (n > 4000) warning(sprintf(paste("n = %d observed cells: each iteration costs O(n^3) and will be slow;",
                                      "use an optimised BLAS (MKL, OpenBLAS or Accelerate)"), n), call. = FALSE)
  if (control$progress)
    message(sprintf("fitting %d observed cells (%d units x %d time steps): %d chain%s x %d iterations%s",
                    n, frame$m, T, control$n_chains, if (control$n_chains > 1) "s" else "", control$n_iter,
                    if (control$parallel && control$n_chains > 1 && control$n_cores > 1)
                      sprintf(", in parallel on %d workers", min(control$n_chains, control$n_cores)) else ""))

  res <- with_seed(control$seed, run_chains(setup, control))

  to_array <- function(l) {
    a <- array(NA_real_, c(nrow(l[[1]]), length(l), ncol(l[[1]])),
               dimnames = list(NULL, NULL, colnames(l[[1]])))
    for (i in seq_along(l)) a[, i, ] <- l[[i]]
    a
  }
  xs <- lapply(res, `[[`, "x")
  xs <- lapply(xs, function(x) { colnames(x) <- cells$cell_id; x })
  trace <- do.call(rbind, lapply(res, function(r) cbind(chain = r$chain, r$theta)))
  rownames(trace) <- NULL
  setup$eig <- NULL                                   # it is on the frame
  structure(list(trace = trace,
                 draws = list(theta = to_array(lapply(res, `[[`, "draw_theta")),
                              beta = to_array(lapply(res, `[[`, "beta")), x = to_array(xs)),
                 stored_iter = res[[1]]$stored_iter,
                 seconds = vapply(res, `[[`, numeric(1), "seconds"),
                 adaptation = lapply(res, function(r) list(centred = r$am_centred, surrogate = r$am_surrogate)),
                 setup = setup, data = data, priors = priors, control = control, nugget = isTRUE(nugget)),
            class = "stareal_fit")
}

# ---- draws ------------------------------------------------------------------

# the hyperparameter trace after burn-in, natural scale, as an
# iterations x chains x variables array (the full trace, not only the stored
# draws: more draws for theta's R-hat and ESS)
trace_array <- function(fit) {
  tn <- theta_names(fit$nugget)
  post <- fit$trace[fit$trace$iteration > fit$control$n_burn, ]
  ch <- sort(unique(post$chain))
  l <- lapply(ch, function(c) theta_natural(as.matrix(post[post$chain == c, tn])))
  a <- array(NA_real_, c(nrow(l[[1]]), length(l), ncol(l[[1]])), dimnames = list(NULL, NULL, colnames(l[[1]])))
  for (i in seq_along(l)) a[, i, ] <- l[[i]]
  a
}

# the stored draws: hyperparameters (natural scale) and fixed effects, and
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
#' The stored draws as a `posterior::draws_array`: the hyperparameters
#' (`sigma`, `rho_s`, `rho_t`, `sigma_n`), the fixed effects, and, if
#' `latent`, the logit prevalence of each observed cell, nugget included,
#' named `logit[<cell_id>]`.
#'
#' @param x A `stareal_fit`.
#' @param latent Include the observed cells' latent values.
#' @param ... Not used.
#' @return A `posterior::draws_array`.
#' @exportS3Method posterior::as_draws_array
as_draws_array.stareal_fit <- function(x, latent = TRUE, ...) {
  a <- stored_array(x, latent)
  dimnames(a) <- list(iteration = NULL, chain = NULL, variable = dimnames(a)[[3]])
  posterior::as_draws_array(a)
}

# ---- summaries --------------------------------------------------------------

# one row per variable of an iterations x chains x variables array:
# posterior summaries plus rank-normalised R-hat, bulk/tail ESS and MCSE
# (from Africa_areal/R/diagnostics.R, diagnostics_table)
draws_table <- function(arr, type) {
  s <- posterior::summarise_draws(posterior::as_draws_array(arr), "mean", "sd",
                                  ~quantile(.x, c(0.025, 0.5, 0.975)), "rhat", "ess_bulk",
                                  "ess_tail", "mcse_mean")
  s <- as.data.frame(s)
  names(s) <- c("variable", "mean", "sd", "q2.5", "q50", "q97.5", "rhat", "ess_bulk", "ess_tail", "mcse_mean")
  cbind(type = type, s, n_draws = prod(dim(arr)[1:2]))
}

# acceptance of the two theta moves after burn-in
acceptance <- function(fit) {
  post <- fit$trace[fit$trace$iteration > fit$control$n_burn, ]
  c(centred = mean(post$accept_centred), surrogate = mean(post$accept_surrogate))
}

#' @export
print.stareal_fit <- function(x, ...) {
  ct <- x$control
  cat(sprintf("<stareal_fit> %d observed cells on %d units x %d time steps%s\n", nrow(x$data$cells),
              x$data$frame$m, x$data$time$T, if (x$nugget) ", with nugget" else ""))
  cat(sprintf("  %d chain%s x %d iterations (%d burn-in); %d stored draws in all; %s per chain\n",
              ct$n_chains, if (ct$n_chains > 1) "s" else "", ct$n_iter, ct$n_burn, ct$n_store * ct$n_chains,
              format_seconds(mean(x$seconds))))
  ta <- trace_array(x)
  med <- apply(ta, 3, stats::median)
  cat("  posterior medians: ", paste(sprintf("%s %.3g", names(med), med), collapse = ", "), "\n", sep = "")
  cat("  see summary(), diagnose() and plot()\n")
  invisible(x)
}

#' Summarise a fit
#'
#' Posterior summaries of the hyperparameters and fixed effects, with
#' R-hat and effective sample sizes. For the full check, including every
#' observed cell, use [diagnose()].
#'
#' @param object A `stareal_fit`.
#' @param ... Not used.
#' @return An object of class `summary.stareal_fit`: `table` (one row per
#'   parameter), `acceptance` (of the two hyperparameter moves, after
#'   burn-in) and `seconds_per_iteration`.
#' @export
summary.stareal_fit <- function(object, ...) {
  sa <- stored_array(object, latent = FALSE)
  beta <- sa[, , object$setup$beta_names, drop = FALSE]
  tab <- rbind(draws_table(trace_array(object), "theta"), draws_table(beta, "beta"))
  structure(list(table = tab, acceptance = acceptance(object),
                 seconds_per_iteration = mean(object$seconds) / object$control$n_iter,
                 n_cells = nrow(object$data$cells), control = object$control),
            class = "summary.stareal_fit")
}

#' @export
print.summary.stareal_fit <- function(x, digits = 3, ...) {
  cat(sprintf("<stareal_fit summary> %d observed cells; %d chain%s, %d iterations after burn-in\n\n",
              x$n_cells, x$control$n_chains, if (x$control$n_chains > 1) "s" else "",
              x$control$n_iter - x$control$n_burn))
  tab <- x$table[, c("variable", "mean", "sd", "q2.5", "q50", "q97.5", "rhat", "ess_bulk", "ess_tail")]
  tab$ess_bulk <- round(tab$ess_bulk); tab$ess_tail <- round(tab$ess_tail)
  num <- c("mean", "sd", "q2.5", "q50", "q97.5")
  tab[num] <- lapply(tab[num], signif, digits)
  tab$rhat <- round(tab$rhat, 3)
  print(tab, row.names = FALSE)
  cat(sprintf("\nacceptance after burn-in: centred %.2f, surrogate %.2f; %.3g s per iteration\n",
              x$acceptance[["centred"]], x$acceptance[["surrogate"]], x$seconds_per_iteration))
  invisible(x)
}

#' Trace plots of a fit
#'
#' The hyperparameter traces on the natural scale, one colour per chain,
#' burn-in included (left of the dashed line).
#'
#' @param x A `stareal_fit`.
#' @param ... Passed to [graphics::matplot()].
#' @return `x`, invisibly.
#' @export
plot.stareal_fit <- function(x, ...) {
  tn <- theta_names(x$nugget)
  ch <- sort(unique(x$trace$chain))
  nat <- theta_natural(as.matrix(x$trace[tn]))
  it <- x$trace$iteration[x$trace$chain == ch[1]]
  op <- graphics::par(mfrow = c(ceiling(ncol(nat) / 2), 2), mar = c(4, 4, 1.5, 1))
  on.exit(graphics::par(op))
  cols <- grDevices::hcl.colors(max(2, length(ch)), "Dark 3")[seq_along(ch)]
  for (v in colnames(nat)) {
    m <- vapply(ch, function(c) nat[x$trace$chain == c, v], numeric(length(it)))
    graphics::matplot(it, m, type = "l", lty = 1, col = cols, xlab = "iteration", ylab = v, main = v, ...)
    graphics::abline(v = x$control$n_burn, lty = 2, col = "grey40")
  }
  invisible(x)
}
