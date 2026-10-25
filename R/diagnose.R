# diagnose.R
#
# summary() and diagnose() for a fit: one table, one class.
# summary() covers the hyperparameters and fixed effects; diagnose() adds
# every observed cell's latent value and a pass/fail verdict. The lesson of
# 2026-10-11 in Africa_areal: a 40k run whose hyperparameters looked
# converged (R-hat <= 1.04) had latent cells at median ESS 14; a run is only
# converged when the latent line is clean too.
#
# R-hat, ESS and MCSE are the rank-normalised versions from the `posterior`
# package. ESS does not depend on the scale; R-hat does slightly, so the
# latent cells are checked on the scale they are reported on (prevalence).
#
# Ported from Africa_areal/R/diagnostics.R (2026-10-12). Changes: reads a
# stareal_fit; thresholds are arguments (the defaults are Africa_areal's
# settings.yml); the latent rows carry the cell's area id and time label;
# the fixed effects' ESS is checked too; summary() and diagnose() share one
# class.
#
# ------------------------------------------------------------------

#' Summarise a fit
#'
#' Posterior summaries of the hyperparameters and fixed effects, with R-hat
#' and effective sample sizes (ESS). The hyperparameters are summarised from
#' the full trace after burn-in; the fixed effects (and, with `latent =
#' TRUE`, the observed cells) from the stored draws. For the full check of
#' convergence, use [diagnose()].
#'
#' @param object A `stareal_fit`.
#' @param latent Also summarise each observed cell's latent prevalence
#'   (nugget included).
#' @param ... Not used.
#' @return An object of class `summary.stareal_fit`: `table` (one row per
#'   parameter: `type` (`theta`, `beta` or `latent`), `variable`, `mean`,
#'   `sd`, `q2.5`, `q50`, `q97.5`, `rhat`, `ess_bulk`, `ess_tail`,
#'   `mcse_mean`, `n_draws`, and for latent rows `cell_id`, `id` and
#'   `time`), `acceptance` (of the two hyperparameter moves, after burn-in),
#'   `seconds_per_iteration`, and, from [diagnose()], `checks` and `pass`.
#' @export
summary.stareal_fit <- function(object, latent = FALSE, ...) {
  check_dots(...)
  fit <- object
  tab <- rbind(draws_table(trace_array(fit), "theta"), draws_table(fit$draws$beta, "beta"))
  tab$cell_id <- NA_real_; tab$id <- NA_character_; tab$time <- NA_character_
  if (latent) {
    cells <- fit$data$cells
    prev <- stats::plogis(fit$draws$x)
    dimnames(prev)[[3]] <- sprintf("prev[%s]", cells$cell_id)
    lt <- draws_table(prev, "latent")
    lt$cell_id <- cells$cell_id; lt$id <- cells$id; lt$time <- cells$time
    tab <- rbind(tab, lt)
  }
  rownames(tab) <- NULL
  structure(list(table = tab, acceptance = acceptance(fit),
                 seconds_per_iteration = mean(fit$seconds) / fit$control$n_iter,
                 n_cells = nrow(fit$data$cells), control = fit$control, design = fit$design,
                 checks = NULL, pass = NULL),
            class = "summary.stareal_fit")
}

#' Check convergence and mixing of a fit
#'
#' [summary()][summary.stareal_fit] with every observed cell included and a
#' verdict. Rank-normalised R-hat, bulk and tail effective sample sizes
#' (ESS) and Monte Carlo standard errors are computed for every parameter:
#' the hyperparameters (from the full trace after burn-in), the fixed
#' effects, and the latent prevalence of every observed cell, nugget
#' included (from the stored draws; this is not the nugget-free surface of
#' [predict()][predict.stareal_fit]). A fit passes only if all are
#' acceptable. Hyperparameters can look converged while the latent cells
#' mix poorly, so the latent cells are always checked.
#'
#' @param fit A [fit_areal()] object.
#' @param rhat_max Largest acceptable R-hat, for every parameter.
#' @param ess_theta_min Smallest acceptable bulk ESS of each hyperparameter
#'   and each fixed effect.
#' @param ess_latent_median_min Smallest acceptable median bulk ESS over the
#'   observed cells.
#' @return A `summary.stareal_fit` (see [summary()][summary.stareal_fit])
#'   with the latent rows, and `checks` (a data frame: check, value,
#'   threshold, pass), `pass` (`TRUE` if every check passes) and `latent`
#'   (a few statistics over the cells).
#' @examples
#' fr <- areal_lattice(4, 4)
#' tf <- time_frame(n_steps = 5)
#' sim <- simulate_areal(fr, tf, n_cells = 30, seed = 1)
#' fit <- fit_areal(sim, control = list(n_iter = 100, n_burn = 50, n_store = 25, n_chains = 2,
#'                                      n_cores = 1, ess_sweeps = 2, verbose = FALSE))
#' suppressWarnings(diagnose(fit))   # far too short to pass
#' @export
diagnose <- function(fit, rhat_max = 1.05, ess_theta_min = 400, ess_latent_median_min = 100) {
  if (!inherits(fit, "stareal_fit")) stop("`fit` must be made by fit_areal()", call. = FALSE)
  s <- summary(fit, latent = TRUE)
  tab <- s$table
  th <- tab[tab$type == "theta", ]; be <- tab[tab$type == "beta", ]; la <- tab[tab$type == "latent", ]
  checks <- data.frame(
    check = c("theta R-hat (max)", "beta R-hat (max)", "latent R-hat (max)", "theta bulk ESS (min)",
              "beta bulk ESS (min)", "latent bulk ESS (median)"),
    value = c(max(th$rhat), max(be$rhat), max(la$rhat), min(th$ess_bulk), min(be$ess_bulk),
              stats::median(la$ess_bulk)),
    threshold = c(rhat_max, rhat_max, rhat_max, ess_theta_min, ess_theta_min, ess_latent_median_min),
    direction = c("<=", "<=", "<=", ">=", ">=", ">="),
    stringsAsFactors = FALSE)
  checks$pass <- ifelse(checks$direction == "<=", checks$value <= checks$threshold, checks$value >= checks$threshold)
  checks$pass[is.na(checks$pass)] <- FALSE
  s$checks <- checks
  s$pass <- all(checks$pass)
  s$latent <- list(n = nrow(la), rhat_over = sum(la$rhat > rhat_max, na.rm = TRUE),
                   ess_min = min(la$ess_bulk), mcse_mean_max = max(la$mcse_mean))
  s
}

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
print.summary.stareal_fit <- function(x, digits = 3, ...) {
  ct <- x$control
  cat(sprintf("<stareal_fit summary> %d observed cells; %d chain%s, %d iterations after burn-in\n",
              x$n_cells, ct$n_chains, if (ct$n_chains > 1) "s" else "", ct$n_iter - ct$n_burn))
  cat(sprintf("  hyperparameters from the trace (%d draws), fixed effects%s from the %d stored draws\n\n",
              x$table$n_draws[1], if (!is.null(x$checks)) " and cells" else "", ct$n_chains * ct$n_store))
  tab <- x$table[x$table$type != "latent", c("variable", "mean", "sd", "q2.5", "q50", "q97.5", "rhat",
                                             "ess_bulk", "ess_tail")]
  tab$ess_bulk <- round(tab$ess_bulk); tab$ess_tail <- round(tab$ess_tail)
  num <- c("mean", "sd", "q2.5", "q50", "q97.5")
  tab[num] <- lapply(tab[num], signif, digits)
  tab$rhat <- round(tab$rhat, 3)
  print(tab, row.names = FALSE)
  cat(sprintf("\nacceptance after burn-in: centred %.2f, surrogate %.2f; %.3g s per iteration\n",
              x$acceptance[["centred"]], x$acceptance[["surrogate"]], x$seconds_per_iteration))
  sc <- Filter(function(cl) cl$scale != 1 || cl$centre != 0, x$design$cols)
  if (length(sc)) cat(sprintf("standardised covariates (effects per sd): %s\n",
                              paste(sprintf("%s (centre %.3g, sd %.3g)", vapply(sc, `[[`, "", "name"),
                                            vapply(sc, `[[`, 0, "centre"), vapply(sc, `[[`, 0, "scale")),
                                    collapse = ", ")))
  if (!is.null(x$checks)) {
    cat(sprintf("latent cells: %d; R-hat over threshold in %d; bulk ESS median %.0f, min %.0f; max MCSE of the mean %.2g (prevalence)\n\n",
                x$latent$n, x$latent$rhat_over, x$checks$value[6], x$latent$ess_min, x$latent$mcse_mean_max))
    for (i in seq_len(nrow(x$checks)))
      cat(sprintf("  %-26s %8.3g  %s %-6g %s\n", x$checks$check[i], x$checks$value[i], x$checks$direction[i],
                  x$checks$threshold[i], if (x$checks$pass[i]) "ok" else "FAIL"))
    cat(sprintf("\n%s\n", if (x$pass) "PASS" else "FAIL: run longer chains (see ?mcmc_control)"))
  } else {
    cat("for every observed cell and a verdict, use diagnose()\n")
  }
  invisible(x)
}
