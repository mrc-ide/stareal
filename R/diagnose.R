# diagnose.R
#
# Convergence and mixing for EVERY parameter of a fit (the hyperparameters,
# the fixed effects and each observed cell's latent value) and a pass/fail
# verdict. The lesson of 2026-10-11 in Africa_areal: a 40k run whose
# hyperparameters looked converged (R-hat <= 1.04) had latent cells at
# median ESS 14; a run is only converged when the latent line is clean too.
#
# R-hat and ESS are the rank-normalised versions from the `posterior`
# package, so they do not depend on the scale (logit or prevalence).
#
# Ported from Africa_areal/R/diagnostics.R (2026-10-12). Changes: reads a
# stareal_fit, not an archive on disk; the thresholds are arguments (the
# defaults are Africa_areal's settings.yml); the latent rows carry the
# cell's area id and time label; the verdict lists the checks that failed.
#
# ------------------------------------------------------------------

#' Check convergence and mixing of a fit
#'
#' Computes rank-normalised R-hat, bulk and tail effective sample sizes
#' (ESS) and Monte Carlo standard errors for every parameter: the
#' hyperparameters (from the full trace after burn-in), the fixed effects,
#' and the prevalence of every observed cell (from the stored draws). A fit
#' passes only if all of them are acceptable. Hyperparameters can look
#' converged while the latent cells mix poorly, so the latent cells are
#' always checked.
#'
#' @param fit A [fit_areal()] object.
#' @param rhat_max Largest acceptable R-hat, for every parameter.
#' @param ess_theta_min Smallest acceptable bulk ESS of each hyperparameter.
#' @param ess_latent_median_min Smallest acceptable median bulk ESS over the
#'   observed cells.
#' @return An object of class `stareal_diagnostics`:
#'   \describe{
#'     \item{table}{one row per parameter: `type` (`theta`, `beta` or
#'       `latent`), `variable`, posterior mean, sd and quantiles (latent
#'       cells on the prevalence scale), `rhat`, `ess_bulk`, `ess_tail`,
#'       `mcse_mean`, `n_draws`, and for latent cells `cell_id`, `id` and
#'       `time`}
#'     \item{checks}{a data frame of the checks: value, threshold, pass}
#'     \item{pass}{`TRUE` if every check passes}
#'   }
#' @examples
#' fr <- areal_lattice(4, 4)
#' tf <- time_frame(T = 5)
#' sim <- simulate_areal(fr, tf, n_cells = 30, seed = 1)
#' d <- areal_data(sim$data, fr, tf, y = "positive", k = "tested", area = "id", t = "t",
#'                 verbose = FALSE)
#' fit <- fit_areal(d, control = mcmc_control(n_iter = 100, n_burn = 50, n_store = 25,
#'                                            n_chains = 2, parallel = FALSE,
#'                                            ess_sweeps = 2, progress = FALSE))
#' diagnose(fit)   # far too short to pass
#' @export
diagnose <- function(fit, rhat_max = 1.05, ess_theta_min = 400, ess_latent_median_min = 100) {
  if (!inherits(fit, "stareal_fit")) stop("`fit` must be made by fit_areal()", call. = FALSE)
  sa <- stored_array(fit, latent = TRUE)
  cells <- fit$data$cells
  beta <- sa[, , fit$setup$beta_names, drop = FALSE]
  prev <- stats::plogis(sa[, , sprintf("logit[%s]", cells$cell_id), drop = FALSE])
  dimnames(prev)[[3]] <- sprintf("prev[%s]", cells$cell_id)
  tab <- rbind(draws_table(trace_array(fit), "theta"), draws_table(beta, "beta"), draws_table(prev, "latent"))
  is_lat <- tab$type == "latent"
  tab$cell_id <- NA_integer_; tab$id <- NA_character_; tab$time <- NA_character_
  tab$cell_id[is_lat] <- cells$cell_id; tab$id[is_lat] <- cells$id; tab$time[is_lat] <- cells$time

  th <- tab[tab$type == "theta", ]; be <- tab[tab$type == "beta", ]; la <- tab[is_lat, ]
  checks <- data.frame(
    check = c("theta R-hat (max)", "beta R-hat (max)", "latent R-hat (max)", "theta bulk ESS (min)",
              "latent bulk ESS (median)"),
    value = c(max(th$rhat), max(be$rhat), max(la$rhat), min(th$ess_bulk), stats::median(la$ess_bulk)),
    threshold = c(rhat_max, rhat_max, rhat_max, ess_theta_min, ess_latent_median_min),
    stringsAsFactors = FALSE)
  checks$pass <- ifelse(seq_len(nrow(checks)) <= 3, checks$value <= checks$threshold,
                        checks$value >= checks$threshold)
  checks$pass[is.na(checks$pass)] <- FALSE
  structure(list(table = tab, checks = checks, pass = all(checks$pass),
                 latent = list(n = nrow(la), rhat_over = sum(la$rhat > rhat_max, na.rm = TRUE),
                               ess_min = min(la$ess_bulk), mcse_mean_max = max(la$mcse_mean)),
                 acceptance = acceptance(fit),
                 seconds_per_iteration = mean(fit$seconds) / fit$control$n_iter),
            class = "stareal_diagnostics")
}

#' @export
print.stareal_diagnostics <- function(x, ...) {
  cat("<stareal_diagnostics>\n")
  tb <- x$table[x$table$type != "latent", c("variable", "mean", "q2.5", "q97.5", "rhat", "ess_bulk", "ess_tail")]
  tb[c("mean", "q2.5", "q97.5")] <- lapply(tb[c("mean", "q2.5", "q97.5")], signif, 3)
  tb$rhat <- round(tb$rhat, 3); tb$ess_bulk <- round(tb$ess_bulk); tb$ess_tail <- round(tb$ess_tail)
  print(tb, row.names = FALSE)
  cat(sprintf("latent cells: %d; R-hat over threshold in %d; bulk ESS median %.0f, min %.0f; max MCSE of the mean %.2g (prevalence)\n",
              x$latent$n, x$latent$rhat_over, x$checks$value[5], x$latent$ess_min, x$latent$mcse_mean_max))
  cat(sprintf("acceptance after burn-in: centred %.2f, surrogate %.2f\n\n",
              x$acceptance[["centred"]], x$acceptance[["surrogate"]]))
  for (i in seq_len(nrow(x$checks)))
    cat(sprintf("  %-26s %8.3g  %s %-6g %s\n", x$checks$check[i], x$checks$value[i],
                if (i <= 3) "<=" else ">=", x$checks$threshold[i], if (x$checks$pass[i]) "ok" else "FAIL"))
  cat(sprintf("\n%s\n", if (x$pass) "PASS" else "FAIL: run longer chains (see ?mcmc_control)"))
  invisible(x)
}
