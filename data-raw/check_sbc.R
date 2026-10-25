# check_sbc.R
#
# Simulation-based calibration (Talts et al. 2018, arXiv:1804.06788) of
# fit_areal(): draw the parameters from the default priors, simulate data
# from the model, fit it, and record the rank of each true value among 99
# thinned posterior draws. If the sampler and the model code are right,
# every parameter's ranks are uniform on 0..99, whatever the design. A
# check of the software, not of the model's fit to real data.
#
# A small lattice (6 x 6, 10 steps, 60 surveyed cells of 50) keeps each fit
# to seconds. Short chains are enough for the ranks: most fits do not pass
# diagnose()'s ESS thresholds, but 99 draws thinned from 500 are close to
# independent.
#
# Not part of the package. From the package root:
#   Rscript data-raw/check_sbc.R
# About 15 minutes for 200 replicates on 4 cores.
#
# First run, 2026-10-25: 200 replicates; chi-square p-values for uniform
# ranks 0.32-0.61, Kolmogorov-Smirnov 0.08-0.71, for all six parameters.
#
# ------------------------------------------------------------------

n_blocks <- 4                            # worker processes
n_rep <- 50                              # replicates per block

one_block <- function(block, n_rep, root) {
  Sys.setenv(VECLIB_MAXIMUM_THREADS = "1")
  suppressMessages(pkgload::load_all(root, quiet = TRUE))
  fr <- areal_lattice(6, 6)
  tf <- time_frame(n_steps = 10)
  ctl <- list(n_iter = 4000, n_burn = 1500, n_store = 250, n_chains = 2, n_cores = 1, theta_thin = 2,
              verbose = FALSE)
  res <- list()
  for (r in seq_len(n_rep)) {
    seed <- 1000 * block + r
    set.seed(seed)
    th <- c(sigma = abs(rnorm(1, 0, 2)), rho_s = rbeta(1, 2, 1), rho_t = rbeta(1, 1, 1),
            sigma_n = abs(rnorm(1, 0, 1)))                  # the defaults of areal_priors()
    b <- rnorm(2, 0, 3)
    sim <- simulate_areal(fr, tf, sigma = th[["sigma"]], rho_s = th[["rho_s"]], rho_t = th[["rho_t"]],
                          sigma_n = th[["sigma_n"]], beta = b, n_cells = 60, k = 50, seed = seed)
    fit <- tryCatch(suppressWarnings(fit_areal(sim, control = c(ctl, seed = seed))), error = function(e) NULL)
    if (is.null(fit)) next
    truth <- c(th, intercept = b[1], trend = b[2])
    dm <- posterior::as_draws_matrix(fit)
    dm <- dm[round(seq(1, nrow(dm), length.out = 99)), names(truth)]
    res[[length(res) + 1]] <- data.frame(block = block, rep = r, variable = names(truth), truth = truth,
                                         rank = colSums(sweep(unclass(dm), 2, truth, "<")))
  }
  do.call(rbind, res)
}

cl <- parallel::makePSOCKcluster(n_blocks)
r <- do.call(rbind, parallel::clusterApply(cl, seq_len(n_blocks), one_block, n_rep = n_rep, root = getwd()))
parallel::stopCluster(cl)

cat(sprintf("%d replicates\n", length(unique(paste(r$block, r$rep)))))
out <- do.call(rbind, lapply(split(r, r$variable), function(z) {
  bins <- table(cut(z$rank, breaks = seq(-0.5, 99.5, length.out = 11)))
  data.frame(variable = z$variable[1], n = nrow(z),
             chisq_p = round(suppressWarnings(stats::chisq.test(bins))$p.value, 3),
             ks_p = round(suppressWarnings(stats::ks.test((z$rank + stats::runif(nrow(z))) / 100, "punif"))$p.value, 3),
             bins = paste(as.integer(bins), collapse = " "))
}))
print(out, row.names = FALSE)
