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
#   Rscript data-raw/check_sbc.R [model] [n_blocks] [n_rep] [first_block]
# first_block (default 1) numbers the blocks, and so sets the seeds: a
# second batch with first_block = n_blocks + 1 uses new seeds.
# model: leroux (the default), bym2, car, icar (simulated with the Leroux
# prior at rho_s = 1, fitted with spatial = "icar", so rho_s is not
# sampled) or fixed (Leroux, simulated and fitted with rho_t fixed at 0.5);
# or an additive model (main effects and interaction types, 2026-11-14):
#   time_rw2       Leroux, main = "time", temporal = "rw2" (interaction IV)
#   both_rw1_III   Leroux, main = "both", temporal = "rw1", interaction = "III"
#   space_bym2_II  BYM2, main = "space", interaction = "II"
# with the main effects' hyperparameters also drawn from the default priors.
# n_blocks worker processes run n_rep replicates each (default 4 x 50).
# About 15 minutes for 200 replicates on 4 cores.
#
# First run, 2026-10-25: 200 replicates; chi-square p-values for uniform
# ranks 0.32-0.61, Kolmogorov-Smirnov 0.08-0.71, for all six parameters.
# The new priors and fixed hyperparameters, 2026-10-28 (5 x 40 replicates
# each, on Windows with MKL): uniform for every parameter of bym2, car and
# fixed (chi-square p 0.12-0.94, KS 0.20-0.96) and icar (all p >= 0.23 at
# 10 x 40; at 200 replicates, rho_t had chi-square p 0.017, KS 0.44).
#
# Additive models, 2026-11-14 (Windows, MKL; 6 x 34 replicates each):
#   time_rw2 and space_bym2_II: uniform for every parameter (chi-square p
#   0.11-0.99, KS 0.28-0.95). both_rw1_III: uniform except sigma_space
#   (chi-square p 0.045, KS 0.49). A second batch of new seeds (first_block
#   7) flagged sigma (chi-square 0.005, KS 0.09), sigma_n (0.038, 0.031)
#   and sigma_space (0.041, 0.61). A third, with chains three times as long
#   (SBC_LONG=3, 16 x 13, first_block 20): those three uniform (0.77, 0.995,
#   0.58) and rho_s flagged instead (chi-square 0.017, KS 0.15), with only
#   57 of 208 fits passing diagnose(). The flag moves from parameter to
#   parameter between batches, as poor mixing does and a bias would not: on
#   60 cells this model (seven hyperparameters; a type III interaction
#   turns into the nugget as rho_s -> 0) is weakly identified. ?fit_areal
#   says so.
#
# ------------------------------------------------------------------

args <- commandArgs(TRUE)
model <- if (length(args) >= 1) args[1] else "leroux"
additive <- list(time_rw2 = list(spatial = "leroux", main = "time", temporal = "rw2", interaction = "IV"),
                 both_rw1_III = list(spatial = "leroux", main = "both", temporal = "rw1", interaction = "III"),
                 space_bym2_II = list(spatial = "bym2", main = "space", temporal = "ar1", interaction = "II"))
stopifnot(model %in% c("leroux", "bym2", "car", "icar", "fixed", names(additive)))
n_blocks <- if (length(args) >= 2) as.integer(args[2]) else 4L     # worker processes
n_rep <- if (length(args) >= 3) as.integer(args[3]) else 50L       # replicates per block
first <- if (length(args) >= 4) as.integer(args[4]) else 1L       # number of the first block (seeds)
long <- as.numeric(Sys.getenv("SBC_LONG", "1"))   # chain-length multiplier (SBC_LONG=3: three times as long)
Sys.setenv(VECLIB_MAXIMUM_THREADS = "1", MKL_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1")   # one BLAS thread per worker

one_block <- function(block, n_rep, root, model, additive, long) {
  suppressMessages(pkgload::load_all(root, quiet = TRUE))
  fr <- areal_lattice(6, 6)
  tf <- time_frame(n_steps = 10)
  ctl <- list(n_iter = 4000 * long, n_burn = 1500 * long, n_store = 250, n_chains = 2, n_cores = 1, theta_thin = 2,
              verbose = FALSE)
  res <- list()
  for (r in seq_len(n_rep)) {
    seed <- 1000 * block + r
    set.seed(seed)
    th <- c(sigma = abs(rnorm(1, 0, 2)), rho_s = rbeta(1, 2, 1), rho_t = rbeta(1, 1, 1),
            sigma_n = abs(rnorm(1, 0, 1)))                  # the defaults of areal_priors()
    b <- rnorm(2, 0, 3)
    th_main <- c(sigma_space = abs(rnorm(1, 0, 1)), rho_space = rbeta(1, 2, 1), sigma_time = abs(rnorm(1, 0, 1)),
                 rho_time = rbeta(1, 1, 1))              # (drawn after the others: the four above are unchanged)
    if (model %in% names(additive)) {
      a <- additive[[model]]
      sim <- do.call(simulate_areal, c(list(fr, tf, sigma = th[["sigma"]], rho_s = th[["rho_s"]], rho_t = th[["rho_t"]],
                                            sigma_n = th[["sigma_n"]], beta = b, n_cells = 60, k = 50, seed = seed),
                                       as.list(th_main), a))
      fit <- tryCatch(suppressWarnings(do.call(fit_areal, c(list(sim, control = c(ctl, seed = seed)), a))),
                      error = function(e) NULL)
      if (is.null(fit)) next
      truth <- c(sim$truth$theta, intercept = b[1], trend = b[2])
      dm <- posterior::as_draws_matrix(fit)
      dm <- dm[round(seq(1, nrow(dm), length.out = 99)), names(truth)]
      res[[length(res) + 1]] <- data.frame(block = block, rep = r, variable = names(truth), truth = truth,
                                           rank = colSums(sweep(unclass(dm), 2, truth, "<")),
                                           pass = diagnose(fit)$pass)
      next
    }
    if (model == "icar") th[["rho_s"]] <- 1
    if (model == "fixed") th[["rho_t"]] <- 0.5
    sim <- simulate_areal(fr, tf, sigma = th[["sigma"]], rho_s = th[["rho_s"]], rho_t = th[["rho_t"]],
                          sigma_n = th[["sigma_n"]], beta = b, n_cells = 60, k = 50, seed = seed,
                          spatial = if (model %in% c("bym2", "car")) model else "leroux")
    args <- list(sim, control = c(ctl, seed = seed))
    if (model %in% c("bym2", "car", "icar")) args$spatial <- model
    if (model == "fixed") args$fixed <- list(rho_t = 0.5)
    fit <- tryCatch(suppressWarnings(do.call(fit_areal, args)), error = function(e) NULL)
    if (is.null(fit)) next
    truth <- c(th, intercept = b[1], trend = b[2])
    truth <- truth[setdiff(names(truth), switch(model, icar = "rho_s", fixed = "rho_t", character(0)))]
    dm <- posterior::as_draws_matrix(fit)
    dm <- dm[round(seq(1, nrow(dm), length.out = 99)), names(truth)]
    res[[length(res) + 1]] <- data.frame(block = block, rep = r, variable = names(truth), truth = truth,
                                         rank = colSums(sweep(unclass(dm), 2, truth, "<")))
  }
  do.call(rbind, res)
}

cl <- parallel::makePSOCKcluster(n_blocks)
r <- do.call(rbind, parallel::clusterApply(cl, first - 1L + seq_len(n_blocks), one_block, n_rep = n_rep, root = getwd(),
                                           model = model, additive = additive, long = long))
parallel::stopCluster(cl)

cat(sprintf("model %s: %d replicates of %d\n", model, length(unique(paste(r$block, r$rep))), n_blocks * n_rep))
out <- do.call(rbind, lapply(split(r, r$variable), function(z) {
  bins <- table(cut(z$rank, breaks = seq(-0.5, 99.5, length.out = 11)))
  data.frame(variable = z$variable[1], n = nrow(z),
             chisq_p = round(suppressWarnings(stats::chisq.test(bins))$p.value, 3),
             ks_p = round(suppressWarnings(stats::ks.test((z$rank + stats::runif(nrow(z))) / 100, "punif"))$p.value, 3),
             bins = paste(as.integer(bins), collapse = " "))
}))
print(out, row.names = FALSE)
if ("pass" %in% names(r))
  cat(sprintf("fits passing diagnose(): %d of %d\n", sum(r$pass[!duplicated(paste(r$block, r$rep))]),
              length(unique(paste(r$block, r$rep)))))
