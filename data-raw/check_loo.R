# check_loo.R
#
# loo() (integrated PSIS-LOO, from one fit) against exact leave-one-out:
# one refit per observed cell, each held-out cell scored by crossval()'s
# exact quadrature log score. On simulated data (6 x 6 lattice, 8 steps,
# 40 surveyed cells), for several priors, sample sizes and with and without
# the nugget. A check of the approximation, not of the model.
#
# Not part of the package. From the package root:
#   Rscript data-raw/check_loo.R
# About 25 minutes on 4 cores (the exact refits; loo() takes seconds).
#
# First run, 2026-10-28 (4 chains x 3,000 iterations):
#   case            elpd exact  elpd loo   diff   max |cell diff|  cor     k max
#   leroux, k 50       -99.35     -98.99   0.35       0.31        0.9995   0.68
#   bym2, k 50        -137.68    -137.75  -0.07       0.11        0.9994   0.65
#   leroux, k 300     -165.98    -165.53   0.44       0.32        0.9993   0.65
#   car, no nugget    -119.88    -119.69   0.18       0.07        0.9999   0.78
# The exact values carry their own Monte Carlo error (one short refit per
# cell).
#
# ------------------------------------------------------------------

Sys.setenv(VECLIB_MAXIMUM_THREADS = "1", MKL_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1")
pkgload::load_all(".", quiet = TRUE)
fr <- areal_lattice(6, 6)
tf <- time_frame(n_steps = 8)
ctl <- list(n_iter = 3000, n_burn = 1000, n_store = 250, n_chains = 4, n_cores = 4, verbose = FALSE)
cases <- list(
  leroux_k50 = list(spatial = "leroux", k = 50, nugget = TRUE),
  bym2_k50 = list(spatial = "bym2", k = 50, nugget = TRUE),
  leroux_k300 = list(spatial = "leroux", k = 300, nugget = TRUE),
  car_no_nugget = list(spatial = "car", k = 50, nugget = FALSE))
out <- list()
for (nm in names(cases)) {
  cs <- cases[[nm]]
  sim <- simulate_areal(fr, tf, n_cells = 40, k = cs$k, seed = 11, spatial = cs$spatial,
                        sigma_n = if (cs$nugget) 0.4 else 0)
  fit <- fit_areal(sim, spatial = cs$spatial, nugget = cs$nugget, control = ctl)
  l <- loo(fit)
  ex <- crossval(fit, folds = seq_len(nrow(fit$data$cells)))       # each cell its own fold
  e <- ex$cells$log_score[match(fit$data$cells$cell_id, ex$cells$cell_id)]
  a <- l$cells$elpd_loo
  out[[nm]] <- data.frame(case = nm, pass = diagnose(fit)$pass, elpd_exact = sum(e), elpd_loo = sum(a),
                          diff = sum(a) - sum(e), max_abs_cell = max(abs(a - e)), cor = cor(a, e),
                          k_max = max(l$cells$pareto_k))
}
print(do.call(rbind, out), row.names = FALSE)
