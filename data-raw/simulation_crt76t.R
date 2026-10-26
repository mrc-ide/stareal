# simulation_crt76t.R
#
# Makes inst/extdata/simulation_crt76t.rds: the repeated-simulation study
# of the tutorial "Simulate, fit, check", which imports it. Twenty
# simulations, each fitted and predicted, are too long to run inside a
# tutorial, so they are run once here.
#
# Each simulation copies the survey design of the crt K76T data of tutorial
# 1 (setup_crt76t.R: the same area-months and numbers tested), with the
# parameters the tutorial uses, close to that fit's posterior. Each is fitted
# with the default sampler settings, the length of the tutorial's fit. Kept:
# each parameter's posterior mean and 95% interval, and the coverage and
# width of the surface's 95% intervals, for surveyed and unsurveyed
# area-months.
#
# Rerun when the sampler, the model or the data change. From the package
# root:
#   Rscript data-raw/simulation_crt76t.R
#
# ------------------------------------------------------------------

source("data-raw/setup_crt76t.R")

truth <- list(sigma = 3, rho_s = 0.97, rho_t = 0.997, sigma_n = 0.6, beta = c(intercept = 0.5, trend = -2))
study <- lapply(1:20, function(r) {
  s <- simulate_areal(fr, tf, cells = d$cells, sigma = truth$sigma, rho_s = truth$rho_s, rho_t = truth$rho_t,
                      sigma_n = truth$sigma_n, beta = truth$beta, seed = 1000 + r)
  f <- fit_areal(s, control = list(verbose = FALSE))
  tab <- summary(f)$table
  pr <- predict(f, verbose = FALSE)
  tv <- s$truth$prevalence[cbind(match(pr$area$id, fr$ids), pr$area$t)]
  message(sprintf("simulation %d done", r))
  list(params = data.frame(rep = r, variable = tab$variable, mean = tab$mean, q2.5 = tab$q2.5, q97.5 = tab$q97.5,
                           truth = c(unlist(truth[c("sigma", "rho_s", "rho_t", "sigma_n")]), truth$beta)[tab$variable],
                           pass = diagnose(f)$pass),
       surface = aggregate(cbind(covered = pr$area$q2.5 <= tv & tv <= pr$area$q97.5,
                                 width = pr$area$q97.5 - pr$area$q2.5) ~ observed,
                           data.frame(observed = pr$area$k > 0), mean))
})

dir.create("inst/extdata", showWarnings = FALSE, recursive = TRUE)
saveRDS(list(truth = truth,
             params = do.call(rbind, lapply(study, `[[`, "params")),
             surface = do.call(rbind, lapply(seq_along(study), function(r) cbind(rep = r, study[[r]]$surface))),
             built = list(date = Sys.Date(), stareal = as.character(utils::packageVersion("stareal")),
                          sampler = SAMPLER_VERSION)),
        "inst/extdata/simulation_crt76t.rds", compress = "xz")
message("saved inst/extdata/simulation_crt76t.rds")
