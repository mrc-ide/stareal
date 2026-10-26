# crossval_crt76t.R
#
# Makes inst/extdata/crossval_crt76t.rds: the cross-validation results of
# the tutorial "Cross-validation and model comparison", which imports them.
# Cross-validating five models with ten folds, plus country folds for one,
# is 58 refits: too long to run inside a tutorial, so it is run once here.
#
# The models are fitted to the crt K76T data of tutorial 1 (setup_crt76t.R):
#   default        ~ trend()                   (tutorial 1)
#   no_nugget      ~ trend(), nugget = FALSE
#   no_trend       ~ 1
#   smooth_trend   ~ s1 + s2 + s3             (tutorial 3: a spline in time)
#   kenya_uganda   ~ trend() + ku + ku_trend  (tutorial 3: Kenya and Uganda)
# Every refit uses the default sampler settings (mcmc_control()), the length
# of the tutorials' fits, and the folds come from seed 1, so compare_cv()
# can line the models up. The default model is also cross-validated with
# one fold per country.
#
# Rerun when the sampler, the model or the data change. From the package
# root:
#   Rscript data-raw/crossval_crt76t.R
#
# ------------------------------------------------------------------

source("data-raw/setup_crt76t.R")

# the covariates of tutorial 3, built as there
year <- 1989 + (seq_len(tf$n_steps) - 1) / 12
basis <- splines::ns(year, df = 3)
by_time <- data.frame(time = tf$labels, s1 = basis[, 1], s2 = basis[, 2], s3 = basis[, 3])
by_area <- data.frame(id = fr$ids, ku = as.numeric(fr$groups$country %in% c("Kenya", "Uganda")))
by_cell <- expand.grid(id = fr$ids, time = tf$labels, stringsAsFactors = FALSE)
by_cell$ku_trend <- by_area$ku[match(by_cell$id, by_area$id)] *
  seq(-1, 1, length.out = tf$n_steps)[match(by_cell$time, tf$labels)]

quiet <- list(verbose = FALSE)
models <- list(
  default      = fit_areal(d, control = quiet),
  no_nugget    = fit_areal(d, nugget = FALSE, control = quiet),
  no_trend     = fit_areal(d, formula = ~ 1, control = quiet),
  smooth_trend = fit_areal(d, ~ s1 + s2 + s3, covariates = by_time, scale_covariates = FALSE, control = quiet),
  kenya_uganda = fit_areal(d, ~ trend() + ku + ku_trend, covariates = list(by_area, by_cell),
                           scale_covariates = FALSE, control = quiet)
)
checks <- lapply(models, function(f) diagnose(f)$checks)
cv <- lapply(names(models), function(m) {
  message("cross-validating ", m)
  crossval(models[[m]], folds = 10, seed = 1)
})
names(cv) <- names(models)
country <- fr$groups$country[d$cells$area]
cv_country <- crossval(models$default, folds = country, seed = 1)

dir.create("inst/extdata", showWarnings = FALSE, recursive = TRUE)
saveRDS(list(cv = cv, cv_country = cv_country, checks = checks,
             built = list(date = Sys.Date(), stareal = as.character(utils::packageVersion("stareal")),
                          sampler = SAMPLER_VERSION)),
        "inst/extdata/crossval_crt76t.rds", compress = "xz")
message(sprintf("saved inst/extdata/crossval_crt76t.rds (%.0f kB)", file.size("inst/extdata/crossval_crt76t.rds") / 1e3))
