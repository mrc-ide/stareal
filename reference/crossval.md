# Cross-validate a fit

Splits the observed cells at random into `folds` groups. Each group is
held out in turn: the model is refitted to the other cells, with the
fit's formula, covariates, priors, nugget and sampler settings, and then
predicts the held-out cells. Use it to check calibration out of sample
and to compare models fitted to the same cells. The same `seed` gives
the same folds, so models are compared on the same splits. The mean log
score is higher for the better model; about 5% of observations should
fall outside their 95% predictive interval.

## Usage

``` r
crossval(
  fit,
  folds = 10,
  fold = NULL,
  control = fit$control,
  checkpoint_dir = NULL,
  seed = 1
)
```

## Arguments

- fit:

  A
  [`fit_areal()`](https://mrc-ide.github.io/stareal/reference/fit_areal.md)
  object.

- folds:

  The number of folds, assigned at random; or a fold label for each
  observed cell, in the order of `fit$data$cells`, for folds of your
  own: by site, by country (spatial cross-validation, which tests
  extrapolation to unsurveyed places), or by period.

- fold:

  Which folds to run (default all), as numbers or labels, e.g. to run
  them in separate sessions.

- control:

  Sampler settings for the refits:
  [`mcmc_control()`](https://mrc-ide.github.io/stareal/reference/mcmc_control.md)
  or a list of some of its arguments (default: the fit's).

- checkpoint_dir:

  A folder for the refits' checkpoints: fold `k` of `K` with split seed
  `s` checkpoints in `<checkpoint_dir>/cv_seed<s>_fold_<k>_of_<K>`.
  `NULL`: none.

- seed:

  Seed of the random split into folds.

## Value

An object of class `stareal_cv`:

- cells:

  one row per held-out cell: `fold`, `cell_id`, `id`, `t`, `time`, `y`,
  `k`, `obs_prev` (`100 y / k`), `pred_mean`, `pred_lo`, `pred_hi` (mean
  and 95% interval of the cell's prevalence, with the nugget), `rep_lo`,
  `rep_hi` (95% posterior-predictive interval of `100 y / k`); all
  prevalences in percent; `pit` (randomised probability integral
  transform of `y`), `outside` (`pit` outside 0.025-0.975), and
  `log_score`, the log predictive density of `y`

- folds:

  one row per fold, and `overall`: the number held out, RMSE of
  `pred_mean` against `obs_prev` (percentage points), percent outside
  the 95% interval, mean log score and correlation

- fold_of:

  the fold of every observed cell

Compare models with
[`compare_cv()`](https://mrc-ide.github.io/stareal/reference/compare_cv.md).

## Details

Held-out cells are predicted as observations, with the cell nugget, so
`pred_*` describe the prevalence a new survey of the cell would see, not
the smoother surface of
[predict()](https://mrc-ide.github.io/stareal/reference/predict.stareal_fit.md).

Each fold is a full refit, so this costs `folds` (or `length(fold)`)
times the fit.

## Examples

``` r
# \donttest{
fr <- areal_lattice(4, 4)
tf <- time_frame(n_steps = 5)
sim <- simulate_areal(fr, tf, n_cells = 40, seed = 1)
fit <- fit_areal(sim, control = list(n_iter = 100, n_burn = 50, n_store = 25, n_chains = 1,
                                     ess_sweeps = 2, verbose = FALSE))
cv <- crossval(fit, folds = 4)
cv
#> <stareal_cv> 4-fold cross-validation (seed 1), 4 folds run, 40 cells held out
#>   RMSE 13.75 points; 2.5% outside the 95% predictive interval; mean log score -3.254; correlation 0.51
# }
```
