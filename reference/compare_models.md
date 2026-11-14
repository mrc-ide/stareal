# Compare models fitted to the same data

One row per fit: its model, whether it passes
[`diagnose()`](https://mrc-ide.github.io/stareal/reference/diagnose.md),
its approximate leave-one-out score from
[loo()](https://mrc-ide.github.io/stareal/reference/loo.stareal_fit.md),
and, if given, its
[`crossval()`](https://mrc-ide.github.io/stareal/reference/crossval.md)
scores. All fits must be made to the same observed cells (the same data,
whatever the formula, spatial prior or nugget).

## Usage

``` r
compare_models(..., cv = NULL, loo = NULL)
```

## Arguments

- ...:

  Two or more
  [`fit_areal()`](https://mrc-ide.github.io/stareal/reference/fit_areal.md)
  objects, named (e.g. `compare_models(leroux = fit1, bym2 = fit2)`);
  unnamed ones are named by their expressions.

- cv:

  Optional:
  [`crossval()`](https://mrc-ide.github.io/stareal/reference/crossval.md)
  results, a list named as the fits, with the same folds and seed for
  every model (see
  [`compare_cv()`](https://mrc-ide.github.io/stareal/reference/compare_cv.md)).

- loo:

  Optional:
  [loo()](https://mrc-ide.github.io/stareal/reference/loo.stareal_fit.md)
  results already computed, a list named as the fits (any fit left out
  is computed here).

## Value

A data frame, best first (by `elpd_loo`): `model`, `spatial` (the
spatial prior and any fixed hyperparameters), `nugget`, `n_hyper` (free
hyperparameters), `n_beta` (fixed effects), `pass`
([`diagnose()`](https://mrc-ide.github.io/stareal/reference/diagnose.md)),
`elpd_loo`, `se_elpd_loo`, `elpd_diff`, `se_diff`, `p_loo`, `n_k_high`
(cells with a high Pareto k, see
[loo()](https://mrc-ide.github.io/stareal/reference/loo.stareal_fit.md));
and with `cv`: `cv_log_score`, `cv_diff`, `cv_se_diff`,
`cv_pct_outside_95`, `cv_rmse`. The `loo` results are attached as
attribute `loo`.

## Details

Scores are expected log predictive densities: higher is better.
`elpd_diff` is each model's difference from the best (0 for the best,
negative for the others), summed over cells, and `se_diff` its standard
error from the cell-by-cell differences, as in
[`loo::loo_compare()`](https://mc-stan.org/loo/reference/loo_compare.html).
A difference smaller than about two standard errors is not a clear
preference. The cross-validation columns are per cell (`cv_log_score` is
the mean log score of
[`compare_cv()`](https://mrc-ide.github.io/stareal/reference/compare_cv.md)),
as there.

Why not AIC, BIC or DIC: AIC needs a number of parameters, which a model
with a latent value per cell does not have (the effective number,
`p_loo` here, depends on the hyperparameters); BIC approximates a
marginal likelihood that, in a hierarchical model, is driven by the
priors; and DIC is unreliable for latent Gaussian models. Predictive
scores avoid all three problems.

## Examples

``` r
# \donttest{
sim <- simulate_areal(areal_lattice(4, 4), time_frame(n_steps = 5), n_cells = 40, seed = 1)
ctl <- list(n_iter = 200, n_burn = 100, n_store = 50, n_chains = 2, n_cores = 1, verbose = FALSE)
compare_models(leroux = fit_areal(sim, control = ctl),
               bym2 = fit_areal(sim, spatial = "bym2", control = ctl),
               no_nugget = fit_areal(sim, nugget = FALSE, control = ctl))
#>       model spatial nugget n_hyper n_beta  pass  elpd_loo se_elpd_loo
#> 1      bym2    bym2   TRUE       4      2 FALSE -125.4016    4.186919
#> 2    leroux  leroux   TRUE       4      2 FALSE -126.1897    3.992879
#> 3 no_nugget  leroux  FALSE       3      2 FALSE -126.5311    4.311237
#>    elpd_diff  se_diff    p_loo n_k_high
#> 1  0.0000000 0.000000 4.134368        4
#> 2 -0.7880569 1.502263 4.249251        8
#> 3 -1.1295210 1.500586 6.349915       10
# }
```
