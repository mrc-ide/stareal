# Approximate leave-one-out cross-validation of a fit

Estimates how well the model predicts each observed cell when that cell
is left out, from the one fit: no refits. Each observed cell's own
latent value is integrated out exactly before Pareto-smoothed importance
sampling (PSIS; Vehtari et al. 2016, 2017), which makes the
approximation stable where plain PSIS-LOO on this model is not. The
result is a `loo` object, so
[`loo::loo_compare()`](https://mc-stan.org/loo/reference/loo_compare.html)
and the `loo` package's other tools work on it;
[`compare_models()`](https://mrc-ide.github.io/stareal/reference/compare_models.md)
puts several fits side by side.

## Usage

``` r
# S3 method for class 'stareal_fit'
loo(x, n_draws = NULL, verbose = FALSE, ...)
```

## Arguments

- x:

  A
  [`fit_areal()`](https://mrc-ide.github.io/stareal/reference/fit_areal.md)
  object.

- n_draws:

  Number of stored draws to use, taken evenly from every chain (default:
  all).

- verbose:

  Report progress.

- ...:

  Not used.

## Value

A `psis_loo` object from the `loo` package (see
[`loo::loo()`](https://mc-stan.org/loo/reference/loo.html)), with class
`stareal_loo` added and an extra table `cells`: one row per observed
cell, with `cell_id`, `id`, `t`, `time`, `y`, `k`, `elpd_loo` (its log
predictive density when left out) and `pareto_k`.

## Details

The unit left out is an observed cell (all the surveys of one area and
time step), predicted as a new survey of the cell would be, nugget
included, as in
[`crossval()`](https://mrc-ide.github.io/stareal/reference/crossval.md).
Leaving one cell out keeps its neighbours in, so this measures how well
the model fills gaps between data; it says little about extrapolation to
areas with no data. For that, use
[`crossval()`](https://mrc-ide.github.io/stareal/reference/crossval.md)
with folds by area or region.

Cells with a high Pareto \\\hat k\\ (above 0.7, or less with fewer than
about 2,200 draws: the `loo` package's threshold \\\min(1 - 1/\log\_{10}
S, 0.7)\\ for \\S\\ draws) have unreliable estimates; the print-out
counts them. A few are normal for cells far from the rest of the data,
or with very large samples. If there are many, check against
[`crossval()`](https://mrc-ide.github.io/stareal/reference/crossval.md),
or use more draws.

## References

Vehtari, A., Mononen, T., Tolvanen, V., Sivula, T. and Winther, O.
(2016). Bayesian leave-one-out cross-validation approximations for
Gaussian latent variable models. *Journal of Machine Learning Research*
17(103), 1-38.

Vehtari, A., Gelman, A. and Gabry, J. (2017). Practical Bayesian model
evaluation using leave-one-out cross-validation and WAIC. *Statistics
and Computing* 27, 1413-1432.

## Examples

``` r
sim <- simulate_areal(areal_lattice(4, 4), time_frame(n_steps = 5), n_cells = 40, seed = 1)
fit <- fit_areal(sim, control = list(n_iter = 200, n_burn = 100, n_store = 50, n_chains = 2, n_cores = 1,
                                     verbose = FALSE))
loo(fit)
#> <stareal_loo> integrated PSIS-LOO over 40 observed cells, 100 draws
#>   elpd_loo -126.2 (se 4.0); p_loo 4.2; 8 cells with Pareto k > 0.50
```
