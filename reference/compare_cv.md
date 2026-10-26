# Compare cross-validated models

Compares models cross-validated on the same cells with the same folds
(the same `folds` and `seed` in
[`crossval()`](https://mrc-ide.github.io/stareal/reference/crossval.md)).
The log scores are compared cell by cell, which removes most of the
variation between cells: the difference from the best model is given
with its standard error, the standard deviation of the per-cell
differences over the square root of the number of cells. A difference of
less than about two standard errors is not a clear preference.

## Usage

``` r
compare_cv(...)
```

## Arguments

- ...:

  Two or more `stareal_cv` objects, named (e.g.
  `compare_cv(default = cv1, no_nugget = cv2)`); unnamed ones are named
  by their expressions.

## Value

A data frame, one row per model, best first: `model`, `n_cells`,
`mean_log_score`, `diff` (from the best model), `se_diff`,
`pct_outside_95`, `rmse`, `correlation`.

## Examples

``` r
# \donttest{
sim <- simulate_areal(areal_lattice(4, 4), time_frame(n_steps = 5), n_cells = 40, seed = 1)
ctl <- list(n_iter = 100, n_burn = 50, n_store = 25, n_chains = 1, ess_sweeps = 2, verbose = FALSE)
fit1 <- fit_areal(sim, control = ctl)
fit2 <- fit_areal(sim, nugget = FALSE, control = ctl)
compare_cv(nugget = crossval(fit1, folds = 4), no_nugget = crossval(fit2, folds = 4))
#>       model n_cells mean_log_score        diff    se_diff pct_outside_95
#> 1 no_nugget      40      -3.213514  0.00000000 0.00000000            5.0
#> 2    nugget      40      -3.254170 -0.04065535 0.03941771            2.5
#>       rmse correlation
#> 1 13.44538   0.5446642
#> 2 13.75298   0.5074830
# }
```
