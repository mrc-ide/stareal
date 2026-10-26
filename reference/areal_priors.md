# Priors for an areal model

Priors are on the natural scale of each parameter. The field has
precision \\\tau = 1/\sigma^2\\; the fixed effects are integrated out
analytically, so their prior must be Gaussian.

## Usage

``` r
areal_priors(
  sigma_sd = 2,
  rho_s = c(2, 1),
  rho_t = c(1, 1),
  sigma_n_sd = 1,
  beta_mean = 0,
  beta_sd = 3
)
```

## Arguments

- sigma_sd:

  Scale of the half-normal prior on `sigma`, the scale of the space-time
  field: \\\sigma \sim N^+(0, \code{sigma\\sd}^2)\\.

- rho_s:

  Shape parameters of the Beta prior on `rho_s`, the spatial dependence
  of the Leroux prior (0: independent areas, 1: intrinsic CAR). The
  default Beta(2, 1) favours strong spatial dependence.

- rho_t:

  Shape parameters of the Beta prior on `rho_t`, the lag-one correlation
  of the AR1 in time. Default uniform. With one time step the data say
  nothing about `rho_t`, and its posterior is this prior.

- sigma_n_sd:

  Scale of the half-normal prior on `sigma_n`, the standard deviation of
  the cell-level nugget.

- beta_mean, beta_sd:

  Means and standard deviations of the independent normal priors on the
  fixed effects (logit scale). Give one value (used for every fixed
  effect), one value per fixed effect in the order of the design, or a
  vector named by design column, e.g. `c(intercept = 0, rain = 0)`,
  where unnamed columns take the default. The design columns are the
  intercept, `trend` (if
  [`trend()`](https://mrc-ide.github.io/stareal/reference/trend.md) is
  in the formula), then the covariates in the order of the formula, a
  factor giving one column per level after the first (named by the
  covariate and the level, e.g. `zonehigh`). They are listed in
  `fit$design$names`.

## Value

A list of the settings.

## Details

[`fit_areal()`](https://mrc-ide.github.io/stareal/reference/fit_areal.md)
also accepts a plain list of any of these arguments, e.g.
`priors = list(sigma_sd = 1)`; the others take their defaults.

## Examples

``` r
areal_priors()
#> $sigma_sd
#> [1] 2
#> 
#> $rho_s
#> [1] 2 1
#> 
#> $rho_t
#> [1] 1 1
#> 
#> $sigma_n_sd
#> [1] 1
#> 
#> $beta_mean
#> [1] 0
#> 
#> $beta_sd
#> [1] 3
#> 
areal_priors(sigma_sd = 1, rho_t = c(4, 1))
#> $sigma_sd
#> [1] 1
#> 
#> $rho_s
#> [1] 2 1
#> 
#> $rho_t
#> [1] 4 1
#> 
#> $sigma_n_sd
#> [1] 1
#> 
#> $beta_mean
#> [1] 0
#> 
#> $beta_sd
#> [1] 3
#> 
```
