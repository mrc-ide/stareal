# Simulate survey data from the model

Draws one realisation of the latent surface on a frame and time frame,
then binomial surveys at randomly chosen area-time cells (one survey per
cell). The logit prevalence is the fixed effects set by `formula` (by
default an intercept and a linear trend running from -1 at the first
time step to 1 at the last) + a Leroux x AR1 field, with a cell-level
nugget added at each surveyed cell.

## Usage

``` r
simulate_areal(
  frame,
  time,
  formula = ~trend(),
  covariates = NULL,
  scale_covariates = TRUE,
  beta = c(-1, 1),
  sigma = 1,
  rho_s = 0.9,
  rho_t = 0.9,
  sigma_n = 0.3,
  n_cells = 100,
  k = 50,
  seed = NULL,
  cells = NULL
)
```

## Arguments

- frame:

  A
  [`areal_frame()`](https://mrc-ide.github.io/stareal/reference/areal_frame.md).

- time:

  A
  [`time_frame()`](https://mrc-ide.github.io/stareal/reference/time_frame.md).

- formula, covariates, scale_covariates:

  The fixed effects, as in
  [`fit_areal()`](https://mrc-ide.github.io/stareal/reference/fit_areal.md).
  Simulating and fitting with the same settings gives effects on the
  same scale.

- beta:

  Fixed effects on the logit scale: one per design column, in order
  (intercept, `trend` if in the formula, then the covariates; see
  [`areal_priors()`](https://mrc-ide.github.io/stareal/reference/areal_priors.md)),
  or named by design column.

- sigma:

  Scale of the field: the standard deviation of the Leroux prior is
  `sigma` times a factor set by the graph (\\\tau = 1/\sigma^2\\).

- rho_s:

  Spatial dependence of the Leroux prior, in (0, 1).

- rho_t:

  Lag-one temporal correlation of the AR1, in (0, 1).

- sigma_n:

  Standard deviation of the cell-level nugget (0 for none).

- n_cells:

  Number of distinct area-time cells surveyed, chosen at random.

- k:

  Number tested per survey: one number, or one per cell.

- seed:

  Random seed. If given, the caller's random-number stream is left
  untouched.

- cells:

  Instead of `n_cells` random cells: the cells to survey, as a data
  frame with an area `id` column, a `t` (step number) or `time` (label)
  column, and optionally `k` (number tested). The `$cells` of an
  [`areal_data()`](https://mrc-ide.github.io/stareal/reference/areal_data.md)
  object qualifies, so a simulation can copy a real survey design: its
  areas, time steps and sample sizes.

## Value

A
[`areal_data()`](https://mrc-ide.github.io/stareal/reference/areal_data.md)
object (`stareal_data`), ready for
[`fit_areal()`](https://mrc-ide.github.io/stareal/reference/fit_areal.md),
with `$truth`: `theta` (natural scale), `beta` (named by design column),
`logit` and `prevalence` (areas x time steps matrices of the surface
without the nugget; prevalence in percent), `x` (the logit at each
surveyed cell with its nugget, in the order of `$cells`: what a fit
estimates for the observed cells) and `design`. `$surveys` has columns
`id`, `t`, `time`, `date` (calendar frames), `y`, `k`.

## Examples

``` r
fr <- areal_lattice(6, 6)
tf <- time_frame(2001, 2020, "year")
sim <- simulate_areal(fr, tf, n_cells = 60, seed = 1)
sim
#> <stareal_data> 60 of 60 surveys used, pooled into 60 cells (simulated)
#>   29 of 36 units and 20 of 20 time steps have data (8.33% of the 720 area-time cells)
#>   tested per cell: median 50, range 50-50; total 3,000 tested, 710 positive
#>   truth: sigma 1, rho_s 0.9, rho_t 0.9, sigma_n 0.3; intercept -1, trend 1
```
