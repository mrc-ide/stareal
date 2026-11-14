# Fit the space-time areal model

Fits the model to the observed area-time cells of `data` by MCMC. The
logit prevalence of cell (area \\a\\, time \\t\\) is \$\$x\_{at} =
z\_{at}^\top \beta + f\_{at} + \epsilon\_{at},\$\$ where \\z\_{at}\\
holds the fixed effects set by `formula` (by default an intercept and a
linear trend in time), \\f\\ is a Gaussian field with separable
covariance (a conditional autoregressive prior on the frame's graph in
space, with scale `sigma` and dependence `rho_s`; a unit-variance AR1
with correlation `rho_t` in time), and \\\epsilon\\ is a cell-level
nugget with standard deviation `sigma_n`. Positives are binomial given
\\x\\.

## Usage

``` r
fit_areal(
  data,
  formula = ~trend(),
  covariates = NULL,
  scale_covariates = TRUE,
  nugget = TRUE,
  spatial = "leroux",
  main = "none",
  interaction = "IV",
  temporal = "ar1",
  fixed = NULL,
  priors = areal_priors(),
  control = mcmc_control(),
  checkpoint_dir = NULL
)
```

## Arguments

- data:

  A
  [`areal_data()`](https://mrc-ide.github.io/stareal/reference/areal_data.md)
  object (or a simulation from
  [`simulate_areal()`](https://mrc-ide.github.io/stareal/reference/simulate_areal.md)).

- formula:

  A one-sided formula for the fixed effects. It always has an intercept;
  [`trend()`](https://mrc-ide.github.io/stareal/reference/trend.md) adds
  a linear trend in time (from -1 at the first time step to 1 at the
  last); other terms are names of columns of `covariates`. Factor or
  character covariates get one dummy per level after the first.
  Transforms and interactions must be made as columns.

- covariates:

  A data frame, or a list of data frames, of covariates. Each is keyed
  by an `id` column (areas, for covariates that do not change over
  time), a `t` (step number) or `time` (label) column (for covariates
  the same in every area), or both (covariates by area and time). Each
  must cover every area, time step or area-time cell it is keyed by,
  since predictions are made everywhere.

- scale_covariates:

  Standardise numeric covariates to mean 0 and sd 1, over the areas,
  time steps or cells they are given for, so that the default prior on
  their effects (sd 3, see
  [`areal_priors()`](https://mrc-ide.github.io/stareal/reference/areal_priors.md))
  is weakly informative whatever their units. Effects are then per sd of
  the covariate; the centres and scales are in `fit$design` and in
  [`summary()`](https://rdrr.io/r/base/summary.html). Unscaled
  covariates in large units (e.g. population counts) can make the prior
  variance of the mean so large that the field is lost in rounding; the
  fit then stops and says so.

- nugget:

  Include the cell-level nugget.

- spatial:

  The spatial prior: `"leroux"`, `"bym2"`, `"car"`, `"iid"` or `"icar"`
  (see "Spatial priors"), for every component that is structured in
  space. `"iid"` and `"icar"` fix `rho_space` too.

- main:

  Main effects: `"none"`, `"space"`, `"time"` or `"both"` (see "Main
  effects and interaction types").

- interaction:

  The space-time interaction: `"IV"`, `"III"`, `"II"`, `"I"` or
  `"none"`.

- temporal:

  The prior of the time main effect: `"ar1"`, `"rw1"` or `"rw2"`.

- fixed:

  Hyperparameters held at given values instead of sampled: a named list
  (or vector) with any of the model's hyperparameters: `sigma`, `rho_s`
  (between 0 and 1; below 1 for `"car"`), `rho_t` (at least 0 and below
  1), `sigma_n` (with `nugget = TRUE`), and with main effects
  `sigma_space`, `rho_space`, `sigma_time` and `rho_time`. At least one
  hyperparameter must stay free.

- priors:

  [`areal_priors()`](https://mrc-ide.github.io/stareal/reference/areal_priors.md),
  or a list of some of its arguments.

- control:

  [`mcmc_control()`](https://mrc-ide.github.io/stareal/reference/mcmc_control.md),
  or a list of some of its arguments.

- checkpoint_dir:

  A folder to save each chain's state in. If a fit is interrupted,
  running it again with the same folder resumes every chain where it
  stopped, and gives exactly the draws of an uninterrupted run; a
  finished fit run again is read back from the folder. A folder left by
  a different fit (other data, frame, model or sampler settings) is
  refused. Use a folder outside synced storage (Dropbox, OneDrive),
  which can lock files. `NULL`: no checkpoints.

## Value

An object of class `stareal_fit`:

- trace:

  the hyperparameter trace, burn-in included: chain, iteration, the free
  hyperparameters on the unconstrained scale (`log_sigma`,
  `logit_rho_s`, `logit_rho_t`, `log_sigma_n`, and with main effects
  `log_sigma_space`, `logit_rho_space`, `log_sigma_time`,
  `logit_rho_time`), `lp` (the log prior of the hyperparameters plus the
  log density of the latent values given them; `lp + ll` is the
  unnormalised log joint), the log-likelihood `ll`, and the acceptance
  rates of the two hyperparameter moves

- draws:

  the stored posterior draws, as arrays (draw x chain x variable):
  `theta` (unconstrained), `beta`, and `x` (logit prevalence of each
  observed cell, nugget included, in the order of `data$cells`)

- stored_iter, seconds:

  the iterations the draws were stored at; run time per chain

- design:

  the fixed-effect design: formula, names, and the centre and scale of
  each standardised covariate

- model:

  the model: `spatial` (as given), `family` (the prior it is:
  `"leroux"`, `"bym2"` or `"car"`), `fixed` (the fixed hyperparameters,
  natural scale, named), `main`, `interaction` and `temporal`

- data, priors, control, nugget, checkpoint_dir, call:

  as given (priors with one prior per fixed effect)

- version:

  what made the fit: the *stareal* version, its commit (`sha`, when
  installed from GitHub), the sampler version, R and the BLAS library

Use
[summary()](https://mrc-ide.github.io/stareal/reference/summary.stareal_fit.md),
[`diagnose()`](https://mrc-ide.github.io/stareal/reference/diagnose.md),
[plot()](https://mrc-ide.github.io/stareal/reference/plot.stareal_fit.md),
[predict()](https://mrc-ide.github.io/stareal/reference/predict.stareal_fit.md),
[loo()](https://mrc-ide.github.io/stareal/reference/loo.stareal_fit.md),
[`compare_models()`](https://mrc-ide.github.io/stareal/reference/compare_models.md)
and `posterior::as_draws_df(fit)` (or `as_draws_array()`, ...).

## Spatial priors

`spatial` chooses the prior on the graph. In each, `rho_s = 0` makes the
areas independent and `rho_s` near 1 makes neighbours alike (the
intrinsic CAR, ICAR, at 1). \\D\\ is the diagonal matrix of each area's
number of neighbours and \\W\\ the adjacency, so \\R = D - W\\ is the
graph Laplacian.

- `"leroux"`:

  (the default; Leroux et al. 2000) precision \\\tau(\rho_s R + (1 -
  \rho_s) I)\\, \\\sigma^2 = 1/\tau\\.

- `"bym2"`:

  (Riebler et al. 2016) covariance \\\sigma^2(\rho_s R\_\*^+ + (1 -
  \rho_s) I)\\: a mix of a scaled ICAR (\\R\_\*^+\\, the generalised
  inverse of \\R\\ scaled so that the geometric mean of its variances
  is 1) and independent noise. `sigma` is then roughly the marginal sd
  and `rho_s` the share of the variance that is spatially structured.
  Needs a connected graph.

- `"car"`:

  the proper CAR, precision \\\tau(D - \rho_s W)\\. Its marginal
  variances differ between areas (larger where there are fewer
  neighbours). `rho_s` must stay below 1, and every area needs a
  neighbour.

- `"iid"`:

  independent areas: the Leroux prior with `rho_s` fixed at 0.

- `"icar"`:

  the intrinsic CAR: the Leroux prior with `rho_s` fixed at 1. The field
  then sums to zero over the areas at every time step, so the regional
  level and its changes over time come from the fixed effects alone.
  Needs a connected graph.

Any hyperparameter can also be held at a value with `fixed`, e.g.
`fixed = list(rho_t = 0)` for a field with no memory from one time step
to the next. Fixed hyperparameters are not sampled and do not appear in
[summary()](https://mrc-ide.github.io/stareal/reference/summary.stareal_fit.md)
or
[`diagnose()`](https://mrc-ide.github.io/stareal/reference/diagnose.md).

## Main effects and interaction types

By default the field \\f\\ is one separable space-time component. `main`
and `interaction` make it a sum of independent components (Knorr-Held
2000, *Statistics in Medicine* 19, 2555-2567): \$\$f\_{at} =
f^{ST}\_{at} + f^S_a + f^T_t.\$\$

- `main`:

  `"none"` (the default), `"space"` (a spatial effect \\f^S\\, the same
  at every time step, under the prior `spatial`, with its own scale
  `sigma_space` and dependence `rho_space`), `"time"` (a temporal effect
  \\f^T\\, the same in every area, with scale `sigma_time`) or `"both"`.

- `temporal`:

  the prior of the time main effect: `"ar1"` (the default; lag-one
  correlation `rho_time`), `"rw1"` or `"rw2"` (first- or second-order
  random walks). The random walks are scaled so that `sigma_time` is
  roughly the marginal standard deviation whatever the number of time
  steps (Sorbye and Rue 2014). An RW1 sums to zero over time and an RW2
  also has zero slope, so the intercept and
  [`trend()`](https://mrc-ide.github.io/stareal/reference/trend.md)
  carry the level and the straight line, and the RW2 the curve about it:
  the natural pair. An RW1 or AR1 time effect overlaps
  [`trend()`](https://mrc-ide.github.io/stareal/reference/trend.md) in
  part, which is not an error but weakens what the trend coefficient
  means. Random walks are for the time main effect only.

- `interaction`:

  the space-time interaction \\f^{ST}\\, with scale `sigma`, by whether
  its time and space factors are structured: `"IV"` (the default; AR1 in
  time with `rho_t` x the spatial prior with `rho_s`), `"III"`
  (independent in time x the spatial prior), `"II"` (AR1 in time x
  independent areas), `"I"` (independent cells) or `"none"`. Type I
  needs `nugget = FALSE`: on the observed cells an independent cell
  effect is the nugget again. Type III with a nugget comes close to that
  when `rho_s` is small, and then mixes slowly: check
  [`diagnose()`](https://mrc-ide.github.io/stareal/reference/diagnose.md)
  and run longer.

With the defaults the model is exactly the separable one above. With the
ICAR, BYM2 or a random walk a component sums to zero, so the regional
level always comes from the fixed effects. More components mean more
hyperparameters, each less well identified than in a simpler model:
expect to need longer runs, and compare models with
[loo()](https://mrc-ide.github.io/stareal/reference/loo.stareal_fit.md)
or
[`crossval()`](https://mrc-ide.github.io/stareal/reference/crossval.md).

The unobserved cells and the fixed effects are integrated out exactly,
so the sampler works on the observed cells only; its cost per iteration
grows as the cube of the number of observed cells. The full area-time
surface is reconstructed afterwards by
[predict()](https://mrc-ide.github.io/stareal/reference/predict.stareal_fit.md).

## Examples

``` r
fr <- areal_lattice(4, 4)
tf <- time_frame(n_steps = 5)
sim <- simulate_areal(fr, tf, n_cells = 30, seed = 1)
fit <- fit_areal(sim, control = list(n_iter = 100, n_burn = 50, n_store = 25, n_chains = 1,
                                     ess_sweeps = 2, verbose = FALSE))
fit
#> <stareal_fit> 30 observed cells on 16 units x 5 time steps, with nugget
#>   spatial prior: leroux
#>   formula ~trend(): intercept, trend
#>   1 chain x 100 iterations (50 burn-in); 25 stored draws in all; 0.08 s per chain
#>   posterior medians: sigma 1.41, rho_s 0.69, rho_t 0.929, sigma_n 0.241
#>   see summary(), diagnose(), plot() and predict()

# the same data under a BYM2 prior, and with no memory in time
fit_bym2 <- fit_areal(sim, spatial = "bym2", control = list(n_iter = 100, n_burn = 50, n_store = 25,
                                                           n_chains = 1, ess_sweeps = 2, verbose = FALSE))
fit_iid_t <- fit_areal(sim, fixed = list(rho_t = 0), control = list(n_iter = 100, n_burn = 50, n_store = 25,
                                                                   n_chains = 1, ess_sweeps = 2, verbose = FALSE))

# a regional curve over time (an RW2 about the trend) plus the space-time field
fit_rw2 <- fit_areal(sim, main = "time", temporal = "rw2",
                     control = list(n_iter = 100, n_burn = 50, n_store = 25, n_chains = 1, ess_sweeps = 2,
                                    verbose = FALSE))
```
