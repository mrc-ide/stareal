# Fit the space-time areal model

Fits the model to the observed area-time cells of `data` by MCMC. The
logit prevalence of cell (area \\a\\, time \\t\\) is \$\$x\_{at} =
z\_{at}^\top \beta + f\_{at} + \epsilon\_{at},\$\$ where \\z\_{at}\\
holds the fixed effects set by `formula` (by default an intercept and a
linear trend in time), \\f\\ is a Gaussian field with separable
covariance (a Leroux prior on the frame's graph in space, with scale
`sigma` and dependence `rho_s`; a unit-variance AR1 with correlation
`rho_t` in time), and \\\epsilon\\ is a cell-level nugget with standard
deviation `sigma_n`. Positives are binomial given \\x\\.

## Usage

``` r
fit_areal(
  data,
  formula = ~trend(),
  covariates = NULL,
  scale_covariates = TRUE,
  nugget = TRUE,
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

  the hyperparameter trace, burn-in included: chain, iteration, the
  hyperparameters on the unconstrained scale (`log_sigma`,
  `logit_rho_s`, `logit_rho_t`, `log_sigma_n`), `lp` (the log prior of
  the hyperparameters plus the log density of the latent values given
  them; `lp + ll` is the unnormalised log joint), the log-likelihood
  `ll`, and the acceptance rates of the two hyperparameter moves

- draws:

  the stored posterior draws, as arrays (draw x chain x variable):
  `theta` (unconstrained), `beta`, and `x` (logit prevalence of each
  observed cell, nugget included, in the order of `data$cells`)

- stored_iter, seconds:

  the iterations the draws were stored at; run time per chain

- design:

  the fixed-effect design: formula, names, and the centre and scale of
  each standardised covariate

- data, priors, control, nugget, checkpoint_dir, call:

  as given (priors with one prior per fixed effect)

Use
[summary()](https://mrc-ide.github.io/stareal/reference/summary.stareal_fit.md),
[`diagnose()`](https://mrc-ide.github.io/stareal/reference/diagnose.md),
[plot()](https://mrc-ide.github.io/stareal/reference/plot.stareal_fit.md),
[predict()](https://mrc-ide.github.io/stareal/reference/predict.stareal_fit.md)
and `posterior::as_draws_df(fit)` (or `as_draws_array()`, ...).

## Details

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
#>   formula ~trend(): intercept, trend
#>   1 chain x 100 iterations (50 burn-in); 25 stored draws in all; 0.064 s per chain
#>   posterior medians: sigma 1.33, rho_s 0.603, rho_t 0.905, sigma_n 0.252
#>   see summary(), diagnose(), plot() and predict()
```
