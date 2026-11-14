# Check convergence and mixing of a fit

[summary()](https://mrc-ide.github.io/stareal/reference/summary.stareal_fit.md)
with every observed cell included and a verdict. Rank-normalised R-hat,
bulk and tail effective sample sizes (ESS) and Monte Carlo standard
errors are computed for every parameter: the hyperparameters (from the
full trace after burn-in), the fixed effects, and the latent prevalence
of every observed cell, nugget included (from the stored draws; this is
not the nugget-free surface of
[predict()](https://mrc-ide.github.io/stareal/reference/predict.stareal_fit.md)).
A fit passes only if all are acceptable. Hyperparameters can look
converged while the latent cells mix poorly, so the latent cells are
always checked.

## Usage

``` r
diagnose(
  fit,
  rhat_max = 1.05,
  ess_theta_min = 400,
  ess_latent_median_min = 100
)
```

## Arguments

- fit:

  A
  [`fit_areal()`](https://mrc-ide.github.io/stareal/reference/fit_areal.md)
  object.

- rhat_max:

  Largest acceptable R-hat, for every parameter.

- ess_theta_min:

  Smallest acceptable bulk ESS of each hyperparameter and each fixed
  effect.

- ess_latent_median_min:

  Smallest acceptable median bulk ESS over the observed cells.

## Value

A `summary.stareal_fit` (see
[summary()](https://mrc-ide.github.io/stareal/reference/summary.stareal_fit.md))
with the latent rows, and `checks` (a data frame: check, value,
threshold, pass), `pass` (`TRUE` if every check passes) and `latent` (a
few statistics over the cells).

## Examples

``` r
fr <- areal_lattice(4, 4)
tf <- time_frame(n_steps = 5)
sim <- simulate_areal(fr, tf, n_cells = 30, seed = 1)
fit <- fit_areal(sim, control = list(n_iter = 100, n_burn = 50, n_store = 25, n_chains = 2,
                                     n_cores = 1, ess_sweeps = 2, verbose = FALSE))
suppressWarnings(diagnose(fit))   # far too short to pass
#> <stareal_fit summary> 30 observed cells; 2 chains, 50 iterations after burn-in
#>   hyperparameters from the trace (100 draws), fixed effects and cells from the 50 stored draws
#> 
#>   variable   mean    sd     q2.5    q50 q97.5  rhat ess_bulk ess_tail
#>      sigma  1.200 0.297  0.55200  1.170 1.610 1.122       12       46
#>      rho_s  0.648 0.186  0.23600  0.679 0.926 1.045       43       55
#>      rho_t  0.818 0.201  0.19500  0.920 0.985 1.438       18       23
#>    sigma_n  0.319 0.186  0.06350  0.248 0.765 1.084       24       23
#>  intercept -0.539 0.455 -1.32000 -0.520 0.298 1.073       81       56
#>      trend  0.791 0.345  0.00924  0.870 1.340 1.009       46       18
#> 
#> acceptance after burn-in: centred 0.12, surrogate 0.20; 0.0006 s per iteration
#> latent cells: 30; R-hat over threshold in 10; bulk ESS median 56, min 11; max MCSE of the mean 2.2 (prevalence, %)
#> 
#>   theta R-hat (max)              1.44  <= 1.05   FAIL
#>   beta R-hat (max)               1.07  <= 1.05   FAIL
#>   latent R-hat (max)             1.19  <= 1.05   FAIL
#>   theta bulk ESS (min)             12  >= 400    FAIL
#>   beta bulk ESS (min)            45.8  >= 400    FAIL
#>   latent bulk ESS (median)         56  >= 100    FAIL
#> 
#> FAIL: run longer chains (see ?mcmc_control)
```
