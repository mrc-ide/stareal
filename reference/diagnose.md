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
#>   variable   mean    sd   q2.5    q50 q97.5  rhat ess_bulk ess_tail
#>      sigma  1.240 0.268  0.618  1.220 1.790 1.209        8       53
#>      rho_s  0.609 0.206  0.196  0.603 0.972 1.058       24       14
#>      rho_t  0.859 0.091  0.623  0.874 0.977 1.431        5       19
#>    sigma_n  0.309 0.132  0.129  0.277 0.579 1.070       26       53
#>  intercept -0.549 0.681 -1.950 -0.556 0.768 1.027       38       34
#>      trend  0.956 0.368  0.446  0.966 1.600 1.114       20       56
#> 
#> acceptance after burn-in: centred 0.17, surrogate 0.20; 0.00071 s per iteration
#> latent cells: 30; R-hat over threshold in 10; bulk ESS median 47, min 16; max MCSE of the mean 1.5 (prevalence, %)
#> 
#>   theta R-hat (max)              1.43  <= 1.05   FAIL
#>   beta R-hat (max)               1.11  <= 1.05   FAIL
#>   latent R-hat (max)              1.2  <= 1.05   FAIL
#>   theta bulk ESS (min)           5.11  >= 400    FAIL
#>   beta bulk ESS (min)            19.9  >= 400    FAIL
#>   latent bulk ESS (median)       46.8  >= 100    FAIL
#> 
#> FAIL: run longer chains (see ?mcmc_control)
```
