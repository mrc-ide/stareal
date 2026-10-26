# Summarise a fit

Posterior summaries of the hyperparameters and fixed effects, with R-hat
and effective sample sizes (ESS). The hyperparameters are summarised
from the full trace after burn-in; the fixed effects (and, with
`latent = TRUE`, the observed cells) from the stored draws. For the full
check of convergence, use
[`diagnose()`](https://mrc-ide.github.io/stareal/reference/diagnose.md).

## Usage

``` r
# S3 method for class 'stareal_fit'
summary(object, latent = FALSE, ...)
```

## Arguments

- object:

  A `stareal_fit`.

- latent:

  Also summarise each observed cell's latent prevalence, in percent
  (nugget included).

- ...:

  Not used.

## Value

An object of class `summary.stareal_fit`: `table` (one row per
parameter: `type` (`theta`, `beta` or `latent`), `variable`, `mean`,
`sd`, `q2.5`, `q50`, `q97.5`, `rhat`, `ess_bulk`, `ess_tail`,
`mcse_mean`, `n_draws`, and for latent rows `cell_id`, `id` and `time`),
`acceptance` (of the two hyperparameter moves, after burn-in),
`seconds_per_iteration`, and, from
[`diagnose()`](https://mrc-ide.github.io/stareal/reference/diagnose.md),
`checks` and `pass`.
