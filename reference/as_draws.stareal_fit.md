# Posterior draws of a fit

The stored draws as a `posterior` draws object: the hyperparameters
(`sigma`, `rho_s`, `rho_t`, `sigma_n`), the fixed effects, and, if
`latent`, the logit prevalence of each observed cell, nugget included,
named `logit[<cell_id>]`. All the `posterior` conversions work:
`as_draws_df()`, `as_draws_array()`, `as_draws_matrix()`,
`as_draws_list()`, `as_draws_rvars()`.

## Usage

``` r
# S3 method for class 'stareal_fit'
as_draws(x, latent = TRUE, ...)

# S3 method for class 'stareal_fit'
as_draws_array(x, latent = TRUE, ...)

# S3 method for class 'stareal_fit'
as_draws_df(x, latent = TRUE, ...)

# S3 method for class 'stareal_fit'
as_draws_matrix(x, latent = TRUE, ...)

# S3 method for class 'stareal_fit'
as_draws_list(x, latent = TRUE, ...)

# S3 method for class 'stareal_fit'
as_draws_rvars(x, latent = TRUE, ...)
```

## Arguments

- x:

  A `stareal_fit`.

- latent:

  Include the observed cells' latent values.

- ...:

  Not used.

## Value

A `posterior` draws object.
