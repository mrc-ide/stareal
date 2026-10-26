# Settings for the MCMC sampler

Each iteration makes `ess_sweeps` updates of the latent values of the
observed cells (elliptical slice sampling) and `theta_moves` pairs of
Metropolis moves on the hyperparameters. During burn-in, the proposals
adapt and the Gaussian approximation that steers the latent updates is
refitted to the chain (twice, at half and three-quarters of burn-in);
both are fixed afterwards, so the draws after burn-in come from an exact
MCMC. Chains start from dispersed hyperparameters and latent values, so
that R-hat in
[`diagnose()`](https://mrc-ide.github.io/stareal/reference/diagnose.md)
can reveal chains that have not mixed.

## Usage

``` r
mcmc_control(
  n_iter = 3000,
  n_burn = 1000,
  n_store = 250,
  n_chains = 4,
  ess_sweeps = 10,
  theta_moves = 2,
  theta_thin = 1,
  seed = 1,
  n_cores = NULL,
  blas_threads = 1,
  checkpoint_every = 500,
  verbose = TRUE
)
```

## Arguments

- n_iter:

  Iterations per chain, burn-in included.

- n_burn:

  Burn-in iterations (adaptation; discarded).

- n_store:

  Posterior draws kept per chain, spread evenly over the iterations
  after burn-in (the last iteration is always kept). These are the draws
  used for prediction.

- n_chains:

  Number of chains.

- ess_sweeps:

  Latent updates per iteration.

- theta_moves:

  Hyperparameter move pairs per iteration.

- theta_thin:

  Keep every `theta_thin`-th iteration of the hyperparameter trace
  (burn-in included). At least 4 traced iterations must fall after
  burn-in.

- seed:

  Random seed, any whole number \>= 0. Each chain gets its own seed,
  drawn from this one, and R's default generator is used whatever the
  session's [`RNGkind()`](https://rdrr.io/r/base/Random.html). The
  caller's random-number stream is left untouched.

- n_cores:

  Number of worker processes; chains run one per worker. Default: one
  per chain, at most the number of cores less one. Use 1 to run the
  chains one after another in this R session.

- blas_threads:

  BLAS threads per worker. One is safest; with an optimised BLAS (MKL,
  OpenBLAS, Accelerate) and large fits, more may be faster, as long as
  `n_cores * blas_threads` does not exceed the cores.

- checkpoint_every:

  Iterations between checkpoints, when
  [`fit_areal()`](https://mrc-ide.github.io/stareal/reference/fit_areal.md)
  is given a `checkpoint_dir`.

- verbose:

  Print progress: a line per chain after 10 iterations, then every tenth
  of the run (at most every 500 iterations), with the time left.

## Value

A list of the settings.

## Details

The defaults are a pilot run. A run for publication is typically much
longer (e.g. 40,000 iterations with 8,000 burn-in and
`theta_thin = 10`). Use
[`diagnose()`](https://mrc-ide.github.io/stareal/reference/diagnose.md)
to judge whether a run is long enough.

[`fit_areal()`](https://mrc-ide.github.io/stareal/reference/fit_areal.md)
also accepts a plain list of any of these arguments, e.g.
`control = list(n_iter = 500)`; the others take their defaults.

## Examples

``` r
mcmc_control(n_iter = 500, n_burn = 200, n_store = 100, n_chains = 2)
#> $n_iter
#> [1] 500
#> 
#> $n_burn
#> [1] 200
#> 
#> $n_store
#> [1] 100
#> 
#> $n_chains
#> [1] 2
#> 
#> $ess_sweeps
#> [1] 10
#> 
#> $theta_moves
#> [1] 2
#> 
#> $theta_thin
#> [1] 1
#> 
#> $seed
#> [1] 1
#> 
#> $n_cores
#> [1] 2
#> 
#> $blas_threads
#> [1] 1
#> 
#> $checkpoint_every
#> [1] 500
#> 
#> $verbose
#> [1] TRUE
#> 
```
