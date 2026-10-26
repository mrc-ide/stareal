# Predict the full area-time surface from a fit

Reconstructs the latent prevalence surface for every area and time step,
in percent (as all prevalences in *stareal* are), observed or not, once
per posterior draw. It does this by conditional simulation (Matheron's
rule), which is exact given each draw. It then summarises the draws per
area and time step. Each grouping in `groups` is aggregated draw by
draw, as a weighted mean of its areas' prevalence, and then summarised.
Its uncertainty therefore accounts for the correlation between areas.

## Usage

``` r
# S3 method for class 'stareal_fit'
predict(
  object,
  groups = NULL,
  weights = NULL,
  n_draws = NULL,
  probs = c(0.025, 0.5, 0.975),
  exceedance = NULL,
  keep_draws = NULL,
  seed = 1,
  n_cores = NULL,
  blas_threads = 1,
  block_bytes = 5e+08,
  verbose = TRUE,
  ...
)
```

## Arguments

- object:

  A
  [`fit_areal()`](https://mrc-ide.github.io/stareal/reference/fit_areal.md)
  object.

- groups:

  Names of groupings of the frame (see
  [`areal_frame()`](https://mrc-ide.github.io/stareal/reference/areal_frame.md))
  to aggregate to. Areas are always summarised.

- weights:

  Weights of the areas in group aggregates, for example population or
  number of infections. The aggregate of a group is \\\sum_a w\_{at}
  p\_{at} / \sum_a w\_{at}\\ for each draw. Give one of: `NULL` (equal
  weights); a vector with one value per area (named by area id, or in
  frame order); an areas x time steps matrix (rows named by area id or
  in frame order; columns named by time label or in order); or a data
  frame with columns `id`, `weight` and optionally `time` (label) or `t`
  (step number), with every area-time combination present once.

- n_draws:

  Number of posterior draws to use, taken evenly from the stored draws
  of all chains. `NULL`: all of them.

- probs:

  Quantiles to report.

- exceedance:

  Prevalence thresholds, in percent (e.g. `5` for 5%); for each, the
  posterior probability that prevalence exceeds it is reported, in a
  column `p_gt_<threshold>` (e.g. `p_gt_5`).

- keep_draws:

  Ids of areas whose prevalence draws are returned in full.

- seed:

  Random seed, a whole number. Each draw gets its own seed, drawn from
  this one. The caller's random-number stream is left untouched.

- n_cores:

  Number of worker processes. Default: 1 (this session) if the whole
  surface fits in one block of `block_bytes`, otherwise the cores less
  one. Every worker holds its own copy of the frame's eigenvectors (8
  m^2 bytes).

- blas_threads:

  BLAS threads per worker.

- block_bytes:

  Approximate memory for one block's draws (bytes). Areas are processed
  in blocks of at most this size (up to `n_cores` at once, each also
  needing about half as much again per exceedance threshold). Each block
  regenerates every draw's noise, so fewer, larger blocks are faster.

- verbose:

  Report progress: the draws conditioned and predicted so far, with an
  estimate of the time left.

- ...:

  Not used.

## Value

An object of class `stareal_pred`:

- area:

  one row per area and time step, in the order of area, then time step:
  `id`, the frame's groupings, `t`, `time`, `date` (start of the step;
  calendar time frames), the pooled data `y` and `k` (0 where there are
  none), then `mean`, `sd` and one column per quantile (`q2.5`, `q50`,
  ...) of prevalence, in percent (0 to 100), and one column per
  threshold (`p_gt_5`, ...), the probability (0 to 1) that prevalence
  exceeds it

- groups:

  a list with one table per grouping in `groups`: `id` (the group), `t`,
  `time`, `date` and the same statistics, ordered by group, then time
  step

- draws:

  if `keep_draws`: prevalence draws, in percent, an array (area x time
  step x draw)

- frame, time, settings:

  for plotting and the record (the frame without its eigenvectors)

[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) gives the
area table (or a group's, with `group =`);
`sf::st_as_sf(pr, time = ...)` gives it as polygons for one time step.

## Details

The surface is the field plus fixed effects, without the cell-level
nugget: it is the prevalence of the area-time cell, not of one survey in
it. At observed cells it therefore differs from the latent values that
[`diagnose()`](https://mrc-ide.github.io/stareal/reference/diagnose.md)
and
[`crossval()`](https://mrc-ide.github.io/stareal/reference/crossval.md)
report, which include the nugget.

## Examples

``` r
fr <- areal_lattice(4, 4, blocks = c(2, 2))
tf <- time_frame(n_steps = 5)
sim <- simulate_areal(fr, tf, n_cells = 30, seed = 1)
fit <- fit_areal(sim, control = list(n_iter = 100, n_burn = 50, n_store = 25, n_chains = 1,
                                     ess_sweeps = 2, verbose = FALSE))
pr <- predict(fit, groups = "block", verbose = FALSE)
head(pr$area)
#>     id block t time  y  k     mean        sd      q2.5      q50    q97.5
#> 1 r1c1  b1_1 1    1 15 50 27.76358  5.244571 17.475060 27.50821 36.50347
#> 2 r1c1  b1_1 2    2  0  0 38.19425  9.548869 23.791310 38.46915 57.77599
#> 3 r1c1  b1_1 3    3 24 50 51.40362  9.589957 34.541165 51.02984 68.07607
#> 4 r1c1  b1_1 4    4  0  0 57.20519 14.301829 28.487682 60.11959 79.10870
#> 5 r1c1  b1_1 5    5  0  0 62.75272 15.872400 32.238574 64.81357 84.13614
#> 6 r1c2  b1_1 1    1  0  0 19.58383  6.093364  9.670767 18.63618 31.13343
pr$groups$block
#>      id t time     mean       sd     q2.5      q50    q97.5
#> 1  b1_1 1    1 24.03205 3.339114 17.35945 23.82604 29.27107
#> 2  b1_1 2    2 35.36298 5.107788 28.11229 35.39810 44.00332
#> 3  b1_1 3    3 48.02901 5.219416 38.02603 48.01526 55.32425
#> 4  b1_1 4    4 54.69573 7.097510 41.77794 54.45067 65.19264
#> 5  b1_1 5    5 61.53742 7.767997 46.79119 61.43908 73.61319
#> 6  b1_2 1    1 18.90888 4.012992 10.62747 18.21439 25.70471
#> 7  b1_2 2    2 28.94552 4.544896 21.08419 29.11678 36.34366
#> 8  b1_2 3    3 40.27671 4.483327 31.00523 40.63563 47.14644
#> 9  b1_2 4    4 49.06061 5.120647 38.14857 49.90036 57.80228
#> 10 b1_2 5    5 58.39946 8.131232 41.76937 59.84926 69.39045
#> 11 b2_1 1    1 13.39723 3.590305  7.06265 13.56779 19.36011
#> 12 b2_1 2    2 20.73469 5.461339 12.25727 20.16163 30.96478
#> 13 b2_1 3    3 31.83823 6.241591 20.92492 32.29420 42.53589
#> 14 b2_1 4    4 39.93082 5.033122 30.40158 40.06794 46.27097
#> 15 b2_1 5    5 51.20371 7.962376 33.12924 53.09598 63.34420
#> 16 b2_2 1    1 20.87254 5.869345 12.44359 20.81967 31.38286
#> 17 b2_2 2    2 32.90281 6.357740 22.19187 31.84416 44.45471
#> 18 b2_2 3    3 42.92875 4.035860 36.92262 41.96534 51.29063
#> 19 b2_2 4    4 50.30957 5.148825 42.72004 50.06091 59.25782
#> 20 b2_2 5    5 57.20326 8.332660 41.83662 57.89208 71.48372
```
