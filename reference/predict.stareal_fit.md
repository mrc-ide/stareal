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
  group_draws = NULL,
  keep_times = NULL,
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

  Groupings to aggregate to (areas are always summarised). Give names of
  the frame's groupings (see
  [`areal_frame()`](https://mrc-ide.github.io/stareal/reference/areal_frame.md)),
  or a list whose elements are such names or new groupings, named, each
  a group label per area (in frame order, or named by area id). For
  example `groups = list("country", zone = zones)` aggregates to the
  frame's countries and to `zones`. A new grouping's name may not be one
  the package uses (see
  [`areal_frame()`](https://mrc-ide.github.io/stareal/reference/areal_frame.md))
  or another grouping of the frame.

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

- group_draws:

  Names of groupings in `groups` whose aggregates' draws are returned in
  full (`TRUE`: all of them). Quantities derived from an aggregate need
  these: the change between two years, the slope over a window, or a
  contrast between groups. Draws of two separate fits are independent of
  each other, so a contrast between fits treats them as independent.

- keep_times:

  Time steps (labels or step numbers) at which the draws of `keep_draws`
  and `group_draws` are kept (default: all). The summaries always cover
  every step.

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
  step x draw), at the steps of `keep_times`

- group_draws:

  if `group_draws`: a list with one array per grouping (group x time
  step x draw) of the aggregates' draws, in percent, at the steps of
  `keep_times`. Draw `d` is the same posterior draw in every area, group
  and time step, and in `draws`, so differences and contrasts can be
  taken draw by draw.

- frame, time, settings:

  for plotting and the record (the frame without its eigenvectors;
  groupings given in `groups` are added to its groupings)

- version:

  what made the prediction (see
  [`fit_areal()`](https://mrc-ide.github.io/stareal/reference/fit_areal.md))

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
#>     id block t time  y  k     mean        sd     q2.5      q50    q97.5
#> 1 r1c1  b1_1 1    1 10 50 19.37656  5.528601 10.35617 20.73413 28.16670
#> 2 r1c1  b1_1 2    2  0  0 35.43921 13.130501 10.27090 36.38516 57.78506
#> 3 r1c1  b1_1 3    3 27 50 47.93018  8.129873 33.11751 47.69561 59.96458
#> 4 r1c1  b1_1 4    4  0  0 51.18328 12.987212 23.71910 50.03455 73.21236
#> 5 r1c1  b1_1 5    5  0  0 59.61234 14.380280 35.54940 57.83576 83.86135
#> 6 r1c2  b1_1 1    1  0  0 20.65590  8.796724  8.41993 17.82740 36.45238
pr$groups$block
#>      id t time     mean        sd      q2.5      q50    q97.5
#> 1  b1_1 1    1 18.51268  5.612897 10.233253 18.44443 29.60909
#> 2  b1_1 2    2 33.81141  7.871984 23.997584 31.31148 54.03462
#> 3  b1_1 3    3 43.10237  4.541246 36.199974 41.97948 51.70674
#> 4  b1_1 4    4 48.08474  6.181675 36.158179 47.75567 58.71610
#> 5  b1_1 5    5 58.19363  6.926333 45.065057 57.78694 70.59609
#> 6  b1_2 1    1 28.89787  8.517588 12.842407 29.33691 41.81237
#> 7  b1_2 2    2 41.33818  7.534992 24.077741 43.03764 53.67862
#> 8  b1_2 3    3 54.90501  5.257541 45.837494 54.47138 64.36265
#> 9  b1_2 4    4 57.76655  7.888980 38.719592 60.85025 65.79321
#> 10 b1_2 5    5 68.85309  6.035256 57.839158 68.94146 77.70193
#> 11 b2_1 1    1 11.79531  5.299858  6.332922 10.09583 24.54373
#> 12 b2_1 2    2 21.95405  9.179031 12.489411 19.13452 42.72574
#> 13 b2_1 3    3 28.09970  8.071707 17.641259 27.38193 44.44286
#> 14 b2_1 4    4 33.46455  7.127400 24.833365 31.63830 48.25212
#> 15 b2_1 5    5 46.31201 10.765997 33.582365 43.12296 71.51349
#> 16 b2_2 1    1 23.16693  7.419873 10.771692 23.65972 37.80891
#> 17 b2_2 2    2 35.27516  7.803052 24.565714 35.04836 50.97115
#> 18 b2_2 3    3 43.05841  5.527102 32.749415 43.42534 52.98384
#> 19 b2_2 4    4 48.15162  7.589381 32.657491 48.47052 59.38774
#> 20 b2_2 5    5 60.23039  7.779010 48.263709 59.44103 75.26021

# a grouping given here, with the draws of its aggregates at two steps
halves <- ifelse(as.integer(sub("r([0-9]+)c.*", "\\1", fr$ids)) <= 2, "north", "south")
pr2 <- predict(fit, groups = list("block", half = halves), group_draws = "half", keep_times = c(1, 5),
               verbose = FALSE)
# the posterior probability that the north's prevalence fell between steps 1 and 5
mean(pr2$group_draws$half["north", "5", ] < pr2$group_draws$half["north", "1", ])
#> [1] 0
```
