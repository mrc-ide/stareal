# Changelog

## stareal 0.1.0

The first numbered version: everything below, plus the 2026-10-28
changes.

### 2026-10-28: spatial priors, model comparison, group draws

- **More spatial priors.** `fit_areal(spatial = )` chooses the prior on
  the graph: `"leroux"` (the default, as before), `"bym2"` (Riebler et
  al. 2016: a scaled ICAR mixed with independent noise, with `sigma`
  roughly the marginal sd and `rho_s` the structured share of the
  variance), `"car"` (the proper CAR), and the shortcuts `"iid"` (Leroux
  with `rho_s = 0`) and `"icar"` (Leroux with `rho_s = 1`). All use the
  observed-cell machinery unchanged, so they cost the same to fit and
  predict. `simulate_areal(spatial = )` draws from them too.
- **Fixed hyperparameters.** `fit_areal(fixed = list(rho_t = 0))` holds
  any of `sigma`, `rho_s`, `rho_t` and `sigma_n` at a value instead of
  sampling it, e.g. a field with no memory between time steps. Fixed
  values are left out of the draws,
  [`summary()`](https://rdrr.io/r/base/summary.html) and
  [`diagnose()`](https://mrc-ide.github.io/stareal/reference/diagnose.md).
- **[`loo()`](https://mc-stan.org/loo/reference/loo.html)**: approximate
  leave-one-out cross-validation from one fit, with no refits. Each
  observed cell’s own latent value is integrated out before
  Pareto-smoothed importance sampling (Vehtari et al. 2016), which makes
  it stable where plain PSIS-LOO on this model is not. It returns a
  `loo` object, so
  [`loo_compare()`](https://mc-stan.org/loo/reference/loo_compare.html)
  works on it.
- **[`compare_models()`](https://mrc-ide.github.io/stareal/reference/compare_models.md)**:
  one row per fit of the same data, best first: the model, whether it
  passes
  [`diagnose()`](https://mrc-ide.github.io/stareal/reference/diagnose.md),
  the leave-one-out score with the difference from the best and its
  standard error, and, if given,
  [`crossval()`](https://mrc-ide.github.io/stareal/reference/crossval.md)
  scores. Its help page says why it does not use AIC, BIC or DIC.
- **Draws of aggregates.** `predict(group_draws = )` returns the draws
  of chosen groupings’ aggregates (group x time step x draw), the draws
  the summaries are made from, so changes between years, slopes over a
  window and contrasts can be computed draw by draw. `keep_times =`
  keeps those draws, and `keep_draws`, at chosen time steps only.
- **Groupings given at prediction time.**
  `predict(groups = list("country", zone = labels))` aggregates to a
  grouping that is not in the frame, without rebuilding the frame (or
  refitting).
- **What made a result.** Fits, predictions, cross-validations and
  [`loo()`](https://mc-stan.org/loo/reference/loo.html) results carry
  `$version`: the *stareal* version, its commit when installed from
  GitHub, the sampler version, R and the BLAS.
- The default model’s draws are unchanged, bit for bit, and its
  checkpoint folders still resume.
- **Website.** A new tutorial 4, “Choosing a model”: the five spatial
  priors, switching the nugget and the trend off, a fixed `rho_t`, and
  [`loo()`](https://mc-stan.org/loo/reference/loo.html) and
  [`compare_models()`](https://mrc-ide.github.io/stareal/reference/compare_models.md),
  checked against tutorial 5’s cross-validation. Tutorial 1 gains
  changes between years and a grouping made at prediction time, from
  `group_draws`. The simulation tutorial (now 6) asks whether
  [`compare_models()`](https://mrc-ide.github.io/stareal/reference/compare_models.md)
  finds the prior that made the data. The Method page covers every prior
  and the leave-one-out. The cross-validation tutorial is now 5.
- **Checked.** Simulation-based calibration (`data-raw/check_sbc.R`, now
  for any prior) gives uniform ranks for every parameter under BYM2, the
  proper CAR, ICAR and a fixed `rho_t`. On simulated data,
  [`loo()`](https://mc-stan.org/loo/reference/loo.html) agrees with
  exact leave-one-out (40 refits) to within 0.5 in total elpd, with
  cell-by-cell correlation above 0.999, in seconds instead of minutes
  (`data-raw/check_loo.R`).

### 2026-10-26: website

- The simulated-data “Getting started” vignette is gone: the tutorials
  cover all of it with real data. “Get started” on the website now opens
  tutorial 1, and the notes on how cost scales moved to the Method page.
  For an offline tour, run `demo("walkthrough", package = "stareal")`.

### 2026-10-26: percent, panel labels, progress

- **Prevalence is in percent (0 to 100) everywhere**:
  [`predict()`](https://rdrr.io/r/stats/predict.html)’s summaries and
  kept draws, [`plot()`](https://rdrr.io/r/graphics/plot.default.html),
  [`crossval()`](https://mrc-ide.github.io/stareal/reference/crossval.md)’s
  observed and predicted prevalences (and so its RMSE, in percentage
  points),
  [`simulate_areal()`](https://mrc-ide.github.io/stareal/reference/simulate_areal.md)’s
  true prevalence, and `summary(fit, latent = TRUE)`. Exceedance
  thresholds are given in percent too: `exceedance = 5` reports
  `p_gt_5`, the probability (0 to 1) that prevalence exceeds 5%.
- **Panel labels.** `plot(pr, ids = , labels = )` titles time-series
  panels with your own labels, one per id or named by id.
- **Progress from [`predict()`](https://rdrr.io/r/stats/predict.html).**
  It reports the draws conditioned and predicted so far, with an
  estimate of the time left, in both phases and from parallel workers
  too.

### 2026-10-25: example data and the website

- **Example data.** `igad_crt76t` holds 328 surveys of the
  chloroquine-resistance marker crt K76T in the eight IGAD countries,
  from the public STAVE collation. The surveys were selected by
  location, not by STAVE’s country names, which are unreliable. No
  boundaries ship with the package: the examples and tutorials download
  GADM with
  [`geodata::gadm()`](https://rspatial.github.io/geodata/reference/gadm.html).
- **Plots.** `plot(fit)` shows every parameter of
  [`summary()`](https://rdrr.io/r/base/summary.html), fixed effects
  included. Time-series plots of predictions (`plot(pr)`, by area or
  group) add the observed prevalence of each surveyed area and time
  step, with its exact 95% binomial interval (`data = TRUE`).
- **Title.** “Spatial-Temporal Areal Estimation from Survey Data”.
- **Borders that nearly touch.** `areal_frame(adjacency_snap = )` makes
  units whose borders come within that many metres neighbours. This is
  for boundary files assembled from different sources, which leave small
  gaps along national borders.
- **Folds of your own.** `crossval(folds = )` takes a fold label per
  observed cell, e.g. a country, for spatial cross-validation.
- **Comparing models.**
  [`compare_cv()`](https://mrc-ide.github.io/stareal/reference/compare_cv.md)
  compares models cross-validated on the same folds, cell by cell, with
  the standard error of the difference.
- **Simulating at a real design.** `simulate_areal(cells = )` surveys
  the given area-time cells, with their numbers tested, so a simulation
  can copy a real survey design.
- **Fix: cross-validation log scores.** The quadrature behind
  [`crossval()`](https://mrc-ide.github.io/stareal/reference/crossval.md)’s
  log score found the integrand’s mode by plain Newton steps. When a
  held-out survey was far from its predictive mean, these overshot, and
  the log score came out wrong by thousands (e.g. -18,698 for a survey
  inside its 95% interval). The steps are now damped. PIT, intervals and
  RMSE were not affected.
- **Website.** Five tutorials (crt K76T by month, k13 R561H in East
  Africa, covariates, cross-validation, simulation), each running in
  minutes. The two long computations, the cross-validation of five
  models and a 20-replicate simulation study, ship with the package in
  `inst/extdata/` and are imported by their tutorials.

### 2026-10-25: review and restructure

#### Interface

- **Simulations.**
  [`simulate_areal()`](https://mrc-ide.github.io/stareal/reference/simulate_areal.md)
  returns data ready to fit (a `stareal_data` with `$truth`), so
  `fit_areal(sim)` works directly.
- **Fit arguments.** `fit_areal(data, formula = ~ trend() + x, ...)`:
  `formula` replaces `mean`, and `checkpoint_dir` is now an argument of
  [`fit_areal()`](https://mrc-ide.github.io/stareal/reference/fit_areal.md)
  and
  [`crossval()`](https://mrc-ide.github.io/stareal/reference/crossval.md).
  `priors` and `control` accept plain lists of just the settings to
  change, and fixed-effect priors can be named by design column.
- **Settings.**
  [`areal_priors()`](https://mrc-ide.github.io/stareal/reference/areal_priors.md)
  and
  [`mcmc_control()`](https://mrc-ide.github.io/stareal/reference/mcmc_control.md)
  return plain validated lists. In
  [`mcmc_control()`](https://mrc-ide.github.io/stareal/reference/mcmc_control.md),
  `n_cores = 1` replaces `parallel = FALSE`, and `verbose` replaces
  `progress`.
- **Summaries.** [`summary()`](https://rdrr.io/r/base/summary.html) and
  [`diagnose()`](https://mrc-ide.github.io/stareal/reference/diagnose.md)
  return one class.
  [`diagnose()`](https://mrc-ide.github.io/stareal/reference/diagnose.md)
  adds every observed cell, a check of the fixed effects’ ESS, and the
  verdict.
- **Prediction.** `predict(fit, groups = ...)` replaces `levels`, and
  group tables are now under `pr$groups`. New: `as.data.frame(pr)` and
  `sf::st_as_sf(pr, time = ...)`.
- **Draws.** All `posterior::as_draws_*()` conversions work on a fit.
- **Time.** `time_frame(n_steps = )` replaces `T`, and its field is
  `$n_steps`.
  [`trend()`](https://mrc-ide.github.io/stareal/reference/trend.md) runs
  from -1 at the first time step to 1 at the last, and is dropped when
  there is one step.
- [`areal_data()`](https://mrc-ide.github.io/stareal/reference/areal_data.md)
  no longer prints on assignment, and its surveys carry `.id` (the area
  id).

#### Sampler

- The pseudo-data that steer the latent and hyperparameter moves are
  refitted at the chain’s mean during burn-in, then frozen. For rare
  outcomes this mixes the latent cells about 8 times better.
- Chains start from dispersed points, so R-hat can reveal chains that
  have not mixed.
- Seeds:
  - each chain and each prediction draw gets its own seed, drawn from
    the user’s, so any seed works;
  - results do not depend on the session’s
    [`RNGkind()`](https://rdrr.io/r/base/Random.html).

#### Running

- **Live progress.** Worker processes’ progress now appears in RStudio,
  not only in a terminal.
- **Interrupts.** On interrupt or error, the workers are stopped instead
  of running on.
- **Checkpoint folders.** A folder is tied to its data, graph, model,
  settings and sampler version, but not to the R version.

#### Cross-validation

- The log score integrates the held-out Gaussian exactly.
- “Outside the 95% interval” uses a randomised PIT, so a calibrated
  model shows about 5% whatever the number tested.

#### Fixes

About 80 smaller fixes: input validation (factors, date-times, numeric
ids, non-finite values, name clashes), memory bounds in prediction,
messages, and documentation. See `review_stareal-v0.1_2026-10-24.md` in
the source repository.
