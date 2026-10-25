# stareal (development version)

## 2026-10-25: example data and the website

- **Example data.** `igad_adm1` holds the 268 administrative units of the
  eight IGAD countries, from geoBoundaries (Uganda by district).
  `igad_crt76t` (272 surveys) and `igad_k13_622i` (423 surveys) hold two
  drug-resistance markers from the public STAVE collation.
- **Borders that nearly touch.** `areal_frame(adjacency_snap = )` makes
  units whose borders come within that many metres neighbours. This is
  for boundary files assembled from different sources, which leave small
  gaps along national borders.
- **Folds of your own.** `crossval(folds = )` takes a fold label per
  observed cell, e.g. a country, for spatial cross-validation.
- **Comparing models.** `compare_cv()` compares models cross-validated
  on the same folds, cell by cell, with the standard error of the
  difference.
- **Simulating at a real design.** `simulate_areal(cells = )` surveys
  the given area-time cells, with their numbers tested, so a simulation
  can copy a real survey design.
- **Fix: cross-validation log scores.** The quadrature behind
  `crossval()`'s log score found the integrand's mode by plain Newton
  steps. When a held-out survey was far from its predictive mean, these
  overshot, and the log score came out wrong by thousands (e.g. -18,698
  for a survey inside its 95% interval). The steps are now damped. PIT,
  intervals and RMSE were not affected.
- **Website.** Five tutorials (crt K76T, k13 R622I, covariates,
  cross-validation, simulation), with their long fits precomputed by the
  scripts in `vignettes/articles/_build/`.

## 2026-10-25: review and restructure

### Interface

- **Simulations.** `simulate_areal()` returns data ready to fit (a
  `stareal_data` with `$truth`), so `fit_areal(sim)` works directly.
- **Fit arguments.** `fit_areal(data, formula = ~ trend() + x, ...)`:
  `formula` replaces `mean`, and `checkpoint_dir` is now an argument of
  `fit_areal()` and `crossval()`. `priors` and `control` accept plain
  lists of just the settings to change, and fixed-effect priors can be
  named by design column.
- **Settings.** `areal_priors()` and `mcmc_control()` return plain
  validated lists. In `mcmc_control()`, `n_cores = 1` replaces
  `parallel = FALSE`, and `verbose` replaces `progress`.
- **Summaries.** `summary()` and `diagnose()` return one class.
  `diagnose()` adds every observed cell, a check of the fixed effects'
  ESS, and the verdict.
- **Prediction.** `predict(fit, groups = ...)` replaces `levels`, and
  group tables are now under `pr$groups`. New: `as.data.frame(pr)` and
  `sf::st_as_sf(pr, time = ...)`.
- **Draws.** All `posterior::as_draws_*()` conversions work on a fit.
- **Time.** `time_frame(n_steps = )` replaces `T`, and its field is
  `$n_steps`. `trend()` runs from -1 at the first time step to 1 at the
  last, and is dropped when there is one step.
- `areal_data()` no longer prints on assignment, and its surveys carry
  `.id` (the area id).

### Sampler

- The pseudo-data that steer the latent and hyperparameter moves are
  refitted at the chain's mean during burn-in, then frozen. For rare
  outcomes this mixes the latent cells about 8 times better.
- Chains start from dispersed points, so R-hat can reveal chains that
  have not mixed.
- Seeds:
  - each chain and each prediction draw gets its own seed, drawn from
    the user's, so any seed works;
  - results do not depend on the session's `RNGkind()`.

### Running

- **Live progress.** Worker processes' progress now appears in RStudio,
  not only in a terminal.
- **Interrupts.** On interrupt or error, the workers are stopped instead
  of running on.
- **Checkpoint folders.** A folder is tied to its data, graph, model,
  settings and sampler version, but not to the R version.

### Cross-validation

- The log score integrates the held-out Gaussian exactly.
- "Outside the 95% interval" uses a randomised PIT, so a calibrated model
  shows about 5% whatever the number tested.

### Fixes

About 80 smaller fixes: input validation (factors, date-times, numeric
ids, non-finite values, name clashes), memory bounds in prediction,
messages, and documentation. See `review_stareal-v0.1_2026-10-24.md` in
the source repository.
