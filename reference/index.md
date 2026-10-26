# Package index

## Space and time

- [`areal_frame()`](https://mrc-ide.github.io/stareal/reference/areal_frame.md)
  : Build a spatial frame
- [`areal_lattice()`](https://mrc-ide.github.io/stareal/reference/areal_lattice.md)
  : A regular lattice frame
- [`plot(`*`<stareal_frame>`*`)`](https://mrc-ide.github.io/stareal/reference/plot.stareal_frame.md)
  : Map a frame and its graph
- [`time_frame()`](https://mrc-ide.github.io/stareal/reference/time_frame.md)
  : Build a time frame
- [`time_index()`](https://mrc-ide.github.io/stareal/reference/time_index.md)
  : Map dates to time steps

## Data

- [`areal_data()`](https://mrc-ide.github.io/stareal/reference/areal_data.md)
  : Attach survey data to a frame
- [`plot(`*`<stareal_data>`*`)`](https://mrc-ide.github.io/stareal/reference/plot.stareal_data.md)
  : Plot survey data
- [`simulate_areal()`](https://mrc-ide.github.io/stareal/reference/simulate_areal.md)
  : Simulate survey data from the model

## Fitting

- [`fit_areal()`](https://mrc-ide.github.io/stareal/reference/fit_areal.md)
  : Fit the space-time areal model
- [`trend()`](https://mrc-ide.github.io/stareal/reference/trend.md) : A
  linear trend in time
- [`areal_priors()`](https://mrc-ide.github.io/stareal/reference/areal_priors.md)
  : Priors for an areal model
- [`mcmc_control()`](https://mrc-ide.github.io/stareal/reference/mcmc_control.md)
  : Settings for the MCMC sampler
- [`summary(`*`<stareal_fit>`*`)`](https://mrc-ide.github.io/stareal/reference/summary.stareal_fit.md)
  : Summarise a fit
- [`diagnose()`](https://mrc-ide.github.io/stareal/reference/diagnose.md)
  : Check convergence and mixing of a fit
- [`plot(`*`<stareal_fit>`*`)`](https://mrc-ide.github.io/stareal/reference/plot.stareal_fit.md)
  : Trace plots of a fit
- [`as_draws(`*`<stareal_fit>`*`)`](https://mrc-ide.github.io/stareal/reference/as_draws.stareal_fit.md)
  [`as_draws_array(`*`<stareal_fit>`*`)`](https://mrc-ide.github.io/stareal/reference/as_draws.stareal_fit.md)
  [`as_draws_df(`*`<stareal_fit>`*`)`](https://mrc-ide.github.io/stareal/reference/as_draws.stareal_fit.md)
  [`as_draws_matrix(`*`<stareal_fit>`*`)`](https://mrc-ide.github.io/stareal/reference/as_draws.stareal_fit.md)
  [`as_draws_list(`*`<stareal_fit>`*`)`](https://mrc-ide.github.io/stareal/reference/as_draws.stareal_fit.md)
  [`as_draws_rvars(`*`<stareal_fit>`*`)`](https://mrc-ide.github.io/stareal/reference/as_draws.stareal_fit.md)
  : Posterior draws of a fit

## Prediction and checking

- [`predict(`*`<stareal_fit>`*`)`](https://mrc-ide.github.io/stareal/reference/predict.stareal_fit.md)
  : Predict the full area-time surface from a fit
- [`plot(`*`<stareal_pred>`*`)`](https://mrc-ide.github.io/stareal/reference/plot.stareal_pred.md)
  : Plot predictions
- [`as.data.frame(`*`<stareal_pred>`*`)`](https://mrc-ide.github.io/stareal/reference/as.data.frame.stareal_pred.md)
  : Predictions as a data frame
- [`st_as_sf(`*`<stareal_pred>`*`)`](https://mrc-ide.github.io/stareal/reference/st_as_sf.stareal_pred.md)
  : Predictions as polygons
- [`crossval()`](https://mrc-ide.github.io/stareal/reference/crossval.md)
  : Cross-validate a fit
- [`compare_cv()`](https://mrc-ide.github.io/stareal/reference/compare_cv.md)
  : Compare cross-validated models

## Example data

Drug-resistance surveys in the eight IGAD countries.

- [`igad_crt76t`](https://mrc-ide.github.io/stareal/reference/igad_crt76t.md)
  : Surveys of the chloroquine-resistance marker crt K76T in IGAD
  countries

## Package

- [`stareal`](https://mrc-ide.github.io/stareal/reference/stareal-package.md)
  [`stareal-package`](https://mrc-ide.github.io/stareal/reference/stareal-package.md)
  :

  *stareal*: spatial-temporal areal estimation from survey data
