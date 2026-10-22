
<!-- README.md is generated from README.Rmd. Please edit that file -->

# stareal

Space–time areal models for binomial survey data, fitted on the observed
cells only.

**Status: in development (v0.0.0.9000).**

## The idea

Survey counts (positives / tested) are pooled into area × time cells and
modelled as binomial. The logit prevalence is fixed effects (intercept,
trend, covariates) plus a separable Gaussian field (Leroux in space, AR1
in time) plus a cell-level nugget. Almost all area × time cells have no
data. The prior is Gaussian, so those cells can be integrated out
exactly, and the model is fitted on the n observed cells only instead of
on all m × T cells. The full table is then reconstructed by conditional
simulation (Matheron’s rule), one surface per posterior draw. Coarser
groupings (e.g. regions, countries) are aggregated draw by draw, so
their uncertainty accounts for the correlation between areas.

The method was developed for continental maps of antimalarial
drug-resistance markers in Africa: about 4,300 districts by 360 months,
fitted on about 1,800 observed cells.

## Installation

``` r
# install.packages("remotes")
remotes::install_github("mrc-ide/stareal")
```

Use an optimised BLAS/LAPACK (Intel MKL, OpenBLAS or Apple Accelerate).
The sampler is almost entirely dense linear algebra.

## A first run

``` r
library(stareal)
fr <- areal_frame(polygons, id = "district_id", groups = c(country = "iso3"))  # an sf object
tf <- time_frame(2000, 2024, "year")
d <- areal_data(surveys, fr, tf, y = "positive", k = "tested",
                coords = c("lon", "lat"), date = "date")
fit <- fit_areal(d, control = mcmc_control(checkpoint_dir = "chains"))
diagnose(fit)                                   # every parameter, latent cells included
pr <- predict(fit, levels = "country")
pr$area; pr$country                             # tidy tables of estimates
plot(pr, time = "2024")
```

`vignette("stareal")` runs a complete analysis on simulated data, and
`demo("walkthrough", package = "stareal")` runs every function.

## Licence

MIT.
