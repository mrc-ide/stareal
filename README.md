
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
drug-resistance markers in Africa: thousands of districts by hundreds of
months, fitted on a couple of thousand observed cells.

## Installation

``` r
# install.packages("remotes")
remotes::install_github("mrc-ide/stareal")
```

Use an optimised BLAS/LAPACK (Intel MKL, OpenBLAS or Apple Accelerate).
The sampler is almost entirely dense linear algebra.

## A first map

``` r
library(stareal)
fr <- areal_frame(igad_adm1, id = "id", groups = c(country = "country"),
                  link_within = "country", adjacency_snap = 2000) # areas and their neighbours
tf <- time_frame(1989, 2024, "year")                             # the time grid
d  <- areal_data(igad_crt76t, fr, tf, y = "positive", k = "tested",
                 coords = c("longitude", "latitude"), date = "collection_day")
fit <- fit_areal(d)                                              # MCMC, chains in parallel
diagnose(fit)                                                    # convergence of every parameter
pr <- predict(fit, groups = "country")                           # every area-year, and countries
plot(pr, time = "2024")
```

The data are the chloroquine-resistance marker crt K76T in the eight
IGAD countries, from the public
[STAVE](https://github.com/mrc-ide/STAVE) collation. The [package
website](https://mrc-ide.github.io/stareal/) has tutorials, from this
analysis to rare markers, covariates, cross-validation and simulation
studies.

## Licence

MIT.
