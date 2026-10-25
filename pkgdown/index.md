# stareal

**Space-time maps of prevalence from survey data, on any set of areas and
any time grid.**

stareal fits a Bayesian space-time model to binomial survey data, such as the
prevalence of a drug-resistance marker, on administrative units or grid
cells over months, quarters or years. It estimates every area in every time
step, surveyed or not, with uncertainty, and aggregates the estimates to
regions or countries draw by draw.

![Predicted prevalence of the chloroquine-resistance marker crt K76T across the eight IGAD countries, from the first tutorial.](man/figures/hero-crt76t.png)

- **Any polygons, any time step.** Your own boundaries (or an adjacency
  matrix), months, quarters or years, covariates by area, by time or both.
- **Fitted on the observed cells only.** The unsurveyed area-time cells are
  integrated out exactly, so the cost depends on the number of surveyed
  cells, not on the size of the map: continental maps over decades are
  feasible.
- **Honest group summaries.** Regional and national estimates are taken from
  joint posterior draws, so their uncertainty accounts for the correlation
  between areas.

## Installation

```r
# install.packages("remotes")
remotes::install_github("mrc-ide/stareal")
```

The sampler is almost entirely dense linear algebra: use an optimised BLAS
(Intel MKL, OpenBLAS or Apple Accelerate).

## A first map

```r
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

## Tutorials

| | |
|---|---|
| [**1. Mapping crt K76T in the Greater Horn of Africa**](articles/crt76t.html) | From a STAVE data object to maps, country trends and probability statements. |
| [**2. Rare and emerging markers: k13 R622I**](articles/k13-622i.html) | Mostly-zero data, exceedance probabilities, and when to doubt the map. |
| [**3. Covariates**](articles/covariates.html) | Covariates by area, by time and by both: building the tables and reading the coefficients. |
| [**4. Cross-validation and model comparison**](articles/crossval.html) | Which model predicts held-out surveys best, and is it calibrated? |
| [**5. Simulate, fit, check**](articles/simulation.html) | Data simulated from known parameters: does the fit recover the truth? |

[**Get started**](articles/stareal.html) runs through every step on simulated
data in a few minutes. The [**Method**](articles/method.html) page explains
the model and how it is fitted.

## Citation

Verity, R. stareal: space-time areal models fitted on the observed cells
only. R package, <https://github.com/mrc-ide/stareal>.
