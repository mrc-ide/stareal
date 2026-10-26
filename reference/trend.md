# A linear trend in time

For use in the `formula` of
[`fit_areal()`](https://mrc-ide.github.io/stareal/reference/fit_areal.md)
and
[`simulate_areal()`](https://mrc-ide.github.io/stareal/reference/simulate_areal.md)
only: a linear effect of time, running from -1 at the first time step to
1 at the last, so its coefficient is half the change in logit prevalence
over the window. With one time step it is dropped.

## Usage

``` r
trend()
```

## Value

Not meant to be called.
