# Trace plots of a fit

One panel per parameter in
[summary()](https://mrc-ide.github.io/stareal/reference/summary.stareal_fit.md),
one colour per chain. The hyperparameters are drawn from the full trace
on the natural scale, burn-in included (left of the dashed line). The
fixed effects are drawn from the stored draws, which start after
burn-in.

## Usage

``` r
# S3 method for class 'stareal_fit'
plot(x, ...)
```

## Arguments

- x:

  A `stareal_fit`.

- ...:

  Passed to
  [`graphics::matplot()`](https://rdrr.io/r/graphics/matplot.html) (e.g.
  `col`, `lty`, `lwd`).

## Value

`x`, invisibly.
