# Build a time frame

The model's time axis is a regular grid of `n_steps` steps. A step is a
calendar month, quarter or year, running from the step containing
`start` to the step containing `end`. Alternatively, give only `n_steps`
for a plain index `1..n_steps`; survey data then carry the time step as
a whole number.

## Usage

``` r
time_frame(
  start = NULL,
  end = NULL,
  step = c("month", "quarter", "year"),
  n_steps = NULL
)
```

## Arguments

- start, end:

  The first and last dates of the window: `Date`s, or strings `"YYYY"`,
  `"YYYY-MM"` or `"YYYY-MM-DD"`, or whole years as numbers. A year or
  month stands for the whole period: as `start` its first day, as `end`
  its last (so `time_frame(2001, 2020, "month")` runs from January 2001
  to December 2020). The window then covers every step that contains a
  day between `start` and `end`.

- step:

  `"month"`, `"quarter"` or `"year"`.

- n_steps:

  Instead of `start` and `end`: the number of steps of a plain index
  frame.

## Value

An object of class `stareal_time`: `n_steps`, `step`, `labels` (one per
step, e.g. `"2020-03"`, `"2020-Q1"`, `"2020"`), and for calendar steps
`starts`, the first day of each step.

## Details

The temporal prior is a first-order autoregression between consecutive
steps, so the step sets the scale on which the temporal correlation
\\\rho_t\\ is read.

## See also

[`time_index()`](https://mrc-ide.github.io/stareal/reference/time_index.md)
to map dates to steps.

## Examples

``` r
time_frame("1995-01", "2024-12", "month")
#> <stareal_time> monthly steps, 1995-01 to 2024-12 (360 steps)
time_frame(2000, 2009, "year")
#> <stareal_time> yearly steps, 2000 to 2009 (10 steps)
time_frame(n_steps = 12)
#> <stareal_time> plain index, 12 steps
```
