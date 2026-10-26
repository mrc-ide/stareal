# Map dates to time steps

Map dates to time steps

## Usage

``` r
time_index(time, dates)
```

## Arguments

- time:

  A `stareal_time` with calendar steps.

- dates:

  Dates: `Date`, date-times (read in their own time zone), or strings
  `"YYYY-MM-DD"`.

## Value

Integer step indices (1..n_steps), `NA` for missing or unreadable dates
and for dates outside the window.

## Examples

``` r
tf <- time_frame("2020-01", "2020-12")
time_index(tf, as.Date(c("2020-01-15", "2020-12-31", "2021-01-01")))
#> [1]  1 12 NA
```
