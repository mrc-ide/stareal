# Attach survey data to a frame

Assigns each survey to an areal unit and a time step, then pools surveys
that fall in the same unit and step: positives add, and so do the
numbers tested. The pooled cells are the only data the model sees.
Surveys that cannot be used are kept, with the reason, so that every
exclusion can be accounted for: print the result to see them.

## Usage

``` r
areal_data(
  data,
  frame,
  time,
  y,
  k,
  area = NULL,
  coords = NULL,
  crs = 4326,
  date = NULL,
  t = NULL,
  snap_max_m = 10000,
  verbose = TRUE
)
```

## Arguments

- data:

  A data frame, one row per survey.

- frame:

  A
  [`areal_frame()`](https://mrc-ide.github.io/stareal/reference/areal_frame.md).

- time:

  A
  [`time_frame()`](https://mrc-ide.github.io/stareal/reference/time_frame.md).

- y, k:

  Names of the columns holding the number of positives and the number
  tested.

- area:

  Name of a column of unit ids matching `frame$ids`. Give this or
  `coords`.

- coords:

  Names of the x and y (longitude and latitude) columns. Needs a frame
  built from polygons.

- crs:

  Coordinate reference system of `coords` (default WGS84
  longitude/latitude).

- date:

  Name of a column of dates (`Date`, date-times, or `"YYYY-MM-DD"`
  text). Give this or `t`.

- t:

  Name of a column of time steps: step numbers (1..n_steps), or the time
  frame's labels (e.g. `"2020-03"`). Required for a plain index time
  frame.

- snap_max_m:

  Maximum distance (metres) for snapping points outside every polygon to
  the nearest unit. Use 0 for no snapping.

- verbose:

  Report, in one line, how many surveys were used and how many cells
  they make.

## Value

An object of class `stareal_data`:

- cells:

  one row per observed cell, ordered by `cell_id`
  (`(area - 1) * n_steps + t`): `cell_id`, `area` (row of the frame),
  `id`, `t`, `time` (label), `y`, `k`, `n_surveys`

- surveys:

  `data` with `.id` (the area id, `NA` if none), `.t`, `.snap_m` and
  `.excluded` (`NA` if used, otherwise the reason) added

- frame, time:

  the frame and time frame

- truth:

  `NULL` (see
  [`simulate_areal()`](https://mrc-ide.github.io/stareal/reference/simulate_areal.md))

Print it for the exclusions;
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) shows where and
when the data are.

## Details

Reasons for exclusion, checked in this order:

- missing or invalid counts: `y` or `k` missing or unreadable, not whole
  numbers, negative, `k = 0`, or `y > k`. These come with a warning,
  because they are usually data errors.

- missing area, or missing or unreadable coordinates;

- area id not in the frame. This also comes with a warning.

- outside the frame: coordinates more than `snap_max_m` from every unit;

- missing date or time step;

- unreadable date (not `YYYY-MM-DD`), or a time step that is not a step
  number or label of the time frame. These come with a warning.

- outside the time window.

Points that fall outside every polygon but within `snap_max_m` of one
are assigned to the nearest unit. Survey coordinates are often rounded,
and coastal sites fall just offshore. Snap distances are recorded in
`$surveys$.snap_m`. A point inside two overlapping units goes to the
first, with a warning.

## Examples

``` r
fr <- areal_lattice(3, 3)
tf <- time_frame(n_steps = 4)
surveys <- data.frame(unit = c("r1c1", "r1c1", "r2c3", "r9c9"),
                      step = c(1, 1, 4, 2), pos = c(3, 1, 0, 2), n = c(10, 5, 8, 9))
d <- areal_data(surveys, fr, tf, y = "pos", k = "n", area = "unit", t = "step")
#> Warning: 1 survey with area ids not in the frame
#> 3 of 4 surveys used, pooled into 2 cells (print the result for the exclusions)
d
#> <stareal_data> 3 of 4 surveys used, pooled into 2 cells
#>   2 of 9 units and 2 of 4 time steps have data (5.56% of the 36 area-time cells)
#>   tested per cell: median 11.5, range 8-15; total 23 tested, 4 positive
#>   excluded:
#>          1  area not in frame
d$cells
#>    cell_id area   id t time y  k n_surveys
#> 1        1    1 r1c1 1    1 4 15         2
#> 24      24    6 r2c3 4    4 0  8         1
```
