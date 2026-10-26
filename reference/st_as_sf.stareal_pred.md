# Predictions as polygons

The area (or group) table for one time step joined to the frame's
polygons, ready for mapping with your own tools (e.g.
[`ggplot2::geom_sf()`](https://ggplot2.tidyverse.org/reference/ggsf.html)).
Groups are drawn by dissolving their areas.

## Usage

``` r
# S3 method for class 'stareal_pred'
st_as_sf(x, time, group = "area", ...)
```

## Arguments

- x:

  A `stareal_pred` whose frame has polygons.

- time:

  A time label or step number.

- group:

  `"area"` or one of the groupings predicted.

- ...:

  Not used.

## Value

An `sf` object, one row per area (or group).
