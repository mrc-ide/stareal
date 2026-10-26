# Map a frame and its graph

Draws the polygons, the adjacency edges between unit centres, and any
island links (in red). Needs a frame built from polygons (or a lattice
with `sf` installed).

## Usage

``` r
# S3 method for class 'stareal_frame'
plot(x, edges = TRUE, ...)
```

## Arguments

- x:

  A `stareal_frame`.

- edges:

  Draw the adjacency edges.

- ...:

  Passed to the polygon plot (e.g. `border`, `lwd`, `main`).

## Value

`x`, invisibly.
