# A regular lattice frame

A frame of `nrow` x `ncol` square cells, for examples, simulation
studies and testing. Units are numbered row by row from the top left.
Square polygons are attached if `sf` is installed, so the frame can be
mapped and survey coordinates assigned to it.

## Usage

``` r
areal_lattice(nrow, ncol, queen = TRUE, blocks = NULL, cellsize = 10000)
```

## Arguments

- nrow, ncol:

  Number of rows and columns.

- queen:

  Cells sharing only a corner are neighbours (`TRUE`) or not (`FALSE`).

- blocks:

  Optional `c(rows, cols)`: adds a grouping `block` of rectangular
  blocks of this size (e.g. `c(2, 5)`; one number gives square blocks),
  for aggregation.

- cellsize:

  Side of a cell in metres; the polygons are placed in a Lambert
  azimuthal equal-area projection centred on (0, 0) longitude and
  latitude.

## Value

A `stareal_frame` (see
[`areal_frame()`](https://mrc-ide.github.io/stareal/reference/areal_frame.md)).
Ids are `"r<i>c<j>"`.

## Examples

``` r
fr <- areal_lattice(4, 5, blocks = c(2, 5))
fr
#> <stareal_frame> 20 units, 55 edges (mean 5.5 neighbours), 1 connected component
#>   built from: 4 x 5 lattice
#>   groupings: block (2) 
```
