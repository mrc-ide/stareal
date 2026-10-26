# Build a spatial frame

A frame is the set of areal units the model is defined on, with their
adjacency graph and any groupings to aggregate predictions to (regions,
countries, ...). Build it either from polygons, in which case adjacency
is computed and survey coordinates can later be assigned to units, or
from a ready-made adjacency matrix.

## Usage

``` r
areal_frame(
  polygons = NULL,
  id = NULL,
  adjacency = NULL,
  ids = NULL,
  groups = NULL,
  queen = TRUE,
  link_islands = TRUE,
  link_within = NULL,
  link_crs = NULL,
  adjacency_snap = 0,
  link_simplify = NULL,
  verbose = TRUE
)
```

## Arguments

- polygons:

  An `sf` object with one row per areal unit (POLYGON or MULTIPOLYGON).
  Units are taken in row order. Needs the `sf` and `spdep` packages.

- id:

  For `polygons`: the name of the column holding unique unit ids.

- adjacency:

  Instead of `polygons`: a square, symmetric 0/1 matrix (base or
  `Matrix`) with zero diagonal, one row per unit.

- ids:

  For `adjacency`: unique unit ids, in row order. Defaults to the row
  names of `adjacency`, or `1..m`.

- groups:

  Groupings of the units, for aggregating predictions. For `polygons`:
  names of columns of `polygons`. If the vector is named, the names
  become the group names (e.g.
  `c(region = "NAME_1", country = "GID_0")`). For `adjacency`: a data
  frame with one row per unit in row order, one column per grouping.
  Grouping names become column names of the prediction tables, so names
  the package uses there are refused (`id`, `area`, `t`, `time`, `y`,
  `k`, `mean`, `sd`, ...).

- queen:

  For `polygons`: units sharing only a corner are neighbours (`TRUE`,
  queen) or not (`FALSE`, rook).

- link_islands:

  For `polygons`: join disconnected components to the largest one (see
  Details).

- link_within:

  For `polygons`: the name of a grouping inside which island links are
  made where possible, or `NULL`.

- link_crs:

  The projected coordinate reference system in which distances are
  measured. It is used for island links here, and for snapping survey
  points in
  [`areal_data()`](https://mrc-ide.github.io/stareal/reference/areal_data.md).
  Its units must be metres. The default is the projection of `polygons`
  if they are projected in metres (and not Web Mercator), and otherwise
  a Lambert azimuthal equal-area projection centred on the units
  (allowing for frames that cross the 180th meridian). For a frame
  spanning a large part of the globe, give a suitable projection here.

- adjacency_snap:

  For `polygons`: units whose borders come within this distance (metres,
  measured in `link_crs`) count as neighbours. The default, 0, needs
  shared borders to match exactly, as they do within one boundary file.
  Boundaries assembled from several sources (e.g. one file per country)
  rarely meet exactly at national borders: a snap of 1-2 km joins them.

- link_simplify:

  Tolerance (in `link_crs` units) used to simplify polygons before
  measuring island-link distances, which is much faster. A unit that
  vanishes when simplified keeps its full detail. The default is
  1/10,000 of the diagonal of the bounding box. Use 0 for full detail.

- verbose:

  Print progress messages.

## Value

An object of class `stareal_frame`:

- ids, m:

  unit ids (in unit order) and their number

- W:

  sparse symmetric 0/1 adjacency, island links included

- groups:

  data frame: `id` and one column per grouping

- links:

  data frame of added island edges, or `NULL`

- n_components:

  connected components of the final graph

- eig:

  eigendecomposition of the Laplacian: `values`, `vectors`

- polygons:

  the `sf` polygons (id and grouping columns), or `NULL`

- crs_planar:

  the projection used for distances, or `NULL`

- source:

  what the frame was built from

## Details

The spatial prior is a Leroux conditional autoregressive model on the
graph, which uses the eigendecomposition of the graph Laplacian \\R =
D - W\\. It is computed once here. It is a dense \\m \times m\\ problem,
and takes about a minute for a few thousand units with an optimised
BLAS.

**Islands.** A graph in several pieces leaves each piece without
neighbours to borrow strength from. With polygons and
`link_islands = TRUE` (the default), every component other than the
largest is joined to it by one edge. The edge runs from the component's
unit closest to the main component to the closest unit of the main
component, by polygon distance. If `link_within` names a grouping, that
unit is looked for first among main-component units in the same group
(e.g. the same country), and elsewhere only if the group has none. The
added edges are listed in `$links`. A frame built from an adjacency
matrix cannot be linked like this, and a disconnected one is accepted
with a warning.

**Save the frame.** Build a frame once and save it
([`saveRDS()`](https://rdrr.io/r/base/readRDS.html)): the
eigendecomposition is the slow step, and its eigenvectors depend on the
BLAS library, so a frame rebuilt on another machine gives slightly
different random draws for the same seed.

## See also

[`areal_lattice()`](https://mrc-ide.github.io/stareal/reference/areal_lattice.md)
for a regular grid,
[`areal_data()`](https://mrc-ide.github.io/stareal/reference/areal_data.md)
to attach survey data.

## Examples

``` r
# from an adjacency matrix: four units in a line, two regions
W <- matrix(0, 4, 4)
W[cbind(1:3, 2:4)] <- 1
W <- W + t(W)
fr <- areal_frame(adjacency = W, ids = c("a", "b", "c", "d"),
                  groups = data.frame(region = c("N", "N", "S", "S")))
#> eigendecomposition of the graph Laplacian (m = 4)
fr
#> <stareal_frame> 4 units, 3 edges (mean 1.5 neighbours), 1 connected component
#>   built from: adjacency matrix (no polygons)
#>   groupings: region (2) 
```
