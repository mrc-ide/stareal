# frame.R
#
# The spatial frame: the areal units, their adjacency graph, optional
# groupings for aggregation, and the one-off eigendecomposition of the graph
# Laplacian that the sampler and prediction rely on.
#
# Polygon frames generalise Africa_areal/pipeline/01_frame.R (2026-10-12):
# queen (or rook) adjacency via spdep::poly2nb on planar geometry; each
# island component linked to the main component by one edge, from the
# island unit closest to the mainland to the closest mainland unit, by
# polygon distance in an equal-area projection, staying inside the island's
# own group (`link_within`, e.g. country) where that group has mainland
# units. data-raw/check_frame_africa.R checks that the two agree on the
# Africa ADMIN1 frame.
#
# ------------------------------------------------------------------

#' Build a spatial frame
#'
#' A frame is the set of areal units the model is defined on, with their
#' adjacency graph and any groupings to aggregate predictions to (regions,
#' countries, ...). Build it either from polygons, in which case adjacency
#' is computed and survey coordinates can later be assigned to units, or
#' from a ready-made adjacency matrix.
#'
#' The spatial prior is a Leroux conditional autoregressive model on the
#' graph, which uses the eigendecomposition of the graph Laplacian
#' \eqn{R = D - W}. It is computed once here. It is a dense \eqn{m \times m}
#' problem, and takes about a minute for a few thousand units with an
#' optimised BLAS.
#'
#' **Islands.** A graph in several pieces leaves each piece without
#' neighbours to borrow strength from. With polygons and
#' `link_islands = TRUE` (the default), every component other than the
#' largest is joined to it by one edge. The edge runs from the component's
#' unit closest to the main component to the closest unit of the main
#' component, by polygon distance. If `link_within` names a grouping, that
#' unit is looked for first among main-component units in the same group
#' (e.g. the same country), and elsewhere only if the group has none. The
#' added edges are listed in `$links`. A frame built from an adjacency
#' matrix cannot be linked like this, and a disconnected one is accepted
#' with a warning.
#'
#' @param polygons An `sf` object with one row per areal unit (POLYGON or
#'   MULTIPOLYGON). Units are taken in row order. Needs the `sf` and
#'   `spdep` packages.
#' @param id For `polygons`: the name of the column holding unique unit ids.
#' @param adjacency Instead of `polygons`: a square, symmetric 0/1 matrix
#'   (base or `Matrix`) with zero diagonal, one row per unit.
#' @param ids For `adjacency`: unique unit ids, in row order. Defaults to the
#'   row names of `adjacency`, or `1..m`.
#' @param groups Groupings of the units, for aggregating predictions. For
#'   `polygons`: names of columns of `polygons`. If the vector is named,
#'   the names become the group names (e.g.
#'   `c(region = "NAME_1", country = "GID_0")`). For `adjacency`: a data
#'   frame with one row per unit in row order, one column per grouping.
#' @param queen For `polygons`: units sharing only a corner are neighbours
#'   (`TRUE`, queen) or not (`FALSE`, rook).
#' @param link_islands For `polygons`: join disconnected components to the
#'   largest one (see Details).
#' @param link_within For `polygons`: the name of a grouping inside which
#'   island links are made where possible, or `NULL`.
#' @param link_crs The projected coordinate reference system in which
#'   distances are measured. It is used for island links here, and for
#'   snapping survey points in [areal_data()]. Its units must be metres. The
#'   default is the projection of `polygons` if they are projected, and
#'   otherwise a Lambert azimuthal equal-area projection centred on their
#'   bounding box.
#' @param link_simplify Tolerance (in `link_crs` units) used to simplify
#'   polygons before measuring island-link distances, which is much faster.
#'   A unit that vanishes when simplified keeps its full detail. The default
#'   is 1/10,000 of the diagonal of the bounding box. Use 0 for full detail.
#' @param verbose Print progress messages.
#'
#' @return An object of class `stareal_frame`:
#'   \describe{
#'     \item{ids, m}{unit ids (in unit order) and their number}
#'     \item{W}{sparse symmetric 0/1 adjacency, island links included}
#'     \item{groups}{data frame: `id` and one column per grouping}
#'     \item{links}{data frame of added island edges, or `NULL`}
#'     \item{n_components}{connected components of the final graph}
#'     \item{eig}{eigendecomposition of the Laplacian: `values`, `vectors`}
#'     \item{polygons}{the `sf` polygons (id and grouping columns), or `NULL`}
#'     \item{crs_planar}{the projection used for distances, or `NULL`}
#'   }
#' @seealso [areal_lattice()] for a regular grid, [areal_data()] to attach
#'   survey data.
#' @examples
#' # from an adjacency matrix: four units in a line, two regions
#' W <- matrix(0, 4, 4)
#' W[cbind(1:3, 2:4)] <- 1
#' W <- W + t(W)
#' fr <- areal_frame(adjacency = W, ids = c("a", "b", "c", "d"),
#'                   groups = data.frame(region = c("N", "N", "S", "S")))
#' fr
#' @export
areal_frame <- function(polygons = NULL, id = NULL, adjacency = NULL, ids = NULL, groups = NULL,
                        queen = TRUE, link_islands = TRUE, link_within = NULL,
                        link_crs = NULL, link_simplify = NULL, verbose = TRUE) {
  if (is.null(polygons) == is.null(adjacency))
    stop("give exactly one of `polygons` and `adjacency`", call. = FALSE)
  if (!is.null(polygons)) {
    fr <- frame_from_polygons(polygons, id, groups, queen, link_islands, link_within,
                              link_crs, link_simplify, verbose)
  } else {
    fr <- frame_from_adjacency(adjacency, ids, groups)
  }
  if (fr$n_components > 1)
    warning(sprintf(paste("the graph has %d connected components; units in the smaller ones borrow",
                          "no strength from the rest (component sizes: %s)"),
                    fr$n_components, paste(sort(tabulate(fr$component), decreasing = TRUE), collapse = ", ")),
            call. = FALSE)
  if (any(Matrix::rowSums(fr$W) == 0))
    warning(sprintf("%d units have no neighbours", sum(Matrix::rowSums(fr$W) == 0)), call. = FALSE)
  if (fr$m > 8000)
    warning(sprintf("m = %d units: the eigendecomposition is dense (m x m) and will be slow and large",
                    fr$m), call. = FALSE)
  say(verbose, "eigendecomposition of the graph Laplacian (m = %d)", fr$m)
  fr$eig <- laplacian_eigen(fr$W, fr$n_components)
  fr$component <- NULL
  structure(fr, class = "stareal_frame")
}

# the one-off eigendecomposition; a Laplacian has eigenvalues >= 0, with one
# zero per connected component
laplacian_eigen <- function(W, n_components) {
  R <- as.matrix(Matrix::Diagonal(x = Matrix::rowSums(W)) - W)
  eig <- eigen(R, symmetric = TRUE)
  tol <- 1e-8 * max(1, max(abs(eig$values)))
  stopifnot(min(eig$values) > -tol, sum(abs(eig$values) < tol) == n_components)
  eig$values[abs(eig$values) < tol] <- 0
  list(values = eig$values, vectors = eig$vectors)
}

check_ids <- function(ids) {
  ids <- as.character(ids)
  if (anyNA(ids) || any(ids == "")) stop("unit ids must not be missing or empty", call. = FALSE)
  if (anyDuplicated(ids)) stop(sprintf("unit ids must be unique; duplicated: %s",
                                       paste(utils::head(unique(ids[duplicated(ids)]), 5), collapse = ", ")),
                               call. = FALSE)
  ids
}

check_groups <- function(groups, m) {
  if (is.null(groups)) return(NULL)
  groups <- as.data.frame(groups, stringsAsFactors = FALSE)
  if (nrow(groups) != m) stop(sprintf("`groups` has %d rows; the frame has %d units", nrow(groups), m),
                              call. = FALSE)
  if (any(names(groups) %in% c("id", "")) || anyDuplicated(names(groups)))
    stop("grouping names must be unique, non-empty and not 'id'", call. = FALSE)
  groups[] <- lapply(groups, as.character)
  bad <- vapply(groups, anyNA, logical(1))
  if (any(bad)) stop(sprintf("grouping(s) with missing values: %s", paste(names(groups)[bad], collapse = ", ")),
                     call. = FALSE)
  groups
}

frame_from_adjacency <- function(adjacency, ids, groups) {
  m <- nrow(adjacency)
  if (is.null(m) || ncol(adjacency) != m) stop("`adjacency` must be a square matrix", call. = FALSE)
  if (is.null(ids)) ids <- if (!is.null(rownames(adjacency))) rownames(adjacency) else seq_len(m)
  if (length(ids) != m) stop(sprintf("%d ids for %d units", length(ids), m), call. = FALSE)
  ids <- check_ids(ids)
  W <- as_dgC(adjacency)
  if (any(is.na(W@x)) || !all(W@x == 1)) stop("`adjacency` must contain only 0 and 1", call. = FALSE)
  if (any(Matrix::diag(W) != 0)) stop("`adjacency` must have a zero diagonal", call. = FALSE)
  if (!Matrix::isSymmetric(W)) stop("`adjacency` must be symmetric", call. = FALSE)
  dimnames(W) <- NULL
  groups <- check_groups(groups, m)
  comp <- graph_components(adjacency_list(W))
  gtab <- data.frame(id = ids, stringsAsFactors = FALSE)
  if (!is.null(groups)) gtab <- cbind(gtab, groups)
  list(ids = ids, m = m, W = W, groups = gtab,
       links = NULL, n_components = max(comp), component = comp,
       polygons = NULL, crs_planar = NULL)
}

frame_from_polygons <- function(polygons, id, groups, queen, link_islands, link_within,
                                link_crs, link_simplify, verbose) {
  need_pkg("sf", "for frames built from polygons")
  need_pkg("spdep", "for frames built from polygons")
  if (!inherits(polygons, "sf")) stop("`polygons` must be an sf object", call. = FALSE)
  if (is.null(id) || !id %in% names(polygons)) stop("`id` must name a column of `polygons`", call. = FALSE)
  ids <- check_ids(polygons[[id]])
  m <- length(ids)
  gt <- unique(as.character(sf::st_geometry_type(polygons)))
  if (!all(gt %in% c("POLYGON", "MULTIPOLYGON")))
    stop(sprintf("`polygons` must be POLYGON/MULTIPOLYGON, not %s", paste(gt, collapse = ", ")), call. = FALSE)
  if (is.na(sf::st_crs(polygons))) stop("`polygons` has no coordinate reference system", call. = FALSE)

  if (!is.null(groups)) {
    gcols <- unname(groups)
    gnames <- if (is.null(names(groups))) gcols else ifelse(names(groups) == "", gcols, names(groups))
    miss <- setdiff(gcols, names(polygons))
    if (length(miss)) stop(sprintf("grouping column(s) not in `polygons`: %s", paste(miss, collapse = ", ")),
                           call. = FALSE)
    groups <- stats::setNames(sf::st_drop_geometry(polygons)[gcols], gnames)
  }
  groups <- check_groups(groups, m)
  if (!is.null(link_within) && !link_within %in% names(groups))
    stop("`link_within` must be the name of one of the groupings", call. = FALSE)

  geom <- sf::st_geometry(polygons)
  gtab <- data.frame(id = ids, stringsAsFactors = FALSE)
  if (!is.null(groups)) gtab <- cbind(gtab, groups)
  poly <- sf::st_sf(gtab, geometry = geom)

  with_s2_off({
    if (!all(sf::st_is_valid(poly))) {
      say(verbose, "repairing %d invalid polygons (sf::st_make_valid)", sum(!sf::st_is_valid(poly)))
      poly <- sf::st_make_valid(poly)
    }
    say(verbose, "adjacency (%s) of %d units", if (queen) "queen" else "rook", m)
    nb <- suppressMessages(spdep::poly2nb(poly, queen = queen))     # (sf's planar-on-lon/lat note: deliberate)
    nb <- lapply(seq_len(m), function(i) setdiff(as.integer(nb[[i]]), 0L))   # 0L: spdep's "no neighbours"
    for (i in seq_len(m)) for (j in nb[[i]]) nb[[j]] <- union(nb[[j]], i)    # symmetrise (poly2nb nearly always is)

    crs_planar <- planar_crs(poly, link_crs)
    comp <- graph_components(nb)
    links <- NULL
    if (link_islands && max(comp) > 1) {
      res <- link_components(poly, nb, comp, crs_planar, link_within, link_simplify, verbose)
      nb <- res$nb; links <- res$links
      comp <- graph_components(nb)
    }
  })
  W <- Matrix::sparseMatrix(i = rep(seq_len(m), lengths(nb)), j = unlist(nb), x = 1, dims = c(m, m))
  W <- as_dgC(W)
  say(verbose, "%d units, %d edges, mean %.1f neighbours, %d component%s", m, length(W@x) / 2,
      mean(lengths(nb)), max(comp), if (max(comp) > 1) "s" else "")
  list(ids = ids, m = m, W = W, groups = sf::st_drop_geometry(poly), links = links,
       n_components = max(comp), component = comp, polygons = poly, crs_planar = crs_planar)
}

# the projection used for distances: the polygons' own if projected,
# otherwise Lambert azimuthal equal-area centred on their bounding box
planar_crs <- function(poly, link_crs) {
  if (!is.null(link_crs)) return(sf::st_crs(link_crs))
  if (!isTRUE(sf::st_is_longlat(poly))) return(sf::st_crs(poly))
  bb <- sf::st_bbox(sf::st_transform(poly, 4326))
  sf::st_crs(sprintf("+proj=laea +lat_0=%.4f +lon_0=%.4f +datum=WGS84 +units=m +no_defs",
                     (bb[["ymin"]] + bb[["ymax"]]) / 2, (bb[["xmin"]] + bb[["xmax"]]) / 2))
}

# one edge from every non-main component to the main one (see areal_frame)
link_components <- function(poly, nb, comp, crs_planar, link_within, link_simplify, verbose) {
  proj <- sf::st_transform(poly, crs_planar)
  if (is.null(link_simplify)) {
    bb <- sf::st_bbox(proj)
    link_simplify <- 1e-4 * sqrt((bb[["xmax"]] - bb[["xmin"]])^2 + (bb[["ymax"]] - bb[["ymin"]])^2)
  }
  simp <- proj
  if (link_simplify > 0) {
    simp <- sf::st_simplify(proj, dTolerance = link_simplify)
    empty <- sf::st_is_empty(simp)                 # small islands that vanish keep full detail
    sf::st_geometry(simp)[empty] <- sf::st_geometry(proj)[empty]
  }
  grp <- if (is.null(link_within)) NULL else poly[[link_within]]
  links <- list()
  for (cid in setdiff(sort(unique(comp)), 1L)) {
    isl <- which(comp == cid)
    mland <- which(comp == 1L)
    if (!is.null(grp)) {
      same <- mland[grp[mland] %in% grp[isl]]
      if (length(same)) mland <- same
    }
    d <- sf::st_distance(simp[isl, ], simp[mland, ])
    d <- matrix(as.numeric(d), nrow(d))
    ij <- arrayInd(which.min(d), dim(d))[1, ]
    i <- isl[ij[1]]; j <- mland[ij[2]]
    nb[[i]] <- union(nb[[i]], j); nb[[j]] <- union(nb[[j]], i)
    links[[length(links) + 1]] <- data.frame(from = poly$id[i], to = poly$id[j], n_units = length(isl),
                                             distance_km = round(d[ij[1], ij[2]] / 1e3, 1),
                                             stringsAsFactors = FALSE)
  }
  links <- do.call(rbind, links)
  say(verbose, "linked %d island component%s to the main component", nrow(links),
      if (nrow(links) > 1) "s" else "")
  list(nb = nb, links = links)
}

#' A regular lattice frame
#'
#' A frame of `nrow` x `ncol` square cells, for examples, simulation
#' studies and testing. Units are numbered row by row from the top left.
#' Square polygons are attached if `sf` is installed, so the frame can be
#' mapped and survey coordinates assigned to it.
#'
#' @param nrow,ncol Number of rows and columns.
#' @param queen Cells sharing only a corner are neighbours (`TRUE`) or not
#'   (`FALSE`).
#' @param blocks Optional `c(rows, cols)`: adds a grouping `block` of
#'   rectangular blocks of this size (e.g. `c(2, 2)`), for aggregation.
#' @param cellsize Side of a cell in metres; the polygons are placed in a
#'   Lambert azimuthal equal-area projection centred on (0, 0) longitude
#'   and latitude.
#' @param verbose Print progress messages.
#' @return A `stareal_frame` (see [areal_frame()]). Ids are `"r<i>c<j>"`.
#' @examples
#' fr <- areal_lattice(4, 5, blocks = c(2, 5))
#' fr
#' @export
areal_lattice <- function(nrow, ncol, queen = TRUE, blocks = NULL, cellsize = 10000, verbose = FALSE) {
  stopifnot(nrow >= 1, ncol >= 1, nrow * ncol >= 2)
  r <- rep(seq_len(nrow), each = ncol); cc <- rep(seq_len(ncol), times = nrow)
  m <- nrow * ncol
  ids <- sprintf("r%dc%d", r, cc)
  d_r <- abs(outer(r, r, "-")); d_c <- abs(outer(cc, cc, "-"))
  W <- 1 * ((d_r + d_c == 1) | (queen & d_r == 1 & d_c == 1))
  groups <- if (is.null(blocks)) NULL else
    data.frame(block = sprintf("b%d_%d", (r - 1) %/% blocks[1] + 1, (cc - 1) %/% blocks[2] + 1))
  fr <- areal_frame(adjacency = W, ids = ids, groups = groups, verbose = verbose)
  if (requireNamespace("sf", quietly = TRUE)) {
    crs <- sf::st_crs("+proj=laea +lat_0=0 +lon_0=0 +datum=WGS84 +units=m +no_defs")
    sq <- lapply(seq_len(m), function(a) {
      x0 <- (cc[a] - 1) * cellsize; y0 <- (nrow - r[a]) * cellsize
      sf::st_polygon(list(cbind(c(x0, x0 + cellsize, x0 + cellsize, x0, x0),
                                c(y0, y0, y0 + cellsize, y0 + cellsize, y0))))
    })
    fr$polygons <- sf::st_sf(fr$groups, geometry = sf::st_sfc(sq, crs = crs))
    fr$crs_planar <- crs
  }
  fr
}

#' @export
print.stareal_frame <- function(x, ...) {
  deg <- Matrix::rowSums(x$W)
  cat(sprintf("<stareal_frame> %d units, %d edges (mean %.1f neighbours), %d connected component%s\n",
              x$m, length(x$W@x) / 2, mean(deg), x$n_components, if (x$n_components > 1) "s" else ""))
  cat(sprintf("  built from: %s\n", if (is.null(x$polygons)) "adjacency matrix" else "polygons"))
  if (!is.null(x$links)) cat(sprintf("  island links added: %d\n", nrow(x$links)))
  g <- setdiff(names(x$groups), "id")
  if (length(g)) {
    cat("  groupings:", paste(sprintf("%s (%d)", g, vapply(x$groups[g], function(v) length(unique(v)), 1L)),
                              collapse = ", "), "\n")
  } else cat("  groupings: none\n")
  invisible(x)
}

#' Map a frame and its graph
#'
#' Draws the polygons, the adjacency edges between unit centres, and any
#' island links (in red). Needs a frame built from polygons (or a lattice
#' with `sf` installed).
#' @param x A `stareal_frame`.
#' @param edges Draw the adjacency edges.
#' @param ... Passed to the polygon plot.
#' @return `x`, invisibly.
#' @export
plot.stareal_frame <- function(x, edges = TRUE, ...) {
  if (is.null(x$polygons)) stop("this frame has no polygons to plot", call. = FALSE)
  need_pkg("sf", "to plot a frame")
  geom <- sf::st_geometry(x$polygons)
  graphics::plot(geom, border = "grey60", lwd = 0.3, ...)
  if (edges || !is.null(x$links)) {
    xy <- with_s2_off(suppressWarnings(sf::st_coordinates(sf::st_point_on_surface(geom))))
    if (edges) {
      e <- Matrix::which(Matrix::triu(x$W) != 0, arr.ind = TRUE)
      graphics::segments(xy[e[, 1], 1], xy[e[, 1], 2], xy[e[, 2], 1], xy[e[, 2], 2], col = "grey40", lwd = 0.4)
    }
    if (!is.null(x$links)) {
      i <- match(x$links$from, x$ids); j <- match(x$links$to, x$ids)
      graphics::segments(xy[i, 1], xy[i, 2], xy[j, 1], xy[j, 2], col = "red", lwd = 1.5)
    }
  }
  invisible(x)
}
