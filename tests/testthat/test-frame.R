# Checks of frame.R: frames from adjacency matrices and from polygons agree,
# the eigendecomposition is that of the Laplacian, island linking picks the
# nearest unit (within the group where asked), and bad input is refused.

line_W <- function(m) {
  W <- matrix(0, m, m); W[cbind(1:(m - 1), 2:m)] <- 1; W + t(W)
}

test_that("an adjacency frame stores the graph, ids, groups and the Laplacian eigendecomposition", {
  W <- line_W(5); rownames(W) <- letters[1:5]
  fr <- areal_frame(adjacency = W, groups = data.frame(region = c("N", "N", "N", "S", "S")), verbose = FALSE)
  expect_s3_class(fr, "stareal_frame")
  expect_identical(fr$ids, letters[1:5])
  expect_equal(as.matrix(fr$W), unname(W))
  expect_identical(fr$groups$region, c("N", "N", "N", "S", "S"))
  R <- diag(rowSums(W)) - W
  expect_lt(max(abs(fr$eig$vectors %*% diag(fr$eig$values) %*% t(fr$eig$vectors) - R)), 1e-10)
  expect_equal(sum(fr$eig$values == 0), 1)
  expect_output(print(fr), "5 units, 4 edges")
})

test_that("bad adjacency input is refused", {
  W <- line_W(4)
  expect_error(areal_frame(adjacency = W, polygons = 1), "exactly one")
  expect_error(areal_frame(adjacency = W[, -1]), "square")
  A <- W; A[1, 2] <- 0; expect_error(areal_frame(adjacency = A, verbose = FALSE), "symmetric")
  A <- W; A[1, 1] <- 1; expect_error(areal_frame(adjacency = A, verbose = FALSE), "diagonal")
  A <- W * 2; expect_error(areal_frame(adjacency = A, verbose = FALSE), "0 and 1")
  expect_error(areal_frame(adjacency = W, ids = c("a", "a", "b", "c"), verbose = FALSE), "unique")
  expect_error(areal_frame(adjacency = W, ids = 1:3, verbose = FALSE), "3 ids for 4")
  expect_error(areal_frame(adjacency = W, groups = data.frame(g = 1:3), verbose = FALSE), "3 rows")
  expect_error(areal_frame(adjacency = W, groups = data.frame(g = c(1, NA, 1, 1)), verbose = FALSE), "missing")
})

test_that("a disconnected adjacency frame is accepted with a warning", {
  W <- matrix(0, 5, 5); W[1, 2] <- W[2, 1] <- W[2, 3] <- W[3, 2] <- W[4, 5] <- W[5, 4] <- 1
  expect_warning(fr <- areal_frame(adjacency = W, verbose = FALSE), "2 connected components")
  expect_equal(fr$n_components, 2)
  expect_equal(sum(fr$eig$values == 0), 2)
  W[4, 5] <- W[5, 4] <- 0
  expect_warning(expect_warning(areal_frame(adjacency = W, verbose = FALSE), "components"), "no neighbours")
})

test_that("a lattice has the right queen and rook graphs", {
  fq <- areal_lattice(3, 4); fr <- areal_lattice(3, 4, queen = FALSE)
  expect_equal(fq$m, 12)
  expect_equal(length(fr$W@x) / 2, 3 * 3 + 2 * 4)            # rook: horizontal + vertical edges
  expect_equal(length(fq$W@x) / 2, 17 + 2 * 3 * 2)           # queen adds two diagonals per 2 x 2 square
  expect_identical(fq$ids[c(1, 4, 5)], c("r1c1", "r1c4", "r2c1"))
  fb <- areal_lattice(4, 4, blocks = c(2, 2))
  expect_equal(length(unique(fb$groups$block)), 4)
})

test_that("frames from polygons reproduce the lattice graph (queen and rook, projected and lon/lat)", {
  skip_if_not_installed("sf"); skip_if_not_installed("spdep")
  for (queen in c(TRUE, FALSE)) {
    lat <- areal_lattice(3, 4, queen = queen)
    fp <- areal_frame(lat$polygons, id = "id", queen = queen, verbose = FALSE)
    expect_equal(as.matrix(fp$W), as.matrix(lat$W))
    expect_identical(fp$ids, lat$ids)
    expect_lt(max(abs(fp$eig$values - lat$eig$values)), 1e-10)
  }
  ll <- sf::st_transform(areal_lattice(3, 4)$polygons, 4326)
  fl <- areal_frame(ll, id = "id", verbose = FALSE)
  expect_equal(as.matrix(fl$W), as.matrix(areal_lattice(3, 4)$W))
  expect_match(fl$crs_planar$proj4string, "laea")
})

test_that("islands are linked to the nearest mainland unit, within their group where asked", {
  skip_if_not_installed("sf"); skip_if_not_installed("spdep")
  lat <- areal_lattice(3, 3)                        # 10 km cells, x 0-30 km, y 0-30 km; row 1 at the top
  sq <- function(x0, y0, s = 1e4) sf::st_polygon(list(cbind(c(x0, x0 + s, x0 + s, x0, x0),
                                                            c(y0, y0, y0 + s, y0 + s, y0))))
  isl <- sf::st_sf(id = "island", geometry = sf::st_sfc(sq(5e4, 2e4), crs = sf::st_crs(lat$polygons)))
  poly <- rbind(lat$polygons[, "id"], isl)
  poly$country <- ifelse(poly$id %in% c("r1c1", "r2c1", "r3c1", "island"), "A", "B")

  fr <- areal_frame(poly, id = "id", groups = c(country = "country"), verbose = FALSE)
  expect_equal(fr$n_components, 1)
  expect_equal(nrow(fr$links), 1)
  expect_identical(unlist(fr$links[c("from", "to")]), c(from = "island", to = "r1c3"))   # 20 km away
  expect_equal(fr$links$distance_km, 20)

  fr_w <- areal_frame(poly, id = "id", groups = c(country = "country"), link_within = "country", verbose = FALSE)
  expect_identical(fr_w$links$to, "r1c1")                                                   # 40 km, same country
  expect_equal(fr_w$links$distance_km, 40)

  expect_warning(expect_warning(f0 <- areal_frame(poly, id = "id", link_islands = FALSE, verbose = FALSE),
                                "components"), "no neighbours")
  expect_equal(f0$n_components, 2)
})

test_that("bad polygon input is refused", {
  skip_if_not_installed("sf"); skip_if_not_installed("spdep")
  p <- areal_lattice(2, 2)$polygons
  expect_error(areal_frame(p, id = "nope"), "`id`")
  expect_error(areal_frame(p, id = "id", groups = "nope", verbose = FALSE), "not in `polygons`")
  expect_error(areal_frame(p, id = "id", link_within = "nope", verbose = FALSE), "link_within")
  p2 <- p; p2$id[2] <- p2$id[1]
  expect_error(areal_frame(p2, id = "id", verbose = FALSE), "unique")
  expect_error(areal_frame(data.frame(id = 1:2), id = "id"), "sf object")
})
