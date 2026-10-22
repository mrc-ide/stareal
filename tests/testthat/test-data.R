# Checks of data.R: pooling, cell ids, every exclusion reason, and
# assignment of coordinates (inside, snapped, too far, lon/lat input).

fr <- areal_lattice(3, 3, blocks = c(3, 1))
ti <- time_frame(T = 4)

test_that("surveys by id and step are pooled into cells with the right ids", {
  s <- data.frame(unit = c("r1c1", "r1c1", "r2c3", "r1c1"), step = c(1, 1, 4, 2),
                  pos = c(3, 1, 0, 2), n = c(10, 5, 8, 9))
  d <- areal_data(s, fr, ti, y = "pos", k = "n", area = "unit", t = "step", verbose = FALSE)
  expect_s3_class(d, "stareal_data")
  cl <- d$cells
  expect_identical(cl$id, c("r1c1", "r1c1", "r2c3"))
  expect_identical(cl$t, c(1L, 2L, 4L))
  expect_identical(cl$y, c(4L, 2L, 0L)); expect_identical(cl$k, c(15L, 9L, 8L))
  expect_identical(cl$n_surveys, c(2L, 1L, 1L))
  expect_identical(cl$cell_id, (cl$area - 1L) * ti$T + cl$t)
  expect_true(all(is.na(d$surveys$.excluded)))
  expect_identical(as.data.frame(d), cl)
  expect_output(print(d), "4 of 4 surveys used, pooled into 3 cells")
})

test_that("every exclusion reason is recorded, with warnings for data errors", {
  s <- data.frame(unit = c("r1c1", "r1c1", "r1c1", "r1c1", "r1c1", NA, "zz", "r2c2", "r2c2", "r3c3"),
                  step = c(1, 1, 1, 1, 1, 1, 1, NA, 9, 2),
                  pos = c(NA, 1.5, 0, 6, 1, 1, 1, 1, 1, 1), n = c(5, 5, 0, 5, 5, 5, 5, 5, 5, 5))
  expect_warning(expect_warning(
    d <- areal_data(s, fr, ti, y = "pos", k = "n", area = "unit", t = "step", verbose = FALSE),
    "invalid counts"), "not in the frame")
  expect_identical(d$surveys$.excluded,
                   c("missing count", "invalid count (not a whole number >= 0)", "none tested (k = 0)",
                     "positives exceed tested (y > k)", NA, "missing area", "area not in frame",
                     "missing time step", "outside time window", NA))
  expect_equal(nrow(d$cells), 2)
  expect_error(areal_data(s[1:4, ], fr, ti, y = "pos", k = "n", area = "unit", t = "step", verbose = FALSE),
               "no survey could be used") |> suppressWarnings()
})

test_that("dates are mapped to calendar steps", {
  tm <- time_frame("2020-01", "2020-06", "month")
  s <- data.frame(unit = "r1c1", d = as.Date(c("2020-01-31", "2020-02-01", "2020-07-01", NA)), pos = 1, n = 2)
  d <- areal_data(s, fr, tm, y = "pos", k = "n", area = "unit", date = "d", verbose = FALSE)
  expect_identical(d$cells$time, c("2020-01", "2020-02"))
  expect_identical(d$surveys$.excluded[3:4], c("outside time window", "missing date"))
  expect_error(areal_data(s, fr, ti, y = "pos", k = "n", area = "unit", date = "d"), "plain index")
})

test_that("argument errors", {
  s <- data.frame(unit = "r1c1", step = 1, pos = 1, n = 2)
  expect_error(areal_data(s, fr, ti, y = "pos", k = "n", t = "step"), "`area` and `coords`")
  expect_error(areal_data(s, fr, ti, y = "pos", k = "n", area = "unit"), "`date` and `t`")
  expect_error(areal_data(s, fr, ti, y = "pos", k = "nope", area = "unit", t = "step"), "nope")
  expect_error(areal_data(s, list(), ti, y = "pos", k = "n", area = "unit", t = "step"), "areal_frame")
})

test_that("coordinates are assigned by point-in-polygon, with snapping", {
  skip_if_not_installed("sf")
  crs <- sf::st_crs(fr$polygons)                     # metres; x 0-30 km, y 0-30 km, row 1 at the top
  s <- data.frame(x = c(5e3, 25e3, 31e3, 80e3, NA), y = c(25e3, 5e3, 15e3, 15e3, 1),
                  step = 1, pos = 1, n = 2)
  d <- areal_data(s, fr, ti, y = "pos", k = "n", coords = c("x", "y"), crs = crs, t = "step",
                  snap_max_m = 5000, verbose = FALSE)
  expect_identical(fr$ids[d$surveys$.area[1:3]], c("r1c1", "r3c3", "r2c3"))
  expect_equal(d$surveys$.snap_m[1:4], c(0, 0, 1000, 50000))
  expect_identical(d$surveys$.excluded[4:5], c("outside frame (> 5 km from every unit)", "missing coordinates"))
  # the same points given as longitude / latitude
  ll <- sf::st_coordinates(sf::st_transform(sf::st_as_sf(s[1:4, ], coords = c("x", "y"), crs = crs), 4326))
  s_ll <- cbind(s[1:4, ], lon = ll[, 1], lat = ll[, 2])
  d2 <- areal_data(s_ll, fr, ti, y = "pos", k = "n", coords = c("lon", "lat"), t = "step",
                   snap_max_m = 5000, verbose = FALSE)
  expect_identical(d2$surveys$.area, d$surveys$.area[1:4])
  # coordinates need polygons
  fa <- areal_frame(adjacency = as.matrix(fr$W), verbose = FALSE)
  expect_error(areal_data(s, fa, ti, y = "pos", k = "n", coords = c("x", "y"), t = "step"), "polygons")
})
