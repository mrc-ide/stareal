# Checks of data.R: pooling, cell ids, every exclusion reason, input
# forms (factors, labels, date-times, numeric ids), and assignment of
# coordinates (inside, snapped, too far, lon/lat input, unreadable text).

fr <- areal_lattice(3, 3, blocks = c(3, 1))
ti <- time_frame(n_steps = 4)

test_that("surveys by id and step are pooled into cells with the right ids", {
  s <- data.frame(unit = c("r1c1", "r1c1", "r2c3", "r1c1"), step = c(1, 1, 4, 2),
                  pos = c(3, 1, 0, 2), n = c(10, 5, 8, 9))
  expect_message(d <- areal_data(s, fr, ti, y = "pos", k = "n", area = "unit", t = "step"),
                 "4 of 4 surveys used, pooled into 3 cells")
  expect_s3_class(d, "stareal_data")
  cl <- d$cells
  expect_identical(cl$id, c("r1c1", "r1c1", "r2c3"))
  expect_identical(cl$t, c(1L, 2L, 4L))
  expect_equal(cl$y, c(4, 2, 0)); expect_equal(cl$k, c(15, 9, 8))
  expect_identical(cl$n_surveys, c(2L, 1L, 1L))
  expect_identical(cl$cell_id, (cl$area - 1L) * ti$n_steps + cl$t)
  expect_true(all(is.na(d$surveys$.excluded)))
  expect_identical(d$surveys$.id, c("r1c1", "r1c1", "r2c3", "r1c1"))
  expect_identical(as.data.frame(d), cl)
  expect_output(print(d), "4 of 4 surveys used, pooled into 3 cells")
  expect_silent(areal_data(s, fr, ti, y = "pos", k = "n", area = "unit", t = "step", verbose = FALSE))
})

test_that("every exclusion reason is recorded, with warnings for data errors", {
  s <- data.frame(unit = c("r1c1", "r1c1", "r1c1", "r1c1", "r1c1", NA, "zz", "r2c2", "r2c2", "r3c3", "r3c3"),
                  step = c(1, 1, 1, 1, 1, 1, 1, NA, 9, 2, "two"),
                  pos = c(NA, 1.5, 0, 6, 1, 1, 1, 1, 1, 1, 1), n = c(5, 5, 0, 5, 5, 5, 5, 5, 5, 5, 5))
  w <- character(0)
  d <- withCallingHandlers(areal_data(s, fr, ti, y = "pos", k = "n", area = "unit", t = "step", verbose = FALSE),
                           warning = function(e) { w <<- c(w, conditionMessage(e)); invokeRestart("muffleWarning") })
  expect_true(any(grepl("invalid counts", w)) && any(grepl("not in the frame", w)) && any(grepl("unreadable", w)))
  expect_identical(d$surveys$.excluded,
                   c("missing or unreadable count", "invalid count (not a whole number >= 0)", "none tested (k = 0)",
                     "positives exceed tested (y > k)", NA, "missing area", "area not in frame",
                     "missing time step", "outside time window", NA,
                     "time step not a step number or label of the time frame"))
  expect_equal(nrow(d$cells), 2)
  expect_error(areal_data(s[1:4, ], fr, ti, y = "pos", k = "n", area = "unit", t = "step", verbose = FALSE),
               "no survey could be used") |> suppressWarnings()
})

test_that("factor columns are read by their labels, not their codes", {
  s <- data.frame(unit = factor(c("r1c1", "r2c2")), step = factor(c("3", "10")),
                  pos = factor(c("3", "10")), n = factor(c("20", "20")))
  d <- areal_data(s, fr, time_frame(n_steps = 12), y = "pos", k = "n", area = "unit", t = "step", verbose = FALSE)
  expect_equal(d$cells$y, c(3, 10)); expect_equal(d$cells$t, c(3L, 10L))
})

test_that("time steps can be given as labels; dates as Date, text or date-times", {
  tm <- time_frame("2020-01", "2020-06", "month")
  s <- data.frame(unit = "r1c1", lab = c("2020-01", "2020-03", "2021-01"), pos = 1, n = 2)
  d <- suppressWarnings(areal_data(s, fr, tm, y = "pos", k = "n", area = "unit", t = "lab", verbose = FALSE))
  expect_identical(d$cells$time, c("2020-01", "2020-03"))
  s <- data.frame(unit = "r1c1", d = as.Date(c("2020-01-31", "2020-02-01", "2020-07-01", NA)), pos = 1, n = 2)
  d <- areal_data(s, fr, tm, y = "pos", k = "n", area = "unit", date = "d", verbose = FALSE)
  expect_identical(d$cells$time, c("2020-01", "2020-02"))
  expect_identical(d$surveys$.excluded[3:4], c("outside time window", "missing date"))
  s$d <- c("2020-01-31", "2020-02-01", "not a date", NA)
  d <- suppressWarnings(areal_data(s, fr, tm, y = "pos", k = "n", area = "unit", date = "d", verbose = FALSE))
  expect_identical(d$surveys$.excluded[3], "unreadable date (expected YYYY-MM-DD)")
  s2 <- data.frame(unit = "r1c1", d = as.POSIXct("2020-02-01 00:00", tz = "Africa/Nairobi"), pos = 1, n = 2)
  expect_identical(areal_data(s2, fr, tm, y = "pos", k = "n", area = "unit", date = "d", verbose = FALSE)$cells$time,
                   "2020-02")
  s3 <- data.frame(unit = "r1c1", d = 2020, pos = 1, n = 2)
  expect_error(areal_data(s3, fr, tm, y = "pos", k = "n", area = "unit", date = "d"), "numeric")
  expect_error(areal_data(s, fr, ti, y = "pos", k = "n", area = "unit", date = "d"), "plain index")
  # years in a step column of a yearly frame: excluded, with a hint
  ty <- time_frame(2001, 2005, "year")
  s4 <- data.frame(unit = "r1c1", yr = c(2001, 2002), pos = 1, n = 2)
  expect_warning(expect_error(areal_data(s4, fr, ty, y = "pos", k = "n", area = "unit", t = "yr", verbose = FALSE),
                              "no survey"), "not years")
})

test_that("numeric ids match whatever their storage", {
  W <- matrix(0, 2, 2); W[1, 2] <- W[2, 1] <- 1
  f2 <- areal_frame(adjacency = W, ids = c(100000, 2500000), verbose = FALSE)
  expect_identical(f2$ids, c("100000", "2500000"))
  s <- data.frame(unit = c(100000L, 2500000L), step = 1, pos = 1, n = 2)
  expect_equal(nrow(areal_data(s, f2, time_frame(n_steps = 1), y = "pos", k = "n", area = "unit", t = "step",
                               verbose = FALSE)$cells), 2)
})

test_that("pooled counts beyond the integer range stay exact", {
  s <- data.frame(unit = "r1c1", step = 1, pos = c(1.5e9, 1e9), n = c(2e9, 2e9))
  d <- areal_data(s, fr, ti, y = "pos", k = "n", area = "unit", t = "step", verbose = FALSE)
  expect_equal(d$cells$k, 4e9); expect_equal(d$cells$y, 2.5e9)
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
  s <- data.frame(x = c(5e3, 25e3, 31e3, 80e3, NA, "north"), y = c(25e3, 5e3, 15e3, 15e3, 1, 1),
                  step = 1, pos = 1, n = 2)
  d <- areal_data(s, fr, ti, y = "pos", k = "n", coords = c("x", "y"), crs = crs, t = "step",
                  snap_max_m = 5000, verbose = FALSE)
  expect_identical(d$surveys$.id[1:3], c("r1c1", "r3c3", "r2c3"))
  expect_equal(d$surveys$.snap_m[1:4], c(0, 0, 1000, 50000))
  expect_identical(d$surveys$.excluded[4:6], c("outside frame (> 5 km from every unit)",
                                               "missing or unreadable coordinates", "missing or unreadable coordinates"))
  # the same points given as longitude / latitude
  s <- s[1:4, ]; s$x <- as.numeric(s$x)
  ll <- sf::st_coordinates(sf::st_transform(sf::st_as_sf(s, coords = c("x", "y"), crs = crs), 4326))
  s_ll <- cbind(s, lon = ll[, 1], lat = ll[, 2])
  d2 <- areal_data(s_ll, fr, ti, y = "pos", k = "n", coords = c("lon", "lat"), t = "step",
                   snap_max_m = 5000, verbose = FALSE)
  expect_identical(d2$surveys$.id, d$surveys$.id[1:4])
  # coordinates need polygons
  fa <- areal_frame(adjacency = as.matrix(fr$W), verbose = FALSE)
  expect_error(areal_data(s, fa, ti, y = "pos", k = "n", coords = c("x", "y"), t = "step"), "polygons")
})
