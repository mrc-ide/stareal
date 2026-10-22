# Checks of time.R: step counts, labels, date -> step mapping at the window
# edges, and input forms.

test_that("monthly, quarterly and yearly frames have the right steps and labels", {
  tm <- time_frame("1995-01", "2024-12", "month")
  expect_equal(tm$T, 360)
  expect_identical(tm$labels[c(1, 360)], c("1995-01", "2024-12"))
  tq <- time_frame("2019-11-20", "2021-02-01", "quarter")      # moved to the starts of their quarters
  expect_identical(tq$labels, c("2019-Q4", paste0("2020-Q", 1:4), "2021-Q1"))
  ty <- time_frame(2001, 2020, "year")
  expect_equal(ty$T, 20)
  expect_identical(ty$starts[1], as.Date("2001-01-01"))
  expect_equal(time_frame(as.Date("2020-02-29"), "2020", "month")$T, 11)   # "2020" as end = all of 2020
  expect_equal(time_frame(2001, 2020, "month")$T, 240)
  expect_equal(time_frame("2020-03", "2020-03", "month")$T, 1)
  expect_error(time_frame(2020, 2019, "year"), "before")
  expect_output(print(tm), "monthly steps, 1995-01 to 2024-12 \\(T = 360\\)")
})

test_that("dates map to the right step, and NA outside the window", {
  tm <- time_frame("1995-01", "2024-12", "month")
  d <- as.Date(c("1994-12-31", "1995-01-01", "1995-01-31", "1995-02-01", "2024-12-31", "2025-01-01", NA))
  expect_identical(time_index(tm, d), c(NA, 1L, 1L, 2L, 360L, NA, NA))
  tq <- time_frame("2020", "2021", "quarter")
  expect_identical(time_index(tq, c("2020-03-31", "2020-04-01", "2021-12-31")), c(1L, 2L, 8L))
  ty <- time_frame(2001, 2005, "year")
  expect_identical(time_index(ty, c("2001-06-01", "2005-12-31", "2006-01-01")), c(1L, 5L, NA))
})

test_that("a plain index frame", {
  ti <- time_frame(T = 12)
  expect_equal(ti$T, 12)
  expect_identical(ti$step, "index")
  expect_error(time_index(ti, Sys.Date()), "no dates")
  expect_error(time_frame(T = 0), "positive")
  expect_error(time_frame(2000, T = 3), "either")
  expect_error(time_frame("2020-01"), "give `start` and `end`")
  expect_error(time_frame("not a date", "2020"), "not a date")
})
