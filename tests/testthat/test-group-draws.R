# Checks of predict()'s group draws, groupings given at prediction time and
# kept time steps:
#   1. The group draws' means and quantiles equal the group summaries, for
#      a grouping nested in the blocks and for one assembled from partial
#      sums, in this session and in parallel with many blocks.
#   2. Draw d is the same posterior draw everywhere: with equal weights, a
#      group's draw is the mean of its areas' kept draws.
#   3. A grouping given at prediction time with the same labels as one of
#      the frame's gives the same results; new groupings can be named by id;
#      bad ones are refused.
#   4. keep_times restricts the kept draws, not the summaries.

gd_frame <- function() {
  k <- 5; m <- k^2
  r <- rep(seq_len(k), each = k); cc <- rep(seq_len(k), times = k)
  W <- 1 * (abs(outer(r, r, "-")) + abs(outer(cc, cc, "-")) == 1)
  areal_frame(adjacency = W, ids = sprintf("u%02d", seq_len(m)),
              groups = data.frame(row = paste0("R", r), half = ifelse(r <= 2, "top", "bottom"),
                                  col = paste0("C", cc)), verbose = FALSE)
}
fr_gd <- gd_frame()
sim <- simulate_areal(fr_gd, time_frame(n_steps = 6), n_cells = 40, k = 30, seed = 4)
fit <- fit_areal(sim, control = list(n_iter = 60, n_burn = 20, n_store = 15, n_chains = 2, ess_sweeps = 2,
                                     theta_moves = 1, n_cores = 1, verbose = FALSE))

draws_vs_table <- function(dr, tab, stat = "mean") {
  f <- switch(stat, mean = mean, q50 = function(v) stats::quantile(v, 0.5, names = FALSE, type = 7))
  m <- apply(dr, c(1, 2), f)
  max(abs(m[cbind(match(tab$id, rownames(m)), match(tab$time, colnames(m)))] - tab[[stat]]))
}

test_that("group draws agree with the group summaries, nested or not, in any blocks", {
  for (setting in list(list(n_cores = 1, block_bytes = 5e8), list(n_cores = 2, block_bytes = 3000))) {
    pr <- do.call(predict, c(list(fit, groups = c("half", "col", "row"), group_draws = TRUE, n_draws = 20,
                                  verbose = FALSE), setting))
    expect_identical(names(pr$group_draws), c("half", "col", "row"))
    for (g in names(pr$group_draws)) {
      expect_identical(dim(pr$group_draws[[g]]), c(length(unique(fr_gd$groups[[g]])), 6L, 20L))
      expect_lt(draws_vs_table(pr$group_draws[[g]], pr$groups[[g]]), 1e-10)
      expect_lt(draws_vs_table(pr$group_draws[[g]], pr$groups[[g]], "q50"), 1e-10)
    }
  }
})

test_that("group draws agree with the summaries under weights that vary by area and time", {
  set.seed(3)
  W <- matrix(stats::runif(25 * 6, 0.1, 5), 25, 6)          # (equal weights would hide a wrong denominator)
  for (setting in list(list(n_cores = 1, block_bytes = 5e8), list(n_cores = 2, block_bytes = 3000))) {
    pr <- do.call(predict, c(list(fit, groups = c("half", "col"), group_draws = TRUE, weights = W, n_draws = 20,
                                  keep_times = c(2, 3, 6), verbose = FALSE), setting))
    for (g in c("half", "col"))
      expect_lt(draws_vs_table(pr$group_draws[[g]], pr$groups[[g]][pr$groups[[g]]$t %in% c(2, 3, 6), ]), 1e-10)
  }
})

test_that("a group draw is the mean of its areas' draws (the same posterior draw)", {
  top <- fr_gd$ids[fr_gd$groups$half == "top"]
  pr <- predict(fit, groups = "half", group_draws = "half", keep_draws = top, n_draws = 10, verbose = FALSE)
  expect_lt(max(abs(apply(pr$draws, c(2, 3), mean) - pr$group_draws$half["top", , ])), 1e-10)
})

test_that("groupings given at prediction time work like the frame's", {
  a <- predict(fit, groups = "half", group_draws = TRUE, n_draws = 10, verbose = FALSE)
  b <- predict(fit, groups = list(half2 = fr_gd$groups$half), group_draws = TRUE, n_draws = 10, verbose = FALSE)
  expect_equal(b$groups$half2, a$groups$half)
  expect_equal(unname(b$group_draws$half2), unname(a$group_draws$half))
  expect_identical(b$frame$groups$half2, fr_gd$groups$half)
  expect_identical(b$settings$new_groupings, "half2")
  # named by id, in another order; mixed with a frame grouping
  lab <- stats::setNames(fr_gd$groups$half, fr_gd$ids)[rev(fr_gd$ids)]
  cc <- predict(fit, groups = list("row", h = lab), n_draws = 10, verbose = FALSE)
  expect_equal(cc$groups$h, a$groups$half)
  expect_identical(names(cc$groups), c("row", "h"))
  # refused
  expect_error(predict(fit, groups = list(h = 1:3), verbose = FALSE), "one label per area")
  expect_error(predict(fit, groups = list(half = rep("x", 25)), verbose = FALSE), "already a grouping")
  expect_error(predict(fit, groups = list(mean = rep("x", 25)), verbose = FALSE), "used by the package")
  expect_error(predict(fit, groups = "nope", verbose = FALSE), "not among the frame's groupings")
  expect_error(predict(fit, groups = "half", group_draws = "row", verbose = FALSE), "must name groupings")
})

test_that("keep_times restricts the kept draws, not the summaries", {
  pr <- predict(fit, groups = "half", group_draws = "half", keep_draws = "u01", keep_times = c("2", "5"),
                n_draws = 10, verbose = FALSE)
  expect_identical(dimnames(pr$group_draws$half)[[2]], c("2", "5"))
  expect_identical(dim(pr$draws), c(1L, 2L, 10L))
  expect_identical(nrow(pr$groups$half), 12L)
  expect_lt(draws_vs_table(pr$draws, pr$area[pr$area$id == "u01" & pr$area$t %in% c(2, 5), ]), 1e-10)
  expect_output(print(pr), "aggregate draws kept for half, at 2 time step")
  expect_error(predict(fit, keep_times = 9, verbose = FALSE), "time label")
  expect_false(is.null(pr$version$stareal))
})
