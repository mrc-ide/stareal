# Checks of predict() (two-phase prediction, block-wise aggregation) and of
# the weights. Ported from Africa_areal/tests/test_predict.R check 4, plus
# non-nested groupings, the choice of blocks, kept draws and weights:
#   1. predict() gives the same area and group summaries as brute force:
#      whole surfaces per draw (predict_surface), weighted means per group,
#      then summaries. The groupings are rows (nested in halves), halves
#      (the blocks) and columns (not nested: partial sums). Run in parallel
#      and one block at a time, with several blocks forced by memory.
#   2. Kept draws equal the brute-force surfaces.
#   3. Weights given as a vector, a matrix and a long data frame agree.

# a 5 x 5 lattice with three groupings
pred_frame <- function() {
  k <- 5; m <- k^2
  r <- rep(seq_len(k), each = k); cc <- rep(seq_len(k), times = k)
  W <- 1 * (abs(outer(r, r, "-")) + abs(outer(cc, cc, "-")) == 1)
  areal_frame(adjacency = W, ids = sprintf("u%02d", seq_len(m)),
              groups = data.frame(row = paste0("R", r), half = ifelse(r <= 2, "top", "bottom"),
                                  col = paste0("C", cc)), verbose = FALSE)
}

pred_fit <- function() {
  fr <- pred_frame(); tf <- time_frame(T = 6)
  sim <- simulate_areal(fr, tf, n_cells = 40, k = 30, seed = 4)
  d <- areal_data(sim$data, fr, tf, y = "positive", k = "tested", area = "id", t = "t", verbose = FALSE)
  fit_areal(d, control = mcmc_control(n_iter = 60, n_burn = 20, n_store = 15, n_chains = 2, ess_sweeps = 2,
                                      theta_moves = 1, parallel = FALSE, progress = FALSE))
}

brute_force <- function(fit, idx, seed, W, probs, exceedance) {
  fr <- fit$data$frame; T <- fit$data$time$T
  ctx <- list(eig = fr$eig, obs = fit$setup$obs, b0 = fit$setup$b0, B = fit$setup$B, T = T)
  pool <- function(a) matrix(a, prod(dim(a)[1:2]), dim(a)[3])
  th <- pool(fit$draws$theta); xo <- pool(fit$draws$x)
  surf <- lapply(seq_along(idx), function(d) plogis(predict_surface(xo[idx[d], ], th[idx[d], ], seed + d, ctx)))
  area <- summarise_draw_matrix(sapply(surf, function(P) as.vector(t(P))), probs, exceedance)  # area, then t
  group <- function(gv) {
    lv <- sort(unique(gv)); G <- as.matrix(group_matrix(gv, lv))
    A <- sapply(surf, function(P) as.vector(t((G %*% (W * P)) / (G %*% W))))                 # group, then t
    cbind(id = rep(lv, each = T), summarise_draw_matrix(A, probs, exceedance))
  }
  list(surf = surf, area = area, row = group(fr$groups$row), half = group(fr$groups$half),
       col = group(fr$groups$col))
}

stat_cols <- c("mean", "sd", "q2.5", "q50", "q97.5", "p_gt_0.2", "p_gt_0.5")
max_err <- function(a, b) max(abs(as.matrix(a[stat_cols]) - as.matrix(b[stat_cols])))

test_that("predict() equals brute force, for nested and non-nested groupings", {
  fit <- pred_fit()
  fr <- fit$data$frame
  set.seed(31)
  Wt <- matrix(runif(fr$m * 6, 0.5, 2), fr$m, 6, dimnames = list(fr$ids, NULL))
  args <- list(fit, levels = c("row", "half", "col"), weights = Wt, n_draws = 20,
               exceedance = c(0.2, 0.5), keep_draws = c("u03", "u17"), seed = 100, progress = FALSE)
  # two blocks forced by memory (the surface is 24 kB), one at a time in
  # this session: blocks are the halves, rows are summarised inside them,
  # columns come from partial sums
  seq_pr <- do.call(predict, c(args, list(parallel = FALSE, block_bytes = 15000)))
  expect_identical(seq_pr$settings$n_blocks, 2L)
  expect_identical(seq_pr$settings$block_grouping, "half")
  idx <- seq_pr$settings$draws_used
  expect_length(idx, 20)
  bf <- brute_force(fit, idx, 100, Wt, c(0.025, 0.5, 0.975), c(0.2, 0.5))
  expect_lt(max_err(seq_pr$area, bf$area), 1e-10)
  for (lv in c("row", "half", "col")) {
    expect_identical(seq_pr[[lv]]$id, bf[[lv]]$id)
    expect_lt(max_err(seq_pr[[lv]], bf[[lv]]), 1e-10)
  }
  expect_identical(seq_pr$area$id, rep(fr$ids, each = 6))
  # observed data carried onto the area table
  cells <- fit$data$cells
  at <- (match(cells$id, fr$ids) - 1) * 6 + cells$t
  expect_identical(seq_pr$area$k[at], cells$k)
  expect_equal(sum(seq_pr$area$k), sum(cells$k))
  # many blocks: more than there are halves, so single areas, and every
  # grouping from partial sums
  small <- do.call(predict, c(args, list(parallel = FALSE, block_bytes = 2000)))
  expect_null(small$settings$block_grouping)
  expect_gt(small$settings$n_blocks, 2)
  for (lv in c("area", "row", "half", "col")) expect_lt(max_err(small[[lv]], seq_pr[[lv]]), 1e-12)
  expect_lt(max(abs(small$draws - seq_pr$draws)), 1e-12)   # rows computed in other groupings: rounding only
  # kept draws
  for (a in c("u03", "u17")) {
    i <- match(a, fr$ids)
    expect_lt(max(abs(seq_pr$draws[a, , ] - sapply(bf$surf, function(P) P[i, ]))), 1e-12)
  }
  # in parallel (two workers): the same, up to BLAS rounding
  skip_on_cran()
  par_pr <- do.call(predict, c(args, list(parallel = TRUE, n_cores = 2)))
  for (lv in c("area", "row", "half", "col")) expect_lt(max_err(par_pr[[lv]], seq_pr[[lv]]), 1e-8)
})

test_that("predict() leaves the user's random-number stream untouched and checks its inputs", {
  fit <- pred_fit()
  set.seed(5); before <- runif(3)
  set.seed(5); pr <- predict(fit, n_draws = 5, parallel = FALSE, progress = FALSE); after <- runif(3)
  expect_identical(after, before)
  expect_s3_class(pr, "stareal_pred")
  expect_output(print(pr), "stareal_pred")
  expect_error(predict(fit, levels = "nope", parallel = FALSE), "nope")
  expect_error(predict(fit, n_draws = 1000, parallel = FALSE), "n_draws")
  expect_error(predict(fit, keep_draws = "zz", parallel = FALSE), "zz")
  expect_error(predict(fit, exceedance = 2, parallel = FALSE), "exceedance")
})

test_that("weights as a vector, a matrix and a data frame agree", {
  fr <- pred_frame(); tf <- time_frame(T = 3)
  v <- setNames(seq_len(fr$m), rev(fr$ids))                 # named, in reverse order
  Wv <- weights_matrix(v, fr, tf)
  expect_equal(unname(Wv[, 1]), rev(seq_len(fr$m)))
  Wm <- weights_matrix(matrix(rev(seq_len(fr$m)), fr$m, 3), fr, tf)
  expect_identical(Wm, Wv)
  long <- data.frame(id = rep(fr$ids, 3), t = rep(1:3, each = fr$m), weight = rep(rev(seq_len(fr$m)), 3))
  expect_identical(weights_matrix(long[sample(nrow(long)), ], fr, tf), Wv)
  expect_identical(weights_matrix(data.frame(id = names(v), weight = unname(v)), fr, tf), Wv)
  expect_error(weights_matrix(long[-1, ], fr, tf), "no value")
  expect_error(weights_matrix(-Wv, fr, tf), ">= 0")
  expect_error(weights_matrix(1:3, fr, tf), "one value per area")
})

test_that("block_grouping prefers the grouping that nests the most levels", {
  g <- pred_frame()$groups
  expect_identical(block_grouping(g, c("row", "half", "col"), 2), "half")
  expect_identical(block_grouping(g, c("row", "col"), 2), "row")      # a tie (5 groups each): the first
  expect_null(block_grouping(g, "half", 4))                           # too few groups for 4 blocks
  expect_true(nests_in(g$row, g$half))
  expect_false(nests_in(g$col, g$half))
})
