# Checks of predict() (two-phase prediction, block-wise aggregation), its
# inputs and its methods. Ported from Africa_areal/tests/test_predict.R
# check 4, plus non-nested groupings, the choice of blocks, kept draws,
# weights and the output methods:
#   1. predict() gives the same area and group summaries as brute force:
#      whole surfaces per draw (predict_surface), weighted means per group,
#      then summaries. The groupings are rows (nested in halves), halves
#      (the blocks) and columns (not nested: partial sums). Run in parallel
#      and in this session, with two blocks and with many forced by memory.
#   2. Kept draws equal the brute-force surfaces.
#   3. Weights given as a vector, a matrix and a long data frame agree.
#   4. Blocks: a grouping is used only if its largest group fits; blocks are
#      contiguous runs.
#   5. Inputs are checked before any work, and the methods work.

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
  sim <- simulate_areal(pred_frame(), time_frame(n_steps = 6), n_cells = 40, k = 30, seed = 4)
  fit_areal(sim, control = list(n_iter = 60, n_burn = 20, n_store = 15, n_chains = 2, ess_sweeps = 2,
                                theta_moves = 1, n_cores = 1, verbose = FALSE))
}

brute_force <- function(fit, idx, seed, W, probs, exceedance) {
  fr <- fit$data$frame; T <- fit$data$time$n_steps
  st <- fit_setup_of(fit)
  ctx <- list(eig = fr$eig, obs = st$obs, b0 = st$b0, B = st$B, T = T)
  pool <- function(a) matrix(a, prod(dim(a)[1:2]), dim(a)[3])
  th <- pool(fit$draws$theta); xo <- pool(fit$draws$x)
  seeds <- with_seed(seed, sample.int(.Machine$integer.max - 1L, length(idx)))
  surf <- lapply(seq_along(idx), function(d) plogis(predict_surface(xo[idx[d], ], th[idx[d], ], seeds[d], ctx)))
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
  args <- list(fit, groups = c("row", "half", "col"), weights = Wt, n_draws = 20,
               probs = c(0.975, 0.025, 0.5), exceedance = c(0.5, 0.2), keep_draws = c("u03", "u17", "u03"),
               seed = 100, verbose = FALSE)
  # two blocks forced by memory (the surface is 24 kB), one at a time in
  # this session: blocks are the halves, rows are summarised inside them,
  # columns come from partial sums
  pr <- do.call(predict, c(args, list(n_cores = 1, block_bytes = 15000)))
  expect_identical(pr$settings$n_blocks, 2L)
  expect_identical(pr$settings$block_grouping, "half")
  expect_identical(pr$settings$probs, c(0.025, 0.5, 0.975))      # sorted
  idx <- pr$settings$draws_used
  expect_length(idx, 20)
  bf <- brute_force(fit, idx, 100, Wt, c(0.025, 0.5, 0.975), c(0.2, 0.5))
  expect_lt(max_err(pr$area, bf$area), 1e-10)
  for (g in c("row", "half", "col")) {
    expect_identical(pr$groups[[g]]$id, bf[[g]]$id)
    expect_lt(max_err(pr$groups[[g]], bf[[g]]), 1e-10)
  }
  expect_identical(pr$area$id, rep(fr$ids, each = 6))
  expect_null(pr$frame$eig)
  # observed data carried onto the area table
  cells <- fit$data$cells
  at <- (match(cells$id, fr$ids) - 1) * 6 + cells$t
  expect_equal(pr$area$k[at], cells$k)
  expect_equal(sum(pr$area$k), sum(cells$k))
  # many blocks: more than there are halves, so contiguous runs of areas,
  # and every grouping from partial sums
  small <- do.call(predict, c(args, list(n_cores = 1, block_bytes = 2000)))
  expect_null(small$settings$block_grouping)
  expect_gt(small$settings$n_blocks, 2)
  expect_lt(max_err(small$area, pr$area), 1e-12)
  for (g in c("row", "half", "col")) expect_lt(max_err(small$groups[[g]], pr$groups[[g]]), 1e-12)
  # kept draws (a repeated id is kept once)
  expect_equal(dim(pr$draws), c(2, 6, 20))
  for (a in c("u03", "u17")) {
    i <- match(a, fr$ids)
    expect_lt(max(abs(pr$draws[a, , ] - sapply(bf$surf, function(P) P[i, ]))), 1e-12)
  }
  # in parallel (two workers): the same, up to BLAS rounding
  skip_on_cran()
  par_pr <- do.call(predict, c(args, list(n_cores = 2, block_bytes = 15000)))
  expect_lt(max_err(par_pr$area, pr$area), 1e-8)
  for (g in c("row", "half", "col")) expect_lt(max_err(par_pr$groups[[g]], pr$groups[[g]]), 1e-8)
})

test_that("predict() checks its inputs first, handles one draw, and leaves the stream untouched", {
  fit <- pred_fit()
  set.seed(5); before <- runif(3)
  set.seed(5); pr <- predict(fit, n_draws = 5, verbose = FALSE); after <- runif(3)
  expect_identical(after, before)
  expect_s3_class(pr, "stareal_pred")
  expect_output(print(pr), "stareal_pred")
  one <- predict(fit, n_draws = 1, groups = "col", verbose = FALSE)
  expect_equal(one$settings$n_draws, 1)
  expect_true(all(one$area$sd == 0 | is.na(one$area$sd)))
  expect_error(predict(fit, groups = "nope"), "nope")
  expect_error(predict(fit, n_draws = 1000), "n_draws")
  expect_error(predict(fit, n_draws = 2.5), "n_draws")
  expect_error(predict(fit, keep_draws = "zz"), "zz")
  expect_error(predict(fit, exceedance = 2), "exceedance")
  expect_error(predict(fit, probs = numeric(0)), "probs")
  expect_error(predict(fit, seed = NULL), "seed")
  expect_error(predict(fit, block_bytes = 0), "block_bytes")
  expect_error(predict(fit, levels = "row"), "unused argument")
})

test_that("weights as a vector, a matrix and a data frame agree", {
  fr <- pred_frame(); tf <- time_frame(n_steps = 3)
  v <- setNames(seq_len(fr$m), rev(fr$ids))                 # named, in reverse order
  Wv <- weights_matrix(v, fr, tf)
  expect_equal(unname(Wv[, 1]), rev(seq_len(fr$m)))
  Wm <- weights_matrix(matrix(rev(seq_len(fr$m)), fr$m, 3), fr, tf)
  expect_identical(Wm, Wv)
  long <- data.frame(id = rep(fr$ids, 3), t = rep(1:3, each = fr$m), weight = rep(rev(seq_len(fr$m)), 3))
  expect_identical(weights_matrix(long[sample(nrow(long)), ], fr, tf), Wv)
  long_lab <- data.frame(id = long$id, time = as.character(long$t), weight = long$weight)
  expect_identical(weights_matrix(long_lab, fr, tf), Wv)
  expect_identical(weights_matrix(data.frame(id = names(v), weight = unname(v)), fr, tf), Wv)
  expect_error(weights_matrix(long[-1, ], fr, tf), "no value")
  expect_error(weights_matrix(rbind(long, long[1, ]), fr, tf), "more than one row")
  expect_error(weights_matrix(-Wv, fr, tf), ">= 0")
  expect_error(weights_matrix(1:3, fr, tf), "one value per area")
  # matrix columns named by time label are matched, not taken in order
  M <- matrix(1:(fr$m * 3), fr$m, 3, dimnames = list(NULL, c("3", "1", "2")))
  expect_equal(weights_matrix(M, fr, tf)[, 1], M[, 2], ignore_attr = TRUE)
})

test_that("blocks: a grouping is used only if its largest group fits, and blocks are contiguous", {
  g <- pred_frame()$groups
  expect_identical(block_grouping(g, c("row", "half", "col"), 2), "half")
  expect_identical(block_grouping(g, c("row", "col"), 2), "row")      # a tie (5 groups each): the first
  expect_null(block_grouping(g, "half", 4))                           # too few groups for 4 blocks
  expect_identical(block_grouping(g, c("row", "half"), 2, bytes_per_area = 1, block_bytes = 12), "row")  # half too big
  expect_true(nests_in(g$row, g$half))
  expect_false(nests_in(g$col, g$half))
  b <- make_blocks(g, c("row", "half"), "half", 2)
  expect_length(b, 2)
  expect_true(all(vapply(b, function(i) length(unique(g$half[i])), 1) == 1))        # whole halves
  b5 <- make_blocks(g, "col", NULL, 5)
  expect_true(all(vapply(b5, function(i) length(unique(g$col[i])), 1) == 1))        # contiguous: one column each
  expect_setequal(unlist(b5), seq_len(nrow(g)))
})

test_that("as.data.frame, st_as_sf and plot", {
  skip_if_not_installed("sf"); skip_if_not_installed("ggplot2")
  sim <- simulate_areal(areal_lattice(4, 4, blocks = 2), time_frame(2001, 2005, "year"), n_cells = 30, seed = 1)
  fit <- fit_areal(sim, control = list(n_iter = 60, n_burn = 20, n_store = 10, n_chains = 1, ess_sweeps = 2,
                                       verbose = FALSE))
  pr <- predict(fit, groups = "block", verbose = FALSE)
  expect_identical(as.data.frame(pr), pr$area)
  expect_identical(as.data.frame(pr, group = "block"), pr$groups$block)
  expect_error(as.data.frame(pr, group = "region"), "group")
  s <- sf::st_as_sf(pr, time = "2003")
  expect_s3_class(s, "sf"); expect_equal(nrow(s), 16)
  sb <- sf::st_as_sf(pr, time = 3, group = "block")
  expect_equal(nrow(sb), 4)
  expect_s3_class(pr$area$date, "Date")
  expect_s3_class(plot(pr, time = "2003"), "ggplot")
  expect_s3_class(plot(pr, group = "block"), "ggplot")
  expect_s3_class(plot(pr), "ggplot")
  expect_error(plot(pr, time = "1999"), "time label")
  expect_error(plot(pr, levels = "block"), "unused argument")
  one <- fit_areal(simulate_areal(areal_lattice(4, 4), time_frame(n_steps = 1), formula = ~ 1, beta = -1,
                                  n_cells = 10, seed = 2),
                   formula = ~ 1, control = list(n_iter = 60, n_burn = 20, n_store = 10, n_chains = 1,
                                                 ess_sweeps = 2, verbose = FALSE))
  expect_s3_class(plot(predict(one, verbose = FALSE)), "ggplot")       # one time step: point ranges
})
