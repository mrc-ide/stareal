# predict.R
#
# predict() for a fit: the full area x time surface, one surface per
# posterior draw by Matheron's rule (predict_core.R), reduced to summaries
# per area and time step; and aggregation of each joint draw to the frame's
# groupings before summarising, never from area-level summaries, which
# would ignore the correlations.
#
# Ported from Africa_areal/R/predict.R (run_prediction, predict_block,
# summarise_draw_matrix, make_blocks, group_matrix), 2026-10-12. Changes:
# - any number of named groupings instead of admin1 and admin0. Blocks are
#   whole groups of one grouping (block_grouping(), below); a grouping
#   nested in it is complete inside each block and is summarised there, as
#   admin1 was; any other is assembled from per-block partial sums, as
#   admin0 was, added up as each block's result arrives;
# - blocks are contiguous runs of areas ordered by the groupings, so a
#   block holds few groups of any grouping, and the partial sums stay small;
# - the number of blocks bounds memory (block_bytes per block), and a
#   grouping whose largest group would not fit is not used for blocks;
# - weights from a matrix, a vector or a data frame (weights_matrix());
# - quantiles named q2.5, q50, ..., exceedance p_gt_<threshold>; no
#   rounding; areas in cell_id order (area, then time step);
# - optionally keeps the draws of chosen areas;
# - draws are selected evenly from the stored draws of all chains, and
#   each draw's seed is drawn from the user's seed.
#
# How it scales (Africa ADMIN2: m = 4,319, T up to 600, 1,000 draws):
# holding all draws would need ~70 GB, so the work is split in two phases.
#   Phase A (parallel over draws): the conditioning step, which needs only
#     the observed areas. Per draw it keeps w (length n), beta* and the seed.
#   Phase B (parallel over blocks of areas): for each draw, regenerate the
#     same Z from the seed and compute the block's rows of the surface; keep
#     the block's draws (rows x T x draws), summarise them, and aggregate
#     within the block. Each block regenerates the m x T noise of every
#     draw, so fewer, larger blocks are faster.
#
# ------------------------------------------------------------------

#' Predict the full area-time surface from a fit
#'
#' Reconstructs the latent prevalence surface for every area and time step,
#' in percent (as all prevalences in *stareal* are),
#' observed or not, once per posterior draw. It does this by conditional
#' simulation (Matheron's rule), which is exact given each draw. It then
#' summarises the draws per area and time step. Each grouping in `groups` is
#' aggregated draw by draw, as a weighted mean of its areas' prevalence,
#' and then summarised. Its uncertainty therefore accounts for the
#' correlation between areas.
#'
#' The surface is the field plus fixed effects, without the cell-level
#' nugget: it is the prevalence of the area-time cell, not of one survey in
#' it. At observed cells it therefore differs from the latent values that
#' [diagnose()] and [crossval()] report, which include the nugget.
#'
#' @param object A [fit_areal()] object.
#' @param groups Names of groupings of the frame (see [areal_frame()]) to
#'   aggregate to. Areas are always summarised.
#' @param weights Weights of the areas in group aggregates, for example
#'   population or number of infections. The aggregate of a group is
#'   \eqn{\sum_a w_{at} p_{at} / \sum_a w_{at}} for each draw. Give one of:
#'   `NULL` (equal weights); a vector with one value per area (named by area
#'   id, or in frame order); an areas x time steps matrix (rows named by
#'   area id or in frame order; columns named by time label or in order);
#'   or a data frame with columns `id`, `weight` and optionally `time`
#'   (label) or `t` (step number), with every area-time combination present
#'   once.
#' @param n_draws Number of posterior draws to use, taken evenly from the
#'   stored draws of all chains. `NULL`: all of them.
#' @param probs Quantiles to report.
#' @param exceedance Prevalence thresholds, in percent (e.g. `5` for 5%);
#'   for each, the posterior probability that prevalence exceeds it is
#'   reported, in a column `p_gt_<threshold>` (e.g. `p_gt_5`).
#' @param keep_draws Ids of areas whose prevalence draws are returned in full.
#' @param seed Random seed, a whole number. Each draw gets its own seed,
#'   drawn from this one. The caller's random-number stream is left
#'   untouched.
#' @param n_cores Number of worker processes. Default: 1 (this session) if
#'   the whole surface fits in one block of `block_bytes`, otherwise the
#'   cores less one. Every worker holds its own copy of the frame's
#'   eigenvectors (8 m^2 bytes).
#' @param blas_threads BLAS threads per worker.
#' @param block_bytes Approximate memory for one block's draws (bytes).
#'   Areas are processed in blocks of at most this size (up to `n_cores`
#'   at once, each also needing about half as much again per exceedance
#'   threshold). Each block regenerates every draw's noise, so fewer,
#'   larger blocks are faster.
#' @param verbose Report progress: the draws conditioned and predicted so
#'   far, with an estimate of the time left.
#' @param ... Not used.
#' @return An object of class `stareal_pred`:
#'   \describe{
#'     \item{area}{one row per area and time step, in the order of area, then
#'       time step: `id`, the frame's groupings, `t`, `time`, `date` (start
#'       of the step; calendar time frames), the pooled data `y` and `k` (0
#'       where there are none), then `mean`, `sd` and one column per
#'       quantile (`q2.5`, `q50`, ...) of prevalence, in percent (0 to 100),
#'       and one column per threshold (`p_gt_5`, ...), the probability
#'       (0 to 1) that prevalence exceeds it}
#'     \item{groups}{a list with one table per grouping in `groups`: `id`
#'       (the group), `t`, `time`, `date` and the same statistics, ordered
#'       by group, then time step}
#'     \item{draws}{if `keep_draws`: prevalence draws, in percent, an array
#'       (area x time step x draw)}
#'     \item{frame, time, settings}{for plotting and the record (the frame
#'       without its eigenvectors)}
#'   }
#'   `as.data.frame()` gives the area table (or a group's, with `group =`);
#'   `sf::st_as_sf(pr, time = ...)` gives it as polygons for one time step.
#' @examples
#' fr <- areal_lattice(4, 4, blocks = c(2, 2))
#' tf <- time_frame(n_steps = 5)
#' sim <- simulate_areal(fr, tf, n_cells = 30, seed = 1)
#' fit <- fit_areal(sim, control = list(n_iter = 100, n_burn = 50, n_store = 25, n_chains = 1,
#'                                      ess_sweeps = 2, verbose = FALSE))
#' pr <- predict(fit, groups = "block", verbose = FALSE)
#' head(pr$area)
#' pr$groups$block
#' @export
predict.stareal_fit <- function(object, groups = NULL, weights = NULL, n_draws = NULL,
                                probs = c(0.025, 0.5, 0.975), exceedance = NULL, keep_draws = NULL,
                                seed = 1, n_cores = NULL, blas_threads = 1, block_bytes = 5e8, verbose = TRUE, ...) {
  check_dots(...)
  fit <- object
  fr <- fit$data$frame; tf <- fit$data$time; m <- fr$m; T <- tf$n_steps
  # ---- every argument checked before any work ----
  gnames <- setdiff(names(fr$groups), "id")
  groups <- unique(groups)
  if (length(setdiff(groups, gnames)))
    stop(sprintf("group(s) not among the frame's groupings (%s): %s",
                 if (length(gnames)) paste(gnames, collapse = ", ") else "none",
                 paste(setdiff(groups, gnames), collapse = ", ")), call. = FALSE)
  if (!is.numeric(probs) || !length(probs) || anyNA(probs) || any(probs <= 0 | probs >= 1))
    stop("`probs` must be one or more probabilities in (0, 1)", call. = FALSE)
  probs <- sort(unique(probs))
  if (!is.null(exceedance) && (!is.numeric(exceedance) || anyNA(exceedance) || any(exceedance <= 0 | exceedance >= 100)))
    stop("`exceedance` must be prevalences in percent, between 0 and 100 (e.g. 5 for 5%)", call. = FALSE)
  exceedance <- sort(unique(exceedance))
  keep <- if (is.null(keep_draws)) integer(0) else unique(match(id_chr(keep_draws), fr$ids))
  if (anyNA(keep)) stop(sprintf("`keep_draws` ids not in the frame: %s",
                                paste(keep_draws[is.na(match(id_chr(keep_draws), fr$ids))], collapse = ", ")),
                        call. = FALSE)
  seed <- check_whole(seed, "seed", 0)
  blas_threads <- check_whole(blas_threads, "blas_threads")
  n_cores <- check_whole(n_cores, "n_cores", null_ok = TRUE)
  if (!is.numeric(block_bytes) || length(block_bytes) != 1 || !is.finite(block_bytes) || block_bytes <= 0)
    stop("`block_bytes` must be one positive number (bytes)", call. = FALSE)
  pool <- function(a) matrix(a, prod(dim(a)[1:2]), dim(a)[3])
  th <- pool(fit$draws$theta); xo <- pool(fit$draws$x)
  D_all <- nrow(xo)
  if (is.null(n_draws)) n_draws <- D_all
  n_draws <- check_whole(n_draws, "n_draws")
  if (n_draws > D_all) stop(sprintf("`n_draws` must be at most the %d stored draws", D_all), call. = FALSE)
  W <- weights_matrix(weights, fr, tf)
  for (g in groups) {
    den <- rowsum(W, fr$groups[[g]])
    if (any(den <= 0)) stop(sprintf("some groups of `%s` have zero total weight in some time steps", g),
                            call. = FALSE)
  }

  idx <- unique(round(seq(1, D_all, length.out = n_draws)))
  D <- length(idx)
  seeds <- with_seed(seed, sample.int(.Machine$integer.max - 1L, D))
  bytes <- 8 * m * T * D
  n_blocks_mem <- ceiling(bytes / block_bytes)
  if (is.null(n_cores)) n_cores <- if (n_blocks_mem <= 1) 1L else default_cores(32)
  ctx <- list(eig = fr$eig, obs = make_obs(fit$data$cells$area, fit$data$cells$t, T,
                                           X = design_eval(fit$design, fit$data$cells$area, fit$data$cells$t)),
              b0 = fit$priors$beta_mean, B = diag(fit$priors$beta_sd^2, length(fit$priors$beta_sd)), T = T,
              design = design_rows(fit$design))
  job_fail <- function(what) function(e) stop(sprintf("prediction %s failed (job %d): %s", what, e$job,
                                                      conditionMessage(e)), call. = FALSE)

  # ---- phase A: conditioning, in chunks of draws ----
  say(verbose, "phase A: conditioning %d draws on %d observed cells", D, ctx$obs$n)
  chunks <- parallel::splitIndices(D, min(D, max(10L, 4L * n_cores)))
  t0 <- proc.time()[["elapsed"]]; done <- 0
  progress_A <- function(i, v) {                       # in the session, as each chunk of draws returns
    done <<- done + length(v)
    el <- proc.time()[["elapsed"]] - t0
    say(verbose, "phase A: %d of %d draws conditioned (%s left)", done, D, format_seconds(el / done * (D - done)))
    v
  }
  jobA <- function(ds, sh) lapply(ds, function(d)
    utils::getFromNamespace("prepare_draw", "stareal")(sh$x[d, ], sh$theta[d, ], sh$seeds[d], sh$ctx))
  prep <- tryCatch(run_jobs(chunks, jobA, shared = list(x = xo[idx, , drop = FALSE], theta = th[idx, , drop = FALSE],
                                                       seeds = seeds, ctx = ctx),
                            n_workers = n_cores, blas_threads = blas_threads, verbose = verbose,
                            collect = progress_A),
                   stareal_job_error = job_fail("phase A"))
  prep <- unlist(prep, recursive = FALSE)

  # ---- phase B: blocks of areas ----
  n_blocks <- max(min(n_cores, m), n_blocks_mem)
  bg <- block_grouping(fr$groups, groups, n_blocks, 8 * T * D, block_bytes)
  blocks <- make_blocks(fr$groups, groups, bg, n_blocks)
  nested <- vapply(groups, function(g) !is.null(bg) && nests_in(fr$groups[[g]], fr$groups[[bg]]), logical(1))
  say(verbose, "phase B: %d block%s of areas (%d-%d areas each)%s", length(blocks), if (length(blocks) > 1) "s" else "",
      min(lengths(blocks)), max(lengths(blocks)), if (is.null(bg)) "" else sprintf(", whole groups of `%s`", bg))
  stats <- list(probs = probs, exceedance = exceedance)
  gtab <- fr$groups[groups]
  # partial sums of the groupings not nested in the blocks, added up as each
  # block arrives
  acc <- list()
  for (g in groups[!nested]) {
    lvl <- sort(unique(gtab[[g]]))
    acc[[g]] <- list(id = lvl, num = matrix(0, length(lvl) * T, D), den = matrix(0, length(lvl), T))
  }
  collect <- function(i, r) {
    for (g in names(acc)) {
      ps <- r$groups[[g]]
      G <- length(acc[[g]]$id)
      i0 <- match(ps$id, acc[[g]]$id)
      rows <- as.vector(outer(i0, (seq_len(T) - 1) * G, `+`))       # (group, t) at (t - 1) * G + i
      acc[[g]]$num[rows, ] <<- acc[[g]]$num[rows, ] + ps$num
      acc[[g]]$den[i0, ] <<- acc[[g]]$den[i0, ] + ps$den
      r$groups[[g]] <- NULL
    }
    r
  }
  jobB <- function(b, sh) utils::getFromNamespace("predict_block", "stareal")(b, sh)
  res <- tryCatch(run_jobs(seq_along(blocks), jobB,
                           shared = list(ctx = c(ctx, list(prep = prep)), blocks = blocks, W = W, groups = gtab,
                                         nested = nested, stats = stats, keep = keep, verbose = verbose),
                           n_workers = n_cores, blas_threads = blas_threads, verbose = verbose, collect = collect),
                  stareal_job_error = job_fail("phase B"))

  # ---- assemble ----
  label <- function(df, ...) {
    out <- data.frame(..., t = df$t, time = tf$labels[df$t], stringsAsFactors = FALSE, check.names = FALSE)
    if (!is.null(tf$starts)) out$date <- tf$starts[df$t]
    out
  }
  out <- list()
  u <- do.call(rbind, lapply(res, `[[`, "area"))
  u <- u[order(u$area, u$t), ]
  cells <- fit$data$cells
  yk <- match((u$area - 1L) * T + u$t, cells$cell_id)
  lab <- label(u, id = fr$ids[u$area], fr$groups[u$area, gnames, drop = FALSE])
  lab$y <- ifelse(is.na(yk), 0, cells$y[yk]); lab$k <- ifelse(is.na(yk), 0, cells$k[yk])
  out$area <- cbind(lab, u[setdiff(names(u), c("area", "t"))])
  rownames(out$area) <- NULL
  out$groups <- list()
  for (g in groups) {
    if (nested[[g]]) {
      tb <- do.call(rbind, lapply(res, function(r) r$groups[[g]]))
    } else {
      a <- acc[[g]]; G <- length(a$id)
      tb <- cbind(data.frame(id = rep(a$id, times = T), t = rep(seq_len(T), each = G), stringsAsFactors = FALSE),
                  summarise_draw_matrix(a$num / as.vector(a$den), probs, exceedance))
    }
    tb <- tb[order(tb$id, tb$t), ]
    out$groups[[g]] <- cbind(label(tb, id = tb$id), tb[setdiff(names(tb), c("id", "t"))])
    rownames(out$groups[[g]]) <- NULL
  }
  if (length(keep)) {
    dr <- array(NA_real_, c(length(keep), T, D), dimnames = list(fr$ids[keep], tf$labels, NULL))
    for (r in res) if (length(r$draws)) for (a in names(r$draws)) dr[a, , ] <- r$draws[[a]]
    out$draws <- dr
  }
  fr$eig <- NULL                                       # 8 m^2 bytes that nothing downstream uses
  out$frame <- fr; out$time <- tf
  out$settings <- list(groups = groups, n_draws = D, draws_used = idx, probs = probs, exceedance = exceedance,
                       seed = seed, weights = if (is.null(weights)) "equal" else "given", block_grouping = bg,
                       n_blocks = length(blocks))
  structure(out, class = "stareal_pred")
}

# ---- phase B for one block ----------------------------------------------------

# sh: ctx (with prep, the prepare_draw() outputs), blocks, W (m x T weights,
# rows named by id), groups (data frame of the groupings' labels, frame
# order), nested (per grouping: complete inside every block), stats, keep
# (areas whose draws are returned), verbose (report progress about ten times
# per block, as messages, which reach the session from workers too).
# Prevalence is in percent from here on. Returns the block's area summaries; per
# grouping either its summaries (nested) or partial sums (num per group,
# time and draw; den per group and time); and the kept draws.
predict_block <- function(b, sh) {
  ctx <- sh$ctx; rows <- sh$blocks[[b]]; nr <- length(rows); T <- ctx$T; D <- length(ctx$prep)
  P <- matrix(NA_real_, nr * T, D)                 # cell (area i, step t) at (t - 1) * nr + i
  parts <- rows_parts(rows, ctx)
  t0 <- proc.time()[["elapsed"]]; every <- max(1L, ceiling(D / 10))
  for (d in seq_len(D)) {
    P[, d] <- 100 * stats::plogis(as.vector(predict_rows(rows, ctx$prep[[d]], ctx, parts)))
    if (isTRUE(sh$verbose) && (d %% every == 0 || d == D)) {
      el <- proc.time()[["elapsed"]] - t0
      message(sprintf("phase B: block %d of %d, %d of %d draws predicted (%s left)", b, length(sh$blocks), d, D,
                      format_seconds(el / d * (D - d))))
    }
  }
  st <- sh$stats
  out <- list(area = cbind(data.frame(area = rep(rows, times = T), t = rep(seq_len(T), each = nr)),
                           summarise_draw_matrix(P, st$probs, st$exceedance)))
  Wb <- sh$W[rows, , drop = FALSE]                 # nr x T
  out$groups <- list()
  for (g in names(sh$groups)) {
    gv <- sh$groups[[g]][rows]
    lvl <- sort(unique(gv))
    G <- group_matrix(gv, lvl)
    den <- as.matrix(G %*% Wb)                     # groups x T
    num <- matrix(NA_real_, length(lvl) * T, D)
    for (d in seq_len(D)) num[, d] <- as.vector(as.matrix(G %*% (Wb * matrix(P[, d], nr, T))))
    out$groups[[g]] <- if (sh$nested[[g]]) {
      cbind(data.frame(id = rep(lvl, times = T), t = rep(seq_len(T), each = length(lvl)),
                       stringsAsFactors = FALSE),
            summarise_draw_matrix(num / as.vector(den), st$probs, st$exceedance))
    } else list(id = lvl, num = num, den = den)
  }
  kb <- intersect(sh$keep, rows)
  out$draws <- lapply(stats::setNames(kb, rownames(sh$W)[kb]), function(a) {
    i <- match(a, rows)
    P[(seq_len(T) - 1) * nr + i, , drop = FALSE]   # T x D
  })
  out
}

# rows = cells, columns = draws (prevalence, percent) -> mean, sd, quantiles
# and exceedance probabilities (thresholds in percent)
summarise_draw_matrix <- function(P, probs, exceedance) {
  q <- matrixStats::rowQuantiles(P, probs = probs, drop = FALSE)
  colnames(q) <- paste0("q", 100 * probs)
  out <- data.frame(mean = rowMeans(P), sd = matrixStats::rowSds(P), q, check.names = FALSE)
  for (e in exceedance) out[[paste0("p_gt_", e)]] <- rowMeans(P > e)
  out
}

# ---- blocks and groupings ----------------------------------------------------

# is every group of `g` inside one group of `h`?
nests_in <- function(g, h) all(tapply(h, g, function(v) length(unique(v))) == 1)

# The grouping whose whole groups make the blocks: among the requested
# groupings with at least n_blocks groups and whose largest group fits in
# block_bytes (bytes_per_area each), the one with the most groupings nested
# in it (ties: the one with more groups, for balance). NULL: blocks of
# single areas, and every grouping is assembled from partial sums.
block_grouping <- function(groups, levels, n_blocks, bytes_per_area = 0, block_bytes = Inf) {
  fits <- function(g) length(unique(groups[[g]])) >= n_blocks &&
    max(table(groups[[g]])) * bytes_per_area <= block_bytes
  cand <- levels[vapply(levels, fits, logical(1))]
  if (!length(cand)) return(NULL)
  score <- vapply(cand, function(c) sum(vapply(levels, function(g) nests_in(groups[[g]], groups[[c]]),
                                               logical(1))), numeric(1))
  n_groups <- vapply(cand, function(c) length(unique(groups[[c]])), numeric(1))
  cand[order(-score, -n_groups)][1]
}

# Blocks of areas for phase B: areas ordered by the requested groupings
# (coarsest first) and then frame order, cut into about n_blocks contiguous
# runs of similar size, never splitting a group of `bg`. Contiguity keeps
# the groups of any grouping in few blocks.
make_blocks <- function(gtab, levels, bg, n_blocks) {
  m <- nrow(gtab)
  lev <- levels[order(vapply(levels, function(g) length(unique(gtab[[g]])), numeric(1)))]   # coarsest first
  ord <- if (length(lev)) do.call(order, c(unname(as.list(gtab[lev])), list(seq_len(m)))) else seq_len(m)
  unit <- if (is.null(bg)) seq_len(m) else match(gtab[[bg]], unique(gtab[[bg]][ord]))   # indivisible units
  runs <- split(ord, factor(unit[ord], levels = unique(unit[ord])))                     # in order
  target <- m / n_blocks
  blocks <- list(); cur <- integer(0)
  for (i in seq_along(runs)) {
    cur <- c(cur, runs[[i]])
    if (length(cur) >= target && length(blocks) < n_blocks - 1) { blocks[[length(blocks) + 1]] <- cur; cur <- integer(0) }
  }
  if (length(cur)) blocks[[length(blocks) + 1]] <- cur
  lapply(blocks, sort)
}

# group-membership matrix: one row per group in `levels`, one column per area
group_matrix <- function(group, levels) {
  Matrix::sparseMatrix(i = match(group, levels), j = seq_along(group), x = 1,
                       dims = c(length(levels), length(group)))
}

# ---- weights ------------------------------------------------------------------

# the user's weights as an m x T matrix in frame order, rows named by id
weights_matrix <- function(weights, frame, time) {
  m <- frame$m; T <- time$n_steps; ids <- frame$ids
  by_id <- function(nm) {                          # row order of named rows/entries
    if (is.null(nm)) return(seq_len(m))
    i <- match(ids, id_chr(nm))
    if (anyNA(i)) stop(sprintf("`weights` missing for %d area(s), e.g. %s", sum(is.na(i)), ids[which(is.na(i))[1]]),
                       call. = FALSE)
    i
  }
  if (is.null(weights)) {
    W <- matrix(1, m, T)
  } else if (is.data.frame(weights)) {
    if (!all(c("id", "weight") %in% names(weights)))
      stop("a `weights` data frame needs columns `id` and `weight`", call. = FALSE)
    a <- match(id_chr(weights$id), ids)
    if (anyNA(a)) stop("`weights` has ids not in the frame", call. = FALSE)
    if ("t" %in% names(weights) || "time" %in% names(weights)) {
      t <- if ("t" %in% names(weights)) time_steps(weights$t, time, "`weights`", "t") else
        time_steps(weights$time, time, "`weights`", "time")
      if (anyDuplicated(cbind(a, t)))
        stop("`weights` has more than one row for some area-time combinations", call. = FALSE)
      W <- matrix(NA_real_, m, T)
      W[cbind(a, t)] <- weights$weight
      if (anyNA(W)) stop(sprintf("`weights` has no value for %d area-time combinations", sum(is.na(W))),
                         call. = FALSE)
    } else {
      if (anyDuplicated(a)) stop("`weights` has duplicated ids (add a `time` or `t` column for weights over time)",
                                 call. = FALSE)
      v <- rep(NA_real_, m); v[a] <- weights$weight
      if (anyNA(v)) stop(sprintf("`weights` has no value for %d areas", sum(is.na(v))), call. = FALSE)
      W <- matrix(v, m, T)
    }
  } else if (is.matrix(weights)) {
    if (!identical(dim(weights), c(m, T)))
      stop(sprintf("a `weights` matrix must be %d areas x %d time steps", m, T), call. = FALSE)
    cols <- seq_len(T)
    if (!is.null(colnames(weights))) {
      cols <- match(time$labels, colnames(weights))
      if (anyNA(cols)) stop("the column names of a `weights` matrix must be the time labels", call. = FALSE)
    }
    W <- weights[by_id(rownames(weights)), cols, drop = FALSE]
  } else if (is.numeric(weights)) {
    if (length(weights) != m) stop(sprintf("a `weights` vector needs one value per area (%d)", m), call. = FALSE)
    W <- matrix(weights[by_id(names(weights))], m, T)
  } else {
    stop("`weights` must be NULL, a vector, a matrix or a data frame", call. = FALSE)
  }
  W <- matrix(as.numeric(W), m, T, dimnames = list(ids, NULL))
  if (any(!is.finite(W)) || any(W < 0)) stop("`weights` must be finite and >= 0", call. = FALSE)
  W
}

# ---- methods -------------------------------------------------------------------

#' @export
print.stareal_pred <- function(x, ...) {
  s <- x$settings
  cat(sprintf("<stareal_pred> %d areas x %d time steps from %d posterior draws\n", x$frame$m, x$time$n_steps,
              s$n_draws))
  if (length(s$groups)) cat(sprintf("  aggregated to: %s (%s weights)\n",
                                    paste(sprintf("%s (%d)", s$groups, vapply(s$groups, function(g)
                                      length(unique(x$groups[[g]]$id)), 1L)), collapse = ", "), s$weights))
  stat <- c("mean", "sd", paste0("q", 100 * s$probs), if (length(s$exceedance)) paste0("p_gt_", s$exceedance))
  cat(sprintf("  statistics: %s\n", paste(stat, collapse = ", ")))
  if (!is.null(x$draws)) cat(sprintf("  draws kept for %d area(s)\n", dim(x$draws)[1]))
  cat(sprintf("  tables: $area%s\n", if (length(s$groups)) paste0(", ", paste0("$groups$", s$groups, collapse = ", ")) else ""))
  invisible(x)
}

#' Predictions as a data frame
#'
#' @param x A `stareal_pred`.
#' @param row.names,optional Not used.
#' @param group `"area"` or one of the groupings predicted.
#' @param ... Not used.
#' @return The area table, or the group's table.
#' @export
as.data.frame.stareal_pred <- function(x, row.names = NULL, optional = FALSE, group = "area", ...) {
  pred_table(x, group)
}

pred_table <- function(x, group) {
  if (identical(group, "area")) return(x$area)
  if (!group %in% names(x$groups))
    stop(sprintf("`group` must be \"area\"%s", if (length(x$groups))
      paste0(" or one of ", paste(names(x$groups), collapse = ", ")) else " (no groupings were predicted)"),
      call. = FALSE)
  x$groups[[group]]
}

# a time step from a label or a step number
pred_time <- function(x, time) {
  if (length(time) != 1) stop("`time` must be one time label or step number", call. = FALSE)
  t <- match(as.character(time), x$time$labels)
  if (is.na(t) && is.numeric(time) && time == round(time) && time >= 1 && time <= x$time$n_steps) t <- as.integer(time)
  if (is.na(t)) stop(sprintf("`time` must be a time label (%s to %s) or a step number (1 to %d)",
                             x$time$labels[1], x$time$labels[x$time$n_steps], x$time$n_steps), call. = FALSE)
  t
}

#' Predictions as polygons
#'
#' The area (or group) table for one time step joined to the frame's
#' polygons, ready for mapping with your own tools (e.g.
#' `ggplot2::geom_sf()`). Groups are drawn by dissolving their areas.
#'
#' @param x A `stareal_pred` whose frame has polygons.
#' @param time A time label or step number.
#' @param group `"area"` or one of the groupings predicted.
#' @param ... Not used.
#' @return An `sf` object, one row per area (or group).
#' @exportS3Method sf::st_as_sf
st_as_sf.stareal_pred <- function(x, time, group = "area", ...) {
  if (is.null(x$frame$polygons)) stop("this frame has no polygons", call. = FALSE)
  need_pkg("sf", "to make polygons")
  t <- pred_time(x, time)
  tab <- pred_table(x, group)
  tab <- tab[tab$t == t, ]
  poly <- x$frame$polygons
  if (group == "area") {
    geom <- sf::st_geometry(poly)[match(tab$id, poly$id)]
  } else {
    geom <- with_s2_off(do.call(c, lapply(tab$id, function(g) sf::st_union(sf::st_geometry(poly)[poly[[group]] == g]))))
  }
  sf::st_sf(tab, geometry = geom)
}

#' Plot predictions
#'
#' With `time`, a map of one statistic for one time step. This needs a
#' frame with polygons. Groups are drawn by colouring their areas. Without
#' `time`, time series of chosen areas or groups, with the outer quantiles
#' as a band. With `data = TRUE`, the observed prevalence (positives /
#' tested) of every surveyed area and time step in the panel is added, with
#' its exact (Clopper-Pearson) 95% binomial interval. In a group's panel
#' these are the surveyed areas of the group, while the band is for the
#' group as a whole, so the two need not agree.
#'
#' @param x A `stareal_pred`.
#' @param group `"area"` or one of the groupings predicted.
#' @param time For a map: a time label (as in `x$time$labels`) or a time
#'   step number.
#' @param ids For time series: the area or group ids to show (default: up
#'   to six, those with the most people tested first).
#' @param stat For a map: the statistic to show (a column of the table,
#'   e.g. `"mean"`, `"q50"`, `"p_gt_5"`).
#' @param data For time series: add the observed prevalence and its exact
#'   binomial interval.
#' @param labels For time series: panel titles in place of the ids. Either
#'   one per id, in the order of `ids`, or a vector named by id (ids not
#'   named keep their id). Labels must be unique.
#' @param ... Not used.
#' @return A `ggplot` object (needs the `ggplot2` package). Prevalence is in
#'   percent.
#' @examples
#' fr <- areal_lattice(4, 4)
#' sim <- simulate_areal(fr, time_frame(n_steps = 5), n_cells = 30, seed = 1)
#' fit <- fit_areal(sim, control = list(n_iter = 100, n_burn = 50, n_store = 25, n_chains = 1,
#'                                      ess_sweeps = 2, verbose = FALSE))
#' pr <- predict(fit, verbose = FALSE)
#' plot(pr, ids = c("r1c1", "r1c2"), labels = c("first unit", "second unit"))
#' @export
plot.stareal_pred <- function(x, group = "area", time = NULL, ids = NULL, stat = "mean", data = TRUE,
                              labels = NULL, ...) {
  check_dots(...)
  need_pkg("ggplot2", "to plot predictions")
  .data <- ggplot2::.data
  tab <- pred_table(x, group)
  qcols <- grep("^q[0-9.]+$", names(tab), value = TRUE)
  qv <- as.numeric(sub("^q", "", qcols))
  if (!is.null(time)) {
    if (is.null(x$frame$polygons)) stop("this frame has no polygons to map", call. = FALSE)
    if (!stat %in% names(tab)) stop(sprintf("`stat` must be a column of the table, e.g. %s",
                                            paste(c("mean", qcols), collapse = ", ")), call. = FALSE)
    t <- pred_time(x, time)
    rows <- tab[tab$t == t, ]
    key <- if (group == "area") x$frame$ids else x$frame$groups[[group]]
    poly <- x$frame$polygons
    poly$value <- rows[[stat]][match(key, rows$id)]
    return(ggplot2::ggplot(poly) +
             ggplot2::geom_sf(ggplot2::aes(fill = .data$value), colour = "grey70", linewidth = 0.1) +
             ggplot2::scale_fill_viridis_c(name = stat) +
             ggplot2::labs(title = sprintf("%s, %s", if (group == "area") "areas" else group, x$time$labels[t])) +
             ggplot2::theme_minimal())
  }
  if (is.null(ids)) {
    tested <- if (group == "area") tapply(tab$k, tab$id, sum) else
      tapply(x$area$k, x$frame$groups[[group]][match(x$area$id, x$frame$ids)], sum)
    ids <- names(sort(tested, decreasing = TRUE))[seq_len(min(6, length(tested)))]
  }
  sub <- tab[tab$id %in% ids, ]
  if (!nrow(sub)) stop("none of `ids` are in the predictions", call. = FALSE)
  lab <- panel_labels(ids, labels)
  sub$id <- factor(sub$id, levels = ids, labels = lab)
  sub$x <- if (is.null(x$time$starts)) sub$t else x$time$starts[sub$t]   # dates on calendar frames
  mid <- if ("q50" %in% names(sub)) "q50" else "mean"
  lo <- qcols[which.min(qv)]; hi <- qcols[which.max(qv)]
  p <- ggplot2::ggplot(sub, ggplot2::aes(x = .data$x))
  if (x$time$n_steps > 1) {
    p <- p + ggplot2::geom_ribbon(ggplot2::aes(ymin = .data[[lo]], ymax = .data[[hi]]), fill = "#2563eb", alpha = 0.2) +
      ggplot2::geom_line(ggplot2::aes(y = .data[[mid]]), colour = "#1e40af")
  } else {
    p <- p + ggplot2::geom_pointrange(ggplot2::aes(y = .data[[mid]], ymin = .data[[lo]], ymax = .data[[hi]]),
                                      colour = "#1e40af")
  }
  obs <- if (data) observed_in(x, group, ids) else NULL
  if (!is.null(obs)) {
    obs$id <- factor(obs$id, levels = ids, labels = lab)
    obs$x <- if (is.null(x$time$starts)) obs$t else x$time$starts[obs$t]
    p <- p + ggplot2::geom_linerange(data = obs, ggplot2::aes(ymin = .data$ci_lo, ymax = .data$ci_hi),
                                     colour = "grey45", alpha = 0.6) +
      ggplot2::geom_point(data = obs, ggplot2::aes(y = .data$prev), size = 0.9, colour = "grey15")
  }
  p + ggplot2::facet_wrap(~id) +
    ggplot2::labs(x = if (is.null(x$time$starts)) "time step" else NULL, y = "prevalence (%)",
                  subtitle = paste0(sprintf("%s with %s-%s %s", mid, lo, hi, if (x$time$n_steps > 1) "band" else "interval"),
                                    if (!is.null(obs)) "; points: observed, with exact 95% binomial intervals")) +
    ggplot2::theme_minimal()
}

# Panel titles for a time-series plot: the ids, or `labels` (one per id in
# order, or named by id)
panel_labels <- function(ids, labels) {
  if (is.null(labels)) return(ids)
  labels <- stats::setNames(as.character(labels), names(labels))      # keeping the names
  if (!is.null(names(labels)) && all(nzchar(names(labels)))) {
    out <- ifelse(ids %in% names(labels), labels[match(ids, names(labels))], ids)
  } else {
    if (length(labels) != length(ids))
      stop(sprintf("`labels` must give one label per id (%d), or be named by id", length(ids)), call. = FALSE)
    out <- labels
  }
  if (anyDuplicated(out)) stop("`labels` must be unique, or panels would merge", call. = FALSE)
  out
}

# The surveyed area-time cells behind the panels `ids` of a time-series
# plot (areas, or the areas of groups), with the observed prevalence and
# its exact 95% binomial interval; NULL if there are none. `id` is the
# panel (area or group) id.
observed_in <- function(x, group, ids) {
  a <- x$area[x$area$k > 0, c("id", "t", "y", "k")]
  if (group != "area") a$id <- x$frame$groups[[group]][match(a$id, x$frame$ids)]
  a <- a[a$id %in% ids, ]
  if (!nrow(a)) return(NULL)
  ci <- binom_ci(a$y, a$k)
  cbind(a, prev = 100 * a$y / a$k, ci_lo = 100 * ci$lo, ci_hi = 100 * ci$hi)       # percent
}

# Exact (Clopper-Pearson) binomial interval for y positive of k
binom_ci <- function(y, k, level = 0.95) {
  a <- (1 - level) / 2
  data.frame(lo = ifelse(y == 0, 0, stats::qbeta(a, y, k - y + 1)),
             hi = ifelse(y == k, 1, stats::qbeta(1 - a, y + 1, k - y)))
}
