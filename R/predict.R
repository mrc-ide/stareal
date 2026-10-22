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
#   admin0 was;
# - the number of blocks also bounds memory (block_bytes per block);
# - weights from a matrix, a vector or a data frame (weights_matrix());
# - quantiles named q2.5, q50, ..., exceedance p_gt_<threshold>; no
#   rounding; areas in cell_id order (area, then time step);
# - optionally keeps the draws of chosen areas;
# - draws are selected evenly from the stored draws of all chains (the
#   Africa code used them all), and draw d uses seed `seed + d`.
#
# How it scales (Africa ADMIN2: m = 4,319, T up to 600, 1,000 draws):
# holding all draws would need ~70 GB, so the work is split in two phases.
#   Phase A (parallel over draws): the conditioning step, which needs only
#     the observed areas. Per draw it keeps w (length n), beta* and the seed.
#   Phase B (parallel over blocks of areas): for each draw, regenerate the
#     same Z from the seed and compute the block's rows of the surface; keep
#     the block's draws (rows x T x draws), summarise them, and aggregate
#     within the block.
#
# ------------------------------------------------------------------

#' Predict the full area-time surface from a fit
#'
#' Reconstructs the latent prevalence surface for every area and time step,
#' observed or not, once per posterior draw. It does this by conditional
#' simulation (Matheron's rule), which is exact given each draw. It then
#' summarises the draws per area and time step. Each grouping in `levels` is
#' aggregated draw by draw, as a weighted mean of its areas' prevalence,
#' and then summarised. Its uncertainty therefore accounts for the
#' correlation between areas.
#'
#' The surface is the field plus fixed effects, without the cell-level
#' nugget: it is the prevalence of the area-time cell, not of one survey in
#' it.
#'
#' @param object A [fit_areal()] object.
#' @param levels Names of groupings of the frame (see [areal_frame()]) to
#'   aggregate to. Areas are always summarised.
#' @param weights Weights of the areas in group aggregates, for example
#'   population or number of infections. The aggregate of a group is
#'   \eqn{\sum_a w_{at} p_{at} / \sum_a w_{at}} for each draw. Give one of:
#'   `NULL` (equal weights); a vector with one value per area (named by area
#'   id, or in frame order); an areas x time steps matrix (rows named by
#'   area id, or in frame order); or a data frame with columns `id`,
#'   `weight` and optionally `time` (time label) or `t` (time step), with
#'   every area-time combination present.
#' @param n_draws Number of posterior draws to use, taken evenly from the
#'   stored draws of all chains. `NULL`: all of them.
#' @param probs Quantiles to report.
#' @param exceedance Prevalence thresholds; for each, the posterior
#'   probability that prevalence exceeds it is reported.
#' @param keep_draws Ids of areas whose prevalence draws are returned in full.
#' @param seed Draw `d` uses seed `seed + d`. The caller's random-number
#'   stream is left untouched.
#' @param parallel,n_cores,blas_threads Run on worker processes, as in
#'   [mcmc_control()]. `n_cores` defaults to the number of cores less one.
#' @param block_bytes Approximate memory for one block's draws (bytes). More
#'   blocks are used if the surface does not fit.
#' @param progress Print progress messages.
#' @param ... Not used.
#' @return An object of class `stareal_pred`:
#'   \describe{
#'     \item{area}{one row per area and time step, in the order of area, then
#'       time step: `id`, the frame's groupings, `t`, `time`, the pooled
#'       data `y` and `k` (0 where there are none), then `mean`, `sd`, one
#'       column per quantile (`q2.5`, `q50`, ...) and per threshold
#'       (`p_gt_0.1`, ...), all on the prevalence scale}
#'     \item{<level>}{for each grouping in `levels`: `id` (the group), `t`,
#'       `time` and the same statistics}
#'     \item{draws}{if `keep_draws`: prevalence draws, an array (area x time
#'       step x draw)}
#'     \item{frame, time, settings}{for plotting and the record}
#'   }
#' @examples
#' fr <- areal_lattice(4, 4, blocks = c(2, 2))
#' tf <- time_frame(T = 5)
#' sim <- simulate_areal(fr, tf, n_cells = 30, seed = 1)
#' d <- areal_data(sim$data, fr, tf, y = "positive", k = "tested", area = "id", t = "t",
#'                 verbose = FALSE)
#' fit <- fit_areal(d, control = mcmc_control(n_iter = 100, n_burn = 50, n_store = 25,
#'                                            n_chains = 1, ess_sweeps = 2, progress = FALSE))
#' pr <- predict(fit, levels = "block", parallel = FALSE, progress = FALSE)
#' head(pr$area)
#' pr$block
#' @export
predict.stareal_fit <- function(object, levels = NULL, weights = NULL, n_draws = NULL,
                                probs = c(0.025, 0.5, 0.975), exceedance = NULL, keep_draws = NULL,
                                seed = 1, parallel = TRUE, n_cores = NULL, blas_threads = 1,
                                block_bytes = 5e8, progress = TRUE, ...) {
  fit <- object
  fr <- fit$data$frame; tf <- fit$data$time; m <- fr$m; T <- tf$T
  gnames <- setdiff(names(fr$groups), "id")
  if (length(setdiff(levels, gnames)))
    stop(sprintf("level(s) not among the frame's groupings (%s): %s",
                 if (length(gnames)) paste(gnames, collapse = ", ") else "none",
                 paste(setdiff(levels, gnames), collapse = ", ")), call. = FALSE)
  if (!is.numeric(probs) || any(probs <= 0 | probs >= 1)) stop("`probs` must be in (0, 1)", call. = FALSE)
  if (!is.null(exceedance) && (!is.numeric(exceedance) || any(exceedance <= 0 | exceedance >= 1)))
    stop("`exceedance` must be prevalences in (0, 1)", call. = FALSE)
  keep <- if (is.null(keep_draws)) integer(0) else match(as.character(keep_draws), fr$ids)
  if (anyNA(keep)) stop(sprintf("`keep_draws` ids not in the frame: %s",
                                paste(keep_draws[is.na(keep)], collapse = ", ")), call. = FALSE)
  W <- weights_matrix(weights, fr, tf)
  for (lv in levels) {
    den <- rowsum(W, fr$groups[[lv]])
    if (any(den <= 0)) stop(sprintf("some groups of `%s` have zero total weight in some time steps", lv),
                            call. = FALSE)
  }

  # posterior draws, pooled chain by chain, then thinned evenly
  pool <- function(a) matrix(a, prod(dim(a)[1:2]), dim(a)[3])
  th <- pool(fit$draws$theta); xo <- pool(fit$draws$x)
  D_all <- nrow(xo)
  if (is.null(n_draws)) n_draws <- D_all
  if (n_draws < 1 || n_draws > D_all)
    stop(sprintf("`n_draws` must be between 1 and the %d stored draws", D_all), call. = FALSE)
  idx <- unique(round(seq(1, D_all, length.out = n_draws)))
  D <- length(idx)
  seeds <- seed + seq_len(D)
  if (is.null(n_cores)) n_cores <- max(1L, parallel::detectCores() - 1L, na.rm = TRUE)
  n_workers <- if (parallel) as.integer(n_cores) else 1L

  ctx <- list(eig = fr$eig, obs = fit$setup$obs, b0 = fit$setup$b0, B = fit$setup$B, T = T,
              design = design_rows(fit))
  check <- function(res, what) {
    bad <- which(vapply(res, inherits, logical(1), "try-error"))
    if (length(bad)) stop(sprintf("prediction %s failed: %s", what,
                                  conditionMessage(attr(res[[bad[1]]], "condition"))), call. = FALSE)
    res
  }

  # ---- phase A: conditioning, in chunks of draws ----
  say(progress, "phase A: conditioning %d draws on %d observed cells", D, ctx$obs$n)
  chunks <- split(seq_len(D), cut(seq_len(D), min(D, max(1, n_workers) * 4), labels = FALSE))
  jobA <- function(ds, sh) lapply(ds, function(d)
    utils::getFromNamespace("prepare_draw", "stareal")(sh$x[d, ], sh$theta[d, ], sh$seeds[d], sh$ctx))
  prep <- check(par_lapply(chunks, jobA, shared = list(x = xo[idx, , drop = FALSE], theta = th[idx, , drop = FALSE],
                                                     seeds = seeds, ctx = ctx),
                           n_workers = n_workers, blas_threads = blas_threads, quiet = !progress), "phase A")
  prep <- unlist(prep, recursive = FALSE)

  # ---- phase B: blocks of areas ----
  bytes <- 8 * m * T * D
  n_blocks <- max(n_workers, ceiling(bytes / block_bytes))
  bg <- block_grouping(fr$groups, levels, n_blocks)
  blocks <- make_blocks(if (is.null(bg)) fr$ids else fr$groups[[bg]], n_blocks)
  nested <- vapply(levels, function(lv) !is.null(bg) && nests_in(fr$groups[[lv]], fr$groups[[bg]]), logical(1))
  say(progress, "phase B: %d block%s of areas (%d-%d areas each)%s", length(blocks), if (length(blocks) > 1) "s" else "",
      min(lengths(blocks)), max(lengths(blocks)), if (is.null(bg)) "" else sprintf(", whole groups of `%s`", bg))
  stats <- list(probs = probs, exceedance = exceedance)
  groups <- fr$groups[levels]
  jobB <- function(b, sh) utils::getFromNamespace("predict_block", "stareal")(b, sh)
  res <- check(par_lapply(seq_along(blocks), jobB,
                          shared = list(ctx = c(ctx, list(prep = prep)), blocks = blocks, W = W, groups = groups,
                                        nested = nested, stats = stats, keep = keep),
                          n_workers = n_workers, blas_threads = blas_threads, quiet = !progress), "phase B")

  # ---- assemble ----
  out <- list()
  u <- do.call(rbind, lapply(res, `[[`, "area"))
  u <- u[order(u$area, u$t), ]
  cells <- fit$data$cells
  ci <- (u$area - 1L) * T + u$t
  yk <- match(ci, cells$cell_id)
  lab <- data.frame(id = fr$ids[u$area], fr$groups[u$area, gnames, drop = FALSE], t = u$t,
                    time = tf$labels[u$t], y = ifelse(is.na(yk), 0L, cells$y[yk]),
                    k = ifelse(is.na(yk), 0L, cells$k[yk]), stringsAsFactors = FALSE)
  out$area <- cbind(lab, u[setdiff(names(u), c("area", "t"))])
  rownames(out$area) <- NULL
  for (lv in levels) {
    if (nested[[lv]]) {
      g <- do.call(rbind, lapply(res, function(r) r$levels[[lv]]))
    } else {
      lvl <- sort(unique(groups[[lv]]))
      G <- length(lvl)
      num <- matrix(0, G * T, D); den <- matrix(0, G, T)
      for (r in res) {
        ps <- r$levels[[lv]]
        i0 <- match(ps$id, lvl)
        rows <- as.vector(outer(i0, (seq_len(T) - 1) * G, `+`))       # (group, t) at (t - 1) * G + i
        num[rows, ] <- num[rows, ] + ps$num
        den[i0, ] <- den[i0, ] + ps$den
      }
      g <- cbind(data.frame(id = rep(lvl, times = T), t = rep(seq_len(T), each = G), stringsAsFactors = FALSE),
                 summarise_draw_matrix(num / as.vector(den), probs, exceedance))
    }
    g <- g[order(g$id, g$t), ]
    out[[lv]] <- cbind(data.frame(id = g$id, t = g$t, time = tf$labels[g$t], stringsAsFactors = FALSE),
                       g[setdiff(names(g), c("id", "t"))])
    rownames(out[[lv]]) <- NULL
  }
  if (length(keep)) {
    dr <- array(NA_real_, c(length(keep), T, D), dimnames = list(fr$ids[keep], tf$labels, NULL))
    for (r in res) if (length(r$draws)) for (a in names(r$draws)) dr[a, , ] <- r$draws[[a]]
    out$draws <- dr
  }
  out$frame <- fr; out$time <- tf
  out$settings <- list(levels = levels, n_draws = D, draws_used = idx, probs = probs, exceedance = exceedance,
                       seed = seed, weights = if (is.null(weights)) "equal" else "given", block_grouping = bg,
                       n_blocks = length(blocks))
  structure(out, class = "stareal_pred")
}

# ---- phase B for one block ----------------------------------------------------

# sh: ctx (with prep, the prepare_draw() outputs), blocks, W (m x T weights),
# groups (data frame of the levels' group vectors, frame order), nested
# (per level: complete inside every block), stats, keep (areas whose draws
# are returned). Returns the block's area summaries; per level either its
# summaries (nested) or partial sums (num per group, time and draw; den per
# group and time); and the kept draws.
predict_block <- function(b, sh) {
  ctx <- sh$ctx; rows <- sh$blocks[[b]]; nr <- length(rows); T <- ctx$T; D <- length(ctx$prep)
  P <- matrix(NA_real_, nr * T, D)                 # cell (area i, step t) at (t - 1) * nr + i
  parts <- rows_parts(rows, ctx)
  for (d in seq_len(D)) P[, d] <- stats::plogis(as.vector(predict_rows(rows, ctx$prep[[d]], ctx, parts)))
  st <- sh$stats
  out <- list(area = cbind(data.frame(area = rep(rows, times = T), t = rep(seq_len(T), each = nr)),
                           summarise_draw_matrix(P, st$probs, st$exceedance)))
  Wb <- sh$W[rows, , drop = FALSE]                 # nr x T
  out$levels <- list()
  for (lv in names(sh$groups)) {
    gv <- sh$groups[[lv]][rows]
    lvl <- sort(unique(gv))
    G <- group_matrix(gv, lvl)
    den <- as.matrix(G %*% Wb)                     # groups x T
    num <- matrix(NA_real_, length(lvl) * T, D)
    for (d in seq_len(D)) num[, d] <- as.vector(as.matrix(G %*% (Wb * matrix(P[, d], nr, T))))
    out$levels[[lv]] <- if (sh$nested[[lv]]) {
      cbind(data.frame(id = rep(lvl, times = T), t = rep(seq_len(T), each = length(lvl)),
                       stringsAsFactors = FALSE),
            summarise_draw_matrix(num / as.vector(den), st$probs, st$exceedance))
    } else list(id = lvl, num = num, den = den)
  }
  kb <- intersect(sh$keep, rows)
  out$draws <- lapply(stats::setNames(kb, ctx_ids(sh, kb)), function(a) {
    i <- match(a, rows)
    P[(seq_len(T) - 1) * nr + i, , drop = FALSE]   # T x D
  })
  out
}

ctx_ids <- function(sh, a) rownames(sh$W)[a]

# rows = cells, columns = draws (prevalence) -> mean, sd, quantiles and
# exceedance probabilities
summarise_draw_matrix <- function(P, probs, exceedance) {
  q <- matrixStats::rowQuantiles(P, probs = probs, drop = FALSE)
  colnames(q) <- paste0("q", 100 * probs)
  out <- data.frame(mean = rowMeans(P), sd = matrixStats::rowSds(P), q, check.names = FALSE)
  for (e in exceedance) out[[paste0("p_gt_", e)]] <- rowMeans(P > e)
  out
}

# ---- blocks and groupings ----------------------------------------------------

# Blocks of areas for phase B: whole groups (one label per area), packed
# largest-first into n_blocks blocks of roughly equal size
make_blocks <- function(parent, n_blocks) {
  sizes <- sort(table(parent), decreasing = TRUE)
  load <- numeric(n_blocks); members <- vector("list", n_blocks)
  for (p in names(sizes)) {
    b <- which.min(load)
    members[[b]] <- c(members[[b]], p); load[b] <- load[b] + sizes[[p]]
  }
  members <- members[lengths(members) > 0]
  lapply(members, function(ps) which(parent %in% ps))
}

# is every group of `g` inside one group of `h`?
nests_in <- function(g, h) all(tapply(h, g, function(v) length(unique(v))) == 1)

# The grouping whose whole groups make the blocks: among the levels with at
# least n_blocks groups, the one with the most levels nested in it (ties:
# the one with more groups, for balance). NULL: blocks of single areas,
# and every level is assembled from partial sums.
block_grouping <- function(groups, levels, n_blocks) {
  cand <- levels[vapply(levels, function(lv) length(unique(groups[[lv]])) >= n_blocks, logical(1))]
  if (!length(cand)) return(NULL)
  score <- vapply(cand, function(c) sum(vapply(levels, function(lv) nests_in(groups[[lv]], groups[[c]]),
                                               logical(1))), numeric(1))
  n_groups <- vapply(cand, function(c) length(unique(groups[[c]])), numeric(1))
  cand[order(-score, -n_groups)][1]
}

# group-membership matrix: one row per group in `levels`, one column per area
group_matrix <- function(group, levels) {
  Matrix::sparseMatrix(i = match(group, levels), j = seq_along(group), x = 1,
                       dims = c(length(levels), length(group)))
}

# ---- weights ------------------------------------------------------------------

# the user's weights as an m x T matrix in frame order, rows named by id
weights_matrix <- function(weights, frame, time) {
  m <- frame$m; T <- time$T; ids <- frame$ids
  by_id <- function(v, nm) {                       # reorder named rows/entries to frame order
    if (is.null(nm)) return(v)
    i <- match(ids, nm)
    if (anyNA(i)) stop(sprintf("`weights` missing for %d area(s), e.g. %s", sum(is.na(i)), ids[which(is.na(i))[1]]),
                       call. = FALSE)
    i
  }
  if (is.null(weights)) {
    W <- matrix(1, m, T)
  } else if (is.data.frame(weights)) {
    if (!all(c("id", "weight") %in% names(weights)))
      stop("a `weights` data frame needs columns `id` and `weight`", call. = FALSE)
    a <- match(as.character(weights$id), ids)
    if (anyNA(a)) stop("`weights` has ids not in the frame", call. = FALSE)
    if ("t" %in% names(weights) || "time" %in% names(weights)) {
      t <- if ("t" %in% names(weights)) as.integer(weights$t) else match(as.character(weights$time), time$labels)
      if (anyNA(t) || any(t < 1 | t > T)) stop("`weights` has time steps outside the time frame", call. = FALSE)
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
    W <- weights[by_id(seq_len(m), rownames(weights)), , drop = FALSE]
  } else if (is.numeric(weights)) {
    if (length(weights) != m) stop(sprintf("a `weights` vector needs one value per area (%d)", m), call. = FALSE)
    W <- matrix(weights[by_id(seq_len(m), names(weights))], m, T)
  } else {
    stop("`weights` must be NULL, a vector, a matrix or a data frame", call. = FALSE)
  }
  W <- matrix(as.numeric(W), m, T, dimnames = list(ids, NULL))
  if (any(!is.finite(W)) || any(W < 0)) stop("`weights` must be finite and >= 0", call. = FALSE)
  W
}

# ---- print and plot ----------------------------------------------------------

#' @export
print.stareal_pred <- function(x, ...) {
  s <- x$settings
  cat(sprintf("<stareal_pred> %d areas x %d time steps from %d posterior draws\n", x$frame$m, x$time$T, s$n_draws))
  if (length(s$levels)) cat(sprintf("  aggregated to: %s (%s weights)\n",
                                    paste(sprintf("%s (%d)", s$levels, vapply(s$levels, function(l)
                                      length(unique(x[[l]]$id)), 1L)), collapse = ", "), s$weights))
  stat <- setdiff(names(x$area), c("id", setdiff(names(x$frame$groups), "id"), "t", "time", "y", "k"))
  cat(sprintf("  statistics: %s\n", paste(stat, collapse = ", ")))
  if (!is.null(x$draws)) cat(sprintf("  draws kept for %d area(s)\n", dim(x$draws)[1]))
  cat(sprintf("  tables: %s\n", paste(c("area", s$levels), collapse = ", ")))
  invisible(x)
}

#' Plot predictions
#'
#' With `time`, a map of one statistic for one time step. This needs a
#' frame with polygons. Groups of a level are drawn by colouring their areas.
#' Without `time`, time series of chosen units or groups, with the outer
#' quantiles as a band and the observed prevalence (positives / tested) as
#' points sized by the number tested.
#'
#' @param x A `stareal_pred`.
#' @param level `"area"` or one of the levels predicted.
#' @param time For a map: a time label (as in `x$time$labels`) or a time
#'   step number.
#' @param ids For time series: the area or group ids to show (default: up
#'   to six, those with the most data first).
#' @param stat For a map: the statistic to show (a column of the table,
#'   e.g. `"mean"`, `"q50"`, `"p_gt_0.1"`).
#' @param ... Not used.
#' @return A `ggplot` object (needs the `ggplot2` package).
#' @export
plot.stareal_pred <- function(x, level = "area", time = NULL, ids = NULL, stat = "mean", ...) {
  need_pkg("ggplot2", "to plot predictions")
  .data <- ggplot2::.data
  tab <- x[[level]]
  if (!is.data.frame(tab)) stop(sprintf("no predictions at level `%s`", level), call. = FALSE)
  qcols <- grep("^q[0-9.]+$", names(tab), value = TRUE)
  if (!is.null(time)) {
    if (is.null(x$frame$polygons)) stop("this frame has no polygons to map", call. = FALSE)
    if (!stat %in% names(tab)) stop(sprintf("`stat` must be a column of the table, e.g. %s",
                                            paste(c("mean", qcols), collapse = ", ")), call. = FALSE)
    t <- if (is.numeric(time)) as.integer(time) else match(as.character(time), x$time$labels)
    if (is.na(t) || t < 1 || t > x$time$T) stop("`time` is not a time step of the predictions", call. = FALSE)
    rows <- tab[tab$t == t, ]
    key <- if (level == "area") x$frame$ids else x$frame$groups[[level]]
    poly <- x$frame$polygons
    poly$value <- rows[[stat]][match(key, rows$id)]
    return(ggplot2::ggplot(poly) +
             ggplot2::geom_sf(ggplot2::aes(fill = .data$value), colour = "grey70", linewidth = 0.1) +
             ggplot2::scale_fill_viridis_c(name = stat) +
             ggplot2::labs(title = sprintf("%s, %s", if (level == "area") "areas" else level, x$time$labels[t])) +
             ggplot2::theme_minimal())
  }
  if (is.null(ids)) {
    if (level == "area") {
      tot <- tapply(tab$k, tab$id, sum)
      ids <- names(sort(tot, decreasing = TRUE))[seq_len(min(6, length(tot)))]
    } else ids <- utils::head(unique(tab$id), 6)
  }
  sub <- tab[tab$id %in% ids, ]
  if (!nrow(sub)) stop("none of `ids` are in the predictions", call. = FALSE)
  sub$id <- factor(sub$id, levels = ids)
  sub$x <- if (is.null(x$time$starts)) sub$t else x$time$starts[sub$t]   # dates on calendar frames
  mid <- if ("q50" %in% names(sub)) "q50" else "mean"
  lo <- qcols[1]; hi <- qcols[length(qcols)]
  p <- ggplot2::ggplot(sub, ggplot2::aes(x = .data$x)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = .data[[lo]], ymax = .data[[hi]]), fill = "#2563eb", alpha = 0.2) +
    ggplot2::geom_line(ggplot2::aes(y = .data[[mid]]), colour = "#1e40af") +
    ggplot2::facet_wrap(~id) +
    ggplot2::labs(x = if (is.null(x$time$starts)) "time step" else NULL, y = "prevalence",
                  subtitle = sprintf("%s with %s-%s band", mid, lo, hi)) +
    ggplot2::theme_minimal()
  if (level == "area" && any(sub$k > 0))
    p <- p + ggplot2::geom_point(data = sub[sub$k > 0, ], ggplot2::aes(y = .data$y / .data$k, size = .data$k),
                                 shape = 21, fill = "white") +
      ggplot2::scale_size_area(name = "tested", max_size = 3)
  p
}
