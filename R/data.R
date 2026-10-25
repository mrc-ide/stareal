# data.R
#
# Survey data to pooled area-time cells. Generalises
# Africa_areal/pipeline/04_cells.R (2026-10-12): units from ids or from
# coordinates (point-in-polygon, with a snap to the nearest unit for points
# just outside every polygon); time from dates, step numbers or labels;
# surveys in the same cell pooled (y and k add, since they are binomial with
# the same probability). Every exclusion is counted here, with its reason,
# and nowhere else.
#
# ------------------------------------------------------------------

#' Attach survey data to a frame
#'
#' Assigns each survey to an areal unit and a time step, then pools surveys
#' that fall in the same unit and step: positives add, and so do the
#' numbers tested. The pooled cells are the only data the model sees.
#' Surveys that cannot be used are kept, with the reason, so that every
#' exclusion can be accounted for: print the result to see them.
#'
#' Reasons for exclusion, checked in this order:
#' \itemize{
#'   \item missing or invalid counts: `y` or `k` missing or unreadable, not
#'     whole numbers, negative, `k = 0`, or `y > k`. These come with a
#'     warning, because they are usually data errors.
#'   \item missing area, or missing or unreadable coordinates;
#'   \item area id not in the frame. This also comes with a warning.
#'   \item outside the frame: coordinates more than `snap_max_m` from every
#'     unit;
#'   \item missing date or time step;
#'   \item unreadable date (not `YYYY-MM-DD`), or a time step that is not a
#'     step number or label of the time frame. These come with a warning.
#'   \item outside the time window.
#' }
#'
#' Points that fall outside every polygon but within `snap_max_m` of one
#' are assigned to the nearest unit. Survey coordinates are often rounded,
#' and coastal sites fall just offshore. Snap distances are recorded in
#' `$surveys$.snap_m`. A point inside two overlapping units goes to the
#' first, with a warning.
#'
#' @param data A data frame, one row per survey.
#' @param frame A [areal_frame()].
#' @param time A [time_frame()].
#' @param y,k Names of the columns holding the number of positives and the
#'   number tested.
#' @param area Name of a column of unit ids matching `frame$ids`. Give this
#'   or `coords`.
#' @param coords Names of the x and y (longitude and latitude) columns.
#'   Needs a frame built from polygons.
#' @param crs Coordinate reference system of `coords` (default WGS84
#'   longitude/latitude).
#' @param date Name of a column of dates (`Date`, date-times, or
#'   `"YYYY-MM-DD"` text). Give this or `t`.
#' @param t Name of a column of time steps: step numbers (1..n_steps), or
#'   the time frame's labels (e.g. `"2020-03"`). Required for a plain
#'   index time frame.
#' @param snap_max_m Maximum distance (metres) for snapping points outside
#'   every polygon to the nearest unit. Use 0 for no snapping.
#' @param verbose Report, in one line, how many surveys were used and how
#'   many cells they make.
#'
#' @return An object of class `stareal_data`:
#'   \describe{
#'     \item{cells}{one row per observed cell, ordered by `cell_id`
#'       (`(area - 1) * n_steps + t`): `cell_id`, `area` (row of the
#'       frame), `id`, `t`, `time` (label), `y`, `k`, `n_surveys`}
#'     \item{surveys}{`data` with `.id` (the area id, `NA` if none), `.t`,
#'       `.snap_m` and `.excluded` (`NA` if used, otherwise the reason)
#'       added}
#'     \item{frame, time}{the frame and time frame}
#'     \item{truth}{`NULL` (see [simulate_areal()])}
#'   }
#'   Print it for the exclusions; `plot()` shows where and when the data are.
#' @examples
#' fr <- areal_lattice(3, 3)
#' tf <- time_frame(n_steps = 4)
#' surveys <- data.frame(unit = c("r1c1", "r1c1", "r2c3", "r9c9"),
#'                       step = c(1, 1, 4, 2), pos = c(3, 1, 0, 2), n = c(10, 5, 8, 9))
#' d <- areal_data(surveys, fr, tf, y = "pos", k = "n", area = "unit", t = "step")
#' d
#' d$cells
#' @export
areal_data <- function(data, frame, time, y, k, area = NULL, coords = NULL, crs = 4326,
                       date = NULL, t = NULL, snap_max_m = 10000, verbose = TRUE) {
  if (!inherits(frame, "stareal_frame")) stop("`frame` must be made by areal_frame()", call. = FALSE)
  if (!inherits(time, "stareal_time")) stop("`time` must be made by time_frame()", call. = FALSE)
  if (inherits(data, "sf")) data <- sf::st_drop_geometry(data)
  data <- as.data.frame(data)
  if (is.null(area) == is.null(coords)) stop("give exactly one of `area` and `coords`", call. = FALSE)
  if (is.null(date) == is.null(t)) stop("give exactly one of `date` and `t`", call. = FALSE)
  if (time$step == "index" && !is.null(date))
    stop("`time` is a plain index frame: give the time step as a column of step numbers (`t`)", call. = FALSE)
  if (!is.null(coords) && length(coords) != 2) stop("`coords` must name two columns (x, y)", call. = FALSE)
  cols <- c(y, k, area, coords, date, t)
  miss <- setdiff(cols, names(data))
  if (length(miss)) stop(sprintf("column(s) not in `data`: %s", paste(miss, collapse = ", ")), call. = FALSE)
  if (!nrow(data)) stop("`data` has no rows", call. = FALSE)

  N <- nrow(data)
  why <- rep(NA_character_, N)
  flag <- function(cond, reason) {
    cond[is.na(cond)] <- FALSE                     # e.g. comparisons with a missing count
    cond <- cond & is.na(why)
    why[cond] <<- reason
    sum(cond)
  }

  # ---- counts ----
  yy <- as_num(data[[y]]); kk <- as_num(data[[k]])
  n_bad_counts <- flag(is.na(yy) | is.na(kk), "missing or unreadable count") +
    flag(yy != round(yy) | kk != round(kk) | yy < 0 | kk < 0, "invalid count (not a whole number >= 0)") +
    flag(kk == 0, "none tested (k = 0)") +
    flag(yy > kk, "positives exceed tested (y > k)")

  # ---- area ----
  snap <- rep(NA_real_, N)
  n_unknown <- 0; n_overlap <- 0
  if (!is.null(area)) {
    aid <- id_chr(data[[area]])
    flag(is.na(aid) | trimws(aid) == "", "missing area")
    a <- match(aid, frame$ids)
    n_unknown <- flag(is.na(a), "area not in frame")
  } else {
    if (is.null(frame$polygons)) stop("assigning coordinates needs a frame built from polygons", call. = FALSE)
    xn <- as_num(data[[coords[1]]]); yn <- as_num(data[[coords[2]]])
    flag(!is.finite(xn) | !is.finite(yn), "missing or unreadable coordinates")
    res <- assign_points(xn, yn, crs, frame, snap_max_m, ok = is.na(why))
    a <- res$area; snap <- res$snap_m; n_overlap <- res$n_overlap
    flag(is.na(a), sprintf("outside frame (> %g km from every unit)", snap_max_m / 1e3))
  }

  # ---- time ----
  n_bad_time <- 0
  if (!is.null(date)) {
    raw <- data[[date]]
    dd <- to_date(raw, sprintf("column `%s`", date))
    flag(is.na(raw) | (is.character(raw) & trimws(as.character(raw)) == ""), "missing date")
    n_bad_time <- flag(is.na(dd), "unreadable date (expected YYYY-MM-DD)")
    tt <- time_index(time, dd)
  } else {
    raw <- data[[t]]
    flag(is.na(raw), "missing time step")
    tt <- time_steps(raw, time, "", "t", na_ok = TRUE)
    num <- as_num(raw)
    outside <- is.na(tt) & !is.na(num) & num == round(num)          # a whole number, but not 1..T
    n_bad_time <- flag(is.na(tt) & !outside, "time step not a step number or label of the time frame")
    if (sum(outside & is.na(why)) && all(num[outside & is.na(why)] > time$n_steps) && time$step != "index")
      warning(sprintf("`%s` holds numbers above %d, the number of steps: give step numbers (1..%d) or labels (e.g. \"%s\"), not years",
                      t, time$n_steps, time$n_steps, time$labels[1]), call. = FALSE)
  }
  flag(is.na(tt), "outside time window")

  # ---- pool ----
  use <- is.na(why)
  if (!any(use)) {
    ex <- table(why)
    stop("no survey could be used: ", paste(sprintf("%d %s", ex, names(ex)), collapse = "; "), call. = FALSE)
  }
  T <- time$n_steps
  cid <- (a[use] - 1L) * T + tt[use]
  pooled_y <- rowsum(yy[use], cid); pooled_k <- rowsum(kk[use], cid)
  ns <- rowsum(rep(1L, sum(use)), cid)
  cell_id <- as.integer(rownames(pooled_y))                   # rowsum sorts the groups
  ca <- (cell_id - 1L) %/% T + 1L; ct <- (cell_id - 1L) %% T + 1L
  # counts stay doubles: pooled totals can exceed the integer range
  cells <- data.frame(cell_id = cell_id, area = ca, id = frame$ids[ca], t = ct, time = time$labels[ct],
                      y = pooled_y[, 1], k = pooled_k[, 1], n_surveys = as.integer(ns[, 1]),
                      stringsAsFactors = FALSE)

  surveys <- data
  surveys$.id <- frame$ids[a]
  surveys$.t <- tt
  surveys$.snap_m <- snap
  surveys$.excluded <- why

  plural <- function(n) if (n > 1) "s" else ""
  if (n_bad_counts) warning(sprintf("%d survey%s excluded for missing or invalid counts", n_bad_counts,
                                    plural(n_bad_counts)), call. = FALSE)
  if (n_unknown) warning(sprintf("%d survey%s with area ids not in the frame", n_unknown, plural(n_unknown)),
                         call. = FALSE)
  if (n_bad_time) warning(sprintf("%d survey%s excluded for an unreadable date or time step", n_bad_time,
                                  plural(n_bad_time)), call. = FALSE)
  if (n_overlap) warning(sprintf("%d survey point%s fall in more than one unit (overlapping polygons); assigned to the first",
                                 n_overlap, plural(n_overlap)), call. = FALSE)
  say(verbose, "%d of %d surveys used, pooled into %d cells (print the result for the exclusions)",
      sum(use), N, nrow(cells))
  structure(list(cells = cells, surveys = surveys, frame = frame, time = time, truth = NULL),
            class = "stareal_data")
}

# point-in-polygon on the plane, then a snap to the nearest unit (by
# distance in the frame's planar projection) for points within snap_max_m
assign_points <- function(xn, yn, crs, frame, snap_max_m, ok) {
  N <- length(xn)
  area <- rep(NA_integer_, N); snap_m <- rep(NA_real_, N)
  idx <- which(ok)
  if (!length(idx)) return(list(area = area, snap_m = snap_m, n_overlap = 0L))
  with_s2_off({
    pts <- sf::st_as_sf(data.frame(x = xn[idx], y = yn[idx]), coords = c("x", "y"), crs = crs)
    pts <- sf::st_transform(pts, sf::st_crs(frame$polygons))
    hit <- suppressMessages(sf::st_within(pts, frame$polygons))   # planar on purpose (with_s2_off)
    n_overlap <- sum(lengths(hit) > 1)
    a <- vapply(hit, function(h) if (length(h)) h[1] else NA_integer_, integer(1))
    s <- ifelse(is.na(a), NA_real_, 0)
    miss <- which(is.na(a))
    if (length(miss) && snap_max_m > 0) {
      pp <- sf::st_transform(pts[miss, ], frame$crs_planar)
      ap <- sf::st_transform(frame$polygons, frame$crs_planar)
      near <- sf::st_nearest_feature(pp, ap)
      dist <- as.numeric(sf::st_distance(pp, ap[near, ], by_element = TRUE))
      s[miss] <- dist
      a[miss[dist <= snap_max_m]] <- near[dist <= snap_max_m]
    }
  })
  area[idx] <- a; snap_m[idx] <- s
  list(area = area, snap_m = snap_m, n_overlap = n_overlap)
}

#' @export
print.stareal_data <- function(x, ...) {
  cl <- x$cells; s <- x$surveys
  big <- function(v) format(v, big.mark = ",", scientific = FALSE, trim = TRUE)
  cat(sprintf("<stareal_data> %s of %s surveys used, pooled into %s cells%s\n",
              big(sum(is.na(s$.excluded))), big(nrow(s)), big(nrow(cl)),
              if (!is.null(x$truth)) " (simulated)" else ""))
  cat(sprintf("  %d of %d units and %d of %d time steps have data (%.2f%% of the %s area-time cells)\n",
              length(unique(cl$area)), x$frame$m, length(unique(cl$t)), x$time$n_steps,
              100 * nrow(cl) / (x$frame$m * x$time$n_steps), big(x$frame$m * x$time$n_steps)))
  cat(sprintf("  tested per cell: median %s, range %s-%s; total %s tested, %s positive\n",
              big(stats::median(cl$k)), big(min(cl$k)), big(max(cl$k)), big(sum(cl$k)), big(sum(cl$y))))
  snapped <- sum(!is.na(s$.snap_m) & s$.snap_m > 0 & is.na(s$.excluded))
  if (snapped) cat(sprintf("  %d survey%s snapped to the nearest unit (max %.1f km)\n", snapped,
                           if (snapped > 1) "s" else "",
                           max(s$.snap_m[is.na(s$.excluded)], na.rm = TRUE) / 1e3))
  ex <- table(s$.excluded)
  if (length(ex)) {
    cat("  excluded:\n")
    for (r in names(ex)) cat(sprintf("    %6d  %s\n", ex[[r]], r))
  }
  if (!is.null(x$truth)) {
    th <- x$truth$theta
    cat(sprintf("  truth: sigma %.3g, rho_s %.3g, rho_t %.3g, sigma_n %.3g; %s\n",
                th[["sigma"]], th[["rho_s"]], th[["rho_t"]], th[["sigma_n"]],
                paste(sprintf("%s %.3g", names(x$truth$beta), x$truth$beta), collapse = ", ")))
  }
  invisible(x)
}

#' @export
as.data.frame.stareal_data <- function(x, ...) x$cells

#' Plot survey data
#'
#' Two panels. The first maps the total tested in each unit, if the frame
#' has polygons. The second shows the total tested in each time step.
#' @param x A `stareal_data`.
#' @param ... Unused.
#' @return `x`, invisibly.
#' @export
plot.stareal_data <- function(x, ...) {
  cl <- x$cells
  per <- function(index, n) {                      # total tested per unit / per time step
    out <- numeric(n); s <- rowsum(cl$k, index); out[as.integer(rownames(s))] <- s[, 1]; out
  }
  has_map <- !is.null(x$frame$polygons) && requireNamespace("sf", quietly = TRUE)
  op <- graphics::par(mfrow = c(1, if (has_map) 2 else 1), mar = c(4, 4, 2, 1))
  on.exit(graphics::par(op))
  if (has_map) {
    tested <- per(cl$area, x$frame$m)
    br <- unique(c(0, 1, stats::quantile(tested[tested > 0], c(0.25, 0.5, 0.75, 1))))
    pal <- c("grey95", grDevices::hcl.colors(length(br) - 2, "Blues 3", rev = TRUE))
    cls <- cut(tested, br, right = FALSE, include.lowest = TRUE)
    graphics::plot(sf::st_geometry(x$frame$polygons), col = pal[as.integer(cls)], border = "grey70",
                   lwd = 0.2, main = "tested per unit")
    graphics::legend("bottomleft", legend = levels(cls), fill = pal, title = "tested", bty = "n", cex = 0.7)
  }
  tk <- per(cl$t, x$time$n_steps)
  graphics::barplot(tk, names.arg = x$time$labels, border = NA, col = "#2563eb",
                    main = "tested per time step", las = 2, cex.names = 0.6)
  invisible(x)
}
