# data.R
#
# Survey data to pooled area-time cells. Generalises
# Africa_areal/pipeline/04_cells.R (2026-10-12): units from ids or from
# coordinates (point-in-polygon, with a snap to the nearest unit for points
# just outside every polygon); time from dates or from a step index;
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
#' exclusion can be accounted for.
#'
#' Reasons for exclusion, checked in this order:
#' \itemize{
#'   \item missing or invalid counts: `y` or `k` missing, not whole numbers,
#'     negative, `k = 0`, or `y > k`. These come with a warning, because
#'     they are usually data errors.
#'   \item missing area or coordinates;
#'   \item area id not in the frame. This also comes with a warning.
#'   \item outside the frame: coordinates more than `snap_max_m` from every
#'     unit;
#'   \item missing date or time step;
#'   \item outside the time window.
#' }
#'
#' Points that fall outside every polygon but within `snap_max_m` of one
#' are assigned to the nearest unit. Survey coordinates are often rounded,
#' and coastal sites fall just offshore. Snap distances are recorded in
#' `$surveys$.snap_m`.
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
#' @param date Name of a column of dates (`Date`, or `"YYYY-MM-DD"`
#'   strings). Give this or `t`.
#' @param t Name of a column of integer time steps (1..T). Required for a
#'   plain index time frame.
#' @param snap_max_m Maximum distance (metres) for snapping points outside
#'   every polygon to the nearest unit. Use 0 for no snapping.
#' @param verbose Print a summary of assignment and exclusions.
#'
#' @return An object of class `stareal_data`:
#'   \describe{
#'     \item{cells}{one row per observed cell: `cell_id` (=
#'       `(area - 1) * T + t`), `area` (unit index), `id`, `t`, `time`
#'       (label), `y`, `k`, `n_surveys`; ordered by `cell_id`}
#'     \item{surveys}{`data` with `.area`, `.t`, `.snap_m` and `.excluded`
#'       (`NA` if used, otherwise the reason) added}
#'     \item{frame, time}{the frame and time frame}
#'   }
#' @examples
#' fr <- areal_lattice(3, 3)
#' tf <- time_frame(T = 4)
#' surveys <- data.frame(unit = c("r1c1", "r1c1", "r2c3", "r9c9"),
#'                       step = c(1, 1, 4, 2), pos = c(3, 1, 0, 2), n = c(10, 5, 8, 9))
#' d <- areal_data(surveys, fr, tf, y = "pos", k = "n", area = "unit", t = "step")
#' d$cells
#' @export
areal_data <- function(data, frame, time, y, k, area = NULL, coords = NULL, crs = 4326,
                       date = NULL, t = NULL, snap_max_m = 10000, verbose = TRUE) {
  if (!inherits(frame, "stareal_frame")) stop("`frame` must be made by areal_frame()", call. = FALSE)
  if (!inherits(time, "stareal_time")) stop("`time` must be made by time_frame()", call. = FALSE)
  data <- as.data.frame(data)
  if (is.null(area) == is.null(coords)) stop("give exactly one of `area` and `coords`", call. = FALSE)
  if (is.null(date) == is.null(t)) stop("give exactly one of `date` and `t`", call. = FALSE)
  if (time$step == "index" && !is.null(date))
    stop("`time` is a plain index frame: give the time step as an integer column (`t`)", call. = FALSE)
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
  }

  # ---- counts ----
  yy <- suppressWarnings(as.numeric(data[[y]])); kk <- suppressWarnings(as.numeric(data[[k]]))
  flag(is.na(yy) | is.na(kk), "missing count")
  flag(yy != round(yy) | kk != round(kk) | yy < 0 | kk < 0, "invalid count (not a whole number >= 0)")
  flag(kk == 0, "none tested (k = 0)")
  flag(yy > kk, "positives exceed tested (y > k)")
  n_bad_counts <- sum(!is.na(why))

  # ---- area ----
  snap <- rep(NA_real_, N)
  if (!is.null(area)) {
    aid <- as.character(data[[area]])
    flag(is.na(aid) | aid == "", "missing area")
    a <- match(aid, frame$ids)
    n_unknown <- sum(is.na(a) & is.na(why))
    flag(is.na(a), "area not in frame")
  } else {
    if (is.null(frame$polygons)) stop("assigning coordinates needs a frame built from polygons", call. = FALSE)
    xy <- data[coords]
    flag(is.na(xy[[1]]) | is.na(xy[[2]]), "missing coordinates")
    res <- assign_points(xy, crs, frame, snap_max_m, ok = is.na(why))
    a <- res$area; snap <- res$snap_m
    n_unknown <- 0
    flag(is.na(a), sprintf("outside frame (> %g km from every unit)", snap_max_m / 1e3))
  }

  # ---- time ----
  if (!is.null(date)) {
    dd <- tryCatch(as.Date(data[[date]]), error = function(e)
      stop(sprintf("column `%s` cannot be read as dates", date), call. = FALSE))
    flag(is.na(dd), "missing date")
    tt <- time_index(time, dd)
  } else {
    tv <- suppressWarnings(as.numeric(data[[t]]))
    flag(is.na(tv), "missing time step")
    tt <- ifelse(!is.na(tv) & tv == round(tv) & tv >= 1 & tv <= time$T, as.integer(tv), NA_integer_)
  }
  flag(is.na(tt), "outside time window")

  # ---- pool ----
  use <- is.na(why)
  if (!any(use)) {
    ex <- table(why)
    stop("no survey could be used: ", paste(sprintf("%d %s", ex, names(ex)), collapse = "; "), call. = FALSE)
  }
  cid <- (a[use] - 1L) * time$T + tt[use]
  pooled_y <- rowsum(yy[use], cid); pooled_k <- rowsum(kk[use], cid)
  ns <- rowsum(rep(1L, sum(use)), cid)
  cell_id <- as.integer(rownames(pooled_y))                   # rowsum sorts the groups
  ca <- (cell_id - 1L) %/% time$T + 1L; ct <- (cell_id - 1L) %% time$T + 1L
  cells <- data.frame(cell_id = cell_id, area = ca, id = frame$ids[ca], t = ct, time = time$labels[ct],
                      y = as.integer(pooled_y[, 1]), k = as.integer(pooled_k[, 1]),
                      n_surveys = as.integer(ns[, 1]), stringsAsFactors = FALSE)

  surveys <- data
  surveys$.area <- a
  surveys$.t <- tt
  surveys$.snap_m <- snap
  surveys$.excluded <- why

  if (n_bad_counts) warning(sprintf("%d survey%s excluded for missing or invalid counts", n_bad_counts,
                                    if (n_bad_counts > 1) "s" else ""),
                            call. = FALSE)
  if (n_unknown) warning(sprintf("%d survey%s with area ids not in the frame", n_unknown, if (n_unknown > 1) "s" else ""), call. = FALSE)
  out <- structure(list(cells = cells, surveys = surveys, frame = frame, time = time), class = "stareal_data")
  if (verbose) print(out)
  out
}

# point-in-polygon on the plane, then a snap to the nearest unit (by
# distance in the frame's planar projection) for points within snap_max_m
assign_points <- function(xy, crs, frame, snap_max_m, ok) {
  N <- nrow(xy)
  area <- rep(NA_integer_, N); snap_m <- rep(NA_real_, N)
  idx <- which(ok)
  if (!length(idx)) return(list(area = area, snap_m = snap_m))
  with_s2_off({
    pts <- sf::st_as_sf(data.frame(x = as.numeric(xy[[1]][idx]), y = as.numeric(xy[[2]][idx])),
                        coords = c("x", "y"), crs = crs)
    pts <- sf::st_transform(pts, sf::st_crs(frame$polygons))
    hit <- suppressMessages(sf::st_within(pts, frame$polygons))   # planar on purpose (with_s2_off)
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
  list(area = area, snap_m = snap_m)
}

#' @export
print.stareal_data <- function(x, ...) {
  cl <- x$cells; s <- x$surveys
  cat(sprintf("<stareal_data> %d of %d surveys used, pooled into %d cells\n",
              sum(is.na(s$.excluded)), nrow(s), nrow(cl)))
  cat(sprintf("  %d of %d units and %d of %d time steps have data (%.2f%% of the %d area-time cells)\n",
              length(unique(cl$area)), x$frame$m, length(unique(cl$t)), x$time$T,
              100 * nrow(cl) / (x$frame$m * x$time$T), x$frame$m * x$time$T))
  cat(sprintf("  tested per cell: median %g, range %g-%g; total %g tested, %g positive\n",
              stats::median(cl$k), min(cl$k), max(cl$k), sum(cl$k), sum(cl$y)))
  snapped <- sum(!is.na(s$.snap_m) & s$.snap_m > 0 & is.na(s$.excluded))
  if (snapped) cat(sprintf("  %d survey%s snapped to the nearest unit (max %.1f km)\n", snapped,
                           if (snapped > 1) "s" else "",
                           max(s$.snap_m[is.na(s$.excluded)], na.rm = TRUE) / 1e3))
  ex <- table(s$.excluded)
  if (length(ex)) {
    cat("  excluded:\n")
    for (r in names(ex)) cat(sprintf("    %6d  %s\n", ex[[r]], r))
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
  }
  tk <- per(cl$t, x$time$T)
  graphics::barplot(tk, names.arg = x$time$labels, border = NA, col = "#2563eb",
                    main = "tested per time step", las = 2, cex.names = 0.6)
  invisible(x)
}
