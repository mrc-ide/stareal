# time.R
#
# The time frame: a regular grid of T steps (months, quarters, years, or a
# plain index 1..T). Africa_areal hard-wired months from January of a start
# year; here the step is a choice, and the labels travel with the outputs.
#
# ------------------------------------------------------------------

#' Build a time frame
#'
#' The model's time axis is a regular grid of `T` steps. A step is a
#' calendar month, quarter or year, running from the step containing
#' `start` to the step containing `end`. Alternatively, give only `T` for a
#' plain index `1..T`. Survey data must then carry the time step as an
#' integer.
#'
#' The temporal prior is a first-order autoregression between consecutive
#' steps, so the step sets the scale on which the temporal correlation
#' \eqn{\rho_t} is read.
#'
#' @param start,end The first and last dates of the window: `Date`s, or
#'   strings `"YYYY"`, `"YYYY-MM"` or `"YYYY-MM-DD"`, or whole years as
#'   numbers. A year or month stands for the whole period: as `start` its
#'   first day, as `end` its last (so `time_frame(2001, 2020, "month")`
#'   runs from January 2001 to December 2020). The window then covers every
#'   step that contains a day between `start` and `end`.
#' @param step `"month"`, `"quarter"` or `"year"`.
#' @param T Instead of `start` and `end`: the number of steps of a plain
#'   index frame.
#' @return An object of class `stareal_time`: `T`, `step`, `labels` (one per
#'   step, e.g. `"2020-03"`, `"2020-Q1"`, `"2020"`), and for calendar steps
#'   `starts`, the first day of each step.
#' @seealso [time_index()] to map dates to steps.
#' @examples
#' time_frame("1995-01", "2024-12", "month")
#' time_frame(2000, 2009, "year")
#' time_frame(T = 12)
#' @export
time_frame <- function(start = NULL, end = NULL, step = c("month", "quarter", "year"), T = NULL) {
  if (!is.null(T)) {
    if (!is.null(start) || !is.null(end)) stop("give either `start` and `end`, or `T`", call. = FALSE)
    T <- as.integer(T)
    if (length(T) != 1 || is.na(T) || T < 1) stop("`T` must be a positive whole number", call. = FALSE)
    return(structure(list(T = T, step = "index", labels = as.character(seq_len(T)), starts = NULL),
                     class = "stareal_time"))
  }
  if (is.null(start) || is.null(end)) stop("give `start` and `end` (or `T` for a plain index)", call. = FALSE)
  step <- match.arg(step)
  s <- step_floor(as_date_like(start, "start"), step)
  e <- step_floor(as_date_like(end, "end", at_end = TRUE), step)
  if (e < s) stop("`end` is before `start`", call. = FALSE)
  by <- c(month = "month", quarter = "3 months", year = "year")[[step]]
  starts <- seq(s, e, by = by)
  structure(list(T = length(starts), step = step, labels = step_label(starts, step), starts = starts),
            class = "stareal_time")
}

#' Map dates to time steps
#'
#' @param time A `stareal_time` with calendar steps.
#' @param dates Dates (`Date`, or strings `"YYYY-MM-DD"`).
#' @return Integer step indices (1..T), `NA` for missing dates or dates
#'   outside the window.
#' @examples
#' tf <- time_frame("2020-01", "2020-12")
#' time_index(tf, as.Date(c("2020-01-15", "2020-12-31", "2021-01-01")))
#' @export
time_index <- function(time, dates) {
  if (time$step == "index") stop("an index time frame has no dates; give the time step as an integer",
                                 call. = FALSE)
  d <- as.Date(dates)
  i <- step_number(d, time$step) - step_number(time$starts[1], time$step) + 1L
  i[!is.na(i) & (i < 1L | i > time$T)] <- NA_integer_
  as.integer(i)
}

# a running count of steps since year 0, so that differences give indices
step_number <- function(d, step) {
  y <- as.integer(format(d, "%Y")); mo <- as.integer(format(d, "%m"))
  switch(step, month = 12L * y + mo - 1L, quarter = 4L * y + (mo - 1L) %/% 3L, year = y)
}

step_floor <- function(d, step) {
  y <- format(d, "%Y"); mo <- as.integer(format(d, "%m"))
  switch(step,
         month = as.Date(sprintf("%s-%02d-01", y, mo)),
         quarter = as.Date(sprintf("%s-%02d-01", y, 3L * ((mo - 1L) %/% 3L) + 1L)),
         year = as.Date(sprintf("%s-01-01", y)))
}

step_label <- function(starts, step) {
  switch(step,
         month = format(starts, "%Y-%m"),
         quarter = sprintf("%s-Q%d", format(starts, "%Y"), (as.integer(format(starts, "%m")) - 1L) %/% 3L + 1L),
         year = format(starts, "%Y"))
}

# Date from a Date, a whole year, or "YYYY" / "YYYY-MM" / "YYYY-MM-DD". A
# year or month is its first day, or its last if at_end.
as_date_like <- function(x, what, at_end = FALSE) {
  if (length(x) != 1 || is.na(x)) stop(sprintf("`%s` must be a single date", what), call. = FALSE)
  if (inherits(x, "Date")) return(x)
  if (inherits(x, "POSIXt")) return(as.Date(x))
  period <- function(first, by) if (at_end) seq(first, by = by, length.out = 2)[2] - 1 else first
  if (is.numeric(x) && x == round(x)) return(period(as.Date(sprintf("%d-01-01", as.integer(x))), "year"))
  x <- as.character(x)
  if (grepl("^[0-9]{4}$", x)) return(period(as.Date(paste0(x, "-01-01")), "year"))
  if (grepl("^[0-9]{4}-[0-9]{1,2}$", x)) return(period(as.Date(paste0(x, "-01")), "month"))
  d <- tryCatch(as.Date(x), error = function(e) NA)
  if (is.na(d)) stop(sprintf("`%s` (\"%s\") is not a date", what, x), call. = FALSE)
  d
}

#' @export
print.stareal_time <- function(x, ...) {
  if (x$step == "index") {
    cat(sprintf("<stareal_time> plain index, T = %d steps\n", x$T))
  } else {
    cat(sprintf("<stareal_time> %s steps, %s to %s (T = %d)\n", paste0(x$step, "ly"),
                x$labels[1], x$labels[x$T], x$T))
  }
  invisible(x)
}
