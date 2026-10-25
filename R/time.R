# time.R
#
# The time frame: a regular grid of n_steps steps (months, quarters, years,
# or a plain index 1..n_steps). Africa_areal hard-wired months from January
# of a start year; here the step is a choice, and the labels travel with
# the outputs. (Inside the package's code the number of steps is written T,
# as in the derivation.)
#
# ------------------------------------------------------------------

#' Build a time frame
#'
#' The model's time axis is a regular grid of `n_steps` steps. A step is a
#' calendar month, quarter or year, running from the step containing
#' `start` to the step containing `end`. Alternatively, give only `n_steps`
#' for a plain index `1..n_steps`; survey data then carry the time step as
#' a whole number.
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
#' @param n_steps Instead of `start` and `end`: the number of steps of a
#'   plain index frame.
#' @return An object of class `stareal_time`: `n_steps`, `step`, `labels`
#'   (one per step, e.g. `"2020-03"`, `"2020-Q1"`, `"2020"`), and for
#'   calendar steps `starts`, the first day of each step.
#' @seealso [time_index()] to map dates to steps.
#' @examples
#' time_frame("1995-01", "2024-12", "month")
#' time_frame(2000, 2009, "year")
#' time_frame(n_steps = 12)
#' @export
time_frame <- function(start = NULL, end = NULL, step = c("month", "quarter", "year"), n_steps = NULL) {
  if (!is.null(n_steps)) {
    if (!is.null(start) || !is.null(end)) stop("give either `start` and `end`, or `n_steps`", call. = FALSE)
    n_steps <- check_whole(n_steps, "n_steps")
    return(structure(list(n_steps = n_steps, step = "index", labels = as.character(seq_len(n_steps)), starts = NULL),
                     class = "stareal_time"))
  }
  if (is.null(start) || is.null(end)) stop("give `start` and `end` (or `n_steps` for a plain index)", call. = FALSE)
  step <- match.arg(step)
  s <- step_floor(as_date_like(start, "start"), step)
  e <- step_floor(as_date_like(end, "end", at_end = TRUE), step)
  if (e < s) stop("`end` is before `start`", call. = FALSE)
  by <- c(month = "month", quarter = "3 months", year = "year")[[step]]
  starts <- seq(s, e, by = by)
  structure(list(n_steps = length(starts), step = step, labels = step_label(starts, step), starts = starts),
            class = "stareal_time")
}

#' Map dates to time steps
#'
#' @param time A `stareal_time` with calendar steps.
#' @param dates Dates: `Date`, date-times (read in their own time zone), or
#'   strings `"YYYY-MM-DD"`.
#' @return Integer step indices (1..n_steps), `NA` for missing or
#'   unreadable dates and for dates outside the window.
#' @examples
#' tf <- time_frame("2020-01", "2020-12")
#' time_index(tf, as.Date(c("2020-01-15", "2020-12-31", "2021-01-01")))
#' @export
time_index <- function(time, dates) {
  if (!inherits(time, "stareal_time")) stop("`time` must be made by time_frame()", call. = FALSE)
  if (time$step == "index") stop("an index time frame has no dates; give the time step as a whole number",
                                 call. = FALSE)
  d <- to_date(dates, "`dates`")
  i <- step_number(d, time$step) - step_number(time$starts[1], time$step) + 1L
  i[!is.na(i) & (i < 1L | i > time$n_steps)] <- NA_integer_
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

# Date from a Date, a whole year (4 digits), or "YYYY" / "YYYY-MM" /
# "YYYY-MM-DD". A year or month is its first day, or its last if at_end.
as_date_like <- function(x, what, at_end = FALSE) {
  bad <- function() stop(sprintf("`%s` (%s) is not a date: give a Date, a year such as 2020, or \"2020\", \"2020-06\" or \"2020-06-15\"",
                                 what, paste(format(x), collapse = ", ")), call. = FALSE)
  if (length(x) != 1 || is.na(x)) stop(sprintf("`%s` must be a single date", what), call. = FALSE)
  if (inherits(x, "Date")) return(x)
  if (inherits(x, "POSIXt")) return(to_date(x, what))
  period <- function(first, by) if (at_end) seq(first, by = by, length.out = 2)[2] - 1 else first
  if (is.numeric(x)) {
    if (x != round(x) || x < 1000 || x > 9999) bad()
    return(period(as.Date(sprintf("%d-01-01", as.integer(x))), "year"))
  }
  x <- trimws(as.character(x))
  d <- tryCatch({
    if (grepl("^[0-9]{4}$", x)) period(as.Date(paste0(x, "-01-01")), "year")
    else if (grepl("^[0-9]{4}-[0-9]{1,2}$", x)) period(as.Date(paste0(x, "-01")), "month")
    else if (grepl("^[0-9]{4}-[0-9]{1,2}-[0-9]{1,2}$", x)) as.Date(x)
    else NA
  }, error = function(e) NA)
  if (is.na(d)) bad()
  d
}

#' @export
print.stareal_time <- function(x, ...) {
  if (x$step == "index") {
    cat(sprintf("<stareal_time> plain index, %d steps\n", x$n_steps))
  } else {
    cat(sprintf("<stareal_time> %s steps, %s to %s (%d steps)\n", paste0(x$step, "ly"),
                x$labels[1], x$labels[x$n_steps], x$n_steps))
  }
  invisible(x)
}
