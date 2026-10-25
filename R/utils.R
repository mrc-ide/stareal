# utils.R
#
# Small helpers shared across the package. Internal.
#
# ------------------------------------------------------------------

# Evaluate `code` with the random-number generator set to `seed`, then put
# the caller's state back, so package functions given a seed never move the
# user's stream. seed = NULL: evaluate with the current stream.
# The generator is pinned to R's defaults (Mersenne-Twister, Inversion,
# Rejection) whatever kind the session uses, so a seed gives the same
# numbers in this session and in worker processes. The caller's kind is
# restored with its state.
with_seed <- function(seed, code) {
  if (is.null(seed)) return(code)
  old <- get0(".Random.seed", envir = globalenv(), inherits = FALSE)
  old_kind <- RNGkind()
  on.exit({
    if (is.null(old)) {
      if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) rm(".Random.seed", envir = globalenv())
      do.call(RNGkind, as.list(old_kind))
    } else assign(".Random.seed", old, envir = globalenv())
  })
  set_seed(seed)
  code
}

# set.seed() with R's default generators, whatever the session's RNGkind
set_seed <- function(seed) {
  set.seed(seed, kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
}

# Evaluate `code` with sf's spherical geometry (s2) switched off.
# Administrative boundaries (e.g. GADM) are usually valid on the plane but
# not on the sphere, and adjacency and point-in-polygon are planar questions
# at this scale.
with_s2_off <- function(code) {
  old <- sf::sf_use_s2()
  on.exit(suppressMessages(sf::sf_use_s2(old)))
  suppressMessages(sf::sf_use_s2(FALSE))
  code
}

need_pkg <- function(pkg, why) {
  if (!requireNamespace(pkg, quietly = TRUE))
    stop(sprintf("package '%s' is needed %s; install it with install.packages(\"%s\")", pkg, why, pkg),
         call. = FALSE)
}

say <- function(verbose, ...) if (verbose) message(sprintf(...))

# Unit ids as character, the same string whatever the storage type: a
# whole-number id is "100000", never "1e+05", whether it came as an integer
# or a double. NA stays NA.
id_chr <- function(x) {
  if (is.factor(x)) x <- as.character(x)
  out <- if (is.numeric(x)) format(x, scientific = FALSE, trim = TRUE, digits = 15) else as.character(x)
  out[is.na(x)] <- NA_character_
  out
}

# numbers from a column that may be a factor (its labels, not its codes) or
# text; anything unreadable becomes NA
as_num <- function(v) suppressWarnings(as.numeric(if (is.factor(v)) as.character(v) else v))

# Dates from Date, POSIXct/lt (in their own time zone, so a local-midnight
# timestamp east of UTC stays on its day) or "YYYY-MM-DD" text. Numbers are
# refused: they are usually years, or days since 1970 by accident.
# Unreadable text becomes NA.
to_date <- function(x, what) {
  if (is.numeric(x))
    stop(sprintf(paste("%s is numeric; give dates as Date or \"YYYY-MM-DD\" text (a year is not a date:",
                       "build one with as.Date(sprintf(\"%%d-%%02d-01\", year, month)))"), what), call. = FALSE)
  if (inherits(x, "POSIXt")) {
    tz <- attr(x, "tzone")
    return(as.Date(format(x, "%Y-%m-%d", tz = if (is.null(tz) || !nzchar(tz[1])) "" else tz[1])))
  }
  if (inherits(x, "Date")) return(x)
  x <- as.character(x)
  out <- as.Date(x, optional = TRUE, tryFormats = "%Y-%m-%d")
  out[!grepl("^\\s*[0-9]{4}-[0-9]{1,2}-[0-9]{1,2}", x)] <- NA          # "2020-1", "20200101": not dates
  out
}

# a symmetric 0/1 adjacency as a general sparse matrix (dgCMatrix)
as_dgC <- function(W) {
  if (!methods::is(W, "Matrix")) W <- Matrix::Matrix(as.matrix(W), sparse = TRUE)
  W <- methods::as(methods::as(methods::as(W, "dMatrix"), "generalMatrix"), "CsparseMatrix")
  Matrix::drop0(W)
}

# neighbour lists from a sparse symmetric adjacency (column i's row indices)
adjacency_list <- function(W) {
  m <- ncol(W)
  split(W@i + 1L, factor(rep(seq_len(m), diff(W@p)), levels = seq_len(m)))
}

# connected components of a graph given as neighbour lists: an integer
# component label per node, 1 = the largest component
graph_components <- function(nb) {
  m <- length(nb); comp <- integer(m); n_comp <- 0L
  for (s in seq_len(m)) {
    if (comp[s]) next
    n_comp <- n_comp + 1L; comp[s] <- n_comp; queue <- s
    while (length(queue)) {
      v <- nb[[queue[1]]]; queue <- queue[-1]
      new <- v[comp[v] == 0L]
      comp[new] <- n_comp; queue <- c(queue, new)
    }
  }
  # relabel by size, largest first (ties: first seen)
  sizes <- tabulate(comp, n_comp)
  match(comp, order(-sizes, seq_len(n_comp)))
}

# Default number of worker processes: the cores less one, at most `cap`
# (the jobs there are to run), and 2 under R CMD check's core limit.
default_cores <- function(cap = Inf) {
  n <- suppressWarnings(parallel::detectCores())
  if (is.na(n)) n <- 1L
  n <- max(1L, n - 1L)
  if (nzchar(Sys.getenv("_R_CHECK_LIMIT_CORES_"))) n <- min(n, 2L)
  as.integer(max(1L, min(n, cap)))
}

# check one whole number >= min (NULL allowed if null_ok)
check_whole <- function(x, nm, min = 1, null_ok = FALSE) {
  if (null_ok && is.null(x)) return(NULL)
  if (!is.numeric(x) || length(x) != 1 || !is.finite(x) || x != round(x) || x < min)
    stop(sprintf("`%s` must be a whole number >= %d", nm, min), call. = FALSE)
  if (x > .Machine$integer.max) stop(sprintf("`%s` is too large", nm), call. = FALSE)
  as.integer(x)
}

# stop on arguments swallowed by `...` (misspelt names)
check_dots <- function(...) {
  if (...length()) {
    nm <- names(list(...))
    stop(sprintf("unused argument(s): %s", paste(ifelse(is.null(nm) | nm == "", "<unnamed>", nm), collapse = ", ")),
         call. = FALSE)
  }
}
