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

#' Plot traces for observed cells posterior
#'
#' @param fit_object A stareal fit object
#'
#' @returns Trace plots of the posterior for all observed cells
#' @examples
#' \donttest{
#' fr<- areal_lattice(4, 4)
#' tf<- time_frame(n_steps = 5)
#' sim<- simulate_areal(fr, tf, n_cells = 40, seed = 1)
#' fit<- fit_areal(sim, control = list(n_iter = 100, n_burn = 50, n_store = 25, n_chains = 1,
#'                                      ess_sweeps = 2, verbose = FALSE))
#' plot_cells_trace(fit)
#' }
#' @export
plot_cells_trace<- function(fit_object){
  #fit_object$draws$x<- plogis(fit_object$draws$x)
  n_iter<- dim(fit_object$draws$x)[1]
  n_chains<- dim(fit_object$draws$x)[2]
  n_cells<- dim(fit_object$draws$x)[3]

  set.seed(0);cols<- randomcoloR::randomColor(n_chains)
  par(mfrow=c(3,3))
  for(i in 1:n_cells){
    if(n_chains==1){
      plot(fit_object$draws$x[,1,i], type="l",xlab ="MCMC iterations", main=paste("cell",i), ylab = "", col = adjustcolor("red", alpha.f = 0.5), ylim = range(as.numeric(fit_object$draws$x[,,i])))
      grid()
    }else{
      plot(fit_object$draws$x[,1,i], type="l",xlab ="MCMC iterations", main=paste("cell",i), ylab = "", col = adjustcolor("red", alpha.f = 0.5), ylim = range(as.numeric(fit_object$draws$x[,,i])))
      grid()
      for(j in 2:n_chains){
        lines(1:n_iter, fit_object$draws$x[,j,i], col = adjustcolor(cols[j], alpha.f = 0.5))
      }
    }
  }
}

#' Conditional loo given observed cell posterior
#'
#' @param fit A stareal fit object
#'
#' @returns conditional loo
#' @examples
#' \donttest{
#' fr<- areal_lattice(4, 4)
#' tf<- time_frame(n_steps = 5)
#' sim<- simulate_areal(fr, tf, n_cells = 40, seed = 1)
#' fit<- fit_areal(sim, control = list(n_iter = 100, n_burn = 50, n_store = 25, n_chains = 1,
#'                                      ess_sweeps = 2, verbose = FALSE))
#' conditional_loo_fit(fit)
#' }
#' @export
conditional_loo_fit<- function(fit){
  x<- fit$draws$x
  d<- dim(x)
  cells<- fit$data$cells
  y<- fit$data$cells$y
  N<- fit$data$cells$k

  p<- plogis(x)
  log_lik<- array(NA_real_, dim = d)

  for(j in seq_len(d[3])) {
    log_lik[, , j] <- dbinom(y[j], size = N[j], prob = p[, , j], log = TRUE)
  }
  # Convert to draws x observations
  log_lik<- matrix(log_lik, nrow = d[1] * d[2], ncol = d[3])
  colnames(log_lik)<- cells$cell_id
  looCV<- loo::loo(log_lik)
  return(looCV)
}

#' Conditional waic given observed cell posterior
#'
#' @param fit A stareal fit object
#'
#' @returns conditional waic
#' @examples
#' \donttest{
#' fr<- areal_lattice(4, 4)
#' tf<- time_frame(n_steps = 5)
#' sim<- simulate_areal(fr, tf, n_cells = 40, seed = 1)
#' fit<- fit_areal(sim, control = list(n_iter = 100, n_burn = 50, n_store = 25, n_chains = 1,
#'                                      ess_sweeps = 2, verbose = FALSE))
#' conditional_waic_fit(fit)
#' }
#' @export
conditional_waic_fit<- function(fit){
  x<- fit$draws$x
  d<- dim(x)
  cells<- fit$data$cells
  y<- fit$data$cells$y
  N<- fit$data$cells$k
  p<- plogis(x)
  log_lik<- array(NA_real_, dim = d)

  for (j in seq_len(d[3])) {
    log_lik[, , j] <- dbinom(y[j], size = N[j], prob = p[, , j], log = TRUE)
  }
  # Convert to draws x observations
  log_lik<- matrix(log_lik, nrow = d[1] * d[2], ncol = d[3])
  colnames(log_lik)<- cells$cell_id
  WAIC<- loo::waic(log_lik)
  return(WAIC)
}


#' Create stareal_data from positive and tested matrices
#'
#' Converts matrices of positive and tested counts into pooled
#' area-time observations.
#'
#' Rows represent spatial units and columns represent time steps.
#'
#' @param positive_matrix Matrix of positive counts.
#' @param tested_matrix Matrix of numbers tested.
#' @param frame A stareal_frame.
#' @param time A stareal_time.
#'
#' @return An object of class stareal_data.
#' @export
stareal_data_from_matrix <- function(positive_matrix,tested_matrix,frame,time){
  nrows<- nrow(positive_matrix)
  ncols<- ncol(positive_matrix)
  if(nrows != length(frame$ids))
    stop("Number of matrix rows must equal number of spatial units", call. = FALSE)
  if (ncols != time$n_steps)
    stop("Number of matrix columns must equal number of time steps",
         call. = FALSE)

  yy<- as_num(as.vector(t(positive_matrix)))
  kk<- as_num(as.vector(t(tested_matrix)))
  cells<- expand.grid(area = seq_len(nrows),  t = seq_len(ncols))
  cells$cell_id<- (cells$area - 1L) * ncols + cells$t
  cells$id <- frame$ids[cells$area]
  cells$time <- time$labels[cells$t]
  cells$y<- yy
  cells$k<- kk

  valid<- is.finite(cells$y) &
    is.finite(cells$k) &
    cells$y >= 0 &
    cells$k > 0 &
    cells$y <= cells$k &
    cells$y == round(cells$y) &
    cells$k == round(cells$k)
  cells<- cells[valid, , drop = FALSE]

  if(!nrow(cells))
    stop("No valid area-time observations", call. = FALSE)

  cells<- cells[, c("cell_id","area","id","t","time", "y", "k")]
  cells$n_surveys<- NA_integer_
  structure(list(cells = cells, surveys = NULL,
      frame = frame, time = time, truth = NULL),
    class = "stareal_data")
}
