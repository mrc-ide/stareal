# utils.R
#
# Small helpers shared across the package. Internal.
#
# ------------------------------------------------------------------

# Evaluate `code` with the random-number generator set to `seed`, then put
# the caller's state back, so package functions given a seed never move the
# user's stream. seed = NULL: evaluate with the current stream.
with_seed <- function(seed, code) {
  if (is.null(seed)) return(code)
  old <- get0(".Random.seed", envir = globalenv(), inherits = FALSE)
  on.exit(if (is.null(old)) rm(".Random.seed", envir = globalenv())
          else assign(".Random.seed", old, envir = globalenv()))
  set.seed(seed)
  code
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
