# workers.R
#
# Running independent jobs (MCMC chains, blocks of prediction) on worker
# processes inside the R session. Internal.
#
# run_jobs(X, FUN, shared) is lapply(X, FUN, shared) on socket workers
# (parallel::makePSOCKcluster), or in this session when one worker is asked
# for. Why it is written out rather than parallel::clusterApplyLB():
# - progress: the session polls the workers (socketSelect with a timeout)
#   and relays their messages from a log file, so progress lines appear in
#   RStudio and Rgui, where a worker's console output does not;
# - interrupts and errors: if the session stops waiting for any reason
#   (Esc, an error in another job, a worker that died), every worker is
#   killed by process id. stopCluster() alone leaves a busy worker running
#   to the end of its job (a chain could run on for hours, writing its
#   checkpoints);
# - errors come back as the job's own message, with the job's index;
# - results can be consumed as they arrive (`collect`), so the session never
#   holds every job's result at once.
# Why worker processes and not forks: on macOS, a forked worker that calls
# Accelerate's BLAS after the parent has used it dies, and setting
# VECLIB_MAXIMUM_THREADS from inside R does not help once BLAS has run
# (checked on this Mac, 2026-10-22); OpenBLAS has the same problem on
# Linux. Fresh processes, started with the BLAS thread limits already in
# their environment, are safe everywhere.
#
# Sending and receiving use parallel's own (unexported) sendCall() and
# recvResult(), the functions clusterApplyLB() is built on.
#
# ------------------------------------------------------------------

# where a worker keeps the objects shared by all its jobs
.worker <- new.env(parent = emptyenv())

# FUN(x, shared) for each element of X; returns the list of results, or of
# collect(i, result) when `collect` is given (called as each job finishes,
# in completion order). A job that fails stops everything with an error of
# class "stareal_job_error" (fields job, message), after the other workers
# have been stopped.
run_jobs <- function(X, FUN, shared = list(), n_workers = 1L, blas_threads = 1L, verbose = TRUE,
                     collect = NULL) {
  n <- length(X)
  out <- vector("list", n)
  keep <- function(i, v) out[[i]] <<- if (is.null(collect)) v else collect(i, v)
  n_workers <- min(as.integer(n_workers), n)
  if (n_workers < 2) {
    for (i in seq_len(n)) {
      v <- tryCatch(FUN(X[[i]], shared), error = function(e) job_error(i, conditionMessage(e)))
      keep(i, v)
    }
    return(out)
  }

  vars <- c("VECLIB_MAXIMUM_THREADS", "MKL_NUM_THREADS", "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS")
  old <- Sys.getenv(vars, unset = NA)
  on.exit({
    for (v in vars) if (is.na(old[[v]])) Sys.unsetenv(v) else do.call(Sys.setenv, as.list(old[v]))
  }, add = TRUE)
  do.call(Sys.setenv, as.list(stats::setNames(rep(as.character(blas_threads), length(vars)), vars)))
  log <- tempfile("stareal_workers_", fileext = ".log")
  file.create(log)
  cl <- parallel::makePSOCKcluster(n_workers, outfile = log)
  pids <- unlist(parallel::clusterCall(cl, Sys.getpid))
  finished <- FALSE
  on.exit({
    if (!finished) kill_processes(pids)
    try(parallel::stopCluster(cl), silent = TRUE)
    unlink(log)
  }, add = TRUE)

  # the functions sent to workers must not carry this namespace with them
  # before the workers have loaded it, hence environment = globalenv()
  load_pkg <- function(dev_path, lib_paths) {
    .libPaths(lib_paths)
    if (is.null(dev_path)) loadNamespace("stareal", lib.loc = lib_paths) else
      pkgload::load_all(dev_path, quiet = TRUE, helpers = FALSE, attach_testthat = FALSE)
    NULL
  }
  put_shared <- function(sh) {
    assign("shared", sh, envir = utils::getFromNamespace(".worker", "stareal"))
    NULL
  }
  wrap <- function(x, fun) {
    sh <- get("shared", envir = utils::getFromNamespace(".worker", "stareal"))
    tryCatch(list(value = fun(x, sh)), error = function(e) list(error = conditionMessage(e)))
  }
  environment(load_pkg) <- environment(put_shared) <- environment(wrap) <- environment(FUN) <- globalenv()
  parallel::clusterCall(cl, load_pkg, dev_path(), lib_paths())
  parallel::clusterCall(cl, put_shared, shared)

  send <- utils::getFromNamespace("sendCall", "parallel")
  recv <- utils::getFromNamespace("recvResult", "parallel")
  busy <- rep(NA_integer_, n_workers)                    # the job each worker is running
  next_job <- 1L
  start <- function(w) {
    send(cl[[w]], wrap, list(X[[next_job]], FUN))
    busy[w] <<- next_job
    next_job <<- next_job + 1L
  }
  for (w in seq_len(n_workers)) start(w)
  relay <- log_relay(log, verbose)
  n_done <- 0L
  while (n_done < n) {
    relay()
    active <- which(!is.na(busy))
    ready <- socketSelect(lapply(cl[active], `[[`, "con"), timeout = 0.5)
    for (w in active[ready]) {
      i <- busy[w]
      r <- tryCatch(recv(cl[[w]]), error = function(e) NULL)
      relay()
      if (is.null(r)) job_error(i, "its worker process died (out of memory?)")
      if (inherits(r, "try-error")) job_error(i, as.character(r))
      if (!is.null(r$error)) job_error(i, r$error)
      keep(i, r$value)
      n_done <- n_done + 1L
      busy[w] <- NA_integer_
      if (next_job <= n) start(w)
    }
  }
  relay()
  finished <- TRUE
  out
}

job_error <- function(i, message) {
  stop(structure(class = c("stareal_job_error", "error", "condition"),
                 list(message = message, call = NULL, job = i)))
}

# a function that copies new complete lines of the worker log to this
# session's messages (dropping the cluster's own start-up lines)
log_relay <- function(log, verbose) {
  pos <- 0
  partial <- ""
  function() {
    if (!verbose) return(invisible())
    size <- file.size(log)
    if (is.na(size) || size <= pos) return(invisible())
    con <- file(log, "rb")
    on.exit(close(con))
    seek(con, pos)
    txt <- paste0(partial, readChar(con, size - pos, useBytes = TRUE))
    pos <<- size
    lines <- strsplit(txt, "\n", fixed = TRUE)[[1]]
    if (!endsWith(txt, "\n")) { partial <<- lines[length(lines)]; lines <- lines[-length(lines)] } else partial <<- ""
    lines <- lines[nzchar(trimws(lines)) & !grepl("^starting worker pid=", lines)]
    for (l in lines) message(l)
    invisible()
  }
}

kill_processes <- function(pids) {
  for (p in pids) try(tools::pskill(p, tools::SIGTERM), silent = TRUE)
}

# the source folder when this package was loaded by pkgload::load_all()
# (as when testing), else NULL (installed)
dev_path <- function() {
  if (!requireNamespace("pkgload", quietly = TRUE) || !pkgload::is_dev_package("stareal")) return(NULL)
  getNamespaceInfo(asNamespace("stareal"), "path")
}

# library paths for workers: the one this session's stareal came from first
lib_paths <- function() {
  own <- if (is.null(dev_path())) dirname(getNamespaceInfo(asNamespace("stareal"), "path")) else NULL
  unique(c(own, .libPaths()))
}
