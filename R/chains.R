# chains.R
#
# Running MCMC chains: one chain (run_chain, with checkpoint and resume) and
# all chains of a fit (run_chains: one chain per worker process, in
# parallel, inside the R session). The moves themselves are in sampler.R.
# Internal.
#
# run_chain() is ported from Africa_areal/R/chains.R (2026-10-12). Changes:
# - `setup` (fit_setup(), below) carries the likelihood as one object, so
#   lik$loglik, lik$pseudo and lik$pseudo_var stand where ll_fun, g_star and
#   S_aux stood; `profile` is now `control` (mcmc_control());
# - any number of fixed effects (the Africa code stored two: intercept and
#   trend);
# - no heartbeat file (there is no supervisor: chains that outlive the R
#   session were dropped, Bob 2026-10-22);
# - the state is also checkpointed at the last iteration, so a finished
#   chain is read back, not rerun;
# - progress messages say how long is left, and the first comes after 10
#   iterations so the cost of a fit is known early (plan section 6). Messages
#   use no random numbers, so the chain is unchanged.
# The draws are those of the Africa code for the same seed.
#
# Why worker processes and not forks: on macOS, a forked worker that calls
# Accelerate's BLAS after the parent has already used it dies, and setting
# VECLIB_MAXIMUM_THREADS from inside R does not help once BLAS has run
# (checked on this Mac, 2026-10-22). OpenBLAS has the same problem on
# Linux. Fresh worker processes, started with the BLAS thread limits already
# in their environment, are safe everywhere.
#
# ------------------------------------------------------------------

# ---- everything the sampler needs about one data set -------------------------

# cells: areal_data()$cells (area, t, y, k). eig: the frame's
# eigendecomposition. X: the fixed-effect design of the observed cells
# (intercept + trend if NULL). Returns the fixed inputs of the MCMC.
fit_setup <- function(cells, eig, T, priors, nugget, X = NULL) {
  if (is.null(X)) X <- design_X(cells$t, T)
  if (is.null(colnames(X))) colnames(X) <- if (ncol(X) == 2) c("intercept", "trend") else paste0("beta", seq_len(ncol(X)))
  obs <- make_obs(cells$area, cells$t, T, X = X)
  p <- ncol(X)
  b0 <- rep_len(priors$beta_mean, p); B <- diag(rep_len(priors$beta_sd, p)^2, p)
  list(obs = obs, eig = eig, lik = lik_binomial(cells$y, cells$k), b0 = b0, B = B,
       mu_o = as.vector(obs$X %*% b0), nugget = isTRUE(nugget), log_prior = make_log_prior(priors),
       beta_names = colnames(X))
}

# ---- one chain --------------------------------------------------------------

# setup: fit_setup() output. control: mcmc_control() output. Returns
#   theta       data frame: iteration, theta (unconstrained), lp, ll, and the
#               acceptance of the centred and surrogate moves (mean over the
#               thinning window), every theta_thin-th iteration, burn-in
#               included
#   x           n_store x n latent draws (evenly spaced after burn-in)
#   draw_theta  n_store x d theta at those same iterations
#   beta        n_store x p fixed-effect draws, one from beta | x_o, theta
#               per stored latent draw
#   stored_iter, seconds, final adaptation states
#
# checkpoint: a file path. Every control$checkpoint_every iterations the
#   whole state (latent values, theta, adaptation, outputs so far and the
#   random-number state) is saved there; if the file exists at the start,
#   the chain continues from it and produces exactly the chain an
#   uninterrupted run would have (tests/testthat/test-fit.R).
# stop_after: stop at this iteration, after saving a checkpoint, and return
#   NULL (for the tests).
run_chain <- function(chain, setup, control, checkpoint = NULL, stop_after = NULL) {
  obs <- setup$obs; eig <- setup$eig; B <- setup$B; mu_o <- setup$mu_o
  lik <- setup$lik; prior <- setup$log_prior
  S_aux <- lik$pseudo_var; g_star <- lik$pseudo; ll_fun <- lik$loglik
  n_iter <- control$n_iter; n_burn <- control$n_burn; p_beta <- ncol(obs$X)
  step <- (n_iter - n_burn) %/% control$n_store
  store_at <- n_burn + step * seq_len(control$n_store)
  every <- if (n_iter >= 5000) 500L else max(10L, n_iter %/% 10L)
  ck_every <- as.integer(control$checkpoint_every)
  XBX <- obs$X %*% B %*% t(obs$X)

  if (!is.null(checkpoint) && file.exists(checkpoint)) {
    ck <- readRDS(checkpoint)                             # resume
    list2env(ck$state, environment())
    assign(".Random.seed", ck$seed, envir = globalenv())
    it0 <- ck$it; elapsed0 <- ck$elapsed
    if (control$progress && it0 < n_iter)
      message(sprintf("chain %d: resuming from the checkpoint at iteration %d", chain, it0))
  } else {
    set.seed(control$seed * 1000 + chain)
    # dispersed starting values; latents at the empirical logits
    theta <- c(log(1), qlogis(0.9), qlogis(0.95)) + rnorm(3, 0, 0.1)
    if (setup$nugget) theta <- c(theta, log(0.5) + rnorm(1, 0, 0.1))
    x <- g_star
    ll <- ll_fun(x)
    # starting proposal sds on the unconstrained scale (roughly the posterior
    # sds of the first pilots; the adaptation corrects them). The two moves
    # adapt separately: their acceptance landscapes differ.
    sd0 <- c(0.03, 0.4, 0.12, 0.1)[seq_along(theta)]
    am_c <- am_init(theta, sd0); am_s <- am_init(theta, sd0)
    trace_at <- seq(control$theta_thin, n_iter, by = control$theta_thin)
    tr <- matrix(NA_real_, length(trace_at), length(theta) + 5,
                 dimnames = list(NULL, c("iteration", theta_names(setup$nugget), "lp", "ll",
                                         "accept_centred", "accept_surrogate")))
    out_x <- matrix(NA_real_, control$n_store, obs$n)
    out_beta <- matrix(NA_real_, control$n_store, p_beta, dimnames = list(NULL, setup$beta_names))
    out_theta <- matrix(NA_real_, control$n_store, length(theta),
                        dimnames = list(NULL, theta_names(setup$nugget)))
    acc_c_win <- acc_s_win <- 0; j <- 0; jt <- 0
    acc_c_rep <- acc_s_rep <- 0
    it0 <- 0; elapsed0 <- 0
  }
  # derived from the state, so not checkpointed (and recomputed identically)
  cur <- list(theta = theta, sig = build_sigma_oo(theta, obs, eig, B))
  cur$lp <- prior(theta) + log_dmvnorm_chol(x, mu_o, cur$sig$L, cur$sig$logdet)
  prop_c <- am_prop_chol(am_c); prop_s <- am_prop_chol(am_s)
  t0 <- proc.time()[["elapsed"]]
  # the iteration of the last progress message (for the acceptance rates)
  last_rep <- if (it0 >= every) (it0 %/% every) * every else if (it0 >= 10L && every > 10L) 10L else 0L
  save_checkpoint <- function(it) {
    state <- list(theta = cur$theta, x = x, ll = ll, am_c = am_c, am_s = am_s, tr = tr, out_x = out_x,
                  out_beta = out_beta, out_theta = out_theta, acc_c_win = acc_c_win,
                  acc_s_win = acc_s_win, j = j, jt = jt, acc_c_rep = acc_c_rep, acc_s_rep = acc_s_rep)
    tmp <- paste0(checkpoint, ".tmp")                    # write, then rename: never a half-written file
    saveRDS(list(it = it, elapsed = elapsed0 + proc.time()[["elapsed"]] - t0, state = state,
                 seed = get(".Random.seed", envir = globalenv())), tmp, compress = FALSE)
    file.rename(tmp, checkpoint)
  }

  for (it in seq.int(it0 + 1L, length.out = max(0L, n_iter - it0))) {
    # 1. latents given theta
    cur$sig <- ess_surrogate_prep(cur$sig, mu_o, g_star, S_aux)
    for (sweep in seq_len(control$ess_sweeps)) {
      st <- ess_step_surrogate(x, ll, cur$sig, g_star, S_aux, ll_fun); x <- st$x; ll <- st$ll
    }
    # 2. theta: centred move, then surrogate-data move, theta_moves times
    acc_c <- acc_s <- 0
    for (mv in seq_len(control$theta_moves)) {
      cur$lp <- prior(cur$theta) + log_dmvnorm_chol(x, mu_o, cur$sig$L, cur$sig$logdet)
      cur <- mh_step_theta(cur, x, mu_o, obs, eig, B, prop_c, prior)
      acc_c <- acc_c + cur$accepted
      if (it <= n_burn) { am_c <- am_update(am_c, cur$theta, cur$alpha); prop_c <- am_prop_chol(am_c) }
      st <- mh_step_theta_surrogate(cur, x, ll, mu_o, obs, eig, B, S_aux, prop_s, prior, ll_fun)
      cur <- st$cur; x <- st$x; ll <- st$ll
      acc_s <- acc_s + st$accepted
      if (it <= n_burn) { am_s <- am_update(am_s, cur$theta, st$alpha); prop_s <- am_prop_chol(am_s) }
    }
    acc_c_win <- acc_c_win + acc_c / control$theta_moves
    acc_s_win <- acc_s_win + acc_s / control$theta_moves
    acc_c_rep <- acc_c_rep + acc_c / control$theta_moves
    acc_s_rep <- acc_s_rep + acc_s / control$theta_moves

    if (it %% control$theta_thin == 0) {
      jt <- jt + 1
      tr[jt, ] <- c(it, cur$theta, cur$lp, ll, acc_c_win / control$theta_thin, acc_s_win / control$theta_thin)
      acc_c_win <- acc_s_win <- 0
    }
    if (it %in% store_at) {
      j <- j + 1
      out_x[j, ] <- x; out_theta[j, ] <- cur$theta
      bp <- beta_posterior(x, obs$X, cur$sig$S - XBX, setup$b0, B)
      out_beta[j, ] <- bp$mean + as.vector(t(chol(bp$cov)) %*% rnorm(p_beta))
    }
    first <- it == 10L && every > 10L
    if (control$progress && (it %% every == 0 || first)) {
      p <- theta_to_par(cur$theta)
      spi <- (proc.time()[["elapsed"]] - t0) / (it - it0)
      n_rep <- it - last_rep; last_rep <- it
      message(sprintf("chain %d  iter %d/%d  %.3g s/iter (%s left)  accept centred %.2f surrogate %.2f  sigma %.2f rho_s %.3f rho_t %.3f%s",
                      chain, it, n_iter, spi, format_seconds(spi * (n_iter - it)), acc_c_rep / n_rep,
                      acc_s_rep / n_rep, p$sigma, p$rho_s, p$rho_t,
                      if (setup$nugget) sprintf(" sigma_n %.2f", p$sigma_n) else ""))
      acc_c_rep <- acc_s_rep <- 0
    }
    stop_here <- !is.null(stop_after) && it == stop_after
    if (!is.null(checkpoint) && (it %% ck_every == 0 || it == n_iter || stop_here)) save_checkpoint(it)
    if (stop_here) return(NULL)
  }
  list(chain = chain, theta = as.data.frame(tr), x = out_x, draw_theta = out_theta, beta = out_beta,
       stored_iter = store_at, seconds = elapsed0 + proc.time()[["elapsed"]] - t0,
       am_centred = am_c, am_surrogate = am_s)
}

format_seconds <- function(s) {
  if (s < 10) sprintf("%.2g s", s) else if (s < 90) sprintf("%.0f s", s) else if (s < 5400) sprintf("%.0f min", s / 60) else sprintf("%.1f h", s / 3600)
}

# ---- all chains of a fit ----------------------------------------------------

# Run control$n_chains chains and return their run_chain() outputs, in chain
# order. Chains run on worker processes when control$parallel and there is
# more than one chain; otherwise one after another in this session. Either
# way the draws are the same, since each chain seeds itself: identical up to
# BLAS rounding (~1e-13; workers run BLAS on one thread, this session may
# not).
run_chains <- function(setup, control) {
  chains <- seq_len(control$n_chains)
  ck <- function(c) if (is.null(control$checkpoint_dir)) NULL else
    file.path(control$checkpoint_dir, sprintf("chain_%d.rds", c))
  if (!is.null(control$checkpoint_dir)) checkpoint_guard(control$checkpoint_dir, setup, control)
  n_workers <- min(length(chains), control$n_cores)
  if (!control$parallel || n_workers < 2) {
    return(lapply(chains, function(c) run_chain(c, setup, control, checkpoint = ck(c))))
  }
  res <- par_lapply(lapply(chains, function(c) list(chain = c, checkpoint = ck(c))),
                    function(a, setup, control)
                      utils::getFromNamespace("run_chain", "stareal")(a$chain, setup, control, a$checkpoint),
                    n_workers = n_workers, blas_threads = control$blas_threads, quiet = !control$progress,
                    setup = setup, control = control)
  bad <- which(vapply(res, inherits, logical(1), "try-error"))
  if (length(bad))
    stop(sprintf("chain %d failed: %s%s", bad[1], conditionMessage(attr(res[[bad[1]]], "condition")),
                 if (is.null(control$checkpoint_dir)) "" else
                   "\nRun fit_areal() again with the same `checkpoint_dir` to resume from the last checkpoints."),
         call. = FALSE)
  res
}

# lapply(X, FUN, ...) on fresh worker processes, which load this package
# (from source, under pkgload::load_all, as when testing). The BLAS thread
# limits are put in the environment before the workers start, since that is
# when MKL, OpenBLAS and Accelerate read them; this session's own values are
# restored afterwards. Worker messages reach this session's console unless
# `quiet`. A
# failed element comes back as a try-error rather than stopping the others.
par_lapply <- function(X, FUN, n_workers, blas_threads = 1L, quiet = FALSE, ...) {
  vars <- c("VECLIB_MAXIMUM_THREADS", "MKL_NUM_THREADS", "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS")
  old <- Sys.getenv(vars, unset = NA)
  on.exit({
    for (v in vars) if (is.na(old[[v]])) Sys.unsetenv(v) else do.call(Sys.setenv, as.list(old[v]))
  }, add = TRUE)
  do.call(Sys.setenv, as.list(stats::setNames(rep(as.character(blas_threads), length(vars)), vars)))
  cl <- parallel::makePSOCKcluster(n_workers, outfile = if (quiet) nullfile() else "")
  on.exit(parallel::stopCluster(cl), add = TRUE)
  # the functions sent to workers must not carry this namespace with them
  # before the workers have loaded it, hence environment = globalenv()
  load_pkg <- function(dev_path, lib_paths) {
    .libPaths(lib_paths)                                # libraries this session added at run time
    if (is.null(dev_path)) loadNamespace("stareal") else pkgload::load_all(dev_path, quiet = TRUE, helpers = FALSE, attach_testthat = FALSE)
    NULL
  }
  environment(load_pkg) <- globalenv()
  parallel::clusterCall(cl, load_pkg, dev_path(), .libPaths())
  wrap <- function(x, worker_fun, ...) try(worker_fun(x, ...), silent = TRUE)
  environment(wrap) <- globalenv()
  environment(FUN) <- globalenv()
  parallel::clusterApplyLB(cl, X, wrap, worker_fun = FUN, ...)
}

# the source folder when this package was loaded by pkgload::load_all(),
# else NULL (installed)
dev_path <- function() {
  if (!requireNamespace("pkgload", quietly = TRUE) || !pkgload::is_dev_package("stareal")) return(NULL)
  getNamespaceInfo(asNamespace("stareal"), "path")
}

# ---- checkpoints belong to one fit ------------------------------------------

# A checkpoint folder is tied to one fit by a fingerprint of the data, the
# model and the sampler settings, so a folder left by a different fit is
# never resumed by mistake. Settings that do not change the draws (parallel,
# workers, progress) are left out of it.
checkpoint_guard <- function(dir, setup, control) {
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
  fp <- fingerprint(list(setup$obs[c("a_o", "t_o", "T", "X")], setup$lik$y, setup$lik$k, setup$b0, setup$B,
                         setup$nugget, environment(setup$log_prior)$priors,
                         control[c("n_iter", "n_burn", "n_store", "n_chains", "ess_sweeps", "theta_moves",
                                   "theta_thin", "seed", "checkpoint_every")]))
  fp_file <- file.path(dir, "fingerprint.txt")
  old <- list.files(dir, pattern = "^chain_[0-9]+\\.rds$")
  if (length(old) && (!file.exists(fp_file) || !identical(readLines(fp_file, warn = FALSE), fp)))
    stop(sprintf(paste("`checkpoint_dir` (%s) holds checkpoints from a different fit (other data, model or",
                       "settings). Choose another folder, or delete the chain_*.rds files to start again."), dir),
         call. = FALSE)
  writeLines(fp, fp_file)
  if (length(old)) message(sprintf("checkpoints found in %s: resuming", dir))
  invisible(fp)
}

fingerprint <- function(x) {
  tmp <- tempfile(); on.exit(unlink(tmp))
  saveRDS(x, tmp, compress = FALSE)
  unname(tools::md5sum(tmp))
}
