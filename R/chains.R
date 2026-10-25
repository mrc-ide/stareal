# chains.R
#
# Running MCMC chains: one chain (run_chain, with checkpoint and resume) and
# all chains of a fit (run_chains: one chain per worker process, in
# parallel, inside the R session; workers.R). The moves themselves are in
# sampler.R. Internal.
#
# run_chain() is ported from Africa_areal/R/chains.R (2026-10-12). Changes:
# - `setup` (fit_setup(), below) carries the likelihood as one object, so
#   lik$loglik, lik$pseudo and lik$pseudo_var stand where ll_fun, g_star and
#   S_aux stood; `profile` is now `control` (mcmc_control());
# - any number of fixed effects;
# - the pseudo-data that steer the latent ESS and the surrogate theta move
#   are refitted during burn-in (at 1/2 and 3/4 of it) to a Gaussian
#   approximation of the likelihood at the chain's mean latent values, then
#   frozen (likelihood.R says why). Africa_areal kept the empirical logits
#   throughout, which for cells with y = 0 or few positives put the ellipse
#   in the wrong place;
# - chains start from dispersed points: hyperparameters spread about the
#   same centre (sd 0.5 on log scales, 1 on logit scales) and latent values
#   drawn about the empirical logits. Africa_areal started every chain at
#   nearly the same point, which weakens R-hat;
# - each chain has its own seed, drawn from the user's (any seed works);
# - stored draws are spread over the whole post-burn-in run;
# - no heartbeat file; the state is also checkpointed at the last
#   iteration, so a finished chain is read back, not rerun; a checkpoint
#   that cannot be written stops the chain with a message;
# - progress messages say how long is left (this chain plus any queued
#   behind it), and the first comes after 10 iterations.
#
# ------------------------------------------------------------------

# bump when a change to the sampler means old checkpoints must not be
# continued
SAMPLER_VERSION <- 2L

# ---- everything the sampler needs about one data set -------------------------

# cells: areal_data()$cells (area, t, y, k). eig: the frame's
# eigendecomposition. X: the fixed-effect design of the observed cells.
# priors: areal_priors(), with beta_mean / beta_sd already of length p.
fit_setup <- function(cells, eig, T, priors, nugget, X) {
  obs <- make_obs(cells$area, cells$t, T, X = X)
  p <- ncol(X)
  b0 <- priors$beta_mean; B <- diag(priors$beta_sd^2, p)
  list(obs = obs, eig = eig, lik = lik_binomial(cells$y, cells$k), b0 = b0, B = B,
       mu_o = as.vector(obs$X %*% b0), nugget = isTRUE(nugget), log_prior = make_log_prior(priors),
       beta_names = colnames(X))
}

# ---- one chain --------------------------------------------------------------

# setup: fit_setup() output. control: mcmc_control() output. seed: this
# chain's seed. queued: whole chains waiting behind this one on the same
# worker (for the time-left estimate). Returns
#   theta       data frame: iteration, theta (unconstrained), lp, ll, and the
#               acceptance of the centred and surrogate moves (mean over the
#               thinning window), every theta_thin-th iteration, burn-in
#               included. lp is the log prior of theta plus the log density
#               of the latent values given theta (the centred move's
#               target); lp + ll is the unnormalised log joint
#   x           n_store x n latent draws (spread over the post-burn-in run)
#   draw_theta  n_store x d theta at those same iterations
#   beta        n_store x p fixed-effect draws, one from beta | x_o, theta
#               per stored latent draw
#   stored_iter, seconds
#
# checkpoint: a file path. Every control$checkpoint_every iterations the
#   whole state (latent values, theta, pseudo-data, adaptation, outputs so
#   far and the random-number state) is saved there; if the file exists at
#   the start, the chain continues from it and produces exactly the chain an
#   uninterrupted run would have (tests/testthat/test-fit.R).
# stop_after: stop at this iteration, after saving a checkpoint, and return
#   NULL (for the tests).
run_chain <- function(chain, setup, control, seed, checkpoint = NULL, stop_after = NULL, queued = 0) {
  obs <- setup$obs; eig <- setup$eig; B <- setup$B; mu_o <- setup$mu_o
  lik <- setup$lik; prior <- setup$log_prior; ll_fun <- lik$loglik
  n_iter <- control$n_iter; n_burn <- control$n_burn; p_beta <- ncol(obs$X)
  store_at <- store_iterations(n_iter, n_burn, control$n_store)
  every <- if (n_iter >= 5000) 500L else max(10L, n_iter %/% 10L)
  ck_every <- as.integer(control$checkpoint_every)
  XBX <- obs$X %*% B %*% t(obs$X)
  # refits of the pseudo-data: at the end of iterations refit_at, each from
  # the mean of x over the iterations since the previous refit (or since a
  # quarter of burn-in). None for a burn-in shorter than 20 iterations.
  refit_at <- if (n_burn >= 20) unique(c(n_burn %/% 2, (3 * n_burn) %/% 4)) else integer(0)
  window_start <- n_burn %/% 4

  if (!is.null(checkpoint) && file.exists(checkpoint)) {
    ck <- readRDS(checkpoint)                             # resume
    list2env(ck$state, environment())
    assign(".Random.seed", ck$seed, envir = globalenv())
    it0 <- ck$it; elapsed0 <- ck$elapsed
    if (control$verbose && it0 < n_iter)
      message(sprintf("chain %d: resuming from the checkpoint at iteration %d", chain, it0))
  } else {
    set_seed(seed)
    # dispersed starting values: hyperparameters about a common centre,
    # latents about the empirical logits
    d <- if (setup$nugget) 4L else 3L
    theta <- (c(log(1), stats::qlogis(0.9), stats::qlogis(0.95), log(0.5)) +
                rnorm(4, 0, c(0.5, 1, 1, 0.5)))[seq_len(d)]
    g_star <- lik$pseudo; S_aux <- lik$pseudo_var
    x <- g_star + rnorm(obs$n, 0, sqrt(pmin(S_aux, 4)))
    ll <- ll_fun(x)
    # starting proposal sds on the unconstrained scale (roughly the posterior
    # sds of the first Africa_areal pilots; the adaptation corrects them). The
    # two moves adapt separately: their acceptance landscapes differ.
    sd0 <- c(0.03, 0.4, 0.12, 0.1)[seq_len(d)]
    am_c <- am_init(theta, sd0); am_s <- am_init(theta, sd0)
    trace_at <- seq(control$theta_thin, n_iter, by = control$theta_thin)
    tr <- matrix(NA_real_, length(trace_at), d + 5,
                 dimnames = list(NULL, c("iteration", theta_names(setup$nugget), "lp", "ll",
                                         "accept_centred", "accept_surrogate")))
    out_x <- matrix(NA_real_, control$n_store, obs$n)
    out_beta <- matrix(NA_real_, control$n_store, p_beta, dimnames = list(NULL, setup$beta_names))
    out_theta <- matrix(NA_real_, control$n_store, d, dimnames = list(NULL, theta_names(setup$nugget)))
    acc_c_win <- acc_s_win <- 0; j <- 0; jt <- 0
    acc_c_rep <- acc_s_rep <- 0; n_rep <- 0               # acceptance since the last progress line
    x_sum <- numeric(obs$n); x_cnt <- 0                   # for the next refit of the pseudo-data
    it0 <- 0; elapsed0 <- 0
  }
  # derived from the state, so not checkpointed (and recomputed identically)
  cur <- list(theta = theta, sig = build_sigma_oo(theta, obs, eig, B))
  cur$lp <- prior(theta) + log_dmvnorm_chol(x, mu_o, cur$sig$L, cur$sig$logdet)
  prop_c <- am_prop_chol(am_c); prop_s <- am_prop_chol(am_s)
  t0 <- proc.time()[["elapsed"]]
  save_checkpoint <- function(it) {
    state <- list(theta = cur$theta, x = x, ll = ll, am_c = am_c, am_s = am_s, tr = tr, out_x = out_x,
                  out_beta = out_beta, out_theta = out_theta, acc_c_win = acc_c_win,
                  acc_s_win = acc_s_win, j = j, jt = jt, acc_c_rep = acc_c_rep, acc_s_rep = acc_s_rep,
                  n_rep = n_rep, g_star = g_star, S_aux = S_aux, x_sum = x_sum, x_cnt = x_cnt)
    # write, then rename: never a half-written file (the temporary name is
    # per process, so two processes never share one)
    tmp <- sprintf("%s.%d.tmp", checkpoint, Sys.getpid())
    saveRDS(list(it = it, elapsed = elapsed0 + proc.time()[["elapsed"]] - t0, state = state,
                 seed = get(".Random.seed", envir = globalenv())), tmp, compress = FALSE)
    for (attempt in 1:3) {
      if (isTRUE(suppressWarnings(file.rename(tmp, checkpoint)))) return(invisible())
      Sys.sleep(1)                                        # a sync client or scanner may hold the file
    }
    unlink(tmp)
    stop(sprintf("could not write the checkpoint %s (is the folder synced, e.g. by Dropbox?)", checkpoint),
         call. = FALSE)
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
    # 3. burn-in only: refit the pseudo-data at the chain's mean latent values
    if (length(refit_at) && it > window_start && it <= max(refit_at)) {
      x_sum <- x_sum + x; x_cnt <- x_cnt + 1
      if (it %in% refit_at) {
        ps <- lik$pseudo_at(x_sum / x_cnt)
        g_star <- ps$pseudo; S_aux <- ps$pseudo_var
        cur$sig$sur <- NULL; cur$sig$m_star <- NULL; cur$sig$L_post <- NULL   # rebuilt from the new pseudo-data
        x_sum <- numeric(obs$n); x_cnt <- 0
      }
    }
    acc_c_win <- acc_c_win + acc_c / control$theta_moves
    acc_s_win <- acc_s_win + acc_s / control$theta_moves
    acc_c_rep <- acc_c_rep + acc_c / control$theta_moves
    acc_s_rep <- acc_s_rep + acc_s / control$theta_moves
    n_rep <- n_rep + 1

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
    if (it %% every == 0 || (it == 10L && every > 10L)) {
      if (control$verbose) {
        p <- theta_to_par(cur$theta)
        spi <- (proc.time()[["elapsed"]] - t0) / (it - it0)
        message(sprintf("chain %d  iter %d/%d  %.3g s/iter (%s left)  accept centred %.2f surrogate %.2f  sigma %.2f rho_s %.3f rho_t %.3f%s",
                        chain, it, n_iter, spi, format_seconds(spi * (n_iter - it + queued * n_iter)),
                        acc_c_rep / n_rep, acc_s_rep / n_rep, p$sigma, p$rho_s, p$rho_t,
                        if (setup$nugget) sprintf(" sigma_n %.2f", p$sigma_n) else ""))
      }
      acc_c_rep <- acc_s_rep <- 0; n_rep <- 0
    }
    stop_here <- !is.null(stop_after) && it == stop_after
    if (!is.null(checkpoint) && (it %% ck_every == 0 || it == n_iter || stop_here)) save_checkpoint(it)
    if (stop_here) return(NULL)
  }
  list(chain = chain, theta = as.data.frame(tr), x = out_x, draw_theta = out_theta, beta = out_beta,
       stored_iter = store_at, seconds = elapsed0 + proc.time()[["elapsed"]] - t0)
}

format_seconds <- function(s) {
  if (s < 10) sprintf("%.2g s", s) else if (s < 90) sprintf("%.0f s", s) else if (s < 5400)
    sprintf("%.0f min", s / 60) else sprintf("%.1f h", s / 3600)
}

# the chains' seeds, drawn from the user's seed
chain_seeds <- function(seed, n_chains) with_seed(seed, sample.int(.Machine$integer.max - 1L, n_chains))

# ---- all chains of a fit ----------------------------------------------------

# Run control$n_chains chains and return their run_chain() outputs, in chain
# order. Chains run on control$n_cores worker processes (one chain each at a
# time), or one after another in this session when n_cores is 1. Either
# way the draws are the same, since each chain seeds itself: identical up to
# BLAS rounding (~1e-13; workers run BLAS on blas_threads threads, this
# session may not). W: the frame's adjacency, for the checkpoint guard.
run_chains <- function(setup, control, checkpoint_dir = NULL, W = NULL) {
  chains <- seq_len(control$n_chains)
  seeds <- chain_seeds(control$seed, control$n_chains)
  ck <- function(c) if (is.null(checkpoint_dir)) NULL else file.path(checkpoint_dir, sprintf("chain_%d.rds", c))
  if (!is.null(checkpoint_dir)) checkpoint_guard(checkpoint_dir, setup, control, W)
  n_workers <- min(control$n_chains, control$n_cores)
  # whole chains queued behind chain c on its worker (chains start in waves)
  waves <- ceiling(control$n_chains / n_workers)
  queued <- waves - ceiling(chains / n_workers)
  jobs <- lapply(chains, function(c) list(chain = c, seed = seeds[c], checkpoint = ck(c), queued = queued[c]))
  run_one <- function(a, shared)
    utils::getFromNamespace("run_chain", "stareal")(a$chain, shared$setup, shared$control, a$seed,
                                                    checkpoint = a$checkpoint, queued = a$queued)
  tryCatch(run_jobs(jobs, run_one, shared = list(setup = setup, control = control), n_workers = n_workers,
                    blas_threads = control$blas_threads, verbose = control$verbose),
           stareal_job_error = function(e)
             stop(sprintf("chain %d failed: %s%s", e$job, conditionMessage(e),
                          if (is.null(checkpoint_dir)) "" else
                            "\nRun fit_areal() again with the same `checkpoint_dir` to resume from the last checkpoints."),
                  call. = FALSE))
}

# ---- checkpoints belong to one fit ------------------------------------------

# A checkpoint folder is tied to one fit by a fingerprint of the data, the
# graph, the model, the sampler settings that change the draws, and the
# sampler version, so a folder left by a different fit is never resumed by
# mistake. Settings that do not change the draws (workers, BLAS threads,
# progress, checkpoint_every) are left out, and so is the R version.
checkpoint_guard <- function(dir, setup, control, W = NULL) {
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
  fp <- fingerprint(list(SAMPLER_VERSION, setup$obs[c("a_o", "t_o", "T", "X")], setup$lik$y, setup$lik$k,
                         setup$b0, setup$B, setup$nugget, environment(setup$log_prior)$priors,
                         if (is.null(W)) NULL else list(dim(W), W@i, W@p),
                         control[c("n_iter", "n_burn", "n_store", "n_chains", "ess_sweeps", "theta_moves",
                                   "theta_thin", "seed")]))
  fp_file <- file.path(dir, "fingerprint.txt")
  old <- list.files(dir, pattern = "^chain_[0-9]+\\.rds$")
  if (length(old) && (!file.exists(fp_file) || !identical(readLines(fp_file, warn = FALSE), fp)))
    stop(sprintf(paste("`checkpoint_dir` (%s) holds checkpoints from a different fit (other data, frame, model,",
                       "settings or sampler version). Choose another folder, or delete the chain_*.rds files",
                       "to start again."), dir),
         call. = FALSE)
  writeLines(fp, fp_file)
  if (length(old) && control$verbose) message(sprintf("checkpoints found in %s: resuming", dir))
  invisible(fp)
}

# an md5 fingerprint of an R object, independent of the R version (the
# serialisation header, which records it, is dropped)
fingerprint <- function(x) {
  b <- serialize(x, NULL, version = 2)
  tmp <- tempfile(); on.exit(unlink(tmp))
  writeBin(b[-seq_len(14)], tmp)
  unname(tools::md5sum(tmp))
}
