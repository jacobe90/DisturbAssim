## Run ORJointUpdate for sites and settings, saving one .rds per run.
##
## Ported from the run chunks of OR_joint_update2.Rmd (now legacy/models/). Each run is
## saved as <outdir>/site_<idx>_distprod_<TF>_<source>_<model>_tau<prior>[_bias].rds, which
## or-model/analysis-or.R parses. Whether the process error was propagated through the
## disturbance is stored in the run ($distProcErr), not in the name: keep runs with it on and
## off in different folders.

## Generate the ORJointUpdate inputs for each requested site, run the joint update and
## save each run. Returns the runs invisibly (keyed by site).
## Rparams.source ("GEDI"/"LandTrendr") and Rparams.model ("het"/"con") describe the
## observation-error parameters in use; together with tauPrior ("gamma"/"loguniform") they
## are recorded in each output filename. distProcErr: see ORJointUpdate.
run_joint_update_sites <- function(site.idx, Rparams = Rparams.GEDI.het,
                                   useDistProduct = TRUE,
                                   tauPrior = c("gamma", "loguniform"),
                                   Rparams.source = c("GEDI", "LandTrendr"),
                                   Rparams.model  = c("het", "con"),
                                   outdir = "experiments/or-model",
                                   bias = NULL, distProcErr = TRUE) {
  Rparams.source <- match.arg(Rparams.source)
  Rparams.model  <- match.arg(Rparams.model)
  tauPrior       <- match.arg(tauPrior)
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
  results <- list()
  for (site in site.idx) {
    args <- make_or_args(site, useDistProduct, tauPrior)
    out <- ORJointUpdate(Y = args$Y, Rparams = Rparams,
                         Q = Q, PAR = args$PAR, Xic = ICprior,
                         THETA = theta[1:11], dist = args$disturb.params, prior = args$prior,
                         grid = args$grid, bias = bias, distProcErr = distProcErr)
    out$year <- atime                  # calendar year of each DA timestep
    out$obs  <- args$Y[, 3]            # assimilated obs, aligned with out$year (NA = none)
    out$biasConfig <- bias
    out$distProcErr <- distProcErr
    out$extractFile <- get0("or.extract.file", ifnotfound = or.files$extract)   # the site set
    fname <- run_file(outdir, site, useDistProduct, Rparams.source, Rparams.model, tauPrior, !is.null(bias))
    saveRDS(out, fname)
    results[[as.character(site)]] <- out
  }
  invisible(results)
}

run_file <- function(outdir, site, useDistProduct, source, model, tauPrior, bias)
  file.path(outdir, sprintf("site_%d_distprod_%s_%s_%s_tau%s%s.rds", site, useDistProduct, source,
                            model, tauPrior, ifelse(bias, "_bias", "")))

## The experiment grid: every combination of observation-error model (names of
## Rparams.grid), disturbance product (with / without), tau prior and site, with the file
## each run is saved to. Row i of the grid over ALL sites is run with seed i, so a run's
## result does not depend on which subset of sites is run, in what order or on how many
## workers (and the default settings give the same seeds as the old full grid).
make_exp_grid <- function(R = names(Rparams.grid), useDistProduct = c(TRUE, FALSE), tauPrior = "gamma",
                          sites = seq_len(nrow(dextract)), outdir, bias = NULL) {
  g <- expand.grid(R = R, useDistProduct = useDistProduct, tauPrior = tauPrior,
                   site = seq_len(nrow(dextract)), stringsAsFactors = FALSE)
  g$seed <- seq_len(nrow(g))
  g <- g[g$site %in% sites, , drop = FALSE]
  g$file <- run_file(outdir, g$site, g$useDistProduct, sapply(Rparams.grid[g$R], `[[`, "source"),
                     sapply(Rparams.grid[g$R], `[[`, "model"), g$tauPrior, !is.null(bias))
  rownames(g) <- NULL
  g
}

## Run every row of exp.grid on `cores` forked workers, collecting runs as they finish and
## printing progress (% finished, failures, elapsed time, ETA) after each one. Existing
## files are skipped unless overwrite = TRUE, so an interrupted grid can be restarted. A
## worker still running after `timeout` seconds is killed and its run retried once.
## Forked workers can deadlock inside RStudio; run this with Rscript.
## Returns, per row, the file written or a "try-error" with the failure.
run_or_grid <- function(exp.grid, cores = 1, timeout = 300, overwrite = FALSE, bias = NULL,
                        distProcErr = TRUE) {
  dir.create(unique(dirname(exp.grid$file)), recursive = TRUE, showWarnings = FALSE)
  run_row <- function(i) {
    e <- exp.grid[i, ]
    r <- Rparams.grid[[e$R]]
    set.seed(e$seed)
    ## catch errors per run so one failure doesn't take down the worker. Messages and
    ## warnings are muffled inside the worker: in RStudio a forked child can otherwise block
    ## forever on a lock held by another RStudio thread at fork time.
    tryCatch(suppressWarnings(suppressMessages({
      run_joint_update_sites(e$site, Rparams = r$Rparams, useDistProduct = e$useDistProduct,
                             tauPrior = e$tauPrior, Rparams.source = r$source, Rparams.model = r$model,
                             outdir = dirname(e$file), bias = bias, distProcErr = distProcErr)
      e$file
    })), error = function(err) structure(conditionMessage(err), class = "try-error"))
  }
  label <- function(i) sprintf("site %d, %s, distprod %s, tau %s", exp.grid$site[i], exp.grid$R[i],
                               exp.grid$useDistProduct[i], exp.grid$tauPrior[i])

  n.total <- nrow(exp.grid)
  out     <- vector("list", n.total)
  todo    <- if (overwrite) seq_len(n.total) else which(!file.exists(exp.grid$file))
  done0   <- setdiff(seq_len(n.total), todo)
  out[done0] <- exp.grid$file[done0]
  n.skip  <- length(done0)
  n.done  <- n.skip; n.fail <- 0
  running <- list()   # pid -> list(job, row, start)
  retried <- integer(0)
  t0      <- Sys.time()
  report  <- function() {
    elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    n.new   <- n.done - n.skip
    eta     <- if (n.new > 0) elapsed / n.new * (n.total - n.done) else NA
    say("%d / %d runs finished (%.1f%%) | %d failed | %d running | elapsed %s | ETA %s",
        n.done, n.total, 100 * n.done / n.total, n.fail, length(running), fmt_secs(elapsed), fmt_secs(eta))
  }
  say("%d runs, %d already done (skipped), %d to run on %d worker(s)", n.total, n.skip, length(todo), cores)
  report()
  while (length(todo) || length(running)) {
    while (length(running) < cores && length(todo)) {
      i   <- todo[1]; todo <- todo[-1]
      job <- parallel::mcparallel(run_row(i), silent = TRUE)
      running[[as.character(job$pid)]] <- list(job = job, row = i, start = Sys.time())
    }
    for (pid in names(running)) {
      if (as.numeric(difftime(Sys.time(), running[[pid]]$start, units = "secs")) < timeout) next
      i <- running[[pid]]$row
      tools::pskill(as.integer(pid), tools::SIGKILL)
      parallel::mccollect(running[[pid]]$job, wait = TRUE)   # reap the killed child
      running[[pid]] <- NULL
      if (!(i %in% retried)) {
        retried <- c(retried, i); todo <- c(i, todo)
        say("  TIMEOUT: %s -- killed after %gs, retrying", label(i), timeout)
      } else {
        out[[i]] <- structure(sprintf("timed out twice (> %gs)", timeout), class = "try-error")
        n.done <- n.done + 1; n.fail <- n.fail + 1
        say("  FAILED: %s: timed out twice", label(i))
      }
    }
    if (!length(running)) next
    res <- parallel::mccollect(lapply(running, `[[`, "job"), wait = FALSE, timeout = 5)
    for (pid in names(res)) {
      i <- running[[pid]]$row
      o <- res[[pid]]
      if (is.null(o)) o <- structure("worker died (no result returned)", class = "try-error")
      out[[i]] <- o
      n.done <- n.done + 1
      if (inherits(o, "try-error")) { n.fail <- n.fail + 1; say("  FAILED: %s: %s", label(i), o) }
      running[[pid]] <- NULL
    }
    if (length(res)) report()
  }
  out
}
