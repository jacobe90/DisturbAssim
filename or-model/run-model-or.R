## Run the OR joint-update filter on the loaded data.
##
## Uses the inputs built by load-data-or.R: those already in this R session (when run by
## load-data-and-run-model.R), otherwise the cache data/cache/or-model-inputs.RData.
## Runs ORJointUpdate for every requested (site, observation-error model, disturbance
## product, tau prior) on forked workers, printing progress after each run, and saves each
## run to <outdir>/site_<idx>_distprod_<TF>_<source>_<model>_tau<prior>[_bias].rds.
## Runs already on disk are skipped, so an interrupted run can be restarted with the same
## command. Then runs or-model/analysis-or.R on <outdir>, which writes its summary reports
## there (exps_summary.pdf, distprod_TRUE_vs_FALSE.pdf, elist.RData, stats.RData), and with
## at most --plot-max sites also writes <outdir>/site-plots.pdf.
##
## Run from the repository root:
##   Rscript or-model/run-model-or.R                         # full grid (all sites, 8 settings)
##   Rscript or-model/run-model-or.R --sites 300 --R GEDI_het --distprod TRUE --tau-prior loguniform
##   Rscript or-model/run-model-or.R --sites 1:50 --jobs 4 --bias off
##
## Options (defaults in brackets):
##   --sites       sites, an R expression or "all"                          [all]
##   --R           observation-error models, comma-separated, from
##                 GEDI_het, GEDI_con, LandTrendr_het, LandTrendr_con        [all four]
##   --distprod    TRUE, FALSE or TRUE,FALSE (use the disturbance product)   [TRUE,FALSE]
##   --tau-prior   gamma, loguniform or gamma,loguniform                     [gamma]
##   --bias        on / off: online bias correction of the undisturbed wood
##                 forecast (rho = 0.9, denom = "N", skip = 1)               [on]
##   --dist-proc-err  on / off: propagate the process error through the
##                 disturbance (see ORJointUpdate)                           [on]
##   --outdir      output folder         [experiments/or-model/runs, with "-<site set>" for
##                 a site set other than the tracked one and "-noDistQ" with --dist-proc-err off]
##   --jobs        parallel workers                                          [cores - 2]
##   --timeout     seconds before a stuck run is killed and retried once     [300]
##   --overwrite   rerun runs that are already on disk
##   --plot-max    largest number of sites for which site plots are drawn    [10]
##   --no-analysis skip the analysis reports

source("or-model/R/config.R")
check_repo_root()
suppressPackageStartupMessages(library(BayesianTools))   # VSEM
for (f in c("R/disturbance.R", "or-model/R/or-joint-update.R", "or-model/R/make-or-args.R",
            "or-model/R/run-grid.R", "or-model/R/plots.R")) source(f)

## ---- options ---------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
opt  <- function(name, default) { i <- which(args == name); if (length(i)) args[i + 1] else default }
split_opt <- function(x) strsplit(x, ",", fixed = TRUE)[[1]]
on_off <- function(x, name) { if (!x %in% c("on", "off")) stop(name, " must be on or off"); x == "on" }
cfg <- list(
  sites       = opt("--sites", "all"),
  R           = split_opt(opt("--R", paste(names(Rparams.grid), collapse = ","))),
  distprod    = as.logical(split_opt(opt("--distprod", "TRUE,FALSE"))),
  tauPrior    = split_opt(opt("--tau-prior", "gamma")),
  bias        = on_off(opt("--bias", "on"), "--bias"),
  distProcErr = on_off(opt("--dist-proc-err", "on"), "--dist-proc-err"),
  outdir      = opt("--outdir", NA),
  jobs        = as.integer(opt("--jobs", max(1L, parallel::detectCores() - 2L))),
  timeout     = as.numeric(opt("--timeout", 300)),
  overwrite   = "--overwrite" %in% args,
  plotMax     = as.integer(opt("--plot-max", 10)),
  analysis    = !"--no-analysis" %in% args)
if (length(bad <- setdiff(cfg$R, names(Rparams.grid))))
  stop("unknown --R: ", paste(bad, collapse = ", "), "; choose from ", paste(names(Rparams.grid), collapse = ", "))
if (anyNA(cfg$distprod)) stop("--distprod must be TRUE, FALSE or TRUE,FALSE")
if (length(bad <- setdiff(cfg$tauPrior, c("gamma", "loguniform")))) stop("unknown --tau-prior: ", bad)
bias.cfg <- if (cfg$bias) list(rho = 0.9, denom = "N", skip = 1) else NULL

## ---- inputs ----------------------------------------------------------------------------
if (exists("or.inputs.loaded") && isTRUE(or.inputs.loaded)) {
  say("OR model: using the inputs loaded in this session")
} else {
  if (!file.exists(or.files$orInputs))
    stop("no ", or.files$orInputs, ": run  Rscript or-model/load-data-or.R  first", call. = FALSE)
  say("OR model: reading %s (built %s)", or.files$orInputs,
      format(file.mtime(or.files$orInputs), "%Y-%m-%d %H:%M"))
  load(or.files$orInputs)
}
if (!exists("or.extract.file")) or.extract.file <- or.files$extract   # caches from before --extract
say("OR model: sites from %s (%d sites)", or.extract.file, nrow(dextract))
if (is.na(cfg$outdir)) {
  set.tag <- if (or.extract.file == or.files$extract) "" else
    paste0("-", sub("^extract\\.h03v04-?", "", tools::file_path_sans_ext(basename(or.extract.file))))
  cfg$outdir <- paste0("experiments/or-model/runs", set.tag, if (cfg$distProcErr) "" else "-noDistQ")
}
sites <- if (cfg$sites == "all") seq_len(nrow(dextract)) else eval(parse(text = cfg$sites))
if (any(sites < 1 | sites > nrow(dextract))) stop("--sites must be in 1:", nrow(dextract))

## ---- run -------------------------------------------------------------------------------
exp.grid <- make_exp_grid(R = cfg$R, useDistProduct = cfg$distprod, tauPrior = cfg$tauPrior, sites = sites,
                          outdir = cfg$outdir, bias = bias.cfg)
say("OR model: %d site(s) x %d setting(s) = %d runs", length(sites), nrow(exp.grid) / length(sites), nrow(exp.grid))
message(sprintf("  observation error: %s\n  disturbance product: %s\n  tau prior: %s\n  bias correction: %s\n  process error through disturbance: %s\n  output: %s",
                paste(cfg$R, collapse = ", "), paste(cfg$distprod, collapse = ", "),
                paste(cfg$tauPrior, collapse = ", "), if (cfg$bias) "on" else "off",
                if (cfg$distProcErr) "on" else "off", cfg$outdir))
## runs already in outdir are skipped (and all of them go into the analysis): make sure
## they were made with the same setting and site set
old <- dir(cfg$outdir, pattern = "^site_.*\\.rds$", full.names = TRUE)
if (length(old) && !cfg$overwrite) {
  prev <- readRDS(old[1])
  prevQ <- isTRUE(prev$distProcErr)                             # older runs: off
  prevX <- if (is.null(prev$extractFile)) or.files$extract else prev$extractFile
  if (prevQ != cfg$distProcErr)
    stop(sprintf("%s holds runs made with --dist-proc-err %s; use another --outdir or --overwrite",
                 cfg$outdir, if (prevQ) "on" else "off"), call. = FALSE)
  if (prevX != or.extract.file)
    stop(sprintf("%s holds runs on the sites of %s, not %s; use another --outdir",
                 cfg$outdir, prevX, or.extract.file), call. = FALSE)
}
t.run <- Sys.time()
grid.out <- run_or_grid(exp.grid, cores = cfg$jobs, timeout = cfg$timeout, overwrite = cfg$overwrite,
                        bias = bias.cfg, distProcErr = cfg$distProcErr)
failed <- which(sapply(grid.out, inherits, "try-error"))
say("OR model: %d of %d runs ok, %d failed, in %s", nrow(exp.grid) - length(failed), nrow(exp.grid),
    length(failed), fmt_secs(as.numeric(difftime(Sys.time(), t.run, units = "secs"))))
if (length(failed)) print(cbind(exp.grid[failed, c("R", "useDistProduct", "tauPrior", "site")],
                                error = sapply(grid.out[failed], as.character)))

## ---- site plots ------------------------------------------------------------------------
ok <- setdiff(seq_len(nrow(exp.grid)), failed)
if (length(ok) && length(sites) <= cfg$plotMax) {
  f <- file.path(cfg$outdir, "site-plots.pdf")
  pdf(f, width = 9, height = 4.5)
  for (i in ok) {
    h <- load_hist(exp.grid$file[i])
    par(mfrow = c(1, 2))
    plot_analysis(h, main = sprintf("Site %d: wood-pool analysis", h$site))
    plot_pdist(h, main = "probability of disturbance")
    mtext(sub("^site_\\d+_", "", h$tag), side = 1, line = -1, outer = TRUE, cex = 0.7)
  }
  invisible(dev.off())
  say("site plots: %s", f)
}
## ---- analysis reports --------------------------------------------------------------------
if (cfg$analysis && length(ok)) {
  say("OR model: analysis of the runs in %s", cfg$outdir)
  ae <- new.env()
  ae$analysis.indir <- cfg$outdir
  tryCatch(sys.source("or-model/analysis-or.R", envir = ae),
           error = function(e) say("analysis FAILED: %s", conditionMessage(e)))
}
if (length(failed)) quit(status = 1)
