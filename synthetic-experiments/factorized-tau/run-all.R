## Run every factorized-tau experiment (factorized vs 3-D grid vs known tau) and collect
## their summaries.
##
## Each experiment runs in its own Rscript process (they share global names and set their
## own seed), with its output in experiments/synthetic/factorized-tau/<experiment>/ and its
## console log in experiments/synthetic/factorized-tau/logs/<experiment>.log. The combined
## summaries and the status of every run are written to
## experiments/synthetic/factorized-tau/summary-all.txt.
##
## Run from the repository root:
##   Rscript synthetic-experiments/factorized-tau/run-all.R              # all, one at a time
##   Rscript synthetic-experiments/factorized-tau/run-all.R --jobs 4     # all, 4 in parallel
##   Rscript synthetic-experiments/factorized-tau/run-all.R multisite-singledisturb-multipool

experiments <- c("singlesite-singledisturb-multipool",
                 "singlesite-multidisturb-multipool",
                 "singlesite-multidisturb-multipool-multipft",
                 "multisite-singledisturb-multipool")
script.dir <- "synthetic-experiments/factorized-tau"
out.dir    <- "experiments/synthetic/factorized-tau"
log.dir    <- file.path(out.dir, "logs")

## ---- arguments ---------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
jobs <- 1L
if (length(i <- which(args == "--jobs"))) {
  jobs <- as.integer(args[i + 1]); args <- args[-c(i, i + 1)]
  if (is.na(jobs) || jobs < 1) stop("--jobs needs a positive integer")
}
if (length(args)) {
  bad <- setdiff(args, experiments)
  if (length(bad)) stop("unknown experiment(s): ", paste(bad, collapse = ", "),
                        "\navailable: ", paste(experiments, collapse = ", "))
  experiments <- args
}
if (!file.exists(file.path(script.dir, "or-joint-update-factorized.R")))
  stop("run from the repository root")
dir.create(log.dir, recursive = TRUE, showWarnings = FALSE)

## ---- run ---------------------------------------------------------------------------------
run_one <- function(name) {
  log <- file.path(log.dir, paste0(name, ".log"))
  t0 <- Sys.time()
  status <- system2("Rscript", file.path(script.dir, paste0(name, ".R")), stdout = log, stderr = log)
  mins <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  message(sprintf("%-45s %s in %.1f min", name, if (status == 0) "done" else sprintf("FAILED (status %d)", status), mins))
  list(name = name, status = status, minutes = mins, log = log)
}
message(sprintf("running %d experiment(s), %d at a time; logs in %s", length(experiments), jobs, log.dir))
res <- parallel::mclapply(experiments, run_one, mc.cores = jobs, mc.preschedule = FALSE)

## ---- collect -----------------------------------------------------------------------------
status <- data.frame(experiment = experiments,
                     status = sapply(res, function(r) if (r$status == 0) "ok" else "FAILED"),
                     minutes = round(sapply(res, `[[`, "minutes"), 1))
out <- c(sprintf("Factorized-tau experiments, run %s", format(Sys.time(), "%Y-%m-%d %H:%M")), "",
         capture.output(print(status, row.names = FALSE)))
for (r in res) {
  sm <- file.path(out.dir, r$name, "summary.txt")
  out <- c(out, "", strrep("=", 80), r$name, strrep("=", 80))
  out <- c(out, if (r$status == 0 && file.exists(sm)) readLines(sm)
                else c(sprintf("failed; last lines of %s:", r$log), tail(readLines(r$log), 15)))
}
writeLines(out, file.path(out.dir, "summary-all.txt"))
message("summaries: ", file.path(out.dir, "summary-all.txt"))
if (any(status$status != "ok")) quit(status = 1)
