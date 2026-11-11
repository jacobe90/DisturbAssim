## Dropped forecast mass of the factorized-tau filter across seeds.
##
## For every factorized-tau experiment and every seed, rebuilds the experiment's set-up
## (drivers, ensemble, truth, observations; all drawn from the seed) and its filter
## arguments, then runs only the factorized-tau filter (no 3-D grid or known-tau
## comparison). Records, per step, the forecast probability mass dropped for collapsed
## importance weights and for singular covariances, and the number of grid points left
## with zero mass after dropping (filter-diagnostics.R).
##
## Output, in experiments/synthetic/dropped-across-seeds/:
##   <experiment>.pdf   dropped mass and zero-mass grid points, seeds x years
##   summary.txt        per experiment: totals and the seeds / years with any drops
##   results.rds        the per-run records (runs/: one record and log per run)
##
## Run from the repository root:
##   Rscript synthetic-experiments/factorized-tau/run-seeds-dropped.R               # all
##   Rscript synthetic-experiments/factorized-tau/run-seeds-dropped.R --jobs 8 --seeds 1:5
##   Rscript synthetic-experiments/factorized-tau/run-seeds-dropped.R --plot-only   # redraw

source("synthetic-experiments/factorized-tau/or-joint-update-factorized.R")

experiments <- c("singlesite-singledisturb-multipool",
                 "singlesite-multidisturb-multipool",
                 "singlesite-multidisturb-multipool-multipft",
                 "multisite-singledisturb-multipool")
seeds   <- c(42, 1:20)
jobs    <- max(1L, parallel::detectCores() - 2L)
out.dir <- "experiments/synthetic/dropped-across-seeds"

args <- commandArgs(trailingOnly = TRUE)
opt  <- function(name) { i <- which(args == name); if (length(i)) args[i + 1] else NULL }
if (!is.null(opt("--jobs")))  jobs  <- as.integer(opt("--jobs"))
if (!is.null(opt("--seeds"))) seeds <- eval(parse(text = opt("--seeds")))
plot.only <- "--plot-only" %in% args
dir.create(out.dir, recursive = TRUE, showWarnings = FALSE)

## ---- run ---------------------------------------------------------------------------------
## Each (experiment, seed) runs in its own Rscript process (the set-ups use global
## variables): the experiment's set-up runs up to the comparison line of its factorized-tau
## script, then only the factorized filter is run on its filter arguments.
run.dir <- file.path(out.dir, "runs")
if ("--worker" %in% args) {
  exp <- opt("--worker"); seed <- as.numeric(opt("--seed")); file <- opt("--out")
  t0 <- Sys.time()
  source_until(file.path("synthetic-experiments/factorized-tau", paste0(exp, ".R")),
               "^cmp  <- run_comparison", envir = globalenv())
  set.seed(seed)
  fit <- do.call(ORJointUpdateFactorized,
                 modifyList(filter_args, list(grid = grid, forecastQuantiles = FALSE, keepMix = integer(0))))
  out <- c(fit[c("wDropIS", "wDropCov", "wDrop", "nNone", "nActive", "nBad")],
           list(experiment = exp, seed = seed, minutes = as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  saveRDS(out, file)
  quit(save = "no")
}

run_one <- function(exp, seed) {
  file <- file.path(run.dir, sprintf("%s_seed%02d.rds", exp, seed))
  log  <- sub("\\.rds$", ".log", file)
  st <- system2("Rscript", c("synthetic-experiments/factorized-tau/run-seeds-dropped.R",
                             "--worker", exp, "--seed", seed, "--out", file), stdout = log, stderr = log)
  out <- if (st == 0 && file.exists(file)) readRDS(file)
         else list(experiment = exp, seed = seed, error = paste(tail(readLines(log), 3), collapse = " | "))
  message(sprintf("%-45s seed %2d  %s", exp, seed,
                  if (is.null(out$error)) sprintf("done (%.1f min)", out$minutes) else paste("FAILED:", out$error)))
  out
}

rds <- file.path(out.dir, "results.rds")
if (plot.only) {
  res <- readRDS(rds)
} else {
  dir.create(run.dir, showWarnings = FALSE)
  todo <- expand.grid(seed = seeds, exp = experiments, stringsAsFactors = FALSE)
  message(sprintf("%d runs (%d experiments x %d seeds), %d at a time", nrow(todo), length(experiments), length(seeds), jobs))
  res <- parallel::mclapply(seq_len(nrow(todo)), function(i) run_one(todo$exp[i], todo$seed[i]),
                            mc.cores = jobs, mc.preschedule = FALSE)
  saveRDS(res, rds)
}

## ---- figures and summary -----------------------------------------------------------------
## seeds x years heatmap; zero cells white, labelled cells optional
heat <- function(M, zlim, cols, main, key, labels = NULL) {
  ny <- ncol(M); ns <- nrow(M)
  layout(matrix(1:2, 1), widths = c(1, 0.12))
  par(mar = c(4, 4.5, 3.2, 0.5))
  Z <- M; Z[M == 0 | is.na(M)] <- NA
  image(1:ny, 1:ns, t(Z), zlim = zlim, col = cols, xlab = "year", ylab = "seed", yaxt = "n", main = main)
  axis(2, at = 1:ns, labels = rownames(M), las = 1, cex.axis = 0.65)
  abline(h = 1:ns + 0.5, v = 1:ny + 0.5, col = "grey92")
  box()
  if (!is.null(labels) && any(!is.na(Z))) { ok <- !is.na(Z); text(col(M)[ok], row(M)[ok], labels[ok], cex = 0.55) }
  par(mar = c(4, 0.5, 3.2, 3))
  image(1, seq(zlim[1], zlim[2], length.out = length(cols)), t(seq_along(cols)), col = cols, xaxt = "n", yaxt = "n", xlab = "", ylab = "")
  axis(4, las = 1, cex.axis = 0.75); mtext(key, side = 4, line = 2, cex = 0.85)
  layout(1)
}

ok  <- res[sapply(res, function(r) is.null(r$error))]
bad <- res[sapply(res, function(r) !is.null(r$error))]
mass.cols  <- hcl.colors(40, "YlOrRd", rev = TRUE)
count.cols <- hcl.colors(20, "Purples", rev = TRUE)
lines.out  <- c(sprintf("Factorized-tau filter: dropped forecast components across seeds (%s)",
                        format(Sys.time(), "%Y-%m-%d %H:%M")), "")

for (exp in experiments) {
  rs <- ok[sapply(ok, `[[`, "experiment") == exp]
  if (!length(rs)) next
  rs <- rs[order(sapply(rs, `[[`, "seed"))]
  sd <- sapply(rs, `[[`, "seed")
  get <- function(f) { m <- t(sapply(rs, `[[`, f)); rownames(m) <- sd; m }
  IS <- get("wDropIS"); CV <- get("wDropCov"); NN <- get("nNone"); NA_ <- get("nActive"); NB <- get("nBad")
  TOT <- IS + CV
  lg  <- function(M) { L <- log10(M); L[M <= 0] <- 0; L }   # 0 entries are left blank below
  floorlog <- -12
  zl  <- c(floorlog, 0)
  clamp <- function(L) { L[L != 0] <- pmax(L[L != 0], floorlog); L }
  nsteps <- length(TOT)

  pdf(file.path(out.dir, paste0(exp, ".pdf")), width = 8.5, height = 6.5)
  ## page 1: overview
  par(mar = c(1, 1, 1, 1))
  plot.new()
  txt <- c(exp, "",
           sprintf("factorized-tau filter, %d seeds x %d years (%d filter steps)", length(sd), ncol(TOT), nsteps), "",
           sprintf("steps with any dropped mass:            %d (%d seeds)", sum(TOT > 0), sum(rowSums(TOT) > 0)),
           sprintf("  from collapsed importance weights:    %d", sum(IS > 0)),
           sprintf("  from singular covariances:            %d (%d components in all)", sum(CV > 0), sum(NB)),
           sprintf("largest mass dropped in one step:       %.3g", max(TOT)),
           sprintf("mean mass dropped per step:             %.3g", mean(TOT)), "",
           sprintf("grid points with zero mass after dropping: %d in all, in %d steps (%d seeds)",
                   sum(NN), sum(NN > 0), sum(rowSums(NN) > 0)),
           sprintf("active grid points per step:            %d - %d", min(NA_), max(NA_)), "",
           "Dropped mass: share of the forecast probability (over the tau blocks, history",
           "groups and scenarios) on components given zero weight before the update.",
           "Zero-mass grid points: grid points all of whose components were dropped; they",
           "are restored at their prior weights (on unweighted placeholder moments).")
  text(0, 1, paste(txt, collapse = "\n"), adj = c(0, 1), family = "mono", cex = 0.8)

  ## page 2: total dropped mass
  heat(clamp(lg(TOT)), zl, mass.cols, sprintf("%s\nforecast mass dropped (all causes); blank = none", exp), "log10 mass (bottom: <= -12)")
  ## page 3: by cause
  heat(clamp(lg(IS)), zl, mass.cols, "mass dropped for collapsed importance weights; blank = none", "log10 mass (bottom: <= -12)")
  heat(clamp(lg(CV)), zl, mass.cols, "mass dropped for singular covariances; blank = none", "log10 mass (bottom: <= -12)")
  ## page 4: zero-mass grid points (labels: zero-mass / active)
  heat(NN, c(0.5, max(1.5, max(NN) + 0.5)), count.cols,
       "zero-mass grid points after dropping\n(label: zero-mass / active grid points); blank = none",
       "grid points", labels = matrix(sprintf("%d/%d", NN, NA_), nrow(NN)))
  ## page 5: per-year totals across seeds
  op <- par(mfrow = c(2, 1), mar = c(4, 4.6, 2.6, 1))
  ns <- colSums(TOT > 0)
  barplot(ns, names.arg = seq_len(ncol(TOT)), col = "#fd8d3c", border = NA, ylim = c(0, max(1, ns) * 1.15),
          yaxt = "n", xlab = "year", ylab = "seeds", main = sprintf("Seeds with any dropped mass, by year (of %d seeds)", length(sd)))
  axis(2, at = unique(round(pretty(c(0, max(1, ns))))), las = 1)
  nz <- colSums(NN)
  barplot(nz, names.arg = seq_len(ncol(NN)), col = "#807dba", border = NA, ylim = c(0, max(1, nz) * 1.15),
          yaxt = "n", xlab = "year", ylab = "grid points", main = "Zero-mass grid points summed over seeds, by year")
  axis(2, at = unique(round(pretty(c(0, max(1, nz))))), las = 1)
  par(op)
  dev.off()

  lines.out <- c(lines.out, exp,
    sprintf("  seeds %s", paste(sd, collapse = ",")),
    sprintf("  steps with dropped mass %d of %d (importance %d, covariance %d); max %.3g; seeds affected: %s",
            sum(TOT > 0), nsteps, sum(IS > 0), sum(CV > 0), max(TOT),
            if (any(TOT > 0)) paste(sd[rowSums(TOT) > 0], collapse = ",") else "none"),
    sprintf("  zero-mass grid points %d in %d steps; seeds affected: %s", sum(NN), sum(NN > 0),
            if (any(NN > 0)) paste(sd[rowSums(NN) > 0], collapse = ",") else "none"),
    sprintf("  run time per seed (set-up + filter): %.2f min on average", mean(sapply(rs, `[[`, "minutes"))), "")
}
if (length(bad)) lines.out <- c(lines.out, "FAILED runs:",
  sapply(bad, function(r) sprintf("  %s seed %d: %s", r$experiment, r$seed, r$error)))
writeLines(lines.out, file.path(out.dir, "summary.txt"))
cat(lines.out, sep = "\n")
