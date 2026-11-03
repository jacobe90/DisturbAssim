## Does a smaller observation-error variance R fix the over-wide forecasts?
##
## Reruns the 100-site subset of claude-forecast-width-10.8.26 (bias correction on,
## no disturbance product) for each observation-error model with R scaled by c:
##   R_t = c * (m_R y_t + b_R),   c in c.grid
## and compares one-step-ahead calibration, the learned process variance and
## disturbance detection across c. Each run uses the same seed as its grid run, so
## c = 1 reproduces experiments/fullRunBias exactly (checked below).
##
##   Rscript claude-working-docs/r_scale_experiment.R
##   (EXP_DIR=/path/to/experiments Rscript ... if the experiment folders are elsewhere)

out.dir <- "claude-working-docs"
exp.dir <- Sys.getenv("EXP_DIR", "experiments")
c.grid  <- c(1, 0.5, 0.25, 0.1, 0.05, 0.02, 0.01, 0.005)
n.sub   <- 100

## ---- model setup + grid definitions from OR_joint_update2.Rmd --------------------
pdf(NULL)
rfile <- tempfile(fileext = ".R")
knitr::purl("OR_joint_update2.Rmd", output = rfile, quiet = TRUE)
code <- readLines(rfile)
hdr  <- grep("^## ----", code); rs <- grep("^## ----run-single", code); rg <- grep("^## ----run-grid", code)
g    <- code[rg:(min(hdr[hdr > rg]) - 1)]
g    <- g[seq_len(grep("^fmt_secs", g) - 1)]          # grid definitions, not the worker pool
source(textConnection(c(code[seq_len(rs - 1)], g)))
dev.off()
NS <- nrow(dextract)

## same subset as forecast_width_analysis.R
set.seed(1); sub <- sort(sample(NS, n.sub))
Rmods <- names(Rparams.grid)
Rlab  <- c(GEDI_het = "GEDI heteroskedastic", GEDI_con = "GEDI constant",
           LandTrendr_het = "LandTrendr heteroskedastic", LandTrendr_con = "LandTrendr constant")

## ---- runs ------------------------------------------------------------------------
jobs <- expand.grid(row = which(exp.grid$site %in% sub & !exp.grid$useDistProduct), c = c.grid)
jobs$dir  <- file.path(exp.dir, "rScale", sprintf("c%g", jobs$c))
jobs$file <- with(jobs, file.path(dir, basename(exp.grid$file[row])))
todo <- which(!file.exists(jobs$file))
if (length(todo)) {
  message(sprintf("running %d runs", length(todo)))
  res <- parallel::mclapply(todo, function(j) {
    i <- jobs$row[j]; e <- exp.grid[i, ]; r <- Rparams.grid[[e$R]]
    set.seed(i)                               # same seed as the grid run of this row
    tryCatch({
      suppressMessages(run_joint_update_sites(e$site, Rparams = jobs$c[j] * r$Rparams,
        useDistProduct = e$useDistProduct, tauPrior = e$tauPrior, Rparams.source = r$source,
        Rparams.model = r$model, outdir = jobs$dir[j], bias = grid.bias))
      NULL
    }, error = function(err) conditionMessage(err))
  }, mc.cores = grid.cores)
  bad <- !sapply(res, is.null)
  if (any(bad)) { message(sum(bad), " runs failed:"); print(unique(unlist(res[bad]))) }
}
jobs <- jobs[file.exists(jobs$file), ]

## c = 1 must reproduce the grid
chk <- head(jobs[jobs$c == 1, ], 40)
cat("c = 1 reproduces fullRunBias (max |diff| of predMean over 40 runs):",
    max(sapply(seq_len(nrow(chk)), function(j)
      max(abs(readRDS(chk$file[j])$predMean -
              readRDS(file.path(exp.dir, "fullRunBias", basename(chk$file[j])))$predMean), na.rm = TRUE))), "\n")

## ---- site-years ----------------------------------------------------------------
## all observed site-years t >= 2 (PIT of the full predictive, as in the PIT plots), and
## the undisturbed subset (t >= 3, P(disturbed | y_{1:t}) < 0.5) for the k = 0 scores
sy <- do.call(rbind, lapply(seq_len(nrow(jobs)), function(j) {
  h <- readRDS(jobs$file[j]); e <- exp.grid[jobs$row[j], ]
  data.frame(R.mod = e$R, c = jobs$c[j], site = e$site, t = seq_len(NT), y = h$obs, pit = h$pit,
             m0 = h$predMean0, V0 = h$predVar0, Q0 = h$predProcVar0, R = h$predR,
             pD = h$pDist, procVar = h$procVar, bias = h$bias)
}))
sy  <- subset(sy, t >= 2 & is.finite(y))
und <- subset(sy, t >= 3 & pD < 0.5)
und$e <- und$y - und$m0; und$z <- und$e / sqrt(und$V0); und$P0 <- und$V0 - und$R - und$Q0

## calibration error of a PIT sample: mean |empirical - nominal| coverage of the
## central p-intervals, p = 0.05, ..., 0.95 (0 if calibrated; > 0 either way)
cov.err <- function(u) { p <- seq(0.05, 0.95, 0.05); u <- u[is.finite(u)]
  mean(abs(sapply(p, function(pp) mean(abs(2 * u - 1) <= pp)) - p)) }
## signed version at the central 50% interval (> 0: too wide, < 0: too narrow)
cov50 <- function(u) { u <- u[is.finite(u)]; mean(abs(2 * u - 1) <= 0.5) - 0.5 }

tab <- do.call(rbind, lapply(split(sy, list(sy$R.mod, sy$c)), function(d) {
  u <- und[und$R.mod == d$R.mod[1] & und$c == d$c[1], ]
  last <- d[d$t == NT, ]
  data.frame(R.mod = d$R.mod[1], c = d$c[1],
    calErr = cov.err(d$pit), cov50 = cov50(d$pit), in95 = 100 * mean(abs(2 * d$pit - 1) <= 0.95, na.rm = TRUE),
    meanPIT = mean(d$pit, na.rm = TRUE), pit0.05 = 100 * mean(d$pit < 0.05 | d$pit > 0.95, na.rm = TRUE),
    Ez2 = mean(u$z^2), medz = median(abs(u$z)),
    logS = mean(dnorm(u$e, 0, sqrt(u$V0), log = TRUE)),
    R = mean(u$R), Q = mean(u$Q0), P = mean(u$P0), e2 = mean(u$e^2),
    vT = median(last$procVar), pD = mean(d$pD), pD.gt.5 = 100 * mean(d$pD > 0.5),
    bT = median(last$bias))
}))
tab <- tab[order(match(tab$R.mod, Rmods), -tab$c), ]
cat("\n== calibration vs R scale c (bias-corrected, no product, 100 sites) ==\n")
cat("calErr: mean |coverage error| over central intervals (0 = calibrated); cov50: coverage error of the\n",
    "central 50% interval (> 0 too wide); pit0.05: % of PIT outside [0.05, 0.95] (10 if calibrated);\n",
    "Ez2 / medz / logS: k = 0 predictive over undisturbed site-years (1 / 0.674 / higher is better);\n",
    "R, Q, P, e2: mean variance budget and squared innovation; vT: median v_T in 2017; pD: mean P(disturbed);\n",
    "pD.gt.5: % of site-years flagged disturbed by the filter; bT: median bias b in 2017\n")
print(cbind(tab[, 1:2], signif(tab[, -(1:2)], 3)), row.names = FALSE)
write.csv(tab, file.path(out.dir, "r_scale_scores.csv"), row.names = FALSE)

## ---- figures ---------------------------------------------------------------------
ink <- "#52514e"; grid.col <- "#e6e6e3"; c.bar <- "#2a78d6"
.par <- function(...) par(las = 1, bty = "n", col.axis = ink, fg = ink, font.main = 1, ...)

## (1) PIT histograms: rows = observation-error model, columns = c
br <- seq(0, 1, 0.05); nb <- length(br) - 1
png(file.path(out.dir, "rs_pit.png"), width = 15, height = 11, units = "in", res = 140)
.par(mfrow = c(length(Rmods), length(c.grid)), mar = c(3.2, 3.4, 3, 0.6), oma = c(1.5, 2, 3, 0))
for (m in Rmods) for (cc in c.grid) {
  u <- sy$pit[sy$R.mod == m & sy$c == cc]; u <- u[is.finite(u)]
  dens <- hist(u, br, plot = FALSE)$counts * nb / length(u)
  band <- qbinom(c(0.025, 0.975), length(u), 1 / nb) * nb / length(u)
  k <- tab[tab$R.mod == m & tab$c == cc, ]
  plot(NA, xlim = c(0, 1), ylim = c(0, 4), xlab = "", ylab = "",
       main = sprintf("c = %g\nerror %.3f, %.1f%% in 95%%", cc, k$calErr, k$in95), cex.main = 0.95)
  rect(0, band[1], 1, band[2], col = grid.col, border = NA)
  rect(br[-length(br)] + 0.004, 0, br[-1] - 0.004, pmin(dens, 4), col = c.bar, border = NA)
  abline(h = 1, lty = 2, col = ink)
  if (cc == c.grid[1]) mtext(Rlab[m], side = 2, line = 2.6, las = 0, col = ink, cex = 0.85)
}
mtext("PIT u", side = 1, outer = TRUE, line = 0.2, col = ink)
mtext("One-step-ahead PIT, all site-years 1991-2017, with R scaled by c (bias-corrected, no product, 100 sites)   [dashed = uniform; gray = 95% range per bin]",
      outer = TRUE, line = 1, cex = 1.05, col = ink)
dev.off()

## (2) scores vs c, one line per observation-error model
cols <- c(GEDI_het = "#2a78d6", GEDI_con = "#1a4f99", LandTrendr_het = "#eb6834", LandTrendr_con = "#a8461f")
ltys <- c(GEDI_het = 1, GEDI_con = 2, LandTrendr_het = 1, LandTrendr_con = 2)
pnl <- list(
  list(v = "calErr", lab = "PIT calibration error (0 = calibrated)", ref = 0),
  list(v = "cov50",  lab = "coverage error of central 50% interval\n(> 0 too wide, < 0 too narrow)", ref = 0),
  list(v = "logS",   lab = "mean log predictive density\n(higher is better)", ref = NA),
  list(v = "Q",      lab = "mean process variance Q", ref = NA, log = "y"),
  list(v = "pD.gt.5", lab = "% of site-years with P(disturbed) > 0.5", ref = NA))
png(file.path(out.dir, "rs_scores.png"), width = 16, height = 4.2, units = "in", res = 140)
.par(mfrow = c(1, length(pnl)), mar = c(4.4, 4.6, 4.2, 0.8))
for (p in pnl) {
  yl <- range(tab[[p$v]], p$ref, na.rm = TRUE)
  plot(NA, xlim = rev(range(c.grid)), ylim = yl, log = paste0("x", if (!is.null(p$log)) p$log else ""),
       xlab = "R scale c", ylab = "", main = p$lab, cex.main = 1)
  abline(h = axTicks(2), col = grid.col)
  if (!is.na(p$ref)) abline(h = p$ref, lty = 2, col = ink)
  for (m in Rmods) { k <- tab[tab$R.mod == m, ]
    lines(k$c, k[[p$v]], col = cols[m], lwd = 2, lty = ltys[m]); points(k$c, k[[p$v]], pch = 19, col = cols[m]) }
  if (p$v == "calErr") legend("topright", bty = "n", cex = 0.8, text.col = ink, legend = Rlab[Rmods],
                              col = cols[Rmods], lty = ltys[Rmods], lwd = 2)
}
dev.off()
message("wrote ", file.path(out.dir, c("rs_pit.png", "rs_scores.png", "r_scale_scores.csv")), collapse = " ")
