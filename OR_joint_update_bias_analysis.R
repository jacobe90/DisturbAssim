## Bias-correction analysis for ORJointUpdate (OR_joint_update2.Rmd).
##
## Runs the joint update for one site without bias correction and with the online
## bias correction (posterior-weighted, undisturbed components only; see the `bias`
## argument of ORJointUpdate) under each denominator:
##   denom = "N": b_{t+1} = S_t / N_t   (weighted discounted mean)
##   denom = "t": b_{t+1} = S_t / n_t   (n_t = # obs used; shrinks b toward 0)
## It prints a table of summary scores, saves it as CSV, and writes a comparison figure:
##   (a) one-step-ahead forecast of wood C with 95% band (mean +/- 1.96 sd), plus obs
##   (b) innovations y_t - E[y_t | y_{1:t-1}]
##   (c) bias estimate b_t and the P(z = 0)-weighted residuals r_t
##   (d) posterior probability of disturbance
##   (e) the summary-score table
##
## Run from the repo root:  Rscript OR_joint_update_bias_analysis.R
## Needs the same data files as OR_joint_update2.Rmd (STprior, extract.h03v04,
## dis.postprocessed, fit, daymet.h03v04 .RData).

## ---- settings ---------------------------------------------------------------
site    <- 300
## runs to compare: label -> bias config (NULL = no correction). skip = 1, burnin = 0,
## wMin = 0 by default.
cfgs <- list(
  "no correction" = NULL,
  "denom = N"     = list(rho = 0.9, denom = "N"),
  "denom = t"     = list(rho = 0.9, denom = "t")
)
## per-run output directory under outdir, and plot style
run.dir <- c("no correction" = "none", "denom = N" = "denomN", "denom = t" = "denomT")
run.col <- c("no correction" = "#8a8984", "denom = N" = "#2a78d6", "denom = t" = "#eb6834")
run.lty <- c("no correction" = 2, "denom = N" = 1, "denom = t" = 4)
outdir  <- "experiments/biasCorrection"
figfile <- sprintf("claude-working-docs/bias_correction_site%d.png", site)
tabfile <- sub("\\.png$", "_scores.csv", figfile)
rerun   <- FALSE   # TRUE: rerun even if saved runs exist

## ---- scores ------------------------------------------------------------------
## Summary scores of one run over the years in `keep`:
##   meanInn  mean one-step-ahead innovation y_t - E[y_t | y_{1:t-1}]
##   RMSE     sqrt(mean innovation^2)  (= sqrt(meanInn^2 + sdInn^2))
##   sdInn    sd of the innovations (scatter left after removing the mean)
##   meanPIT  mean PIT of y_t under the one-step-ahead predictive (0.5 = unbiased)
##   MAE_ana  mean |y_t - analysis mean|
##   probD    mean posterior P(disturbed)
score_run <- function(h, keep) {
  i <- (h$obs - h$predMean)[keep]
  c(meanInn = mean(i), RMSE = sqrt(mean(i^2)), sdInn = sqrt(mean((i - mean(i))^2)),
    meanPIT = mean(h$pit[keep]), MAE_ana = mean(abs(h$obs - h$mean)[keep]),
    probD = mean(h$pDist[keep]))
}

## ---- plotting --------------------------------------------------------------
## Comparison figure for a named list of ORJointUpdate runs H (same site; the first
## one is the uncorrected reference), written to figfile as a PNG. scores: the table
## from score_run, one row per run; keep: years it was computed over.
plot_bias_comparison <- function(H, scores, figfile, site, keep) {
  x   <- H[[1]]$year
  nm  <- names(H); nr <- length(H)
  ink <- "#52514e"; grid.col <- "#e6e6e3"; obs.col <- "#0b0b0b"
  col <- run.col[nm]; lty <- run.lty[nm]
  ## side-by-side bar slots for the runs within each year
  bw  <- 0.76 / nr
  bar <- function(v, k, ...) rect(x - 0.38 + (k - 1) * bw + 0.01, 0, x - 0.38 + k * bw - 0.01, v,
                                  col = col[k], border = NA, ...)

  png(figfile, width = 9, height = 12.5, units = "in", res = 150)
  layout(matrix(1:5, 5), heights = c(1.5, 1.1, 1, 0.8, 0.75))
  .panel_par <- function(mar = c(2.2, 5, 2.6, 1))
    par(mar = mar, las = 1, bty = "n", col.axis = ink, fg = ink, cex.axis = 1,
        cex.main = 1.15, font.main = 1, adj = 0)
  grid_y <- function() abline(h = axTicks(2), col = grid.col)
  xl <- range(x) + c(-0.6, 0.6)

  ## (a) one-step-ahead predictive mean +/- 1.96 sd, every run, with obs
  .panel_par()
  yl <- range(c(sapply(H, `[[`, "predMean"), H[[1]]$obs), na.rm = TRUE)
  plot(NA, xlim = xl, ylim = yl, xlab = "", ylab = expression("wood C (kg C m"^-2*")"),
       main = sprintf("Site %d: one-step-ahead forecast of wood C (before assimilating y_t)", site))
  grid_y()
  for (k in seq_len(nr)) {
    h <- H[[k]]; sd <- sqrt(h$predVar)
    polygon(c(x, rev(x)), c(h$predMean - 1.96 * sd, rev(h$predMean + 1.96 * sd)),
            col = adjustcolor(col[k], 0.12), border = NA)
  }
  for (k in seq_len(nr)) lines(x, H[[k]]$predMean, col = col[k], lwd = 2, lty = lty[k])
  points(x, H[[1]]$obs, pch = 21, bg = obs.col, col = "white", cex = 1.2)
  legend("topright", bty = "n", text.col = ink, cex = 1, legend = c(nm, "observation"),
         col = c(col, obs.col), lwd = c(rep(2, nr), NA), lty = c(lty, NA), pch = c(rep(NA, nr), 16))

  ## (b) innovations y_t - forecast mean
  .panel_par()
  plot(NA, xlim = xl, ylim = c(-1.6, 0.8), xlab = "", ylab = "obs - forecast",
       main = "Innovation y_t - E[y_t | y_{1:t-1}]   [bars clipped to axis; scores in table below]")
  grid_y(); abline(h = 0, col = ink)
  clip <- function(v) pmin(pmax(v, -1.6), 0.8)
  for (k in seq_len(nr)) bar(clip(H[[k]]$obs - H[[k]]$predMean), k)

  ## (c) bias traces b_t with residuals r_t (filled = used at full weight, hollow =
  ## down-weighted or skipped), for every corrected run
  .panel_par()
  plot(NA, xlim = xl, ylim = c(-1.2, 0.3), xlab = "", ylab = "bias (kg C m-2)",
       main = "Bias estimate b_t (lines) and P(z=0)-weighted residuals r_t (points; hollow = down-weighted or skipped)")
  grid_y(); abline(h = 0, col = ink)
  off.lab <- character(0)
  for (k in which(!is.na(sapply(H, function(h) h$biasW[2])))) {
    h <- H[[k]]; r <- h$biasResid; w <- h$biasW
    lines(x, h$bias, type = "s", col = col[k], lwd = 2.5, lty = lty[k])
    inr  <- !is.na(r) & r >= -1.2 & r <= 0.3
    full <- seq_along(x) > 1 & w > 0.5 & inr
    points(x[full], r[full], pch = 21, bg = col[k], col = "white", cex = 1.2)
    points(x[!full & inr], r[!full & inr], pch = 21, bg = "white", col = col[k], cex = 1.2)
    off <- which(!inr & !is.na(r))
    off.lab <- sprintf("r=%.1f\nw=%.2f", r[off], w[off]); off.x <- x[off]
  }
  ## off-axis residuals (same in every run: year 1 and the disturbance years)
  if (length(off.lab)) text(off.x, -1.2, off.lab, cex = 0.75, col = ink, pos = 3, xpd = NA)

  ## (d) P(disturbed)
  .panel_par(mar = c(4, 5, 2.6, 1))
  plot(NA, xlim = xl, ylim = c(0, 1), xlab = "", ylab = "P(disturbed)",
       main = "Posterior probability of disturbance")
  abline(h = seq(0, 1, 0.5), col = grid.col)
  for (k in seq_len(nr)) bar(H[[k]]$pDist, k)
  mtext("year", side = 1, line = 2.5, adj = 0.5, col = ink)

  ## (e) score table
  .panel_par(mar = c(0.5, 1, 2.6, 1))
  plot(NA, xlim = c(0, 1), ylim = c(0, nr + 2.2), axes = FALSE, xlab = "", ylab = "",
       main = sprintf("Scores over %d years (excluding %d and P(disturbed) > 0.5 years)",
                      sum(keep), x[1]))
  cols <- colnames(scores)
  cx <- c(0.02, 0.27 + (seq_along(cols) - 1) * 0.12)
  yy <- nr + 1.2 - seq_len(nr)
  text(cx[-1], nr + 1.2, cols, adj = c(1, 0.5), font = 2, col = ink, cex = 1.05)
  segments(0, nr + 0.7, 1, nr + 0.7, col = ink, lwd = 0.8)
  for (k in seq_len(nr)) {
    segments(cx[1], yy[k], cx[1] + 0.035, yy[k], col = col[k], lwd = 2.5, lty = lty[k])
    text(cx[1] + 0.05, yy[k], nm[k], adj = c(0, 0.5), col = ink, cex = 1.05)
    text(cx[-1], yy[k], sprintf("%.3f", scores[k, ]), adj = c(1, 0.5), col = ink, cex = 1.05)
  }
  text(0, -0.1, "meanInn, RMSE, sdInn: innovation vs forecast mean; meanPIT: 0.5 = unbiased; MAE_ana: |obs - analysis mean|; probD: mean P(disturbed)",
       adj = c(0, 0.5), col = ink, cex = 0.8, xpd = NA)
  dev.off()
  invisible(figfile)
}

## ---- load the model setup and ORJointUpdate from OR_joint_update2.Rmd ---------
## Everything before the `run-single` chunk: data, VSEM setup, make_or_args,
## ORJointUpdate and run_joint_update_sites. The run / grid chunks are not executed.
pdf(NULL)   # swallow the diagnostic plots drawn by the setup chunks
rfile <- tempfile(fileext = ".R")
knitr::purl("OR_joint_update2.Rmd", output = rfile, quiet = TRUE)
code  <- readLines(rfile)
stop.at <- grep("^## ----run-single", code)
stopifnot(length(stop.at) == 1)
source(textConnection(code[seq_len(stop.at - 1)]))
dev.off()

## ---- run (or load) every configuration ---------------------------------------
tag <- sprintf("site_%d_distprod_TRUE_GEDI_het_taugamma", site)
H   <- list()
for (nm in names(cfgs)) {
  f <- file.path(outdir, run.dir[[nm]], paste0(tag, if (is.null(cfgs[[nm]])) "" else "_bias", ".rds"))
  if (rerun || !file.exists(f)) {
    set.seed(site)   # same proposal draws for every run
    run_joint_update_sites(site, Rparams = Rparams.GEDI.het, useDistProduct = TRUE,
                           tauPrior = "gamma", Rparams.source = "GEDI", Rparams.model = "het",
                           outdir = dirname(f), bias = cfgs[[nm]])
  }
  H[[nm]] <- readRDS(f)
}

## ---- scores, table, figure -----------------------------------------------------
## Years scored: drop year 1 (forecast from the population IC prior) and any year a
## run flags as disturbed (P(disturbed) > 0.5).
x      <- H[[1]]$year
keep   <- x > x[1] & Reduce(`&`, lapply(H, function(h) h$pDist < 0.5))
scores <- t(sapply(H, score_run, keep = keep))
print(round(scores, 3))
write.csv(data.frame(run = rownames(scores), round(scores, 4), row.names = NULL),
          tabfile, row.names = FALSE)
plot_bias_comparison(H, scores, figfile, site, keep)
message("wrote ", figfile, " and ", tabfile)
