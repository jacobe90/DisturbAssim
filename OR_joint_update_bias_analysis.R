## Bias-correction analysis for ORJointUpdate (OR_joint_update2.Rmd).
##
## Runs the joint update for one site with and without the online bias correction
## (posterior-weighted, undisturbed components only; see the `bias` argument of
## ORJointUpdate), prints summary scores, and writes a four-panel comparison figure:
##   (a) one-step-ahead forecast of wood C with 95% band, both runs, plus observations
##   (b) innovations y_t - E[y_t | y_{1:t-1}]
##   (c) bias estimate b_t and the P(z = 0)-weighted residuals r_t
##   (d) posterior probability of disturbance
##
## Run from the repo root:  Rscript OR_joint_update_bias_analysis.R
## Needs the same data files as OR_joint_update2.Rmd (STprior, extract.h03v04,
## dis.postprocessed, fit, daymet.h03v04 .RData).

## ---- settings ---------------------------------------------------------------
site      <- 300
bias.cfg  <- list(rho = 0.9, denom = "N")   # skip = 1, burnin = 0, wMin = 0 by default
outdir    <- "experiments/biasCorrection"
figfile   <- sprintf("claude-working-docs/bias_correction_site%d.png", site)
rerun     <- FALSE   # TRUE: rerun even if saved runs exist

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

## ---- run (or load) both configurations --------------------------------------
tag   <- sprintf("site_%d_distprod_TRUE_GEDI_het_taugamma", site)
files <- c(none = file.path(outdir, "none", paste0(tag, ".rds")),
           post = file.path(outdir, "post", paste0(tag, "_bias.rds")))
cfgs  <- list(none = NULL, post = bias.cfg)
for (nm in names(files)) {
  if (!rerun && file.exists(files[[nm]])) next
  set.seed(site)   # same proposal draws for both runs
  run_joint_update_sites(site, Rparams = Rparams.GEDI.het, useDistProduct = TRUE,
                         tauPrior = "gamma", Rparams.source = "GEDI", Rparams.model = "het",
                         outdir = file.path(outdir, nm), bias = cfgs[[nm]])
}
h0 <- readRDS(files[["none"]])
h1 <- readRDS(files[["post"]])
x <- h0$year
pal <- list(raw = "#8a8984", cor = "#2a78d6", dist = "#eb6834", obs = "#0b0b0b",
            grid = "#e6e6e3", ink = "#52514e")

## summary: one-step-ahead innovation y_t - E[y_t | y_{1:t-1}], excluding year 1
## (population IC prior) and the flagged disturbance years (pDist > 0.5)
keep <- x > x[1] & h0$pDist < 0.5 & h1$pDist < 0.5
inn0 <- (h0$obs - h0$predMean)[keep]; inn1 <- (h1$obs - h1$predMean)[keep]
pit0 <- h0$pit[keep]; pit1 <- h1$pit[keep]
cat(sprintf("%-14s %8s %8s %8s %8s %8s\n", "", "meanInn", "RMSE", "meanPIT", "MAE_ana", "probD"))
for (nm in c("none", "post")) {
  h <- if (nm == "none") h0 else h1; i <- if (nm == "none") inn0 else inn1
  cat(sprintf("%-14s %8.3f %8.3f %8.3f %8.3f %8.3f\n", nm, mean(i), sqrt(mean(i^2)),
              mean(h$pit[keep]), mean(abs(h$obs - h$mean)[keep]), mean(h$pDist[keep])))
}

png(figfile, width = 9, height = 10.5, units = "in", res = 150)
layout(matrix(1:4, 4), heights = c(1.5, 1.1, 1, 0.8))
base <- function() par(mar = c(2.2, 5, 2.6, 1), las = 1, bty = "n", col.axis = pal$ink,
                       fg = pal$ink, cex.axis = 1, cex.main = 1.15, font.main = 1, adj = 0)
grid_y <- function() abline(h = axTicks(2), col = pal$grid)

## (a) one-step-ahead predictive mean +/- 1.96 sd, both runs, with obs
base()
band <- function(h, col) { sd <- sqrt(h$predVar); polygon(c(x, rev(x)),
  c(h$predMean - 1.96 * sd, rev(h$predMean + 1.96 * sd)), col = adjustcolor(col, 0.15), border = NA) }
yl <- range(c(h0$predMean, h1$predMean, h0$obs), na.rm = TRUE)
plot(NA, xlim = range(x) + c(-0.6, 0.6), ylim = yl, xlab = "", ylab = expression("wood C (kg C m"^-2*")"),
     main = sprintf("Site %d: one-step-ahead forecast of wood C (before assimilating y_t)", site))
grid_y(); band(h0, pal$raw); band(h1, pal$cor)
lines(x, h0$predMean, col = pal$raw, lwd = 2, lty = 2)
lines(x, h1$predMean, col = pal$cor, lwd = 2)
points(x, h0$obs, pch = 21, bg = pal$obs, col = "white", cex = 1.2)
legend("topright", bty = "n", text.col = pal$ink, cex = 1,
       legend = c("no bias correction", "bias-corrected", "observation"),
       col = c(pal$raw, pal$cor, pal$obs), lwd = c(2, 2, NA), lty = c(2, 1, NA), pch = c(NA, NA, 16))

## (b) innovations y_t - forecast mean
base()
i0 <- h0$obs - h0$predMean; i1 <- h1$obs - h1$predMean
plot(NA, xlim = range(x) + c(-0.6, 0.6), ylim = c(-1.6, 0.8), xlab = "", ylab = "obs - forecast",
     main = sprintf("Innovation (outside disturbance years): mean %.2f -> %.2f, RMSE %.2f -> %.2f   [bars clipped to axis]",
                    mean(inn0), mean(inn1), sqrt(mean(inn0^2)), sqrt(mean(inn1^2))))
grid_y(); abline(h = 0, col = pal$ink)
clip <- function(v) pmin(pmax(v, -1.6), 0.8)
rect(x - 0.38, 0, x - 0.02, clip(i0), col = pal$raw, border = NA)
rect(x + 0.02, 0, x + 0.38, clip(i1), col = pal$cor, border = NA)

## (c) bias trace b_t with residuals r_t (filled = full weight, hollow = w_t ~ 0)
base()
r <- h1$biasResid; w <- h1$biasW; used <- seq_along(x) > 1
plot(NA, xlim = range(x) + c(-0.6, 0.6), ylim = c(-1.2, 0.3), xlab = "", ylab = "bias (kg C m-2)",
     main = "Bias estimate b_t (line) and P(z=0)-weighted residuals r_t (points; hollow = down-weighted or skipped)")
grid_y(); abline(h = 0, col = pal$ink)
lines(x, h1$bias, type = "s", col = pal$cor, lwd = 2.5)
inr <- r >= -1.2 & r <= 0.3
full <- used & w > 0.5 & inr
points(x[full], r[full], pch = 21, bg = pal$cor, col = "white", cex = 1.3)
points(x[!full & inr], r[!full & inr], pch = 21, bg = "white", col = pal$cor, cex = 1.3)
off <- which(!inr)
if (length(off)) text(x[off], -1.2, sprintf("r=%.1f\nw=%.2f", r[off], w[off]), cex = 0.75,
                      col = pal$ink, pos = 3, xpd = NA)

## (d) P(disturbed)
par(mar = c(4, 5, 2.6, 1)); base(); par(mar = c(4, 5, 2.6, 1))
plot(NA, xlim = range(x) + c(-0.6, 0.6), ylim = c(0, 1), xlab = "", ylab = "P(disturbed)",
     main = "Posterior probability of disturbance")
abline(h = seq(0, 1, 0.5), col = pal$grid)
rect(x - 0.38, 0, x - 0.02, h0$pDist, col = pal$raw, border = NA)
rect(x + 0.02, 0, x + 0.38, h1$pDist, col = pal$cor, border = NA)
mtext("year", side = 1, line = 2.5, adj = 0.5, col = pal$ink)
dev.off()
message("wrote ", figfile)
