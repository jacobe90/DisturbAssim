## Load saved ORJointUpdate runs and plot one site.
##
## Ported from OR_joint_update2.Rmd (now legacy/models/).

## Load a saved ORJointUpdate run. The site index and run settings are parsed from
## the file name (site_<idx>_distprod_<TF>_<source>_<model>_tau<prior>.rds). Runs
## store their own year axis ($year) and aligned observations ($obs); runs made
## before that (which put 1990 biomass at a 1985 timestep) are rejected.
load_hist <- function(file) {
  h     <- readRDS(file)
  tag   <- sub("\\.rds$", "", basename(file))
  site  <- as.integer(sub("^site_(\\d+)_.*$", "\\1", tag))
  if (is.null(h$year)) stop(file, ": no $year axis -- stale run from before the 1990 alignment fix; rerun it")
  h$site <- site
  h$tag  <- tag
  h
}

## shared look: one palette, light grid, no box
.pal <- list(mean = "#2a78d6", band = "#2a78d6", obs = "#0b0b0b", dist = "#eb6834",
             grid = "#e6e6e3", ink = "#52514e")
.base_par <- function() par(mar = c(4.2, 4.8, 3.2, 1), las = 1, bty = "n",
                            col.axis = .pal$ink, fg = .pal$ink, cex.axis = 0.9)

## Wood-pool analysis: mean (line), 95% interval mean +/- 1.96 sd (band), observations (points)
plot_analysis <- function(h, main = sprintf("Site %d: wood-pool analysis", h$site)) {
  op <- .base_par(); on.exit(par(op))
  x  <- h$year; sd <- sqrt(pmax(h$var, 0))
  lo <- h$mean - 1.96 * sd; hi <- h$mean + 1.96 * sd
  plot(NA, xlim = range(x), ylim = range(c(lo, hi, h$obs), na.rm = TRUE),
       xlab = "year", ylab = expression("wood carbon (kg C m"^-2*")"), main = main)
  abline(h = axTicks(2), col = .pal$grid)
  polygon(c(x, rev(x)), c(lo, rev(hi)), col = adjustcolor(.pal$band, 0.2), border = NA)
  lines(x, h$mean, col = .pal$mean, lwd = 2.5)
  points(x, h$obs, pch = 21, bg = .pal$obs, col = "white", cex = 1.1)
  legend("topright", bty = "n", cex = 0.85, text.col = .pal$ink,
         legend = c("analysis mean", "95% interval", "observation"),
         col = c(.pal$mean, adjustcolor(.pal$band, 0.35), .pal$obs),
         lwd = c(2.5, 10, NA), pch = c(NA, NA, 16))
  invisible(h)
}

## Posterior probability of disturbance per year (bars), with the 0.5 line
plot_pdist <- function(h, main = sprintf("Site %d: probability of disturbance", h$site)) {
  op <- .base_par(); on.exit(par(op))
  x <- h$year; p <- h$pDist
  plot(NA, xlim = range(x) + c(-0.6, 0.6), ylim = c(0, 1), xlab = "year",
       ylab = "P(disturbed)", main = main)
  abline(h = seq(0, 1, 0.25), col = .pal$grid)
  abline(h = 0.5, col = .pal$ink, lty = 2)
  rect(x - 0.38, 0, x + 0.38, p, col = .pal$dist, border = NA)
  big <- which(p >= 0.5)
  if (length(big)) text(x[big], p[big], sprintf("%.2f", p[big]), pos = 3, cex = 0.75,
                        col = .pal$ink, xpd = NA)
  invisible(h)
}
