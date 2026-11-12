## Figures and numbers for claude-filter-explainer-11.12.26.tex.
##
## 1. Synthetic: singlesite-singledisturb-multipool.R (the multi-pool port of ORJointUpdate,
##    known truth, disturbance in year 10) for seeds 42, 1..20.
## 2. OR: ORJointUpdate (or-model/) with the default settings (GEDI heteroskedastic R,
##    bias correction, process error through the disturbance, gamma tau prior) on site 300
##    and on 100 random sites, with and without the disturbance product.
##
## Run from the repository root, after  Rscript or-model/load-data-or.R :
##   Rscript claude-working-docs/claude-filter-explainer-11.12.26-figures.R [--jobs 8]

out.dir <- "claude-working-docs/claude-filter-explainer-11.12.26-files"
run.dir <- "experiments/claude-filter-explainer"           # model runs (not tracked)
dir.create(out.dir, showWarnings = FALSE); dir.create(run.dir, recursive = TRUE, showWarnings = FALSE)
args <- commandArgs(trailingOnly = TRUE)
jobs <- if (length(i <- which(args == "--jobs"))) as.integer(args[i + 1]) else max(1L, parallel::detectCores() - 2L)

ink <- "#3d3d3a"; grid.col <- "#e6e6e3"; c.an <- "#2a78d6"; c.fc <- "#eb6834"; c.tr <- "#0b0b0b"; c.dist <- "#c0392b"
base_par <- function(...) par(las = 1, bty = "n", col.axis = ink, fg = ink, col.lab = ink, col.main = ink, mgp = c(2.4, 0.7, 0), ...)
band <- function(x, lo, hi, col) polygon(c(x, rev(x)), c(lo, rev(hi)), border = NA, col = col)
dmix <- function(x, mix, a) sapply(x, function(v) sum(mix$w * dnorm(v, mix$mu[, a], sqrt(mix$sig[, 4 * a - 3]))))
## quantiles of a discrete distribution (values v, weights w)
wq <- function(v, w, p) { o <- order(v); cw <- cumsum(w[o]) / sum(w); sapply(p, function(pp) v[o][which(cw >= pp)[1]]) }

## ---- 1. synthetic ---------------------------------------------------------------------
syn.file <- file.path(run.dir, "synthetic-seeds.rds")
if (!file.exists(syn.file)) {
  run_seed <- function(s) {
    e <- new.env(); e$seed <- s; e$outdir <- file.path(run.dir, sprintf("synthetic-seed%02d", s))
    capture.output(sys.source("synthetic-experiments/singlesite-singledisturb-multipool.R", envir = e))
    f <- e$or.fit
    keep <- list(seed = s, truth = e$True_states, obs = e$Observed_states, Q = e$Q, R = e$R,
                 dist.at = which(e$disturbance_at == 1), p.dist = e$p.dist,
                 q = f$q, fq = f$fq, pDist = f$pDist, pit = f$pit, procVar = f$procVar,
                 tauAxes = f$tauAxes, tauMarg = f$tauMarg,
                 q.known = e$or.known$q, pDist.known = e$or.known$pDist, an.stats = e$an.stats)
    if (s == 42) { keep$mixF <- f$mixF; keep$mixA <- f$mixA }
    unlink(file.path(e$outdir, "results.rds"))
    keep
  }
  message("synthetic: 21 seeds")
  syn <- parallel::mclapply(c(42, 1:20), run_seed, mc.cores = jobs, mc.preschedule = FALSE)
  saveRDS(syn, syn.file)
}
syn <- readRDS(syn.file)
s42 <- syn[[1]]; stopifnot(s42$seed == 42)
NY <- nrow(s42$truth); yr <- 1:NY; td <- s42$dist.at
pools <- c("leaf", "soil", "wood")

## states, seed 42
pdf(file.path(out.dir, "syn-states.pdf"), width = 10, height = 5.2)
base_par(mfrow = c(2, 2), mar = c(3.6, 4.2, 2.2, 1))
for (a in 1:3) {
  q <- s42$q[, , a]
  plot(NA, xlim = range(yr), ylim = range(c(q[c(1, 5), ], s42$truth[, a], s42$obs[, a]), na.rm = TRUE),
       xlab = "year", ylab = "kg C m-2", main = sprintf("%s carbon", pools[a]))
  abline(h = axTicks(2), col = grid.col); abline(v = td, lty = 2, col = c.dist)
  band(yr, q[1, ], q[5, ], adjustcolor(c.an, 0.22))
  lines(yr, q[3, ], col = c.an, lwd = 2.2)
  lines(yr, s42$truth[, a], col = c.tr, lwd = 1.4, lty = 1)
  points(yr, s42$obs[, a], pch = 21, bg = "white", col = ink, cex = 0.8)
  if (a == 1) legend("topright", bty = "n", cex = 0.8, legend = c("truth", "observation", "analysis median", "95% credible band"),
                     col = c(c.tr, ink, c.an, adjustcolor(c.an, 0.4)), lwd = c(1.4, NA, 2.2, 8), pch = c(NA, 21, NA, NA))
}
plot(NA, xlim = range(yr) + c(-0.6, 0.6), ylim = c(0, 1), xlab = "year", ylab = "P(disturbed | y)", main = "probability of disturbance")
abline(h = c(0, 0.5, 1), col = grid.col)
rect(yr - 0.38, 0, yr + 0.38, s42$pDist, col = c.fc, border = NA)
abline(v = td, lty = 2, col = c.dist); text(td, 1, "true disturbance", pos = 4, cex = 0.8, col = c.dist, xpd = NA)
dev.off()

## forecast and analysis mixtures of wood at the disturbance year, seed 42
pdf(file.path(out.dir, "syn-mixture.pdf"), width = 8, height = 4)
base_par(mar = c(4, 4.2, 1, 1))
mf <- s42$mixF[[td]]; ma <- s42$mixA[[td]]
x <- seq(-1, max(s42$q[5, td - 1, 3], s42$obs[td - 1, 3]) * 1.15, length.out = 800)
dF <- dmix(x, mf, 3); dA <- dmix(x, ma, 3)
plot(NA, xlim = range(x), ylim = c(0, max(dF, dA) * 1.08), xlab = "wood carbon (kg C m-2)", ylab = "density")
abline(h = 0, col = grid.col)
polygon(c(x, rev(x)), c(dF, rep(0, length(x))), col = adjustcolor(c.fc, 0.25), border = c.fc)
polygon(c(x, rev(x)), c(dA, rep(0, length(x))), col = adjustcolor(c.an, 0.25), border = c.an)
abline(v = s42$obs[td, 3], lwd = 2, col = ink); abline(v = s42$truth[td, 3], lty = 2, col = ink)
xd <- x[x < 2][which.max(dF[x < 2])]; xu <- x[x > 2][which.max(dF[x > 2])]
text(xd + 0.15, max(dF[x < 2]) + 0.12, "disturbed mode\n(prior weight 5%)", adj = c(0, 0), cex = 0.8, col = c.fc)
text(xu, max(dF[x > 2]) + 0.12, "undisturbed mode", adj = c(0.5, 0), cex = 0.8, col = c.fc)
legend("topright", bty = "n", cex = 0.85, legend = c("forecast (prior) mixture", "analysis (posterior) mixture", "observation", "truth"),
       fill = c(adjustcolor(c.fc, 0.25), adjustcolor(c.an, 0.25), NA, NA), border = c(c.fc, c.an, NA, NA),
       lty = c(NA, NA, 1, 2), lwd = c(NA, NA, 2, 1), col = c(NA, NA, ink, ink))
dev.off()

## process sd posterior over time, seed 42
pdf(file.path(out.dir, "syn-tau.pdf"), width = 10, height = 3.6)
base_par(mfrow = c(1, 3), mar = c(4, 4.2, 2.2, 1))
for (a in 1:3) {
  sdv <- 1 / sqrt(s42$tauAxes[[a]])
  qq  <- apply(s42$tauMarg[[a]], 2, function(w) wq(sdv, w, c(0.025, 0.5, 0.975)))
  plot(NA, xlim = range(yr), ylim = range(sdv), log = "y", xlab = "year", ylab = "process sd", main = sprintf("%s: process sd posterior", pools[a]))
  abline(h = 10^(-2:1), col = grid.col)
  band(yr, qq[1, ], qq[3, ], adjustcolor(c.an, 0.22)); lines(yr, qq[2, ], col = c.an, lwd = 2.2)
  lines(yr, sqrt(s42$procVar[, a]), col = c.an, lty = 3)
  abline(h = sqrt(s42$Q[a, a]), lwd = 1.5, col = c.tr, lty = 2)
  if (a == 1) legend("bottomleft", bg = "white", box.col = NA, cex = 0.75, legend = c("posterior median", "95% credible band", "sqrt(E[1/tau])", "truth"),
                     col = c(c.an, adjustcolor(c.an, 0.4), c.an, c.tr), lwd = c(2.2, 8, 1, 1.5), lty = c(1, 1, 3, 2))
}
dev.off()

## across seeds
cover <- function(r, q) sapply(1:3, function(a) mean(r$truth[, a] >= r[[q]][1, , a] & r$truth[, a] <= r[[q]][5, , a]))
rmse  <- function(r, q) sapply(1:3, function(a) sqrt(mean((r[[q]][3, , a] - r$truth[, a])^2)))
syn.tab <- data.frame(
  seed = sapply(syn, `[[`, "seed"),
  pDet = sapply(syn, function(r) r$pDist[r$dist.at]),
  pFA  = sapply(syn, function(r) max(r$pDist[-r$dist.at])),
  nFA  = sapply(syn, function(r) sum(r$pDist[-r$dist.at] > 0.5)),
  t(sapply(syn, cover, "q")), t(sapply(syn, function(r) rmse(r, "q") / rmse(r, "q.known"))),
  t(sapply(syn, function(r) sqrt(r$procVar[NY, ]) / sqrt(diag(r$Q)))),
  ## posterior median of the process sd / truth, and whether the 95% band holds the truth, at the last year
  t(sapply(syn, function(r) sapply(1:3, function(a) wq(1 / sqrt(r$tauAxes[[a]]), r$tauMarg[[a]][, NY], 0.5) / sqrt(r$Q[a, a])))),
  t(sapply(syn, function(r) sapply(1:3, function(a) {
    b <- wq(1 / sqrt(r$tauAxes[[a]]), r$tauMarg[[a]][, NY], c(0.025, 0.975)); sd0 <- sqrt(r$Q[a, a]); sd0 >= b[1] & sd0 <= b[2] }))))
names(syn.tab)[5:19] <- c(paste0("cov.", pools), paste0("rmseRatio.", pools), paste0("sdRatio.", pools),
                          paste0("sdMedRatio.", pools), paste0("sdCov.", pools))
write.csv(syn.tab, file.path(out.dir, "synthetic-seeds.csv"), row.names = FALSE)

pdf(file.path(out.dir, "syn-seeds.pdf"), width = 10, height = 3.8)
base_par(mfrow = c(1, 2), mar = c(4, 4.2, 2.2, 1))
P <- t(sapply(syn, `[[`, "pDist"))
plot(NA, xlim = range(yr), ylim = c(0, 1), xlab = "year", ylab = "P(disturbed | y)", main = "P(disturbed), 21 seeds")
abline(h = c(0, 0.5, 1), col = grid.col); abline(v = td, lty = 2, col = c.dist)
for (i in seq_len(nrow(P))) points(yr + runif(NY, -0.15, 0.15), P[i, ], pch = 19, cex = 0.5, col = adjustcolor(ifelse(yr == td, c.dist, ink), 0.6))
text(td + 0.3, 0.75, "true disturbance", pos = 4, cex = 0.8, col = c.dist)
cov <- as.matrix(syn.tab[, paste0("cov.", pools)])
boxplot(cov, names = pools, ylim = c(0, 1), col = adjustcolor(c.an, 0.25), border = c.an, ylab = "share of years",
        main = "truth inside the 95% band, per seed")
abline(h = 0.95, lty = 2, col = ink)
dev.off()

## ---- 2. OR ----------------------------------------------------------------------------
source("or-model/R/load-or-model.R")
bias.cfg <- list(rho = 0.9, denom = "N", skip = 1)
set.seed(1); sub <- sort(sample(nrow(dextract), 100))
eg <- make_exp_grid(R = "GEDI_het", useDistProduct = c(TRUE, FALSE), sites = sort(unique(c(300, sub))),
                    outdir = file.path(run.dir, "or"), bias = bias.cfg)
res <- run_or_grid(eg, cores = jobs, bias = bias.cfg)
if (any(sapply(res, inherits, "try-error"))) stop("OR runs failed")
H <- lapply(eg$file, readRDS)

## site 300, with the product
h <- H[[which(eg$site == 300 & eg$useDistProduct)]]
h0 <- H[[which(eg$site == 300 & !eg$useDistProduct)]]
x <- h$year
pdf(file.path(out.dir, "or-site300.pdf"), width = 10, height = 5.2)
base_par(mfrow = c(2, 2), mar = c(3.6, 4.2, 2.2, 1))
plot(NA, xlim = range(x), ylim = range(c(h$ci95, h$obs), na.rm = TRUE), xlab = "year", ylab = "kg C m-2", main = "wood carbon: analysis")
abline(h = axTicks(2), col = grid.col)
band(x, h$ci95[, "lo"], h$ci95[, "hi"], adjustcolor(c.an, 0.22)); lines(x, h$mean, col = c.an, lwd = 2.2)
points(x, h$obs, pch = 21, bg = "white", col = ink, cex = 0.8)
legend("topright", bty = "n", cex = 0.8, legend = c("observation", "analysis mean", "95% credible band"),
       col = c(ink, c.an, adjustcolor(c.an, 0.4)), lwd = c(NA, 2.2, 8), pch = c(21, NA, NA))
pm <- h$predMean; ps <- sqrt(h$predVar)
plot(NA, xlim = range(x), ylim = range(c(pm - 2 * ps, pm + 2 * ps, h$obs), na.rm = TRUE), xlab = "year", ylab = "kg C m-2",
     main = "wood carbon: one-step-ahead forecast")
abline(h = axTicks(2), col = grid.col)
band(x, pm - 1.96 * ps, pm + 1.96 * ps, adjustcolor(c.fc, 0.22)); lines(x, pm, col = c.fc, lwd = 2.2)
points(x, h$obs, pch = 21, bg = "white", col = ink, cex = 0.8)
plot(NA, xlim = range(x) + c(-0.6, 0.6), ylim = c(0, 1), xlab = "year", ylab = "P(disturbed | y)", main = "probability of disturbance")
abline(h = c(0, 0.5, 1), col = grid.col)
rect(x - 0.4, 0, x, h$pDist, col = c.fc, border = NA); rect(x, 0, x + 0.4, h0$pDist, col = "grey55", border = NA)
fl <- dtime[dextract[300, ] > 0]; if (length(fl)) points(fl, rep(1.04, length(fl)), pch = 25, bg = c.dist, col = c.dist, xpd = NA)
legend("left", bty = "n", cex = 0.8, legend = c("with product", "without product", "product flag"),
       fill = c(c.fc, "grey55", NA), border = NA, pch = c(NA, NA, 25), col = c(NA, NA, c.dist), pt.bg = c(NA, NA, c.dist))
plot(x, sqrt(h$procVar), type = "l", lwd = 2.2, col = c.an, log = "y", xlab = "year", ylab = "process sd (kg C m-2)",
     main = "wood process sd, sqrt(E[1/tau | y])")
abline(h = axTicks(2), col = grid.col)
dev.off()

## 100 random sites: calibration and detection
pit <- lapply(c(TRUE, FALSE), function(p) unlist(lapply(H[eg$useDistProduct == p & eg$site %in% sub], function(h) h$pit[-1])))
loss <- function(s) { d <- -diff(Yw[s, ]); max(d) }            # largest one-year drop in observed wood
pmax <- function(p) sapply(sub, function(s) max(H[[which(eg$site == s & eg$useDistProduct == p)]]$pDist))
or.tab <- data.frame(site = sub, loss = sapply(sub, loss), pmax.prod = pmax(TRUE), pmax.noprod = pmax(FALSE),
                     flagged = sapply(sub, function(s) any(dextract[s, as.character(paste0("t", intersect(dtime, atime)))] > 0)))
write.csv(or.tab, file.path(out.dir, "or-100-sites.csv"), row.names = FALSE)

pdf(file.path(out.dir, "or-100.pdf"), width = 10, height = 3.9)
base_par(mfrow = c(1, 2), mar = c(4, 4.2, 2.2, 1))
br <- seq(0, 1, 0.1)
hp <- hist(pit[[1]], breaks = br, plot = FALSE); hn <- hist(pit[[2]], breaks = br, plot = FALSE)
plot(NA, xlim = c(0, 1), ylim = c(0, max(hp$density, hn$density) * 1.1), xlab = "PIT of y_t under its one-step-ahead forecast",
     ylab = "density", main = "forecast calibration (flat = calibrated)")
abline(h = 1, lty = 2, col = ink)
rect(br[-11], 0, br[-1] - 0.05, hp$density, col = c.fc, border = NA); rect(br[-11] + 0.05, 0, br[-1], hn$density, col = "grey55", border = NA)
legend("topright", bty = "n", cex = 0.8, legend = c("with product", "without product"), fill = c(c.fc, "grey55"), border = NA)
plot(or.tab$loss, or.tab$pmax.noprod, pch = 19, col = adjustcolor(ifelse(or.tab$flagged, c.dist, ink), 0.7), ylim = c(0, 1),
     xlab = "largest one-year drop in observed wood (kg C m-2)", ylab = "max over years of P(disturbed | y)", main = "detection without the product")
abline(h = 0.5, lty = 2, col = grid.col)
legend("bottomright", bty = "n", cex = 0.8, legend = c("flagged by the product", "not flagged"), pch = 19, col = c(c.dist, ink))
dev.off()

## ---- numbers for the slides ---------------------------------------------------------------
f2 <- function(x) sprintf("%.2f", x)
nums <- c(
  sprintf("\\newcommand{\\SynDetMed}{%s}", f2(median(syn.tab$pDet))),
  sprintf("\\newcommand{\\SynDetMin}{%s}", f2(min(syn.tab$pDet))),
  sprintf("\\newcommand{\\SynDetN}{%d}", sum(syn.tab$pDet > 0.5)),
  sprintf("\\newcommand{\\SynFAN}{%d}", sum(syn.tab$nFA)),
  sprintf("\\newcommand{\\SynFASeeds}{%d}", sum(syn.tab$nFA > 0)),
  sprintf("\\newcommand{\\SynYears}{%d}", length(syn) * (NY - 1)),
  sprintf("\\newcommand{\\SynCov}{%s}", paste(sprintf("%s %.0f\\%%", pools, 100 * colMeans(cov)), collapse = ", ")),
  sprintf("\\newcommand{\\SynRMSE}{%s}", paste(sprintf("%s %.2f", pools, apply(syn.tab[, paste0("rmseRatio.", pools)], 2, median)), collapse = ", ")),
  sprintf("\\newcommand{\\SynSDMed}{%s}", paste(sprintf("%s %.2f", pools, apply(syn.tab[, paste0("sdMedRatio.", pools)], 2, median)), collapse = ", ")),
  sprintf("\\newcommand{\\SynSDCov}{%s}", paste(sprintf("%s %d", pools, colSums(syn.tab[, paste0("sdCov.", pools)])), collapse = ", ")),
  sprintf("\\newcommand{\\SynSD}{%s}", paste(sprintf("%s %.2f", pools, apply(syn.tab[, paste0("sdRatio.", pools)], 2, median)), collapse = ", ")),
  sprintf("\\newcommand{\\ORPitIn}{%.0f\\%%}", 100 * mean(abs(2 * pit[[1]] - 1) <= 0.95)),
  sprintf("\\newcommand{\\ORPitInNo}{%.0f\\%%}", 100 * mean(abs(2 * pit[[2]] - 1) <= 0.95)),
  sprintf("\\newcommand{\\ORPitMid}{%.0f\\%%}", 100 * mean(pit[[1]] > 0.25 & pit[[1]] < 0.75)),
  sprintf("\\newcommand{\\ORDetBig}{%d of %d}", sum(or.tab$pmax.noprod[or.tab$loss > 3] > 0.5), sum(or.tab$loss > 3)),
  sprintf("\\newcommand{\\ORDetSmall}{%d of %d}", sum(or.tab$pmax.noprod[or.tab$loss < 1] > 0.5), sum(or.tab$loss < 1)),
  sprintf("\\newcommand{\\ORSiteSD}{%s}", f2(sqrt(h$procVar[length(x)]))),
  sprintf("\\newcommand{\\ORSiteP}{%s}", f2(max(h$pDist))),
  sprintf("\\newcommand{\\ORSitePNo}{%s}", f2(max(h0$pDist))))
writeLines(nums, file.path(out.dir, "numbers.tex"))
cat(nums, sep = "\n")
print(summary(syn.tab))
