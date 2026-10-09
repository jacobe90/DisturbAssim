## Why are the bias-corrected one-step-ahead forecasts too wide?
## Figures and numbers for claude-working-docs/claude-forecast-width-10.8.26.tex.
##
## Inputs (run from the repo root):
##   experiments/fullRunBias    full grid, bias correction on   (OR_joint_update2.Rmd run-grid)
##   experiments/fullRunNoBias  full grid, bias correction off  (same, grid.bias <- NULL)
## Runs a 100-site subset of both grids again into experiments/forecastWidth, with the
## same seeds as the grid, to get the no-disturbance predictive (predMean0 / predVar0 /
## predProcVar0), which the grid runs predate.
##
##   Rscript claude-working-docs/forecast_width_analysis.R
##   (EXP_DIR=/path/to/experiments Rscript ... if the grid outputs are elsewhere)

out.dir <- "claude-working-docs"
exp.dir <- Sys.getenv("EXP_DIR", "experiments")   # override if the runs live elsewhere
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

exps <- expand.grid(R = names(Rparams.grid), useDistProduct = c(FALSE, TRUE), stringsAsFactors = FALSE)
exps$label <- sprintf("%s %s R, %s", sapply(Rparams.grid[exps$R], `[[`, "source"),
                      ifelse(grepl("het", exps$R), "heteroskedastic", "constant"),
                      ifelse(exps$useDistProduct, "product", "no product"))
run_file <- function(dir, site, R, prod, bias)
  file.path(exp.dir, dir, sprintf("site_%d_distprod_%s_%s_%s_taugamma%s.rds", site, prod,
            Rparams.grid[[R]]$source, Rparams.grid[[R]]$model, if (bias) "_bias" else ""))

## product-flagged site-years (prior P(disturbed) > 0.5), NS x NT
flag <- t(sapply(seq_len(nrow(dextract)), function(s) 1 - make_or_args(s, TRUE)$disturb.params$p[1, ] > 0.5))
NS <- nrow(dextract)

ink <- "#52514e"; grid.col <- "#e6e6e3"; c.raw <- "#8a8984"; c.cor <- "#2a78d6"
.par <- function(...) par(las = 1, bty = "n", col.axis = ink, fg = ink, font.main = 1, ...)

## =================================================================================
## 1. learned process variance v_t = E[1/tau | y_{1:t}] over time, both grids
## =================================================================================
procvar <- function(dir, R, prod, bias)
  t(sapply(seq_len(NS), function(s) readRDS(run_file(dir, s, R, prod, bias))$procVar))
PV <- lapply(seq_len(nrow(exps)), function(i) list(
  raw = procvar("fullRunNoBias", exps$R[i], exps$useDistProduct[i], FALSE),
  cor = procvar("fullRunBias",   exps$R[i], exps$useDistProduct[i], TRUE)))
yr <- atime

png(file.path(out.dir, "fw_procvar_time.png"), width = 14, height = 7.5, units = "in", res = 150)
.par(mfrow = c(2, 4), mar = c(4, 4.6, 3.4, 0.8), oma = c(0, 0, 2.6, 0))
yl <- range(sapply(PV, function(p) apply(rbind(p$raw, p$cor), 2, quantile, c(.1, .9))))
for (i in seq_len(nrow(exps))) {
  plot(NA, xlim = range(yr), ylim = yl, log = "y", xlab = "year",
       ylab = if (i %% 4 == 1) expression(v[t] == E*"[1/"*tau*" | "*y[1:t]*"]  (kg C m"^-2*")"^2) else "",
       main = exps$label[i], cex.main = 1.05)
  abline(h = 10^(-4:1), col = grid.col)
  for (k in c("raw", "cor")) {
    q <- apply(PV[[i]][[k]], 2, quantile, c(.25, .5, .75))
    col <- if (k == "raw") c.raw else c.cor
    polygon(c(yr, rev(yr)), c(q[1, ], rev(q[3, ])), col = adjustcolor(col, 0.2), border = NA)
    lines(yr, q[2, ], col = col, lwd = 2.2, lty = if (k == "raw") 2 else 1)
  }
  if (i == 1) legend("topright", bty = "n", cex = 0.9, text.col = ink,
                     legend = c("uncorrected: median, IQR", "bias-corrected: median, IQR"),
                     col = c(c.raw, c.cor), lwd = 2.2, lty = c(2, 1))
}
mtext("Learned process variance over time across 538 sites (posterior mean of 1/tau after assimilating y_t)",
      outer = TRUE, line = 0.8, cex = 1.1, col = ink)
dev.off()

cat("\n== 1. process variance v_T in 2017 (median over sites), and share of sites where corrected < uncorrected ==\n")
pv.tab <- t(sapply(seq_len(nrow(exps)), function(i) {
  a <- PV[[i]]$raw[, NT]; b <- PV[[i]]$cor[, NT]
  c(raw.median = median(a), cor.median = median(b), median.ratio = median(b / a),
    pct.smaller = 100 * mean(b < a), raw.1995 = median(PV[[i]]$raw[, 6]), cor.1995 = median(PV[[i]]$cor[, 6]))
}))
rownames(pv.tab) <- exps$label; print(signif(pv.tab, 3))

## =================================================================================
## 2. variance budget of the no-disturbance predictive (100-site subset rerun)
## =================================================================================
set.seed(1); sub <- sort(sample(NS, n.sub))
rows <- which(exp.grid$site %in% sub)
jobs <- expand.grid(row = rows, bias = c(FALSE, TRUE))
jobs$file <- mapply(function(r, b) run_file(file.path("forecastWidth", if (b) "bias" else "nobias"),
                    exp.grid$site[r], exp.grid$R[r], exp.grid$useDistProduct[r], b), jobs$row, jobs$bias)
todo <- which(!file.exists(jobs$file))
if (length(todo)) {
  message(sprintf("rerunning %d subset runs", length(todo)))
  invisible(parallel::mclapply(todo, function(j) {
    i <- jobs$row[j]; e <- exp.grid[i, ]; r <- Rparams.grid[[e$R]]
    set.seed(i)                               # same seed as the grid run of this row
    suppressMessages(run_joint_update_sites(e$site, Rparams = r$Rparams, useDistProduct = e$useDistProduct,
      tauPrior = e$tauPrior, Rparams.source = r$source, Rparams.model = r$model,
      outdir = dirname(jobs$file[j]), bias = if (jobs$bias[j]) grid.bias else NULL))
    NULL
  }, mc.cores = grid.cores))
}

## long table of every site-year in the subset: one-step-ahead quantities given z_t = 0
sy <- do.call(rbind, lapply(seq_len(nrow(jobs)), function(j) {
  h <- readRDS(jobs$file[j]); e <- exp.grid[jobs$row[j], ]
  data.frame(exp = which(exps$R == e$R & exps$useDistProduct == e$useDistProduct),
             bias = jobs$bias[j], site = e$site, t = seq_len(NT), y = h$obs,
             m0 = h$predMean0, V0 = h$predVar0, Q0 = h$predProcVar0, R = h$predR,
             V = h$predVar, pD = h$pDist, flag = flag[e$site, ], procVar = h$procVar)
}))
## undisturbed site-years: not product-flagged in t-1 or t, P(disturbed | y_{1:t}) < 0.5,
## and t >= 3 (year 1 starts from the population IC prior; year 2 still carries it)
sy$flagPrev <- ave(sy$flag, sy$exp, sy$bias, sy$site, FUN = function(f) c(FALSE, head(f, -1)))
und <- subset(sy, t >= 3 & !flag & !flagPrev & pD < 0.5 & is.finite(y))
und$e0  <- und$y - und$m0                    # innovation vs the no-disturbance forecast
und$P0  <- und$V0 - und$R - und$Q0           # propagated state (analysis) uncertainty
und$z   <- und$e0 / sqrt(und$V0)             # standardised innovation
und$pit0 <- pnorm(und$z)

## consistency: the subset reruns reproduce the grid runs
chk <- head(jobs[jobs$bias, ], 50)
cat("\nsubset rerun reproduces grid procVar (max abs diff over 50 runs):",
    max(sapply(seq_len(nrow(chk)), function(j) {
      e <- exp.grid[chk$row[j], ]
      max(abs(readRDS(chk$file[j])$procVar -
              readRDS(run_file("fullRunBias", e$site, e$R, e$useDistProduct, TRUE))$procVar))
    })), "\n")

budget <- do.call(rbind, lapply(split(und, list(und$exp, und$bias)), function(d) data.frame(
  exp = d$exp[1], bias = d$bias[1], n = nrow(d),
  E.e2  = mean(d$e0^2), med.e2 = median(d$e0^2),
  R = mean(d$R), Q = mean(d$Q0), P = mean(d$P0), V0 = mean(d$V0),
  dist.tail = mean(d$V - d$V0),
  var.z = mean(d$z^2), med.abs.z = median(abs(d$z)),
  in95 = 100 * mean(abs(d$z) <= qnorm(0.975)),
  ## scale c on R that would make E[e0^2] = c R + Q + P
  c.R = (mean(d$e0^2) - mean(d$Q0) - mean(d$P0)) / mean(d$R))))
budget <- budget[order(budget$exp, budget$bias), ]
budget$label <- exps$label[budget$exp]
cat("\n== 2. variance budget of the no-disturbance predictive, undisturbed site-years (subset) ==\n")
print(cbind(budget[, c("label", "bias", "n")], signif(budget[, c("E.e2", "med.e2", "R", "Q", "P", "V0",
      "dist.tail", "var.z", "med.abs.z", "in95", "c.R")], 3)), row.names = FALSE)

## figure: stacked mean variance budget vs mean / median squared innovation
png(file.path(out.dir, "fw_variance_budget.png"), width = 12, height = 5.6, units = "in", res = 150)
.par(mfrow = c(1, 2), mar = c(7.5, 4.8, 3.4, 0.8))
for (b in c(FALSE, TRUE)) {
  B <- budget[budget$bias == b, ]
  M <- rbind(R = B$R, `process 1/tau` = B$Q, `propagated state` = B$P)
  yl <- c(0, max(colSums(M), B$E.e2) * 1.45)
  x <- barplot(M, col = c("#eb6834", "#2a78d6", "#8fb6ea"), border = NA, ylim = yl, names.arg = rep("", ncol(M)),
               ylab = expression("variance  (kg C m"^-2*")"^2),
               main = if (b) "Bias-corrected" else "Uncorrected", space = 0.35)
  abline(h = pretty(yl), col = grid.col); barplot(M, col = c("#eb6834", "#2a78d6", "#8fb6ea"), border = NA,
         add = TRUE, names.arg = rep("", ncol(M)), space = 0.35, axes = FALSE)
  points(x, B$E.e2, pch = 23, bg = "white", col = ink, cex = 1.6, lwd = 1.5)
  points(x, B$med.e2 / qchisq(0.5, 1), pch = 21, bg = ink, col = "white", cex = 1.5)
  text(x, -yl[2] * 0.03, sub(", ", "\n", B$label), srt = 35, adj = c(1, 1), xpd = NA, cex = 0.75, col = ink)
  if (!b) legend("topleft", bty = "n", cex = 0.72, ncol = 2, text.col = ink,
                 legend = c("R (observation error)", "process variance E[1/tau]", "propagated state variance",
                            "mean squared innovation", "median sq. innovation / 0.455"),
                 fill = c("#eb6834", "#2a78d6", "#8fb6ea", NA, NA), border = NA,
                 pch = c(NA, NA, NA, 23, 21), pt.bg = c(NA, NA, NA, "white", ink), col = c(NA, NA, NA, ink, "white"))
}
dev.off()

## =================================================================================
## 3. model-free: year-to-year change of the observations vs 2R
## =================================================================================
## Under iid observation errors with variance R_t, Var(y_t - y_{t-1}) = Var(x_t - x_{t-1})
## + R_t + R_{t-1} >= R_t + R_{t-1}. Use every site-year pair not flagged by the product.
Y <- Yw[, paste0("t", atime)]
dY <- Y[, -1] - Y[, -NT]
okf <- !(flag[, -1] | flag[, -NT])
d  <- dY[okf & is.finite(dY)]
yy <- (Y[, -1])[okf & is.finite(dY)]; yp <- (Y[, -NT])[okf & is.finite(dY)]
robust.var <- function(v) (1.4826 * mad(v))^2
cat("\n== 3. year-to-year change of the observed biomass, unflagged pairs ==\n")
cat(sprintf("n = %d pairs; Var(dy) = %.4f, robust (MAD) Var(dy) = %.4f, median |dy| = %.3f\n",
            length(d), var(d), robust.var(d), median(abs(d))))
okl <- okf[, -1] & okf[, -(NT - 1)] & is.finite(dY[, -1]) & is.finite(dY[, -(NT - 1)])
cat(sprintf("lag-1 autocorrelation of dy (consecutive unflagged pairs): %.3f\n",
            cor(dY[, -1][okl], dY[, -(NT - 1)][okl])))
cat(sprintf("share of pairs with dy == 0: %.1f%%; obs resolution: %s\n", 100 * mean(d == 0),
            paste(head(sort(unique(round(abs(d), 3))), 4), collapse = ", ")))
dy.tab <- t(sapply(names(Rparams.grid), function(k) {
  rp <- Rparams.grid[[k]]$Rparams
  twoR <- (yy * rp[1] + rp[2]) + (yp * rp[1] + rp[2])
  c(mean.2R = mean(twoR), ratio = var(d) / mean(twoR), ratio.robust = robust.var(d) / mean(twoR),
    rho.min = 1 - robust.var(d) / mean(twoR),
    pct.within.1sd = 100 * mean(abs(d) <= sqrt(twoR)))
}))
print(signif(dy.tab, 3))

png(file.path(out.dir, "fw_obs_change.png"), width = 9, height = 5, units = "in", res = 150)
.par(mar = c(4.4, 4.6, 3.4, 1))
br <- seq(-3, 3, by = 0.05)
dd <- pmin(pmax(d, -3), 3)
hst <- hist(dd, breaks = br, plot = FALSE)
plot(NA, xlim = c(-3, 3), ylim = c(0, max(hst$density) * 1.05),
     xlab = expression(Delta*y[t] == y[t] - y[t-1]*"  (kg C m"^-2*", clipped at "*"±"*3*")"),
     ylab = "density", main = "Observed year-to-year biomass change (unflagged years) vs. iid-error prediction N(0, 2R)")
abline(h = pretty(c(0, max(hst$density))), col = grid.col)
rect(br[-length(br)], 0, br[-1], hst$density, col = "#8a8984", border = NA)
xs <- seq(-3, 3, length.out = 400)
cols <- c(GEDI_het = "#2a78d6", GEDI_con = "#1a4f99", LandTrendr_het = "#eb6834", LandTrendr_con = "#a8461f")
ltys <- c(GEDI_het = 1, GEDI_con = 2, LandTrendr_het = 1, LandTrendr_con = 2)
for (k in names(cols)) lines(xs, dnorm(xs, 0, sqrt(dy.tab[k, "mean.2R"])), col = cols[k], lwd = 2, lty = ltys[k])
legend("topright", bty = "n", cex = 0.85, text.col = ink,
       legend = c("observed", sprintf("N(0, 2R): %s", sub("_", " ", names(cols)))),
       fill = c("#8a8984", NA, NA, NA, NA), border = NA, col = c(NA, cols), lwd = c(NA, 2, 2, 2, 2),
       lty = c(NA, ltys))
dev.off()

## =================================================================================
## 4. standardised innovations of the no-disturbance predictive (bias-corrected)
## =================================================================================
png(file.path(out.dir, "fw_std_innov.png"), width = 14, height = 4, units = "in", res = 150)
.par(mfrow = c(1, 4), mar = c(4.2, 4.2, 3.4, 0.8))
zb <- seq(-4, 4, by = 0.25)
for (k in which(!exps$useDistProduct)) {
  z <- und$z[und$exp == k & und$bias]; z <- pmin(pmax(z, -4), 4)
  h <- hist(z, breaks = zb, plot = FALSE)
  plot(NA, xlim = c(-4, 4), ylim = c(0, 1.6), xlab = expression(z[t] == (y[t] - m[t]^0) / sqrt(V[t]^0)),
       ylab = if (k == 1) "density" else "", main = sprintf("%s\nvar(z) = %.2f (1 if calibrated)",
       sub(", no product", "", exps$label[k]), mean(und$z[und$exp == k & und$bias]^2)), cex.main = 1)
  abline(h = c(0.5, 1, 1.5), col = grid.col)
  rect(zb[-length(zb)], 0, zb[-1], h$density, col = c.cor, border = NA)
  lines(seq(-4, 4, length.out = 200), dnorm(seq(-4, 4, length.out = 200)), lwd = 2, lty = 2, col = ink)
}
dev.off()

## where does the process variance sit relative to the tau grid floor 1/tauMax?
cat(sprintf("\n== 5. tau grid: process variance range [%.0e, %.0e]; share of sites with v_T < 1e-3: ",
            1 / 1e4, 1 / 0.1))
cat(paste(sprintf("%s %.0f%%/%.0f%%", exps$label, sapply(PV, function(p) 100 * mean(p$raw[, NT] < 1e-3)),
                  sapply(PV, function(p) 100 * mean(p$cor[, NT] < 1e-3))), collapse = "; "), "(uncorrected/corrected)\n")

## =================================================================================
## 6. PIT over time: does it flatten as the learned process variance shrinks?
## =================================================================================
## Full grids, no disturbance product (no flagged years), all 538 sites. The forecast
## of y_t uses the tau posterior after y_{t-1}, so year t is paired with v_{t-1}.
## Per-year quantities, from each run's stored one-step-ahead predictive:
##   PIT u_t (full predictive), Q_t = predProcVar (process part), R_t = predR,
##   V_t = predVar (total, incl. the disturbance tail).
pit_tab <- function(dir, R, bias) do.call(rbind, lapply(seq_len(NS), function(s) {
  h <- readRDS(run_file(dir, s, R, FALSE, bias))
  data.frame(site = s, t = seq_len(NT), pit = h$pit, vPrev = c(NA, head(h$procVar, -1)),
             Q = h$predProcVar, R = h$predR, V = h$predVar)
}))
## coverage error of the central 50% interval (> 0: too wide) and mean |coverage
## error| over the central p-intervals, p = 0.05, ..., 0.95 (0 = calibrated)
cov50  <- function(u) { u <- u[is.finite(u)]; mean(abs(2 * u - 1) <= 0.5) - 0.5 }
calerr <- function(u) { p <- seq(0.05, 0.95, 0.05); u <- u[is.finite(u)]
  mean(abs(sapply(p, function(pp) mean(abs(2 * u - 1) <= pp)) - p)) }
periods <- list("1991-95" = 1991:1995, "1996-2000" = 1996:2000, "2001-05" = 2001:2005,
                "2006-11" = 2006:2011, "2012-17" = 2012:2017)
Rm.np <- exps$R[!exps$useDistProduct]
PT <- lapply(setNames(Rm.np, Rm.np), function(R) list(
  raw = subset(pit_tab("fullRunNoBias", R, FALSE), t >= 2 & is.finite(pit)),
  cor = subset(pit_tab("fullRunBias",   R, TRUE),  t >= 2 & is.finite(pit))))
Rlab <- c(GEDI_het = "GEDI heteroskedastic R", GEDI_con = "GEDI constant R",
          LandTrendr_het = "LandTrendr heteroskedastic R", LandTrendr_con = "LandTrendr constant R")

## per-year summary
yr.tab <- do.call(rbind, lapply(Rm.np, function(R) do.call(rbind, lapply(c("raw", "cor"), function(k) {
  d <- PT[[R]][[k]]
  do.call(rbind, lapply(split(d, d$t), function(x) data.frame(
    R.mod = R, run = k, year = atime[x$t[1]], v = median(x$vPrev), cov50 = cov50(x$pit),
    calErr = calerr(x$pit), Qshare = median(x$Q / x$V), Rshare = median(x$R / x$V),
    vMean = mean(x$vPrev), vQ25 = quantile(x$vPrev, 0.25, names = FALSE), vQ75 = quantile(x$vPrev, 0.75, names = FALSE),
    pitMean = mean(x$pit), pitMed = median(x$pit),
    pitQ25 = quantile(x$pit, 0.25, names = FALSE), pitQ75 = quantile(x$pit, 0.75, names = FALSE))))
}))))
cat("\n== 6. PIT over time (no product, 538 sites): per period, median v_{t-1}, coverage error of the central 50% interval, calibration error ==\n")
per.tab <- do.call(rbind, lapply(Rm.np, function(R) do.call(rbind, lapply(c("raw", "cor"), function(k) {
  d <- PT[[R]][[k]]; yrs <- atime[d$t]
  do.call(rbind, lapply(names(periods), function(p) { x <- d[yrs %in% periods[[p]], ]
    data.frame(R.mod = R, run = k, period = p, v = median(x$vPrev), cov50 = cov50(x$pit), calErr = calerr(x$pit),
               Qshare = median(x$Q / x$V), Rshare = median(x$R / x$V)) }))
}))))
print(cbind(per.tab[, 1:3], signif(per.tab[, -(1:3)], 3)), row.names = FALSE)
cat("\nSpearman correlation over years of median v_{t-1} with cov50 (corrected / uncorrected):\n")
for (R in Rm.np) { a <- yr.tab[yr.tab$R.mod == R & yr.tab$run == "cor", ]; b <- yr.tab[yr.tab$R.mod == R & yr.tab$run == "raw", ]
  cat(sprintf("  %-28s %5.2f / %5.2f\n", Rlab[R], cor(a$v, a$cov50, method = "spearman"), cor(b$v, b$cov50, method = "spearman"))) }

## figure (a): PIT histograms by period, rows = observation-error model
br <- seq(0, 1, 0.05); nb <- length(br) - 1
png(file.path(out.dir, "fw_pit_periods.png"), width = 15, height = 11, units = "in", res = 140)
.par(mfrow = c(length(Rm.np), length(periods)), mar = c(3.2, 3.4, 3.4, 0.6), oma = c(1.5, 2, 3, 0))
for (R in Rm.np) for (p in names(periods)) {
  dens <- function(d) { u <- d$pit[atime[d$t] %in% periods[[p]]]; hist(u, br, plot = FALSE)$counts * nb / length(u) }
  d1 <- dens(PT[[R]]$cor); d0 <- dens(PT[[R]]$raw)
  n  <- sum(atime[PT[[R]]$cor$t] %in% periods[[p]])
  band <- qbinom(c(0.025, 0.975), n, 1 / nb) * nb / n
  k1 <- per.tab[per.tab$R.mod == R & per.tab$run == "cor" & per.tab$period == p, ]
  k0 <- per.tab[per.tab$R.mod == R & per.tab$run == "raw" & per.tab$period == p, ]
  plot(NA, xlim = c(0, 1), ylim = c(0, 4), xlab = "", ylab = "", cex.main = 0.95,
       main = sprintf("%s\nmedian v: %.3f corr. / %.3f uncorr.", p, k1$v, k0$v))
  rect(0, band[1], 1, band[2], col = grid.col, border = NA)
  rect(br[-length(br)] + 0.004, 0, br[-1] - 0.004, pmin(d1, 4), col = c.cor, border = NA)
  lines(rep(br, each = 2), c(0, rep(pmin(d0, 4), each = 2), 0), col = ink, lwd = 1.5)
  abline(h = 1, lty = 2, col = ink)
  if (p == names(periods)[1]) mtext(Rlab[R], side = 2, line = 2.6, las = 0, col = ink, cex = 0.85)
}
mtext("PIT u", side = 1, outer = TRUE, line = 0.2, col = ink)
mtext("One-step-ahead PIT by period, no product, 538 sites:  bars = bias-corrected, outline = uncorrected   [dashed = uniform; gray = 95% range per bin]",
      outer = TRUE, line = 1, cex = 1.05, col = ink)
dev.off()

## figure (b): per year t (x axis), one row per observation-error model, bias-corrected
## (blue) vs uncorrected (gray):
##   left:   v_{t-1}, the process variance the forecast of y_t uses -- median and IQR
##           over sites (line, band) and mean over sites (thin line)
##   middle: PIT of y_t over sites -- median and IQR (line, band) and mean (thin line);
##           calibrated: median 0.5, IQR 0.25-0.75
##   right:  coverage error of the central 50% interval (> 0: too wide; 0: calibrated)
cat("\n== 6b. per-year process variance and PIT (no product, 538 sites) ==\n")
print(cbind(yr.tab[, 1:3], signif(yr.tab[, c("v", "vMean", "vQ25", "vQ75", "pitMean", "pitMed", "pitQ25", "pitQ75", "cov50")], 3)),
      row.names = FALSE)
write.csv(yr.tab, file.path(out.dir, "fw_pit_time.csv"), row.names = FALSE)
png(file.path(out.dir, "fw_pit_time.png"), width = 14, height = 13, units = "in", res = 140)
.par(mfrow = c(length(Rm.np), 3), mar = c(3.6, 4.6, 3.2, 0.8), oma = c(1.2, 2, 3.4, 0))
vyl <- range(yr.tab[, c("vQ25", "vQ75", "vMean")])
for (R in Rm.np) {
  X <- list(raw = yr.tab[yr.tab$R.mod == R & yr.tab$run == "raw", ],
            cor = yr.tab[yr.tab$R.mod == R & yr.tab$run == "cor", ])
  sty <- list(raw = list(col = c.raw, lty = 2), cor = list(col = c.cor, lty = 1))
  band_line <- function(x, lo, mid, hi, mn, k) {
    polygon(c(x$year, rev(x$year)), c(x[[lo]], rev(x[[hi]])), col = adjustcolor(sty[[k]]$col, 0.18), border = NA)
    lines(x$year, x[[mid]], col = sty[[k]]$col, lwd = 2.2, lty = sty[[k]]$lty)
    lines(x$year, x[[mn]], col = sty[[k]]$col, lwd = 1, lty = 3)
  }
  ## process variance
  plot(NA, xlim = range(atime[-1]), ylim = vyl, log = "y", xlab = "", ylab = "", cex.main = 1,
       main = expression("process variance "*v[t-1]*"  (kg C m"^-2*")"^2))
  abline(h = 10^(-3:0) * rep(c(1, 2, 5), each = 4), col = grid.col)
  for (k in c("raw", "cor")) band_line(X[[k]], "vQ25", "v", "vQ75", "vMean", k)
  mtext(Rlab[R], side = 2, line = 3.6, las = 0, col = ink, cex = 0.9)
  if (R == Rm.np[1]) legend("topright", bty = "n", cex = 0.85, text.col = ink,
    legend = c("bias-corrected: median, IQR", "uncorrected: median, IQR", "mean over sites"),
    col = c(c.cor, c.raw, ink), lwd = c(2.2, 2.2, 1), lty = c(1, 2, 3))
  ## PIT
  plot(NA, xlim = range(atime[-1]), ylim = c(0, 1), xlab = "", ylab = "", cex.main = 1,
       main = "PIT of the one-step-ahead forecast")
  rect(min(atime) - 1, 0.25, max(atime) + 1, 0.75, col = adjustcolor(ink, 0.06), border = NA)
  abline(h = c(0.25, 0.5, 0.75), col = ink, lty = c(3, 2, 3))
  for (k in c("raw", "cor")) band_line(X[[k]], "pitQ25", "pitMed", "pitQ75", "pitMean", k)
  if (R == Rm.np[1]) legend("bottomright", bty = "n", cex = 0.8, text.col = ink,
    legend = c("calibrated median (0.5) and IQR (0.25-0.75)"), fill = adjustcolor(ink, 0.06), border = NA)
  ## coverage error
  plot(NA, xlim = range(atime[-1]), ylim = range(yr.tab$cov50, 0), xlab = "", ylab = "", cex.main = 1,
       main = "coverage error of central 50% interval (> 0 too wide)")
  abline(h = axTicks(2), col = grid.col); abline(h = 0, lty = 2, col = ink)
  for (k in c("raw", "cor")) lines(X[[k]]$year, X[[k]]$cov50, col = sty[[k]]$col, lwd = 2.2, lty = sty[[k]]$lty)
}
mtext("year", side = 1, outer = TRUE, line = 0, col = ink)
mtext("Process variance and PIT over time, no disturbance product, 538 sites", outer = TRUE, line = 1.2, cex = 1.15, col = ink)
dev.off()

## =================================================================================
## 7. decomposition of the one-step-ahead predictive variance over time
## =================================================================================
## Needs the no-disturbance predictive (predVar0 / predProcVar0), which only the
## 100-site subset reruns of section 2 store, so this uses that subset (no product).
## Exact decomposition of the full predictive variance V_t of y_t (see
## claude-forecast-predictive-10.9.26 for the derivation):
##   V_t = R_t + beta0 Q0_t + beta0 P0_t + D_t
##   R_t       observation-error variance                           (predR)
##   beta0     P(z_t = 0 | y_{1:t-1}) = predProcVar / predProcVar0
##   Q0_t      process variance of the no-disturbance predictive    (predProcVar0)
##   P0_t      propagated state variance = predVar0 - R_t - Q0_t
##   D_t       disturbance part = V_t - R_t - beta0 (predVar0 - R_t)
dec <- do.call(rbind, lapply(which(!exp.grid$useDistProduct[jobs$row]), function(j) {
  h <- readRDS(jobs$file[j]); e <- exp.grid[jobs$row[j], ]
  b0 <- h$predProcVar / h$predProcVar0
  data.frame(R.mod = e$R, run = if (jobs$bias[j]) "cor" else "raw", site = e$site, t = seq_len(NT),
             V = h$predVar, R = h$predR, beta0 = b0, Q = b0 * h$predProcVar0,
             P = b0 * (h$predVar0 - h$predR - h$predProcVar0),
             D = h$predVar - h$predR - b0 * (h$predVar0 - h$predR))
}))
dec <- subset(dec, t >= 2 & is.finite(V))
stopifnot(max(abs(with(dec, R + Q + P + D - V))) < 1e-8)       # the four parts add up to V
comp <- c(R = "observation error R", Q = "process error", P = "propagated state", D = "disturbance")
comp.col <- c(R = "#eb6834", Q = "#2a78d6", P = "#8fb6ea", D = "#8a8984")
dec.yr <- do.call(rbind, lapply(split(dec, list(dec$R.mod, dec$run, dec$t), drop = TRUE), function(x) {
  sh <- sapply(names(comp), function(k) median(x[[k]] / x$V))
  data.frame(R.mod = x$R.mod[1], run = x$run[1], year = atime[x$t[1]], n = nrow(x),
             V = mean(x$V), R = mean(x$R), Q = mean(x$Q), P = mean(x$P), D = mean(x$D),
             medV = median(x$V), beta0 = mean(x$beta0),
             shR = sh[["R"]], shQ = sh[["Q"]], shP = sh[["P"]], shD = sh[["D"]])
}))
dec.yr <- dec.yr[order(match(dec.yr$R.mod, Rm.np), dec.yr$run, dec.yr$year), ]
write.csv(dec.yr, file.path(out.dir, "fw_var_decomp.csv"), row.names = FALSE)
cat("\n== 7. predictive-variance decomposition (100-site subset, no product): mean component per year ==\n")
print(cbind(dec.yr[dec.yr$year %in% c(1991, 1995, 2000, 2005, 2010, 2017), 1:3],
            signif(dec.yr[dec.yr$year %in% c(1991, 1995, 2000, 2005, 2010, 2017), -(1:4)], 3)), row.names = FALSE)

## figure (a): stacked mean components per year; rows = model, columns = uncorrected / corrected
png(file.path(out.dir, "fw_var_decomp.png"), width = 12, height = 13, units = "in", res = 140)
.par(mfrow = c(length(Rm.np), 2), mar = c(3.4, 4.6, 3.2, 0.8), oma = c(1.2, 2, 3.4, 0))
for (R in Rm.np) {
  ymax <- max(dec.yr$V[dec.yr$R.mod == R]) * 1.05
  for (k in c("raw", "cor")) {
    x <- dec.yr[dec.yr$R.mod == R & dec.yr$run == k, ]
    plot(NA, xlim = range(x$year), ylim = c(0, ymax), xlab = "", ylab = "", cex.main = 1,
         main = sprintf("%s, %s", Rlab[R], if (k == "cor") "bias-corrected" else "uncorrected"))
    abline(h = axTicks(2), col = grid.col)
    lo <- rep(0, nrow(x))
    for (cm in names(comp)) {
      hi <- lo + x[[cm]]
      polygon(c(x$year, rev(x$year)), c(lo, rev(hi)), col = comp.col[cm], border = "white", lwd = 0.5)
      lo <- hi
    }
    if (k == "raw") mtext(expression("mean variance  (kg C m"^-2*")"^2), side = 2, line = 3, las = 0, col = ink, cex = 0.75)
    if (R == Rm.np[1] && k == "raw") legend("topright", bty = "n", cex = 0.85, text.col = ink,
      legend = comp, fill = comp.col, border = NA)
  }
}
mtext("year", side = 1, outer = TRUE, line = 0, col = ink)
mtext("One-step-ahead predictive variance V = R + process + propagated state + disturbance (mean over 100 sites, no product)",
      outer = TRUE, line = 1.2, cex = 1.1, col = ink)
dev.off()

## figure (b): median per-site share of V from each component
png(file.path(out.dir, "fw_var_share.png"), width = 16, height = 4.4, units = "in", res = 140)
.par(mfrow = c(1, length(Rm.np)), mar = c(4.2, 4.4, 3.4, 0.8))
for (R in Rm.np) {
  plot(NA, xlim = range(atime[-1]), ylim = c(0, 1), xlab = "year", ylab = if (R == Rm.np[1]) "median share of V over sites" else "",
       main = Rlab[R], cex.main = 1)
  abline(h = seq(0, 1, 0.2), col = grid.col)
  for (k in c("raw", "cor")) { x <- dec.yr[dec.yr$R.mod == R & dec.yr$run == k, ]
    for (cm in names(comp)) lines(x$year, x[[paste0("sh", cm)]], col = comp.col[cm], lwd = if (k == "cor") 2.2 else 1.3,
                                  lty = if (k == "cor") 1 else 3) }
  if (R == Rm.np[1]) legend("topleft", bty = "n", cex = 0.8, text.col = ink, ncol = 2,
    legend = c(comp, "bias-corrected", "uncorrected"), col = c(comp.col, ink, ink),
    lwd = c(rep(2.2, 4), 2.2, 1.3), lty = c(rep(1, 4), 1, 3))
}
dev.off()

## check of the undisturbed-year PIT approximation (claude-forecast-predictive, Sec. 3.5):
##   u_t ~ (1 - beta0) + beta0 * Phi((y_t - m0) / sqrt(V0)), against the stored exact PIT,
## bias-corrected subset, no product, undisturbed site-years (t >= 3, P(disturbed) < 0.5)
ua <- do.call(rbind, lapply(which(!exp.grid$useDistProduct[jobs$row] & jobs$bias), function(j) {
  h <- readRDS(jobs$file[j]); b0 <- h$predProcVar / h$predProcVar0
  data.frame(t = seq_len(NT), u = h$pit, pD = h$pDist,
             ua = (1 - b0) + b0 * pnorm((h$obs - h$predMean0) / sqrt(h$predVar0)))
}))
ua <- subset(ua, t >= 3 & is.finite(u) & pD < 0.5)
d <- ua$u - ua$ua
cat(sprintf("\n== 7b. undisturbed-year PIT approximation: n = %d, mean |u - u~| = %.4f, 95%% |.| = %.4f, max = %.3f, mean(u - u~) = %+.4f, cor = %.4f, IQR u = %.3f, IQR u~ = %.3f ==\n",
            nrow(ua), mean(abs(d)), quantile(abs(d), 0.95), max(abs(d)), mean(d), cor(ua$u, ua$ua), IQR(ua$u), IQR(ua$ua)))
