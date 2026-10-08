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
