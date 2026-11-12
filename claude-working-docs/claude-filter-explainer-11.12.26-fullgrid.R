## Full-grid figures and numbers for claude-filter-explainer-11.12.26.tex.
##
## Reads the full OR grid run with the default settings (run-model-or.R: all 538 sites x 4
## observation-error models x with / without the disturbance product; bias correction on,
## process error propagated through the disturbance), and, for comparison, the earlier grid
## with bias correction but without the disturbed-pool process error
## (experiments/old/fullRunBias).
##
## Run from the repository root, after  Rscript or-model/run-model-or.R :
##   Rscript claude-working-docs/claude-filter-explainer-11.12.26-fullgrid.R [--new DIR] [--old DIR]

source("or-model/R/config.R")
out.dir <- "claude-working-docs/claude-filter-explainer-11.12.26-files"
new.dir <- cli_opt("--new", "experiments/or-model/runs")
old.dir <- cli_opt("--old", Filter(dir.exists, c("experiments/old/fullRunBias", "../../../experiments/old/fullRunBias"))[1])
load(or.files$extract)
colnames(dextract) <- paste0("t", dtime); dextract[is.na(dextract)] <- 0
Yw <- bextract * B2C * 1000 * 1e-4; colnames(Yw) <- paste0("t", ctime)
NS <- nrow(dextract)
Rmods <- c("GEDI_het", "GEDI_con", "LandTrendr_het", "LandTrendr_con")
Rlab  <- c(GEDI_het = "GEDI, heterosk.", GEDI_con = "GEDI, constant", LandTrendr_het = "LandTrendr, heterosk.",
           LandTrendr_con = "LandTrendr, constant")
src <- function(R) sub("_.*", "", R); mdl <- function(R) sub(".*_", "", R)
fname <- function(dir, s, R, p) file.path(dir, sprintf("site_%d_distprod_%s_%s_%s_taugamma_bias.rds", s, p, src(R), mdl(R)))

## the product's flags in the DA window, and the observed one-year wood losses
yrs   <- intersect(dtime, atime)
flagY <- dextract[, paste0("t", yrs)] > 0                         # NS x years (1990-2012)
lossY <- sapply(yrs, function(y) if (y == min(ctime)) rep(NA_real_, NS) else Yw[, paste0("t", y - 1)] - Yw[, paste0("t", y)])
colnames(lossY) <- yrs                                              # wood loss from year-1 to year (kg C m-2)

## per run: summaries
summ <- function(h, s) {
  T <- length(h$year); tT <- max(which(!is.na(h$obs)))
  iy <- match(yrs, h$year)
  list(pit = h$pit[-1],                                              # skip year 1 (IC prior)
       div = (h$obs[tT] - h$mean[tT]) / ((h$ci95[tT, "hi"] - h$ci95[tT, "lo"]) / 2),
       pD  = h$pDist[iy], sd = sqrt(h$procVar[T]))
}
read_grid <- function(dir) {
  out <- list()
  for (R in Rmods) for (p in c(TRUE, FALSE)) {
    f <- fname(dir, seq_len(NS), R, p)
    ok <- file.exists(f)
    message(sprintf("  %s: %s, product %s: %d / %d runs", basename(dir), R, p, sum(ok), NS))
    out[[paste(R, p)]] <- list(R = R, p = p, sites = which(ok), s = lapply(which(ok), function(s) summ(readRDS(f[s]), s)))
  }
  out
}
cache <- file.path("experiments/or-model", "fullgrid-summaries.rds")
if (file.exists(cache) && file.mtime(cache) > max(file.mtime(dir(new.dir, full.names = TRUE)))) {
  G <- readRDS(cache)
} else {
  G <- list(new = read_grid(new.dir), old = if (!is.na(old.dir)) read_grid(old.dir))
  saveRDS(G, cache)
}

## metrics per setting
metrics <- function(e) {
  pit <- unlist(lapply(e$s, `[[`, "pit")); div <- sapply(e$s, `[[`, "div")
  pD  <- do.call(rbind, lapply(e$s, `[[`, "pD"))                    # sites x years
  fl  <- flagY[e$sites, ]; ls <- lossY[e$sites, ]
  big   <- !is.na(ls) & ls > 3                                       # wood loss > 3 kg C m-2 in one year
  quiet <- !is.na(ls) & abs(ls) < 0.5 & !fl                          # small change, not flagged
  c(n = length(e$s),
    pit95 = mean(pit > 0.025 & pit < 0.975), pit50 = mean(pit > 0.25 & pit < 0.75),
    div1 = mean(abs(div) <= 1), divMed = median(abs(div)),
    det = mean(pD[big] > 0.5), nBig = sum(big), fa = mean(pD[quiet] > 0.5), nQuiet = sum(quiet),
    detFlag = mean(pD[fl] > 0.5), sd = median(sapply(e$s, `[[`, "sd")))
}
M <- lapply(G[!sapply(G, is.null)], function(g) t(sapply(g, metrics)))
print(lapply(M, round, 3))
write.csv(do.call(rbind, lapply(names(M), function(n) data.frame(grid = n, setting = rownames(M[[n]]), M[[n]]))),
          file.path(out.dir, "fullgrid-metrics.csv"), row.names = FALSE)

## LaTeX table rows: one per setting (new grid)
pc <- function(x) sprintf("%.0f\\%%", 100 * x)
m <- M$new
rows <- sprintf("%s & %s & %s & %s & %s & %s & %s & %s & %.3f \\\\", Rlab[sub(" .*", "", rownames(m))],
                ifelse(grepl("TRUE", rownames(m)), "yes", "no"), pc(m[, "pit95"]), pc(m[, "pit50"]), pc(m[, "div1"]),
                pc(m[, "det"]), pc(m[, "detFlag"]), pc(m[, "fa"]), m[, "sd"])
writeLines(rows, file.path(out.dir, "rows-fullgrid.tex"))
if (!is.null(M$old)) {
  o <- M$old
  one <- function(lab, x) sprintf("%s & %s & %s & %s & %s & %s & %s & %.3f \\\\", lab, pc(mean(x[, "pit95"])), pc(mean(x[, "pit50"])),
                                   pc(mean(x[, "div1"])), pc(mean(x[, "det"])), pc(mean(x[, "detFlag"])), pc(mean(x[, "fa"])), median(x[, "sd"]))
  cmp <- c(one("propagated (this grid)", m), one("none (earlier grid)", o))
  writeLines(cmp, file.path(out.dir, "rows-fullgrid-compare.tex"))
}
writeLines(c(sprintf("\\newcommand{\\FGNBig}{%d}", m[1, "nBig"]), sprintf("\\newcommand{\\FGNQuiet}{%d}", m[1, "nQuiet"]),
             sprintf("\\newcommand{\\FGRuns}{%d}", sum(m[, "n"]))),
           file.path(out.dir, "numbers-fullgrid.tex"))

## ---- figures -------------------------------------------------------------------------
ink <- "#3d3d3a"; grid.col <- "#e6e6e3"; c.p <- "#eb6834"; c.n <- "#8a8984"; c.dist <- "#c0392b"
base_par <- function(...) par(las = 1, bty = "n", col.axis = ink, fg = ink, col.lab = ink, col.main = ink, mgp = c(2.4, 0.7, 0), ...)

## calibration: PIT histograms per observation-error model, product on / off
pdf(file.path(out.dir, "fullgrid-pit.pdf"), width = 10, height = 4.4)
base_par(mfrow = c(1, 4), mar = c(4, 3.6, 2.4, 0.6), oma = c(0, 1, 0, 0))
br <- seq(0, 1, 0.1)
D <- lapply(Rmods, function(R) lapply(c(TRUE, FALSE), function(p)
  hist(unlist(lapply(G$new[[paste(R, p)]]$s, `[[`, "pit")), breaks = br, plot = FALSE)$density))
names(D) <- Rmods
for (R in Rmods) {
  d <- D[[R]]
  plot(NA, xlim = c(0, 1), ylim = c(0, max(unlist(D)) * 1.05), xlab = "PIT", ylab = if (R == Rmods[1]) "density" else "", main = Rlab[R])
  abline(h = 1, lty = 2, col = ink)
  rect(br[-11], 0, br[-1] - 0.05, d[[1]], col = c.p, border = NA); rect(br[-11] + 0.05, 0, br[-1], d[[2]], col = c.n, border = NA)
  if (R == Rmods[1]) legend("topleft", bty = "n", cex = 0.85, legend = c("with product", "without"), fill = c(c.p, c.n), border = NA)
}
dev.off()

## detection: P(disturbed) per site-year against the observed one-year wood loss (GEDI heteroskedastic)
pdf(file.path(out.dir, "fullgrid-detection.pdf"), width = 10, height = 4.2)
base_par(mfrow = c(1, 2), mar = c(4, 4.2, 2.4, 1))
for (p in c(FALSE, TRUE)) {
  e  <- G$new[[paste("GEDI_het", p)]]
  pD <- do.call(rbind, lapply(e$s, `[[`, "pD")); fl <- flagY[e$sites, ]; ls <- lossY[e$sites, ]
  ok <- !is.na(ls)
  plot(NA, xlim = range(ls[ok]), ylim = c(0, 1), xlab = "observed one-year wood loss (kg C m-2)", ylab = "P(disturbed | y)",
       main = sprintf("GEDI heteroskedastic R, %s product", if (p) "with" else "without"))
  abline(h = c(0, 0.5, 1), col = grid.col); abline(v = 0, col = grid.col)
  points(ls[ok & !fl], pD[ok & !fl], pch = 19, cex = 0.35, col = adjustcolor(ink, 0.25))
  points(ls[ok & fl], pD[ok & fl], pch = 19, cex = 0.5, col = adjustcolor(c.dist, 0.7))
  legend("right", bty = "n", cex = 0.85, legend = c("product flags the year", "not flagged"), pch = 19, col = c(c.dist, ink))
}
dev.off()
message("full-grid figures and numbers in ", out.dir)
