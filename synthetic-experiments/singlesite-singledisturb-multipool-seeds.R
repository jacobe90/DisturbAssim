## Seed sweep of the single site, single disturbance, multi-pool experiment.
##
## Runs singlesite-singledisturb-multipool.R unchanged for several random seeds (each seed
## gives new drivers, parameter ensemble, truth and observations; the true disturbance is
## always in year 10) and compares the posterior probability of disturbance across them:
## detection in year 10, and false alarms in the other years. Seed 42 is the original run.
##
## Run from the repository root:  Rscript synthetic-experiments/singlesite-singledisturb-multipool-seeds.R
## Outputs (in sweep.dir):
##   summary.txt               headline numbers
##   summary_by_seed.csv       per seed: P(disturbed) in year 10, false-alarm counts, max elsewhere
##   false_alarms.csv          every year outside year 10 with P(disturbed) > 0.2, with context
##   pdist_by_seed.png         heatmap of P(disturbed) by seed and year
##   pdist_summary.pdf         heatmap, P(disturbed) by year across seeds, tau learned vs known
##   sweep_results.rds         per-seed P(disturbed), states, observations and quantiles
##   seed_<s>/figures.pdf      the base experiment's figures for each seed

seeds     <- c(42, 1:20)
sweep.dir <- "experiments/synthetic/singlesite-singledisturb-multipool-seeds"
base      <- "synthetic-experiments/singlesite-singledisturb-multipool.R"
n.cores   <- min(length(seeds), max(1, parallel::detectCores() - 2))
fa.level  <- 0.2   # P(disturbed) counted as a false alarm outside the disturbance year
if (!dir.exists(sweep.dir)) dir.create(sweep.dir, recursive = TRUE)

## ---- Run every seed --------------------------------------------------------------------
## Each seed sources the base experiment in its own environment, with its seed and output
## folder set first. The base run's results.rds (tens of MB, mostly stored mixtures) is
## replaced by the compact pieces the comparison needs.
run_seed <- function(s) {
  e <- new.env()
  e$seed   <- s
  e$outdir <- file.path(sweep.dir, sprintf("seed_%02d", s))
  log <- capture.output(sys.source(base, envir = e))
  unlink(file.path(e$outdir, "results.rds"))
  list(seed = s, log = log, p.dist = e$p.dist, R = e$R, Q = e$Q, disturbance_at = e$disturbance_at,
       truth = e$True_states, obs = e$Observed_states,
       pDist = e$or.fit$pDist, pDist.known = e$or.known$pDist,
       q = e$or.fit$q, fq = e$or.fit$fq, q.known = e$or.known$q, procVar = e$or.fit$procVar,
       an.stats = e$an.stats, an0.stats = e$an0.stats)
}
t0  <- Sys.time()
res <- parallel::mclapply(seeds, function(s) tryCatch(run_seed(s), error = function(err)
         list(seed = s, error = conditionMessage(err))), mc.cores = n.cores, mc.preschedule = FALSE)
failed <- sapply(res, function(x) !is.null(x$error))
if (any(failed)) for (x in res[failed]) message(sprintf("seed %d failed: %s", x$seed, x$error))
res   <- res[!failed]
seeds <- sapply(res, `[[`, "seed")
message(sprintf("%d seeds finished in %.1f min", length(res), as.numeric(difftime(Sys.time(), t0, units = "mins"))))

## ---- Compare -----------------------------------------------------------------------------
td <- which(res[[1]]$disturbance_at == 1)
NT <- length(res[[1]]$pDist)
R  <- res[[1]]$R
P  <- t(sapply(res, `[[`, "pDist"))           # seed x year
PK <- t(sapply(res, `[[`, "pDist.known"))
rownames(P) <- rownames(PK) <- seeds
clip <- t(sapply(res, function(x) apply(x$obs == 0, 1, any)))   # an observation clipped to 0

## standardised one-step-ahead innovation per pool, (y - forecast median) / sqrt(forecast var + R),
## with the forecast variance from its 95% interval
z <- lapply(res, function(x) {
  sdF <- (x$fq[5, , ] - x$fq[1, , ]) / 3.92
  (x$obs - x$fq[3, , ]) / sqrt(sdF^2 + rep(diag(R), each = NT))
})

other <- setdiff(seq_len(NT), td)
summ <- data.frame(
  seed = seeds, p.year10 = P[, td], p.year10.tauknown = PK[, td],
  n.false.gt0.2 = rowSums(P[, other] > 0.2), n.false.gt0.5 = rowSums(P[, other] > 0.5),
  max.other = apply(P[, other], 1, max),
  year.of.max = other[apply(P[, other], 1, which.max)],
  max.other.tauknown = apply(PK[, other], 1, max))
write.csv(summ, file.path(sweep.dir, "summary_by_seed.csv"), row.names = FALSE)

fa <- which(P > fa.level & col(P) != td, arr.ind = TRUE)
fa <- fa[order(fa[, 1], fa[, 2]), , drop = FALSE]
fat <- data.frame(
  seed = seeds[fa[, 1]], year = fa[, 2], years.since.dist = fa[, 2] - td,
  p = P[fa], p.tauknown = PK[fa], obs.clipped = clip[fa],
  obs.leaf = sapply(seq_len(nrow(fa)), function(i) res[[fa[i, 1]]]$obs[fa[i, 2], 1]),
  obs.stem = sapply(seq_len(nrow(fa)), function(i) res[[fa[i, 1]]]$obs[fa[i, 2], 3]),
  true.leaf = sapply(seq_len(nrow(fa)), function(i) res[[fa[i, 1]]]$truth[fa[i, 2], 1]),
  true.stem = sapply(seq_len(nrow(fa)), function(i) res[[fa[i, 1]]]$truth[fa[i, 2], 3]),
  z.leaf = sapply(seq_len(nrow(fa)), function(i) z[[fa[i, 1]]][fa[i, 2], 1]),
  z.stem = sapply(seq_len(nrow(fa)), function(i) z[[fa[i, 1]]][fa[i, 2], 3]))
write.csv(fat, file.path(sweep.dir, "false_alarms.csv"), row.names = FALSE)

und  <- col(P) != td
pre  <- col(P) < td
post <- col(P) > td
lines.out <- c(
  sprintf("Seed sweep of %s: %d seeds (%s), prior p.dist = %.3f, true disturbance in year %d",
          base, length(seeds), paste(seeds, collapse = ", "), res[[1]]$p.dist, td),
  sprintf("year %d: P(disturbed) median %.3f, min %.3f (seed %d); tau known: median %.3f, min %.3f",
          td, median(P[, td]), min(P[, td]), seeds[which.min(P[, td])], median(PK[, td]), min(PK[, td])),
  sprintf("other years: %d of %d seed-years with P > 0.2, %d with P > 0.5 (tau known: %d, %d)",
          sum(P[und] > 0.2), sum(und), sum(P[und] > 0.5), sum(PK[und] > 0.2), sum(PK[und] > 0.5)),
  sprintf("  before year %d: mean P %.4f, P > 0.2 in %d of %d", td, mean(P[pre]), sum(P[pre] > 0.2), sum(pre)),
  sprintf("  after  year %d: mean P %.4f, P > 0.2 in %d of %d", td, mean(P[post]), sum(P[post] > 0.2), sum(post)),
  sprintf("  seeds with any year P > 0.2: %d of %d; largest false alarm %.3f (seed %d, year %d)",
          sum(summ$n.false.gt0.2 > 0), length(seeds), max(summ$max.other),
          summ$seed[which.max(summ$max.other)], summ$year.of.max[which.max(summ$max.other)]),
  sprintf("  years with an observation clipped to 0: %d, of which P > 0.2: %d; unclipped years with P > 0.2: %d of %d",
          sum(clip & und), sum(P[clip & und] > 0.2), sum(P[!clip & und] > 0.2), sum(!clip & und)),
  sprintf("seed 42 (original): P(year %d) = %.3f, largest elsewhere %.3f (year %d)", td,
          P["42", td], summ$max.other[summ$seed == 42], summ$year.of.max[summ$seed == 42]))
writeLines(lines.out, file.path(sweep.dir, "summary.txt"))
cat(lines.out, sep = "\n")

## ---- Figures -----------------------------------------------------------------------------
heatmap_pdist <- function() {
  par(mar = c(4.2, 4.2, 3, 6.5))
  cols <- hcl.colors(50, "YlOrRd", rev = TRUE)
  ny <- length(seeds)
  image(1:NT, 1:ny, t(P), col = cols, zlim = c(0, 1), xlab = "year", ylab = "seed", axes = FALSE,
        main = sprintf("Posterior P(disturbed) by seed (true disturbance in year %d)", td))
  axis(1, at = 1:NT); axis(2, at = 1:ny, labels = seeds, las = 1, cex.axis = 0.7)
  abline(v = td + c(-0.5, 0.5), col = "grey40")
  cl <- which(clip, arr.ind = TRUE)
  points(cl[, 2], cl[, 1] - 0.22, pch = 4, cex = 0.6)
  big <- which(P > fa.level & col(P) != td, arr.ind = TRUE)
  text(big[, 2], big[, 1] + 0.18, sprintf("%.2f", P[big]), cex = 0.5)
  legend("topright", inset = c(-0.16, 0), xpd = NA, bty = "n", cex = 0.75,
         legend = c("obs clipped to 0", sprintf("P > %.1f labelled", fa.level)), pch = c(4, NA))
  ## colour key
  usr <- par("usr"); x0 <- usr[2] + 0.6; yk <- seq(usr[3], usr[3] + 0.5 * (usr[4] - usr[3]), length.out = 51)
  rect(x0, yk[-51], x0 + 0.5, yk[-1], col = cols, border = NA, xpd = NA)
  text(x0 + 0.6, yk[c(1, 26, 51)], c("0", "0.5", "1"), adj = 0, cex = 0.7, xpd = NA)
}

png(file.path(sweep.dir, "pdist_by_seed.png"), width = 1200, height = 760, res = 115)
heatmap_pdist()
dev.off()

pdf(file.path(sweep.dir, "pdist_summary.pdf"), width = 9, height = 6)
heatmap_pdist()

## P(disturbed) by year across seeds
par(mar = c(4.2, 4.2, 3, 1))
qq <- apply(P, 2, quantile, c(0.5, 0.9, 1))
plot(1:NT, colMeans(P), type = 'n', ylim = c(0, 1), xlab = "year", ylab = "P(disturbed)",
     main = "P(disturbed) across seeds")
for (s in seq_along(seeds)) lines(1:NT, P[s, ], col = adjustcolor("grey50", 0.4))
lines(1:NT, colMeans(P), lwd = 2.5, col = "#eb6834")
lines(1:NT, qq[2, ], lwd = 1.5, lty = 2, col = "#eb6834")
lines(1:NT, P["42", ], lwd = 2, col = "#2a78d6")
abline(h = res[[1]]$p.dist, lty = 3); abline(v = td, lty = 2)
legend("topright", bty = "n", cex = 0.85,
       legend = c("each seed", "mean", "90th percentile", "seed 42", "prior p.dist"),
       col = c(adjustcolor("grey50", 0.6), "#eb6834", "#eb6834", "#2a78d6", 1),
       lwd = c(1, 2.5, 1.5, 2, 1), lty = c(1, 1, 2, 1, 3))

## tau learned vs tau known, every seed-year outside the disturbance
plot(PK[und], P[und], pch = 19, cex = 0.6, col = ifelse(clip[und], "#e41a1c", adjustcolor(1, 0.5)),
     xlim = c(0, 1), ylim = c(0, 1), xlab = "P(disturbed), tau known", ylab = "P(disturbed), tau learned",
     main = "False-alarm probabilities: tau learned vs known")
abline(0, 1, lty = 2)
legend("topleft", bty = "n", pch = 19, col = c("#e41a1c", adjustcolor(1, 0.5)),
       legend = c("an observation clipped to 0", "no clipped observation"))
dev.off()

saveRDS(list(seeds = seeds, P = P, PK = PK, clip = clip, z = z, summary = summ, false.alarms = fat,
             runs = lapply(res, function(x) x[setdiff(names(x), "log")])),
        file.path(sweep.dir, "sweep_results.rds"))
