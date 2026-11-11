## Numbers and figures for claude-factorized-tau-11.11.26.tex.
##
## For each factorized-tau experiment and seeds 42, 1..20: runs the 3-D grid filter
## (ORJointUpdateGeneral, which stores the joint tau posterior) and the factorized filter on
## the same set-up and importance-sample seed, and records
##   - the dependence of the joint tau posterior at every step: correlation of
##     (log tau_a, log tau_b) and mutual information of (tau_a, tau_b), for each pair
##   - max over years of |P_fact(disturbed) - P_joint(disturbed)|
##   - process sd sqrt(E[1/tau]) at the last year under each, and the posterior sd of
##     log sigma under the 3-D grid
##   - the same for a second 3-D grid run with other importance samples (seed + 1000), as a
##     baseline for the filter's own Monte Carlo noise
##   - run times
## Each run is its own Rscript process (the set-ups use global variables).
##
## Run from the repository root:
##   Rscript claude-working-docs/claude-factorized-tau-11.11.26-analysis.R [--jobs 8]

source("synthetic-experiments/factorized-tau/or-joint-update-factorized.R")
experiments <- c("singlesite-singledisturb-multipool", "singlesite-multidisturb-multipool",
                 "singlesite-multidisturb-multipool-multipft", "multisite-singledisturb-multipool")
short <- c("1 site, 1 dist.", "1 site, 2 dist.", "1 site, PFT switch", "2 sites, 1 dist.")
seeds <- c(42, 1:20)
out.dir <- "claude-working-docs/claude-factorized-tau-11.11.26-files"
args <- commandArgs(trailingOnly = TRUE)
opt  <- function(name) { i <- which(args == name); if (length(i)) args[i + 1] else NULL }
jobs <- if (!is.null(opt("--jobs"))) as.integer(opt("--jobs")) else 8L

## dependence of a joint posterior on the grid: idx (G x 3 axis indices), w (weights)
tau_dep <- function(tj, tauAxes) {
  lt <- sapply(1:3, function(a) log(tauAxes[[a]][tj$idx[, a]]))
  w  <- tj$w / sum(tj$w)
  pr <- rbind(c(1, 2), c(1, 3), c(2, 3))
  t(apply(pr, 1, function(p) {
    a <- p[1]; b <- p[2]
    m <- colSums(w * lt[, p]); cv <- sum(w * (lt[, a] - m[1]) * (lt[, b] - m[2]))
    va <- sum(w * (lt[, a] - m[1])^2); vb <- sum(w * (lt[, b] - m[2])^2)
    cr <- if (va > 0 && vb > 0) cv / sqrt(va * vb) else 0
    J  <- tapply(w, list(tj$idx[, a], tj$idx[, b]), sum); J[is.na(J)] <- 0
    pa <- rowSums(J); pb <- colSums(J); nz <- J > 0
    mi <- sum(J[nz] * log(J[nz] / outer(pa, pb)[nz]))
    c(cor = cr, mi = mi)
  }))
}

if ("--worker" %in% args) {
  exp <- opt("--worker"); seed <- as.numeric(opt("--seed")); file <- opt("--out")
  source_until(file.path("synthetic-experiments/factorized-tau", paste0(exp, ".R")),
               "^cmp  <- run_comparison", envir = globalenv())
  fa <- modifyList(filter_args, list(forecastQuantiles = FALSE, keepMix = integer(0)))
  set.seed(seed); tf <- system.time(ff <- do.call(ORJointUpdateFactorized, c(fa, list(grid = grid))))[3]
  set.seed(seed); tj <- system.time(fj <- do.call(ORJointUpdateGeneral,    c(fa, list(grid = grid))))[3]
  ## baseline: the 3-D grid filter again with other importance samples (its own Monte Carlo noise)
  set.seed(seed + 1000); fj2 <- do.call(ORJointUpdateGeneral, c(fa, list(grid = grid)))
  dep <- lapply(fj$tauJoint, tau_dep, tauAxes = fj$tauAxes)
  pd  <- function(f) 1 - f$pK[, 1]
  NT  <- nrow(fj$procVar)
  ## posterior sd of log sigma_a = -log(tau_a)/2 at the last year, from the 3-D grid marginals
  sdlog <- sapply(1:3, function(a) { w <- fj$tauMarg[[a]][, NT]; l <- -0.5 * log(fj$tauAxes[[a]])
                                     sqrt(sum(w * l^2) - sum(w * l)^2) })
  saveRDS(list(exp = exp, seed = seed, time = c(fact = unname(tf), joint = unname(tj)),
               cor = sapply(dep, function(d) d[, "cor"]), mi = sapply(dep, function(d) d[, "mi"]),
               dP = max(abs(pd(ff) - pd(fj))), dP.mc = max(abs(pd(fj2) - pd(fj))),
               pFact = pd(ff), pJoint = pd(fj), pJoint2 = pd(fj2),
               sd = rbind(fact = sqrt(ff$procVar[NT, ]), joint = sqrt(fj$procVar[NT, ]), joint2 = sqrt(fj2$procVar[NT, ])),
               sdlog = sdlog, Gfact = sum(grid$G), Gjoint = prod(grid$G)), file)
  quit(save = "no")
}

dir.create(file.path(out.dir, "runs"), recursive = TRUE, showWarnings = FALSE)
todo <- expand.grid(seed = seeds, exp = experiments, stringsAsFactors = FALSE)
res <- parallel::mclapply(seq_len(nrow(todo)), function(i) {
  f <- file.path(out.dir, "runs", sprintf("%s_seed%02d.rds", todo$exp[i], todo$seed[i]))
  if (!file.exists(f)) {
    log <- sub("\\.rds$", ".log", f)
    system2("Rscript", c("claude-working-docs/claude-factorized-tau-11.11.26-analysis.R", "--worker", todo$exp[i],
                         "--seed", todo$seed[i], "--out", f), stdout = log, stderr = log)
  }
  if (file.exists(f)) readRDS(f) else NULL
}, mc.cores = jobs, mc.preschedule = FALSE)
res <- Filter(Negate(is.null), res)
message(sprintf("%d of %d runs", length(res), nrow(todo)))

## ---- tables --------------------------------------------------------------------------------
lr <- function(r, m) abs(log(r$sd[m, ] / r$sd["joint", ]))          # per pool
by_exp <- split(res, sapply(res, `[[`, "exp"))[experiments]
tab <- t(sapply(by_exp, function(rs) {
  cr <- unlist(lapply(rs, `[[`, "cor")); mi <- unlist(lapply(rs, `[[`, "mi"))
  lf <- unlist(lapply(rs, lr, "fact")); lj <- unlist(lapply(rs, lr, "joint2"))
  z  <- unlist(lapply(rs, function(r) lr(r, "fact") / r$sdlog))
  c(G.fact = rs[[1]]$Gfact, G.joint = rs[[1]]$Gjoint,
    speedup = median(sapply(rs, function(r) r$time[["joint"]] / r$time[["fact"]])),
    abs.cor.med = median(abs(cr)), abs.cor.max = max(abs(cr)),
    mi.med = median(mi), mi.max = max(mi),
    dP.med = median(sapply(rs, `[[`, "dP")), dP.max = max(sapply(rs, `[[`, "dP")),
    dPmc.med = median(sapply(rs, `[[`, "dP.mc")), dPmc.max = max(sapply(rs, `[[`, "dP.mc")),
    dsd.med = median(lf), dsd.max = max(lf), dsdmc.med = median(lj), dsdmc.max = max(lj),
    z.med = median(z), z.max = max(z), z.gt1 = mean(z > 1))
}))
rownames(tab) <- short
write.csv(round(tab, 4), file.path(out.dir, "table.csv"))
print(round(tab, 3))

## LaTeX rows
f2 <- function(x) sprintf("%.2f", x); f3 <- function(x) sprintf("%.3f", x)
writeLines(sprintf("%s & %s & %s & %s & %s \\\\", short, f2(tab[, "abs.cor.med"]), f2(tab[, "abs.cor.max"]),
                   f3(tab[, "mi.med"]), f3(tab[, "mi.max"])), file.path(out.dir, "rows-dep.tex"))
writeLines(sprintf("%s & %d & %d & %.0f$\\times$ \\\\", short, tab[, "G.fact"], tab[, "G.joint"], tab[, "speedup"]),
           file.path(out.dir, "rows-grid.tex"))
writeLines(sprintf("%s & %s & %s & %s & %s \\\\", short, f3(tab[, "dP.med"]), f3(tab[, "dP.max"]),
                   f3(tab[, "dPmc.med"]), f3(tab[, "dPmc.max"])), file.path(out.dir, "rows-agreeP.tex"))
writeLines(sprintf("%s & %s & %s & %s & %s & %s & %.0f\\%% \\\\", short, f2(tab[, "dsd.med"]), f2(tab[, "dsd.max"]),
                   f2(tab[, "dsdmc.med"]), f2(tab[, "dsdmc.max"]), f2(tab[, "z.med"]), 100 * tab[, "z.gt1"]),
           file.path(out.dir, "rows-agreeSD.tex"))

## ---- figures -------------------------------------------------------------------------------
cols <- c("#1b9e77", "#d95f02", "#7570b3", "#e7298a")
pdf(file.path(out.dir, "dependence.pdf"), width = 9, height = 4)
par(mfrow = c(1, 2), mar = c(4.2, 4.4, 2.4, 1), cex.lab = 1.2, mgp = c(2.6, 1.6, 0))
for (what in c("cor", "mi")) {
  vals <- lapply(by_exp, function(rs) as.vector(sapply(rs, function(r) if (what == "cor") abs(r$cor) else r$mi)))
  boxplot(vals, names = sub(", ", "\n", short), col = adjustcolor(cols, 0.5), border = cols, outcex = 0.4, las = 1, cex.axis = 0.75,
          ylab = if (what == "cor") expression("|" * rho[ab] * "|") else "MI (nats)",
          main = if (what == "cor") "|correlation| of log tau pairs" else "mutual information of tau pairs")
}
dev.off()

pdf(file.path(out.dir, "agreement.pdf"), width = 9, height = 4)
par(mfrow = c(1, 2), mar = c(4.2, 4.4, 2.4, 1), cex.lab = 1.2)
plot(NA, xlim = c(0, 1), ylim = c(0, 1), xlab = "P(disturbed), 3-D grid", ylab = "P(disturbed), factorized",
     main = "Per-year event probability, all seeds")
abline(0, 1, col = "grey60")
for (i in seq_along(by_exp)) for (r in by_exp[[i]]) points(r$pJoint, r$pFact, pch = 19, cex = 0.4, col = adjustcolor(cols[i], 0.6))
legend("topleft", legend = short, col = cols, pch = 19, bty = "n", cex = 0.8)
sdj <- unlist(lapply(by_exp, function(rs) sapply(rs, function(r) r$sd["joint", ])))
sdj2 <- unlist(lapply(by_exp, function(rs) sapply(rs, function(r) r$sd["joint2", ])))
sdf <- unlist(lapply(by_exp, function(rs) sapply(rs, function(r) r$sd["fact", ])))
ce  <- rep(cols, sapply(by_exp, function(rs) 3 * length(rs)))
rng <- range(c(sdj, sdf, sdj2))
plot(sdj, sdf, log = "xy", xlim = rng, ylim = rng, pch = 19, cex = 0.6, col = adjustcolor(ce, 0.7),
     xlab = "process sd, 3-D grid", ylab = "process sd, other run",
     main = "sqrt(E[1/tau]) at year 25, all seeds and pools")
points(sdj, sdj2, pch = 4, cex = 0.6, col = adjustcolor("grey30", 0.6))
abline(0, 1, col = "grey60")
legend("topleft", legend = c("factorized", "3-D grid, other samples"), pch = c(19, 4), col = c("grey30", "grey30"), bty = "n", cex = 0.8)
dev.off()
