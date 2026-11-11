## Meteorology: daymet shortwave for every site -> daily PAR.
##
## Ported from OR_joint_update2.Rmd (now legacy/models/).

## Daymet series for the sites in latlon (lon, lat columns), in the same order: rad[i, ]
## is used as site i's driver. Sites already in the cache (matched by coordinates) are read
## from it; the others are downloaded and added to it, so the cache can hold the sites of
## several site selections.
load_met <- function(latlon, metFile = or.files$met) {
  dm <- list()
  if (file.exists(metFile)) {
    message(sprintf("  reading cached daymet: %s", metFile))
    load(metFile)
  }
  key <- function(lon, lat) sprintf("%.3f_%.3f", lon, lat)
  have <- if (length(dm)) key(sapply(dm, function(d) d$longitude[1]), sapply(dm, function(d) d$latitude[1])) else character(0)
  want <- key(latlon[, 1], latlon[, 2])
  miss <- which(!want %in% have)
  message(sprintf("  %d of %d sites in the cache, %d to download", length(latlon[, 1]) - length(miss),
                  nrow(latlon), length(miss)))
  if (length(miss)) {
    message("  downloading daymet 1985-2017 (slow; the cache is saved every 25 sites)")
    t0 <- Sys.time()
    dir.create(dirname(metFile), recursive = TRUE, showWarnings = FALSE)
    for (m in seq_along(miss)) {
      i <- miss[m]
      dm[[length(dm) + 1]] <- daymetr::download_daymet(site = "h03v04", latlon[i, 2], latlon[i, 1], start = 1985,
                                                       end = 2017, internal = TRUE, simplify = TRUE, silent = TRUE)
      have <- c(have, want[i])
      if (m %% 25 == 0 || m == length(miss)) {
        save(dm, file = metFile)
        el <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
        message(sprintf("  downloaded %d / %d sites (elapsed %s, ETA %s)", m, length(miss),
                        fmt_secs(el), fmt_secs(el / m * (length(miss) - m))))
      }
    }
    message(sprintf("  saved %s (%d sites)", metFile, length(dm)))
  }
  dm[match(want, have)]
}

## Daily PAR (MJ m-2 day-1) for every site: rad (sites x days, 1985-2017, used for the
## spin-up) and radDA (sites x days, restricted to the DA years; daymet has 365 days
## every year, so each DA timestep is one contiguous 365-day block).
met_to_par <- function(dm, years = atime) {
  srad    <- subset(dm[[1]], measurement == "srad..W.m.2.")
  metTime <- as.POSIXct(paste(srad$year, srad$yday, sep = "-"), format = "%Y-%j")
  rad     <- matrix(NA, length(dm), length(metTime))   # [sites x days]
  for (i in seq_along(dm)) {
    s <- subset(dm[[i]], measurement == "srad..W.m.2.")
    rad[i, ] <- s$value * 1e-6 * 86400 * SW2PAR   # W m-2 -> MJ m-2 day-1 PAR
    ## QUESTION -- should the changing day length, from the daymet day length values, be
    ## used (instead of a fixed 24-hr day length)?
  }
  radDA <- rad[, srad$year %in% years, drop = FALSE]
  stopifnot(ncol(radDA) == 365 * length(years))
  list(rad = rad, radDA = radDA, metTime = metTime)
}
