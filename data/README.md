# data

Only the small inputs in `data/inputs/` and `data/paper/` are tracked. Everything else is
large or can be rebuilt. It stays local and is ignored by git (see `.gitignore`).

## Tracked

| file | contents | made by |
|---|---|---|
| `inputs/STprior.RData` | `d.prior` (background disturbance-type rates), `TRANS`, `PBBAR`, `PBVAR` | `legacy/models/StateTransition.Rmd` |
| `inputs/extract.h03v04.RData` | `points` (site coordinates), `bextract` (biomass 1990–2017), `cextract`, `dextract` (disturbance product 1985–2012) | `legacy/models/StateTransition.Rmd` |
| `inputs/dis.postprocessed.RData` | `dist.priors` (disturbance-size mixtures per class), `distlabels` | `legacy/models/RelabelDisturbances.Rmd` |
| `inputs/fit.RData` | VSEM calibration: `theta`, process sds `Ql`, `Qw`, `Qs` | `legacy/models/calibrate.Rmd` |
| `paper/rat.RData`, `paper/SSSD.RData` | inputs to `legacy/paper/methods.Rmd` | `StateTransition.Rmd`, `SingleSite_SingleDisturb.Rmd` |

The OR model (`or-model/`) needs only `data/inputs/`. `or-model/select-sites-or.R` writes
new site sets here as `inputs/extract.h03v04-n<nsamp>-seed<seed>.RData`. Track one of those
if runs on it should be reproducible from a clone.

## Local only (ignored)

| folder | contents |
|---|---|
| `cache/` | `daymet.h03v04.RData` (daymet per site, ~100 MB for 540 sites; `or-model/load-data-or.R` downloads the sites it is missing and adds them) and `or-model-inputs.RData` (everything `load-data-or.R` builds) |
| `raw/` | source rasters for `StateTransition.Rmd` and `or-model/select-sites-or.R`: `Cmonster_disturbance_class_ARD_tile_h03v04.tif`, `conus_biomass_ARD_tile_h03v04.tif`, `conus_land_cover_ARD_tile_h03v04.tif` (0.7–1.4 GB each) |
| `derived/` | rasters and intermediate files made from `raw/` (`disR.tif` is needed by `select-sites-or.R`): `bio04.tif`, `cov04.tif`, `dis04.tif`, `disR.tif`, `dis.processed.tif`, `disR.resamp.RData`, `dis.h03v04.RData`, `stateTrain.RData`, `state-transition-saved.RData` |
| `legacy-results/` | outputs of the legacy notebooks: `exp2–exp5trim.RData`, `elist.RData`, `stats.RData`, `Analytical*.RData`, `multidisturb-ex-values.RData`, `test.RData`, `debugging-or.RData` |

`bio04.tif`, `cov04.tif`, `exp2.RData`, `exp4.RData` and `exp5trim.RData` used to be tracked.
They are 7–67 MB each and are made by the legacy notebooks, so they are no longer tracked.
They are still in the git history.
