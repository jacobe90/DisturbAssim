# DisturbAssim

This project explores how to apply iterative data assimilation to disturbance and land-use change.

## Layout

| folder | contents |
|---|---|
| `or-model/` | the OR joint-update model on the real data (tile h03v04): data loading, the filter, run scripts and analyses |
| `synthetic-experiments/` | the joint-update filters on synthetic VSEM experiments, including the factorized-tau filter |
| `R/` | helpers shared by both: the disturbance models (`disturbance.R`) and `StatsAuxiliary.R` |
| `data/` | inputs; only the small ones are tracked (see `data/README.md`) |
| `experiments/` | run outputs (mostly ignored) |
| `claude-working-docs/` | write-ups and the scripts that made them |
| `legacy/` | earlier model notebooks (`legacy/models`) and the paper (`legacy/paper`) |

## Running the OR model

Run every command from the repository root. Each script prints its progress (step timings,
runs finished, ETA). The header of each script lists all of its options.

### Scripts

| script | what it does |
|---|---|
| `or-model/select-sites-or.R` | draws a new set of sites from the rasters (needs the local rasters, see `data/README.md`) |
| `or-model/load-data-or.R` | builds the model inputs for a site set -> `data/cache/or-model-inputs.RData` |
| `or-model/run-model-or.R` | runs the filter on those inputs, then writes the analysis reports into the run folder |
| `or-model/load-data-and-run-model.R` | `load-data-or.R` then `run-model-or.R` in one R session (takes the options of both) |
| `or-model/analysis-or.R` | the analysis reports on their own, for any run folder (`--indir`) |

A run folder gets one `.rds` per run, plus the analysis reports: `exps_summary.pdf`,
`distprod_TRUE_vs_FALSE.pdf` (when runs exist both with and without the disturbance
product), `elist.RData` and `stats.RData`. With 10 or fewer sites it also gets
`site-plots.pdf`. Runs already in the folder are skipped, so an interrupted command can be
restarted as it is. The analysis covers every run in the folder.

Bias correction and propagating the process error through the disturbance are on by
default. `--bias off` and `--dist-proc-err off` turn them off.

### Regenerating the data

With the tracked sites (538 sites, `data/inputs/extract.h03v04.RData`; about 30 s once
daymet is cached):

```
Rscript or-model/load-data-or.R
```

With a new, larger site selection. `--nsamp` is the number of cells drawn per disturbance
class (5 classes), and about 70% of them are forest and kept. `--nsamp 150` gives about
540 sites, and `--nsamp 300` about 1,100.

```
Rscript or-model/select-sites-or.R --nsamp 300 --seed 1
Rscript or-model/load-data-or.R --extract data/inputs/extract.h03v04-n300-seed1.RData
```

The first `load-data-or.R` for a new site set downloads daymet for the sites that are not
cached yet, at about 4 s per site (about 40 min for 600 new sites). Runs on a site set other
than the tracked one go to `experiments/or-model/runs-<site set>` (e.g. `runs-n300-seed1`).

### Full grid

This covers all sites × 4 observation-error models × with and without the disturbance
product: 4,304 runs, about 18 min on 8 workers.

```
Rscript or-model/run-model-or.R
```

Or rebuild the inputs and run in one go:

```
Rscript or-model/load-data-and-run-model.R
```

Output goes to `experiments/or-model/runs`.

### Smaller grid

For example, the first 50 sites with every setting (400 runs, about 2 min):

```
Rscript or-model/run-model-or.R --sites 1:50 --outdir experiments/or-model/runs-sites1-50
```

Or a random 100 sites, with only the GEDI observation-error models:

```
Rscript or-model/run-model-or.R --sites "sort(sample(538, 100))" --R GEDI_het,GEDI_con --outdir experiments/or-model/runs-100
```

### Single test site

One site and one setting, a few seconds. This writes `site-plots.pdf` and `exps_summary.pdf`:

```
Rscript or-model/run-model-or.R --sites 300 --R GEDI_het --distprod TRUE --outdir experiments/or-model/test-site300
```

With and without the disturbance product, which also writes `distprod_TRUE_vs_FALSE.pdf`:

```
Rscript or-model/run-model-or.R --sites 300 --R GEDI_het --outdir experiments/or-model/test-site300
```

Give a single-site or small-grid run its own `--outdir`. The analysis reports cover every
run in the folder, so a test run in `experiments/or-model/runs` would be mixed into the
full-grid reports.

### Other options

`--tau-prior gamma,loguniform`, `--jobs 4`, `--overwrite` (rerun runs already on disk) and
`--no-analysis`. To write the reports again for an existing run folder:

```
Rscript or-model/analysis-or.R --indir experiments/or-model/runs
```
