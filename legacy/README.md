# legacy

These are earlier notebooks and scripts, kept for reference and not maintained.

- `models/`: the older models and experiments, including `OR_joint_update2.Rmd` and
  `OR_joint_update_analysis.Rmd`, which are now ported to `or-model/`. It also holds the
  data-preparation notebooks (`StateTransition.Rmd`, `RelabelDisturbances.Rmd`,
  `calibrate.Rmd`) and their helpers (`multi.ens.adj.R` and others).
- `paper/`: the paper sources (`Paper.Rmd` and its sections, `bibliography.bibtex`).

These files were written to run from the repository root, with every helper and data file
there. They still use those paths:

- `disturbance.R` and `StatsAuxiliary.R` are now in `R/`.
- The data files are in `data/` (see `data/README.md`).
- Helpers used only by these notebooks sit next to them in `models/`.
