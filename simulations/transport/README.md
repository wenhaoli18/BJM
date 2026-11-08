# Simulation study for `transportBJM()` (experimental)

Kept outside the package build (see `.Rbuildignore`).

- `pilot.R` — first small study (3 shift sizes x 10 replicates).
- `full.R` — full study: 13 scenarios (shift size x target size, label-shift
  violations, visit density, competing risks), 10 methods, landmarks 0 and 2.
  Run from a directory containing an `rlib/` library with this branch of BJM
  installed: `Rscript full.R <first_rep> <last_rep>`; set
  `OPENBLAS_NUM_THREADS=1` when running several processes in parallel.
- `summarize.R` — tables (`sim/out/table_*.csv`) and figures from the
  per-replicate results in `sim/full/`.

## Results (40 replicates, complete)

- `results/replicates/rep_XXX.rds` — per-replicate results (rep_000 = pilot run, excluded).
- `results/pilot/` — pilot study (k = 0.5, 1, 2; 10 replicates each).
- `results/tables_figures/` — table_alpha.csv, table_prediction.csv, fig1_alpha.png,
  fig2_rmse.png, summary.txt.

Write-up of the simulation section: https://claude.ai/artifact/5tPtRiPuQqB8ZFW7FuLMiR
