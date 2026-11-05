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
