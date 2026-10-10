# Common setup, sourced at the top of every numbered script. Run the
# scripts from the package root, e.g.
#   Rscript simulations/interval-censoring/02-imputation-vs-midpoint.R
# Environment variables:
#   N_REPS   replicates per scenario (default: each script's own; 2 = smoke test)
#   N_CORES  parallel workers (default: all cores but one)

suppressMessages(pkgload::load_all(".", quiet = TRUE))
SIM_DIR <- file.path("simulations", "interval-censoring")
RESULTS_DIR <- file.path(SIM_DIR, "results")
source(file.path(SIM_DIR, "R", "data-generating.R"))
source(file.path(SIM_DIR, "R", "helpers.R"))

n_reps_for <- function(default) as.integer(Sys.getenv("N_REPS", default))
N_CORES <- as.integer(Sys.getenv("N_CORES", max(1, parallel::detectCores() - 1)))
GIT_COMMIT <- tryCatch(system("git rev-parse --short HEAD", intern = TRUE), error = function(e) NA)
