# Run every interval-censoring simulation in order, from the package root:
#   Rscript simulations/interval-censoring/run-all.R
#   N_REPS=2 Rscript simulations/interval-censoring/run-all.R   # smoke test
# Each script writes results/<name>.rds and results/<name>.txt.
scripts <- sort(list.files(file.path("simulations", "interval-censoring"), pattern = "^[0-9]{2}-.*\\.R$",
                           full.names = TRUE))
for (f in scripts) {
  cat("\n#####", basename(f), "#####\n")
  t0 <- Sys.time()
  status <- system2(file.path(R.home("bin"), "Rscript"), f)
  cat(sprintf("##### %s: %s, %.1f min\n", basename(f), if (status == 0) "done" else "FAILED",
              as.numeric(Sys.time() - t0, units = "mins")))
}
