# =============================================================================
# run_all.R — Master pipeline runner
# =============================================================================
# Sources all scripts in order. Each script sources 00_config.R internally.
#
# Pipeline order:
#   01  Wide -> long oxygen data
#   02  Spline trimming + diagnostics
#   03  Oxygen model fits (nlsLM) -> carbon units + descriptive plots
#   04  Bayesian thermal models (Arrhenius, Sharpe-Schoolfield, hierarchical)
#   05  All plots, effect sizes, synergy analysis, publication figures
#   06  Dose-response Hill + log-linear models (reads fig1_extra CSV from 05)
#
# You can also run any script individually — each one loads its own
# dependencies via source("00_config.R").
# =============================================================================

script_dir <- if (requireNamespace("rstudioapi", quietly = TRUE) &&
                   rstudioapi::isAvailable() &&
                   nzchar(rstudioapi::getActiveDocumentContext()$path)) {
  dirname(rstudioapi::getActiveDocumentContext()$path)
} else {
  tryCatch(dirname(sys.frame(1)$ofile), error = function(e) getwd())
}

run_script <- function(name) {
  path <- file.path(script_dir, name)
  message("\n", strrep("=", 70))
  message("Running: ", name)
  message(strrep("=", 70), "\n")
  source(path, local = FALSE)
}

run_script("01_longdata.R")
run_script("02_trimming.R")
run_script("03_oxygen_fits.R")
run_script("04_bayesian_models.R")
run_script("05_plots_and_effects.R")
run_script("06_dose_response_brms.R")

message("\n", strrep("=", 70))
message("Pipeline complete.")
message(strrep("=", 70))
