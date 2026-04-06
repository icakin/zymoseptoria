# =============================================================================
# 00_config.R — Shared paths, constants, and helper functions
# =============================================================================
# Sourced by every script in the pipeline.
# Edit base_dir (and any settings below) in this one place only.
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
})

# ===== Project root and output directories ====================================

base_dir    <- "/Users/ilgazcakin/Desktop/Projects/Zymoseptoria"
data_dir    <- file.path(base_dir, "data")
tables_dir  <- file.path(base_dir, "tables")
figures_dir <- file.path(base_dir, "figures")
models_dir  <- file.path(base_dir, "models")

for (d in c(data_dir, tables_dir, figures_dir, models_dir)) {
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
}

# ===== Input files ============================================================

IN_CSV        <- file.path(tables_dir, "Oxygen_Data_Filtered.csv")
TRIM_META_CSV <- file.path(tables_dir, "Oxygen_Trimmed_Series_Metadata.csv")

# ===== Stage 1 output paths ==================================================

pdf_path         <- file.path(figures_dir, "per_series_fits.pdf")
coef_csv         <- file.path(tables_dir, "fit_coefficients_long.csv")
fit_metrics_csv  <- file.path(tables_dir, "fit_metrics.csv")
coef_wide_csv    <- file.path(tables_dir, "fit_coefficients_wide.csv")
derived_csv      <- file.path(tables_dir, "derived_N0_R_results_with_carbon.csv")
group_lookup_csv <- file.path(tables_dir, "group_lookup_with_delta.csv")
dose_key_csv     <- file.path(tables_dir, "dose_key_lookup.csv")

replication_summary_csv <- file.path(tables_dir,
                                     "replication_summary_by_temperature_and_dose.csv")

# ===== Bayesian Arrhenius output paths ========================================

brms_growth_rds  <- file.path(models_dir, "brms_arrhenius_growth_fgC_h_by_dose.rds")
brms_resp_rds    <- file.path(models_dir, "brms_arrhenius_resp_fgC_h_by_dose.rds")
brms_ratio_rds   <- file.path(models_dir, "brms_arrhenius_log_resp_over_growth_by_dose.rds")

post_growth_csv  <- file.path(tables_dir, "posterior_arrhenius_growth_fgC_h_by_dose.csv")
post_resp_csv    <- file.path(tables_dir, "posterior_arrhenius_resp_fgC_h_by_dose.csv")
post_ratio_csv   <- file.path(tables_dir, "posterior_arrhenius_log_resp_over_growth_by_dose.csv")

summary_growth_csv <- file.path(tables_dir, "summary_arrhenius_growth_fgC_h_by_dose.csv")
summary_resp_csv   <- file.path(tables_dir, "summary_arrhenius_resp_fgC_h_by_dose.csv")
summary_ratio_csv  <- file.path(tables_dir, "summary_arrhenius_log_resp_over_growth_by_dose.csv")

# CUE common-slope Arrhenius outputs
brms_cue_rds    <- file.path(models_dir, "brms_arrhenius_log_CUE_common_slope.rds")
post_cue_csv    <- file.path(tables_dir, "posterior_arrhenius_log_CUE_common_slope.csv")
summary_cue_csv <- file.path(tables_dir, "summary_arrhenius_log_CUE_common_slope.csv")

# Biomass-corrected Arrhenius outputs
brms_growth_biomass_rds   <- file.path(models_dir, "brms_arrhenius_growth_C_per_C_h_by_dose.rds")
brms_resp_biomass_rds     <- file.path(models_dir, "brms_arrhenius_respiration_C_per_C_h_by_dose.rds")

post_growth_biomass_csv   <- file.path(tables_dir, "posterior_arrhenius_growth_C_per_C_h_by_dose.csv")
post_resp_biomass_csv     <- file.path(tables_dir, "posterior_arrhenius_respiration_C_per_C_h_by_dose.csv")

summary_growth_biomass_csv <- file.path(tables_dir, "summary_arrhenius_growth_C_per_C_h_by_dose.csv")
summary_resp_biomass_csv   <- file.path(tables_dir, "summary_arrhenius_respiration_C_per_C_h_by_dose.csv")

# ===== Bayesian Arrhenius plot outputs ========================================

bayes_growth_png       <- file.path(figures_dir, "bayesian_arrhenius_growth_fgC_h_by_dose.png")
bayes_resp_png         <- file.path(figures_dir, "bayesian_arrhenius_respiration_fgC_h_by_dose.png")
bayes_ratio_png        <- file.path(figures_dir, "bayesian_arrhenius_resp_over_growth_by_dose.png")
bayes_boltz_growth_png <- file.path(figures_dir, "bayesian_boltzmann_growth_fgC_h_by_dose.png")
bayes_boltz_resp_png   <- file.path(figures_dir, "bayesian_boltzmann_respiration_fgC_h_by_dose.png")
bayes_boltz_ratio_png  <- file.path(figures_dir, "bayesian_boltzmann_resp_over_growth_by_dose.png")
bayes_boltz_cue_png    <- file.path(figures_dir, "bayesian_boltzmann_log_CUE_common_slope.png")

# Biomass-corrected Arrhenius/Boltzmann plot outputs
bayes_growth_biomass_png       <- file.path(figures_dir, "bayesian_arrhenius_growth_C_per_C_h_by_dose.png")
bayes_resp_biomass_png         <- file.path(figures_dir, "bayesian_arrhenius_respiration_C_per_C_h_by_dose.png")
bayes_boltz_growth_biomass_png <- file.path(figures_dir, "bayesian_boltzmann_growth_C_per_C_h_by_dose.png")
bayes_boltz_resp_biomass_png   <- file.path(figures_dir, "bayesian_boltzmann_respiration_C_per_C_h_by_dose.png")

# ===== Bayesian Sharpe-Schoolfield full TPC outputs ===========================

brms_tpc_growth_rds    <- file.path(models_dir, "brms_sharpe_schoolfield_tpc_growth_fgC_h_by_dose.rds")
brms_tpc_resp_rds      <- file.path(models_dir, "brms_sharpe_schoolfield_tpc_respiration_fgC_h_by_dose.rds")

post_tpc_growth_csv    <- file.path(tables_dir, "posterior_sharpe_schoolfield_tpc_growth_fgC_h_by_dose.csv")
post_tpc_resp_csv      <- file.path(tables_dir, "posterior_sharpe_schoolfield_tpc_respiration_fgC_h_by_dose.csv")

summary_tpc_growth_csv <- file.path(tables_dir, "summary_sharpe_schoolfield_tpc_growth_fgC_h_by_dose.csv")
summary_tpc_resp_csv   <- file.path(tables_dir, "summary_sharpe_schoolfield_tpc_respiration_fgC_h_by_dose.csv")

bayes_tpc_growth_png   <- file.path(figures_dir, "bayesian_sharpe_schoolfield_tpc_growth_fgC_h_by_dose.png")
bayes_tpc_resp_png     <- file.path(figures_dir, "bayesian_sharpe_schoolfield_tpc_respiration_fgC_h_by_dose.png")

# Biomass-corrected Sharpe-Schoolfield outputs
brms_tpc_growth_biomass_rds    <- file.path(models_dir, "brms_sharpe_schoolfield_tpc_growth_C_per_C_h_by_dose.rds")
brms_tpc_resp_biomass_rds      <- file.path(models_dir, "brms_sharpe_schoolfield_tpc_respiration_C_per_C_h_by_dose.rds")

post_tpc_growth_biomass_csv    <- file.path(tables_dir, "posterior_sharpe_schoolfield_tpc_growth_C_per_C_h_by_dose.csv")
post_tpc_resp_biomass_csv      <- file.path(tables_dir, "posterior_sharpe_schoolfield_tpc_respiration_C_per_C_h_by_dose.csv")

summary_tpc_growth_biomass_csv <- file.path(tables_dir, "summary_sharpe_schoolfield_tpc_growth_C_per_C_h_by_dose.csv")
summary_tpc_resp_biomass_csv   <- file.path(tables_dir, "summary_sharpe_schoolfield_tpc_respiration_C_per_C_h_by_dose.csv")

bayes_tpc_growth_biomass_png   <- file.path(figures_dir, "bayesian_sharpe_schoolfield_tpc_growth_C_per_C_h_by_dose.png")
bayes_tpc_resp_biomass_png     <- file.path(figures_dir, "bayesian_sharpe_schoolfield_tpc_respiration_C_per_C_h_by_dose.png")

# Posterior distribution plots for SS parameters
postdist_tpc_growth_png         <- file.path(figures_dir, "posterior_distributions_sharpe_schoolfield_growth.png")
postdist_tpc_resp_png           <- file.path(figures_dir, "posterior_distributions_sharpe_schoolfield_respiration.png")
postdist_tpc_growth_biomass_png <- file.path(figures_dir, "posterior_distributions_sharpe_schoolfield_growth_C_per_C_h.png")
postdist_tpc_resp_biomass_png   <- file.path(figures_dir, "posterior_distributions_sharpe_schoolfield_respiration_C_per_C_h.png")

# ===== Other plot outputs =====================================================

box_growth_png             <- file.path(figures_dir, "boxplot_growth_fgC_h.png")
box_resp_png               <- file.path(figures_dir, "boxplot_respiration_fgC_h.png")
scatter_growth_png         <- file.path(figures_dir, "scatter_growth_fgC_h_vs_T.png")
scatter_resp_png           <- file.path(figures_dir, "scatter_respiration_fgC_h_vs_T.png")
resp_over_growth_png       <- file.path(figures_dir, "resp_over_growth_vs_T.png")
cue_vs_T_png               <- file.path(figures_dir, "CUE_vs_T.png")
box_growth_biomass_png     <- file.path(figures_dir, "boxplot_growth_C_per_C_h.png")
box_resp_biomass_png       <- file.path(figures_dir, "boxplot_respiration_C_per_C_h.png")
scatter_growth_biomass_png <- file.path(figures_dir, "scatter_growth_C_per_C_h_vs_T.png")
scatter_resp_biomass_png   <- file.path(figures_dir, "scatter_respiration_C_per_C_h_vs_T.png")

# Per-dose TPC check plots
tpc_check_growth_png         <- file.path(figures_dir, "tpc_check_growth_fgC_h_by_dose_faceted.png")
tpc_check_resp_png           <- file.path(figures_dir, "tpc_check_respiration_fgC_h_by_dose_faceted.png")
tpc_check_growth_biomass_png <- file.path(figures_dir, "tpc_check_growth_C_per_C_h_by_dose_faceted.png")
tpc_check_resp_biomass_png   <- file.path(figures_dir, "tpc_check_respiration_C_per_C_h_by_dose_faceted.png")

# Log-scale TPC check plots
tpc_check_growth_log_png         <- file.path(figures_dir, "tpc_check_growth_fgC_h_log_by_dose_faceted.png")
tpc_check_resp_log_png           <- file.path(figures_dir, "tpc_check_respiration_fgC_h_log_by_dose_faceted.png")
tpc_check_growth_biomass_log_png <- file.path(figures_dir, "tpc_check_growth_C_per_C_h_log_by_dose_faceted.png")
tpc_check_resp_biomass_log_png   <- file.path(figures_dir, "tpc_check_respiration_C_per_C_h_log_by_dose_faceted.png")

# Temperature-specific fungicide effect sizes
fungicide_temp_effect_summary_csv <- file.path(tables_dir, "fungicide_effect_sizes_by_temperature_summary.csv")
fungicide_temp_effect_draws_csv   <- file.path(tables_dir, "fungicide_effect_sizes_by_temperature_draws.csv")
fungicide_temp_effect_png         <- file.path(figures_dir, "fungicide_effect_sizes_by_temperature.png")
fungicide_slope_test_summary_csv  <- file.path(tables_dir, "fungicide_effect_slope_tests_by_temperature_summary.csv")
fungicide_slope_test_draws_csv    <- file.path(tables_dir, "fungicide_effect_slope_tests_by_temperature_draws.csv")

# ===== Hierarchical Arrhenius sensitivity analysis outputs ====================

brms_growth_hier_rds <- file.path(models_dir, "brms_arrhenius_growth_fgC_h_by_dose_hier_cond.rds")
brms_resp_hier_rds   <- file.path(models_dir, "brms_arrhenius_resp_fgC_h_by_dose_hier_cond.rds")
brms_ratio_hier_rds  <- file.path(models_dir, "brms_arrhenius_log_resp_over_growth_by_dose_hier_cond.rds")

brms_growth_biomass_hier_rds <- file.path(models_dir, "brms_arrhenius_growth_C_per_C_h_by_dose_hier_cond.rds")
brms_resp_biomass_hier_rds   <- file.path(models_dir, "brms_arrhenius_respiration_C_per_C_h_by_dose_hier_cond.rds")

summary_growth_hier_csv <- file.path(tables_dir, "summary_arrhenius_growth_fgC_h_by_dose_hier_cond.csv")
summary_resp_hier_csv   <- file.path(tables_dir, "summary_arrhenius_resp_fgC_h_by_dose_hier_cond.csv")
summary_ratio_hier_csv  <- file.path(tables_dir, "summary_arrhenius_log_resp_over_growth_by_dose_hier_cond.csv")

summary_growth_biomass_hier_csv <- file.path(tables_dir, "summary_arrhenius_growth_C_per_C_h_by_dose_hier_cond.csv")
summary_resp_biomass_hier_csv   <- file.path(tables_dir, "summary_arrhenius_respiration_C_per_C_h_by_dose_hier_cond.csv")

post_growth_hier_csv <- file.path(tables_dir, "posterior_arrhenius_growth_fgC_h_by_dose_hier_cond.csv")
post_resp_hier_csv   <- file.path(tables_dir, "posterior_arrhenius_resp_fgC_h_by_dose_hier_cond.csv")
post_ratio_hier_csv  <- file.path(tables_dir, "posterior_arrhenius_log_resp_over_growth_by_dose_hier_cond.csv")

post_growth_biomass_hier_csv <- file.path(tables_dir, "posterior_arrhenius_growth_C_per_C_h_by_dose_hier_cond.csv")
post_resp_biomass_hier_csv   <- file.path(tables_dir, "posterior_arrhenius_respiration_C_per_C_h_by_dose_hier_cond.csv")

bayes_growth_hier_png <- file.path(figures_dir, "bayesian_arrhenius_growth_fgC_h_by_dose_hier_cond.png")
bayes_resp_hier_png   <- file.path(figures_dir, "bayesian_arrhenius_resp_fgC_h_by_dose_hier_cond.png")
bayes_ratio_hier_png  <- file.path(figures_dir, "bayesian_arrhenius_log_resp_over_growth_by_dose_hier_cond.png")

bayes_growth_biomass_hier_png <- file.path(figures_dir, "bayesian_arrhenius_growth_C_per_C_h_by_dose_hier_cond.png")
bayes_resp_biomass_hier_png   <- file.path(figures_dir, "bayesian_arrhenius_respiration_C_per_C_h_by_dose_hier_cond.png")

bayes_boltz_growth_hier_png <- file.path(figures_dir, "bayesian_boltzmann_growth_fgC_h_by_dose_hier_cond.png")
bayes_boltz_resp_hier_png   <- file.path(figures_dir, "bayesian_boltzmann_resp_fgC_h_by_dose_hier_cond.png")
bayes_boltz_ratio_hier_png  <- file.path(figures_dir, "bayesian_boltzmann_log_resp_over_growth_by_dose_hier_cond.png")

bayes_boltz_growth_biomass_hier_png <- file.path(figures_dir, "bayesian_boltzmann_growth_C_per_C_h_by_dose_hier_cond.png")
bayes_boltz_resp_biomass_hier_png   <- file.path(figures_dir, "bayesian_boltzmann_respiration_C_per_C_h_by_dose_hier_cond.png")

# ===== User settings ==========================================================

RMSE_KEEP_THRESHOLD <- 0.5
allowed_doses <- NULL

# ===== USER INPUT =============================================================

N_inoculation_cells_per_L <- 1e7

# ===== Thermal constants ======================================================

k_B    <- 0.00008617
T_ref  <- 293.15
INV_KB <- 1 / k_B

# ===== Zymoseptoria carbon conversion =========================================

CELL_WIDTH_UM  <- 2.5
CELL_LENGTH_UM <- 70
CARBON_DENSITY_FG_PER_UM3 <- 100
RESPIRATORY_QUOTIENT <- 1

cell_radius_um <- CELL_WIDTH_UM / 2
cell_cyl_length_um <- max(CELL_LENGTH_UM - CELL_WIDTH_UM, 0)

CELL_VOLUME_UM3 <-
  pi * (cell_radius_um^2) * cell_cyl_length_um +
  (4 / 3) * pi * (cell_radius_um^3)

CELL_CARBON_FG_PER_CELL <- CELL_VOLUME_UM3 * CARBON_DENSITY_FG_PER_UM3
MG_TO_FG <- 1e12
MIN_TO_H <- 60

# ===== Bayesian settings ======================================================

BAYES_ITER   <- 4000
BAYES_WARMUP <- 2000
BAYES_CHAINS <- 4
BAYES_SEED   <- 123
BAYES_ADAPT  <- 0.95
BAYES_MAX_TD <- 15

# ===== Helper functions =======================================================

norm_dose <- function(x) {
  x_chr  <- as.character(x)
  x_trim <- stringr::str_trim(x_chr)
  ifelse(stringr::str_to_lower(x_trim) == "control", "Control", x_trim)
}

dose_levels_from_data <- function(x) {
  d <- unique(norm_dose(x))
  d <- d[!is.na(d)]

  tibble(Dose = d) %>%
    mutate(
      Dose_num  = suppressWarnings(as.numeric(Dose)),
      Dose_sort = if_else(is.na(Dose_num), -Inf, Dose_num)
    ) %>%
    arrange(Dose_sort, Dose) %>%
    pull(Dose)
}

make_dose_colors <- function(dose_levels) {
  n <- length(dose_levels)
  if (n <= 0) return(setNames(character(0), character(0)))
  cols <- grDevices::hcl.colors(n, palette = "Dark 3")
  names(cols) <- dose_levels
  if ("Control" %in% names(cols)) cols["Control"] <- "#1f77b4"
  cols
}

get_baseline <- function(x, head_max = 30, min_valid = 5) {
  x <- as.numeric(x)
  n <- length(x)
  if (n == 0) return(NA_real_)

  k <- min(head_max, n)
  hv <- x[seq_len(k)]
  hv <- hv[is.finite(hv)]
  if (length(hv) >= min_valid) {
    b <- median(hv, na.rm = TRUE)
    if (is.finite(b)) return(b)
  }
  suppressWarnings(mean(x, na.rm = TRUE))
}

resp_model <- function(r, K, t, O2_0) {
  O2_0 + (K / r) * (1 - exp(r * t))
}

posterior_summary_df <- function(fit, pars = NULL, probs = c(0.025, 0.5, 0.975)) {
  as.data.frame(
    posterior::summarise_draws(
      posterior::as_draws_df(fit, variable = pars),
      mean,
      sd,
      ~posterior::quantile2(.x, probs = probs)
    )
  )
}

safe_logit <- function(p, eps = 1e-8) {
  p2 <- pmin(pmax(p, eps), 1 - eps)
  log(p2 / (1 - p2))
}

pred_ss_log <- function(TK, lnB0, E, Eh, Th) {
  lnB0 -
    (E * 11604.51812) * ((1 / TK) - (1 / 293.15)) -
    log(1 + exp((Eh * 11604.51812) * ((1 / Th) - (1 / TK))))
}

pred_arr_log <- function(TK, alpha, E, k_B = 0.00008617, T_ref = 293.15) {
  boltz_shift <- (1 / (k_B * T_ref)) - (1 / (k_B * TK))
  alpha + E * boltz_shift
}

downsample_draws <- function(df, max_draws = NULL) {
  if (is.null(max_draws)) return(df)
  df %>%
    group_by(Dose) %>%
    slice_head(n = max_draws) %>%
    ungroup()
}

assign_temp_region <- function(T) {
  cut(
    T,
    breaks = c(-Inf, 19, 24, Inf),
    labels = c("Cool (15\u201319\u00b0C)", "Mid (20\u201324\u00b0C)", "Hot (25\u201328\u00b0C)"),
    include.lowest = TRUE,
    right = TRUE
  )
}

make_dir <- function(path) {
  dir.create(path, showWarnings = FALSE, recursive = TRUE)
  path
}

message("00_config.R loaded: base_dir = ", base_dir)
