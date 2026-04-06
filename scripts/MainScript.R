# ─────────────────────────────────────────────────────────────────────────────
# CLEAN K-MODEL PIPELINE
# STAGE 1: oxygen model kept separate (frequentist nlsLM fits)
# STAGE 2A: Bayesian Arrhenius/Boltzmann thermal models
#   - Bayesian Arrhenius for log(growth_fgC_h) by Dose
#   - Bayesian Arrhenius for log(resp_fgC_h) by Dose
#   - Bayesian Arrhenius for log(resp / growth) by Dose
#   - Bayesian Arrhenius for log(growth_C_per_C_h) by Dose
#   - Bayesian Arrhenius for log(respiration_C_per_C_h) by Dose
#   - Bayesian Arrhenius for log(CUE) with common slope
# STAGE 2B: Bayesian full TPC models
#   - Bayesian Sharpe–Schoolfield TPC for log(growth_fgC_h) by Dose
#   - Bayesian Sharpe–Schoolfield TPC for log(resp_fgC_h) by Dose
#   - Bayesian Sharpe–Schoolfield TPC for log(growth_C_per_C_h) by Dose
#   - Bayesian Sharpe–Schoolfield TPC for log(respiration_C_per_C_h) by Dose
# STAGE 2C: Separate hierarchical Arrhenius sensitivity analysis
#   - Same fixed-effects Arrhenius structure as Stage 2A
#   - Additional random intercept on alpha for cond_id = interaction(Dose, T)
#
# Additions implemented from feedback:
# 1) Per-dose faceted TPC check plots
# 2) Explicit replication summary showing all replicate-level Stage 1 rows
#    are retained for Stage 2 (no averaging by condition)
# 3) Log-scale TPC check plots for cleaner visualisation
# 4) Temperature-specific fungicide effect-size summaries and plots
# 5) Bayesian Boltzmann Arrhenius plot for log(CUE) with a common slope
# 6) Bayesian slope/intercept tests for Growth and CUE effect-size curves
# 7) CUE descriptive plot moved to Boltzmann x-axis
#
# Notes:
# - Arrhenius models describe the rising limb.
# - Sharpe–Schoolfield models provide a biologically interpretable full peaked curve.
# - Arrhenius predictor is centred and scaled as: (1 / (k_B * T_ref)) - (1 / (k_B * T))
# - Original fg C h^-1 outputs are kept unchanged.
# - Additional biomass-corrected outputs are added as C per C per h.
# - Main analysis remains non-hierarchical.
# - Hierarchical Arrhenius is a separate sensitivity analysis.
# ─────────────────────────────────────────────────────────────────────────────

rm(list = ls())
gc()

suppressPackageStartupMessages({
  library(tidyverse)
  library(minpack.lm)
  library(tibble)
  library(brms)
  library(posterior)
  library(grid)
})

options(mc.cores = max(1, parallel::detectCores() - 1))

# ===== Files ==================================================================

# Project root and output directories
base_dir    <- "/Users/ilgazcakin/Desktop/ZymoTest"
tables_dir  <- file.path(base_dir, "tables")
figures_dir <- file.path(base_dir, "figures")
models_dir  <- file.path(base_dir, "models")

IN_CSV        <- file.path(tables_dir, "Oxygen_Data_Filtered.csv")
TRIM_META_CSV <- file.path(tables_dir, "Oxygen_Trimmed_Series_Metadata.csv")

pdf_path         <- file.path(figures_dir, "per_series_fits.pdf")
coef_csv         <- file.path(tables_dir, "fit_coefficients_long.csv")
fit_metrics_csv  <- file.path(tables_dir, "fit_metrics.csv")
coef_wide_csv    <- file.path(tables_dir, "fit_coefficients_wide.csv")
derived_csv      <- file.path(tables_dir, "derived_N0_R_results_with_carbon.csv")
group_lookup_csv <- file.path(tables_dir, "group_lookup_with_delta.csv")
dose_key_csv     <- file.path(tables_dir, "dose_key_lookup.csv")

# Replication summary
replication_summary_csv <- file.path(tables_dir, "replication_summary_by_temperature_and_dose.csv")

# Bayesian Arrhenius outputs
brms_growth_rds       <- file.path(models_dir, "brms_arrhenius_growth_fgC_h_by_dose.rds")
brms_resp_rds         <- file.path(models_dir, "brms_arrhenius_resp_fgC_h_by_dose.rds")
brms_ratio_rds        <- file.path(models_dir, "brms_arrhenius_log_resp_over_growth_by_dose.rds")

post_growth_csv       <- file.path(tables_dir, "posterior_arrhenius_growth_fgC_h_by_dose.csv")
post_resp_csv         <- file.path(tables_dir, "posterior_arrhenius_resp_fgC_h_by_dose.csv")
post_ratio_csv        <- file.path(tables_dir, "posterior_arrhenius_log_resp_over_growth_by_dose.csv")

summary_growth_csv    <- file.path(tables_dir, "summary_arrhenius_growth_fgC_h_by_dose.csv")
summary_resp_csv      <- file.path(tables_dir, "summary_arrhenius_resp_fgC_h_by_dose.csv")
summary_ratio_csv     <- file.path(tables_dir, "summary_arrhenius_log_resp_over_growth_by_dose.csv")

# CUE common-slope Arrhenius outputs
brms_cue_rds        <- file.path(models_dir, "brms_arrhenius_log_CUE_common_slope.rds")
post_cue_csv        <- file.path(tables_dir, "posterior_arrhenius_log_CUE_common_slope.csv")
summary_cue_csv     <- file.path(tables_dir, "summary_arrhenius_log_CUE_common_slope.csv")

# Biomass-corrected Arrhenius outputs
brms_growth_biomass_rds       <- file.path(models_dir, "brms_arrhenius_growth_C_per_C_h_by_dose.rds")
brms_resp_biomass_rds         <- file.path(models_dir, "brms_arrhenius_respiration_C_per_C_h_by_dose.rds")

post_growth_biomass_csv       <- file.path(tables_dir, "posterior_arrhenius_growth_C_per_C_h_by_dose.csv")
post_resp_biomass_csv         <- file.path(tables_dir, "posterior_arrhenius_respiration_C_per_C_h_by_dose.csv")

summary_growth_biomass_csv    <- file.path(tables_dir, "summary_arrhenius_growth_C_per_C_h_by_dose.csv")
summary_resp_biomass_csv      <- file.path(tables_dir, "summary_arrhenius_respiration_C_per_C_h_by_dose.csv")

# Bayesian Arrhenius plot outputs
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

# Bayesian Sharpe–Schoolfield full TPC outputs
brms_tpc_growth_rds    <- file.path(models_dir, "brms_sharpe_schoolfield_tpc_growth_fgC_h_by_dose.rds")
brms_tpc_resp_rds      <- file.path(models_dir, "brms_sharpe_schoolfield_tpc_respiration_fgC_h_by_dose.rds")

post_tpc_growth_csv    <- file.path(tables_dir, "posterior_sharpe_schoolfield_tpc_growth_fgC_h_by_dose.csv")
post_tpc_resp_csv      <- file.path(tables_dir, "posterior_sharpe_schoolfield_tpc_respiration_fgC_h_by_dose.csv")

summary_tpc_growth_csv <- file.path(tables_dir, "summary_sharpe_schoolfield_tpc_growth_fgC_h_by_dose.csv")
summary_tpc_resp_csv   <- file.path(tables_dir, "summary_sharpe_schoolfield_tpc_respiration_fgC_h_by_dose.csv")

bayes_tpc_growth_png   <- file.path(figures_dir, "bayesian_sharpe_schoolfield_tpc_growth_fgC_h_by_dose.png")
bayes_tpc_resp_png     <- file.path(figures_dir, "bayesian_sharpe_schoolfield_tpc_respiration_fgC_h_by_dose.png")

# Biomass-corrected Sharpe–Schoolfield outputs
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

# Other plot outputs
box_growth_png             <- file.path(figures_dir, "boxplot_growth_fgC_h.png")
box_resp_png               <- file.path(figures_dir, "boxplot_respiration_fgC_h.png")
scatter_growth_png         <- file.path(figures_dir, "scatter_growth_fgC_h_vs_T.png")
scatter_resp_png           <- file.path(figures_dir, "scatter_respiration_fgC_h_vs_T.png")
resp_over_growth_png       <- file.path(figures_dir, "resp_over_growth_vs_T.png")
cue_vs_T_png               <- file.path(figures_dir, "CUE_vs_T.png")  # now saved as Boltzmann-x CUE descriptive plot
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

# Hierarchical Arrhenius sensitivity analysis outputs
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

# ===== Helpers ================================================================
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
    labels = c("Cool (15–19°C)", "Mid (20–24°C)", "Hot (25–28°C)"),
    include.lowest = TRUE,
    right = TRUE
  )
}

compute_vs_control_effect <- function(pred_df, value_col, dose_ref) {
  val_sym <- rlang::sym(value_col)
  
  base_df <- pred_df %>%
    transmute(Dose, .draw, T, val = !!val_sym)
  
  ctrl_df <- base_df %>%
    filter(Dose == dose_ref) %>%
    select(.draw, T, val_ctrl = val)
  
  base_df %>%
    left_join(ctrl_df, by = c(".draw", "T")) %>%
    mutate(effect = val - val_ctrl)
}

compute_region_deviation <- function(df) {
  overall_df <- df %>%
    group_by(.draw) %>%
    summarise(overall_mean_effect = mean(effect, na.rm = TRUE), .groups = "drop")
  
  region_df <- df %>%
    group_by(temp_region, .draw) %>%
    summarise(region_mean_effect = mean(effect, na.rm = TRUE), .groups = "drop")
  
  region_df %>%
    left_join(overall_df, by = ".draw") %>%
    mutate(deviation = region_mean_effect - overall_mean_effect)
}

summarise_region_deviation <- function(df, trait_name, effect_type = c("negative", "positive")) {
  effect_type <- match.arg(effect_type)
  
  out <- df %>%
    group_by(temp_region) %>%
    summarise(
      median_deviation = median(deviation, na.rm = TRUE),
      q2.5 = quantile(deviation, 0.025, na.rm = TRUE),
      q97.5 = quantile(deviation, 0.975, na.rm = TRUE),
      p_gt0 = mean(deviation > 0, na.rm = TRUE),
      p_lt0 = mean(deviation < 0, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(Trait = trait_name)
  
  if (effect_type == "negative") {
    out %>%
      mutate(
        p_synergy = p_lt0,
        p_antagonism = p_gt0
      ) %>%
      select(-p_gt0, -p_lt0)
  } else {
    out %>%
      mutate(
        p_synergy = p_gt0,
        p_antagonism = p_lt0
      ) %>%
      select(-p_gt0, -p_lt0)
  }
}

label_direction <- function(trait, q2.5, q97.5) {
  if (trait %in% c("Growth", "CUE")) {
    if (q97.5 < 0) return("Synergistic")
    if (q2.5 > 0) return("Antagonistic")
    return("Uncertain")
  } else {
    if (q2.5 > 0) return("Synergistic")
    if (q97.5 < 0) return("Antagonistic")
    return("Uncertain")
  }
}

make_plot_subtitle <- function(mode) {
  if (mode == "observed") {
    paste0(
      "Points show pooled regional deviations from the overall mean fungicide effect, ",
      "averaged over observed experimental temperatures only. Dashed line = average fungicide effect."
    )
  } else {
    paste0(
      "Points show pooled regional deviations from the overall mean fungicide effect, ",
      "averaged over the full posterior temperature grid. Dashed line = average fungicide effect."
    )
  }
}

make_tpc_check_plot <- function(dat, grid, title_txt, ylab_txt, log_scale = FALSE) {
  p <- ggplot(dat, aes(T, y_raw, color = Dose, fill = Dose)) +
    geom_point(size = 2.0, alpha = 0.9) +
    {
      if (nrow(grid) > 0) {
        geom_ribbon(
          data = grid,
          aes(x = T, ymin = raw_q2.5, ymax = raw_q97.5, fill = Dose),
          inherit.aes = FALSE,
          alpha = 0.18,
          colour = NA
        )
      }
    } +
    {
      if (nrow(grid) > 0) {
        geom_line(
          data = grid,
          aes(x = T, y = raw_q50, color = Dose),
          inherit.aes = FALSE,
          linewidth = 1
        )
      }
    } +
    facet_wrap(~Dose, scales = "free_y") +
    scale_color_manual(values = dose_cols, limits = dose_levels, drop = FALSE) +
    scale_fill_manual(values = dose_cols, limits = dose_levels, drop = FALSE) +
    labs(
      title = title_txt,
      x = "Temperature (°C)",
      y = ylab_txt,
      color = "Dose",
      fill = "Dose"
    ) +
    theme_classic(12) +
    theme(legend.position = "none")
  
  if (isTRUE(log_scale)) p <- p + scale_y_log10()
  p
}

summarise_temp_effect <- function(df, trait_name, effect_type = c("negative", "positive")) {
  effect_type <- match.arg(effect_type)
  
  out <- df %>%
    group_by(Dose, T) %>%
    summarise(
      median_effect = median(effect, na.rm = TRUE),
      q2.5 = quantile(effect, 0.025, na.rm = TRUE),
      q97.5 = quantile(effect, 0.975, na.rm = TRUE),
      p_gt0 = mean(effect > 0, na.rm = TRUE),
      p_lt0 = mean(effect < 0, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(Trait = trait_name)
  
  if (effect_type == "negative") {
    out %>%
      mutate(
        p_synergy = p_lt0,
        p_antagonism = p_gt0
      ) %>%
      select(-p_gt0, -p_lt0)
  } else {
    out %>%
      mutate(
        p_synergy = p_gt0,
        p_antagonism = p_lt0
      ) %>%
      select(-p_gt0, -p_lt0)
  }
}

make_temp_effect_plot <- function(plot_df) {
  ggplot(plot_df, aes(x = T, y = median_effect, colour = Dose, fill = Dose)) +
    geom_hline(yintercept = 0, linetype = 2, colour = "grey40", linewidth = 0.5) +
    geom_ribbon(aes(ymin = q2.5, ymax = q97.5), alpha = 0.18, colour = NA) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 1.8) +
    facet_wrap(~Trait, scales = "free_y", ncol = 2) +
    scale_color_manual(values = dose_cols, limits = dose_levels, drop = FALSE) +
    scale_fill_manual(values = dose_cols, limits = dose_levels, drop = FALSE) +
    labs(
      title = "Fungicide effect sizes by temperature",
      subtitle = paste0(
        "Effect size is posterior difference from Control at each temperature. ",
        "Growth, respiration, and respiration/growth are on the log scale; CUE is on the logit scale."
      ),
      x = "Temperature (°C)",
      y = "Effect size vs Control",
      colour = "Dose",
      fill = "Dose"
    ) +
    theme_classic(12)
}

fit_effect_slopes_from_draws <- function(df, trait_name) {
  df %>%
    group_by(Dose, .draw) %>%
    group_modify(~{
      dd <- .x %>% filter(is.finite(T), is.finite(effect))
      if (nrow(dd) < 2 || dplyr::n_distinct(dd$T) < 2) {
        return(tibble(intercept = NA_real_, slope = NA_real_))
      }
      mod <- lm(effect ~ T, data = dd)
      tibble(
        intercept = unname(coef(mod)[1]),
        slope     = unname(coef(mod)[2])
      )
    }) %>%
    ungroup() %>%
    mutate(Trait = trait_name)
}

summarise_effect_slopes <- function(df) {
  df %>%
    group_by(Trait, Dose) %>%
    summarise(
      intercept_median = median(intercept, na.rm = TRUE),
      intercept_q2.5   = quantile(intercept, 0.025, na.rm = TRUE),
      intercept_q97.5  = quantile(intercept, 0.975, na.rm = TRUE),
      slope_median     = median(slope, na.rm = TRUE),
      slope_q2.5       = quantile(slope, 0.025, na.rm = TRUE),
      slope_q97.5      = quantile(slope, 0.975, na.rm = TRUE),
      p_slope_gt0      = mean(slope > 0, na.rm = TRUE),
      p_slope_lt0      = mean(slope < 0, na.rm = TRUE),
      slope_direction  = case_when(
        slope_q2.5 > 0  ~ "Positive",
        slope_q97.5 < 0 ~ "Negative",
        TRUE            ~ "Uncertain"
      ),
      .groups = "drop"
    )
}

# ===== Read oxygen data =======================================================
colspec <- list(
  T         = readr::col_double(),
  Dose      = readr::col_character(),
  Replicate = readr::col_character(),
  Time      = readr::col_double(),
  Oxygen    = readr::col_double()
)

tmp_cols <- try(readr::read_csv(IN_CSV, n_max = 1, show_col_types = FALSE), silent = TRUE)
has_o2fit <- !inherits(tmp_cols, "try-error") && "O2_fit" %in% names(tmp_cols)
if (has_o2fit) colspec$O2_fit <- readr::col_double()

o2f_raw <- readr::read_csv(IN_CSV, col_types = do.call(readr::cols, colspec)) %>%
  mutate(
    T         = as.numeric(T),
    Dose      = norm_dose(Dose),
    Replicate = toupper(as.character(Replicate)),
    series_id = paste0("T=", T, " | Dose=", Dose, " | Rep=", Replicate),
    Oxygen_used = if ("O2_fit" %in% names(.)) dplyr::coalesce(O2_fit, Oxygen) else Oxygen
  ) %>%
  arrange(T, Dose, Replicate, Time)

stopifnot(all(c("T", "Dose", "Replicate", "Time", "Oxygen_used") %in% names(o2f_raw)))

if (is.null(allowed_doses)) {
  allowed_doses <- dose_levels_from_data(o2f_raw$Dose)
} else {
  allowed_doses <- norm_dose(allowed_doses)
}

o2f <- o2f_raw %>% filter(Dose %in% allowed_doses)
if (nrow(o2f) == 0) stop("No rows remain after filtering by allowed_doses.")

dose_levels <- dose_levels_from_data(o2f$Dose)
dose_levels <- c(intersect("Control", dose_levels), setdiff(dose_levels, "Control"))
dose_cols   <- make_dose_colors(dose_levels)

dose_key_tbl <- tibble(
  Dose = dose_levels,
  Dose_key = paste0("D", seq_along(dose_levels))
)
readr::write_csv(dose_key_tbl, dose_key_csv)

ymin_all <- suppressWarnings(min(o2f$Oxygen_used, na.rm = TRUE))
ymax_all <- suppressWarnings(max(o2f$Oxygen_used, na.rm = TRUE))
pad_all  <- 0.02 * (ymax_all - ymin_all)
if (!is.finite(pad_all)) pad_all <- 0.1
Y_LIMITS_SERIES <- c(ymin_all - pad_all, ymax_all + pad_all)

# ===== Read trimming metadata =================================================
trim_meta <- readr::read_csv(TRIM_META_CSV, show_col_types = FALSE) %>%
  mutate(
    T = as.numeric(T),
    Dose = norm_dose(Dose),
    Replicate = toupper(as.character(Replicate)),
    delta_Ninoc_to_N0_min = as.numeric(delta_Ninoc_to_N0_min)
  ) %>%
  select(T, Dose, Replicate, delta_Ninoc_to_N0_min) %>%
  distinct()

stopifnot(all(c("T", "Dose", "Replicate", "delta_Ninoc_to_N0_min") %in% names(trim_meta)))

group_lookup <- o2f %>%
  distinct(T, Dose, Replicate) %>%
  left_join(trim_meta, by = c("T", "Dose", "Replicate")) %>%
  mutate(N_inoculation_cells_per_L = as.numeric(N_inoculation_cells_per_L)) %>%
  arrange(T, Dose, Replicate)

readr::write_csv(group_lookup, group_lookup_csv)

if (any(is.na(group_lookup$delta_Ninoc_to_N0_min))) {
  warning("Missing delta_Ninoc_to_N0_min for some groups in trimming metadata. N0 and R will be NA for those groups.")
}

# ===== Fit one oxygen series ==================================================
fit_one <- function(df, y_limits = NULL, rmse_keep_threshold = 0.06) {
  df0 <- df %>%
    arrange(Time) %>%
    mutate(Time0 = Time - min(Time, na.rm = TRUE))
  
  y <- df0$Oxygen_used
  
  base_plot <- function(subtitle_txt) {
    p <- ggplot(df0, aes(Time0, Oxygen_used)) +
      geom_point(size = 1.3) +
      labs(
        title = df$series_id[1],
        subtitle = subtitle_txt,
        x = "Time (min, rebased)",
        y = "O₂ (mg/L)"
      ) +
      theme_classic(12)
    
    if (!is.null(y_limits) && all(is.finite(y_limits))) {
      p <- p + coord_cartesian(ylim = y_limits)
    }
    p
  }
  
  if (nrow(df0) < 6 || any(!is.finite(y))) {
    return(list(
      coefs = tibble(
        parameter = c("O2_0", "r", "K"),
        Estimate  = c(NA_real_, NA_real_, NA_real_),
        SE        = NA_real_,
        p_value   = NA_real_
      ),
      metrics = tibble(
        T = df$T[1], Dose = df$Dose[1], Replicate = df$Replicate[1],
        n = nrow(df0), r2 = NA_real_, rmse = NA_real_, rss = NA_real_,
        aic = NA_real_, aicc = NA_real_, T_end_min = NA_real_, keep = FALSE
      ),
      keep = FALSE,
      plot = base_plot("Too few/invalid points")
    ))
  }
  
  O2_0_start <- get_baseline(y)
  if (!is.finite(O2_0_start)) O2_0_start <- y[1]
  
  seg_n  <- max(6, floor(0.25 * nrow(df0)))
  seg    <- head(df0, seg_n)
  slope0 <- suppressWarnings(median(diff(seg$Oxygen_used) / diff(seg$Time0), na.rm = TRUE))
  if (!is.finite(slope0)) slope0 <- -1e-5
  
  r_start <- 1e-3
  K_est   <- abs(slope0)
  if (!is.finite(K_est) || K_est <= 0) K_est <- 1e-6
  K_start <- pmin(pmax(K_est, 1e-8), 1e-2)
  
  r_lower <- 1e-6
  r_upper <- 0.05
  K_lower <- max(1e-10, K_est / 50)
  K_upper <- min(1e-2,  K_est * 50)
  
  yr <- range(y[is.finite(y)], na.rm = TRUE)
  yspan <- diff(yr)
  if (!is.finite(yspan) || yspan <= 0) yspan <- abs(y[1]) + 1
  O2_lower <- yr[1] - 0.5 * yspan
  O2_upper <- yr[2] + 0.5 * yspan
  
  fit <- try(
    nlsLM(
      Oxygen_used ~ resp_model(r, K, Time0, O2_0),
      data = df0,
      start = list(r = r_start, K = K_start, O2_0 = O2_0_start),
      lower = c(r = r_lower, K = K_lower, O2_0 = O2_lower),
      upper = c(r = r_upper, K = K_upper, O2_0 = O2_upper),
      control = nls.lm.control(maxiter = 1500, ftol = 1e-12, ptol = 1e-12)
    ),
    silent = TRUE
  )
  
  if (inherits(fit, "try-error")) {
    return(list(
      coefs = tibble(
        parameter = c("O2_0", "r", "K"),
        Estimate  = c(NA_real_, NA_real_, NA_real_),
        SE        = NA_real_,
        p_value   = NA_real_
      ),
      metrics = tibble(
        T = df$T[1], Dose = df$Dose[1], Replicate = df$Replicate[1],
        n = nrow(df0), r2 = NA_real_, rmse = NA_real_, rss = NA_real_,
        aic = NA_real_, aicc = NA_real_, T_end_min = NA_real_, keep = FALSE
      ),
      keep = FALSE,
      plot = base_plot("Fit failed")
    ))
  }
  
  preds <- as.numeric(predict(fit, df0))
  n <- nrow(df0)
  
  rss <- sum((df0$Oxygen_used - preds)^2, na.rm = TRUE)
  rmse <- sqrt(rss / n)
  
  ss_tot <- sum((df0$Oxygen_used - mean(df0$Oxygen_used, na.rm = TRUE))^2, na.rm = TRUE)
  r2 <- if (is.finite(ss_tot) && ss_tot > 1e-12) 1 - rss / ss_tot else NA_real_
  
  k_param <- 3L
  aic <- n * log(rss / n) + 2 * k_param
  aicc <- if (n > k_param + 1) aic + (2 * k_param * (k_param + 1)) / (n - k_param - 1) else NA_real_
  
  T_end_min <- suppressWarnings(max(df0$Time0, na.rm = TRUE))
  
  co <- coef(fit)
  low_vec <- c(r = r_lower, K = K_lower, O2_0 = O2_lower)
  up_vec  <- c(r = r_upper, K = K_upper, O2_0 = O2_upper)
  on_boundary <- any(abs(co - low_vec[names(co)]) < 1e-10 | abs(co - up_vec[names(co)]) < 1e-10, na.rm = TRUE)
  
  keep <- (rmse < rmse_keep_threshold) && (n >= 6) && !on_boundary
  
  co_sum <- as.data.frame(summary(fit)$parameters) %>%
    tibble::rownames_to_column("parameter") %>%
    as_tibble()
  
  nm <- names(co_sum)
  if ("Std. Error" %in% nm) nm[nm == "Std. Error"] <- "SE"
  if ("Pr(>|t|)"  %in% nm) nm[nm == "Pr(>|t|)"]    <- "p_value"
  names(co_sum) <- nm
  
  if (!"SE" %in% names(co_sum)) co_sum$SE <- NA_real_
  if (!"p_value" %in% names(co_sum)) co_sum$p_value <- NA_real_
  
  metrics <- tibble(
    T = df$T[1], Dose = df$Dose[1], Replicate = df$Replicate[1],
    n = n, r2 = r2, rmse = rmse, rss = rss,
    aic = aic, aicc = aicc, T_end_min = T_end_min, keep = keep
  )
  
  subtitle_txt <- paste0(
    "R²=", ifelse(is.na(r2), "NA", sprintf("%.3f", r2)),
    " | RMSE=", sprintf("%.3g", rmse),
    " | AIC=", sprintf("%.1f", aic),
    if (!keep) " — FILTERED" else ""
  )
  
  p <- ggplot(df0, aes(Time0, Oxygen_used)) +
    geom_point(size = 1.3) +
    geom_line(aes(y = preds), linewidth = 0.9, color = "red") +
    labs(
      title = df$series_id[1],
      subtitle = subtitle_txt,
      x = "Time (min, rebased)",
      y = "O₂ (mg/L)"
    ) +
    theme_classic(12)
  
  if (!is.null(y_limits) && all(is.finite(y_limits))) {
    p <- p + coord_cartesian(ylim = y_limits)
  }
  
  list(coefs = co_sum, metrics = metrics, keep = keep, plot = p)
}

# ===== Run oxygen fits ========================================================
groups <- o2f %>%
  group_by(T, Dose, Replicate) %>%
  group_split()

all_coef_rows <- list()
all_metrics   <- list()

pdf(pdf_path, width = 6.8, height = 4.6)
for (g in groups) {
  res <- fit_one(g, y_limits = Y_LIMITS_SERIES, rmse_keep_threshold = RMSE_KEEP_THRESHOLD)
  all_metrics[[length(all_metrics) + 1]] <- res$metrics
  
  if (isTRUE(res$keep)) {
    print(res$plot)
    all_coef_rows[[length(all_coef_rows) + 1]] <-
      tibble(T = g$T[1], Dose = g$Dose[1], Replicate = g$Replicate[1]) %>%
      bind_cols(res$coefs) %>%
      mutate(T_end_min = res$metrics$T_end_min[1])
  }
}
dev.off()

coef_out <- bind_rows(all_coef_rows)
fit_metrics_out <- bind_rows(all_metrics)

readr::write_csv(coef_out, coef_csv)
readr::write_csv(fit_metrics_out, fit_metrics_csv)

# ===== Wide coefficients ======================================================
coef_wide <- coef_out %>%
  select(T, Dose, Replicate, parameter, Estimate, T_end_min) %>%
  tidyr::pivot_wider(names_from = parameter, values_from = Estimate) %>%
  group_by(T, Dose, Replicate) %>%
  summarise(
    r = first(r),
    K = first(K),
    O2_0 = first(O2_0),
    T_end_min = first(T_end_min),
    .groups = "drop"
  ) %>%
  arrange(T, Dose, Replicate)

readr::write_csv(coef_wide, coef_wide_csv)

# ===== Compute sample-specific N0, respiration, carbon units, CUE ============
results <- coef_wide %>%
  left_join(group_lookup, by = c("T", "Dose", "Replicate")) %>%
  mutate(
    N0_cells_per_L = dplyr::if_else(
      is.finite(N_inoculation_cells_per_L) & N_inoculation_cells_per_L > 0 &
        is.finite(delta_Ninoc_to_N0_min) & delta_Ninoc_to_N0_min >= 0 &
        is.finite(r) & r > 0,
      N_inoculation_cells_per_L * exp(r * delta_Ninoc_to_N0_min),
      NA_real_
    ),
    C_tot_O2_mg_per_L = dplyr::if_else(
      is.finite(K) & is.finite(r) & r > 0 &
        is.finite(T_end_min) & T_end_min > 0,
      (K / r) * (exp(r * T_end_min) - 1),
      NA_real_
    ),
    biomass_integral_cells_min_per_L = dplyr::if_else(
      is.finite(N0_cells_per_L) & N0_cells_per_L > 0 &
        is.finite(r) & r > 0 &
        is.finite(T_end_min) & T_end_min > 0,
      N0_cells_per_L * (exp(r * T_end_min) - 1) / r,
      NA_real_
    ),
    R_O2_mg_cell_min = dplyr::if_else(
      is.finite(C_tot_O2_mg_per_L) & C_tot_O2_mg_per_L > 0 &
        is.finite(biomass_integral_cells_min_per_L) & biomass_integral_cells_min_per_L > 0,
      C_tot_O2_mg_per_L / biomass_integral_cells_min_per_L,
      NA_real_
    ),
    cell_volume_um3 = CELL_VOLUME_UM3,
    cell_carbon_fg = CELL_CARBON_FG_PER_CELL,
    growth_fgC_h = dplyr::if_else(
      is.finite(r) & r > 0,
      r * cell_carbon_fg * MIN_TO_H,
      NA_real_
    ),
    respiration_fgC_h = dplyr::if_else(
      is.finite(R_O2_mg_cell_min) & R_O2_mg_cell_min > 0,
      R_O2_mg_cell_min * MG_TO_FG * RESPIRATORY_QUOTIENT * MIN_TO_H,
      NA_real_
    ),
    growth_C_per_C_h = dplyr::if_else(
      is.finite(growth_fgC_h) & growth_fgC_h > 0 &
        is.finite(cell_carbon_fg) & cell_carbon_fg > 0,
      growth_fgC_h / cell_carbon_fg,
      NA_real_
    ),
    respiration_C_per_C_h = dplyr::if_else(
      is.finite(respiration_fgC_h) & respiration_fgC_h > 0 &
        is.finite(cell_carbon_fg) & cell_carbon_fg > 0,
      respiration_fgC_h / cell_carbon_fg,
      NA_real_
    ),
    CUE = dplyr::if_else(
      is.finite(growth_fgC_h) & growth_fgC_h > 0 &
        is.finite(respiration_fgC_h) & respiration_fgC_h > 0,
      growth_fgC_h / (growth_fgC_h + respiration_fgC_h),
      NA_real_
    ),
    resp_over_growth = dplyr::if_else(
      is.finite(respiration_fgC_h) & respiration_fgC_h > 0 &
        is.finite(growth_fgC_h) & growth_fgC_h > 0,
      respiration_fgC_h / growth_fgC_h,
      NA_real_
    )
  ) %>%
  arrange(T, Dose, Replicate)

# ===== Exclude specific outlier points =======================================
EXCLUDE_POINTS <- tibble(
  T         = c(26,     27,     28,     28,    26,    28, 28, 28),
  Dose      = c("0.06", "0.06", "0.06", "1",   "0.5", "0.25", "0.06", "0.12"),
  Replicate = c("R3",   "R2",   "R1",   "R1",  "R3",  "R1", "R2", "R3")
)

results <- results %>%
  anti_join(EXCLUDE_POINTS, by = c("T", "Dose", "Replicate"))

message(sprintf(
  "Exclusion filter applied: %d point(s) removed. %d rows retained.",
  nrow(EXCLUDE_POINTS), nrow(results)
))
# =============================================================================

readr::write_csv(results, derived_csv)

# ===== Descriptive plots ======================================================
results_plot <- results %>%
  filter(Dose %in% allowed_doses) %>%
  mutate(Dose = factor(Dose, levels = dose_levels))

# ===== Replication summary ====================================================
replication_summary <- results_plot %>%
  count(T, Dose, name = "n_replicates") %>%
  arrange(Dose, T)

readr::write_csv(replication_summary, replication_summary_csv)

p_box_growth <- ggplot(results_plot %>% filter(is.finite(growth_fgC_h)), aes(Dose, growth_fgC_h, fill = Dose)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.8) +
  facet_wrap(~T, scales = "free_y") +
  scale_fill_manual(values = dose_cols, limits = dose_levels, drop = FALSE) +
  labs(title = "Growth in carbon units", x = NULL, y = "Growth (fg C h^-1)") +
  theme_classic(12) +
  theme(legend.position = "none")

p_box_resp <- ggplot(results_plot %>% filter(is.finite(respiration_fgC_h)), aes(Dose, respiration_fgC_h, fill = Dose)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.8) +
  facet_wrap(~T, scales = "free_y") +
  scale_fill_manual(values = dose_cols, limits = dose_levels, drop = FALSE) +
  labs(title = "Respiration in carbon units", x = NULL, y = "Respiration (fg C h^-1)") +
  theme_classic(12) +
  theme(legend.position = "none")

ratio_dat <- results_plot %>%
  filter(
    is.finite(growth_fgC_h), growth_fgC_h > 0,
    is.finite(respiration_fgC_h), respiration_fgC_h > 0
  ) %>%
  mutate(
    log_resp_over_growth = log(resp_over_growth),
    TK = T + 273.15
  )

p_rg_t <- ggplot(ratio_dat, aes(T, resp_over_growth, color = Dose)) +
  geom_point(size = 2.1, alpha = 0.9) +
  scale_color_manual(values = dose_cols, limits = dose_levels, drop = FALSE) +
  labs(title = "Respiration / Growth vs Temperature", x = "Temperature (°C)", y = "Respiration / Growth", color = "Dose") +
  theme_classic(12)

p_growth_vs_T <- ggplot(results_plot %>% filter(is.finite(growth_fgC_h), growth_fgC_h > 0), aes(T, growth_fgC_h, color = Dose)) +
  geom_point(size = 2.1, alpha = 0.9) +
  scale_color_manual(values = dose_cols, limits = dose_levels, drop = FALSE) +
  labs(title = "Growth vs Temperature", x = "Temperature (°C)", y = "Growth (fg C h^-1)", color = "Dose") +
  theme_classic(12)

p_resp_vs_T <- ggplot(results_plot %>% filter(is.finite(respiration_fgC_h), respiration_fgC_h > 0), aes(T, respiration_fgC_h, color = Dose)) +
  geom_point(size = 2.1, alpha = 0.9) +
  scale_color_manual(values = dose_cols, limits = dose_levels, drop = FALSE) +
  labs(title = "Respiration vs Temperature", x = "Temperature (°C)", y = "Respiration (fg C h^-1)", color = "Dose") +
  theme_classic(12)

p_box_growth_biomass <- ggplot(
  results_plot %>% filter(is.finite(growth_C_per_C_h), growth_C_per_C_h > 0),
  aes(Dose, growth_C_per_C_h, fill = Dose)
) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.8) +
  facet_wrap(~T, scales = "free_y") +
  scale_fill_manual(values = dose_cols, limits = dose_levels, drop = FALSE) +
  labs(title = "Biomass-corrected growth", x = NULL, y = "Growth (C per C per h)") +
  theme_classic(12) +
  theme(legend.position = "none")

p_box_resp_biomass <- ggplot(
  results_plot %>% filter(is.finite(respiration_C_per_C_h), respiration_C_per_C_h > 0),
  aes(Dose, respiration_C_per_C_h, fill = Dose)
) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.8) +
  facet_wrap(~T, scales = "free_y") +
  scale_fill_manual(values = dose_cols, limits = dose_levels, drop = FALSE) +
  labs(title = "Biomass-corrected respiration", x = NULL, y = "Respiration (C per C per h)") +
  theme_classic(12) +
  theme(legend.position = "none")

p_growth_biomass_vs_T <- ggplot(
  results_plot %>% filter(is.finite(growth_C_per_C_h), growth_C_per_C_h > 0),
  aes(T, growth_C_per_C_h, color = Dose)
) +
  geom_point(size = 2.1, alpha = 0.9) +
  scale_color_manual(values = dose_cols, limits = dose_levels, drop = FALSE) +
  labs(title = "Biomass-corrected growth vs Temperature", x = "Temperature (°C)", y = "Growth (C per C per h)", color = "Dose") +
  theme_classic(12)

p_resp_biomass_vs_T <- ggplot(
  results_plot %>% filter(is.finite(respiration_C_per_C_h), respiration_C_per_C_h > 0),
  aes(T, respiration_C_per_C_h, color = Dose)
) +
  geom_point(size = 2.1, alpha = 0.9) +
  scale_color_manual(values = dose_cols, limits = dose_levels, drop = FALSE) +
  labs(title = "Biomass-corrected respiration vs Temperature", x = "Temperature (°C)", y = "Respiration (C per C per h)", color = "Dose") +
  theme_classic(12)

# ===== Prepare data for Bayesian thermal models ===============================
growth_dat <- results_plot %>%
  filter(is.finite(growth_fgC_h), growth_fgC_h > 0) %>%
  mutate(Dose = as.character(Dose)) %>%
  left_join(dose_key_tbl, by = "Dose") %>%
  transmute(
    T = as.numeric(T),
    TK = T + 273.15,
    boltz_shift = (1 / (k_B * T_ref)) - (1 / (k_B * (T + 273.15))),
    Dose = factor(Dose, levels = dose_levels),
    Dose_key = factor(Dose_key, levels = dose_key_tbl$Dose_key),
    Replicate = factor(Replicate),
    cond_id = factor(interaction(Dose, T, drop = TRUE)),
    y = log(as.numeric(growth_fgC_h)),
    y_raw = as.numeric(growth_fgC_h)
  )

resp_dat <- results_plot %>%
  filter(is.finite(respiration_fgC_h), respiration_fgC_h > 0) %>%
  mutate(Dose = as.character(Dose)) %>%
  left_join(dose_key_tbl, by = "Dose") %>%
  transmute(
    T = as.numeric(T),
    TK = T + 273.15,
    boltz_shift = (1 / (k_B * T_ref)) - (1 / (k_B * (T + 273.15))),
    Dose = factor(Dose, levels = dose_levels),
    Dose_key = factor(Dose_key, levels = dose_key_tbl$Dose_key),
    Replicate = factor(Replicate),
    cond_id = factor(interaction(Dose, T, drop = TRUE)),
    y = log(as.numeric(respiration_fgC_h)),
    y_raw = as.numeric(respiration_fgC_h)
  )

growth_biomass_dat <- results_plot %>%
  filter(is.finite(growth_C_per_C_h), growth_C_per_C_h > 0) %>%
  mutate(Dose = as.character(Dose)) %>%
  left_join(dose_key_tbl, by = "Dose") %>%
  transmute(
    T = as.numeric(T),
    TK = T + 273.15,
    boltz_shift = (1 / (k_B * T_ref)) - (1 / (k_B * (T + 273.15))),
    Dose = factor(Dose, levels = dose_levels),
    Dose_key = factor(Dose_key, levels = dose_key_tbl$Dose_key),
    Replicate = factor(Replicate),
    cond_id = factor(interaction(Dose, T, drop = TRUE)),
    y = log(as.numeric(growth_C_per_C_h)),
    y_raw = as.numeric(growth_C_per_C_h)
  )

resp_biomass_dat <- results_plot %>%
  filter(is.finite(respiration_C_per_C_h), respiration_C_per_C_h > 0) %>%
  mutate(Dose = as.character(Dose)) %>%
  left_join(dose_key_tbl, by = "Dose") %>%
  transmute(
    T = as.numeric(T),
    TK = T + 273.15,
    boltz_shift = (1 / (k_B * T_ref)) - (1 / (k_B * (T + 273.15))),
    Dose = factor(Dose, levels = dose_levels),
    Dose_key = factor(Dose_key, levels = dose_key_tbl$Dose_key),
    Replicate = factor(Replicate),
    cond_id = factor(interaction(Dose, T, drop = TRUE)),
    y = log(as.numeric(respiration_C_per_C_h)),
    y_raw = as.numeric(respiration_C_per_C_h)
  )

ratio_dat_arr <- ratio_dat %>%
  mutate(Dose = as.character(Dose)) %>%
  left_join(dose_key_tbl, by = "Dose") %>%
  filter(is.finite(log_resp_over_growth)) %>%
  transmute(
    T = as.numeric(T),
    TK = as.numeric(TK),
    boltz_shift = (1 / (k_B * T_ref)) - (1 / (k_B * TK)),
    Dose = factor(Dose, levels = dose_levels),
    Dose_key = factor(Dose_key, levels = dose_key_tbl$Dose_key),
    Replicate = factor(Replicate),
    cond_id = factor(interaction(Dose, T, drop = TRUE)),
    y = as.numeric(log_resp_over_growth),
    y_raw = as.numeric(resp_over_growth)
  )

cue_dat <- ratio_dat %>%
  mutate(Dose = as.character(Dose)) %>%
  left_join(dose_key_tbl, by = "Dose") %>%
  filter(is.finite(CUE), CUE > 0) %>%
  transmute(
    T = as.numeric(T),
    TK = as.numeric(TK),
    boltz_shift = (1 / (k_B * T_ref)) - (1 / (k_B * TK)),
    Dose = factor(Dose, levels = dose_levels),
    Dose_key = factor(Dose_key, levels = dose_key_tbl$Dose_key),
    Replicate = factor(Replicate),
    cond_id = factor(interaction(Dose, T, drop = TRUE)),
    y = log(as.numeric(CUE)),
    y_raw = as.numeric(CUE)
  )

# ===== Descriptive CUE plot on Boltzmann x-axis ===============================
p_cue_boltz_scatter <- ggplot(cue_dat, aes(boltz_shift, y_raw, color = Dose)) +
  geom_point(size = 2.1, alpha = 0.9) +
  scale_y_log10() +
  scale_color_manual(values = dose_cols, limits = dose_levels, drop = FALSE) +
  labs(
    title = "CUE vs Boltzmann temperature scale",
    x = expression(frac(1, k[T[ref]]) - frac(1, k*T)),
    y = "CUE (log scale)",
    color = "Dose"
  ) +
  theme_classic(12)

# ===== Bayesian Arrhenius model ==============================================
fit_bayes_arr_by_dose <- function(dat) {
  if (nrow(dat) < 8) stop("Too few rows for Bayesian Arrhenius fit.")
  if (dplyr::n_distinct(dat$Dose_key) < 1) stop("No dose levels found.")
  if (dplyr::n_distinct(dat$TK) < 3) stop("Too few distinct temperatures for Arrhenius fit.")
  
  bf_arr <- brms::bf(
    y ~ alpha + E * boltz_shift,
    alpha ~ 0 + Dose_key,
    E     ~ 0 + Dose_key,
    nl = TRUE
  )
  
  pri <- c(
    brms::prior("normal(0, 5)", nlpar = "alpha"),
    brms::prior("normal(0.65, 0.35)", nlpar = "E", lb = 0, ub = 3),
    brms::prior("exponential(1)", class = "sigma")
  )
  
  brms::brm(
    formula = bf_arr,
    data = dat,
    family = gaussian(),
    prior = pri,
    iter = BAYES_ITER,
    warmup = BAYES_WARMUP,
    chains = BAYES_CHAINS,
    seed = BAYES_SEED,
    control = list(adapt_delta = BAYES_ADAPT, max_treedepth = BAYES_MAX_TD),
    backend = "rstan",
    refresh = 0,
    init = 0
  )
}

fit_bayes_arr_cue_common_slope <- function(dat) {
  if (nrow(dat) < 8) stop("Too few rows for Bayesian CUE Arrhenius fit.")
  if (dplyr::n_distinct(dat$Dose_key) < 1) stop("No dose levels found.")
  if (dplyr::n_distinct(dat$TK) < 3) stop("Too few distinct temperatures for CUE Arrhenius fit.")
  
  bf_arr_cue <- brms::bf(
    y ~ alpha + E * boltz_shift,
    alpha ~ 0 + Dose_key,
    E ~ 1,
    nl = TRUE
  )
  
  pri <- c(
    brms::prior("normal(0, 5)", nlpar = "alpha"),
    brms::prior("normal(0, 1)", nlpar = "E"),
    brms::prior("exponential(1)", class = "sigma")
  )
  
  brms::brm(
    formula = bf_arr_cue,
    data = dat,
    family = gaussian(),
    prior = pri,
    iter = BAYES_ITER,
    warmup = BAYES_WARMUP,
    chains = BAYES_CHAINS,
    seed = BAYES_SEED,
    control = list(adapt_delta = BAYES_ADAPT, max_treedepth = BAYES_MAX_TD),
    backend = "rstan",
    refresh = 0,
    init = 0
  )
}

# ===== Separate hierarchical Arrhenius sensitivity analysis ===================
fit_bayes_arr_by_dose_hier <- function(dat) {
  if (nrow(dat) < 8) stop("Too few rows for hierarchical Bayesian Arrhenius fit.")
  if (dplyr::n_distinct(dat$Dose_key) < 1) stop("No dose levels found.")
  if (dplyr::n_distinct(dat$TK) < 3) stop("Too few distinct temperatures for Arrhenius fit.")
  if (!"cond_id" %in% names(dat)) stop("cond_id column not found.")
  
  bf_arr_hier <- brms::bf(
    y ~ alpha + E * boltz_shift,
    alpha ~ 0 + Dose_key + (1 | cond_id),
    E     ~ 0 + Dose_key,
    nl = TRUE
  )
  
  pri <- c(
    brms::prior("normal(0, 5)", nlpar = "alpha"),
    brms::prior("normal(0.65, 0.35)", nlpar = "E", lb = 0, ub = 3),
    brms::prior("student_t(3, 0, 2.5)", class = "sd", nlpar = "alpha", group = "cond_id"),
    brms::prior("exponential(1)", class = "sigma")
  )
  
  brms::brm(
    formula = bf_arr_hier,
    data = dat,
    family = gaussian(),
    prior = pri,
    iter = BAYES_ITER,
    warmup = BAYES_WARMUP,
    chains = BAYES_CHAINS,
    seed = BAYES_SEED,
    control = list(adapt_delta = 0.99, max_treedepth = BAYES_MAX_TD),
    backend = "rstan",
    refresh = 0,
    init = 0
  )
}

# ===== Bayesian full TPC: Sharpe–Schoolfield =================================
fit_bayes_ss_tpc_by_dose <- function(dat) {
  if (nrow(dat) < 8) stop("Too few rows for Bayesian Sharpe–Schoolfield TPC fit.")
  if (dplyr::n_distinct(dat$Dose_key) < 1) stop("No dose levels found.")
  if (dplyr::n_distinct(dat$T) < 4) stop("Too few distinct temperatures for TPC fit.")
  
  max_TK_obs <- max(dat$TK, na.rm = TRUE)
  
  lnB0_prior_mean <- log(stats::median(dat$y_raw, na.rm = TRUE))
  Th_prior_mean <- max_TK_obs + 1.5
  
  prior_lnB0 <- sprintf("normal(%0.8f, 2)", lnB0_prior_mean)
  prior_Th   <- sprintf("normal(%0.8f, 3)", Th_prior_mean)
  
  bf_ss <- brms::bf(
    y ~ lnB0 -
      (E * 11604.51812) * ((1 / TK) - (1 / 293.15)) -
      log(1 + exp((Eh * 11604.51812) * ((1 / Th) - (1 / TK)))),
    lnB0 ~ 0 + Dose_key,
    E    ~ 0 + Dose_key,
    Eh   ~ 0 + Dose_key,
    Th   ~ 0 + Dose_key,
    nl = TRUE
  )
  
  pri <- c(
    do.call(brms::prior, list(prior = prior_lnB0, nlpar = "lnB0")),
    brms::prior("normal(0.65, 0.35)", nlpar = "E",  lb = 0, ub = 3),
    brms::prior("normal(4, 1.5)",     nlpar = "Eh", lb = 0, ub = 10),
    do.call(brms::prior, list(prior = prior_Th, nlpar = "Th", lb = max_TK_obs - 1, ub = 340)),
    brms::prior("exponential(1)", class = "sigma")
  )
  
  brms::brm(
    formula = bf_ss,
    data = dat,
    family = gaussian(),
    prior = pri,
    iter = BAYES_ITER,
    warmup = BAYES_WARMUP,
    chains = BAYES_CHAINS,
    seed = BAYES_SEED,
    control = list(adapt_delta = 0.99, max_treedepth = BAYES_MAX_TD),
    backend = "rstan",
    refresh = 0
  )
}

# ===== Fit Bayesian thermal models ===========================================
fit_growth_ok <- nrow(growth_dat) >= 8 && dplyr::n_distinct(growth_dat$TK) >= 3 && dplyr::n_distinct(growth_dat$Dose_key) >= 1
fit_resp_ok   <- nrow(resp_dat)   >= 8 && dplyr::n_distinct(resp_dat$TK)   >= 3 && dplyr::n_distinct(resp_dat$Dose_key)   >= 1
fit_ratio_ok  <- nrow(ratio_dat_arr) >= 8 && dplyr::n_distinct(ratio_dat_arr$TK) >= 3 && dplyr::n_distinct(ratio_dat_arr$Dose_key) >= 1
fit_cue_ok    <- nrow(cue_dat) >= 8 && dplyr::n_distinct(cue_dat$TK) >= 3 && dplyr::n_distinct(cue_dat$Dose_key) >= 1

fit_growth_biomass_ok <- nrow(growth_biomass_dat) >= 8 && dplyr::n_distinct(growth_biomass_dat$TK) >= 3 && dplyr::n_distinct(growth_biomass_dat$Dose_key) >= 1
fit_resp_biomass_ok   <- nrow(resp_biomass_dat)   >= 8 && dplyr::n_distinct(resp_biomass_dat$TK)   >= 3 && dplyr::n_distinct(resp_biomass_dat$Dose_key)   >= 1

fit_growth_brm <- if (fit_growth_ok) fit_bayes_arr_by_dose(growth_dat) else NULL
fit_resp_brm   <- if (fit_resp_ok)   fit_bayes_arr_by_dose(resp_dat) else NULL
fit_ratio_brm  <- if (fit_ratio_ok)  fit_bayes_arr_by_dose(ratio_dat_arr) else NULL
fit_cue_brm    <- if (fit_cue_ok)    fit_bayes_arr_cue_common_slope(cue_dat) else NULL

fit_growth_biomass_brm <- if (fit_growth_biomass_ok) fit_bayes_arr_by_dose(growth_biomass_dat) else NULL
fit_resp_biomass_brm   <- if (fit_resp_biomass_ok)   fit_bayes_arr_by_dose(resp_biomass_dat) else NULL

fit_growth_tpc_brm <- if (fit_growth_ok) fit_bayes_ss_tpc_by_dose(growth_dat) else NULL
fit_resp_tpc_brm   <- if (fit_resp_ok)   fit_bayes_ss_tpc_by_dose(resp_dat) else NULL

fit_growth_tpc_biomass_brm <- if (fit_growth_biomass_ok) fit_bayes_ss_tpc_by_dose(growth_biomass_dat) else NULL
fit_resp_tpc_biomass_brm   <- if (fit_resp_biomass_ok)   fit_bayes_ss_tpc_by_dose(resp_biomass_dat) else NULL

fit_growth_hier_brm <- if (fit_growth_ok) fit_bayes_arr_by_dose_hier(growth_dat) else NULL
fit_resp_hier_brm   <- if (fit_resp_ok)   fit_bayes_arr_by_dose_hier(resp_dat) else NULL
fit_ratio_hier_brm  <- if (fit_ratio_ok)  fit_bayes_arr_by_dose_hier(ratio_dat_arr) else NULL

fit_growth_biomass_hier_brm <- if (fit_growth_biomass_ok) fit_bayes_arr_by_dose_hier(growth_biomass_dat) else NULL
fit_resp_biomass_hier_brm   <- if (fit_resp_biomass_ok)   fit_bayes_arr_by_dose_hier(resp_biomass_dat) else NULL

if (!is.null(fit_growth_brm)) saveRDS(fit_growth_brm, brms_growth_rds)
if (!is.null(fit_resp_brm)) saveRDS(fit_resp_brm, brms_resp_rds)
if (!is.null(fit_ratio_brm)) saveRDS(fit_ratio_brm, brms_ratio_rds)
if (!is.null(fit_cue_brm)) saveRDS(fit_cue_brm, brms_cue_rds)
if (!is.null(fit_growth_biomass_brm)) saveRDS(fit_growth_biomass_brm, brms_growth_biomass_rds)
if (!is.null(fit_resp_biomass_brm)) saveRDS(fit_resp_biomass_brm, brms_resp_biomass_rds)

if (!is.null(fit_growth_tpc_brm)) saveRDS(fit_growth_tpc_brm, brms_tpc_growth_rds)
if (!is.null(fit_resp_tpc_brm)) saveRDS(fit_resp_tpc_brm, brms_tpc_resp_rds)
if (!is.null(fit_growth_tpc_biomass_brm)) saveRDS(fit_growth_tpc_biomass_brm, brms_tpc_growth_biomass_rds)
if (!is.null(fit_resp_tpc_biomass_brm)) saveRDS(fit_resp_tpc_biomass_brm, brms_tpc_resp_biomass_rds)

if (!is.null(fit_growth_hier_brm)) saveRDS(fit_growth_hier_brm, brms_growth_hier_rds)
if (!is.null(fit_resp_hier_brm)) saveRDS(fit_resp_hier_brm, brms_resp_hier_rds)
if (!is.null(fit_ratio_hier_brm)) saveRDS(fit_ratio_hier_brm, brms_ratio_hier_rds)
if (!is.null(fit_growth_biomass_hier_brm)) saveRDS(fit_growth_biomass_hier_brm, brms_growth_biomass_hier_rds)
if (!is.null(fit_resp_biomass_hier_brm)) saveRDS(fit_resp_biomass_hier_brm, brms_resp_biomass_hier_rds)

# ===== Posterior summaries ====================================================
get_key_par_names_arr <- function(keys) {
  c(
    paste0("b_alpha_Dose_key", keys),
    paste0("b_E_Dose_key", keys)
  )
}

get_key_par_names_ss <- function(keys) {
  c(
    paste0("b_lnB0_Dose_key", keys),
    paste0("b_E_Dose_key", keys),
    paste0("b_Eh_Dose_key", keys),
    paste0("b_Th_Dose_key", keys)
  )
}

extract_arr_posterior <- function(fit, dose_key_tbl) {
  draws <- posterior::as_draws_df(fit) %>% as_tibble()
  
  purrr::map_dfr(seq_len(nrow(dose_key_tbl)), function(i) {
    d  <- dose_key_tbl$Dose[i]
    dk <- dose_key_tbl$Dose_key[i]
    tibble(
      Dose = d,
      Dose_key = dk,
      .draw = seq_len(nrow(draws)),
      alpha = draws[[paste0("b_alpha_Dose_key", dk)]],
      E     = draws[[paste0("b_E_Dose_key", dk)]],
      A     = exp(alpha)
    )
  })
}

extract_arr_cue_common_posterior <- function(fit, dose_key_tbl) {
  draws <- posterior::as_draws_df(fit) %>% as_tibble()
  
  purrr::map_dfr(seq_len(nrow(dose_key_tbl)), function(i) {
    d  <- dose_key_tbl$Dose[i]
    dk <- dose_key_tbl$Dose_key[i]
    tibble(
      Dose = d,
      Dose_key = dk,
      .draw = seq_len(nrow(draws)),
      alpha = draws[[paste0("b_alpha_Dose_key", dk)]],
      E     = draws[["b_E_Intercept"]],
      A     = exp(alpha)
    )
  })
}

extract_ss_posterior <- function(fit, dose_key_tbl) {
  draws <- posterior::as_draws_df(fit) %>% as_tibble()
  
  purrr::map_dfr(seq_len(nrow(dose_key_tbl)), function(i) {
    d  <- dose_key_tbl$Dose[i]
    dk <- dose_key_tbl$Dose_key[i]
    
    tibble(
      Dose = d,
      Dose_key = dk,
      .draw = seq_len(nrow(draws)),
      lnB0 = draws[[paste0("b_lnB0_Dose_key", dk)]],
      E    = draws[[paste0("b_E_Dose_key", dk)]],
      Eh   = draws[[paste0("b_Eh_Dose_key", dk)]],
      Th   = draws[[paste0("b_Th_Dose_key", dk)]],
      Th_C = draws[[paste0("b_Th_Dose_key", dk)]] - 273.15
    )
  })
}

if (!is.null(fit_growth_brm)) {
  pars_growth <- c(get_key_par_names_arr(dose_key_tbl$Dose_key), "sigma")
  summ_growth <- posterior_summary_df(fit_growth_brm, pars = pars_growth)
  readr::write_csv(summ_growth, summary_growth_csv)
  readr::write_csv(extract_arr_posterior(fit_growth_brm, dose_key_tbl), post_growth_csv)
}

if (!is.null(fit_resp_brm)) {
  pars_resp <- c(get_key_par_names_arr(dose_key_tbl$Dose_key), "sigma")
  summ_resp <- posterior_summary_df(fit_resp_brm, pars = pars_resp)
  readr::write_csv(summ_resp, summary_resp_csv)
  readr::write_csv(extract_arr_posterior(fit_resp_brm, dose_key_tbl), post_resp_csv)
}

if (!is.null(fit_ratio_brm)) {
  pars_ratio <- c(get_key_par_names_arr(dose_key_tbl$Dose_key), "sigma")
  summ_ratio <- posterior_summary_df(fit_ratio_brm, pars = pars_ratio)
  readr::write_csv(summ_ratio, summary_ratio_csv)
  readr::write_csv(extract_arr_posterior(fit_ratio_brm, dose_key_tbl), post_ratio_csv)
}

if (!is.null(fit_cue_brm)) {
  pars_cue <- c(paste0("b_alpha_Dose_key", dose_key_tbl$Dose_key), "b_E_Intercept", "sigma")
  summ_cue <- posterior_summary_df(fit_cue_brm, pars = pars_cue)
  readr::write_csv(summ_cue, summary_cue_csv)
  readr::write_csv(extract_arr_cue_common_posterior(fit_cue_brm, dose_key_tbl), post_cue_csv)
}

if (!is.null(fit_growth_biomass_brm)) {
  pars_growth_biomass <- c(get_key_par_names_arr(dose_key_tbl$Dose_key), "sigma")
  summ_growth_biomass <- posterior_summary_df(fit_growth_biomass_brm, pars = pars_growth_biomass)
  readr::write_csv(summ_growth_biomass, summary_growth_biomass_csv)
  readr::write_csv(extract_arr_posterior(fit_growth_biomass_brm, dose_key_tbl), post_growth_biomass_csv)
}

if (!is.null(fit_resp_biomass_brm)) {
  pars_resp_biomass <- c(get_key_par_names_arr(dose_key_tbl$Dose_key), "sigma")
  summ_resp_biomass <- posterior_summary_df(fit_resp_biomass_brm, pars = pars_resp_biomass)
  readr::write_csv(summ_resp_biomass, summary_resp_biomass_csv)
  readr::write_csv(extract_arr_posterior(fit_resp_biomass_brm, dose_key_tbl), post_resp_biomass_csv)
}

if (!is.null(fit_growth_tpc_brm)) {
  pars_tpc_g <- c(get_key_par_names_ss(dose_key_tbl$Dose_key), "sigma")
  summ_tpc_g <- posterior_summary_df(fit_growth_tpc_brm, pars = pars_tpc_g)
  readr::write_csv(summ_tpc_g, summary_tpc_growth_csv)
  readr::write_csv(extract_ss_posterior(fit_growth_tpc_brm, dose_key_tbl), post_tpc_growth_csv)
}

if (!is.null(fit_resp_tpc_brm)) {
  pars_tpc_r <- c(get_key_par_names_ss(dose_key_tbl$Dose_key), "sigma")
  summ_tpc_r <- posterior_summary_df(fit_resp_tpc_brm, pars = pars_tpc_r)
  readr::write_csv(summ_tpc_r, summary_tpc_resp_csv)
  readr::write_csv(extract_ss_posterior(fit_resp_tpc_brm, dose_key_tbl), post_tpc_resp_csv)
}

if (!is.null(fit_growth_tpc_biomass_brm)) {
  pars_tpc_g_bio <- c(get_key_par_names_ss(dose_key_tbl$Dose_key), "sigma")
  summ_tpc_g_bio <- posterior_summary_df(fit_growth_tpc_biomass_brm, pars = pars_tpc_g_bio)
  readr::write_csv(summ_tpc_g_bio, summary_tpc_growth_biomass_csv)
  readr::write_csv(extract_ss_posterior(fit_growth_tpc_biomass_brm, dose_key_tbl), post_tpc_growth_biomass_csv)
}

if (!is.null(fit_resp_tpc_biomass_brm)) {
  pars_tpc_r_bio <- c(get_key_par_names_ss(dose_key_tbl$Dose_key), "sigma")
  summ_tpc_r_bio <- posterior_summary_df(fit_resp_tpc_biomass_brm, pars = pars_tpc_r_bio)
  readr::write_csv(summ_tpc_r_bio, summary_tpc_resp_biomass_csv)
  readr::write_csv(extract_ss_posterior(fit_resp_tpc_biomass_brm, dose_key_tbl), post_tpc_resp_biomass_csv)
}

if (!is.null(fit_growth_hier_brm)) {
  pars <- c(get_key_par_names_arr(dose_key_tbl$Dose_key), "sigma")
  readr::write_csv(posterior_summary_df(fit_growth_hier_brm, pars = pars), summary_growth_hier_csv)
  readr::write_csv(extract_arr_posterior(fit_growth_hier_brm, dose_key_tbl), post_growth_hier_csv)
}

if (!is.null(fit_resp_hier_brm)) {
  pars <- c(get_key_par_names_arr(dose_key_tbl$Dose_key), "sigma")
  readr::write_csv(posterior_summary_df(fit_resp_hier_brm, pars = pars), summary_resp_hier_csv)
  readr::write_csv(extract_arr_posterior(fit_resp_hier_brm, dose_key_tbl), post_resp_hier_csv)
}

if (!is.null(fit_ratio_hier_brm)) {
  pars <- c(get_key_par_names_arr(dose_key_tbl$Dose_key), "sigma")
  readr::write_csv(posterior_summary_df(fit_ratio_hier_brm, pars = pars), summary_ratio_hier_csv)
  readr::write_csv(extract_arr_posterior(fit_ratio_hier_brm, dose_key_tbl), post_ratio_hier_csv)
}

if (!is.null(fit_growth_biomass_hier_brm)) {
  pars <- c(get_key_par_names_arr(dose_key_tbl$Dose_key), "sigma")
  readr::write_csv(posterior_summary_df(fit_growth_biomass_hier_brm, pars = pars), summary_growth_biomass_hier_csv)
  readr::write_csv(extract_arr_posterior(fit_growth_biomass_hier_brm, dose_key_tbl), post_growth_biomass_hier_csv)
}

if (!is.null(fit_resp_biomass_hier_brm)) {
  pars <- c(get_key_par_names_arr(dose_key_tbl$Dose_key), "sigma")
  readr::write_csv(posterior_summary_df(fit_resp_biomass_hier_brm, pars = pars), summary_resp_biomass_hier_csv)
  readr::write_csv(extract_arr_posterior(fit_resp_biomass_hier_brm, dose_key_tbl), post_resp_biomass_hier_csv)
}

# ===== Posterior prediction grids ============================================
make_arr_grid_from_posterior <- function(fit, dat, dose_key_tbl) {
  dr <- posterior::as_draws_df(fit) %>% as_tibble()
  
  grid_base <- dat %>%
    group_by(Dose, Dose_key) %>%
    summarise(T_min = min(T), T_max = max(T), .groups = "drop") %>%
    right_join(dose_key_tbl, by = c("Dose", "Dose_key")) %>%
    mutate(
      T_min = ifelse(is.na(T_min), min(dat$T, na.rm = TRUE), T_min),
      T_max = ifelse(is.na(T_max), max(dat$T, na.rm = TRUE), T_max)
    )
  
  purrr::map_dfr(seq_len(nrow(dose_key_tbl)), function(i) {
    d  <- dose_key_tbl$Dose[i]
    dk <- dose_key_tbl$Dose_key[i]
    
    T_seq <- seq(
      grid_base$T_min[grid_base$Dose_key == dk],
      grid_base$T_max[grid_base$Dose_key == dk],
      length.out = 200
    )
    TK_seq <- T_seq + 273.15
    boltz_shift <- (1 / (k_B * T_ref)) - (1 / (k_B * TK_seq))
    
    alpha <- dr[[paste0("b_alpha_Dose_key", dk)]]
    E     <- dr[[paste0("b_E_Dose_key", dk)]]
    
    pred_log_mat <- sapply(seq_along(boltz_shift), function(j) alpha + E * boltz_shift[j])
    pred_raw_mat <- exp(pred_log_mat)
    
    tibble(
      Dose = d,
      Dose_key = dk,
      T = T_seq,
      boltz_shift = boltz_shift,
      log_q2.5  = apply(pred_log_mat, 2, quantile, probs = 0.025, na.rm = TRUE),
      log_q50   = apply(pred_log_mat, 2, quantile, probs = 0.5,   na.rm = TRUE),
      log_q97.5 = apply(pred_log_mat, 2, quantile, probs = 0.975, na.rm = TRUE),
      raw_q2.5  = apply(pred_raw_mat, 2, quantile, probs = 0.025, na.rm = TRUE),
      raw_q50   = apply(pred_raw_mat, 2, quantile, probs = 0.5,   na.rm = TRUE),
      raw_q97.5 = apply(pred_raw_mat, 2, quantile, probs = 0.975, na.rm = TRUE)
    )
  })
}

make_arr_grid_from_posterior_common_slope <- function(fit, dat, dose_key_tbl) {
  dr <- posterior::as_draws_df(fit) %>% as_tibble()
  
  grid_base <- dat %>%
    group_by(Dose, Dose_key) %>%
    summarise(T_min = min(T), T_max = max(T), .groups = "drop") %>%
    right_join(dose_key_tbl, by = c("Dose", "Dose_key")) %>%
    mutate(
      T_min = ifelse(is.na(T_min), min(dat$T, na.rm = TRUE), T_min),
      T_max = ifelse(is.na(T_max), max(dat$T, na.rm = TRUE), T_max)
    )
  
  purrr::map_dfr(seq_len(nrow(dose_key_tbl)), function(i) {
    d  <- dose_key_tbl$Dose[i]
    dk <- dose_key_tbl$Dose_key[i]
    
    T_seq <- seq(
      grid_base$T_min[grid_base$Dose_key == dk],
      grid_base$T_max[grid_base$Dose_key == dk],
      length.out = 200
    )
    TK_seq <- T_seq + 273.15
    boltz_shift <- (1 / (k_B * T_ref)) - (1 / (k_B * TK_seq))
    
    alpha <- dr[[paste0("b_alpha_Dose_key", dk)]]
    E     <- dr[["b_E_Intercept"]]
    
    pred_log_mat <- sapply(seq_along(boltz_shift), function(j) alpha + E * boltz_shift[j])
    pred_raw_mat <- exp(pred_log_mat)
    
    tibble(
      Dose = d,
      Dose_key = dk,
      T = T_seq,
      boltz_shift = boltz_shift,
      log_q2.5  = apply(pred_log_mat, 2, quantile, probs = 0.025, na.rm = TRUE),
      log_q50   = apply(pred_log_mat, 2, quantile, probs = 0.5,   na.rm = TRUE),
      log_q97.5 = apply(pred_log_mat, 2, quantile, probs = 0.975, na.rm = TRUE),
      raw_q2.5  = apply(pred_raw_mat, 2, quantile, probs = 0.025, na.rm = TRUE),
      raw_q50   = apply(pred_raw_mat, 2, quantile, probs = 0.5,   na.rm = TRUE),
      raw_q97.5 = apply(pred_raw_mat, 2, quantile, probs = 0.975, na.rm = TRUE)
    )
  })
}

make_ss_grid_from_posterior <- function(fit, dat, dose_key_tbl) {
  dr <- posterior::as_draws_df(fit) %>% as_tibble()
  
  grid_base <- dat %>%
    group_by(Dose, Dose_key) %>%
    summarise(T_min = min(T), T_max = max(T), .groups = "drop") %>%
    right_join(dose_key_tbl, by = c("Dose", "Dose_key")) %>%
    mutate(
      T_min = ifelse(is.na(T_min), min(dat$T, na.rm = TRUE), T_min),
      T_max = ifelse(is.na(T_max), max(dat$T, na.rm = TRUE), T_max)
    )
  
  purrr::map_dfr(seq_len(nrow(dose_key_tbl)), function(i) {
    d  <- dose_key_tbl$Dose[i]
    dk <- dose_key_tbl$Dose_key[i]
    
    T_seq <- seq(
      grid_base$T_min[grid_base$Dose_key == dk],
      grid_base$T_max[grid_base$Dose_key == dk],
      length.out = 200
    )
    TK_seq <- T_seq + 273.15
    
    lnB0 <- dr[[paste0("b_lnB0_Dose_key", dk)]]
    E    <- dr[[paste0("b_E_Dose_key", dk)]]
    Eh   <- dr[[paste0("b_Eh_Dose_key", dk)]]
    Th   <- dr[[paste0("b_Th_Dose_key", dk)]]
    
    pred_log_mat <- sapply(seq_along(TK_seq), function(j) {
      lnB0 -
        (E * 11604.51812) * ((1 / TK_seq[j]) - (1 / 293.15)) -
        log(1 + exp((Eh * 11604.51812) * ((1 / Th) - (1 / TK_seq[j]))))
    })
    pred_raw_mat <- exp(pred_log_mat)
    
    tibble(
      Dose = d,
      Dose_key = dk,
      T = T_seq,
      log_q2.5  = apply(pred_log_mat, 2, quantile, probs = 0.025, na.rm = TRUE),
      log_q50   = apply(pred_log_mat, 2, quantile, probs = 0.5,   na.rm = TRUE),
      log_q97.5 = apply(pred_log_mat, 2, quantile, probs = 0.975, na.rm = TRUE),
      raw_q2.5  = apply(pred_raw_mat, 2, quantile, probs = 0.025, na.rm = TRUE),
      raw_q50   = apply(pred_raw_mat, 2, quantile, probs = 0.5,   na.rm = TRUE),
      raw_q97.5 = apply(pred_raw_mat, 2, quantile, probs = 0.975, na.rm = TRUE)
    )
  })
}

arr_grid_growth <- if (!is.null(fit_growth_brm)) make_arr_grid_from_posterior(fit_growth_brm, growth_dat, dose_key_tbl) else tibble()
arr_grid_resp   <- if (!is.null(fit_resp_brm))   make_arr_grid_from_posterior(fit_resp_brm, resp_dat, dose_key_tbl) else tibble()
arr_grid_ratio  <- if (!is.null(fit_ratio_brm))  make_arr_grid_from_posterior(fit_ratio_brm, ratio_dat_arr, dose_key_tbl) else tibble()
arr_grid_cue    <- if (!is.null(fit_cue_brm))    make_arr_grid_from_posterior_common_slope(fit_cue_brm, cue_dat, dose_key_tbl) else tibble()

arr_grid_growth_biomass <- if (!is.null(fit_growth_biomass_brm)) make_arr_grid_from_posterior(fit_growth_biomass_brm, growth_biomass_dat, dose_key_tbl) else tibble()
arr_grid_resp_biomass   <- if (!is.null(fit_resp_biomass_brm))   make_arr_grid_from_posterior(fit_resp_biomass_brm, resp_biomass_dat, dose_key_tbl) else tibble()

tpc_grid_growth <- if (!is.null(fit_growth_tpc_brm)) make_ss_grid_from_posterior(fit_growth_tpc_brm, growth_dat, dose_key_tbl) else tibble()
tpc_grid_resp   <- if (!is.null(fit_resp_tpc_brm))   make_ss_grid_from_posterior(fit_resp_tpc_brm, resp_dat, dose_key_tbl) else tibble()

tpc_grid_growth_biomass <- if (!is.null(fit_growth_tpc_biomass_brm)) make_ss_grid_from_posterior(fit_growth_tpc_biomass_brm, growth_biomass_dat, dose_key_tbl) else tibble()
tpc_grid_resp_biomass   <- if (!is.null(fit_resp_tpc_biomass_brm))   make_ss_grid_from_posterior(fit_resp_tpc_biomass_brm, resp_biomass_dat, dose_key_tbl) else tibble()

arr_grid_growth_hier <- if (!is.null(fit_growth_hier_brm)) make_arr_grid_from_posterior(fit_growth_hier_brm, growth_dat, dose_key_tbl) else tibble()
arr_grid_resp_hier   <- if (!is.null(fit_resp_hier_brm))   make_arr_grid_from_posterior(fit_resp_hier_brm, resp_dat, dose_key_tbl) else tibble()
arr_grid_ratio_hier  <- if (!is.null(fit_ratio_hier_brm))  make_arr_grid_from_posterior(fit_ratio_hier_brm, ratio_dat_arr, dose_key_tbl) else tibble()

arr_grid_growth_biomass_hier <- if (!is.null(fit_growth_biomass_hier_brm)) make_arr_grid_from_posterior(fit_growth_biomass_hier_brm, growth_biomass_dat, dose_key_tbl) else tibble()
arr_grid_resp_biomass_hier   <- if (!is.null(fit_resp_biomass_hier_brm))   make_arr_grid_from_posterior(fit_resp_biomass_hier_brm, resp_biomass_dat, dose_key_tbl) else tibble()

# ===== Posterior distribution plots for SS parameters =========================
make_ss_posterior_plot <- function(post_df, title_txt) {
  post_long <- post_df %>%
    select(Dose, lnB0, E, Eh, Th) %>%
    pivot_longer(cols = c(lnB0, E, Eh, Th), names_to = "parameter", values_to = "value") %>%
    mutate(parameter = factor(parameter, levels = c("E", "Eh", "lnB0", "Th")))
  
  ggplot(post_long, aes(x = value, color = Dose, fill = Dose)) +
    geom_density(alpha = 0.18, linewidth = 0.9) +
    facet_wrap(~parameter, scales = "free", ncol = 2) +
    scale_color_manual(values = dose_cols, limits = dose_levels, drop = FALSE) +
    scale_fill_manual(values = dose_cols, limits = dose_levels, drop = FALSE) +
    labs(
      title = title_txt,
      x = "Posterior value",
      y = "Density",
      color = "Dose",
      fill = "Dose"
    ) +
    theme_classic(12)
}

p_postdist_growth <- if (file.exists(post_tpc_growth_csv)) {
  make_ss_posterior_plot(
    readr::read_csv(post_tpc_growth_csv, show_col_types = FALSE),
    "Sharpe–Schoolfield posterior distributions: growth"
  )
} else NULL

p_postdist_resp <- if (file.exists(post_tpc_resp_csv)) {
  make_ss_posterior_plot(
    readr::read_csv(post_tpc_resp_csv, show_col_types = FALSE),
    "Sharpe–Schoolfield posterior distributions: respiration"
  )
} else NULL

p_postdist_growth_biomass <- if (file.exists(post_tpc_growth_biomass_csv)) {
  make_ss_posterior_plot(
    readr::read_csv(post_tpc_growth_biomass_csv, show_col_types = FALSE),
    "Sharpe–Schoolfield posterior distributions: biomass-corrected growth"
  )
} else NULL

p_postdist_resp_biomass <- if (file.exists(post_tpc_resp_biomass_csv)) {
  make_ss_posterior_plot(
    readr::read_csv(post_tpc_resp_biomass_csv, show_col_types = FALSE),
    "Sharpe–Schoolfield posterior distributions: biomass-corrected respiration"
  )
} else NULL

# ===== Bayesian plots =========================================================
p_growth_bayes <- ggplot(growth_dat, aes(T, y_raw, color = Dose, fill = Dose)) +
  geom_point(size = 2.1, alpha = 0.9) +
  { if (nrow(arr_grid_growth) > 0) geom_ribbon(data = arr_grid_growth, aes(x = T, ymin = raw_q2.5, ymax = raw_q97.5, fill = Dose), inherit.aes = FALSE, alpha = 0.18, colour = NA) } +
  { if (nrow(arr_grid_growth) > 0) geom_line(data = arr_grid_growth, aes(x = T, y = raw_q50, color = Dose), inherit.aes = FALSE, linewidth = 1) } +
  scale_color_manual(values = dose_cols, limits = dose_levels, drop = FALSE) +
  scale_fill_manual(values = dose_cols, limits = dose_levels, drop = FALSE) +
  labs(title = "Bayesian Arrhenius thermal fit: growth (by Dose)", x = "Temperature (°C)", y = "Growth (fg C h^-1)", color = "Dose", fill = "Dose") +
  theme_classic(12)

p_resp_bayes <- ggplot(resp_dat, aes(T, y_raw, color = Dose, fill = Dose)) +
  geom_point(size = 2.1, alpha = 0.9) +
  { if (nrow(arr_grid_resp) > 0) geom_ribbon(data = arr_grid_resp, aes(x = T, ymin = raw_q2.5, ymax = raw_q97.5, fill = Dose), inherit.aes = FALSE, alpha = 0.18, colour = NA) } +
  { if (nrow(arr_grid_resp) > 0) geom_line(data = arr_grid_resp, aes(x = T, y = raw_q50, color = Dose), inherit.aes = FALSE, linewidth = 1) } +
  scale_color_manual(values = dose_cols, limits = dose_levels, drop = FALSE) +
  scale_fill_manual(values = dose_cols, limits = dose_levels, drop = FALSE) +
  labs(title = "Bayesian Arrhenius thermal fit: respiration (by Dose)", x = "Temperature (°C)", y = "Respiration (fg C h^-1)", color = "Dose", fill = "Dose") +
  theme_classic(12)

p_ratio_bayes <- ggplot(ratio_dat_arr, aes(T, y_raw, color = Dose, fill = Dose)) +
  geom_point(size = 2.1, alpha = 0.9) +
  { if (nrow(arr_grid_ratio) > 0) geom_ribbon(data = arr_grid_ratio, aes(x = T, ymin = raw_q2.5, ymax = raw_q97.5, fill = Dose), inherit.aes = FALSE, alpha = 0.18, colour = NA) } +
  { if (nrow(arr_grid_ratio) > 0) geom_line(data = arr_grid_ratio, aes(x = T, y = raw_q50, color = Dose), inherit.aes = FALSE, linewidth = 1) } +
  scale_color_manual(values = dose_cols, limits = dose_levels, drop = FALSE) +
  scale_fill_manual(values = dose_cols, limits = dose_levels, drop = FALSE) +
  labs(title = "Bayesian Arrhenius thermal fit: respiration / growth (by Dose)", x = "Temperature (°C)", y = "Respiration / Growth", color = "Dose", fill = "Dose") +
  theme_classic(12)

p_cue_arrhenius_bayes <- ggplot(cue_dat, aes(T, y_raw, color = Dose, fill = Dose)) +
  geom_point(size = 2.1, alpha = 0.9) +
  { if (nrow(arr_grid_cue) > 0) geom_ribbon(
    data = arr_grid_cue,
    aes(x = T, ymin = raw_q2.5, ymax = raw_q97.5, fill = Dose),
    inherit.aes = FALSE,
    alpha = 0.18,
    colour = NA
  ) } +
  { if (nrow(arr_grid_cue) > 0) geom_line(
    data = arr_grid_cue,
    aes(x = T, y = raw_q50, color = Dose),
    inherit.aes = FALSE,
    linewidth = 1
  ) } +
  scale_y_log10() +
  scale_color_manual(values = dose_cols, limits = dose_levels, drop = FALSE) +
  scale_fill_manual(values = dose_cols, limits = dose_levels, drop = FALSE) +
  labs(
    title = "Bayesian Arrhenius thermal fit: CUE vs Temperature (common slope)",
    subtitle = if (!is.null(fit_cue_brm)) {
      cue_draws_tmp <- posterior::as_draws_df(fit_cue_brm)
      cue_slope_med <- stats::median(cue_draws_tmp$b_E_Intercept)
      paste0("Common slope E (posterior median) = ", signif(cue_slope_med, 4))
    } else {
      NULL
    },
    x = "Temperature (°C)",
    y = "CUE (log scale)",
    color = "Dose",
    fill = "Dose"
  ) +
  theme_classic(12)

p_growth_biomass_bayes <- ggplot(growth_biomass_dat, aes(T, y_raw, color = Dose, fill = Dose)) +
  geom_point(size = 2.1, alpha = 0.9) +
  { if (nrow(arr_grid_growth_biomass) > 0) geom_ribbon(data = arr_grid_growth_biomass, aes(x = T, ymin = raw_q2.5, ymax = raw_q97.5, fill = Dose), inherit.aes = FALSE, alpha = 0.18, colour = NA) } +
  { if (nrow(arr_grid_growth_biomass) > 0) geom_line(data = arr_grid_growth_biomass, aes(x = T, y = raw_q50, color = Dose), inherit.aes = FALSE, linewidth = 1) } +
  scale_color_manual(values = dose_cols, limits = dose_levels, drop = FALSE) +
  scale_fill_manual(values = dose_cols, limits = dose_levels, drop = FALSE) +
  labs(title = "Bayesian Arrhenius thermal fit: biomass-corrected growth (by Dose)", x = "Temperature (°C)", y = "Growth (C per C per h)", color = "Dose", fill = "Dose") +
  theme_classic(12)

p_resp_biomass_bayes <- ggplot(resp_biomass_dat, aes(T, y_raw, color = Dose, fill = Dose)) +
  geom_point(size = 2.1, alpha = 0.9) +
  { if (nrow(arr_grid_resp_biomass) > 0) geom_ribbon(data = arr_grid_resp_biomass, aes(x = T, ymin = raw_q2.5, ymax = raw_q97.5, fill = Dose), inherit.aes = FALSE, alpha = 0.18, colour = NA) } +
  { if (nrow(arr_grid_resp_biomass) > 0) geom_line(data = arr_grid_resp_biomass, aes(x = T, y = raw_q50, color = Dose), inherit.aes = FALSE, linewidth = 1) } +
  scale_color_manual(values = dose_cols, limits = dose_levels, drop = FALSE) +
  scale_fill_manual(values = dose_cols, limits = dose_levels, drop = FALSE) +
  labs(title = "Bayesian Arrhenius thermal fit: biomass-corrected respiration (by Dose)", x = "Temperature (°C)", y = "Respiration (C per C per h)", color = "Dose", fill = "Dose") +
  theme_classic(12)

p_growth_tpc_bayes <- ggplot(growth_dat, aes(T, y_raw, color = Dose, fill = Dose)) +
  geom_point(size = 2.1, alpha = 0.9) +
  { if (nrow(tpc_grid_growth) > 0) geom_ribbon(data = tpc_grid_growth, aes(x = T, ymin = raw_q2.5, ymax = raw_q97.5, fill = Dose), inherit.aes = FALSE, alpha = 0.18, colour = NA) } +
  { if (nrow(tpc_grid_growth) > 0) geom_line(data = tpc_grid_growth, aes(x = T, y = raw_q50, color = Dose), inherit.aes = FALSE, linewidth = 1) } +
  scale_color_manual(values = dose_cols, limits = dose_levels, drop = FALSE) +
  scale_fill_manual(values = dose_cols, limits = dose_levels, drop = FALSE) +
  scale_y_log10() +
  labs(title = "Bayesian full TPC: growth (Sharpe–Schoolfield, by Dose)", x = "Temperature (°C)", y = "Growth (fg C h^-1, log scale)", color = "Dose", fill = "Dose") +
  theme_classic(12)

p_resp_tpc_bayes <- ggplot(resp_dat, aes(T, y_raw, color = Dose, fill = Dose)) +
  geom_point(size = 2.1, alpha = 0.9) +
  { if (nrow(tpc_grid_resp) > 0) geom_ribbon(data = tpc_grid_resp, aes(x = T, ymin = raw_q2.5, ymax = raw_q97.5, fill = Dose), inherit.aes = FALSE, alpha = 0.18, colour = NA) } +
  { if (nrow(tpc_grid_resp) > 0) geom_line(data = tpc_grid_resp, aes(x = T, y = raw_q50, color = Dose), inherit.aes = FALSE, linewidth = 1) } +
  scale_color_manual(values = dose_cols, limits = dose_levels, drop = FALSE) +
  scale_fill_manual(values = dose_cols, limits = dose_levels, drop = FALSE) +
  scale_y_log10() +
  labs(title = "Bayesian full TPC: respiration (Sharpe–Schoolfield, by Dose)", x = "Temperature (°C)", y = "Respiration (fg C h^-1, log scale)", color = "Dose", fill = "Dose") +
  theme_classic(12)

p_growth_tpc_biomass_bayes <- ggplot(growth_biomass_dat, aes(T, y_raw, color = Dose, fill = Dose)) +
  geom_point(size = 2.1, alpha = 0.9) +
  { if (nrow(tpc_grid_growth_biomass) > 0) geom_ribbon(data = tpc_grid_growth_biomass, aes(x = T, ymin = raw_q2.5, ymax = raw_q97.5, fill = Dose), inherit.aes = FALSE, alpha = 0.18, colour = NA) } +
  { if (nrow(tpc_grid_growth_biomass) > 0) geom_line(data = tpc_grid_growth_biomass, aes(x = T, y = raw_q50, color = Dose), inherit.aes = FALSE, linewidth = 1) } +
  scale_color_manual(values = dose_cols, limits = dose_levels, drop = FALSE) +
  scale_fill_manual(values = dose_cols, limits = dose_levels, drop = FALSE) +
  scale_y_log10() +
  labs(title = "Bayesian full TPC: biomass-corrected growth (Sharpe–Schoolfield, by Dose)", x = "Temperature (°C)", y = "Growth (C per C per h, log scale)", color = "Dose", fill = "Dose") +
  theme_classic(12)

p_resp_tpc_biomass_bayes <- ggplot(resp_biomass_dat, aes(T, y_raw, color = Dose, fill = Dose)) +
  geom_point(size = 2.1, alpha = 0.9) +
  { if (nrow(tpc_grid_resp_biomass) > 0) geom_ribbon(data = tpc_grid_resp_biomass, aes(x = T, ymin = raw_q2.5, ymax = raw_q97.5, fill = Dose), inherit.aes = FALSE, alpha = 0.18, colour = NA) } +
  { if (nrow(tpc_grid_resp_biomass) > 0) geom_line(data = tpc_grid_resp_biomass, aes(x = T, y = raw_q50, color = Dose), inherit.aes = FALSE, linewidth = 1) } +
  scale_color_manual(values = dose_cols, limits = dose_levels, drop = FALSE) +
  scale_fill_manual(values = dose_cols, limits = dose_levels, drop = FALSE) +
  scale_y_log10() +
  labs(title = "Bayesian full TPC: biomass-corrected respiration (Sharpe–Schoolfield, by Dose)", x = "Temperature (°C)", y = "Respiration (C per C per h, log scale)", color = "Dose", fill = "Dose") +
  theme_classic(12)

# ===== Separate TPC check plots by dose ======================================
p_tpc_check_growth <- if (nrow(growth_dat) > 0) {
  make_tpc_check_plot(
    dat = growth_dat,
    grid = tpc_grid_growth,
    title_txt = "TPC check: growth curves by dose",
    ylab_txt = "Growth (fg C h^-1)",
    log_scale = FALSE
  )
} else NULL

p_tpc_check_resp <- if (nrow(resp_dat) > 0) {
  make_tpc_check_plot(
    dat = resp_dat,
    grid = tpc_grid_resp,
    title_txt = "TPC check: respiration curves by dose",
    ylab_txt = "Respiration (fg C h^-1)",
    log_scale = FALSE
  )
} else NULL

p_tpc_check_growth_biomass <- if (nrow(growth_biomass_dat) > 0) {
  make_tpc_check_plot(
    dat = growth_biomass_dat,
    grid = tpc_grid_growth_biomass,
    title_txt = "TPC check: biomass-corrected growth curves by dose",
    ylab_txt = "Growth (C per C per h)",
    log_scale = FALSE
  )
} else NULL

p_tpc_check_resp_biomass <- if (nrow(resp_biomass_dat) > 0) {
  make_tpc_check_plot(
    dat = resp_biomass_dat,
    grid = tpc_grid_resp_biomass,
    title_txt = "TPC check: biomass-corrected respiration curves by dose",
    ylab_txt = "Respiration (C per C per h)",
    log_scale = FALSE
  )
} else NULL

p_tpc_check_growth_log <- if (nrow(growth_dat) > 0) {
  make_tpc_check_plot(
    dat = growth_dat,
    grid = tpc_grid_growth,
    title_txt = "TPC check: growth curves by dose (log scale)",
    ylab_txt = "Growth (fg C h^-1, log scale)",
    log_scale = TRUE
  )
} else NULL

p_tpc_check_resp_log <- if (nrow(resp_dat) > 0) {
  make_tpc_check_plot(
    dat = resp_dat,
    grid = tpc_grid_resp,
    title_txt = "TPC check: respiration curves by dose (log scale)",
    ylab_txt = "Respiration (fg C h^-1, log scale)",
    log_scale = TRUE
  )
} else NULL

p_tpc_check_growth_biomass_log <- if (nrow(growth_biomass_dat) > 0) {
  make_tpc_check_plot(
    dat = growth_biomass_dat,
    grid = tpc_grid_growth_biomass,
    title_txt = "TPC check: biomass-corrected growth curves by dose (log scale)",
    ylab_txt = "Growth (C per C per h, log scale)",
    log_scale = TRUE
  )
} else NULL

p_tpc_check_resp_biomass_log <- if (nrow(resp_biomass_dat) > 0) {
  make_tpc_check_plot(
    dat = resp_biomass_dat,
    grid = tpc_grid_resp_biomass,
    title_txt = "TPC check: biomass-corrected respiration curves by dose (log scale)",
    ylab_txt = "Respiration (C per C per h, log scale)",
    log_scale = TRUE
  )
} else NULL

# ===== Hierarchical Arrhenius plots ==========================================
p_growth_hier_bayes <- if (!is.null(fit_growth_hier_brm)) {
  ggplot(growth_dat, aes(T, y_raw, color = Dose, fill = Dose)) +
    geom_point(size = 2.1, alpha = 0.9) +
    { if (nrow(arr_grid_growth_hier) > 0) geom_ribbon(data = arr_grid_growth_hier, aes(x = T, ymin = raw_q2.5, ymax = raw_q97.5, fill = Dose), inherit.aes = FALSE, alpha = 0.18, colour = NA) } +
    { if (nrow(arr_grid_growth_hier) > 0) geom_line(data = arr_grid_growth_hier, aes(x = T, y = raw_q50, color = Dose), inherit.aes = FALSE, linewidth = 1) } +
    scale_color_manual(values = dose_cols, limits = dose_levels, drop = FALSE) +
    scale_fill_manual(values = dose_cols, limits = dose_levels, drop = FALSE) +
    labs(title = "Hierarchical Bayesian Arrhenius: growth (by Dose)", x = "Temperature (°C)", y = "Growth (fg C h^-1)", color = "Dose", fill = "Dose") +
    theme_classic(12)
} else NULL

p_resp_hier_bayes <- if (!is.null(fit_resp_hier_brm)) {
  ggplot(resp_dat, aes(T, y_raw, color = Dose, fill = Dose)) +
    geom_point(size = 2.1, alpha = 0.9) +
    { if (nrow(arr_grid_resp_hier) > 0) geom_ribbon(data = arr_grid_resp_hier, aes(x = T, ymin = raw_q2.5, ymax = raw_q97.5, fill = Dose), inherit.aes = FALSE, alpha = 0.18, colour = NA) } +
    { if (nrow(arr_grid_resp_hier) > 0) geom_line(data = arr_grid_resp_hier, aes(x = T, y = raw_q50, color = Dose), inherit.aes = FALSE, linewidth = 1) } +
    scale_color_manual(values = dose_cols, limits = dose_levels, drop = FALSE) +
    scale_fill_manual(values = dose_cols, limits = dose_levels, drop = FALSE) +
    labs(title = "Hierarchical Bayesian Arrhenius: respiration (by Dose)", x = "Temperature (°C)", y = "Respiration (fg C h^-1)", color = "Dose", fill = "Dose") +
    theme_classic(12)
} else NULL

p_ratio_hier_bayes <- if (!is.null(fit_ratio_hier_brm)) {
  ggplot(ratio_dat_arr, aes(T, y_raw, color = Dose, fill = Dose)) +
    geom_point(size = 2.1, alpha = 0.9) +
    { if (nrow(arr_grid_ratio_hier) > 0) geom_ribbon(data = arr_grid_ratio_hier, aes(x = T, ymin = raw_q2.5, ymax = raw_q97.5, fill = Dose), inherit.aes = FALSE, alpha = 0.18, colour = NA) } +
    { if (nrow(arr_grid_ratio_hier) > 0) geom_line(data = arr_grid_ratio_hier, aes(x = T, y = raw_q50, color = Dose), inherit.aes = FALSE, linewidth = 1) } +
    scale_color_manual(values = dose_cols, limits = dose_levels, drop = FALSE) +
    scale_fill_manual(values = dose_cols, limits = dose_levels, drop = FALSE) +
    labs(title = "Hierarchical Bayesian Arrhenius: respiration / growth (by Dose)", x = "Temperature (°C)", y = "Respiration / Growth", color = "Dose", fill = "Dose") +
    theme_classic(12)
} else NULL

p_growth_biomass_hier_bayes <- if (!is.null(fit_growth_biomass_hier_brm)) {
  ggplot(growth_biomass_dat, aes(T, y_raw, color = Dose, fill = Dose)) +
    geom_point(size = 2.1, alpha = 0.9) +
    { if (nrow(arr_grid_growth_biomass_hier) > 0) geom_ribbon(data = arr_grid_growth_biomass_hier, aes(x = T, ymin = raw_q2.5, ymax = raw_q97.5, fill = Dose), inherit.aes = FALSE, alpha = 0.18, colour = NA) } +
    { if (nrow(arr_grid_growth_biomass_hier) > 0) geom_line(data = arr_grid_growth_biomass_hier, aes(x = T, y = raw_q50, color = Dose), inherit.aes = FALSE, linewidth = 1) } +
    scale_color_manual(values = dose_cols, limits = dose_levels, drop = FALSE) +
    scale_fill_manual(values = dose_cols, limits = dose_levels, drop = FALSE) +
    labs(title = "Hierarchical Bayesian Arrhenius: biomass-corrected growth (by Dose)", x = "Temperature (°C)", y = "Growth (C per C per h)", color = "Dose", fill = "Dose") +
    theme_classic(12)
} else NULL

p_resp_biomass_hier_bayes <- if (!is.null(fit_resp_biomass_hier_brm)) {
  ggplot(resp_biomass_dat, aes(T, y_raw, color = Dose, fill = Dose)) +
    geom_point(size = 2.1, alpha = 0.9) +
    { if (nrow(arr_grid_resp_biomass_hier) > 0) geom_ribbon(data = arr_grid_resp_biomass_hier, aes(x = T, ymin = raw_q2.5, ymax = raw_q97.5, fill = Dose), inherit.aes = FALSE, alpha = 0.18, colour = NA) } +
    { if (nrow(arr_grid_resp_biomass_hier) > 0) geom_line(data = arr_grid_resp_biomass_hier, aes(x = T, y = raw_q50, color = Dose), inherit.aes = FALSE, linewidth = 1) } +
    scale_color_manual(values = dose_cols, limits = dose_levels, drop = FALSE) +
    scale_fill_manual(values = dose_cols, limits = dose_levels, drop = FALSE) +
    labs(title = "Hierarchical Bayesian Arrhenius: biomass-corrected respiration (by Dose)", x = "Temperature (°C)", y = "Respiration (C per C per h)", color = "Dose", fill = "Dose") +
    theme_classic(12)
} else NULL

# ===== Boltzmann plots ========================================================
p_boltz_growth <- ggplot(growth_dat, aes(boltz_shift, y, color = Dose, fill = Dose)) +
  geom_point(size = 2.1, alpha = 0.9) +
  { if (nrow(arr_grid_growth) > 0) geom_ribbon(data = arr_grid_growth, aes(x = boltz_shift, ymin = log_q2.5, ymax = log_q97.5, fill = Dose), inherit.aes = FALSE, alpha = 0.18, colour = NA) } +
  { if (nrow(arr_grid_growth) > 0) geom_line(data = arr_grid_growth, aes(x = boltz_shift, y = log_q50, color = Dose), inherit.aes = FALSE, linewidth = 1) } +
  scale_color_manual(values = dose_cols, limits = dose_levels, drop = FALSE) +
  scale_fill_manual(values = dose_cols, limits = dose_levels, drop = FALSE) +
  labs(title = "Bayesian Boltzmann plot: log(growth) vs 1/kTref - 1/kT (by Dose)", x = expression(frac(1, k[T[ref]]) - frac(1, k*T)), y = "log(growth fg C h^-1)", color = "Dose", fill = "Dose") +
  theme_classic(12)

p_boltz_resp <- ggplot(resp_dat, aes(boltz_shift, y, color = Dose, fill = Dose)) +
  geom_point(size = 2.1, alpha = 0.9) +
  { if (nrow(arr_grid_resp) > 0) geom_ribbon(data = arr_grid_resp, aes(x = boltz_shift, ymin = log_q2.5, ymax = log_q97.5, fill = Dose), inherit.aes = FALSE, alpha = 0.18, colour = NA) } +
  { if (nrow(arr_grid_resp) > 0) geom_line(data = arr_grid_resp, aes(x = boltz_shift, y = log_q50, color = Dose), inherit.aes = FALSE, linewidth = 1) } +
  scale_color_manual(values = dose_cols, limits = dose_levels, drop = FALSE) +
  scale_fill_manual(values = dose_cols, limits = dose_levels, drop = FALSE) +
  labs(title = "Bayesian Boltzmann plot: log(respiration) vs 1/kTref - 1/kT (by Dose)", x = expression(frac(1, k[T[ref]]) - frac(1, k*T)), y = "log(respiration fg C h^-1)", color = "Dose", fill = "Dose") +
  theme_classic(12)

p_boltz_ratio <- ggplot(ratio_dat_arr, aes(boltz_shift, y, color = Dose, fill = Dose)) +
  geom_point(size = 2.1, alpha = 0.9) +
  { if (nrow(arr_grid_ratio) > 0) geom_ribbon(data = arr_grid_ratio, aes(x = boltz_shift, ymin = log_q2.5, ymax = log_q97.5, fill = Dose), inherit.aes = FALSE, alpha = 0.18, colour = NA) } +
  { if (nrow(arr_grid_ratio) > 0) geom_line(data = arr_grid_ratio, aes(x = boltz_shift, y = log_q50, color = Dose), inherit.aes = FALSE, linewidth = 1) } +
  scale_color_manual(values = dose_cols, limits = dose_levels, drop = FALSE) +
  scale_fill_manual(values = dose_cols, limits = dose_levels, drop = FALSE) +
  labs(title = "Bayesian Boltzmann plot: log(respiration / growth) vs 1/kTref - 1/kT (by Dose)", x = expression(frac(1, k[T[ref]]) - frac(1, k*T)), y = "log(respiration / growth)", color = "Dose", fill = "Dose") +
  theme_classic(12)

p_boltz_cue <- ggplot(cue_dat, aes(boltz_shift, y_raw, color = Dose, fill = Dose)) +
  geom_point(size = 2.1, alpha = 0.9) +
  { if (nrow(arr_grid_cue) > 0) geom_ribbon(
    data = arr_grid_cue,
    aes(x = boltz_shift, ymin = raw_q2.5, ymax = raw_q97.5, fill = Dose),
    inherit.aes = FALSE,
    alpha = 0.18,
    colour = NA
  ) } +
  { if (nrow(arr_grid_cue) > 0) geom_line(
    data = arr_grid_cue,
    aes(x = boltz_shift, y = raw_q50, color = Dose),
    inherit.aes = FALSE,
    linewidth = 1
  ) } +
  scale_y_log10() +
  scale_color_manual(values = dose_cols, limits = dose_levels, drop = FALSE) +
  scale_fill_manual(values = dose_cols, limits = dose_levels, drop = FALSE) +
  labs(
    title = "Bayesian Boltzmann plot: CUE vs 1/kTref - 1/kT (common slope)",
    subtitle = if (!is.null(fit_cue_brm)) {
      cue_draws_tmp <- posterior::as_draws_df(fit_cue_brm)
      cue_slope_med <- stats::median(cue_draws_tmp$b_E_Intercept)
      paste0("Common slope E (posterior median) = ", signif(cue_slope_med, 4))
    } else {
      NULL
    },
    x = expression(frac(1, k[T[ref]]) - frac(1, k*T)),
    y = "CUE (log scale)",
    color = "Dose",
    fill = "Dose"
  ) +
  theme_classic(12)

p_boltz_growth_biomass <- ggplot(growth_biomass_dat, aes(boltz_shift, y, color = Dose, fill = Dose)) +
  geom_point(size = 2.1, alpha = 0.9) +
  { if (nrow(arr_grid_growth_biomass) > 0) geom_ribbon(data = arr_grid_growth_biomass, aes(x = boltz_shift, ymin = log_q2.5, ymax = log_q97.5, fill = Dose), inherit.aes = FALSE, alpha = 0.18, colour = NA) } +
  { if (nrow(arr_grid_growth_biomass) > 0) geom_line(data = arr_grid_growth_biomass, aes(x = boltz_shift, y = log_q50, color = Dose), inherit.aes = FALSE, linewidth = 1) } +
  scale_color_manual(values = dose_cols, limits = dose_levels, drop = FALSE) +
  scale_fill_manual(values = dose_cols, limits = dose_levels, drop = FALSE) +
  labs(title = "Bayesian Boltzmann plot: log(biomass-corrected growth) vs 1/kTref - 1/kT (by Dose)", x = expression(frac(1, k[T[ref]]) - frac(1, k*T)), y = "log(growth C per C per h)", color = "Dose", fill = "Dose") +
  theme_classic(12)

p_boltz_resp_biomass <- ggplot(resp_biomass_dat, aes(boltz_shift, y, color = Dose, fill = Dose)) +
  geom_point(size = 2.1, alpha = 0.9) +
  { if (nrow(arr_grid_resp_biomass) > 0) geom_ribbon(data = arr_grid_resp_biomass, aes(x = boltz_shift, ymin = log_q2.5, ymax = log_q97.5, fill = Dose), inherit.aes = FALSE, alpha = 0.18, colour = NA) } +
  { if (nrow(arr_grid_resp_biomass) > 0) geom_line(data = arr_grid_resp_biomass, aes(x = boltz_shift, y = log_q50, color = Dose), inherit.aes = FALSE, linewidth = 1) } +
  scale_color_manual(values = dose_cols, limits = dose_levels, drop = FALSE) +
  scale_fill_manual(values = dose_cols, limits = dose_levels, drop = FALSE) +
  labs(title = "Bayesian Boltzmann plot: log(biomass-corrected respiration) vs 1/kTref - 1/kT (by Dose)", x = expression(frac(1, k[T[ref]]) - frac(1, k*T)), y = "log(respiration C per C per h)", color = "Dose", fill = "Dose") +
  theme_classic(12)

p_boltz_growth_hier <- if (!is.null(fit_growth_hier_brm)) {
  ggplot(growth_dat, aes(boltz_shift, y, color = Dose, fill = Dose)) +
    geom_point(size = 2.1, alpha = 0.9) +
    { if (nrow(arr_grid_growth_hier) > 0) geom_ribbon(data = arr_grid_growth_hier, aes(x = boltz_shift, ymin = log_q2.5, ymax = log_q97.5, fill = Dose), inherit.aes = FALSE, alpha = 0.18, colour = NA) } +
    { if (nrow(arr_grid_growth_hier) > 0) geom_line(data = arr_grid_growth_hier, aes(x = boltz_shift, y = log_q50, color = Dose), inherit.aes = FALSE, linewidth = 1) } +
    scale_color_manual(values = dose_cols, limits = dose_levels, drop = FALSE) +
    scale_fill_manual(values = dose_cols, limits = dose_levels, drop = FALSE) +
    labs(title = "Hierarchical Bayesian Boltzmann plot: log(growth) vs 1/kTref - 1/kT", x = expression(frac(1, k[T[ref]]) - frac(1, k*T)), y = "log(growth fg C h^-1)", color = "Dose", fill = "Dose") +
    theme_classic(12)
} else NULL

p_boltz_resp_hier <- if (!is.null(fit_resp_hier_brm)) {
  ggplot(resp_dat, aes(boltz_shift, y, color = Dose, fill = Dose)) +
    geom_point(size = 2.1, alpha = 0.9) +
    { if (nrow(arr_grid_resp_hier) > 0) geom_ribbon(data = arr_grid_resp_hier, aes(x = boltz_shift, ymin = log_q2.5, ymax = log_q97.5, fill = Dose), inherit.aes = FALSE, alpha = 0.18, colour = NA) } +
    { if (nrow(arr_grid_resp_hier) > 0) geom_line(data = arr_grid_resp_hier, aes(x = boltz_shift, y = log_q50, color = Dose), inherit.aes = FALSE, linewidth = 1) } +
    scale_color_manual(values = dose_cols, limits = dose_levels, drop = FALSE) +
    scale_fill_manual(values = dose_cols, limits = dose_levels, drop = FALSE) +
    labs(title = "Hierarchical Bayesian Boltzmann plot: log(respiration) vs 1/kTref - 1/kT", x = expression(frac(1, k[T[ref]]) - frac(1, k*T)), y = "log(respiration fg C h^-1)", color = "Dose", fill = "Dose") +
    theme_classic(12)
} else NULL

p_boltz_ratio_hier <- if (!is.null(fit_ratio_hier_brm)) {
  ggplot(ratio_dat_arr, aes(boltz_shift, y, color = Dose, fill = Dose)) +
    geom_point(size = 2.1, alpha = 0.9) +
    { if (nrow(arr_grid_ratio_hier) > 0) geom_ribbon(data = arr_grid_ratio_hier, aes(x = boltz_shift, ymin = log_q2.5, ymax = log_q97.5, fill = Dose), inherit.aes = FALSE, alpha = 0.18, colour = NA) } +
    { if (nrow(arr_grid_ratio_hier) > 0) geom_line(data = arr_grid_ratio_hier, aes(x = boltz_shift, y = log_q50, color = Dose), inherit.aes = FALSE, linewidth = 1) } +
    scale_color_manual(values = dose_cols, limits = dose_levels, drop = FALSE) +
    scale_fill_manual(values = dose_cols, limits = dose_levels, drop = FALSE) +
    labs(title = "Hierarchical Bayesian Boltzmann plot: log(respiration / growth) vs 1/kTref - 1/kT", x = expression(frac(1, k[T[ref]]) - frac(1, k*T)), y = "log(respiration / growth)", color = "Dose", fill = "Dose") +
    theme_classic(12)
} else NULL

p_boltz_growth_biomass_hier <- if (!is.null(fit_growth_biomass_hier_brm)) {
  ggplot(growth_biomass_dat, aes(boltz_shift, y, color = Dose, fill = Dose)) +
    geom_point(size = 2.1, alpha = 0.9) +
    { if (nrow(arr_grid_growth_biomass_hier) > 0) geom_ribbon(data = arr_grid_growth_biomass_hier, aes(x = boltz_shift, ymin = log_q2.5, ymax = log_q97.5, fill = Dose), inherit.aes = FALSE, alpha = 0.18, colour = NA) } +
    { if (nrow(arr_grid_growth_biomass_hier) > 0) geom_line(data = arr_grid_growth_biomass_hier, aes(x = boltz_shift, y = log_q50, color = Dose), inherit.aes = FALSE, linewidth = 1) } +
    scale_color_manual(values = dose_cols, limits = dose_levels, drop = FALSE) +
    scale_fill_manual(values = dose_cols, limits = dose_levels, drop = FALSE) +
    labs(title = "Hierarchical Bayesian Boltzmann plot: log(biomass-corrected growth) vs 1/kTref - 1/kT", x = expression(frac(1, k[T[ref]]) - frac(1, k*T)), y = "log(growth C per C per h)", color = "Dose", fill = "Dose") +
    theme_classic(12)
} else NULL

p_boltz_resp_biomass_hier <- if (!is.null(fit_resp_biomass_hier_brm)) {
  ggplot(resp_biomass_dat, aes(boltz_shift, y, color = Dose, fill = Dose)) +
    geom_point(size = 2.1, alpha = 0.9) +
    { if (nrow(arr_grid_resp_biomass_hier) > 0) geom_ribbon(data = arr_grid_resp_biomass_hier, aes(x = boltz_shift, ymin = log_q2.5, ymax = log_q97.5, fill = Dose), inherit.aes = FALSE, alpha = 0.18, colour = NA) } +
    { if (nrow(arr_grid_resp_biomass_hier) > 0) geom_line(data = arr_grid_resp_biomass_hier, aes(x = boltz_shift, y = log_q50, color = Dose), inherit.aes = FALSE, linewidth = 1) } +
    scale_color_manual(values = dose_cols, limits = dose_levels, drop = FALSE) +
    scale_fill_manual(values = dose_cols, limits = dose_levels, drop = FALSE) +
    labs(title = "Hierarchical Bayesian Boltzmann plot: log(biomass-corrected respiration) vs 1/kTref - 1/kT", x = expression(frac(1, k[T[ref]]) - frac(1, k*T)), y = "log(respiration C per C per h)", color = "Dose", fill = "Dose") +
    theme_classic(12)
} else NULL

# ===== Temperature-specific fungicide effect sizes ============================
DOSE_REF <- "Control"

has_effect_inputs <- file.exists(post_tpc_growth_csv) &&
  file.exists(post_tpc_resp_csv) &&
  file.exists(post_ratio_csv)

p_temp_effect <- NULL
slope_test_summary <- NULL
slope_test_draws   <- NULL

if (has_effect_inputs) {
  post_growth_eff <- readr::read_csv(post_tpc_growth_csv, show_col_types = FALSE) %>%
    mutate(Dose = as.character(Dose))
  post_resp_eff <- readr::read_csv(post_tpc_resp_csv, show_col_types = FALSE) %>%
    mutate(Dose = as.character(Dose))
  post_ratio_eff <- readr::read_csv(post_ratio_csv, show_col_types = FALSE) %>%
    mutate(Dose = as.character(Dose))
  
  growth_pred <- crossing(T = sort(unique(growth_dat$T)), Dose = dose_levels) %>%
    left_join(post_growth_eff, by = "Dose", relationship = "many-to-many") %>%
    mutate(
      TK = T + 273.15,
      pred_log = pred_ss_log(TK, lnB0, E, Eh, Th),
      pred_raw = exp(pred_log)
    ) %>%
    select(Dose, .draw, T, pred_log, pred_raw)
  
  resp_pred <- crossing(T = sort(unique(resp_dat$T)), Dose = dose_levels) %>%
    left_join(post_resp_eff, by = "Dose", relationship = "many-to-many") %>%
    mutate(
      TK = T + 273.15,
      pred_log = pred_ss_log(TK, lnB0, E, Eh, Th),
      pred_raw = exp(pred_log)
    ) %>%
    select(Dose, .draw, T, pred_log, pred_raw)
  
  ratio_pred <- crossing(T = sort(unique(ratio_dat_arr$T)), Dose = dose_levels) %>%
    left_join(post_ratio_eff, by = "Dose", relationship = "many-to-many") %>%
    mutate(
      TK = T + 273.15,
      pred_log = pred_arr_log(TK, alpha, E),
      pred_raw = exp(pred_log)
    ) %>%
    select(Dose, .draw, T, pred_log, pred_raw)
  
  cue_pred <- growth_pred %>%
    rename(growth_log = pred_log, growth_raw = pred_raw) %>%
    inner_join(
      resp_pred %>% rename(resp_log = pred_log, resp_raw = pred_raw),
      by = c("Dose", ".draw", "T")
    ) %>%
    mutate(
      CUE_raw = growth_raw / (growth_raw + resp_raw),
      CUE_logit = safe_logit(CUE_raw)
    ) %>%
    select(Dose, .draw, T, CUE_raw, CUE_logit)
  
  growth_eff <- compute_vs_control_effect(growth_pred, value_col = "pred_log", dose_ref = DOSE_REF) %>%
    filter(Dose != DOSE_REF) %>%
    mutate(Trait = "Growth")
  
  resp_eff <- compute_vs_control_effect(resp_pred, value_col = "pred_log", dose_ref = DOSE_REF) %>%
    filter(Dose != DOSE_REF) %>%
    mutate(Trait = "Respiration")
  
  ratio_eff <- compute_vs_control_effect(ratio_pred, value_col = "pred_log", dose_ref = DOSE_REF) %>%
    filter(Dose != DOSE_REF) %>%
    mutate(Trait = "Respiration / Growth")
  
  cue_eff <- compute_vs_control_effect(cue_pred, value_col = "CUE_logit", dose_ref = DOSE_REF) %>%
    filter(Dose != DOSE_REF) %>%
    mutate(Trait = "CUE")
  
  temp_effect_growth <- summarise_temp_effect(growth_eff, "Growth", effect_type = "negative")
  temp_effect_resp   <- summarise_temp_effect(resp_eff, "Respiration", effect_type = "positive")
  temp_effect_ratio  <- summarise_temp_effect(ratio_eff, "Respiration / Growth", effect_type = "positive")
  temp_effect_cue    <- summarise_temp_effect(cue_eff, "CUE", effect_type = "negative")
  
  temp_effect_summary <- bind_rows(
    temp_effect_growth,
    temp_effect_resp,
    temp_effect_ratio,
    temp_effect_cue
  ) %>%
    mutate(
      Trait = factor(Trait, levels = c("Growth", "Respiration", "Respiration / Growth", "CUE")),
      Dose = factor(Dose, levels = setdiff(dose_levels, "Control"))
    ) %>%
    arrange(Trait, Dose, T)
  
  temp_effect_draws <- bind_rows(
    growth_eff,
    resp_eff,
    ratio_eff,
    cue_eff
  ) %>%
    mutate(
      Trait = factor(Trait, levels = c("Growth", "Respiration", "Respiration / Growth", "CUE")),
      Dose = factor(Dose, levels = setdiff(dose_levels, "Control"))
    ) %>%
    arrange(Trait, Dose, T, .draw)
  
  readr::write_csv(temp_effect_summary, fungicide_temp_effect_summary_csv)
  readr::write_csv(temp_effect_draws, fungicide_temp_effect_draws_csv)
  
  p_temp_effect <- make_temp_effect_plot(temp_effect_summary)
  
  growth_slope_draws <- fit_effect_slopes_from_draws(growth_eff, "Growth")
  cue_slope_draws    <- fit_effect_slopes_from_draws(cue_eff, "CUE")
  
  slope_test_draws <- bind_rows(
    growth_slope_draws,
    cue_slope_draws
  ) %>%
    mutate(
      Trait = factor(Trait, levels = c("Growth", "CUE")),
      Dose  = factor(Dose, levels = setdiff(dose_levels, "Control"))
    ) %>%
    arrange(Trait, Dose, .draw)
  
  slope_test_summary <- summarise_effect_slopes(slope_test_draws) %>%
    mutate(
      Trait = factor(Trait, levels = c("Growth", "CUE")),
      Dose  = factor(Dose, levels = setdiff(dose_levels, "Control"))
    ) %>%
    arrange(Trait, Dose)
  
  readr::write_csv(slope_test_draws, fungicide_slope_test_draws_csv)
  readr::write_csv(slope_test_summary, fungicide_slope_test_summary_csv)
}

# ===== Save plots =============================================================
ggsave(box_growth_biomass_png,     p_box_growth_biomass,      width = 7.5, height = 4.8, dpi = 300)
ggsave(box_resp_biomass_png,       p_box_resp_biomass,        width = 7.5, height = 4.8, dpi = 300)

ggsave(bayes_ratio_png,            p_ratio_bayes,             width = 7.5, height = 4.8, dpi = 300)
ggsave(bayes_growth_biomass_png,   p_growth_biomass_bayes,    width = 7.5, height = 4.8, dpi = 300)
ggsave(bayes_resp_biomass_png,     p_resp_biomass_bayes,      width = 7.5, height = 4.8, dpi = 300)

ggsave(bayes_tpc_growth_biomass_png, p_growth_tpc_biomass_bayes, width = 7.5, height = 4.8, dpi = 300)
ggsave(bayes_tpc_resp_biomass_png,   p_resp_tpc_biomass_bayes,   width = 7.5, height = 4.8, dpi = 300)

ggsave(bayes_boltz_ratio_png,      p_boltz_ratio,             width = 7.5, height = 4.8, dpi = 300)
ggsave(bayes_boltz_cue_png,        p_boltz_cue,               width = 7.5, height = 4.8, dpi = 300)
ggsave(bayes_boltz_growth_biomass_png, p_boltz_growth_biomass, width = 7.5, height = 4.8, dpi = 300)
ggsave(bayes_boltz_resp_biomass_png,   p_boltz_resp_biomass,   width = 7.5, height = 4.8, dpi = 300)

if (!is.null(p_postdist_growth_biomass)) {
  ggsave(postdist_tpc_growth_biomass_png, p_postdist_growth_biomass, width = 8.5, height = 6.5, dpi = 300)
}
if (!is.null(p_postdist_resp_biomass)) {
  ggsave(postdist_tpc_resp_biomass_png, p_postdist_resp_biomass, width = 8.5, height = 6.5, dpi = 300)
}

if (!is.null(p_tpc_check_growth_biomass_log)) {
  ggsave(tpc_check_growth_biomass_log_png, p_tpc_check_growth_biomass_log, width = 9.0, height = 6.2, dpi = 300)
}
if (!is.null(p_tpc_check_resp_biomass_log)) {
  ggsave(tpc_check_resp_biomass_log_png, p_tpc_check_resp_biomass_log, width = 9.0, height = 6.2, dpi = 300)
}

if (!is.null(p_temp_effect)) {
  ggsave(fungicide_temp_effect_png, p_temp_effect, width = 8.8, height = 6.4, dpi = 300)
}

# ===== Final messages =========================================================
message("Done.")
message("Stage 1 oxygen model kept separate: O2(t) = O2_0 + (K/r) * (1 - exp(r*t))")
message("delta_Ninoc_to_N0_min read from: ", TRIM_META_CSV)
message("Single global N_inoculation_cells_per_L used for all groups: ", format(N_inoculation_cells_per_L, scientific = FALSE))
message("Per-sample N0 computed using fitted r and sample-specific delta.")
message("Zymoseptoria cell volume (um^3): ", signif(CELL_VOLUME_UM3, 6))
message("Zymoseptoria carbon per cell (fg C): ", signif(CELL_CARBON_FG_PER_CELL, 6))
message("Derived results with carbon units saved to: ", derived_csv)
message("Per-series fit PDF: ", pdf_path)
message("Doses included: ", paste(allowed_doses, collapse = ", "))
message("Replication summary saved to: ", replication_summary_csv)
message("Bayesian thermal models use all replicate-level rates from Stage 1; no averaging by condition was applied.")

if (!is.null(fit_growth_brm)) {
  message("Saved Bayesian growth Arrhenius model: ", brms_growth_rds)
  message("Saved Bayesian growth posterior summaries: ", summary_growth_csv)
  message("Saved Bayesian growth posterior draws: ", post_growth_csv)
}
if (!is.null(fit_resp_brm)) {
  message("Saved Bayesian respiration Arrhenius model: ", brms_resp_rds)
  message("Saved Bayesian respiration posterior summaries: ", summary_resp_csv)
  message("Saved Bayesian respiration posterior draws: ", post_resp_csv)
}
if (!is.null(fit_ratio_brm)) {
  message("Saved Bayesian Arrhenius model for log(resp/growth): ", brms_ratio_rds)
  message("Saved Bayesian ratio posterior summaries: ", summary_ratio_csv)
  message("Saved Bayesian ratio posterior draws: ", post_ratio_csv)
}
if (!is.null(fit_cue_brm)) {
  message("Saved Bayesian Arrhenius model for log(CUE) with common slope: ", brms_cue_rds)
  message("Saved Bayesian CUE posterior summaries: ", summary_cue_csv)
  message("Saved Bayesian CUE posterior draws: ", post_cue_csv)
  message("Saved Bayesian Boltzmann CUE plot: ", bayes_boltz_cue_png)
}
if (!is.null(fit_growth_biomass_brm)) {
  message("Saved Bayesian biomass-corrected growth Arrhenius model: ", brms_growth_biomass_rds)
  message("Saved Bayesian biomass-corrected growth posterior summaries: ", summary_growth_biomass_csv)
  message("Saved Bayesian biomass-corrected growth posterior draws: ", post_growth_biomass_csv)
}
if (!is.null(fit_resp_biomass_brm)) {
  message("Saved Bayesian biomass-corrected respiration Arrhenius model: ", brms_resp_biomass_rds)
  message("Saved Bayesian biomass-corrected respiration posterior summaries: ", summary_resp_biomass_csv)
  message("Saved Bayesian biomass-corrected respiration posterior draws: ", post_resp_biomass_csv)
}
if (!is.null(fit_growth_tpc_brm)) {
  message("Saved Bayesian full Sharpe–Schoolfield growth model: ", brms_tpc_growth_rds)
  message("Saved Bayesian full Sharpe–Schoolfield growth posterior summaries: ", summary_tpc_growth_csv)
  message("Saved Bayesian full Sharpe–Schoolfield growth posterior draws: ", post_tpc_growth_csv)
}
if (!is.null(fit_resp_tpc_brm)) {
  message("Saved Bayesian full Sharpe–Schoolfield respiration model: ", brms_tpc_resp_rds)
  message("Saved Bayesian full Sharpe–Schoolfield respiration posterior summaries: ", summary_tpc_resp_csv)
  message("Saved Bayesian full Sharpe–Schoolfield respiration posterior draws: ", post_tpc_resp_csv)
}
if (!is.null(fit_growth_tpc_biomass_brm)) {
  message("Saved Bayesian full Sharpe–Schoolfield biomass-corrected growth model: ", brms_tpc_growth_biomass_rds)
  message("Saved Bayesian full Sharpe–Schoolfield biomass-corrected growth posterior summaries: ", summary_tpc_growth_biomass_csv)
  message("Saved Bayesian full Sharpe–Schoolfield biomass-corrected growth posterior draws: ", post_tpc_growth_biomass_csv)
}
if (!is.null(fit_resp_tpc_biomass_brm)) {
  message("Saved Bayesian full Sharpe–Schoolfield biomass-corrected respiration model: ", brms_tpc_resp_biomass_rds)
  message("Saved Bayesian full Sharpe–Schoolfield biomass-corrected respiration posterior summaries: ", summary_tpc_resp_biomass_csv)
  message("Saved Bayesian full Sharpe–Schoolfield biomass-corrected respiration posterior draws: ", post_tpc_resp_biomass_csv)
}

if (!is.null(fit_growth_hier_brm)) {
  message("Saved separate hierarchical Arrhenius growth model: ", brms_growth_hier_rds)
}
if (!is.null(fit_resp_hier_brm)) {
  message("Saved separate hierarchical Arrhenius respiration model: ", brms_resp_hier_rds)
}
if (!is.null(fit_ratio_hier_brm)) {
  message("Saved separate hierarchical Arrhenius ratio model: ", brms_ratio_hier_rds)
}
if (!is.null(fit_growth_biomass_hier_brm)) {
  message("Saved separate hierarchical Arrhenius biomass-corrected growth model: ", brms_growth_biomass_hier_rds)
}
if (!is.null(fit_resp_biomass_hier_brm)) {
  message("Saved separate hierarchical Arrhenius biomass-corrected respiration model: ", brms_resp_biomass_hier_rds)
}

if (file.exists(fungicide_temp_effect_summary_csv)) {
  message("Saved temperature-specific fungicide effect summary: ", fungicide_temp_effect_summary_csv)
}
if (file.exists(fungicide_temp_effect_draws_csv)) {
  message("Saved temperature-specific fungicide effect draws: ", fungicide_temp_effect_draws_csv)
}
if (file.exists(fungicide_temp_effect_png)) {
  message("Saved temperature-specific fungicide effect plot: ", fungicide_temp_effect_png)
}
if (file.exists(fungicide_slope_test_summary_csv)) {
  message("Saved Bayesian slope/intercept summary for Growth and CUE effect-size curves: ", fungicide_slope_test_summary_csv)
}
if (file.exists(fungicide_slope_test_draws_csv)) {
  message("Saved Bayesian slope/intercept draws for Growth and CUE effect-size curves: ", fungicide_slope_test_draws_csv)
}

# ─────────────────────────────────────────────────────────────────────────────
# CLEAN ADD-ON: HIGH-LEVEL REGION SYNERGY / ANTAGONISM PLOT
# ─────────────────────────────────────────────────────────────────────────────

OUT_SUBDIR <- file.path(
  figures_dir, "clean_pooled_region_synergy_antagonism_plot"
)
dir.create(OUT_SUBDIR, showWarnings = FALSE, recursive = TRUE)

plot_summary_csv <- file.path(
  tables_dir, "pooled_region_synergy_antagonism_summary.csv"
)
plot_draws_csv   <- file.path(
  tables_dir, "pooled_region_synergy_antagonism_draws.csv"
)
plot_png         <- file.path(
  OUT_SUBDIR, "pooled_region_synergy_antagonism_plot.png"
)
plot_pdf         <- file.path(
  OUT_SUBDIR, "pooled_region_synergy_antagonism_plot.pdf"
)

TEMP_GRID_N <- 200
TEMP_REGION_BREAKS <- c(-Inf, 19, 24, Inf)
TEMP_REGION_LABELS <- c("Cool (15–19°C)", "Mid (20–24°C)", "Hot (25–28°C)")
MAX_DRAWS_PER_DOSE <- 4000
CUE_EPS <- 1e-8
REGION_POOLING_MODE <- "observed"
BASE_SIZE <- 12

has_region_inputs <- file.exists(derived_csv) &&
  file.exists(dose_key_csv) &&
  file.exists(post_tpc_growth_csv) &&
  file.exists(post_tpc_resp_csv) &&
  file.exists(post_ratio_csv)

if (has_region_inputs) {
  results_obs <- readr::read_csv(derived_csv, show_col_types = FALSE) %>%
    mutate(
      T = as.numeric(T),
      Dose = as.character(Dose)
    )
  
  dose_key_tbl2 <- readr::read_csv(dose_key_csv, show_col_types = FALSE) %>%
    mutate(Dose = as.character(Dose))
  
  dose_levels2 <- dose_key_tbl2$Dose
  treated_doses <- setdiff(dose_levels2, DOSE_REF)
  
  if (!DOSE_REF %in% dose_levels2) {
    stop("Reference dose not found: ", DOSE_REF)
  }
  
  observed_temps <- sort(unique(results_obs$T))
  T_min_obs <- min(observed_temps, na.rm = TRUE)
  T_max_obs <- max(observed_temps, na.rm = TRUE)
  
  if (REGION_POOLING_MODE == "observed") {
    T_pred <- observed_temps
  } else if (REGION_POOLING_MODE == "grid") {
    T_pred <- seq(T_min_obs, T_max_obs, length.out = TEMP_GRID_N)
  } else {
    stop("REGION_POOLING_MODE must be either 'observed' or 'grid'.")
  }
  
  post_growth <- readr::read_csv(post_tpc_growth_csv, show_col_types = FALSE) %>%
    mutate(Dose = as.character(Dose)) %>%
    downsample_draws(MAX_DRAWS_PER_DOSE)
  
  post_resp <- readr::read_csv(post_tpc_resp_csv, show_col_types = FALSE) %>%
    mutate(Dose = as.character(Dose)) %>%
    downsample_draws(MAX_DRAWS_PER_DOSE)
  
  post_ratio <- readr::read_csv(post_ratio_csv, show_col_types = FALSE) %>%
    mutate(Dose = as.character(Dose)) %>%
    downsample_draws(MAX_DRAWS_PER_DOSE)
  
  required_growth_cols <- c("Dose", ".draw", "lnB0", "E", "Eh", "Th")
  required_resp_cols   <- c("Dose", ".draw", "lnB0", "E", "Eh", "Th")
  required_ratio_cols  <- c("Dose", ".draw", "alpha", "E")
  
  if (!all(required_growth_cols %in% names(post_growth))) stop("Growth posterior file is missing required columns.")
  if (!all(required_resp_cols %in% names(post_resp))) stop("Respiration posterior file is missing required columns.")
  if (!all(required_ratio_cols %in% names(post_ratio))) stop("Ratio posterior file is missing required columns.")
  
  growth_pred2 <- crossing(T = T_pred, Dose = dose_levels2) %>%
    left_join(post_growth, by = "Dose", relationship = "many-to-many") %>%
    mutate(
      TK = T + 273.15,
      pred_log = pred_ss_log(TK, lnB0, E, Eh, Th),
      pred_raw = exp(pred_log)
    ) %>%
    select(Dose, .draw, T, pred_log, pred_raw)
  
  resp_pred2 <- crossing(T = T_pred, Dose = dose_levels2) %>%
    left_join(post_resp, by = "Dose", relationship = "many-to-many") %>%
    mutate(
      TK = T + 273.15,
      pred_log = pred_ss_log(TK, lnB0, E, Eh, Th),
      pred_raw = exp(pred_log)
    ) %>%
    select(Dose, .draw, T, pred_log, pred_raw)
  
  ratio_pred2 <- crossing(T = T_pred, Dose = dose_levels2) %>%
    left_join(post_ratio, by = "Dose", relationship = "many-to-many") %>%
    mutate(
      TK = T + 273.15,
      pred_log = pred_arr_log(TK, alpha, E),
      pred_raw = exp(pred_log)
    ) %>%
    select(Dose, .draw, T, pred_log, pred_raw)
  
  cue_pred2 <- growth_pred2 %>%
    rename(growth_log = pred_log, growth_raw = pred_raw) %>%
    inner_join(
      resp_pred2 %>% rename(resp_log = pred_log, resp_raw = pred_raw),
      by = c("Dose", ".draw", "T")
    ) %>%
    mutate(
      CUE_raw = growth_raw / (growth_raw + resp_raw),
      CUE_logit = safe_logit(CUE_raw, CUE_EPS)
    ) %>%
    select(Dose, .draw, T, CUE_raw, CUE_logit)
  
  growth_eff2 <- compute_vs_control_effect(growth_pred2, "pred_log", DOSE_REF) %>%
    filter(Dose %in% treated_doses) %>%
    mutate(Trait = "Growth")
  
  resp_eff2 <- compute_vs_control_effect(resp_pred2, "pred_log", DOSE_REF) %>%
    filter(Dose %in% treated_doses) %>%
    mutate(Trait = "Respiration")
  
  ratio_eff2 <- compute_vs_control_effect(ratio_pred2, "pred_log", DOSE_REF) %>%
    filter(Dose %in% treated_doses) %>%
    mutate(Trait = "Respiration / Growth")
  
  cue_eff2 <- compute_vs_control_effect(cue_pred2, "CUE_logit", DOSE_REF) %>%
    filter(Dose %in% treated_doses) %>%
    mutate(Trait = "CUE")
  
  effect_all <- bind_rows(growth_eff2, resp_eff2, ratio_eff2, cue_eff2) %>%
    mutate(
      temp_region = cut(
        T,
        breaks = TEMP_REGION_BREAKS,
        labels = TEMP_REGION_LABELS,
        include.lowest = TRUE,
        right = TRUE
      ),
      temp_region = factor(temp_region, levels = TEMP_REGION_LABELS)
    )
  
  growth_dev <- effect_all %>%
    filter(Trait == "Growth") %>%
    compute_region_deviation() %>%
    mutate(Trait = "Growth")
  
  resp_dev <- effect_all %>%
    filter(Trait == "Respiration") %>%
    compute_region_deviation() %>%
    mutate(Trait = "Respiration")
  
  ratio_dev <- effect_all %>%
    filter(Trait == "Respiration / Growth") %>%
    compute_region_deviation() %>%
    mutate(Trait = "Respiration / Growth")
  
  cue_dev <- effect_all %>%
    filter(Trait == "CUE") %>%
    compute_region_deviation() %>%
    mutate(Trait = "CUE")
  
  deviation_all <- bind_rows(growth_dev, resp_dev, ratio_dev, cue_dev)
  readr::write_csv(deviation_all, plot_draws_csv)
  
  sum_growth <- deviation_all %>%
    filter(Trait == "Growth") %>%
    summarise_region_deviation("Growth", effect_type = "negative")
  
  sum_resp <- deviation_all %>%
    filter(Trait == "Respiration") %>%
    summarise_region_deviation("Respiration", effect_type = "positive")
  
  sum_ratio <- deviation_all %>%
    filter(Trait == "Respiration / Growth") %>%
    summarise_region_deviation("Respiration / Growth", effect_type = "positive")
  
  sum_cue <- deviation_all %>%
    filter(Trait == "CUE") %>%
    summarise_region_deviation("CUE", effect_type = "negative")
  
  plot_df <- bind_rows(sum_growth, sum_resp, sum_ratio, sum_cue) %>%
    mutate(
      temp_region = factor(temp_region, levels = TEMP_REGION_LABELS),
      Trait = factor(Trait, levels = c("Growth", "Respiration", "Respiration / Growth", "CUE"))
    )
  
  plot_df$direction_label <- mapply(
    FUN = label_direction,
    trait = plot_df$Trait,
    q2.5 = plot_df$q2.5,
    q97.5 = plot_df$q97.5,
    SIMPLIFY = TRUE,
    USE.NAMES = FALSE
  )
  
  plot_df <- plot_df %>%
    mutate(
      point_fill = case_when(
        direction_label == "Synergistic" ~ "black",
        direction_label == "Antagonistic" ~ "white",
        TRUE ~ "grey70"
      )
    )
  
  readr::write_csv(plot_df, plot_summary_csv)
  
  p_clean <- ggplot(plot_df, aes(x = temp_region, y = median_deviation)) +
    geom_hline(yintercept = 0, linetype = 2, colour = "grey40", linewidth = 0.5) +
    geom_linerange(aes(ymin = q2.5, ymax = q97.5), linewidth = 0.7, colour = "black") +
    geom_point(aes(fill = point_fill), shape = 21, size = 3.6, stroke = 0.8, colour = "black") +
    scale_fill_identity() +
    facet_wrap(~Trait, scales = "free_y", ncol = 2) +
    labs(
      title = "Temperature-region synergy / antagonism",
      subtitle = make_plot_subtitle(REGION_POOLING_MODE),
      x = NULL,
      y = "Regional deviation from mean fungicide effect"
    ) +
    theme_classic(base_size = BASE_SIZE) +
    theme(
      strip.background = element_blank(),
      strip.text = element_text(face = "bold"),
      axis.text.x = element_text(angle = 20, hjust = 1),
      plot.title = element_text(face = "bold"),
      panel.spacing = unit(1.1, "lines")
    )
  
  ggsave(plot_png, p_clean, width = 8.6, height = 6.8, dpi = 300)
  
  message("Saved region synergy/antagonism summary table: ", plot_summary_csv)
  message("Saved region synergy/antagonism draw-level table: ", plot_draws_csv)
  message("Saved region synergy/antagonism PNG: ", plot_png)
  message("Saved region synergy/antagonism PDF: ", plot_pdf)
} else {
  message("Region synergy/antagonism add-on skipped because required posterior/input files were not available.")
}


if (exists("p_cue_arrhenius_bayes") && !is.null(p_cue_arrhenius_bayes)) {
  ggsave(
    file.path(figures_dir,
              "bayesian_arrhenius_CUE_common_slope_vs_temperature.png"),
    p_cue_arrhenius_bayes,
    width = 7.5, height = 4.8, dpi = 300
  )
}


# ─────────────────────────────────────────────────────────────────────────────
# ADD-ON: TEST WHETHER GROWTH AND CUE EFFECT-SIZE SLOPES DIFFER FROM ZERO
# For the "fungicide effect sizes by temperature" plot
# ─────────────────────────────────────────────────────────────────────────────

growth_cue_slope_summary_csv <- file.path(
  tables_dir,
  "growth_resp_cue_effect_size_slope_test_summary.csv"
)

growth_cue_slope_draws_csv <- file.path(
  tables_dir,
  "growth_resp_cue_effect_size_slope_test_draws.csv"
)

growth_cue_slope_plot_png <- file.path(
  figures_dir,
  "growth_resp_cue_effect_size_slope_test_plot.png"
)

# helper: fit slope within each posterior draw
fit_posterior_slopes <- function(df, trait_name) {
  df %>%
    filter(is.finite(T), is.finite(effect)) %>%
    group_by(Dose, .draw) %>%
    group_modify(~{
      if (nrow(.x) < 2 || dplyr::n_distinct(.x$T) < 2) {
        return(tibble(intercept = NA_real_, slope = NA_real_))
      }
      m <- lm(effect ~ T, data = .x)
      tibble(
        intercept = unname(coef(m)[1]),
        slope = unname(coef(m)[2])
      )
    }) %>%
    ungroup() %>%
    mutate(Trait = trait_name)
}

# run only for the effect-size curves shown in the plot
if (exists("growth_eff") && exists("resp_eff") && exists("cue_eff")) {
  
  slope_draws_growth <- fit_posterior_slopes(growth_eff, "Growth")
  slope_draws_resp   <- fit_posterior_slopes(resp_eff, "Respiration")
  slope_draws_cue    <- fit_posterior_slopes(cue_eff, "CUE")
  
  slope_draws_gc <- bind_rows(slope_draws_growth, slope_draws_resp, slope_draws_cue) %>%
    filter(is.finite(slope)) %>%
    mutate(
      Trait = factor(Trait, levels = c("Growth", "Respiration", "CUE")),
      Dose = factor(Dose, levels = setdiff(dose_levels, "Control"))
    ) %>%
    arrange(Trait, Dose, .draw)
  
  readr::write_csv(slope_draws_gc, growth_cue_slope_draws_csv)
  
  slope_summary_gc <- slope_draws_gc %>%
    group_by(Trait, Dose) %>%
    summarise(
      intercept_median = median(intercept, na.rm = TRUE),
      intercept_q2.5   = quantile(intercept, 0.025, na.rm = TRUE),
      intercept_q97.5  = quantile(intercept, 0.975, na.rm = TRUE),
      slope_median     = median(slope, na.rm = TRUE),
      slope_q2.5       = quantile(slope, 0.025, na.rm = TRUE),
      slope_q97.5      = quantile(slope, 0.975, na.rm = TRUE),
      p_slope_gt0      = mean(slope > 0, na.rm = TRUE),
      p_slope_lt0      = mean(slope < 0, na.rm = TRUE),
      slope_test = case_when(
        slope_q2.5 > 0  ~ "Slope > 0",
        slope_q97.5 < 0 ~ "Slope < 0",
        TRUE            ~ "Slope overlaps 0"
      ),
      .groups = "drop"
    ) %>%
    arrange(Trait, Dose)
  
  readr::write_csv(slope_summary_gc, growth_cue_slope_summary_csv)
  
  p_slope_gc <- ggplot(
    slope_summary_gc,
    aes(x = Dose, y = slope_median, colour = Dose)
  ) +
    geom_hline(yintercept = 0, linetype = 2, colour = "grey40", linewidth = 0.5) +
    geom_linerange(
      aes(ymin = slope_q2.5, ymax = slope_q97.5),
      linewidth = 1
    ) +
    geom_point(size = 3) +
    facet_wrap(~Trait, scales = "free_y", ncol = 3) +
    scale_color_manual(values = dose_cols, limits = dose_levels, drop = FALSE) +
    labs(
      title = "Posterior slope test for Growth, Respiration, and CUE effect-size curves",
      subtitle = "Slope from posterior effect size ~ temperature. Dashed line = zero slope.",
      x = "Dose",
      y = "Posterior slope"
    ) +
    theme_classic(12) +
    theme(legend.position = "none")
  
  ggsave(growth_cue_slope_plot_png, p_slope_gc, width = 12.0, height = 4.8, dpi = 300)
  
  message("Saved Growth/Respiration/CUE slope-test draws: ", growth_cue_slope_draws_csv)
  message("Saved Growth/Respiration/CUE slope-test summary: ", growth_cue_slope_summary_csv)
  message("Saved Growth/Respiration/CUE slope-test plot: ", growth_cue_slope_plot_png)
  
} else {
  message("Growth/Respiration/CUE slope test skipped because growth_eff, resp_eff, and/or cue_eff were not found.")
}













# ─────────────────────────────────────────────────────────────────────────────
# ADD-ON: COMMON TEMPERATURE LINE FOR GROWTH AND CUE EFFECT SIZES
# One pooled line per facet across treated doses
# ─────────────────────────────────────────────────────────────────────────────

common_line_draws_csv <- file.path(
  tables_dir,
  "growth_resp_cue_common_temperature_line_draws.csv"
)

common_line_summary_csv <- file.path(
  tables_dir,
  "growth_resp_cue_common_temperature_line_summary.csv"
)

common_line_plot_png <- file.path(
  figures_dir,
  "growth_resp_cue_common_temperature_line_plot.png"
)

common_line_plot_pdf <- file.path(
  figures_dir,
  "growth_resp_cue_common_temperature_line_plot.pdf"
)

if (exists("growth_eff") && exists("resp_eff") && exists("cue_eff")) {
  
  # pool across treated doses within each posterior draw and temperature
  common_growth <- growth_eff %>%
    filter(Dose != DOSE_REF) %>%
    group_by(.draw, T) %>%
    summarise(
      effect = mean(effect, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(Trait = "Growth")
  
  common_resp <- resp_eff %>%
    filter(Dose != DOSE_REF) %>%
    group_by(.draw, T) %>%
    summarise(
      effect = mean(effect, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(Trait = "Respiration")
  
  common_cue <- cue_eff %>%
    filter(Dose != DOSE_REF) %>%
    group_by(.draw, T) %>%
    summarise(
      effect = mean(effect, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(Trait = "CUE")
  
  common_eff <- bind_rows(common_growth, common_resp, common_cue) %>%
    mutate(Trait = factor(Trait, levels = c("Growth", "Respiration", "CUE")))
  
  # fit one linear trend per posterior draw
  common_line_draws <- common_eff %>%
    group_by(Trait, .draw) %>%
    group_modify(~{
      dd <- .x %>% filter(is.finite(T), is.finite(effect))
      if (nrow(dd) < 2 || dplyr::n_distinct(dd$T) < 2) {
        return(tibble(intercept = NA_real_, slope = NA_real_))
      }
      m <- lm(effect ~ T, data = dd)
      tibble(
        intercept = unname(coef(m)[1]),
        slope = unname(coef(m)[2])
      )
    }) %>%
    ungroup()
  
  readr::write_csv(common_line_draws, common_line_draws_csv)
  
  common_line_summary <- common_line_draws %>%
    group_by(Trait) %>%
    summarise(
      intercept_median = median(intercept, na.rm = TRUE),
      intercept_q2.5   = quantile(intercept, 0.025, na.rm = TRUE),
      intercept_q97.5  = quantile(intercept, 0.975, na.rm = TRUE),
      slope_median     = median(slope, na.rm = TRUE),
      slope_q2.5       = quantile(slope, 0.025, na.rm = TRUE),
      slope_q97.5      = quantile(slope, 0.975, na.rm = TRUE),
      p_slope_gt0      = mean(slope > 0, na.rm = TRUE),
      p_slope_lt0      = mean(slope < 0, na.rm = TRUE),
      slope_direction  = case_when(
        slope_q2.5 > 0  ~ "Positive",
        slope_q97.5 < 0 ~ "Negative",
        TRUE            ~ "Uncertain"
      ),
      .groups = "drop"
    )
  
  readr::write_csv(common_line_summary, common_line_summary_csv)
  
  # make posterior line/ribbon from fitted common slopes
  T_seq <- sort(unique(common_eff$T))
  
  common_pred <- common_line_draws %>%
    filter(is.finite(intercept), is.finite(slope)) %>%
    tidyr::crossing(T = T_seq) %>%
    mutate(pred = intercept + slope * T)
  
  common_pred_sum <- common_pred %>%
    group_by(Trait, T) %>%
    summarise(
      q2.5 = quantile(pred, 0.025, na.rm = TRUE),
      q50  = quantile(pred, 0.5, na.rm = TRUE),
      q97.5 = quantile(pred, 0.975, na.rm = TRUE),
      .groups = "drop"
    )
  
  # raw pooled points for display
  common_points <- common_eff %>%
    group_by(Trait, T) %>%
    summarise(
      median_effect = median(effect, na.rm = TRUE),
      q2.5 = quantile(effect, 0.025, na.rm = TRUE),
      q97.5 = quantile(effect, 0.975, na.rm = TRUE),
      .groups = "drop"
    )
  
  label_df_common <- common_line_summary %>%
    mutate(
      label = paste0(
        "slope = ", signif(slope_median, 3),
        "\n95% CrI [", signif(slope_q2.5, 3), ", ", signif(slope_q97.5, 3), "]",
        "\nP(>0) = ", sprintf("%.3f", p_slope_gt0),
        "\nP(<0) = ", sprintf("%.3f", p_slope_lt0),
        "\n", slope_direction
      )
    )
  
  p_common_line <- ggplot() +
    geom_hline(yintercept = 0, linetype = 2, colour = "grey40", linewidth = 0.5) +
    geom_ribbon(
      data = common_pred_sum,
      aes(x = T, ymin = q2.5, ymax = q97.5),
      fill = "grey70",
      alpha = 0.35
    ) +
    geom_line(
      data = common_pred_sum,
      aes(x = T, y = q50),
      linewidth = 1.1,
      colour = "black"
    ) +
    geom_point(
      data = common_points,
      aes(x = T, y = median_effect),
      size = 2.2,
      colour = "black"
    ) +
    facet_wrap(~Trait, scales = "free_y", ncol = 3) +
    labs(
      title = "Common temperature trend for pooled fungicide effect sizes",
      subtitle = paste0(
        "One pooled line per facet across treated doses only. ",
        "Grey band = 95% credible interval; points = posterior median pooled effect at each temperature."
      ),
      x = "Temperature (°C)",
      y = "Pooled effect size vs Control"
    ) +
    theme_classic(12) +
    geom_text(
      data = label_df_common,
      aes(x = Inf, y = Inf, label = label),
      inherit.aes = FALSE,
      hjust = 1.05,
      vjust = 1.1,
      size = 3.2
    )
  
  ggsave(common_line_plot_png, p_common_line, width = 12.0, height = 4.8, dpi = 300)
  
  message("Saved common-line posterior draws: ", common_line_draws_csv)
  message("Saved common-line summary: ", common_line_summary_csv)
  message("Saved common-line plot PNG: ", common_line_plot_png)
  message("Saved common-line plot PDF: ", common_line_plot_pdf)
  
} else {
  message("Common-line figure skipped because growth_eff, resp_eff, and/or cue_eff were not found.")
}




# ─────────────────────────────────────────────────────────────────────────────
# CLEAN ADD-ONS
# A) Generic posterior-difference helpers
# B) Pairwise posterior differences for all saved posterior files
# C) Half-eye plot: SS biomass-corrected Growth + Respiration vs Control
# D) Dose/Temperature vs Growth/Respiration plots + regressions
# ─────────────────────────────────────────────────────────────────────────────

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggdist)
  library(grid)
})

# =============================================================================
# A) GENERIC HELPERS
# =============================================================================

make_dir <- function(path) {
  dir.create(path, showWarnings = FALSE, recursive = TRUE)
  path
}

make_all_dose_pairs <- function(dose_levels) {
  if (length(dose_levels) < 2) return(tibble())
  
  combn(dose_levels, 2, simplify = FALSE) %>%
    purrr::map_dfr(~tibble(
      Dose_1 = .x[1],
      Dose_2 = .x[2],
      comparison = paste0(.x[2], " - ", .x[1])
    ))
}

read_posterior_csv_checked <- function(path, required_cols, label = NULL) {
  if (!file.exists(path)) {
    message("Skipped: file not found -> ", path)
    return(NULL)
  }
  
  x <- readr::read_csv(path, show_col_types = FALSE)
  missing_cols <- setdiff(required_cols, names(x))
  
  if (length(missing_cols) > 0) {
    message(
      "Skipped ", ifelse(is.null(label), "", paste0(label, " ")),
      "because columns are missing: ",
      paste(missing_cols, collapse = ", ")
    )
    return(NULL)
  }
  
  x
}

compute_pairwise_posterior_differences <- function(post_df, param_cols, dose_levels) {
  stopifnot(all(c("Dose", ".draw") %in% names(post_df)))
  stopifnot(all(param_cols %in% names(post_df)))
  
  dose_pairs <- make_all_dose_pairs(dose_levels)
  if (nrow(dose_pairs) == 0) return(tibble())
  
  post_df <- post_df %>%
    mutate(Dose = as.character(Dose)) %>%
    select(Dose, .draw, all_of(param_cols))
  
  purrr::map_dfr(seq_len(nrow(dose_pairs)), function(i) {
    d1 <- dose_pairs$Dose_1[i]
    d2 <- dose_pairs$Dose_2[i]
    cmp <- dose_pairs$comparison[i]
    
    x1 <- post_df %>%
      filter(Dose == d1) %>%
      rename_with(~paste0(.x, "_1"), all_of(param_cols))
    
    x2 <- post_df %>%
      filter(Dose == d2) %>%
      rename_with(~paste0(.x, "_2"), all_of(param_cols))
    
    joined <- inner_join(x1, x2, by = ".draw")
    if (nrow(joined) == 0) return(tibble())
    
    out <- tibble(
      .draw = joined$.draw,
      Dose_1 = d1,
      Dose_2 = d2,
      comparison = cmp
    )
    
    for (p in param_cols) {
      out[[p]] <- joined[[paste0(p, "_2")]] - joined[[paste0(p, "_1")]]
    }
    
    out
  })
}

summarise_pairwise_posterior_differences <- function(diff_df, param_cols) {
  if (nrow(diff_df) == 0) return(tibble())
  
  diff_df %>%
    pivot_longer(
      cols = all_of(param_cols),
      names_to = "parameter",
      values_to = "difference"
    ) %>%
    group_by(Dose_1, Dose_2, comparison, parameter) %>%
    summarise(
      median_diff = median(difference, na.rm = TRUE),
      mean_diff   = mean(difference, na.rm = TRUE),
      sd_diff     = sd(difference, na.rm = TRUE),
      q2.5        = quantile(difference, 0.025, na.rm = TRUE),
      q17         = quantile(difference, 0.17,  na.rm = TRUE),
      q83         = quantile(difference, 0.83,  na.rm = TRUE),
      q97.5       = quantile(difference, 0.975, na.rm = TRUE),
      p_gt0       = mean(difference > 0, na.rm = TRUE),
      p_lt0       = mean(difference < 0, na.rm = TRUE),
      direction   = case_when(
        q2.5 > 0  ~ "Dose_2 > Dose_1",
        q97.5 < 0 ~ "Dose_2 < Dose_1",
        TRUE      ~ "Uncertain"
      ),
      .groups = "drop"
    ) %>%
    arrange(parameter, Dose_1, Dose_2)
}

plot_pairwise_posterior_differences <- function(summary_df, title_txt) {
  if (nrow(summary_df) == 0) return(NULL)
  
  ggplot(summary_df, aes(x = comparison, y = median_diff)) +
    geom_hline(yintercept = 0, linetype = 2, colour = "grey40", linewidth = 0.5) +
    geom_linerange(aes(ymin = q2.5, ymax = q97.5), linewidth = 0.8) +
    geom_point(size = 2.3) +
    facet_wrap(~parameter, scales = "free_y", ncol = 2) +
    labs(
      title = title_txt,
      subtitle = "Posterior median difference with 95% credible interval. Difference = Dose_2 - Dose_1",
      x = "Dose comparison",
      y = "Posterior difference"
    ) +
    theme_classic(12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

compute_vs_reference_differences <- function(post_df, param_cols, ref_dose, treated_doses, trait_name) {
  stopifnot(all(c("Dose", ".draw") %in% names(post_df)))
  stopifnot(all(param_cols %in% names(post_df)))
  
  post_df <- post_df %>%
    mutate(Dose = as.character(Dose)) %>%
    select(Dose, .draw, all_of(param_cols))
  
  ref_df <- post_df %>%
    filter(Dose == ref_dose) %>%
    rename_with(~paste0(.x, "_ref"), all_of(param_cols))
  
  if (nrow(ref_df) == 0) stop("Reference dose not found: ", ref_dose)
  
  purrr::map_dfr(treated_doses, function(d) {
    trt_df <- post_df %>%
      filter(Dose == d) %>%
      rename_with(~paste0(.x, "_trt"), all_of(param_cols))
    
    joined <- inner_join(trt_df, ref_df, by = ".draw")
    if (nrow(joined) == 0) return(tibble())
    
    out <- tibble(
      Trait = trait_name,
      .draw = joined$.draw,
      Dose = d,
      comparison = paste0(d, " - ", ref_dose)
    )
    
    for (p in param_cols) {
      out[[p]] <- joined[[paste0(p, "_trt")]] - joined[[paste0(p, "_ref")]]
    }
    
    out
  }) %>%
    pivot_longer(
      cols = all_of(param_cols),
      names_to = "parameter",
      values_to = "difference"
    )
}

summarise_reference_differences <- function(draws_df) {
  draws_df %>%
    group_by(Trait, Dose, comparison, parameter) %>%
    summarise(
      median = median(difference, na.rm = TRUE),
      q2.5   = quantile(difference, 0.025, na.rm = TRUE),
      q17    = quantile(difference, 0.17,  na.rm = TRUE),
      q83    = quantile(difference, 0.83,  na.rm = TRUE),
      q97.5  = quantile(difference, 0.975, na.rm = TRUE),
      p_gt0  = mean(difference > 0, na.rm = TRUE),
      p_lt0  = mean(difference < 0, na.rm = TRUE),
      .groups = "drop"
    )
}

extract_lm_summary <- function(dat, formula_obj, response_name, predictor_name) {
  dat2 <- dat %>% filter(stats::complete.cases(model.frame(formula_obj, data = dat)))
  if (nrow(dat2) < 3) return(tibble())
  
  fit <- lm(formula_obj, data = dat2)
  sm  <- summary(fit)
  
  as.data.frame(sm$coefficients) %>%
    tibble::rownames_to_column("term") %>%
    as_tibble() %>%
    setNames(c("term", "estimate", "std_error", "t_value", "p_value")) %>%
    mutate(
      response = response_name,
      predictor = predictor_name,
      n = nrow(dat2),
      r_squared = sm$r.squared,
      adj_r_squared = sm$adj.r.squared,
      aic = AIC(fit),
      bic = BIC(fit)
    ) %>%
    select(response, predictor, term, estimate, std_error, t_value, p_value,
           n, r_squared, adj_r_squared, aic, bic)
}

# =============================================================================
# C) HALF-EYE PLOT: SHARPE–SCHOOLFIELD BIOMASS-CORRECTED VS CONTROL
# =============================================================================

ss_halfeye_dir <- make_dir(
  file.path(figures_dir, "ss_biomass_halfeye_vs_control")
)

ss_halfeye_draws_csv   <- file.path(
  tables_dir, "ss_biomass_halfeye_vs_control_draws.csv"
)
ss_halfeye_summary_csv <- file.path(
  tables_dir, "ss_biomass_halfeye_vs_control_summary.csv"
)
ss_halfeye_plot_png    <- file.path(
  ss_halfeye_dir, "ss_biomass_halfeye_vs_control_plot.png"
)
ss_halfeye_plot_pdf    <- file.path(
  ss_halfeye_dir, "ss_biomass_halfeye_vs_control_plot.pdf"
)

REF_DOSE <- "Control"
SS_PARAMS <- c("lnB0", "E", "Eh", "Th_C")
selected_treated_doses <- NULL

trait_fill_cols <- c("Growth" = "#4C78A8", "Respiration" = "#F58518")
trait_line_cols <- c("Growth" = "#2F5D8A", "Respiration" = "#C96A12")
param_label_map <- c(lnB0 = "ln(B0)", E = "E", Eh = "Eh", Th_C = "Th (°C)")

post_growth_bio <- read_posterior_csv_checked(
  post_tpc_growth_biomass_csv,
  c("Dose", ".draw", SS_PARAMS),
  "SS biomass growth"
)

post_resp_bio <- read_posterior_csv_checked(
  post_tpc_resp_biomass_csv,
  c("Dose", ".draw", SS_PARAMS),
  "SS biomass respiration"
)

if (!is.null(post_growth_bio) && !is.null(post_resp_bio)) {
  treated_doses <- setdiff(dose_levels, REF_DOSE)
  if (!is.null(selected_treated_doses)) {
    treated_doses <- intersect(as.character(selected_treated_doses), treated_doses)
  }
  
  if (length(treated_doses) > 0) {
    diff_growth <- compute_vs_reference_differences(
      post_growth_bio, SS_PARAMS, REF_DOSE, treated_doses, "Growth"
    )
    diff_resp <- compute_vs_reference_differences(
      post_resp_bio, SS_PARAMS, REF_DOSE, treated_doses, "Respiration"
    )
    
    diff_all <- bind_rows(diff_growth, diff_resp)
    
    if (nrow(diff_all) > 0) {
      comparison_levels <- paste0(treated_doses, " - ", REF_DOSE)
      
      diff_all <- diff_all %>%
        mutate(
          Trait = factor(Trait, levels = c("Growth", "Respiration")),
          parameter = factor(parameter, levels = SS_PARAMS),
          comparison = factor(comparison, levels = rev(comparison_levels))
        )
      
      diff_summary <- summarise_reference_differences(diff_all) %>%
        mutate(
          Trait = factor(Trait, levels = c("Growth", "Respiration")),
          parameter = factor(parameter, levels = SS_PARAMS),
          comparison = factor(comparison, levels = rev(comparison_levels))
        )
      
      readr::write_csv(diff_all, ss_halfeye_draws_csv)
      readr::write_csv(diff_summary, ss_halfeye_summary_csv)
      
      p_ss_halfeye <- ggplot(
        diff_all,
        aes(x = difference, y = comparison, fill = Trait, colour = Trait)
      ) +
        geom_vline(xintercept = 0, linetype = 2, linewidth = 0.35, colour = "grey50") +
        ggdist::stat_halfeye(
          adjust = 0.7,
          width = 0.7,
          .width = c(0.66, 0.95),
          justification = -0.2,
          point_interval = median_qi,
          slab_alpha = 0.7,
          interval_size = 0.7,
          point_size = 1.6,
          side = "right"
        ) +
        facet_grid(
          Trait ~ parameter,
          scales = "free_x",
          labeller = labeller(parameter = param_label_map)
        ) +
        scale_fill_manual(values = trait_fill_cols) +
        scale_colour_manual(values = trait_line_cols) +
        labs(
          title = "Posterior differences in Sharpe–Schoolfield parameters",
          subtitle = paste0(
            "Biomass-corrected Growth and Respiration, shown only as treatment - Control contrasts. ",
            "Point = posterior median; thick and thin intervals = 66% and 95% credible intervals."
          ),
          x = "Posterior difference",
          y = NULL,
          fill = NULL,
          colour = NULL
        ) +
        theme_minimal(base_size = 11) +
        theme(
          legend.position = "none",
          strip.text = element_text(face = "bold", size = 11, colour = "black"),
          panel.grid.major.y = element_blank(),
          panel.grid.minor = element_blank(),
          panel.grid.major.x = element_line(colour = "grey88", linewidth = 0.3),
          axis.text.y = element_text(size = 9, colour = "black"),
          axis.text.x = element_text(size = 9, colour = "black"),
          axis.title.x = element_text(size = 11),
          plot.title = element_text(face = "bold", size = 14, colour = "black"),
          plot.subtitle = element_text(size = 10, colour = "black"),
          panel.spacing.x = unit(1, "lines"),
          panel.spacing.y = unit(1, "lines")
        )
      
      ggsave(ss_halfeye_plot_png, p_ss_halfeye, width = 10.5, height = 6.8, dpi = 500, bg = "white")
      
      message("Saved: ", ss_halfeye_draws_csv)
      message("Saved: ", ss_halfeye_summary_csv)
      message("Saved: ", ss_halfeye_plot_png)
      message("Saved: ", ss_halfeye_plot_pdf)
    }
  }
}

# =============================================================================
# D) DOSE/TEMPERATURE VS GROWTH/RESPIRATION PLOTS + REGRESSIONS
# =============================================================================

rate_plot_dir <- make_dir(
  file.path(figures_dir, "dose_temp_vs_rates_plots")
)

rate_files <- list(
  dose_growth_png    = file.path(rate_plot_dir, "dose_vs_growth_rate.png"),
  dose_resp_png      = file.path(
    rate_plot_dir, "dose_vs_respiration_rate.png"
  ),
  temp_growth_png    = file.path(rate_plot_dir, "temp_vs_growth_rate.png"),
  temp_resp_png      = file.path(
    rate_plot_dir, "temp_vs_respiration_rate.png"
  ),
  dose_growth_pdf    = file.path(rate_plot_dir, "dose_vs_growth_rate.pdf"),
  dose_resp_pdf      = file.path(
    rate_plot_dir, "dose_vs_respiration_rate.pdf"
  ),
  temp_growth_pdf    = file.path(rate_plot_dir, "temp_vs_growth_rate.pdf"),
  temp_resp_pdf      = file.path(
    rate_plot_dir, "temp_vs_respiration_rate.pdf"
  ),
  dose_growth_lm_csv = file.path(
    tables_dir, "dose_vs_growth_rate_regression_summary.csv"
  ),
  dose_resp_lm_csv   = file.path(
    tables_dir, "dose_vs_respiration_rate_regression_summary.csv"
  ),
  temp_growth_lm_csv = file.path(
    tables_dir, "temp_vs_growth_rate_regression_summary.csv"
  ),
  temp_resp_lm_csv   = file.path(
    tables_dir, "temp_vs_respiration_rate_regression_summary.csv"
  )
)

plot_dat <- results_plot %>%
  mutate(
    Dose_chr = as.character(Dose),
    Dose_num = suppressWarnings(as.numeric(Dose_chr)),
    Dose_num_for_log = case_when(
      Dose_chr == "Control" ~ 0,
      is.finite(Dose_num) ~ Dose_num,
      TRUE ~ NA_real_
    ),
    log10_dose_plus1 = log10(Dose_num_for_log + 1)
  )

dose_num_levels <- sort(unique(
  plot_dat$Dose_num_for_log[is.finite(plot_dat$Dose_num_for_log)]
))
dose_index_lookup <- setNames(
  seq_along(dose_num_levels) - 1L,
  as.character(dose_num_levels)
)
plot_dat <- plot_dat %>%
  mutate(dose_index = dose_index_lookup[as.character(Dose_num_for_log)])

dose_summary <- plot_dat %>%
  group_by(Dose, Dose_chr, Dose_num_for_log, log10_dose_plus1, dose_index) %>%
  summarise(
    growth_mean = mean(growth_C_per_C_h, na.rm = TRUE),
    resp_mean   = mean(respiration_C_per_C_h, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(Dose = factor(Dose, levels = dose_levels))

temp_summary <- plot_dat %>%
  group_by(T) %>%
  summarise(
    growth_mean = mean(growth_C_per_C_h, na.rm = TRUE),
    resp_mean   = mean(respiration_C_per_C_h, na.rm = TRUE),
    .groups = "drop"
  )

make_rate_plot <- function(dat, xvar, yvar, ylab, title_txt, subtitle_txt,
                           summary_df = NULL, summary_x = NULL, summary_y = NULL,
                           x_breaks = NULL, x_labels = NULL) {
  p <- ggplot(dat, aes(x = .data[[xvar]], y = .data[[yvar]], colour = Dose, fill = Dose)) +
    geom_point(size = 2.2, alpha = 0.9,
               position = if (grepl("dose", xvar, ignore.case = TRUE)) position_jitter(width = 0.02, height = 0) else "identity") +
    geom_smooth(aes(group = 1), method = "loess", se = TRUE, colour = "black", fill = "grey70", linewidth = 0.9) +
    scale_y_log10() +
    scale_colour_manual(values = dose_cols, limits = dose_levels, drop = FALSE) +
    scale_fill_manual(values = dose_cols, limits = dose_levels, drop = FALSE) +
    labs(
      title = title_txt,
      subtitle = subtitle_txt,
      x = ifelse(xvar == "T", "Temperature (°C)", "Dose"),
      y = ylab,
      colour = "Dose",
      fill = "Dose"
    ) +
    theme_classic(12)
  
  if (!is.null(x_breaks) && !is.null(x_labels)) {
    p <- p + scale_x_continuous(breaks = x_breaks, labels = x_labels)
  }
  
  p
}

dose_growth_dat <- plot_dat %>% filter(is.finite(growth_C_per_C_h), growth_C_per_C_h > 0, is.finite(log10_dose_plus1))
dose_resp_dat   <- plot_dat %>% filter(is.finite(respiration_C_per_C_h), respiration_C_per_C_h > 0, is.finite(log10_dose_plus1))
temp_growth_dat <- plot_dat %>% filter(is.finite(T), is.finite(growth_C_per_C_h), growth_C_per_C_h > 0)
temp_resp_dat   <- plot_dat %>% filter(is.finite(T), is.finite(respiration_C_per_C_h), respiration_C_per_C_h > 0)

p_dose_growth <- make_rate_plot(
  dat = dose_growth_dat,
  xvar = "dose_index",
  yvar = "growth_C_per_C_h",
  ylab = "Growth rate (C per C per h, log scale)",
  title_txt = "Dose vs growth rate",
  subtitle_txt = "Points are replicates; black line is loess fit; x-axis equally spaced by dose level.",
  x_breaks = dose_summary$dose_index,
  x_labels = as.character(dose_summary$Dose_num_for_log)
)

p_dose_resp <- make_rate_plot(
  dat = dose_resp_dat,
  xvar = "dose_index",
  yvar = "respiration_C_per_C_h",
  ylab = "Respiration rate (C per C per h, log scale)",
  title_txt = "Dose vs respiration rate",
  subtitle_txt = "Points are replicates; black line is loess fit; x-axis equally spaced by dose level.",
  x_breaks = dose_summary$dose_index,
  x_labels = as.character(dose_summary$Dose_num_for_log)
)

p_temp_growth <- make_rate_plot(
  dat = temp_growth_dat,
  xvar = "T",
  yvar = "growth_C_per_C_h",
  ylab = "Growth rate (C per C per h, log scale)",
  title_txt = "Temperature vs growth rate",
  subtitle_txt = "Points are replicates coloured by dose; black line is global loess fit."
)

p_temp_resp <- make_rate_plot(
  dat = temp_resp_dat,
  xvar = "T",
  yvar = "respiration_C_per_C_h",
  ylab = "Respiration rate (C per C per h, log scale)",
  title_txt = "Temperature vs respiration rate",
  subtitle_txt = "Points are replicates coloured by dose; black line is global loess fit."
)

dose_growth_lm <- extract_lm_summary(dose_growth_dat, log10(growth_C_per_C_h) ~ log10_dose_plus1, "log10(growth_C_per_C_h)", "log10(dose + 1)")
dose_resp_lm   <- extract_lm_summary(dose_resp_dat, log10(respiration_C_per_C_h) ~ log10_dose_plus1, "log10(respiration_C_per_C_h)", "log10(dose + 1)")
temp_growth_lm <- extract_lm_summary(temp_growth_dat, log10(growth_C_per_C_h) ~ T, "log10(growth_C_per_C_h)", "T")
temp_resp_lm   <- extract_lm_summary(temp_resp_dat, log10(respiration_C_per_C_h) ~ T, "log10(respiration_C_per_C_h)", "T")

readr::write_csv(dose_growth_lm, rate_files$dose_growth_lm_csv)
readr::write_csv(dose_resp_lm,   rate_files$dose_resp_lm_csv)
readr::write_csv(temp_growth_lm, rate_files$temp_growth_lm_csv)
readr::write_csv(temp_resp_lm,   rate_files$temp_resp_lm_csv)

ggsave(rate_files$dose_growth_png, p_dose_growth, width = 7.8, height = 5.0, dpi = 300)
ggsave(rate_files$dose_resp_png,   p_dose_resp,   width = 7.8, height = 5.0, dpi = 300)
ggsave(rate_files$temp_growth_png, p_temp_growth, width = 7.8, height = 5.0, dpi = 300)
ggsave(rate_files$temp_resp_png,   p_temp_resp,   width = 7.8, height = 5.0, dpi = 300)

message("Saved: ", rate_files$dose_growth_png)
message("Saved: ", rate_files$dose_resp_png)
message("Saved: ", rate_files$temp_growth_png)
message("Saved: ", rate_files$temp_resp_png)
message("Saved: ", rate_files$dose_growth_lm_csv)
message("Saved: ", rate_files$dose_resp_lm_csv)
message("Saved: ", rate_files$temp_growth_lm_csv)
message("Saved: ", rate_files$temp_resp_lm_csv)

# =============================================================================
# PUBLICATION FIGURES
# =============================================================================

suppressPackageStartupMessages({
  library(patchwork)
  library(mgcv)
})

pub_fig_dir <- make_dir(file.path(figures_dir, "publication_figures"))

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

make_lnB0_halfeye <- function(post_csv) {
  if (!file.exists(post_csv)) return(NULL)
  post_df <- readr::read_csv(post_csv, show_col_types = FALSE) %>%
    mutate(Dose = factor(Dose, levels = rev(dose_levels)))
  ctrl_med <- post_df %>%
    filter(as.character(Dose) == "Control") %>%
    summarise(ref = median(lnB0, na.rm = TRUE)) %>%
    pull(ref)
  ggplot(post_df, aes(x = lnB0, y = Dose, fill = Dose, colour = Dose)) +
    ggdist::stat_halfeye(
      adjust = 0.8, width = 0.65, .width = c(0.66, 0.95),
      point_interval = ggdist::median_qi,
      slab_alpha = 0.65, interval_size = 0.7,
      point_size = 1.5, side = "right", justification = -0.15
    ) +
    geom_vline(xintercept = ctrl_med, linetype = 2,
               colour = "grey40", linewidth = 0.4) +
    scale_fill_manual(values = dose_cols, limits = dose_levels, drop = FALSE) +
    scale_colour_manual(values = dose_cols, limits = dose_levels, drop = FALSE) +
    labs(x = expression(ln(B[0])), y = NULL) +
    theme_classic(11) +
    theme(legend.position = "none", axis.text.y = element_text(size = 9))
}

make_alpha_halfeye <- function(post_csv) {
  if (!file.exists(post_csv)) return(NULL)
  post_df <- readr::read_csv(post_csv, show_col_types = FALSE) %>%
    mutate(Dose = factor(Dose, levels = rev(dose_levels)))
  ctrl_med <- post_df %>%
    filter(as.character(Dose) == "Control") %>%
    summarise(ref = median(alpha, na.rm = TRUE)) %>%
    pull(ref)
  ggplot(post_df, aes(x = alpha, y = Dose, fill = Dose, colour = Dose)) +
    ggdist::stat_halfeye(
      adjust = 0.8, width = 0.65, .width = c(0.66, 0.95),
      point_interval = ggdist::median_qi,
      slab_alpha = 0.65, interval_size = 0.7,
      point_size = 1.5, side = "right", justification = -0.15
    ) +
    geom_vline(xintercept = ctrl_med, linetype = 2,
               colour = "grey40", linewidth = 0.4) +
    scale_fill_manual(values = dose_cols, limits = dose_levels, drop = FALSE) +
    scale_colour_manual(values = dose_cols, limits = dose_levels, drop = FALSE) +
    labs(x = expression(alpha ~ (equiv. ~ ln(B[0]))), y = NULL) +
    theme_classic(11) +
    theme(legend.position = "none", axis.text.y = element_text(size = 9))
}

make_tpc_pub <- function(dat, grid, ylab) {
  p <- ggplot(dat, aes(T, y_raw, colour = Dose, fill = Dose)) +
    geom_point(size = 1.8, alpha = 0.75) +
    {
      if (nrow(grid) > 0)
        geom_ribbon(data = grid,
                    aes(x = T, ymin = raw_q2.5, ymax = raw_q97.5, fill = Dose),
                    inherit.aes = FALSE, alpha = 0.15, colour = NA)
    } +
    {
      if (nrow(grid) > 0)
        geom_line(data = grid,
                  aes(x = T, y = raw_q50, colour = Dose),
                  inherit.aes = FALSE, linewidth = 0.9)
    } +
    scale_y_log10() +
    scale_colour_manual(values = dose_cols, limits = dose_levels, drop = FALSE) +
    scale_fill_manual(values = dose_cols, limits = dose_levels, drop = FALSE) +
    labs(x = "Temperature (\u00b0C)", y = ylab,
         colour = "Dose (mg L\u207b\u00b9)", fill = "Dose (mg L\u207b\u00b9)") +
    theme_classic(11) +
    theme(legend.key.size = unit(0.4, "cm"))
  p
}

make_boltz_cue_pub <- function(dat, grid) {
  p <- ggplot(dat, aes(boltz_shift, y_raw, colour = Dose, fill = Dose)) +
    geom_point(size = 1.8, alpha = 0.75) +
    {
      if (nrow(grid) > 0)
        geom_ribbon(data = grid,
                    aes(x = boltz_shift, ymin = raw_q2.5, ymax = raw_q97.5,
                        fill = Dose),
                    inherit.aes = FALSE, alpha = 0.15, colour = NA)
    } +
    {
      if (nrow(grid) > 0)
        geom_line(data = grid,
                  aes(x = boltz_shift, y = raw_q50, colour = Dose),
                  inherit.aes = FALSE, linewidth = 0.9)
    } +
    scale_y_log10() +
    scale_colour_manual(values = dose_cols, limits = dose_levels, drop = FALSE) +
    scale_fill_manual(values = dose_cols, limits = dose_levels, drop = FALSE) +
    labs(
      x      = expression(frac(1, k[T[ref]]) - frac(1, k * T)),
      y      = "CUE",
      colour = "Dose (mg L\u207b\u00b9)",
      fill   = "Dose (mg L\u207b\u00b9)"
    ) +
    theme_classic(11) +
    theme(legend.key.size = unit(0.4, "cm"))
  p
}

# ---------------------------------------------------------------------------
# Fig. 1: Growth and Respiration TPCs with lnB0 posteriors
# ---------------------------------------------------------------------------

p_tpc_g_pub <- make_tpc_pub(
  growth_biomass_dat, tpc_grid_growth_biomass,
  "Growth rate (C C\u207b\u00b9 h\u207b\u00b9)"
)
p_tpc_r_pub <- make_tpc_pub(
  resp_biomass_dat, tpc_grid_resp_biomass,
  "Respiration rate (C C\u207b\u00b9 h\u207b\u00b9)"
)
p_lnB0_g <- make_lnB0_halfeye(post_tpc_growth_biomass_csv)
p_lnB0_r <- make_lnB0_halfeye(post_tpc_resp_biomass_csv)

if (!is.null(p_lnB0_g) && !is.null(p_lnB0_r)) {
  row_A <- (p_tpc_g_pub + theme(legend.position = "none")) + p_lnB0_g +
    plot_layout(widths = c(2, 1))
  row_B <- (p_tpc_r_pub + theme(legend.position = "none")) + p_lnB0_r +
    plot_layout(widths = c(2, 1))
  fig1 <- (row_A / row_B) +
    plot_annotation(
      tag_levels = "A",
      caption = paste0(
        "Ribbons: 95% posterior credible interval; line: posterior median.",
        "\nOther Sharpe\u2013Schoolfield parameters (E, Eh, Th) show no treatment effect."
      )
    )
  ggsave(
    file.path(pub_fig_dir, "fig1_tpc_lnB0.png"),
    fig1, width = 13, height = 9, dpi = 300, bg = "white"
  )
  message("Saved: ", file.path(pub_fig_dir, "fig1_tpc_lnB0.png"))
}

# ---------------------------------------------------------------------------
# Fig. 2: CUE Boltzmann-Arrhenius + alpha posterior
# ---------------------------------------------------------------------------

p_boltz_cue_pub <- make_boltz_cue_pub(cue_dat, arr_grid_cue)
p_alpha_he      <- make_alpha_halfeye(post_cue_csv)

if (!is.null(p_alpha_he)) {
  fig2 <- (p_boltz_cue_pub + theme(legend.position = "none")) + p_alpha_he +
    plot_layout(widths = c(2, 1)) +
    plot_annotation(
      caption = paste0(
        "CUE modelled as Boltzmann\u2013Arrhenius with a common slope E across doses.",
        "\nOnly \u03b1 (intercept, equiv. ln B\u2080) differs by dose; E shows no treatment effect."
      )
    )
  ggsave(
    file.path(pub_fig_dir, "fig2_cue_alpha.png"),
    fig2, width = 11, height = 5, dpi = 300, bg = "white"
  )
  message("Saved: ", file.path(pub_fig_dir, "fig2_cue_alpha.png"))
}

# ---------------------------------------------------------------------------
# Fig. 3: Effect size vs temperature with GAM smooth
# ---------------------------------------------------------------------------

pub_eff_df <- NULL
if (exists("temp_effect_summary")) {
  pub_eff_df <- temp_effect_summary
} else if (file.exists(fungicide_temp_effect_summary_csv)) {
  pub_eff_df <- readr::read_csv(fungicide_temp_effect_summary_csv,
                                show_col_types = FALSE)
}

if (!is.null(pub_eff_df)) {
  pub_eff_df <- pub_eff_df %>%
    filter(Trait %in% c("Growth", "Respiration", "CUE")) %>%
    mutate(
      Trait = factor(Trait, levels = c("Growth", "Respiration", "CUE")),
      Dose  = factor(as.character(Dose), levels = dose_levels)
    )
  
  gam_preds_list <- pub_eff_df %>%
    group_by(Trait, Dose) %>%
    tidyr::nest() %>%
    mutate(
      preds = purrr::map(data, function(dd) {
        if (nrow(dd) < 4) return(NULL)
        k_val <- min(4L, nrow(dd) - 1L)
        m <- tryCatch(
          mgcv::gam(median_effect ~ s(T, k = k_val, bs = "tp"),
                    data = dd, method = "REML"),
          error = function(e) NULL
        )
        if (is.null(m)) return(NULL)
        T_seq <- seq(min(dd$T), max(dd$T), length.out = 80)
        pred  <- predict(m, newdata = data.frame(T = T_seq), se.fit = TRUE)
        tibble(
          T   = T_seq,
          fit = as.numeric(pred$fit),
          lwr = as.numeric(pred$fit) - 1.96 * as.numeric(pred$se.fit),
          upr = as.numeric(pred$fit) + 1.96 * as.numeric(pred$se.fit)
        )
      })
    )
  
  gam_preds <- gam_preds_list %>%
    filter(!purrr::map_lgl(preds, is.null)) %>%
    select(Trait, Dose, preds) %>%
    tidyr::unnest(preds)
  
  fig3a <- ggplot(pub_eff_df,
                  aes(x = T, y = median_effect, colour = Dose, fill = Dose)) +
    geom_hline(yintercept = 0, linetype = 2, colour = "grey40", linewidth = 0.5) +
    geom_ribbon(aes(ymin = q2.5, ymax = q97.5), alpha = 0.12, colour = NA) +
    geom_point(size = 1.6, alpha = 0.85) +
    {
      if (nrow(gam_preds) > 0)
        geom_ribbon(data = gam_preds,
                    aes(x = T, ymin = lwr, ymax = upr, fill = Dose),
                    inherit.aes = FALSE, alpha = 0.18, colour = NA)
    } +
    {
      if (nrow(gam_preds) > 0)
        geom_line(data = gam_preds,
                  aes(x = T, y = fit, colour = Dose),
                  inherit.aes = FALSE, linewidth = 1.0)
    } +
    facet_wrap(~Trait, scales = "free_y", ncol = 3) +
    scale_colour_manual(values = dose_cols, limits = dose_levels, drop = FALSE) +
    scale_fill_manual(values = dose_cols, limits = dose_levels, drop = FALSE) +
    labs(
      x       = "Temperature (\u00b0C)",
      y       = "Effect size vs control",
      colour  = "Dose (mg L\u207b\u00b9)",
      fill    = "Dose (mg L\u207b\u00b9)",
      caption = "Ribbon (data): 95% posterior CI. Smooth line + ribbon: GAM thin-plate spline (k \u2264 4, REML)."
    ) +
    theme_classic(11) +
    theme(strip.background = element_blank(),
          strip.text = element_text(face = "bold"))
  
  ggsave(
    file.path(pub_fig_dir, "fig3_effect_size_gam.png"),
    fig3a, width = 13, height = 5, dpi = 300, bg = "white"
  )
  message("Saved: ", file.path(pub_fig_dir, "fig3_effect_size_gam.png"))
}

# ---------------------------------------------------------------------------
# Fig. 1 EXTRA: Dose-response curves (dose on x-axis, temperature as color)
# ---------------------------------------------------------------------------

dose_resp_pub_dat <- bind_rows(
  growth_biomass_dat %>%
    mutate(
      Rate     = y_raw,
      Trait    = "Growth rate (C C\u207b\u00b9 h\u207b\u00b9)",
      Dose_num = suppressWarnings(as.numeric(as.character(Dose))),
      Dose_num = ifelse(is.na(Dose_num), 0, Dose_num)
    ),
  resp_biomass_dat %>%
    mutate(
      Rate     = y_raw,
      Trait    = "Respiration rate (C C\u207b\u00b9 h\u207b\u00b9)",
      Dose_num = suppressWarnings(as.numeric(as.character(Dose))),
      Dose_num = ifelse(is.na(Dose_num), 0, Dose_num)
    )
) %>%
  filter(is.finite(Rate), Rate > 0) %>%
  mutate(
    Trait = factor(Trait, levels = c(
      "Growth rate (C C\u207b\u00b9 h\u207b\u00b9)",
      "Respiration rate (C C\u207b\u00b9 h\u207b\u00b9)"
    )),
    T_fac = factor(T)
  )

temp_labels_extra <- sort(unique(dose_resp_pub_dat$T))
n_temps_extra     <- length(temp_labels_extra)
temp_pal_extra    <- grDevices::hcl.colors(n_temps_extra, palette = "Temps")
names(temp_pal_extra) <- as.character(temp_labels_extra)

fig1_extra <- ggplot(dose_resp_pub_dat,
                     aes(x = Dose_num, y = Rate,
                         colour = T_fac, fill = T_fac)) +
  geom_point(size = 1.7, alpha = 0.8) +
  geom_smooth(aes(group = T_fac), method = "gam",
              formula = y ~ s(x, k = 4, bs = "tp"),
              se = TRUE, alpha = 0.15, linewidth = 0.8) +
  scale_y_log10() +
  scale_colour_manual(values = temp_pal_extra,
                      labels = paste0(temp_labels_extra, "\u00b0C")) +
  scale_fill_manual(values = temp_pal_extra,
                    labels = paste0(temp_labels_extra, "\u00b0C")) +
  facet_wrap(~Trait, scales = "free_y", ncol = 2) +
  labs(
    x      = "Prothioconazole concentration (mg L\u207b\u00b9)",
    y      = "Rate (log scale)",
    colour = "Temperature",
    fill   = "Temperature"
  ) +
  theme_classic(11) +
  theme(strip.background = element_blank(),
        strip.text = element_text(face = "bold"))

ggsave(
  file.path(pub_fig_dir, "fig1_extra_dose_response_by_temp.png"),
  fig1_extra, width = 12, height = 5.5, dpi = 300, bg = "white"
)
message("Saved: ", file.path(pub_fig_dir, "fig1_extra_dose_response_by_temp.png"))

# Both axes log scale (log10-transformed values shown directly on axes)
dose_resp_pub_dat_log <- dose_resp_pub_dat %>%
  filter(Dose_num > 0) %>%
  mutate(
    log10_Dose = log10(Dose_num),
    log10_Rate = log10(Rate)
  )

fig1_extra_both_log <- ggplot(dose_resp_pub_dat_log,
                              aes(x = log10_Dose, y = log10_Rate,
                                  colour = T_fac, fill = T_fac)) +
  geom_point(size = 1.7, alpha = 0.8) +
  geom_smooth(aes(group = T_fac), method = "gam",
              formula = y ~ s(x, k = 4, bs = "tp"),
              se = TRUE, alpha = 0.15, linewidth = 0.8) +
  scale_colour_manual(values = temp_pal_extra,
                      labels = paste0(temp_labels_extra, "\u00b0C")) +
  scale_fill_manual(values = temp_pal_extra,
                    labels = paste0(temp_labels_extra, "\u00b0C")) +
  facet_wrap(~Trait, scales = "free_y", ncol = 2) +
  labs(
    x      = "log\u2081\u2080 Prothioconazole concentration (mg L\u207b\u00b9)",
    y      = "log\u2081\u2080 Rate",
    colour = "Temperature",
    fill   = "Temperature"
  ) +
  theme_classic(11) +
  theme(strip.background = element_blank(),
        strip.text = element_text(face = "bold"))

ggsave(
  file.path(pub_fig_dir, "fig1_extra_dose_response_by_temp_both_log.png"),
  fig1_extra_both_log, width = 12, height = 5.5, dpi = 300, bg = "white"
)
message("Saved: ", file.path(pub_fig_dir, "fig1_extra_dose_response_by_temp_both_log.png"))

# CSV for both-log dose-response plot (log10-transformed values, control excluded)
fig1_both_log_csv <- file.path(
  tables_dir, "fig1_extra_dose_response_by_temp_both_log.csv"
)
dose_resp_pub_dat_log %>%
  transmute(
    Trait,
    Temperature_C              = T,
    Replicate                  = as.character(Replicate),
    Prothioconazole_mg_L       = Dose_num,
    log10_Prothioconazole_mg_L = log10_Dose,
    Rate_raw                   = Rate,
    log10_Rate                 = log10_Rate
  ) %>%
  write.csv(fig1_both_log_csv, row.names = FALSE)
message("Saved: ", fig1_both_log_csv)

# Both axes linear (no log scale)
fig1_extra_no_log <- ggplot(dose_resp_pub_dat,
                            aes(x = Dose_num, y = Rate,
                                colour = T_fac, fill = T_fac)) +
  geom_point(size = 1.7, alpha = 0.8) +
  geom_smooth(aes(group = T_fac), method = "gam",
              formula = y ~ s(x, k = 4, bs = "tp"),
              se = TRUE, alpha = 0.15, linewidth = 0.8) +
  scale_colour_manual(values = temp_pal_extra,
                      labels = paste0(temp_labels_extra, "\u00b0C")) +
  scale_fill_manual(values = temp_pal_extra,
                    labels = paste0(temp_labels_extra, "\u00b0C")) +
  facet_wrap(~Trait, scales = "free_y", ncol = 2) +
  labs(
    x      = "Prothioconazole concentration (mg L\u207b\u00b9)",
    y      = "Rate",
    colour = "Temperature",
    fill   = "Temperature"
  ) +
  theme_classic(11) +
  theme(strip.background = element_blank(),
        strip.text = element_text(face = "bold"))

ggsave(
  file.path(pub_fig_dir, "fig1_extra_dose_response_by_temp_no_log.png"),
  fig1_extra_no_log, width = 12, height = 5.5, dpi = 300, bg = "white"
)
message("Saved: ", file.path(pub_fig_dir, "fig1_extra_dose_response_by_temp_no_log.png"))

# CSV for no-log dose-response plot
fig1_no_log_csv <- file.path(
  tables_dir, "fig1_extra_dose_response_by_temp_no_log.csv"
)
dose_resp_pub_dat %>%
  transmute(
    Trait,
    Temperature_C        = T,
    Replicate            = as.character(Replicate),
    Prothioconazole_mg_L = Dose_num,
    Rate_raw             = Rate
  ) %>%
  write.csv(fig1_no_log_csv, row.names = FALSE)
message("Saved: ", fig1_no_log_csv)

# ---------------------------------------------------------------------------
# Fig: CUE vs prothioconazole dose
# Model: logit(CUE) ~ 0 + T_fac + T_fac:log10_dose  (per-temperature
#        intercepts and slopes; fitted on treated doses only)
# Left panel : CUE vs log10(dose) with fitted curves per temperature
# Right panel: posterior of per-temperature dose-response slope (beta_T)
# ---------------------------------------------------------------------------

# --- File paths ---
brms_cue_dose_rds <- file.path(models_dir,  "brms_logit_CUE_dose_slope_by_temp.rds")
post_cue_dose_csv <- file.path(tables_dir,  "posterior_logit_CUE_dose_slope_by_temp.csv")
fig_cue_dose_png  <- file.path(pub_fig_dir, "fig_cue_vs_dose.png")

# --- Data ---
# Treated doses only for the model (log10 of 0 is undefined for Control).
# Control points are kept for visual reference on the left panel only.
cue_dose_dat <- cue_dat %>%
  mutate(
    logit_CUE  = safe_logit(y_raw),
    T_fac      = factor(T),
    Dose_num   = suppressWarnings(as.numeric(as.character(Dose))),
    # Control gets pseudo log10-dose of log10(0.03) = half the lowest dose
    log10_dose = dplyr::if_else(
      is.na(Dose_num),
      log10(0.03),
      log10(Dose_num)
    ),
    is_control = as.character(Dose) == "Control"
  )

cue_dose_model_dat <- cue_dose_dat        # Control now included

# --- Temperature colour palette ---
temp_vals_cue <- sort(unique(cue_dose_dat$T))
temp_pal_cue  <- grDevices::hcl.colors(length(temp_vals_cue), palette = "Temps")
names(temp_pal_cue) <- as.character(temp_vals_cue)

# --- Fit or load model ---
if (file.exists(brms_cue_dose_rds)) {
  fit_cue_dose <- readRDS(brms_cue_dose_rds)
} else {
  fit_cue_dose <- brms::brm(
    logit_CUE ~ 0 + T_fac + T_fac:log10_dose,
    data    = cue_dose_model_dat,
    family  = gaussian(),
    prior   = c(
      brms::prior("normal(0, 3)",    class = "b"),
      brms::prior("exponential(1)", class = "sigma")
    ),
    iter    = BAYES_ITER,
    warmup  = BAYES_WARMUP,
    chains  = BAYES_CHAINS,
    seed    = BAYES_SEED,
    control = list(adapt_delta = BAYES_ADAPT, max_treedepth = BAYES_MAX_TD),
    backend = "rstan",
    refresh = 0,
    init    = 0
  )
  saveRDS(fit_cue_dose, brms_cue_dose_rds)
}

# --- Posterior draws: extract per-temperature slopes ---
# Parameters named "b_T_facXX:log10_dose" in brms output
raw_draws <- brms::as_draws_df(fit_cue_dose)

slope_draws <- raw_draws %>%
  dplyr::select(dplyr::matches("^b_T_fac.*:log10_dose$")) %>%
  tidyr::pivot_longer(
    dplyr::everything(),
    names_to  = "param",
    values_to = "slope"
  ) %>%
  dplyr::mutate(
    T = as.numeric(sub("^b_T_fac(.+):log10_dose$", "\\1", param)),
    T_fac = factor(T)
  )

readr::write_csv(slope_draws, post_cue_dose_csv)
message("Saved: ", post_cue_dose_csv)

# --- Fitted curves: posterior median + 95% CI on a log10-dose grid ---
dose_grid <- expand.grid(
  T_fac      = factor(temp_vals_cue),
  log10_dose = seq(log10(0.03), log10(4), length.out = 60)
) %>%
  dplyr::mutate(T = as.numeric(as.character(T_fac)))

# Use brms posterior_epred for credible ribbon
epred_mat <- brms::posterior_epred(fit_cue_dose, newdata = dose_grid, re_formula = NA)

dose_grid$median_logit <- apply(epred_mat, 2, median)
dose_grid$q2.5_logit   <- apply(epred_mat, 2, quantile, 0.025)
dose_grid$q97.5_logit  <- apply(epred_mat, 2, quantile, 0.975)

# Back-transform to CUE [0, 1]
dose_grid <- dose_grid %>%
  dplyr::mutate(
    median_CUE = plogis(median_logit),
    q2.5_CUE   = plogis(q2.5_logit),
    q97.5_CUE  = plogis(q97.5_logit),
    Dose_num   = 10^log10_dose
  )

# --- Left panel: CUE vs dose (log10 x), fitted curves + raw data ---
p_cue_dose_left <- ggplot(
  dose_grid,
  aes(x = Dose_num, colour = T_fac, fill = T_fac)
) +
  # 95% credible ribbon
  geom_ribbon(
    aes(ymin = q2.5_CUE, ymax = q97.5_CUE),
    alpha = 0.15, colour = NA
  ) +
  # Posterior median line
  geom_line(aes(y = median_CUE), linewidth = 0.9) +
  # Raw data: treated doses
  geom_point(
    data        = cue_dose_model_dat,
    aes(x = Dose_num, y = y_raw, colour = T_fac),
    size        = 1.7,
    alpha       = 0.7,
    inherit.aes = FALSE
  ) +
  # Raw data: Control shown as triangles coloured by temperature
  geom_point(
    data        = cue_dose_dat %>% dplyr::filter(is_control),
    aes(x = 10^log10_dose, y = y_raw, colour = T_fac),
    shape       = 17,
    size        = 1.7,
    alpha       = 0.7,
    inherit.aes = FALSE
  ) +
  scale_x_log10(
    breaks = c(0.06, 0.12, 0.25, 0.5, 1, 2, 4),
    labels = c("0.06", "0.12", "0.25", "0.5", "1", "2", "4")
  ) +
  scale_y_continuous(limits = c(0, 1), expand = expansion(mult = 0.02)) +
  scale_colour_manual(
    values = temp_pal_cue,
    labels = paste0(temp_vals_cue, "\u00b0C"),
    name   = "Temperature"
  ) +
  scale_fill_manual(
    values = temp_pal_cue,
    labels = paste0(temp_vals_cue, "\u00b0C"),
    name   = "Temperature"
  ) +
  labs(
    x = "Prothioconazole (mg L\u207b\u00b9)",
    y = "CUE"
  ) +
  annotate("text", x = 0.035, y = 0.02, label = "Control",
           colour = "grey50", size = 3, hjust = 0.5) +
  theme_classic(11) +
  theme(legend.position = "right")

# --- Right panel: per-temperature slope posterior half-eye ---
# A slope < 0 means CUE decreases with increasing dose (expected inhibition).
p_cue_dose_right <- ggplot(
  slope_draws %>% mutate(T_fac = factor(T, levels = rev(temp_vals_cue))),
  aes(x = slope, y = T_fac, fill = T_fac, colour = T_fac)
) +
  ggdist::stat_halfeye(
    adjust         = 0.8,
    width          = 0.65,
    .width         = c(0.66, 0.95),
    point_interval = ggdist::median_qi,
    slab_alpha     = 0.65,
    interval_size  = 0.7,
    point_size     = 1.5,
    side           = "right",
    justification  = -0.15
  ) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey40", linewidth = 0.4) +
  scale_fill_manual(
    values = temp_pal_cue,
    labels = paste0(temp_vals_cue, "\u00b0C")
  ) +
  scale_colour_manual(
    values = temp_pal_cue,
    labels = paste0(temp_vals_cue, "\u00b0C")
  ) +
  labs(
    x = expression(beta[T] ~ "(slope: \u0394logit(CUE) per log"[10]*" dose)"),
    y = "Temperature (\u00b0C)"
  ) +
  theme_classic(11) +
  theme(legend.position = "none", axis.text.y = element_text(size = 9))

# --- Combine and save ---
fig_cue_dose <- (p_cue_dose_left + theme(legend.position = "right")) +
  p_cue_dose_right +
  plot_layout(widths = c(2, 1)) +
  plot_annotation(
    caption = paste0(
      "logit(CUE) ~ 0 + T_fac + T_fac:log10_dose (Bayesian Gaussian, brms); Control assigned pseudo-dose = 0.03 mg/L.",
      "\nRibbon: 95% posterior credible interval; line: posterior median.",
      "\nTriangles: Control observations (included in model). Dashed line in right panel: slope = 0."
    )
  )

ggsave(
  fig_cue_dose_png,
  fig_cue_dose, width = 13, height = 5.5, dpi = 300, bg = "white"
)
message("Saved: ", fig_cue_dose_png)

