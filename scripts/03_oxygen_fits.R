# ===== Stage 1: Oxygen model fits, carbon unit conversions, and descriptive plots
#
# This script:
#  - Fits exponential-decay respiration models to dissolved oxygen time series data
#  - Computes initial cell counts (N0) and carbon-based metrics from oxygen consumption
#  - Derives growth rate, respiration rate, and CUE (Carbon Use Efficiency)
#  - Filters outlier points and generates descriptive plots
#  - Prepares data for downstream Bayesian thermal models
#
# Output files are saved to tables_dir for use by downstream scripts

# Source shared config (works interactively in RStudio, via source(), or from CLI)
.this_dir <- if (requireNamespace("rstudioapi", quietly = TRUE) &&
                 rstudioapi::isAvailable() &&
                 nzchar(rstudioapi::getActiveDocumentContext()$path)) {
  dirname(rstudioapi::getActiveDocumentContext()$path)
} else {
  tryCatch(dirname(sys.frame(1)$ofile), error = function(e) getwd())
}
source(file.path(.this_dir, "00_config.R"))

suppressPackageStartupMessages({
  library(minpack.lm)
  library(brms)
  library(posterior)
  library(grid)
})

options(mc.cores = max(1, parallel::detectCores() - 1))

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

# Save intermediates for downstream scripts
saveRDS(results_plot, file.path(models_dir, "results_plot.rds"))
saveRDS(ratio_dat,    file.path(models_dir, "ratio_dat.rds"))
saveRDS(dose_levels,  file.path(models_dir, "dose_levels.rds"))
saveRDS(dose_cols,    file.path(models_dir, "dose_cols.rds"))
saveRDS(dose_key_tbl, file.path(models_dir, "dose_key_tbl.rds"))

# Save data prepared for Bayesian models
saveRDS(growth_dat,          file.path(models_dir, "growth_dat.rds"))
saveRDS(resp_dat,            file.path(models_dir, "resp_dat.rds"))
saveRDS(ratio_dat_arr,       file.path(models_dir, "ratio_dat_arr.rds"))
saveRDS(cue_dat,             file.path(models_dir, "cue_dat.rds"))
saveRDS(growth_biomass_dat,  file.path(models_dir, "growth_biomass_dat.rds"))
saveRDS(resp_biomass_dat,    file.path(models_dir, "resp_biomass_dat.rds"))

message("03_oxygen_fits.R complete.")
