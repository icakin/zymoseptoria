# =============================================================================
# Fungal Growth Rate Dose-Response Modelling Across Temperatures
# =============================================================================
# Fits a Hill inhibition model to fungal growth rate vs prothioconazole
# concentration at 7 temperatures (15–28°C) using Bayesian hierarchical
# nonlinear regression in brms. All individual observations used (not means).
# Zero-concentration controls retained throughout.
#
# Reports EC50 and pEC50 following Srinivasan & Lloyd (2024, J. Med. Chem.).
# =============================================================================

library(tidyverse)
library(brms)
library(tidybayes)
library(patchwork)
library(ggdist)

# --- Paths -------------------------------------------------------------------

base_dir <- file.path(
  "/Users/g.yvon-durocher/Library/CloudStorage",
  "OneDrive-UniversityofExeter/Documents/work/Zymo",
  "temperature_prothioconazole/dose_response"
)

data_path <- file.path(base_dir, "data", "fig1_extra_dose_response_by_temp.csv")
# Growth rate outputs
fig_dir   <- file.path(base_dir, "figures", "growth")
mod_dir   <- file.path(base_dir, "models", "growth")
tab_dir   <- file.path(base_dir, "tables", "growth")

# Respiration rate outputs
fig_dir_resp <- file.path(base_dir, "figures", "respiration")
tab_dir_resp <- file.path(base_dir, "tables", "respiration")
mod_dir_resp <- file.path(base_dir, "models", "respiration")

for (d in c(fig_dir, mod_dir, tab_dir, fig_dir_resp, tab_dir_resp, mod_dir_resp)) {
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
}

# --- Colour palette (colour-blind friendly, viridis-inspired) ----------------

temp_pal <- c(
  "15" = "#440154", "18" = "#443A83", "21" = "#31688E",
  "24" = "#21908C", "26" = "#35B779", "27" = "#8FD744", "28" = "#FDE725"
)

# --- Prothioconazole molecular weight (g/mol) --------------------------------
# Used for pEC50 conversion: pEC50 = -log10(EC50_molar)
# where EC50_molar = EC50_mg_L / (MW * 1000)
MW_prothioconazole <- 344.26

# =============================================================================
# 1. DATA WRANGLING
# =============================================================================

dat_all <- read_csv(data_path, show_col_types = FALSE) %>%
  mutate(
    conc = Prothioconazole_mg_L,
    rate = Rate_raw,
    temperature = factor(Temperature_C),
    rep_id = factor(paste(Temperature_C, Replicate, sep = "_"))
  )

# Split by Trait: growth rate for this analysis, respiration kept separately
dat <- dat_all %>% filter(grepl("Growth", Trait))
dat_resp <- dat_all %>% filter(grepl("Respiration", Trait))

cat("\n--- Rows by Trait ---\n")
cat("Growth rate:", nrow(dat), "\n")
cat("Respiration rate:", nrow(dat_resp), "\n")

# For log-scale plotting: zero-conc controls at pseudo-value (half min non-zero)
min_nonzero_conc <- min(dat$conc[dat$conc > 0])
dat <- dat %>%
  mutate(conc_plot = ifelse(conc == 0, min_nonzero_conc / 2, conc))

# Summary statistics
dat_summ <- dat %>%
  group_by(temperature, conc) %>%
  summarise(
    n         = n(),
    mean_rate = mean(rate, na.rm = TRUE),
    se_rate   = sd(rate, na.rm = TRUE) / sqrt(n()),
    .groups   = "drop"
  ) %>%
  mutate(conc_plot = ifelse(conc == 0, min_nonzero_conc / 2, conc))

cat("\n--- Data summary (n per temperature x concentration) ---\n")
print(dat_summ %>% select(temperature, conc, n) %>% pivot_wider(
  names_from = conc, values_from = n
))

write_csv(dat_summ, file.path(tab_dir, "01_data_summary.csv"))

# =============================================================================
# 2. EXPLORATORY PLOTS
# =============================================================================

# --- 01: Rate vs concentration (linear axes, all data) -----------------------

p01 <- ggplot(dat, aes(x = conc, y = rate, colour = temperature)) +
  geom_point(alpha = 0.5, size = 2) +
  geom_point(data = dat_summ, aes(y = mean_rate), size = 3, shape = 18) +
  geom_line(data = dat_summ, aes(y = mean_rate), linewidth = 0.5, alpha = 0.5) +
  scale_colour_manual(values = temp_pal, name = "Temperature (°C)") +
  labs(
    x = "Prothioconazole (mg/L)",
    y = expression("Growth rate (C C"^{-1}~"h"^{-1}*")"),
    title = "Dose-response: growth rate vs prothioconazole concentration",
    subtitle = "All replicates shown; diamonds = means"
  ) +
  theme_bw(base_size = 13)
ggsave(file.path(fig_dir, "01_rate_vs_conc_linear.png"),
       p01, width = 8, height = 5, dpi = 300)

# --- 02: Log-log plot (all data including zero-conc at pseudo-value) ---------

p02 <- ggplot(dat, aes(x = conc_plot, y = rate, colour = temperature)) +
  geom_point(alpha = 0.6, size = 2) +
  scale_x_log10() +
  scale_y_log10() +
  scale_colour_manual(values = temp_pal, name = "Temperature (°C)") +
  labs(
    x = "Prothioconazole concentration (mg/L, log scale)",
    y = expression("Growth rate (C C"^{-1}~"h"^{-1}*", log scale)"),
    title = "Dose-response in log-log space",
    subtitle = "Zero-concentration controls plotted at 0.03 mg/L"
  ) +
  theme_bw(base_size = 13)
ggsave(file.path(fig_dir, "02_rate_vs_conc_loglog.png"),
       p02, width = 8, height = 5, dpi = 300)

# --- 03: Faceted by temperature (linear axes, all data) ----------------------

p03 <- ggplot(dat, aes(x = conc, y = rate)) +
  geom_point(aes(colour = temperature), alpha = 0.5, size = 1.5) +
  geom_point(data = dat_summ, aes(y = mean_rate), size = 2.5, shape = 18) +
  geom_errorbar(data = dat_summ,
                aes(y = mean_rate,
                    ymin = mean_rate - se_rate,
                    ymax = mean_rate + se_rate),
                width = 0.1) +
  facet_wrap(~ temperature, ncol = 4, labeller = label_both) +
  scale_colour_manual(values = temp_pal, guide = "none") +
  labs(
    x = "Prothioconazole (mg/L)",
    y = expression("Growth rate (C C"^{-1}~"h"^{-1}*")"),
    title = "Dose-response faceted by temperature",
    subtitle = "Points = replicates, diamonds = means +/- SE"
  ) +
  theme_bw(base_size = 12)
ggsave(file.path(fig_dir, "03_rate_vs_conc_faceted.png"),
       p03, width = 12, height = 6, dpi = 300)

# =============================================================================
# 3. FIT HILL INHIBITION MODEL
# =============================================================================
# rate = r0 / (1 + (conc / ec50)^n)
# r0 varies by temperature (fixed) + replicate (random)
# ec50 varies by temperature (fixed)
# n (Hill coefficient) varies by temperature (fixed)

formula_hill <- bf(
  rate ~ r0 / (1 + (conc / ec50)^n),
  r0   ~ temperature + (1 | rep_id),
  ec50 ~ temperature,
  n    ~ temperature,
  nl = TRUE
)

priors_hill <- c(
  prior(normal(0.04, 0.02), nlpar = "r0", lb = 0),
  prior(student_t(3, 0, 0.01), class = "sd", nlpar = "r0"),
  prior(lognormal(0, 1), nlpar = "ec50", lb = 0),
  prior(normal(1.5, 1), nlpar = "n", lb = 0.1)
)

cat("\n--- Fitting Hill model ---\n")

fit_hill <- brm(
  formula  = formula_hill,
  data     = dat,
  prior    = priors_hill,
  family   = gaussian(),
  chains   = 4,
  iter     = 4000,
  warmup   = 2000,
  cores    = 4,
  seed     = 42,
  control  = list(adapt_delta = 0.95, max_treedepth = 12),
  file     = file.path(mod_dir, "hill_growth")
)

cat("\n--- Hill Model Summary ---\n")
print(summary(fit_hill))

hill_fixed <- as_tibble(fixef(fit_hill), rownames = "parameter")
write_csv(hill_fixed, file.path(tab_dir, "02_hill_model_parameters.csv"))

# =============================================================================
# 4. EXTRACT HILL PARAMETERS BY TEMPERATURE
# =============================================================================

# Helper to reconstruct temperature-specific values from intercept + contrasts
reconstruct_param <- function(model, param) {
  int_name <- paste0("b_", param, "_Intercept")
  contrast_names <- paste0("b_", param, "_temperature",
                           c("18", "21", "24", "26", "27", "28"))
  all_names <- c(int_name, contrast_names)

  draws <- model %>% spread_draws(!!!syms(all_names))

  draws %>%
    mutate(
      val_15 = !!sym(int_name),
      val_18 = !!sym(int_name) + !!sym(contrast_names[1]),
      val_21 = !!sym(int_name) + !!sym(contrast_names[2]),
      val_24 = !!sym(int_name) + !!sym(contrast_names[3]),
      val_26 = !!sym(int_name) + !!sym(contrast_names[4]),
      val_27 = !!sym(int_name) + !!sym(contrast_names[5]),
      val_28 = !!sym(int_name) + !!sym(contrast_names[6])
    ) %>%
    select(.draw, starts_with("val_")) %>%
    pivot_longer(starts_with("val_"), names_to = "temperature",
                 values_to = "value", names_prefix = "val_") %>%
    mutate(temperature = factor(temperature,
                                levels = as.character(sort(as.numeric(
                                  unique(temperature))))))
}

ec50_draws <- reconstruct_param(fit_hill, "ec50")
n_draws    <- reconstruct_param(fit_hill, "n")
r0_draws   <- reconstruct_param(fit_hill, "r0")

ec50_summary <- ec50_draws %>%
  group_by(temperature) %>%
  median_qi(value, .width = 0.95) %>%
  rename(EC50 = value, EC50_lower = .lower, EC50_upper = .upper)

n_summary <- n_draws %>%
  group_by(temperature) %>%
  median_qi(value, .width = 0.95) %>%
  rename(n = value, n_lower = .lower, n_upper = .upper)

r0_summary <- r0_draws %>%
  group_by(temperature) %>%
  median_qi(value, .width = 0.95) %>%
  rename(r0 = value, r0_lower = .lower, r0_upper = .upper)

hill_param_tbl <- ec50_summary %>%
  select(temperature, EC50, EC50_lower, EC50_upper) %>%
  left_join(n_summary %>% select(temperature, n, n_lower, n_upper),
            by = "temperature") %>%
  left_join(r0_summary %>% select(temperature, r0, r0_lower, r0_upper),
            by = "temperature")

cat("\n--- Hill model parameters by temperature ---\n")
print(hill_param_tbl)
write_csv(hill_param_tbl, file.path(tab_dir, "03_hill_parameters_by_temperature.csv"))

# =============================================================================
# 5. PARAMETER PLOTS
# =============================================================================

# --- 04: EC50 by temperature -------------------------------------------------

p04 <- ggplot(ec50_summary, aes(x = temperature, y = EC50,
                                 colour = temperature)) +
  geom_pointinterval(aes(ymin = EC50_lower, ymax = EC50_upper), size = 3) +
  scale_colour_manual(values = temp_pal, guide = "none") +
  labs(
    x = "Temperature (°C)",
    y = expression(EC[50]~"(mg/L)"),
    title = expression("EC"[50]~"by temperature — Hill model"),
    subtitle = "Posterior median and 95% credible interval"
  ) +
  theme_bw(base_size = 13)
ggsave(file.path(fig_dir, "04_ec50_by_temperature.png"),
       p04, width = 7, height = 5, dpi = 300)

# --- 05: EC50 posterior densities --------------------------------------------

p05 <- ggplot(ec50_draws, aes(x = value, fill = temperature,
                               colour = temperature)) +
  geom_density(alpha = 0.35, linewidth = 0.6) +
  scale_fill_manual(values = temp_pal, name = "Temperature (°C)") +
  scale_colour_manual(values = temp_pal, name = "Temperature (°C)") +
  labs(
    x = expression(EC[50]~"(mg/L)"),
    y = "Posterior density",
    title = expression("Posterior distributions of EC"[50]~"across temperatures")
  ) +
  theme_bw(base_size = 13)
ggsave(file.path(fig_dir, "05_ec50_posterior_density.png"),
       p05, width = 8, height = 5, dpi = 300)

# --- 06: Hill n by temperature -----------------------------------------------

p06 <- ggplot(n_summary, aes(x = temperature, y = n, colour = temperature)) +
  geom_pointinterval(aes(ymin = n_lower, ymax = n_upper), size = 3) +
  scale_colour_manual(values = temp_pal, guide = "none") +
  labs(
    x = "Temperature (°C)",
    y = "Hill coefficient (n)",
    title = "Hill coefficient (n) by temperature",
    subtitle = "Posterior median and 95% credible interval"
  ) +
  theme_bw(base_size = 13)
ggsave(file.path(fig_dir, "06_hill_n_by_temperature.png"),
       p06, width = 7, height = 5, dpi = 300)

# --- 07: Hill n posterior densities ------------------------------------------

p07 <- ggplot(n_draws, aes(x = value, fill = temperature,
                            colour = temperature)) +
  geom_density(alpha = 0.35, linewidth = 0.6) +
  scale_fill_manual(values = temp_pal, name = "Temperature (°C)") +
  scale_colour_manual(values = temp_pal, name = "Temperature (°C)") +
  labs(
    x = "Hill coefficient (n)",
    y = "Posterior density",
    title = "Posterior distributions of Hill coefficient across temperatures"
  ) +
  theme_bw(base_size = 13)
ggsave(file.path(fig_dir, "07_hill_n_posterior_density.png"),
       p07, width = 8, height = 5, dpi = 300)

# =============================================================================
# 6. MODEL PREDICTIONS
# =============================================================================

pred_grid <- expand_grid(
  temperature = factor(sort(unique(as.numeric(levels(dat$temperature)))),
                       levels = levels(dat$temperature)),
  conc = c(0, seq(min_nonzero_conc, 4, length.out = 200))
) %>%
  mutate(conc_plot = ifelse(conc == 0, min_nonzero_conc / 2, conc))

# Population-level predictions (marginalise over replicate random effects)
hill_pred_draws <- pred_grid %>%
  add_epred_draws(fit_hill, ndraws = 500, re_formula = NA)
hill_pred_summary <- hill_pred_draws %>%
  group_by(temperature, conc, conc_plot) %>%
  median_qi(.epred, .width = 0.95)
hill_pred_nz <- hill_pred_summary %>% filter(conc > 0)

# --- 08: Overlay, linear axes ------------------------------------------------

p08 <- ggplot() +
  geom_ribbon(data = hill_pred_summary,
              aes(x = conc, ymin = .lower, ymax = .upper,
                  fill = temperature), alpha = 0.15) +
  geom_line(data = hill_pred_summary,
            aes(x = conc, y = .epred, colour = temperature),
            linewidth = 0.8) +
  geom_point(data = dat,
             aes(x = conc, y = rate, colour = temperature),
             alpha = 0.5, size = 1.8) +
  scale_colour_manual(values = temp_pal, name = "Temperature (°C)") +
  scale_fill_manual(values = temp_pal, name = "Temperature (°C)") +
  labs(
    x = "Prothioconazole (mg/L)",
    y = expression("Growth rate (C C"^{-1}~"h"^{-1}*")"),
    title = expression("Growth rate (C C"^{-1}~"h"^{-1}*")"),
    subtitle = "Hill model — posterior median and 95% credible band"
  ) +
  theme_bw(base_size = 13)
ggsave(file.path(fig_dir, "08_hill_fit_linear.png"),
       p08, width = 8, height = 5, dpi = 300)

# --- 09: Faceted, linear axes ------------------------------------------------

p09 <- ggplot() +
  geom_ribbon(data = hill_pred_summary,
              aes(x = conc, ymin = .lower, ymax = .upper),
              fill = "steelblue", alpha = 0.2) +
  geom_line(data = hill_pred_summary,
            aes(x = conc, y = .epred),
            colour = "steelblue", linewidth = 0.8) +
  geom_point(data = dat, aes(x = conc, y = rate), alpha = 0.6, size = 1.5) +
  facet_wrap(~ temperature, ncol = 4, labeller = label_both) +
  labs(
    x = "Prothioconazole (mg/L)",
    y = expression("Growth rate (C C"^{-1}~"h"^{-1}*")"),
    title = "Hill model fit by temperature (linear scale)",
    subtitle = "Posterior median and 95% credible band"
  ) +
  theme_bw(base_size = 12)
ggsave(file.path(fig_dir, "09_hill_fit_linear_faceted.png"),
       p09, width = 12, height = 6, dpi = 300)

# --- 10: Overlay, log-log axes -----------------------------------------------

p10 <- ggplot() +
  geom_ribbon(data = hill_pred_nz,
              aes(x = conc, ymin = .lower, ymax = .upper,
                  fill = temperature), alpha = 0.15) +
  geom_line(data = hill_pred_nz,
            aes(x = conc, y = .epred, colour = temperature),
            linewidth = 0.8) +
  geom_point(data = dat,
             aes(x = conc_plot, y = rate, colour = temperature),
             alpha = 0.6, size = 2) +
  scale_x_log10() + scale_y_log10() +
  scale_colour_manual(values = temp_pal, name = "Temperature (°C)") +
  scale_fill_manual(values = temp_pal, name = "Temperature (°C)") +
  labs(
    x = "Prothioconazole concentration (mg/L, log scale)",
    y = expression("Rate (log scale)"),
    title = expression("Growth rate (C C"^{-1}~"h"^{-1}*")"),
    subtitle = "Hill model — log-log scale"
  ) +
  theme_bw(base_size = 13)
ggsave(file.path(fig_dir, "10_hill_fit_loglog.png"),
       p10, width = 8, height = 5, dpi = 300)

# --- 11: Faceted, log-log axes -----------------------------------------------

p11 <- ggplot() +
  geom_ribbon(data = hill_pred_nz,
              aes(x = conc, ymin = .lower, ymax = .upper),
              fill = "steelblue", alpha = 0.2) +
  geom_line(data = hill_pred_nz,
            aes(x = conc, y = .epred),
            colour = "steelblue", linewidth = 0.8) +
  geom_point(data = dat, aes(x = conc_plot, y = rate),
             alpha = 0.6, size = 1.5) +
  scale_x_log10() + scale_y_log10() +
  facet_wrap(~ temperature, ncol = 4, labeller = label_both) +
  labs(
    x = "Prothioconazole concentration (mg/L, log scale)",
    y = expression("Rate (log scale)"),
    title = "Hill model fit by temperature (log-log scale)",
    subtitle = "Posterior median and 95% credible band"
  ) +
  theme_bw(base_size = 12)
ggsave(file.path(fig_dir, "11_hill_fit_loglog_faceted.png"),
       p11, width = 12, height = 6, dpi = 300)

# =============================================================================
# 7. POSTERIOR PREDICTIVE CHECK
# =============================================================================

p12 <- pp_check(fit_hill, ndraws = 100) +
  labs(title = "Posterior predictive check — Hill model") +
  theme_bw(base_size = 13)
ggsave(file.path(fig_dir, "12_pp_check.png"),
       p12, width = 7, height = 5, dpi = 300)

# =============================================================================
# 8. EC50 SUMMARY PLOT
# =============================================================================

p13 <- ggplot(ec50_summary, aes(y = temperature, x = EC50,
                                 colour = temperature)) +
  geom_pointinterval(aes(xmin = EC50_lower, xmax = EC50_upper),
                     size = 4, linewidth = 1.2) +
  scale_colour_manual(values = temp_pal, guide = "none") +
  labs(
    y = "Temperature (°C)",
    x = expression(EC[50]~"(mg/L prothioconazole)"),
    title = expression(bold("Antifungal sensitivity: EC"[50]~"across temperatures")),
    subtitle = "Hill model — posterior median and 95% credible interval"
  ) +
  theme_bw(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold")
  )
ggsave(file.path(fig_dir, "13_ec50_summary.png"),
       p13, width = 8, height = 5, dpi = 300)

# =============================================================================
# 9. pEC50 CONVERSION (Srinivasan & Lloyd, 2024)
# =============================================================================
# pEC50 = -log10(EC50_molar) = -log10(EC50_mg_L / (MW * 1000))
# Higher pEC50 = greater potency/sensitivity (analogous to pH).

pec50_draws <- ec50_draws %>%
  mutate(
    ec50_mg_L  = value,
    ec50_molar = ec50_mg_L / (MW_prothioconazole * 1000),
    pec50      = -log10(ec50_molar)
  )

pec50_summary <- pec50_draws %>%
  group_by(temperature) %>%
  median_qi(pec50, .width = 0.95) %>%
  rename(pEC50 = pec50, pEC50_lower = .lower, pEC50_upper = .upper)

pec50_report <- pec50_draws %>%
  group_by(temperature) %>%
  summarise(
    pEC50_mean   = mean(pec50),
    pEC50_sd     = sd(pec50),
    pEC50_median = median(pec50),
    pEC50_lower  = quantile(pec50, 0.025),
    pEC50_upper  = quantile(pec50, 0.975),
    EC50_mg_L_median = median(ec50_mg_L),
    EC50_mg_L_lower  = quantile(ec50_mg_L, 0.025),
    EC50_mg_L_upper  = quantile(ec50_mg_L, 0.975),
    .groups = "drop"
  )

cat("\n--- pEC50 by temperature (Srinivasan & Lloyd 2024 convention) ---\n")
print(pec50_report)
write_csv(pec50_report, file.path(tab_dir, "04_pec50_by_temperature.csv"))

# --- 14: pEC50 by temperature — primary result figure -------------------------

p14 <- ggplot(pec50_summary, aes(y = temperature, x = pEC50,
                                  colour = temperature)) +
  geom_pointinterval(aes(xmin = pEC50_lower, xmax = pEC50_upper),
                     size = 4, linewidth = 1.2) +
  scale_colour_manual(values = temp_pal, guide = "none") +
  annotate("text", x = Inf, y = Inf, label = "More sensitive \u2192",
           hjust = 1.1, vjust = 2, size = 3.5, colour = "grey40",
           fontface = "italic") +
  labs(
    y = "Temperature (°C)",
    x = expression(pEC[50]~"(-log"[10]~"M)"),
    title = expression(bold("Antifungal sensitivity: pEC"[50]~"across temperatures")),
    subtitle = "Hill model — posterior median and 95% CrI"
  ) +
  theme_bw(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold")
  )
ggsave(file.path(fig_dir, "14_pec50_by_temperature.png"),
       p14, width = 8, height = 5, dpi = 300)

# --- 15: pEC50 posterior densities --------------------------------------------

p15 <- ggplot(pec50_draws, aes(x = pec50, fill = temperature,
                                colour = temperature)) +
  geom_density(alpha = 0.35, linewidth = 0.6) +
  scale_fill_manual(values = temp_pal, name = "Temperature (°C)") +
  scale_colour_manual(values = temp_pal, name = "Temperature (°C)") +
  labs(
    x = expression(pEC[50]~"(-log"[10]~"M)"),
    y = "Posterior density",
    title = expression("Posterior distributions of pEC"[50]~"across temperatures")
  ) +
  theme_bw(base_size = 13)
ggsave(file.path(fig_dir, "15_pec50_posterior_density.png"),
       p15, width = 8, height = 5, dpi = 300)

# --- 16: pEC50 forest plot (mean +/- SD) --------------------------------------

p16 <- ggplot(pec50_report, aes(y = temperature, x = pEC50_mean,
                                 colour = temperature)) +
  geom_point(size = 3) +
  geom_errorbar(aes(xmin = pEC50_mean - pEC50_sd,
                    xmax = pEC50_mean + pEC50_sd),
                width = 0.2, linewidth = 0.8, orientation = "y") +
  geom_errorbar(aes(xmin = pEC50_mean - 2 * pEC50_sd,
                    xmax = pEC50_mean + 2 * pEC50_sd),
                width = 0.1, linewidth = 0.4, linetype = "dashed",
                orientation = "y") +
  scale_colour_manual(values = temp_pal, guide = "none") +
  labs(
    y = "Temperature (°C)",
    x = expression(pEC[50]~"(-log"[10]~"M)"),
    title = expression("pEC"[50]~"\u00b1 SD across temperatures"),
    subtitle = "Solid bars = \u00b11 SD; dashed bars = \u00b12 SD"
  ) +
  theme_bw(base_size = 14) +
  theme(panel.grid.minor = element_blank())
ggsave(file.path(fig_dir, "16_pec50_forest_plot.png"),
       p16, width = 8, height = 5, dpi = 300)

# =============================================================================
# 10. PAIRWISE TEMPERATURE COMPARISONS ON pEC50 SCALE
# =============================================================================
# Bayesian pairwise contrasts on the log10 scale as recommended by
# Srinivasan & Lloyd (2024).

temp_levels <- levels(pec50_draws$temperature)
pairs <- combn(temp_levels, 2, simplify = FALSE)

pairwise_results <- map_dfr(pairs, function(pair) {
  draws_a <- pec50_draws %>% filter(temperature == pair[1]) %>% pull(pec50)
  draws_b <- pec50_draws %>% filter(temperature == pair[2]) %>% pull(pec50)
  diff <- draws_a - draws_b

  tibble(
    temp_A      = pair[1],
    temp_B      = pair[2],
    mean_diff   = mean(diff),
    sd_diff     = sd(diff),
    lower_95    = quantile(diff, 0.025),
    upper_95    = quantile(diff, 0.975),
    prob_A_gt_B = mean(diff > 0),
    prob_B_gt_A = mean(diff < 0)
  )
})

cat("\n--- Pairwise pEC50 comparisons ---\n")
cat("prob_A_gt_B = posterior probability that temperature A has higher pEC50\n")
cat("(greater sensitivity) than temperature B.\n")
cat("Values > 0.95 or < 0.05 indicate strong evidence.\n\n")
print(pairwise_results, n = Inf)
write_csv(pairwise_results, file.path(tab_dir, "05_pec50_pairwise_comparisons.csv"))

# --- 17: Pairwise comparison heatmap -----------------------------------------

# Build full matrix with diagonal
all_temps_f17 <- levels(pec50_draws$temperature)
pair_matrix <- bind_rows(
  pairwise_results %>% select(temp_A, temp_B, prob_A_gt_B),
  pairwise_results %>% transmute(temp_A = temp_B, temp_B = temp_A,
                                  prob_A_gt_B = prob_B_gt_A),
  tibble(temp_A = all_temps_f17, temp_B = all_temps_f17,
         prob_A_gt_B = NA_real_)
) %>%
  distinct(temp_A, temp_B, .keep_all = TRUE) %>%
  mutate(
    temp_A = factor(temp_A, levels = all_temps_f17),
    temp_B = factor(temp_B, levels = all_temps_f17),
    label = ifelse(is.na(prob_A_gt_B), "\u2014", sprintf("%.2f", prob_A_gt_B))
  )

p17 <- ggplot(pair_matrix, aes(x = temp_B, y = temp_A, fill = prob_A_gt_B)) +
  geom_tile(colour = "white", linewidth = 0.5) +
  geom_text(aes(label = label), size = 3.5) +
  scale_fill_gradient2(
    low = "#2166AC", mid = "white", high = "#B2182B",
    midpoint = 0.5, limits = c(0, 1),
    na.value = "grey90",
    name = expression("P(row > col)")
  ) +
  labs(
    x = "Temperature (°C)",
    y = "Temperature (°C)",
    title = expression("Pairwise comparison of pEC"[50]~"between temperatures"),
    subtitle = expression("Values show P(row has higher pEC"[50]~"than column)")
  ) +
  theme_minimal(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 0),
    panel.grid = element_blank()
  ) +
  coord_equal()
ggsave(file.path(fig_dir, "17_pec50_pairwise_heatmap.png"),
       p17, width = 8, height = 7, dpi = 300)

# =============================================================================
# 11. PUBLICATION-READY SUMMARY TABLE
# =============================================================================
# Report pEC50 +/- SD as recommended by Srinivasan & Lloyd (2024)

summary_table <- pec50_report %>%
  mutate(
    pEC50_report = sprintf("%.2f \u00b1 %.2f", pEC50_mean, pEC50_sd),
    EC50_mg_L_report = sprintf("%.2f [%.2f, %.2f]",
                                EC50_mg_L_median,
                                EC50_mg_L_lower,
                                EC50_mg_L_upper)
  ) %>%
  select(temperature, pEC50_report, EC50_mg_L_report,
         pEC50_mean, pEC50_sd, pEC50_median, pEC50_lower, pEC50_upper,
         EC50_mg_L_median, EC50_mg_L_lower, EC50_mg_L_upper)

cat("\n--- Publication-ready summary: report pEC50 +/- SD (Srinivasan & Lloyd, 2024) ---\n")
print(summary_table)
write_csv(summary_table, file.path(tab_dir, "06_summary_table_publication.csv"))

# =============================================================================
# 12. HILL COEFFICIENT n SUMMARY
# =============================================================================
# Hill coefficient n interpretation (Srinivasan & Lloyd 2024; Prinz 2010):
#   n ~ 1: simple hyperbolic inhibition
#   n > 1: cooperative binding or threshold effects
#   n < 1: negative cooperativity or heterogeneous target populations

n_report <- n_draws %>%
  group_by(temperature) %>%
  summarise(
    n_mean   = mean(value),
    n_sd     = sd(value),
    n_median = median(value),
    n_lower  = quantile(value, 0.025),
    n_upper  = quantile(value, 0.975),
    .groups  = "drop"
  )

cat("\n--- Hill coefficient n by temperature ---\n")
print(n_report)
write_csv(n_report, file.path(tab_dir, "07_hill_n_by_temperature.csv"))

# =============================================================================
# 13. COMPUTE GLOBAL MEAN pEC50 (AVERAGED ACROSS TEMPERATURES)
# =============================================================================

# Average pEC50 across all 7 temperatures for each posterior draw
global_pec50_draws <- pec50_draws %>%
  group_by(.draw) %>%
  summarise(pec50 = mean(pec50), .groups = "drop")

global_pec50_mean <- mean(global_pec50_draws$pec50)
global_pec50_sd   <- sd(global_pec50_draws$pec50)
global_pec50_qi   <- quantile(global_pec50_draws$pec50, c(0.025, 0.975))

cat("\n--- Global mean pEC50 (averaged across temperatures) ---\n")
cat(sprintf("pEC50 = %.2f +/- %.2f [%.2f, %.2f]\n",
            global_pec50_mean, global_pec50_sd,
            global_pec50_qi[1], global_pec50_qi[2]))

# =============================================================================
# 14. GROWTH COMPOSITE PANELS (saved for combined figure later)
# =============================================================================

# Panel: Hill model fit, linear axes
pa_growth <- ggplot() +
  geom_ribbon(data = hill_pred_summary,
              aes(x = conc, ymin = .lower, ymax = .upper,
                  fill = temperature), alpha = 0.15) +
  geom_line(data = hill_pred_summary,
            aes(x = conc, y = .epred, colour = temperature),
            linewidth = 0.8) +
  geom_point(data = dat,
             aes(x = conc, y = rate, colour = temperature),
             alpha = 0.5, size = 1.8) +
  scale_colour_manual(values = temp_pal, name = "Temperature (°C)") +
  scale_fill_manual(values = temp_pal, name = "Temperature (°C)") +
  labs(
    x = "Prothioconazole (mg/L)",
    y = expression("Growth rate (C C"^{-1}~"h"^{-1}*")")
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    legend.position = "none"
  )

# Panel: pEC50 by temperature + global mean as separate point
# Add global mean as an extra row
pec50_plot_data <- pec50_report %>%
  select(temperature, pEC50_mean, pEC50_sd) %>%
  bind_rows(tibble(
    temperature = "Global",
    pEC50_mean  = global_pec50_mean,
    pEC50_sd    = global_pec50_sd
  )) %>%
  mutate(temperature = factor(temperature,
                              levels = c("15", "18", "21", "24", "26", "27", "28", "Global")))

# Vertical separator position: between 28 and Global (7.5 on numeric scale)
pb_growth <- ggplot(pec50_plot_data, aes(x = temperature, y = pEC50_mean,
                                          colour = temperature)) +
  geom_vline(xintercept = 7.5, linetype = "dashed", colour = "grey70") +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = pEC50_mean - pEC50_sd,
                    ymax = pEC50_mean + pEC50_sd),
                width = 0.2, linewidth = 0.8) +
  geom_errorbar(aes(ymin = pEC50_mean - 2 * pEC50_sd,
                    ymax = pEC50_mean + 2 * pEC50_sd),
                width = 0.1, linewidth = 0.4, linetype = "dashed") +
  scale_colour_manual(
    values = c(temp_pal, "Global" = "black"),
    guide = "none"
  ) +
  annotate("text", x = Inf, y = Inf, label = "More sensitive \u2191",
           hjust = 1.1, vjust = 1.5, size = 3, colour = "grey40",
           fontface = "italic") +
  labs(
    x = "Temperature (°C)",
    y = expression(pEC[50]~"(-log"[10]~"M)")
  ) +
  theme_bw(base_size = 11) +
  theme(panel.grid = element_blank())

cat("\n--- Growth composite panels prepared ---\n")

# =============================================================================
# SAVE GROWTH MODEL OBJECT
# =============================================================================

saveRDS(fit_hill, file.path(mod_dir, "fit_hill.rds"))

cat("\n--- Growth rate analysis complete ---\n")
cat("Figures:", fig_dir, "\n")
cat("Tables:", tab_dir, "\n")
cat("Model:", mod_dir, "\n")

# =============================================================================
# 15. BLISS INDEPENDENCE ANALYSIS — TEMPERATURE × FUNGICIDE INTERACTION
# =============================================================================
# Quantify whether high temperature and antifungal stress act synergistically,
# additively, or antagonistically on fungal growth rate using Bliss independence.
#
# Bliss independence assumes temperature and fungicide stress act through
# independent pathways, so their fractional effects multiply.
# delta < 0 → synergy (fungicide more effective than expected)
# delta > 0 → antagonism (fungicide less effective than expected)

cat("\n\n========== BLISS INDEPENDENCE ANALYSIS ==========\n\n")

# --- 15a. Identify reference temperature (highest baseline growth rate) -------

r0_medians <- r0_draws %>%
  group_by(temperature) %>%
  summarise(r0_med = median(value), .groups = "drop")
T_ref <- r0_medians$temperature[which.max(r0_medians$r0_med)]
cat("Reference temperature (highest r0):", as.character(T_ref), "°C\n")

# --- 15b. Join parameter draws into a single wide frame ----------------------

param_draws <- r0_draws %>%
  rename(r0 = value) %>%
  left_join(ec50_draws %>% rename(ec50 = value),
            by = c(".draw", "temperature")) %>%
  left_join(n_draws %>% rename(n_hill = value),
            by = c(".draw", "temperature"))

# --- 15c. Compute pointwise Bliss deviation on posterior draws ----------------

# Non-zero concentrations from the data
conc_levels <- sort(unique(dat$conc[dat$conc > 0]))

bliss_draws <- param_draws %>%
  # Attach reference temperature parameters for each draw
  left_join(
    param_draws %>%
      filter(temperature == T_ref) %>%
      select(.draw, r0_ref = r0, ec50_ref = ec50, n_ref = n_hill),
    by = ".draw"
  ) %>%
  # Expand across concentrations
  crossing(conc = conc_levels) %>%
  mutate(
    # Fractional effect at this temperature (Hill dose-response)
    f_obs = 1 / (1 + (conc / ec50)^n_hill),
    # Fractional effect at reference temperature
    f_ref = 1 / (1 + (conc / ec50_ref)^n_ref),
    # Bliss deviation: negative = synergy, positive = antagonism
    delta  = f_obs - f_ref
  )

# --- 15d. Summarise Bliss deviation -------------------------------------------

bliss_summary <- bliss_draws %>%
  group_by(temperature, conc) %>%
  summarise(
    delta_mean   = mean(delta),
    delta_median = median(delta),
    delta_lower  = quantile(delta, 0.025),
    delta_upper  = quantile(delta, 0.975),
    prob_synergy     = mean(delta < 0),
    prob_antagonism  = mean(delta > 0),
    .groups = "drop"
  )

cat("\n--- Pointwise Bliss deviation summary ---\n")
print(bliss_summary)
write_csv(bliss_summary, file.path(tab_dir, "08_bliss_deviation_pointwise.csv"))

# --- 15e. Figure 19: Bliss deviation vs temperature, coloured by dose --------

p19 <- ggplot(bliss_summary, aes(x = as.numeric(as.character(temperature)),
                                  y = delta_mean,
                                  colour = factor(conc))) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_ribbon(aes(ymin = delta_lower, ymax = delta_upper,
                  fill = factor(conc)),
              alpha = 0.1, colour = NA) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  scale_colour_viridis_d(name = "Dose (mg/L)", option = "plasma") +
  scale_fill_viridis_d(name = "Dose (mg/L)", option = "plasma", guide = "none") +
  annotate("text", x = 15, y = -Inf, label = "Synergy \u2193",
           hjust = 0, vjust = -0.5, size = 3, colour = "grey40", fontface = "italic") +
  annotate("text", x = 15, y = Inf, label = "Antagonism \u2191",
           hjust = 0, vjust = 1.5, size = 3, colour = "grey40", fontface = "italic") +
  labs(
    x = "Temperature (\u00b0C)",
    y = expression(Delta~"(Bliss deviation)"),
    title = "Bliss Independence: temperature \u00d7 fungicide interaction",
    subtitle = expression(Delta~"= observed fractional effect \u2212 expected under independence. Dashed = additivity.")
  ) +
  theme_bw(base_size = 13)
ggsave(file.path(fig_dir, "19_bliss_deviation_by_temperature.png"),
       p19, width = 8, height = 5, dpi = 300)

# --- 15f. Figure 20: P(synergy) heatmap --------------------------------------

p20 <- ggplot(bliss_summary,
              aes(x = factor(conc), y = temperature, fill = prob_synergy)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = sprintf("%.2f", prob_synergy)), size = 3.5) +
  scale_fill_gradient2(
    low = "#2166AC", mid = "white", high = "#B2182B",
    midpoint = 0.5, limits = c(0, 1),
    name = "P(synergy)"
  ) +
  labs(
    x = "Prothioconazole (mg/L)",
    y = "Temperature (\u00b0C)",
    title = "Posterior probability of synergy",
    subtitle = "P(\u0394 < 0) at each temperature \u00d7 dose combination"
  ) +
  theme_minimal(base_size = 13) +
  theme(panel.grid = element_blank()) +
  coord_equal()
ggsave(file.path(fig_dir, "20_bliss_prob_synergy_heatmap.png"),
       p20, width = 8, height = 7, dpi = 300)

# --- 15g. Figure 21: Mean delta heatmap --------------------------------------

p21 <- ggplot(bliss_summary,
              aes(x = factor(conc), y = temperature, fill = delta_mean)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = sprintf("%.3f", delta_mean)), size = 3.5) +
  scale_fill_gradient2(
    low = "#2166AC", mid = "white", high = "#B2182B",
    midpoint = 0, name = expression(bar(Delta))
  ) +
  labs(
    x = "Prothioconazole (mg/L)",
    y = "Temperature (\u00b0C)",
    title = "Mean Bliss deviation",
    subtitle = "Blue = synergy (\u0394 < 0), Red = antagonism (\u0394 > 0)"
  ) +
  theme_minimal(base_size = 13) +
  theme(panel.grid = element_blank()) +
  coord_equal()
ggsave(file.path(fig_dir, "21_bliss_delta_heatmap.png"),
       p21, width = 8, height = 7, dpi = 300)

# =============================================================================
# 16. MODEL COMPARISON — ADDITIVE NULL VS INTERACTION
# =============================================================================
# Null model: EC50 and n shared across temperatures (only r0 varies)
# Intermediate models: EC50 varies only, or n varies only
# Full model: both EC50 and n vary (= fit_hill, already fitted)

cat("\n--- Fitting additive null model ---\n")

formula_null <- bf(
  rate ~ r0 / (1 + (conc / ec50)^n),
  r0   ~ temperature + (1 | rep_id),
  ec50 ~ 1,
  n    ~ 1,
  nl = TRUE
)

priors_null <- c(
  prior(normal(0.04, 0.02), nlpar = "r0", lb = 0),
  prior(student_t(3, 0, 0.01), class = "sd", nlpar = "r0"),
  prior(lognormal(0, 1), nlpar = "ec50", lb = 0),
  prior(normal(1.5, 1), nlpar = "n", lb = 0.1)
)

fit_null <- brm(
  formula  = formula_null,
  data     = dat,
  prior    = priors_null,
  family   = gaussian(),
  chains   = 4,
  iter     = 4000,
  warmup   = 2000,
  cores    = 4,
  seed     = 42,
  control  = list(adapt_delta = 0.95, max_treedepth = 12),
  file     = file.path(mod_dir, "hill_growth_null_additive")
)

cat("\n--- Fitting EC50-only interaction model ---\n")

formula_ec50only <- bf(
  rate ~ r0 / (1 + (conc / ec50)^n),
  r0   ~ temperature + (1 | rep_id),
  ec50 ~ temperature,
  n    ~ 1,
  nl = TRUE
)

fit_ec50only <- brm(
  formula  = formula_ec50only,
  data     = dat,
  prior    = priors_hill,
  family   = gaussian(),
  chains   = 4,
  iter     = 4000,
  warmup   = 2000,
  cores    = 4,
  seed     = 42,
  control  = list(adapt_delta = 0.95, max_treedepth = 12),
  file     = file.path(mod_dir, "hill_growth_ec50only")
)

cat("\n--- Fitting n-only interaction model ---\n")

formula_nonly <- bf(
  rate ~ r0 / (1 + (conc / ec50)^n),
  r0   ~ temperature + (1 | rep_id),
  ec50 ~ 1,
  n    ~ temperature,
  nl = TRUE
)

fit_nonly <- brm(
  formula  = formula_nonly,
  data     = dat,
  prior    = priors_hill,
  family   = gaussian(),
  chains   = 4,
  iter     = 4000,
  warmup   = 2000,
  cores    = 4,
  seed     = 42,
  control  = list(adapt_delta = 0.95, max_treedepth = 12),
  file     = file.path(mod_dir, "hill_growth_nonly")
)

# --- 16a. Four-way LOO comparison ---------------------------------------------

cat("\n--- LOO comparison: additive null vs interaction models ---\n")

loo_null     <- loo(fit_null)
loo_ec50only <- loo(fit_ec50only)
loo_nonly    <- loo(fit_nonly)
loo_full     <- loo(fit_hill)

comp_interaction <- loo_compare(loo_null, loo_ec50only, loo_nonly, loo_full)
print(comp_interaction)

# Save comparison table
comp_tbl <- as_tibble(comp_interaction, rownames = "model")
write_csv(comp_tbl, file.path(tab_dir, "09_interaction_model_comparison.csv"))

cat("\nInterpretation:\n")
cat("  - If full model (temperature-varying EC50 + n) beats null → evidence for non-additive interaction\n")
cat("  - elpd_diff and SE quantify the strength of evidence\n")
cat("  - Intermediate models reveal whether interaction acts through EC50 (sensitivity), n (steepness), or both\n")

# --- 16b. Save model objects --------------------------------------------------

saveRDS(fit_null, file.path(mod_dir, "fit_null_additive.rds"))
saveRDS(fit_ec50only, file.path(mod_dir, "fit_ec50only.rds"))
saveRDS(fit_nonly, file.path(mod_dir, "fit_nonly.rds"))

# --- 16c. Figure 22: Model comparison point-and-interval plot -----------------

comp_tbl <- comp_tbl %>%
  mutate(
    model_label = case_when(
      grepl("null", model)     ~ "Null (additive)",
      grepl("ec50only", model) ~ "EC50 varies only",
      grepl("nonly", model)    ~ "n varies only",
      TRUE                     ~ "Full (EC50 + n vary)"
    ),
    model_label = fct_reorder(model_label, elpd_diff)
  )

p22 <- ggplot(comp_tbl, aes(x = elpd_diff, y = model_label)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_pointrange(aes(xmin = elpd_diff - 2 * se_diff,
                      xmax = elpd_diff + 2 * se_diff),
                  size = 0.8) +
  labs(
    x = expression(Delta*"ELPD (relative to best model)"),
    y = NULL,
    title = "Model comparison: additive vs interaction",
    subtitle = "LOO-CV; error bars = \u00b12 SE"
  ) +
  theme_bw(base_size = 13)
ggsave(file.path(fig_dir, "22_model_comparison_interaction.png"),
       p22, width = 8, height = 4, dpi = 300)

# =============================================================================
# 17. TEMPERATURE-SPECIFIC INTERACTION INDEX
# =============================================================================
# Integrate Bliss deviation across all concentrations at each temperature
# to get a single number summarising net interaction strength.
# I(T) = (1/c_max) × ∫₀^c_max [f_obs(c,T) - f_ref(c)] dc
# Computed via trapezoidal rule on posterior draws.

cat("\n--- Computing temperature-specific interaction index ---\n")

# Fine concentration grid for numerical integration
conc_fine <- seq(0, 4, length.out = 200)
c_max <- max(conc_fine)

interaction_index_draws <- param_draws %>%
  left_join(
    param_draws %>%
      filter(temperature == T_ref) %>%
      select(.draw, ec50_ref = ec50, n_ref = n_hill),
    by = ".draw"
  ) %>%
  # For each draw × temperature, integrate delta over concentration
  group_by(.draw, temperature) %>%
  summarise(
    I_T = {
      f_obs_vec <- 1 / (1 + (conc_fine / ec50[1])^n_hill[1])
      f_ref_vec <- 1 / (1 + (conc_fine / ec50_ref[1])^n_ref[1])
      delta_vec <- f_obs_vec - f_ref_vec
      # Trapezoidal integration, normalised by c_max
      (1 / c_max) * sum(diff(conc_fine) * (head(delta_vec, -1) + tail(delta_vec, -1)) / 2)
    },
    .groups = "drop"
  )

# --- 17a. Summarise interaction index -----------------------------------------

interaction_summary <- interaction_index_draws %>%
  group_by(temperature) %>%
  summarise(
    I_mean   = mean(I_T),
    I_median = median(I_T),
    I_lower  = quantile(I_T, 0.025),
    I_upper  = quantile(I_T, 0.975),
    prob_synergy    = mean(I_T < 0),
    prob_antagonism = mean(I_T > 0),
    .groups = "drop"
  )

cat("\n--- Interaction index by temperature ---\n")
cat("I(T) < 0 → net synergy (fungicide more effective than expected)\n")
cat("I(T) > 0 → net antagonism. Reference temperature:", as.character(T_ref), "°C\n\n")
print(interaction_summary)
write_csv(interaction_summary, file.path(tab_dir, "10_interaction_index_by_temperature.csv"))

# --- 17b. Figure 23: Interaction index point-and-interval plot ----------------

p23 <- ggplot(interaction_summary, aes(x = temperature, y = I_median,
                                        colour = temperature)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_pointinterval(aes(ymin = I_lower, ymax = I_upper), size = 3) +
  scale_colour_manual(values = temp_pal, guide = "none") +
  annotate("text", x = 0.5, y = -Inf, label = "Synergy \u2193",
           hjust = 0, vjust = -0.5, size = 3, colour = "grey40", fontface = "italic") +
  annotate("text", x = 0.5, y = Inf, label = "Antagonism \u2191",
           hjust = 0, vjust = 1.5, size = 3, colour = "grey40", fontface = "italic") +
  labs(
    x = "Temperature (\u00b0C)",
    y = "Interaction index I(T)",
    title = "Temperature-specific interaction index (Bliss framework)",
    subtitle = paste0("Integrated Bliss deviation across doses. Reference: ", T_ref, "\u00b0C")
  ) +
  theme_bw(base_size = 13)
ggsave(file.path(fig_dir, "23_interaction_index_by_temperature.png"),
       p23, width = 7, height = 5, dpi = 300)

# --- 17c. Figure 24: Posterior density of I(T) by temperature -----------------

p24 <- ggplot(interaction_index_draws, aes(x = I_T, fill = temperature,
                                            colour = temperature)) +
  geom_density(alpha = 0.3) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  scale_fill_manual(values = temp_pal, name = "Temp (\u00b0C)") +
  scale_colour_manual(values = temp_pal, name = "Temp (\u00b0C)") +
  labs(
    x = "Interaction index I(T)",
    y = "Posterior density",
    title = "Posterior distribution of interaction index by temperature",
    subtitle = "Dashed line = additivity (I = 0)"
  ) +
  theme_bw(base_size = 13)
ggsave(file.path(fig_dir, "24_interaction_index_density.png"),
       p24, width = 8, height = 5, dpi = 300)

# --- 17d. Figure 25: P(synergy) vs P(antagonism) stacked bar chart -----------

interaction_bar <- interaction_summary %>%
  select(temperature, prob_synergy, prob_antagonism) %>%
  pivot_longer(cols = c(prob_synergy, prob_antagonism),
               names_to = "direction", values_to = "probability") %>%
  mutate(direction = ifelse(direction == "prob_synergy", "Synergy", "Antagonism"),
         direction = factor(direction, levels = c("Antagonism", "Synergy")))

p25 <- ggplot(interaction_bar, aes(x = temperature, y = probability,
                                    fill = direction)) +
  geom_col(width = 0.7) +
  geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey50") +
  scale_fill_manual(values = c("Synergy" = "#2166AC", "Antagonism" = "#B2182B"),
                    name = "Direction") +
  labs(
    x = "Temperature (\u00b0C)",
    y = "Posterior probability",
    title = "Evidence for synergy vs antagonism by temperature",
    subtitle = "Stacked bars; dashed line = equal probability"
  ) +
  theme_bw(base_size = 13)
ggsave(file.path(fig_dir, "25_interaction_index_prob_barplot.png"),
       p25, width = 7, height = 5, dpi = 300)

# =============================================================================
# 18. COMBINED INTERACTION SUMMARY
# =============================================================================

bliss_full_summary <- interaction_summary %>%
  mutate(
    interpretation = case_when(
      prob_synergy > 0.95    ~ "Strong synergy",
      prob_synergy > 0.80    ~ "Moderate synergy",
      prob_antagonism > 0.95 ~ "Strong antagonism",
      prob_antagonism > 0.80 ~ "Moderate antagonism",
      TRUE                   ~ "Additive / inconclusive"
    )
  )

cat("\n--- Combined interaction summary ---\n")
print(bliss_full_summary)
cat("\nLOO model comparison:\n")
print(comp_interaction)
write_csv(bliss_full_summary, file.path(tab_dir, "11_interaction_summary.csv"))

cat("\n--- Bliss independence analysis complete ---\n")

# #############################################################################
#
#  RESPIRATION RATE ANALYSIS — POWER LAW
#
# #############################################################################
# Model: log(resp) = log(a) + b * log(conc), i.e. resp = a * conc^b
# Fitted as a linear model on the log-log scale using non-zero concentrations.
# The intercept log(a) = log(respiration at 1 mg/L).
# The slope b = power-law exponent (rate of increase with dose).
# Both vary by temperature; replicate random intercepts capture between-
# replicate variation.

cat("\n\n========== RESPIRATION RATE ANALYSIS ==========\n\n")

# =============================================================================
# R1. DATA WRANGLING — RESPIRATION
# =============================================================================

dat_resp <- dat_resp %>%
  mutate(
    conc_plot = ifelse(conc == 0, min_nonzero_conc / 2, conc),
    log_rate  = log(rate)
  )

resp_summ <- dat_resp %>%
  group_by(temperature, conc) %>%
  summarise(
    n         = n(),
    mean_rate = mean(rate, na.rm = TRUE),
    se_rate   = sd(rate, na.rm = TRUE) / sqrt(n()),
    .groups   = "drop"
  ) %>%
  mutate(conc_plot = ifelse(conc == 0, min_nonzero_conc / 2, conc))

cat("\n--- Respiration data summary (n per temperature x concentration) ---\n")
print(resp_summ %>% select(temperature, conc, n) %>% pivot_wider(
  names_from = conc, values_from = n
))

write_csv(resp_summ, file.path(tab_dir_resp, "01_data_summary.csv"))

# =============================================================================
# R2. EXPLORATORY PLOTS — RESPIRATION
# =============================================================================

r01 <- ggplot(dat_resp, aes(x = conc, y = rate, colour = temperature)) +
  geom_point(alpha = 0.5, size = 2) +
  geom_point(data = resp_summ, aes(y = mean_rate), size = 3, shape = 18) +
  geom_line(data = resp_summ, aes(y = mean_rate), linewidth = 0.5, alpha = 0.5) +
  scale_colour_manual(values = temp_pal, name = "Temperature (°C)") +
  labs(
    x = "Prothioconazole (mg/L)",
    y = expression("Respiration rate (C C"^{-1}~"h"^{-1}*")"),
    title = "Dose-response: respiration rate vs prothioconazole concentration",
    subtitle = "All replicates shown; diamonds = means"
  ) +
  theme_bw(base_size = 13)
ggsave(file.path(fig_dir_resp, "01_resp_vs_conc_linear.png"),
       r01, width = 8, height = 5, dpi = 300)

r02 <- ggplot(dat_resp, aes(x = conc_plot, y = rate, colour = temperature)) +
  geom_point(alpha = 0.6, size = 2) +
  scale_x_log10() + scale_y_log10() +
  scale_colour_manual(values = temp_pal, name = "Temperature (°C)") +
  labs(
    x = "Prothioconazole concentration (mg/L, log scale)",
    y = expression("Respiration rate (C C"^{-1}~"h"^{-1}*", log scale)"),
    title = "Respiration dose-response in log-log space",
    subtitle = "Zero-concentration controls plotted at 0.03 mg/L"
  ) +
  theme_bw(base_size = 13)
ggsave(file.path(fig_dir_resp, "02_resp_vs_conc_loglog.png"),
       r02, width = 8, height = 5, dpi = 300)

r03 <- ggplot(dat_resp, aes(x = conc, y = rate)) +
  geom_point(aes(colour = temperature), alpha = 0.5, size = 1.5) +
  geom_point(data = resp_summ, aes(y = mean_rate), size = 2.5, shape = 18) +
  geom_errorbar(data = resp_summ,
                aes(y = mean_rate,
                    ymin = mean_rate - se_rate,
                    ymax = mean_rate + se_rate),
                width = 0.1) +
  facet_wrap(~ temperature, ncol = 4, labeller = label_both) +
  scale_colour_manual(values = temp_pal, guide = "none") +
  labs(
    x = "Prothioconazole (mg/L)",
    y = expression("Respiration rate (C C"^{-1}~"h"^{-1}*")"),
    title = "Respiration dose-response faceted by temperature",
    subtitle = "Points = replicates, diamonds = means +/- SE"
  ) +
  theme_bw(base_size = 12)
ggsave(file.path(fig_dir_resp, "03_resp_vs_conc_faceted.png"),
       r03, width = 12, height = 6, dpi = 300)

r04 <- ggplot(dat_resp, aes(x = conc_plot, y = rate)) +
  geom_point(aes(colour = temperature), alpha = 0.5, size = 1.5) +
  geom_point(data = resp_summ, aes(y = mean_rate), size = 2.5, shape = 18) +
  scale_x_log10() + scale_y_log10() +
  facet_wrap(~ temperature, ncol = 4, labeller = label_both) +
  scale_colour_manual(values = temp_pal, guide = "none") +
  labs(
    x = "Prothioconazole concentration (mg/L, log scale)",
    y = expression("Respiration rate (log scale)"),
    title = "Respiration dose-response faceted by temperature (log-log)",
    subtitle = "Points = replicates, diamonds = means"
  ) +
  theme_bw(base_size = 12)
ggsave(file.path(fig_dir_resp, "04_resp_vs_conc_loglog_faceted.png"),
       r04, width = 12, height = 6, dpi = 300)

# =============================================================================
# R3. FIT LOG-LINEAR MODEL
# =============================================================================
# log(resp) = intercept + slope * conc
# Equivalent to: resp = a * exp(b * conc)
#   a = exp(intercept) = baseline respiration at zero concentration
#   b = slope = rate of change of log-respiration per mg/L
# All data including zero-concentration controls are used.
# slope > 0 means respiration increases with dose.

cat("\n--- Fitting log-linear model for respiration ---\n")

fit_resp <- brm(
  log_rate ~ conc * temperature + (1 | rep_id),
  data     = dat_resp,
  family   = gaussian(),
  chains   = 4,
  iter     = 4000,
  warmup   = 2000,
  cores    = 4,
  seed     = 42,
  file     = file.path(mod_dir_resp, "loglinear_respiration")
)

cat("\n--- Log-linear Model Summary ---\n")
print(summary(fit_resp))

resp_fixed <- as_tibble(fixef(fit_resp), rownames = "parameter")
write_csv(resp_fixed, file.path(tab_dir_resp, "02_loglinear_parameters.csv"))

# =============================================================================
# R4. EXTRACT INTERCEPT AND SLOPE BY TEMPERATURE
# =============================================================================

# Temperature contrasts: 15°C is the reference level
# Intercept at each temperature = b_Intercept + b_temperatureXX
#   → exp(intercept) = baseline respiration at conc = 0
# Slope at each temperature = b_conc + b_conc:temperatureXX
#   → rate of increase of log(resp) per mg/L prothioconazole

resp_draws <- as_draws_df(fit_resp)

intercept_draws <- tibble(
  .draw = seq_len(nrow(resp_draws)),
  int_15 = resp_draws$b_Intercept,
  int_18 = resp_draws$b_Intercept + resp_draws$b_temperature18,
  int_21 = resp_draws$b_Intercept + resp_draws$b_temperature21,
  int_24 = resp_draws$b_Intercept + resp_draws$b_temperature24,
  int_26 = resp_draws$b_Intercept + resp_draws$b_temperature26,
  int_27 = resp_draws$b_Intercept + resp_draws$b_temperature27,
  int_28 = resp_draws$b_Intercept + resp_draws$b_temperature28
) %>%
  pivot_longer(-`.draw`, names_to = "temperature", values_to = "intercept",
               names_prefix = "int_") %>%
  mutate(temperature = factor(temperature,
                              levels = as.character(c(15, 18, 21, 24, 26, 27, 28))))

slope_draws <- tibble(
  .draw = seq_len(nrow(resp_draws)),
  slp_15 = resp_draws$b_conc,
  slp_18 = resp_draws$b_conc + resp_draws$`b_conc:temperature18`,
  slp_21 = resp_draws$b_conc + resp_draws$`b_conc:temperature21`,
  slp_24 = resp_draws$b_conc + resp_draws$`b_conc:temperature24`,
  slp_26 = resp_draws$b_conc + resp_draws$`b_conc:temperature26`,
  slp_27 = resp_draws$b_conc + resp_draws$`b_conc:temperature27`,
  slp_28 = resp_draws$b_conc + resp_draws$`b_conc:temperature28`
) %>%
  pivot_longer(-`.draw`, names_to = "temperature", values_to = "slope",
               names_prefix = "slp_") %>%
  mutate(temperature = factor(temperature,
                              levels = as.character(c(15, 18, 21, 24, 26, 27, 28))))

# Summaries
intercept_summary <- intercept_draws %>%
  group_by(temperature) %>%
  median_qi(intercept, .width = 0.95)

slope_summary <- slope_draws %>%
  group_by(temperature) %>%
  median_qi(slope, .width = 0.95)

# Back-transform intercept to baseline rate: a = exp(intercept)
a_draws <- intercept_draws %>%
  mutate(a = exp(intercept))

a_summary <- a_draws %>%
  group_by(temperature) %>%
  median_qi(a, .width = 0.95)

# Combined parameter table
resp_param_tbl <- slope_summary %>%
  select(temperature, slope, .lower, .upper) %>%
  rename(slope_lower = .lower, slope_upper = .upper) %>%
  left_join(
    a_summary %>%
      select(temperature, a, .lower, .upper) %>%
      rename(a_lower = .lower, a_upper = .upper),
    by = "temperature"
  )

cat("\n--- Log-linear parameters by temperature ---\n")
cat("slope = b (per mg/L); a = exp(intercept) = baseline resp at conc = 0\n\n")
print(resp_param_tbl)
write_csv(resp_param_tbl, file.path(tab_dir_resp, "03_loglinear_params_by_temperature.csv"))

# =============================================================================
# R5. PARAMETER PLOTS
# =============================================================================

# --- 05: Slope by temperature -------------------------------------------------

r05 <- ggplot(slope_summary, aes(x = temperature, y = slope,
                                  colour = temperature)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_pointinterval(aes(ymin = .lower, ymax = .upper), size = 3) +
  scale_colour_manual(values = temp_pal, guide = "none") +
  labs(
    x = "Temperature (°C)",
    y = expression("Slope (per mg/L)"),
    title = "Dose-response slope by temperature",
    subtitle = "b > 0: respiration increases with dose; posterior median and 95% CrI"
  ) +
  theme_bw(base_size = 13)
ggsave(file.path(fig_dir_resp, "05_slope_by_temperature.png"),
       r05, width = 7, height = 5, dpi = 300)

# --- 06: Slope posterior densities -------------------------------------------

r06 <- ggplot(slope_draws, aes(x = slope, fill = temperature,
                                colour = temperature)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_density(alpha = 0.35, linewidth = 0.6) +
  scale_fill_manual(values = temp_pal, name = "Temperature (°C)") +
  scale_colour_manual(values = temp_pal, name = "Temperature (°C)") +
  labs(
    x = expression("Slope (per mg/L)"),
    y = "Posterior density",
    title = "Posterior distributions of dose-response slope across temperatures"
  ) +
  theme_bw(base_size = 13)
ggsave(file.path(fig_dir_resp, "06_slope_posterior_density.png"),
       r06, width = 8, height = 5, dpi = 300)

# --- 07: Baseline respiration by temperature ----------------------------------

r07 <- ggplot(a_summary, aes(x = temperature, y = a, colour = temperature)) +
  geom_pointinterval(aes(ymin = .lower, ymax = .upper), size = 3) +
  scale_colour_manual(values = temp_pal, guide = "none") +
  labs(
    x = "Temperature (°C)",
    y = expression("Baseline respiration (C C"^{-1}~"h"^{-1}*")"),
    title = "Baseline respiration rate by temperature",
    subtitle = "exp(intercept) = respiration at zero concentration; posterior median and 95% CrI"
  ) +
  theme_bw(base_size = 13)
ggsave(file.path(fig_dir_resp, "07_baseline_by_temperature.png"),
       r07, width = 7, height = 5, dpi = 300)

# --- 08: Baseline posterior densities -----------------------------------------

r08 <- ggplot(a_draws, aes(x = a, fill = temperature,
                            colour = temperature)) +
  geom_density(alpha = 0.35, linewidth = 0.6) +
  scale_fill_manual(values = temp_pal, name = "Temperature (°C)") +
  scale_colour_manual(values = temp_pal, name = "Temperature (°C)") +
  labs(
    x = expression("Baseline respiration (C C"^{-1}~"h"^{-1}*")"),
    y = "Posterior density",
    title = "Posterior distributions of baseline respiration"
  ) +
  theme_bw(base_size = 13)
ggsave(file.path(fig_dir_resp, "08_baseline_posterior_density.png"),
       r08, width = 8, height = 5, dpi = 300)

# =============================================================================
# R6. MODEL PREDICTIONS
# =============================================================================

# Prediction grid including zero concentration
resp_pred_grid <- expand_grid(
  temperature = factor(sort(unique(as.numeric(levels(dat_resp$temperature)))),
                       levels = levels(dat_resp$temperature)),
  conc = seq(0, 4, length.out = 200)
)

# Predicted log_rate (model scale) and back-transformed rate
resp_pred_draws <- resp_pred_grid %>%
  add_epred_draws(fit_resp, ndraws = 500, re_formula = NA)

# Summary on log scale (for log-resp plots)
resp_pred_log_summary <- resp_pred_draws %>%
  group_by(temperature, conc) %>%
  median_qi(.epred, .width = 0.95)

# Summary on natural scale (back-transformed)
resp_pred_summary <- resp_pred_draws %>%
  mutate(pred_rate = exp(.epred)) %>%
  group_by(temperature, conc) %>%
  median_qi(pred_rate, .width = 0.95)

# --- 09: Overlay, log(respiration) vs concentration ---------------------------

r09 <- ggplot() +
  geom_ribbon(data = resp_pred_log_summary,
              aes(x = conc, ymin = .lower, ymax = .upper,
                  fill = temperature), alpha = 0.15) +
  geom_line(data = resp_pred_log_summary,
            aes(x = conc, y = .epred, colour = temperature),
            linewidth = 0.8) +
  geom_point(data = dat_resp,
             aes(x = conc, y = log_rate, colour = temperature),
             alpha = 0.5, size = 1.8) +
  scale_colour_manual(values = temp_pal, name = "Temperature (°C)") +
  scale_fill_manual(values = temp_pal, name = "Temperature (°C)") +
  labs(
    x = "Prothioconazole (mg/L)",
    y = expression("log(Respiration rate)"),
    title = expression("log(Respiration rate) vs concentration"),
    subtitle = "Log-linear model — posterior median and 95% CrI"
  ) +
  theme_bw(base_size = 13)
ggsave(file.path(fig_dir_resp, "09_loglinear_fit_logscale.png"),
       r09, width = 8, height = 5, dpi = 300)

# --- 10: Faceted, log(respiration) vs concentration ---------------------------

r10 <- ggplot() +
  geom_ribbon(data = resp_pred_log_summary,
              aes(x = conc, ymin = .lower, ymax = .upper),
              fill = "steelblue", alpha = 0.2) +
  geom_line(data = resp_pred_log_summary,
            aes(x = conc, y = .epred),
            colour = "steelblue", linewidth = 0.8) +
  geom_point(data = dat_resp, aes(x = conc, y = log_rate),
             alpha = 0.6, size = 1.5) +
  facet_wrap(~ temperature, ncol = 4, labeller = label_both) +
  labs(
    x = "Prothioconazole (mg/L)",
    y = expression("log(Respiration rate)"),
    title = "Log-linear fit by temperature",
    subtitle = "Posterior median and 95% CrI"
  ) +
  theme_bw(base_size = 12)
ggsave(file.path(fig_dir_resp, "10_loglinear_fit_logscale_faceted.png"),
       r10, width = 12, height = 6, dpi = 300)

# =============================================================================
# R7. POSTERIOR PREDICTIVE CHECK
# =============================================================================

r11 <- pp_check(fit_resp, ndraws = 100) +
  labs(title = "Posterior predictive check — log-linear respiration model") +
  theme_bw(base_size = 13)
ggsave(file.path(fig_dir_resp, "11_pp_check.png"),
       r11, width = 7, height = 5, dpi = 300)

# =============================================================================
# R8. SUMMARY TABLE
# =============================================================================

slope_report <- slope_draws %>%
  group_by(temperature) %>%
  summarise(
    slope_mean = mean(slope), slope_sd = sd(slope),
    slope_median = median(slope),
    slope_lower = quantile(slope, 0.025),
    slope_upper = quantile(slope, 0.975), .groups = "drop"
  )

a_report <- a_draws %>%
  group_by(temperature) %>%
  summarise(
    baseline_mean = mean(a), baseline_sd = sd(a),
    baseline_median = median(a),
    baseline_lower = quantile(a, 0.025),
    baseline_upper = quantile(a, 0.975), .groups = "drop"
  )

resp_summary_table <- slope_report %>%
  left_join(a_report, by = "temperature")

cat("\n--- Respiration log-linear parameter summary ---\n")
cat("slope = b (per mg/L); baseline = exp(intercept) = resp at conc = 0\n\n")
print(resp_summary_table)
write_csv(resp_summary_table, file.path(tab_dir_resp, "04_loglinear_summary_table.csv"))

# =============================================================================
# R9. COMPUTE GLOBAL MEAN SLOPE (AVERAGED ACROSS TEMPERATURES)
# =============================================================================

global_slope_draws <- slope_draws %>%
  group_by(.draw) %>%
  summarise(slope = mean(slope), .groups = "drop")

global_slope_mean <- mean(global_slope_draws$slope)
global_slope_sd   <- sd(global_slope_draws$slope)
global_slope_qi   <- quantile(global_slope_draws$slope, c(0.025, 0.975))

cat("\n--- Global mean slope (averaged across temperatures) ---\n")
cat(sprintf("slope = %.3f +/- %.3f [%.3f, %.3f]\n",
            global_slope_mean, global_slope_sd,
            global_slope_qi[1], global_slope_qi[2]))

# =============================================================================
# R10. RESPIRATION COMPOSITE PANELS (saved for combined figure)
# =============================================================================

# Panel: Log-linear fit, log(resp) vs conc
pa_resp <- ggplot() +
  geom_ribbon(data = resp_pred_log_summary,
              aes(x = conc, ymin = .lower, ymax = .upper,
                  fill = temperature), alpha = 0.15) +
  geom_line(data = resp_pred_log_summary,
            aes(x = conc, y = .epred, colour = temperature),
            linewidth = 0.8) +
  geom_point(data = dat_resp,
             aes(x = conc, y = log_rate, colour = temperature),
             alpha = 0.5, size = 1.8) +
  scale_colour_manual(values = temp_pal, name = "Temperature (°C)") +
  scale_fill_manual(values = temp_pal, name = "Temperature (°C)") +
  labs(
    x = "Prothioconazole (mg/L)",
    y = expression("log(Respiration rate)")
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    legend.position = "none"
  )

# Panel: Slope by temperature as forest plot (mean +/- SD) + global mean
slope_report_plot <- slope_report %>%
  select(temperature, slope_mean, slope_sd) %>%
  bind_rows(tibble(
    temperature = "Global",
    slope_mean  = global_slope_mean,
    slope_sd    = global_slope_sd
  )) %>%
  mutate(temperature = factor(temperature,
                              levels = c("15", "18", "21", "24", "26", "27", "28", "Global")))

pb_resp <- ggplot(slope_report_plot, aes(x = temperature, y = slope_mean,
                                          colour = temperature)) +
  geom_vline(xintercept = 7.5, linetype = "dashed", colour = "grey70") +
  geom_hline(yintercept = 0, linetype = "dotted", colour = "grey70") +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = slope_mean - slope_sd,
                    ymax = slope_mean + slope_sd),
                width = 0.2, linewidth = 0.8) +
  geom_errorbar(aes(ymin = slope_mean - 2 * slope_sd,
                    ymax = slope_mean + 2 * slope_sd),
                width = 0.1, linewidth = 0.4, linetype = "dashed") +
  scale_colour_manual(
    values = c(temp_pal, "Global" = "black"),
    guide = "none"
  ) +
  labs(
    x = "Temperature (°C)",
    y = expression("Slope (per mg/L)")
  ) +
  theme_bw(base_size = 11) +
  theme(panel.grid = element_blank())

# =============================================================================
# SAVE RESPIRATION MODEL
# =============================================================================

saveRDS(fit_resp, file.path(mod_dir_resp, "fit_loglinear_respiration.rds"))

cat("\n--- Respiration analysis complete ---\n")

# =============================================================================
# COMBINED COMPOSITE FIGURE — GROWTH + RESPIRATION
# =============================================================================

# Add legend back to top-left panel only
pa_growth_leg <- pa_growth +
  theme(
    legend.position = "top",
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 8),
    legend.key.size = unit(0.4, "cm")
  ) +
  guides(fill = "none",
         colour = guide_legend(nrow = 1, override.aes = list(size = 2)))

library(cowplot)

combined <- plot_grid(
  pa_growth_leg, pb_growth,
  pa_resp, pb_resp,
  ncol = 2, rel_widths = c(2, 1),
  labels = c("a", "b", "c", "d"),
  label_size = 14, label_fontface = "bold"
)

ggsave(file.path(base_dir, "figures", "combined_composite.png"),
       combined, width = 12, height = 10, dpi = 300)

cat("\n--- Combined composite figure saved ---\n")
cat("\n--- All analyses done! ---\n")
