# Zymoseptoria Analysis Pipeline

## Overview

This project analyses how temperature modulates fungicide (prothioconazole) sensitivity
in *Zymoseptoria tritici* growth and respiration rates. The pipeline processes raw oxygen
time-series data through Bayesian thermal performance curve modelling and dose-response
analysis.

The analysis is split into **6 numbered scripts** that run sequentially, plus a shared
configuration file and a master runner.

## Script Map

| Script | What it does | Run time | Key outputs |
|---|---|---|---|
| `00_config.R` | Shared paths, constants, helper functions | Instant (sourced by all) | -- |
| `01_longdata.R` | Read raw oxygen CSVs, pivot wide to long | Seconds | `Oxygen_All_Long.csv` |
| `02_trimming.R` | Spline-based trimming with manual overrides | Minutes | `Oxygen_Data_Filtered.csv`, diagnostics PDF |
| `03_oxygen_fits.R` | nlsLM fits, carbon unit conversions, descriptive plots | Minutes | `derived_N0_R_results_with_carbon.csv`, `.rds` intermediates |
| `04_bayesian_models.R` | Arrhenius + Sharpe-Schoolfield + hierarchical Bayesian fits | **Hours** (MCMC) | `.rds` model files, posterior CSVs, `stage2_workspace.RData` |
| `05_plots_and_effects.R` | All Bayesian plots, effect sizes, synergy, publication figures | Minutes | All `.png` figures, `fig1_extra_dose_response_by_temp.csv` |
| `06_dose_response_brms.R` | Hill inhibition + log-linear dose-response models | Hours (MCMC) | EC50/pEC50 tables, dose-response figures |

## How to Run

### Full pipeline:
```r
source("run_all.R")
```

### Or step by step:
```r
source("01_longdata.R")
source("02_trimming.R")
source("03_oxygen_fits.R")
source("04_bayesian_models.R")      # slow: hours of MCMC
source("05_plots_and_effects.R")
source("06_dose_response_brms.R")   # slow: hours of MCMC
```

### Re-make plots without re-fitting models:
```r
source("05_plots_and_effects.R")    # uses saved stage2_workspace.RData
```

### Run only dose-response analysis (requires fig1_extra CSV from 05):
```r
source("06_dose_response_brms.R")
```

## Dependencies Between Scripts

```
00_config.R  (sourced by every script)
     |
01_longdata.R
     |  writes: Oxygen_All_Long.csv
     v
02_trimming.R
     |  writes: Oxygen_Data_Filtered.csv, Oxygen_Trimmed_Series_Metadata.csv
     v
03_oxygen_fits.R
     |  writes: derived_N0_R_results_with_carbon.csv
     |  saves:  results_plot.rds, ratio_dat.rds, dose_levels.rds,
     |          dose_cols.rds, dose_key_tbl.rds, growth_dat.rds,
     |          resp_dat.rds, ratio_dat_arr.rds, cue_dat.rds,
     |          growth_biomass_dat.rds, resp_biomass_dat.rds
     v
04_bayesian_models.R
     |  loads:  .rds intermediates from 03
     |  saves:  stage2_workspace.RData, all model .rds, posterior CSVs
     v
05_plots_and_effects.R
     |  loads:  stage2_workspace.RData
     |  writes: all figures, effect size CSVs,
     |          fig1_extra_dose_response_by_temp.csv
     v
06_dose_response_brms.R
     |  reads:  fig1_extra_dose_response_by_temp.csv
     |  writes: Hill/log-linear model outputs in tables/growth, tables/respiration
```

## What to Edit

- **`base_dir`** in `00_config.R` -- set to your project root
- **`BAYES_*` settings** in `00_config.R` -- MCMC iterations, chains, etc.
- **`RMSE_KEEP_THRESHOLD`** in `00_config.R` -- quality filter for oxygen fits
- **`EXCLUDE_POINTS`** in `03_oxygen_fits.R` -- the outlier exclusion list
- **`manual_trim_overrides`** in `02_trimming.R` -- manual curve trimming overrides
