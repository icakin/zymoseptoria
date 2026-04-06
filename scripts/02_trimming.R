# ──────────────────────────────────────────────────────────────────────────
# Oxygen-time-series trimming by main spline curve shape
#
# What this script does:
#   1. Reads Oxygen_All_Long.csv
#   2. Creates a short code for each curve (C001, C002, ...)
#   3. Automatically trims each curve
#   4. Lets you manually override start and/or end using curve code
#   5. Saves trimmed data and a PDF showing diagnostics
#
# Manual override rule:
#   - If a curve code is NOT entered below, automatic trimming is used.
#   - If a curve code IS entered:
#       * keep_start filled  -> manual start
#       * keep_start = NA    -> automatic start
#       * keep_end filled    -> manual end
#       * keep_end = NA      -> automatic end
# ──────────────────────────────────────────────────────────────────────────

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
  library(zoo)
  library(stringr)
  library(scales)
})

# ── Main knobs ───────────────────────────────────────────────────────────

TRIM_SPAR <- 0.8

MIN_ROWS_SERIES <- 8
MIN_RANGE_O2    <- 0.05

# Main descending-region detection
NEG_SLOPE_EPS     <- 0.00015
MIN_RUN_POINTS    <- 6
MAX_GAP_POINTS    <- 2
SEARCH_START_FRAC <- 0.08
SEARCH_END_FRAC   <- 0.98

# Start selection
START_BACKOFF_POINTS <- 0
PRETRANS_SLOPE_EPS   <- 0.0010
PRETRANS_MIN_RUN     <- 5
SLOPE_TREND_EPS      <- 0.00001
START_SMOOTH_K       <- 5

# Stable-start controls
START_LOOKBACK_POINTS  <- 120
START_STABLE_RUN_LEN   <- 10
START_NEG_FRAC         <- 0.90
START_WIGGLE_TOL       <- 0.00002
START_MAX_SIGN_CHANGES <- 1

# Protection against very steep initial transition
START_LEVEL_RUN_LEN    <- 8
START_SLOPE_RATIO_MAX  <- 1.03

# End detection after steepest drop
END_SMOOTH_K     <- 5
END_MIN_GAP      <- 2
END_RECOV_RUN    <- 2
END_RECOV_FRAC   <- 0.99
END_MONO_FRAC    <- 0.40
END_SLOPE_EPS    <- 0.00005

# PDF axis detail
X_MAJOR_BREAKS_N <- 12
Y_MAJOR_BREAKS_N <- 10

# ── USER-FRIENDLY MANUAL OVERRIDE SECTION ────────────────────────────────
# HOW TO USE:
# 1. Run the script once
# 2. Open Oxygen_Curve_Code_Key.csv and oxygen_trimming_diagnostics.pdf
# 3. Find the curve code you want to adjust, for example C014
# 4. Enter the code below
# 5. Fill only the side you want to change:
#      keep_start = new start time
#      keep_end   = new end time
# 6. Leave NA if the algorithm is already correct for that side
#
# EXAMPLES:
# - Change only start:
#     "C014", 12, NA
# - Change only end:
#     "C027", NA, 80
# - Change both:
#     "C031", 10, 75
#
# IMPORTANT:
# - NA means: keep the automatic trim for that side
# - If the whole table is empty, all curves use automatic trimming

manual_trim_overrides <- tribble(
  ~curve_code, ~keep_start, ~keep_end,
  "C014",      NA,          2500,
  "C041",      NA,          3500,
  "C099",      500,         1250,
  "C101",      250,         900,
  "C103",      NA,          1600,
  "C105",      250,         1250,
  "C108",      NA,          1000,
  "C122",      NA,          1000,
  "C123",      250,         1750,
  "C127",      750,         NA,
  "C128",      NA,          1600,
  "C129",      NA,          1600,
  "C133",      NA,          1600,
  "C134",      NA,          1100,
  "C135",      NA,          1250,
  "C138",      NA,          1500,
  "C139",      500,         1000,
  "C140",      500,         1500,
  "C141",      500,         1300,
  "C145",      NA,          1500,
  "C146",      1000,        2500,
  "C150",      NA,          1000,
  "C151",      NA,          1000,
  "C152",      NA,          2000,
  "C155",      NA,          1250,
  "C157",      NA,          1000,
  "C156",      750,         1500,
  "C158",      NA,          1500
)

# ── Helpers ──────────────────────────────────────────────────────────────

safe_spline_fit <- function(time, oxygen, spar = TRIM_SPAR) {
  fit <- tryCatch(
    smooth.spline(time, oxygen, spar = spar),
    error = function(e) NULL
  )
  if (is.null(fit)) return(rep(NA_real_, length(time)))
  
  tryCatch(
    predict(fit, x = time)$y,
    error = function(e) rep(NA_real_, length(time))
  )
}

get_slopes <- function(time, y) {
  dt <- diff(time)
  dy <- diff(y)
  s <- dy / dt
  s[!is.finite(s)] <- NA_real_
  s
}

smooth_slopes <- function(slopes, k = START_SMOOTH_K) {
  out <- slopes
  if (length(slopes) >= k) {
    sm <- zoo::rollmedian(slopes, k = k, fill = NA, align = "center")
    sm[is.na(sm)] <- slopes[is.na(sm)]
    out <- sm
  }
  out
}

rle_to_runs <- function(flag_vec) {
  r <- rle(flag_vec)
  ends <- cumsum(r$lengths)
  starts <- c(1, head(ends, -1) + 1)
  
  tibble(
    value = r$values,
    start = starts,
    end   = ends,
    len   = r$lengths
  )
}

merge_descending_runs <- function(run_tbl, max_gap = MAX_GAP_POINTS) {
  if (nrow(run_tbl) == 0) return(run_tbl)
  
  desc_runs <- run_tbl %>% filter(value)
  
  if (nrow(desc_runs) <= 1) {
    return(desc_runs %>% select(start, end, len))
  }
  
  out <- list()
  cur_start <- desc_runs$start[1]
  cur_end   <- desc_runs$end[1]
  
  for (i in 2:nrow(desc_runs)) {
    gap <- desc_runs$start[i] - cur_end - 1
    if (gap <= max_gap) {
      cur_end <- desc_runs$end[i]
    } else {
      out[[length(out) + 1]] <- tibble(start = cur_start, end = cur_end)
      cur_start <- desc_runs$start[i]
      cur_end   <- desc_runs$end[i]
    }
  }
  
  out[[length(out) + 1]] <- tibble(start = cur_start, end = cur_end)
  
  bind_rows(out) %>%
    mutate(len = end - start + 1L)
}

score_descending_runs <- function(time, o2_fit, slopes, run_tbl) {
  if (nrow(run_tbl) == 0) return(run_tbl)
  
  run_tbl %>%
    rowwise() %>%
    mutate(
      time_start = time[start],
      time_end   = time[min(end + 1L, length(time))],
      duration   = time_end - time_start,
      y_start    = o2_fit[start],
      y_end      = o2_fit[min(end + 1L, length(o2_fit))],
      total_drop = y_start - y_end,
      mean_slope = mean(slopes[start:end], na.rm = TRUE),
      score      = pmax(duration, 0) * pmax(total_drop, 0)
    ) %>%
    ungroup()
}

find_main_descending_run <- function(time, o2_fit) {
  slopes <- get_slopes(time, o2_fit)
  
  if (length(slopes) < 3 || all(!is.finite(slopes))) {
    return(list(run = NULL, slopes = slopes, all_runs = NULL))
  }
  
  lo <- max(1L, floor(length(slopes) * SEARCH_START_FRAC))
  hi <- min(length(slopes), ceiling(length(slopes) * SEARCH_END_FRAC))
  
  neg_flag <- rep(FALSE, length(slopes))
  idx_search <- lo:hi
  neg_flag[idx_search] <- is.finite(slopes[idx_search]) & (slopes[idx_search] < -NEG_SLOPE_EPS)
  
  run_tbl <- rle_to_runs(neg_flag)
  run_tbl <- merge_descending_runs(run_tbl, max_gap = MAX_GAP_POINTS)
  
  if (nrow(run_tbl) == 0) {
    return(list(run = NULL, slopes = slopes, all_runs = NULL))
  }
  
  run_tbl <- run_tbl %>% filter(len >= MIN_RUN_POINTS)
  if (nrow(run_tbl) == 0) {
    return(list(run = NULL, slopes = slopes, all_runs = NULL))
  }
  
  scored <- score_descending_runs(time, o2_fit, slopes, run_tbl) %>%
    arrange(desc(score), desc(duration), desc(total_drop))
  
  list(
    run      = scored[1, ],
    slopes   = slopes,
    all_runs = scored
  )
}

find_main_curve_start <- function(time, o2_fit, idx_run_start,
                                  pretrans_slope_eps = PRETRANS_SLOPE_EPS,
                                  pretrans_min_run   = PRETRANS_MIN_RUN,
                                  backoff_points     = START_BACKOFF_POINTS,
                                  slope_trend_eps    = SLOPE_TREND_EPS,
                                  smooth_k           = START_SMOOTH_K) {
  slopes <- get_slopes(time, o2_fit)
  slopes_sm <- smooth_slopes(slopes, k = smooth_k)
  n <- length(o2_fit)
  
  idx_run_start <- max(3L, min(idx_run_start, n))
  idx_seed <- max(2L, min(idx_run_start - backoff_points, idx_run_start))
  
  left_bound <- max(2L, idx_seed - START_LOOKBACK_POINTS)
  candidate_points <- seq(idx_seed, left_bound, by = -1L)
  
  best_start <- idx_seed
  
  for (p in candidate_points) {
    j1 <- p
    j2 <- min(length(slopes_sm), p + START_STABLE_RUN_LEN - 1L)
    
    if ((j2 - j1 + 1L) < START_STABLE_RUN_LEN) next
    
    seg_s <- slopes_sm[j1:j2]
    y2 <- min(n, p + START_STABLE_RUN_LEN)
    seg_y <- o2_fit[p:y2]
    
    if (any(!is.finite(seg_s)) || any(!is.finite(seg_y))) next
    
    cond_neg <- mean(seg_s < 0, na.rm = TRUE) >= START_NEG_FRAC
    cond_mono <- all(diff(seg_y) <= START_WIGGLE_TOL)
    
    s_sign <- sign(seg_s)
    s_sign[abs(seg_s) <= START_WIGGLE_TOL] <- 0
    sign_changes <- sum(diff(s_sign) != 0, na.rm = TRUE)
    cond_wiggle <- sign_changes <= START_MAX_SIGN_CHANGES
    
    k <- START_LEVEL_RUN_LEN
    j3 <- min(length(slopes_sm), p + k - 1L)
    j4 <- min(length(slopes_sm), p + 2L * k - 1L)
    
    cond_level <- TRUE
    
    if ((j3 - p + 1L) >= k && (j4 - (p + k) + 1L) >= k) {
      seg_left  <- slopes_sm[p:j3]
      seg_right <- slopes_sm[(p + k):j4]
      
      if (all(is.finite(seg_left)) && all(is.finite(seg_right))) {
        med_left  <- median(seg_left, na.rm = TRUE)
        med_right <- median(seg_right, na.rm = TRUE)
        
        if (is.finite(med_left) && is.finite(med_right) && abs(med_right) > 0) {
          cond_level <- (abs(med_left) / abs(med_right)) <= START_SLOPE_RATIO_MAX
        }
      }
    }
    
    if (cond_neg && cond_mono && cond_wiggle && cond_level) {
      best_start <- p
    } else {
      break
    }
  }
  
  best_start <- max(1L, min(best_start, idx_run_start))
  as.integer(best_start)
}

find_end_after_steepest_drop <- function(time, o2_fit, idx_peak,
                                         min_gap = END_MIN_GAP,
                                         run_len = END_RECOV_RUN,
                                         smooth_k = END_SMOOTH_K,
                                         recovery_frac = END_RECOV_FRAC,
                                         mono_frac = END_MONO_FRAC,
                                         slope_eps = END_SLOPE_EPS) {
  slopes <- diff(o2_fit) / diff(time)
  nsl <- length(slopes)
  
  if (!is.finite(idx_peak) || is.na(idx_peak) || idx_peak >= nsl - run_len) {
    return(list(idx_end = NA_integer_, idx_steepest = NA_integer_, slopes_sm = slopes))
  }
  
  search_idx <- seq(idx_peak, nsl)
  search_idx <- search_idx[is.finite(slopes[search_idx])]
  
  if (!length(search_idx)) {
    return(list(idx_end = NA_integer_, idx_steepest = NA_integer_, slopes_sm = slopes))
  }
  
  slopes_sm <- slopes
  if (length(slopes) >= smooth_k) {
    sm <- zoo::rollmedian(slopes, k = smooth_k, fill = NA, align = "center")
    keep_na <- is.na(sm)
    sm[keep_na] <- slopes[keep_na]
    slopes_sm <- sm
  }
  
  valid_idx <- search_idx[is.finite(slopes_sm[search_idx])]
  if (!length(valid_idx)) {
    valid_idx <- search_idx
  }
  
  idx_steepest <- valid_idx[which.min(slopes_sm[valid_idx])]
  min_slope <- slopes_sm[idx_steepest]
  
  if (!is.finite(min_slope)) {
    return(list(idx_end = NA_integer_, idx_steepest = idx_steepest, slopes_sm = slopes_sm))
  }
  
  target <- min_slope * recovery_frac
  
  start_i <- idx_steepest + min_gap
  end_i   <- nsl - run_len + 1L
  
  if (start_i > end_i) {
    return(list(
      idx_end      = min(length(o2_fit), idx_steepest + 1L),
      idx_steepest = idx_steepest,
      slopes_sm    = slopes_sm
    ))
  }
  
  for (i in seq(start_i, end_i)) {
    seg <- slopes_sm[i:(i + run_len - 1L)]
    if (!all(is.finite(seg))) next
    
    cond_level <- mean(seg >= target, na.rm = TRUE) >= 0.8
    cond_trend <- mean(diff(seg) >= -slope_eps, na.rm = TRUE) >= mono_frac
    
    if (cond_level && cond_trend) {
      return(list(
        idx_end      = as.integer(i + 1L),
        idx_steepest = idx_steepest,
        slopes_sm    = slopes_sm
      ))
    }
  }
  
  tail_idx <- seq(idx_steepest + 1L, nsl)
  tail_idx <- tail_idx[is.finite(slopes_sm[tail_idx])]
  cand <- tail_idx[slopes_sm[tail_idx] >= target]
  
  if (length(cand)) {
    return(list(
      idx_end      = as.integer(cand[1] + 1L),
      idx_steepest = idx_steepest,
      slopes_sm    = slopes_sm
    ))
  }
  
  list(
    idx_end      = min(length(o2_fit), idx_steepest + 1L),
    idx_steepest = idx_steepest,
    slopes_sm    = slopes_sm
  )
}

trim_one_series <- function(df) {
  df <- df %>% arrange(Time)
  
  if (nrow(df) < MIN_ROWS_SERIES) {
    return(list(ok = FALSE, reason = "Too few rows", data = NULL, meta = NULL))
  }
  
  o2_fit <- safe_spline_fit(df$Time, df$Oxygen, spar = TRIM_SPAR)
  
  if (all(!is.finite(o2_fit))) {
    return(list(ok = FALSE, reason = "Spline failed", data = NULL, meta = NULL))
  }
  
  if ((max(o2_fit, na.rm = TRUE) - min(o2_fit, na.rm = TRUE)) < MIN_RANGE_O2) {
    return(list(ok = FALSE, reason = "Too flat", data = NULL, meta = NULL))
  }
  
  df <- df %>% mutate(O2_fit = o2_fit)
  
  run_out <- find_main_descending_run(df$Time, df$O2_fit)
  main_run <- run_out$run
  slopes   <- run_out$slopes
  all_runs <- run_out$all_runs
  
  if (is.null(main_run) || nrow(main_run) == 0) {
    return(list(ok = FALSE, reason = "No main descending run found", data = NULL, meta = NULL))
  }
  
  idx_run_start <- as.integer(main_run$start[1])
  idx_run_end   <- as.integer(main_run$end[1])
  
  idx_peak <- find_main_curve_start(
    time               = df$Time,
    o2_fit             = df$O2_fit,
    idx_run_start      = idx_run_start,
    pretrans_slope_eps = PRETRANS_SLOPE_EPS,
    pretrans_min_run   = PRETRANS_MIN_RUN,
    backoff_points     = START_BACKOFF_POINTS,
    slope_trend_eps    = SLOPE_TREND_EPS,
    smooth_k           = START_SMOOTH_K
  )
  
  if (!is.finite(idx_peak) || is.na(idx_peak)) {
    idx_peak <- idx_run_start
  }
  
  idx_peak_raw <- idx_peak
  
  end_out <- find_end_after_steepest_drop(
    time          = df$Time,
    o2_fit        = df$O2_fit,
    idx_peak      = idx_peak,
    min_gap       = END_MIN_GAP,
    run_len       = END_RECOV_RUN,
    smooth_k      = END_SMOOTH_K,
    recovery_frac = END_RECOV_FRAC,
    mono_frac     = END_MONO_FRAC,
    slope_eps     = END_SLOPE_EPS
  )
  
  idx_end      <- end_out$idx_end
  idx_steepest <- end_out$idx_steepest
  slopes_sm    <- end_out$slopes_sm
  
  if (!is.finite(idx_end) || is.na(idx_end)) {
    idx_end <- min(nrow(df), idx_run_end + 1L)
    end_reason <- "fallback_main_run_end"
  } else {
    end_reason <- "recovery_after_steepest_drop"
  }
  
  idx_peak <- max(1L, min(idx_peak, nrow(df)))
  idx_end  <- max(idx_peak + 1L, min(idx_end, nrow(df)))
  
  if (!is.finite(idx_steepest) || is.na(idx_steepest)) {
    idx_steepest <- max(idx_peak, min(nrow(df) - 1L, idx_run_start))
  }
  
  if (idx_end <= idx_peak) {
    return(list(ok = FALSE, reason = "Invalid peak/end order", data = NULL, meta = NULL))
  }
  
  trimmed <- df[idx_peak:idx_end, ] %>%
    mutate(
      curve_code           = df$curve_code[1],
      series_id            = df$series_id[1],
      peak_time_raw        = df$Time[idx_peak_raw],
      peak_time            = df$Time[idx_peak],
      main_run_start_time  = df$Time[idx_run_start],
      main_run_end_time    = df$Time[min(idx_run_end + 1L, nrow(df))],
      steepest_drop_time   = df$Time[min(idx_steepest + 1L, nrow(df))],
      end_time_chosen      = df$Time[idx_end],
      peak_idx_raw         = idx_peak_raw,
      peak_idx             = idx_peak,
      main_run_start_idx   = idx_run_start,
      main_run_end_idx     = idx_run_end,
      steepest_drop_idx    = idx_steepest,
      end_idx              = idx_end,
      start_shift_points   = idx_run_start - idx_peak,
      start_shift_time     = df$Time[idx_run_start] - df$Time[idx_peak],
      end_reason           = end_reason
    )
  
  meta <- tibble(
    curve_code = df$curve_code[1],
    series_id = df$series_id[1],
    T = df$T[1],
    Dose = df$Dose[1],
    Replicate = df$Replicate[1],
    n_points_total = nrow(df),
    n_points_trimmed = nrow(trimmed),
    peak_time_raw = df$Time[idx_peak_raw],
    peak_time = df$Time[idx_peak],
    main_run_start_time = df$Time[idx_run_start],
    main_run_end_time = df$Time[min(idx_run_end + 1L, nrow(df))],
    steepest_drop_time = df$Time[min(idx_steepest + 1L, nrow(df))],
    chosen_end_time = df$Time[idx_end],
    fit_duration_min = df$Time[idx_end] - df$Time[idx_peak],
    delta_Ninoc_to_N0_min = df$Time[idx_peak],
    start_shift_points = idx_run_start - idx_peak,
    start_shift_time = df$Time[idx_run_start] - df$Time[idx_peak],
    end_reason = end_reason,
    main_run_duration = main_run$duration[1],
    main_run_drop = main_run$total_drop[1],
    main_run_score = main_run$score[1]
  )
  
  list(
    ok            = TRUE,
    reason        = NA_character_,
    data          = trimmed,
    meta          = meta,
    idx_peak_raw  = idx_peak_raw,
    idx_peak      = idx_peak,
    idx_run_start = idx_run_start,
    idx_run_end   = idx_run_end,
    idx_steepest  = idx_steepest,
    idx_end       = idx_end,
    slopes        = slopes,
    slopes_sm     = slopes_sm,
    all_runs      = all_runs,
    full_df       = df,
    manual_override = FALSE,
    manual_keep_start = NA_real_,
    manual_keep_end   = NA_real_
  )
}

# ── Load data ────────────────────────────────────────────────────────────

raw <- readr::read_csv(file.path(tables_dir, "Oxygen_All_Long.csv"), show_col_types = FALSE) %>%
  mutate(
    T         = as.numeric(T),
    Dose      = as.character(Dose),
    Replicate = as.character(Replicate),
    series_id = paste0("T=", T, " | Dose=", Dose, " | Rep=", Replicate)
  )

required_cols <- c("Time", "T", "Dose", "Replicate", "Oxygen")
missing_cols <- setdiff(required_cols, names(raw))

if (length(missing_cols) > 0) {
  stop(
    "Missing required columns in Oxygen_All_Long.csv: ",
    paste(missing_cols, collapse = ", ")
  )
}

# ── Create short code for each curve ─────────────────────────────────────

curve_key <- raw %>%
  distinct(series_id, T, Dose, Replicate) %>%
  mutate(Dose_num = suppressWarnings(as.numeric(Dose))) %>%
  arrange(T, Dose_num, Dose, Replicate) %>%
  mutate(curve_code = paste0("C", str_pad(row_number(), width = 3, pad = "0"))) %>%
  select(curve_code, series_id, T, Dose, Replicate)

raw <- raw %>%
  left_join(curve_key, by = c("series_id", "T", "Dose", "Replicate"))

readr::write_csv(curve_key, file.path(tables_dir, "Oxygen_Curve_Code_Key.csv"))

# ── Validate manual overrides ────────────────────────────────────────────

if (nrow(manual_trim_overrides) > 0) {
  needed_override_cols <- c("curve_code", "keep_start", "keep_end")
  missing_override_cols <- setdiff(needed_override_cols, names(manual_trim_overrides))
  
  if (length(missing_override_cols) > 0) {
    stop(
      "manual_trim_overrides is missing columns: ",
      paste(missing_override_cols, collapse = ", ")
    )
  }
  
  bad_codes <- setdiff(manual_trim_overrides$curve_code, curve_key$curve_code)
  if (length(bad_codes) > 0) {
    stop(
      "These curve_code values were not found: ",
      paste(bad_codes, collapse = ", ")
    )
  }
  
  bad_ranges <- manual_trim_overrides %>%
    filter(
      !is.na(keep_start) & !is.finite(keep_start) |
        !is.na(keep_end)   & !is.finite(keep_end)
    )
  
  if (nrow(bad_ranges) > 0) {
    stop("Some manual overrides have invalid keep_start or keep_end values.")
  }
}

series_ids <- unique(raw$series_id)

trimmed_lst <- vector("list", length(series_ids))
names(trimmed_lst) <- series_ids

meta_lst <- vector("list", length(series_ids))
names(meta_lst) <- series_ids

diag_lst <- vector("list", length(series_ids))
names(diag_lst) <- series_ids

skipped_log <- tibble(
  curve_code = character(),
  series_id  = character(),
  reason     = character()
)

for (sid in series_ids) {
  df <- raw %>%
    filter(series_id == sid) %>%
    arrange(Time)
  
  this_code <- df$curve_code[1]
  
  out <- trim_one_series(df)
  
  if (!isTRUE(out$ok)) {
    skipped_log <- add_row(
      skipped_log,
      curve_code = this_code,
      series_id = sid,
      reason = out$reason
    )
    next
  }
  
  ov <- manual_trim_overrides %>%
    filter(curve_code == this_code)
  
  if (nrow(ov) > 0) {
    auto_start <- df$Time[out$idx_peak]
    auto_end   <- df$Time[out$idx_end]
    
    final_start <- if (!is.na(ov$keep_start[1])) ov$keep_start[1] else auto_start
    final_end   <- if (!is.na(ov$keep_end[1]))   ov$keep_end[1]   else auto_end
    
    if (!is.finite(final_start) || !is.finite(final_end) || final_end <= final_start) {
      skipped_log <- add_row(
        skipped_log,
        curve_code = this_code,
        series_id = sid,
        reason = "Manual override gave invalid final start/end"
      )
      next
    }
    
    manual_trimmed <- df %>%
      filter(Time >= final_start, Time <= final_end)
    
    if (nrow(manual_trimmed) < 2) {
      skipped_log <- add_row(
        skipped_log,
        curve_code = this_code,
        series_id = sid,
        reason = "Manual override produced fewer than 2 rows"
      )
      next
    }
    
    manual_trimmed <- manual_trimmed %>%
      mutate(
        O2_fit = safe_spline_fit(Time, Oxygen, spar = TRIM_SPAR),
        curve_code = this_code,
        series_id = df$series_id[1],
        peak_time_raw = auto_start,
        peak_time = final_start,
        main_run_start_time = df$Time[out$idx_run_start],
        main_run_end_time = df$Time[min(out$idx_run_end + 1L, nrow(df))],
        steepest_drop_time = if (is.finite(out$idx_steepest)) df$Time[min(out$idx_steepest + 1L, nrow(df))] else NA_real_,
        end_time_chosen = final_end,
        peak_idx_raw = out$idx_peak_raw,
        peak_idx = NA_integer_,
        main_run_start_idx = out$idx_run_start,
        main_run_end_idx = out$idx_run_end,
        steepest_drop_idx = out$idx_steepest,
        end_idx = NA_integer_,
        start_shift_points = NA_integer_,
        start_shift_time = final_start - auto_start,
        end_reason = "manual_override_partial_or_full"
      )
    
    manual_meta <- tibble(
      curve_code = this_code,
      series_id = df$series_id[1],
      T = df$T[1],
      Dose = df$Dose[1],
      Replicate = df$Replicate[1],
      n_points_total = nrow(df),
      n_points_trimmed = nrow(manual_trimmed),
      peak_time_raw = auto_start,
      peak_time = final_start,
      main_run_start_time = df$Time[out$idx_run_start],
      main_run_end_time = df$Time[min(out$idx_run_end + 1L, nrow(df))],
      steepest_drop_time = if (is.finite(out$idx_steepest)) df$Time[min(out$idx_steepest + 1L, nrow(df))] else NA_real_,
      chosen_end_time = final_end,
      fit_duration_min = final_end - final_start,
      delta_Ninoc_to_N0_min = final_start,
      start_shift_points = NA_integer_,
      start_shift_time = final_start - auto_start,
      end_reason = "manual_override_partial_or_full",
      main_run_duration = if (!is.null(out$all_runs) && nrow(out$all_runs) > 0) out$all_runs$duration[1] else NA_real_,
      main_run_drop = if (!is.null(out$all_runs) && nrow(out$all_runs) > 0) out$all_runs$total_drop[1] else NA_real_,
      main_run_score = if (!is.null(out$all_runs) && nrow(out$all_runs) > 0) out$all_runs$score[1] else NA_real_
    )
    
    out$data <- manual_trimmed
    out$meta <- manual_meta
    out$manual_override <- TRUE
    out$manual_keep_start <- if (!is.na(ov$keep_start[1])) ov$keep_start[1] else NA_real_
    out$manual_keep_end   <- if (!is.na(ov$keep_end[1]))   ov$keep_end[1]   else NA_real_
    out$final_start <- final_start
    out$final_end   <- final_end
  } else {
    out$final_start <- df$Time[out$idx_peak]
    out$final_end   <- df$Time[out$idx_end]
  }
  
  trimmed_lst[[sid]] <- out$data
  meta_lst[[sid]]    <- out$meta
  diag_lst[[sid]]    <- out
}

trimmed <- bind_rows(trimmed_lst)
trim_meta <- bind_rows(meta_lst)

readr::write_csv(trimmed, file.path(tables_dir, "Oxygen_Data_Smoothed_Trimmed.csv"))

filtered <- trimmed %>%
  select(curve_code, T, Dose, Replicate, Time, Oxygen, O2_fit)

readr::write_csv(filtered, file.path(tables_dir, "Oxygen_Data_Filtered.csv"))
readr::write_csv(skipped_log, file.path(tables_dir, "Skipped_Series_Log.csv"))

trim_meta <- trim_meta %>%
  mutate(
    Dose_num = suppressWarnings(as.numeric(Dose)),
    Dose_sort = if_else(is.na(Dose_num), -Inf, Dose_num)
  ) %>%
  arrange(T, Dose_sort, Replicate) %>%
  select(-Dose_num, -Dose_sort)

readr::write_csv(trim_meta, file.path(tables_dir, "Oxygen_Trimmed_Series_Metadata.csv"))

# ── Diagnostics PDF ──────────────────────────────────────────────────────

pdf(file.path(figures_dir, "oxygen_trimming_diagnostics.pdf"), width = 8.2, height = 6.2)

for (sid in names(diag_lst)) {
  out <- diag_lst[[sid]]
  if (is.null(out)) next
  
  df <- out$full_df
  this_code <- df$curve_code[1]
  
  idx_peak_raw  <- out$idx_peak_raw
  idx_peak      <- out$idx_peak
  idx_run_start <- out$idx_run_start
  idx_steepest  <- out$idx_steepest
  idx_end       <- out$idx_end
  
  xmin_show <- out$final_start
  xmax_show <- out$final_end
  
  rect_df <- tibble(
    xmin = xmin_show,
    xmax = xmax_show,
    ymin = -Inf,
    ymax = Inf
  )
  
  subtitle_text <- if (isTRUE(out$manual_override)) {
    paste0(
      "MANUAL OVERRIDE | final kept region = ",
      xmin_show, " to ", xmax_show, " min | ",
      "entered start = ",
      ifelse(is.na(out$manual_keep_start), "auto", as.character(out$manual_keep_start)),
      " | entered end = ",
      ifelse(is.na(out$manual_keep_end), "auto", as.character(out$manual_keep_end)),
      " | spar = ", TRIM_SPAR
    )
  } else {
    paste0(
      "Grey dotted = chosen main-curve start | ",
      "Black = chosen start | ",
      "Green dashed = main descending run start | ",
      "Purple dashed = steepest drop | ",
      "Magenta = chosen end | ",
      "spar = ", TRIM_SPAR
    )
  }
  
  p <- ggplot(df, aes(Time, Oxygen)) +
    geom_rect(
      data = rect_df,
      inherit.aes = FALSE,
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
      fill = "orange",
      alpha = 0.12
    ) +
    geom_line(colour = "grey60", linewidth = 0.8) +
    geom_line(aes(y = O2_fit), colour = "blue", linewidth = 1.1)
  
  if (!isTRUE(out$manual_override)) {
    p <- p +
      geom_vline(xintercept = df$Time[idx_peak_raw], colour = "grey20", linetype = "dotted", linewidth = 0.8) +
      geom_vline(xintercept = df$Time[idx_peak], colour = "black", linetype = "solid", linewidth = 0.9) +
      geom_vline(xintercept = df$Time[idx_run_start], colour = "darkgreen", linetype = "dashed", linewidth = 0.8) +
      geom_vline(xintercept = df$Time[min(idx_steepest + 1L, nrow(df))], colour = "purple", linetype = "dashed", linewidth = 0.8) +
      geom_vline(xintercept = df$Time[idx_end], colour = "magenta", linetype = "solid", linewidth = 0.9)
  } else {
    p <- p +
      geom_vline(xintercept = xmin_show, colour = "black", linetype = "solid", linewidth = 0.9) +
      geom_vline(xintercept = xmax_show, colour = "magenta", linetype = "solid", linewidth = 0.9)
  }
  
  p <- p +
    scale_x_continuous(
      breaks = pretty_breaks(n = X_MAJOR_BREAKS_N),
      minor_breaks = pretty_breaks(n = X_MAJOR_BREAKS_N * 2)
    ) +
    scale_y_continuous(
      breaks = pretty_breaks(n = Y_MAJOR_BREAKS_N),
      minor_breaks = pretty_breaks(n = Y_MAJOR_BREAKS_N * 2)
    ) +
    labs(
      title = paste0(this_code, " | ", sid),
      subtitle = subtitle_text,
      x = "Time (min)",
      y = "O2 (mg L^-1)"
    ) +
    theme_classic(base_size = 11) +
    theme(
      axis.text.x = element_text(size = 9),
      axis.text.y = element_text(size = 9),
      axis.ticks.length = unit(0.18, "cm"),
      plot.title = element_text(face = "bold"),
      panel.grid.major = element_line(colour = "grey88", linewidth = 0.25),
      panel.grid.minor = element_line(colour = "grey94", linewidth = 0.2)
    )
  
  print(p)
}

dev.off()

message("Saved: ", file.path(tables_dir, "Oxygen_Curve_Code_Key.csv"))
message("Saved: ", file.path(tables_dir, "Oxygen_Data_Smoothed_Trimmed.csv"))
message("Saved: ", file.path(tables_dir, "Oxygen_Data_Filtered.csv"))
message("Saved: ", file.path(tables_dir, "Oxygen_Trimmed_Series_Metadata.csv"))
message("Saved: ", file.path(tables_dir, "Skipped_Series_Log.csv"))
message("Saved: ", file.path(figures_dir, "oxygen_trimming_diagnostics.pdf"))