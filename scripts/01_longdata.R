# =============================================================================
# 01_longdata.R — Convert raw oxygen CSVs from wide to long format
# =============================================================================
# Reads:  data/*.csv  (raw oxygen time-series at each temperature)
# Writes: tables/Oxygen_All_Long.csv
# =============================================================================

# Source shared config (works interactively in RStudio, via source(), or from CLI)
.this_dir <- if (requireNamespace("rstudioapi", quietly = TRUE) &&
                 rstudioapi::isAvailable() &&
                 nzchar(rstudioapi::getActiveDocumentContext()$path)) {
  dirname(rstudioapi::getActiveDocumentContext()$path)
} else {
  tryCatch(dirname(sys.frame(1)$ofile), error = function(e) getwd())
}
source(file.path(.this_dir, "00_config.R"))

# 0) Check where R is looking
message("Working directory: ", getwd())
message("Reading data from: ", data_dir)

# 1) Find files recursively inside data/
all_files <- list.files(
  path = data_dir,
  recursive = TRUE,
  full.names = TRUE
)

# Keep files like:
# 15zymo_Oxygen.csv, 18zymo_Oxygen.CSV, "15 zymo_Oxygen.csv", etc.
files <- all_files[
  stringr::str_detect(
    basename(all_files),
    regex("^\\s*\\d+\\s*zymo_Oxygen\\s*\\.csv\\s*$", ignore_case = TRUE)
  )
]

if (length(files) == 0) {
  stop(
    paste(
      "No matching zymo oxygen CSVs found in", data_dir,
      "\nCheck file names with list.files(data_dir, recursive = TRUE)."
    )
  )
}

message("Found files:")
print(files)

# 2) Read one file, clean names, add File column, derive T from filename
read_one <- function(path) {
  df <- readr::read_csv(path, show_col_types = FALSE)

  # Trim stray spaces in column names
  names(df) <- stringr::str_trim(names(df))

  # Extract leading temperature number from filename
  temp_from_name <- stringr::str_extract(basename(path), "^\\s*\\d+")
  temp_from_name <- suppressWarnings(as.numeric(stringr::str_trim(temp_from_name)))

  # Add T if missing or all NA
  if (!"T" %in% names(df) || all(is.na(df$T))) {
    df <- df %>% mutate(T = temp_from_name)
  }

  # Add source filename
  df %>% mutate(File = basename(path), .before = 1)
}

# 3) Read and bind all files
raw <- purrr::map_dfr(files, read_one)

# 4) Pivot dose/replicate columns from wide -> long
long_data <- raw %>%
  pivot_longer(
    cols = tidyselect::matches("^(?i)(Control|\\d*\\.?\\d+)_R\\d+$"),
    names_to = c("Dose", "Replicate"),
    names_sep = "_",
    values_to = "Oxygen"
  ) %>%
  mutate(
    Dose = if_else(stringr::str_to_lower(Dose) == "control", "Control", Dose),
    Replicate = toupper(Replicate)
  ) %>%
  tidyr::drop_na(Oxygen)

# 5) Order doses (Control first, then numeric ascending)
dose_levels <- long_data %>%
  distinct(Dose) %>%
  mutate(
    Dose_num = suppressWarnings(as.numeric(Dose)),
    Dose_sort = if_else(is.na(Dose_num), -Inf, Dose_num)
  ) %>%
  arrange(Dose_sort) %>%
  pull(Dose)

# 6) Final tidy output
long_data <- long_data %>%
  mutate(Dose_factor = factor(Dose, levels = dose_levels)) %>%
  arrange(File, T, Dose_factor, Replicate, Time) %>%
  select(any_of(c("File", "Time", "T", "Dose", "Replicate", "Oxygen")))

# 7) Save table
out_table <- file.path(tables_dir, "Oxygen_All_Long.csv")
readr::write_csv(long_data, out_table)

# 8) Checks
message("\nOutput written: ", out_table)
glimpse(long_data)
print(long_data %>% count(File, T, Dose, Replicate), n = Inf)
