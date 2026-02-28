# =============================================================================
# NSCA Poster Analysis
# Neuromuscular Performance Assessment: Basketball vs. Rugby
# Dynamic Strength Index & Temporal Monitoring
#
# Authors: Pablo Tadeo Ríos-Gallardo, PhD & Samuel Montalvo, PhD
# =============================================================================

# =============================================================================
# 0. PACKAGES & PATHS
# =============================================================================
suppressPackageStartupMessages({
  library(tidyverse)
  library(lubridate)
  library(patchwork)
  library(scales)
  library(gt)
  library(gtExtras)
  library(viridis)
})

BASE_DIR <- "/Volumes/ADATA HD710 PRO/fod_rugby"
RAW_DIR  <- file.path(BASE_DIR, "Raw_data")
FIG_DIR  <- file.path(BASE_DIR, "Figures")
TBL_DIR  <- file.path(BASE_DIR, "Tables")

dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(TBL_DIR, showWarnings = FALSE, recursive = TRUE)

cat("=== NSCA Poster Analysis — Basketball & Rugby ===\n\n")

# =============================================================================
# 1. CLEAN OUTPUT DIRECTORIES
# =============================================================================
cat("--- Cleaning output directories ---\n")
invisible(file.remove(list.files(FIG_DIR, full.names = TRUE)))
invisible(file.remove(list.files(TBL_DIR, full.names = TRUE)))
cat("Figures/ and Tables/ cleared.\n\n")

# =============================================================================
# 2. LOAD ALL 5 RAW CSV FILES (UTF-8)
# =============================================================================
cat("--- Loading raw data (UTF-8) ---\n")

read_utf8 <- function(filename) {
  read_csv(
    file.path(RAW_DIR, filename),
    locale        = locale(encoding = "UTF-8"),
    show_col_types = FALSE,
    name_repair   = "minimal"
  )
}

bas_cmj_raw  <- read_utf8("CMJ_Basket - CMJ.csv")
bas_imtp_raw <- read_utf8("IMPT_Basket - Hoja 1.csv")
rug_demo_raw <- read_utf8("CMJ_IMTP_RUGBY - Hoja 1.csv")
rug_cmj_raw  <- read_utf8("CMJ_Rugby - Hoja 1.csv")
rug_imtp_raw <- read_utf8("IMPT_Rugby - Hoja 1.csv")

cat(sprintf("Basketball CMJ:  %d rows\n", nrow(bas_cmj_raw)))
cat(sprintf("Basketball IMTP: %d rows\n", nrow(bas_imtp_raw)))
cat(sprintf("Rugby Demo:      %d rows\n", nrow(rug_demo_raw)))
cat(sprintf("Rugby CMJ:       %d rows\n", nrow(rug_cmj_raw)))
cat(sprintf("Rugby IMTP:      %d rows\n\n", nrow(rug_imtp_raw)))

# =============================================================================
# 3. STANDARDIZE & CLEAN EACH DATASET
# =============================================================================
cat("--- Standardizing datasets ---\n")

# ---- Helper: coerce numeric safely ----
num <- function(x) suppressWarnings(as.numeric(x))

# ---- 3a. Basketball CMJ ----
# Note: Basketball CMJ has 'Week' (capital W), column 4
bas_cmj <- bas_cmj_raw |>
  rename(
    athlete_name     = Name,
    athlete_id       = ExternalId,
    week             = Week,            # capital W in basketball
    date_chr         = Date,
    bw_kg            = `BW [KG]`,
    jump_height_cm   = `Jump Height (Imp-Mom) [cm]`,
    peak_power_w     = `Peak Power [W]`,
    rsi_mod          = `RSI-modified [m/s]`,
    cmj_peak_force_n = `Concentric Peak Force [N]`
  ) |>
  select(athlete_id, athlete_name, week, date_chr, bw_kg,
         jump_height_cm, peak_power_w, rsi_mod, cmj_peak_force_n) |>
  mutate(
    sport  = "Basketball",
    week   = as.integer(num(week)),
    date   = suppressWarnings(dmy(date_chr)),
    across(c(bw_kg, jump_height_cm, peak_power_w, rsi_mod, cmj_peak_force_n), num),
    jump_height_cm   = ifelse(jump_height_cm   < 0 | jump_height_cm   > 120, NA, jump_height_cm),
    peak_power_w     = ifelse(peak_power_w     < 0, NA, peak_power_w),
    rsi_mod          = ifelse(rsi_mod          < 0 | rsi_mod > 5,   NA, rsi_mod),
    cmj_peak_force_n = ifelse(cmj_peak_force_n <= 0, NA, cmj_peak_force_n)
  ) |>
  filter(!is.na(athlete_id), !is.na(week))

# ---- 3b. Basketball IMTP ----
# Note: Basketball IMTP has 'Week' (capital W)
bas_imtp <- bas_imtp_raw |>
  rename(
    athlete_id        = ExternalId,
    week              = Week,           # capital W in basketball
    imtp_peak_force_n = `Peak Vertical Force [N]`
  ) |>
  select(athlete_id, week, imtp_peak_force_n) |>
  mutate(
    week              = as.integer(num(week)),
    imtp_peak_force_n = num(imtp_peak_force_n),
    imtp_peak_force_n = ifelse(imtp_peak_force_n <= 0, NA, imtp_peak_force_n)
  ) |>
  filter(!is.na(athlete_id), !is.na(week))

# ---- 3c. Rugby Demographics ----
# Note: 'externalId' (lowercase e), 'posicion' (no trailing space after read_csv)
rug_demo <- rug_demo_raw |>
  rename_with(~ trimws(.x)) |>          # defensive: strip any trailing spaces
  rename(
    athlete_id  = externalId,
    weight_kg   = `peso (kg)`,
    height_cm   = `altura(cm)`,
    birth_date  = fecha_nacimiento
  ) |>
  mutate(
    weight_kg  = num(weight_kg),
    height_cm  = num(height_cm),
    birth_date = suppressWarnings(dmy(birth_date)),
    age_years  = as.numeric(difftime(as.Date("2025-06-01"),
                                     birth_date, units = "days")) / 365.25
  ) |>
  select(athlete_id, weight_kg, height_cm, age_years)

cat(sprintf("Rugby demographics: %d athletes, weight [%.1f–%.1f kg], height [%.0f–%.0f cm]\n",
            nrow(rug_demo),
            min(rug_demo$weight_kg, na.rm = TRUE),
            max(rug_demo$weight_kg, na.rm = TRUE),
            min(rug_demo$height_cm, na.rm = TRUE),
            max(rug_demo$height_cm, na.rm = TRUE)))

# ---- 3d. Rugby CMJ ----
# Note: Rugby CMJ has 'week' (lowercase w), column 7; fewer columns than basketball
rug_cmj <- rug_cmj_raw |>
  rename(
    athlete_name     = Name,
    athlete_id       = ExternalId,
    # week is already lowercase 'week' in this file
    date_chr         = Date,
    bw_kg            = `BW [KG]`,
    jump_height_cm   = `Jump Height (Imp-Mom) [cm]`,
    peak_power_w     = `Peak Power [W]`,
    rsi_mod          = `RSI-modified [m/s]`,
    cmj_peak_force_n = `Concentric Peak Force [N]`
  ) |>
  select(athlete_id, athlete_name, week, date_chr, bw_kg,
         jump_height_cm, peak_power_w, rsi_mod, cmj_peak_force_n) |>
  mutate(
    sport  = "Rugby",
    week   = as.integer(num(week)),
    date   = suppressWarnings(dmy(date_chr)),
    across(c(bw_kg, jump_height_cm, peak_power_w, rsi_mod, cmj_peak_force_n), num),
    jump_height_cm   = ifelse(jump_height_cm   < 0 | jump_height_cm   > 120, NA, jump_height_cm),
    peak_power_w     = ifelse(peak_power_w     < 0, NA, peak_power_w),
    rsi_mod          = ifelse(rsi_mod          < 0 | rsi_mod > 5,   NA, rsi_mod),
    cmj_peak_force_n = ifelse(cmj_peak_force_n <= 0, NA, cmj_peak_force_n)
  ) |>
  filter(!is.na(athlete_id), !is.na(week))

# ---- 3e. Rugby IMTP ----
# Note: Rugby IMTP has 'week' (lowercase w), column 7
rug_imtp <- rug_imtp_raw |>
  rename(
    athlete_id        = ExternalId,
    # week is already lowercase 'week'
    imtp_peak_force_n = `Peak Vertical Force [N]`
  ) |>
  select(athlete_id, week, imtp_peak_force_n) |>
  mutate(
    week              = as.integer(num(week)),
    imtp_peak_force_n = num(imtp_peak_force_n),
    imtp_peak_force_n = ifelse(imtp_peak_force_n <= 0, NA, imtp_peak_force_n)
  ) |>
  filter(!is.na(athlete_id), !is.na(week))

cat(sprintf(
  "Basketball CMJ clean: %d obs | %d athletes | weeks %d-%d\n",
  nrow(bas_cmj), n_distinct(bas_cmj$athlete_id),
  min(bas_cmj$week), max(bas_cmj$week)))
cat(sprintf(
  "Basketball IMTP clean: %d obs | %d athletes | weeks %d-%d\n",
  nrow(bas_imtp), n_distinct(bas_imtp$athlete_id),
  min(bas_imtp$week), max(bas_imtp$week)))
cat(sprintf(
  "Rugby CMJ clean: %d obs | %d athletes | weeks %d-%d\n",
  nrow(rug_cmj), n_distinct(rug_cmj$athlete_id),
  min(rug_cmj$week), max(rug_cmj$week)))
cat(sprintf(
  "Rugby IMTP clean: %d obs | %d athletes | weeks %d-%d\n\n",
  nrow(rug_imtp), n_distinct(rug_imtp$athlete_id),
  min(rug_imtp$week), max(rug_imtp$week)))

# =============================================================================
# 4. AGGREGATE TO ATHLETE × WEEK (best/mean trial per session)
# =============================================================================
cat("--- Aggregating to athlete × week ---\n")

agg_cmj <- function(df) {
  df |>
    group_by(athlete_id, athlete_name, sport, week) |>
    summarise(
      bw_kg            = mean(bw_kg,            na.rm = TRUE),
      jump_height_cm   = mean(jump_height_cm,   na.rm = TRUE),
      peak_power_w     = mean(peak_power_w,     na.rm = TRUE),
      rsi_mod          = mean(rsi_mod,          na.rm = TRUE),
      cmj_peak_force_n = mean(cmj_peak_force_n, na.rm = TRUE),
      .groups = "drop"
    ) |>
    mutate(across(where(is.numeric), ~ ifelse(is.nan(.x), NA_real_, .x)))
}

agg_imtp <- function(df) {
  df |>
    group_by(athlete_id, week) |>
    summarise(
      # IMTP: best trial (maximum peak force per session)
      imtp_peak_force_n = max(imtp_peak_force_n, na.rm = TRUE),
      .groups = "drop"
    ) |>
    mutate(imtp_peak_force_n = ifelse(is.infinite(imtp_peak_force_n),
                                      NA_real_, imtp_peak_force_n))
}

bas_cmj_agg  <- agg_cmj(bas_cmj)
rug_cmj_agg  <- agg_cmj(rug_cmj)
bas_imtp_agg <- agg_imtp(bas_imtp)
rug_imtp_agg <- agg_imtp(rug_imtp)

# =============================================================================
# 5. MERGE WITHIN SPORTS & CALCULATE DSI
# =============================================================================
cat("--- Merging CMJ + IMTP and calculating DSI ---\n")

merge_and_dsi <- function(cmj_df, imtp_df) {
  left_join(cmj_df, imtp_df, by = c("athlete_id", "week")) |>
    mutate(
      dsi = cmj_peak_force_n / imtp_peak_force_n,
      # Filter physiologically implausible DSI values
      dsi = ifelse(dsi < 0.3 | dsi > 2.0, NA_real_, dsi)
    )
}

bas_perf <- merge_and_dsi(bas_cmj_agg, bas_imtp_agg)
rug_perf <- merge_and_dsi(rug_cmj_agg, rug_imtp_agg)

# Attach Rugby demographics (weight_kg, height_cm, age_years from demo file)
rug_perf <- rug_perf |>
  left_join(rug_demo, by = "athlete_id") |>
  # Use demographic weight if bw_kg from CMJ differs (demo file is reference)
  mutate(bw_kg = coalesce(weight_kg, bw_kg)) |>
  select(-weight_kg)    # keep height_cm and age_years from demo join

cat(sprintf("Basketball: %d athlete-week obs | DSI available: %d\n",
            nrow(bas_perf), sum(!is.na(bas_perf$dsi))))
cat(sprintf("Rugby:      %d athlete-week obs | DSI available: %d\n\n",
            nrow(rug_perf), sum(!is.na(rug_perf$dsi))))

# =============================================================================
# 6. MASTER DATASET (both sports combined)
# =============================================================================
cat("--- Building master dataset ---\n")

# Basketball doesn't have height/age columns — add as NA for consistent binding
bas_perf <- bas_perf |>
  mutate(height_cm = NA_real_, age_years = NA_real_)

master <- bind_rows(bas_perf, rug_perf) |>
  mutate(sport = factor(sport, levels = c("Basketball", "Rugby")))

cat(sprintf("Master dataset: %d rows | %d athletes (%d Basketball, %d Rugby)\n",
            nrow(master),
            n_distinct(master$athlete_id),
            n_distinct(master$athlete_id[master$sport == "Basketball"]),
            n_distinct(master$athlete_id[master$sport == "Rugby"])))

cat(sprintf("DSI: Basketball mean=%.3f SD=%.3f | Rugby mean=%.3f SD=%.3f\n\n",
            mean(master$dsi[master$sport == "Basketball"], na.rm = TRUE),
            sd(master$dsi[master$sport == "Basketball"],   na.rm = TRUE),
            mean(master$dsi[master$sport == "Rugby"],     na.rm = TRUE),
            sd(master$dsi[master$sport == "Rugby"],       na.rm = TRUE)))

# =============================================================================
# 7. TEMPORAL SUMMARIES BY SPORT
# =============================================================================
cat("--- Computing temporal summaries ---\n")

weekly <- master |>
  group_by(sport, week) |>
  summarise(
    n = n(),
    across(c(jump_height_cm, peak_power_w, rsi_mod,
             cmj_peak_force_n, imtp_peak_force_n, dsi),
           list(mean = ~ mean(.x, na.rm = TRUE),
                sd   = ~   sd(.x, na.rm = TRUE)),
           .names = "{.col}__{.fn}"),
    .groups = "drop"
  ) |>
  mutate(across(where(is.numeric), ~ ifelse(is.nan(.x), NA_real_, .x)))

cat(sprintf("Weekly summaries: Basketball %d weeks | Rugby %d weeks\n\n",
            n_distinct(weekly$week[weekly$sport == "Basketball"]),
            n_distinct(weekly$week[weekly$sport == "Rugby"])))

# =============================================================================
# 8. DESIGN SYSTEM (NSCA poster style)
# =============================================================================

# Okabe-Ito colorblind-safe palette
SPORT_COLS  <- c(Basketball = "#0072B2", Rugby = "#E69F00")
SPORT_FILLS <- c(Basketball = "#0072B2", Rugby = "#E69F00")
COL_REF     <- "#D55E00"   # reference lines

# NSCA poster theme (theme_classic base)
theme_nsca <- function(base_size = 12) {
  theme_classic(base_size = base_size) %+replace%
    theme(
      plot.title       = element_text(size = base_size + 2, face = "bold",
                                      margin = margin(b = 4)),
      plot.subtitle    = element_text(size = base_size - 1, color = "grey35",
                                      margin = margin(b = 6)),
      plot.caption     = element_text(size = base_size - 3, color = "grey45",
                                      hjust = 0, margin = margin(t = 8)),
      axis.title       = element_text(size = base_size, face = "bold"),
      axis.text        = element_text(size = base_size - 1),
      legend.position  = "bottom",
      legend.title     = element_text(face = "bold", size = base_size - 1),
      legend.text      = element_text(size = base_size - 1),
      strip.text       = element_text(face = "bold", size = base_size - 1),
      strip.background = element_rect(fill = "grey92", color = NA),
      panel.grid.major.y = element_line(color = "grey88", linewidth = 0.35),
      plot.margin      = margin(10, 14, 10, 14)
    )
}

save_fig <- function(p, filename, width = 14, height = 9, dpi = 300) {
  path <- file.path(FIG_DIR, filename)
  ggsave(path, plot = p, width = width, height = height,
         dpi = dpi, bg = "white", units = "in")
  cat(sprintf("  Saved: %s (%.0f × %.0f in, %d dpi)\n",
              basename(path), width, height, dpi))
}

# =============================================================================
# 9. FIGURE 1 — NEUROMUSCULAR COMPARISON: BASKETBALL vs RUGBY (6-panel)
#    Box + Jitter for all key variables
# =============================================================================
cat("--- Figure 1: Neuromuscular comparison box+jitter ---\n")

var_labels <- c(
  jump_height_cm   = "Jump Height (cm)",
  peak_power_w     = "Peak Power (W)",
  rsi_mod          = "RSI-modified (m/s)",
  cmj_peak_force_n = "CMJ Peak Force (N)",
  imtp_peak_force_n = "IMTP Peak Force (N)",
  dsi              = "Dynamic Strength Index"
)

fig1_long <- master |>
  select(sport, all_of(names(var_labels))) |>
  pivot_longer(-sport, names_to = "variable", values_to = "value") |>
  filter(!is.na(value)) |>
  mutate(
    variable = factor(variable, levels = names(var_labels),
                      labels = unname(var_labels))
  )

fig1 <- ggplot(fig1_long, aes(x = sport, y = value,
                              fill = sport, color = sport)) +
  geom_boxplot(
    alpha = 0.50, outlier.shape = NA, width = 0.48,
    color = "grey25", linewidth = 0.55
  ) +
  geom_jitter(
    width = 0.13, alpha = 0.45, size = 1.6, stroke = 0
  ) +
  facet_wrap(~ variable, scales = "free_y", nrow = 2) +
  scale_fill_manual(values  = SPORT_FILLS, name = "Sport") +
  scale_color_manual(values = SPORT_COLS,  name = "Sport") +
  scale_y_continuous(expand = expansion(mult = 0.07)) +
  labs(
    title    = "Neuromuscular Performance Comparison: Basketball vs. Rugby",
    subtitle = "Female university athletes | Each point = one athlete-session observation",
    x        = NULL,
    y        = "Value",
    caption  = paste0(
      "CMJ = Countermovement Jump  |  IMTP = Isometric Mid-Thigh Pull  |  ",
      "RSI-mod = RSI-modified  |  DSI = Dynamic Strength Index (CMJ Peak Force / IMTP Peak Force)\n",
      "Box = IQR; centre line = median; whiskers = 1.5×IQR; points = individual observations."
    )
  ) +
  theme_nsca() +
  theme(
    axis.text.x = element_text(size = 11, face = "bold"),
    legend.position = "none"          # sport shown on x-axis
  )

save_fig(fig1, "Fig1_Neuromuscular_Comparison.png", width = 16, height = 10)

# =============================================================================
# 10. FIGURE 2 — DSI DISTRIBUTION & FORCE-VELOCITY PROFILE
#     Box + Jitter with reference line and zone annotations
# =============================================================================
cat("--- Figure 2: DSI distribution by sport ---\n")

dsi_data <- master |> filter(!is.na(dsi))

# Compute % with DSI < 1.0 per sport (Force-Velocity Deficit)
dsi_pct <- dsi_data |>
  group_by(sport) |>
  summarise(
    pct_deficit = round(mean(dsi < 1.0) * 100, 1),
    pct_adequate = round(mean(dsi >= 1.0) * 100, 1),
    y_lbl = quantile(dsi, 0.98, na.rm = TRUE) + 0.02,
    .groups = "drop"
  ) |>
  mutate(label = paste0("FV Deficit:\n", pct_deficit, "%"))

fig2 <- ggplot(dsi_data, aes(x = sport, y = dsi,
                              fill = sport, color = sport)) +
  # shading below/above 1.0 for reference
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 1.0,
           fill = "#FFF4E0", alpha = 0.6) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 1.0, ymax = Inf,
           fill = "#E8F4F8", alpha = 0.6) +
  geom_boxplot(
    alpha = 0.60, outlier.shape = NA, width = 0.42,
    color = "grey20", linewidth = 0.65
  ) +
  geom_jitter(width = 0.11, alpha = 0.55, size = 2.5, stroke = 0) +
  geom_hline(yintercept = 1.0, linetype = "dashed",
             color = COL_REF, linewidth = 1.0) +
  geom_text(data = dsi_pct,
            aes(y = y_lbl, label = label),
            inherit.aes = FALSE,
            x = c(1, 2),
            size = 4.2, fontface = "bold",
            color = c(SPORT_COLS["Basketball"], SPORT_COLS["Rugby"])) +
  annotate("text", x = 2.48, y = 1.015, label = "DSI = 1.0",
           hjust = 0, color = COL_REF, size = 3.8, fontface = "italic") +
  annotate("text", x = 0.6, y = 0.85,
           label = "Force-Velocity\nDeficit Zone", size = 3.5,
           hjust = 0, color = "#C07800", fontface = "italic") +
  annotate("text", x = 0.6, y = 1.15,
           label = "Adequate\nStrength Zone", size = 3.5,
           hjust = 0, color = "#005A8E", fontface = "italic") +
  scale_fill_manual(values = SPORT_FILLS, guide = "none") +
  scale_color_manual(values = SPORT_COLS, guide = "none") +
  scale_y_continuous(breaks = seq(0.4, 1.6, 0.1),
                     expand = expansion(mult = 0.08)) +
  labs(
    title    = "Dynamic Strength Index: Basketball vs. Rugby",
    subtitle = "DSI = CMJ Peak Force / IMTP Peak Force  |  Dashed line = reference value 1.0",
    x        = "Sport",
    y        = "Dynamic Strength Index (DSI)",
    caption  = paste0(
      "DSI < 1.0 indicates Force-Velocity Deficit (neuromuscular strength prioritized in training).\n",
      "DSI ≥ 1.0 indicates adequate isometric force capacity. ",
      "Box = IQR; centre line = median; points = individual athlete-sessions."
    )
  ) +
  theme_nsca(base_size = 13) +
  theme(axis.text.x = element_text(size = 13, face = "bold"))

save_fig(fig2, "Fig2_DSI_Sport_Comparison.png", width = 9, height = 8)

# =============================================================================
# 11. FIGURE 3 — TEMPORAL CHANGES ACROSS THE SEASON BY SPORT (6-panel)
# =============================================================================
cat("--- Figure 3: Temporal changes by sport ---\n")

temporal_vars <- c(
  jump_height_cm__mean   = "Jump Height (cm)",
  peak_power_w__mean     = "Peak Power (W)",
  rsi_mod__mean          = "RSI-modified (m/s)",
  cmj_peak_force_n__mean = "CMJ Peak Force (N)",
  imtp_peak_force_n__mean = "IMTP Peak Force (N)",
  dsi__mean              = "Dynamic Strength Index"
)
temporal_sd <- sub("__mean$", "__sd", names(temporal_vars))

temporal_long <- weekly |>
  select(sport, week, n,
         all_of(names(temporal_vars)),
         all_of(temporal_sd)) |>
  pivot_longer(
    cols      = all_of(names(temporal_vars)),
    names_to  = "var_mean",
    values_to = "mean_val"
  ) |>
  mutate(
    sd_col  = sub("__mean$", "__sd", var_mean),
    sd_val  = map2_dbl(row_number(), sd_col, function(i, col) {
      # extract corresponding SD from the original wide data
      weekly[[col]][
        weekly$sport == sport[i] & weekly$week == week[i]
      ][1]
    }),
    variable = factor(var_mean, levels = names(temporal_vars),
                      labels = unname(temporal_vars))
  ) |>
  filter(!is.na(mean_val))

fig3 <- ggplot(temporal_long,
               aes(x = week, y = mean_val,
                   color = sport, fill = sport)) +
  geom_ribbon(
    aes(ymin = mean_val - sd_val,
        ymax = mean_val + sd_val),
    alpha = 0.15, color = NA, na.rm = TRUE
  ) +
  geom_line(linewidth = 1.0, na.rm = TRUE) +
  geom_point(size = 2.0, na.rm = TRUE, alpha = 0.85) +
  facet_wrap(~ variable, scales = "free_y", nrow = 2) +
  scale_color_manual(values = SPORT_COLS, name = "Sport") +
  scale_fill_manual(values  = SPORT_FILLS, name = "Sport") +
  scale_x_continuous(breaks = seq(0, 35, 5),
                     expand = expansion(add = 0.5)) +
  labs(
    title    = "Temporal Changes in Neuromuscular Performance by Sport",
    subtitle = "Weekly group mean ± 1 SD  |  Basketball: weeks 1–33  |  Rugby: weeks 1–31",
    x        = "Monitoring Week",
    y        = "Mean Value",
    caption  = paste0(
      "CMJ = Countermovement Jump  |  IMTP = Isometric Mid-Thigh Pull  |  ",
      "RSI-mod = RSI-modified  |  DSI = Dynamic Strength Index\n",
      "Ribbon = ±1 SD. Note: week numbering represents monitoring weeks within each sport's season."
    )
  ) +
  theme_nsca() +
  theme(legend.position = "bottom",
        legend.key.width = unit(1.2, "cm"))

save_fig(fig3, "Fig3_Temporal_Changes.png", width = 17, height = 10)

# =============================================================================
# 12. FIGURE 4 — CMJ PEAK FORCE vs IMTP PEAK FORCE BY SPORT
# =============================================================================
cat("--- Figure 4: CMJ vs IMTP Peak Force scatter ---\n")

force_data <- master |>
  filter(!is.na(cmj_peak_force_n), !is.na(imtp_peak_force_n))

# Per-sport correlation
r_vals <- force_data |>
  group_by(sport) |>
  summarise(
    r   = cor(cmj_peak_force_n, imtp_peak_force_n, use = "complete.obs"),
    x   = quantile(imtp_peak_force_n, 0.05, na.rm = TRUE),
    y   = max(cmj_peak_force_n, na.rm = TRUE) * 0.97,
    lbl = paste0("r = ", round(r, 2)),
    .groups = "drop"
  )

fig4 <- ggplot(force_data,
               aes(x = imtp_peak_force_n, y = cmj_peak_force_n,
                   color = sport)) +
  geom_point(alpha = 0.45, size = 2.2, stroke = 0) +
  geom_smooth(method = "lm", aes(fill = sport),
              se = TRUE, alpha = 0.15, linewidth = 1.2) +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted",
              color = "grey55", linewidth = 0.9) +
  geom_text(data = r_vals, aes(x = x, y = y, label = lbl, color = sport),
            inherit.aes = FALSE, size = 4.5, fontface = "bold", hjust = 0) +
  scale_color_manual(values = SPORT_COLS,  name = "Sport") +
  scale_fill_manual(values  = SPORT_FILLS, name = "Sport") +
  scale_x_continuous(labels = comma, expand = expansion(mult = 0.05)) +
  scale_y_continuous(labels = comma, expand = expansion(mult = 0.05)) +
  labs(
    title    = "CMJ Peak Force vs. IMTP Peak Force by Sport",
    subtitle = "Each point = one athlete-session  |  Dotted line = line of identity (slope = 1)",
    x        = "IMTP Peak Force (N)",
    y        = "CMJ Peak Force (N)",
    caption  = "Solid lines = OLS regression per sport (shaded = 95% CI). r = Pearson correlation."
  ) +
  theme_nsca(base_size = 13) +
  theme(legend.position = "right")

save_fig(fig4, "Fig4_CMJ_vs_IMTP_Force.png", width = 11, height = 8)

# =============================================================================
# 13. FIGURE 5 — INDIVIDUAL DSI TRAJECTORIES BY SPORT
# =============================================================================
cat("--- Figure 5: Individual DSI trajectories ---\n")

# Athletes with >= 4 DSI observations
ath_traj <- master |>
  filter(!is.na(dsi)) |>
  count(sport, athlete_id) |>
  filter(n >= 4) |>
  select(athlete_id)

fig5_data <- master |>
  filter(!is.na(dsi)) |>
  semi_join(ath_traj, by = "athlete_id")

fig5 <- ggplot(fig5_data,
               aes(x = week, y = dsi,
                   group = athlete_id,
                   color = sport)) +
  geom_hline(yintercept = 1.0, linetype = "dashed",
             color = "grey50", linewidth = 0.85) +
  geom_line(alpha = 0.60, linewidth = 0.80, na.rm = TRUE) +
  geom_point(alpha = 0.55, size = 1.6, na.rm = TRUE) +
  facet_wrap(~ sport, nrow = 1) +
  scale_color_manual(values = SPORT_COLS, guide = "none") +
  scale_x_continuous(breaks = seq(0, 35, 5),
                     expand = expansion(add = 0.5)) +
  scale_y_continuous(breaks = seq(0.4, 1.6, 0.1)) +
  labs(
    title    = "Individual DSI Trajectories Across the Season",
    subtitle = paste0("Athletes with ≥ 4 paired CMJ+IMTP sessions  |  ",
                      "Basketball n=", n_distinct(fig5_data$athlete_id[fig5_data$sport == "Basketball"]),
                      " | Rugby n=",   n_distinct(fig5_data$athlete_id[fig5_data$sport == "Rugby"])),
    x        = "Monitoring Week",
    y        = "Dynamic Strength Index (DSI)",
    caption  = "Dashed line = DSI reference 1.0. Each line represents one athlete's monitoring trajectory."
  ) +
  theme_nsca() +
  theme(strip.text = element_text(size = 13, face = "bold"))

save_fig(fig5, "Fig5_Individual_DSI_Trajectories.png", width = 16, height = 7)

# =============================================================================
# 14. TABLE 1 — DESCRIPTIVE STATISTICS BY SPORT
# =============================================================================
cat("--- Table 1: Descriptive statistics by sport ---\n")

perf_vars <- c(
  "Jump Height (cm)"        = "jump_height_cm",
  "Peak Power (W)"          = "peak_power_w",
  "RSI-modified (m/s)"      = "rsi_mod",
  "CMJ Peak Force (N)"      = "cmj_peak_force_n",
  "IMTP Peak Force (N)"     = "imtp_peak_force_n",
  "Dynamic Strength Index"  = "dsi"
)

desc_tbl_data <- master |>
  select(sport, all_of(unname(perf_vars))) |>
  pivot_longer(-sport, names_to = "var_col", values_to = "value") |>
  filter(!is.na(value)) |>
  group_by(sport, var_col) |>
  summarise(
    n      = n(),
    Mean   = mean(value, na.rm = TRUE),
    SD     = sd(value,   na.rm = TRUE),
    Median = median(value, na.rm = TRUE),
    Min    = min(value,  na.rm = TRUE),
    Max    = max(value,  na.rm = TRUE),
    .groups = "drop"
  ) |>
  mutate(
    `Mean ± SD` = sprintf("%.2f ± %.2f", Mean, SD),
    Median      = round(Median, 2),
    Min         = round(Min,    2),
    Max         = round(Max,    2),
    Variable    = names(perf_vars)[match(var_col, perf_vars)]
  ) |>
  select(Variable, Sport = sport, n, `Mean ± SD`, Median, Min, Max) |>
  mutate(Variable = factor(Variable, levels = names(perf_vars))) |>
  arrange(Variable, Sport)

tbl1_gt <- desc_tbl_data |>
  gt(groupname_col = "Variable") |>
  tab_header(
    title    = md("**Descriptive Statistics by Sport**"),
    subtitle = "Female university athletes — Basketball and Rugby"
  ) |>
  cols_label(
    Sport      = "Sport",
    n          = md("*n*"),
    `Mean ± SD` = "Mean ± SD",
    Median     = "Median",
    Min        = "Min",
    Max        = "Max"
  ) |>
  tab_style(
    style     = cell_text(weight = "bold"),
    locations = cells_column_labels()
  ) |>
  tab_style(
    style     = cell_text(weight = "bold", color = "grey30"),
    locations = cells_row_groups()
  ) |>
  tab_style(
    style     = cell_fill(color = "#EEF4FB"),
    locations = cells_body(rows = Sport == "Basketball")
  ) |>
  tab_source_note(
    md("CMJ = Countermovement Jump; IMTP = Isometric Mid-Thigh Pull; RSI-mod = RSI-modified;  \nDSI = Dynamic Strength Index (CMJ Peak Force / IMTP Peak Force). Values represent athlete-session observations.")
  ) |>
  opt_table_font(font = list(google_font("Roboto"), default_fonts())) |>
  tab_options(
    table.font.size                    = px(12),
    heading.title.font.size            = px(15),
    heading.subtitle.font.size         = px(12),
    column_labels.font.weight          = "bold",
    row_group.font.weight              = "bold",
    table.border.top.color             = "black",
    table.border.top.width             = px(2),
    column_labels.border.bottom.color  = "black",
    column_labels.border.bottom.width  = px(1),
    table.border.bottom.color          = "black",
    table.border.bottom.width          = px(2),
    source_notes.font.size             = px(10)
  )

tbl1_path <- file.path(TBL_DIR, "Table1_Descriptive_Statistics.html")
gtsave(tbl1_gt, tbl1_path)
cat(sprintf("  Saved: %s\n", basename(tbl1_path)))

# =============================================================================
# 15. TABLE 2 — ATHLETE DEMOGRAPHICS BY SPORT
# =============================================================================
cat("--- Table 2: Demographics by sport ---\n")

# Basketball: weight from CMJ file (mean per athlete), height not available
bas_demo_tbl <- master |>
  filter(sport == "Basketball") |>
  distinct(athlete_id, .keep_all = TRUE) |>
  summarise(
    Sport         = "Basketball",
    n             = n(),
    `BM Mean ± SD (kg)` = sprintf("%.1f ± %.1f",
                                   mean(bw_kg, na.rm = TRUE),
                                   sd(bw_kg,   na.rm = TRUE)),
    `Height Mean ± SD (cm)` = "Not available",
    `Age Mean ± SD (yr)`    = "Not available"
  )

# Rugby: weight & height & age from demographics file (one row per athlete)
rug_demo_tbl <- rug_demo |>
  summarise(
    Sport         = "Rugby",
    n             = n(),
    `BM Mean ± SD (kg)` = sprintf("%.1f ± %.1f",
                                   mean(weight_kg, na.rm = TRUE),
                                   sd(weight_kg,   na.rm = TRUE)),
    `Height Mean ± SD (cm)` = sprintf("%.1f ± %.1f",
                                       mean(height_cm, na.rm = TRUE),
                                       sd(height_cm,   na.rm = TRUE)),
    `Age Mean ± SD (yr)` = sprintf("%.1f ± %.1f",
                                    mean(age_years, na.rm = TRUE),
                                    sd(age_years,   na.rm = TRUE))
  )

demo_combined <- bind_rows(bas_demo_tbl, rug_demo_tbl)

tbl2_gt <- demo_combined |>
  gt() |>
  tab_header(
    title    = md("**Athlete Demographic Characteristics**"),
    subtitle = "Female university athletes — Basketball and Rugby"
  ) |>
  cols_label(
    Sport                    = "Sport",
    n                        = md("*n*"),
    `BM Mean ± SD (kg)`      = "Body Mass (kg)",
    `Height Mean ± SD (cm)`  = "Height (cm)",
    `Age Mean ± SD (yr)`     = "Age (yr)"
  ) |>
  tab_style(
    style     = cell_text(weight = "bold"),
    locations = cells_column_labels()
  ) |>
  tab_style(
    style     = cell_fill(color = "#EEF4FB"),
    locations = cells_body(rows = 1)
  ) |>
  tab_source_note(
    md("Values = mean ± SD. BM = Body Mass. Basketball height and age were not available in the provided data files.")
  ) |>
  opt_table_font(font = list(google_font("Roboto"), default_fonts())) |>
  tab_options(
    table.font.size                   = px(13),
    heading.title.font.size           = px(15),
    heading.subtitle.font.size        = px(12),
    column_labels.font.weight         = "bold",
    table.border.top.color            = "black",
    table.border.top.width            = px(2),
    column_labels.border.bottom.color = "black",
    column_labels.border.bottom.width = px(1),
    table.border.bottom.color         = "black",
    table.border.bottom.width         = px(2),
    source_notes.font.size            = px(10)
  )

tbl2_path <- file.path(TBL_DIR, "Table2_Demographics.html")
gtsave(tbl2_gt, tbl2_path)
cat(sprintf("  Saved: %s\n", basename(tbl2_path)))

# =============================================================================
# 16. TABLE 3 — WEEKLY NEUROMUSCULAR SUMMARY (both sports side-by-side)
# =============================================================================
cat("--- Table 3: Weekly neuromuscular summary ---\n")

fmt_ms <- function(m, s) ifelse(is.na(m), "—",
                                 sprintf("%.2f ± %.2f", m, s))
fmt_n0 <- function(m, s) ifelse(is.na(m), "—",
                                  sprintf("%.0f ± %.0f", m, s))
fmt_n3 <- function(m, s) ifelse(is.na(m), "—",
                                  sprintf("%.3f ± %.3f", m, s))

build_weekly_tbl <- function(sport_name) {
  weekly |>
    filter(sport == sport_name,
           !is.na(jump_height_cm__mean)) |>
    transmute(
      Week                        = week,
      N                           = n,
      `Jump Height (cm)`          = fmt_ms(jump_height_cm__mean,    jump_height_cm__sd),
      `RSI-mod (m/s)`             = fmt_ms(rsi_mod__mean,            rsi_mod__sd),
      `CMJ Peak Force (N)`        = fmt_n0(cmj_peak_force_n__mean,  cmj_peak_force_n__sd),
      `IMTP Peak Force (N)`       = fmt_n0(imtp_peak_force_n__mean, imtp_peak_force_n__sd),
      `DSI`                       = fmt_n3(dsi__mean,               dsi__sd)
    )
}

make_weekly_gt <- function(sport_name) {
  build_weekly_tbl(sport_name) |>
    gt() |>
    tab_header(
      title    = md(sprintf("**Weekly Neuromuscular Summary — %s**", sport_name)),
      subtitle = "Group mean ± SD per monitoring week"
    ) |>
    cols_label(
      Week                 = "Week",
      N                    = md("*n*"),
      `Jump Height (cm)`   = "Jump Height (cm)",
      `RSI-mod (m/s)`      = "RSI-mod (m/s)",
      `CMJ Peak Force (N)` = "CMJ PF (N)",
      `IMTP Peak Force (N)`= "IMTP PF (N)",
      DSI                  = "DSI"
    ) |>
    tab_style(
      style     = cell_text(weight = "bold"),
      locations = cells_column_labels()
    ) |>
    tab_style(
      style     = cell_fill(color = "#EEF4FB"),
      locations = cells_body(rows = seq(1, nrow(build_weekly_tbl(sport_name)), 2))
    ) |>
    tab_source_note(
      md("DSI = Dynamic Strength Index; CMJ PF = CMJ Peak Force; IMTP PF = IMTP Peak Force; RSI-mod = RSI-modified.")
    ) |>
    opt_table_font(font = list(google_font("Roboto"), default_fonts())) |>
    tab_options(
      table.font.size                   = px(11),
      heading.title.font.size           = px(14),
      heading.subtitle.font.size        = px(11),
      column_labels.font.weight         = "bold",
      table.border.top.color            = "black",
      table.border.top.width            = px(2),
      column_labels.border.bottom.color = "black",
      column_labels.border.bottom.width = px(1),
      table.border.bottom.color         = "black",
      table.border.bottom.width         = px(2),
      source_notes.font.size            = px(9)
    )
}

tbl3a_path <- file.path(TBL_DIR, "Table3a_Weekly_Summary_Basketball.html")
tbl3b_path <- file.path(TBL_DIR, "Table3b_Weekly_Summary_Rugby.html")
gtsave(make_weekly_gt("Basketball"), tbl3a_path)
gtsave(make_weekly_gt("Rugby"),      tbl3b_path)
cat(sprintf("  Saved: %s\n", basename(tbl3a_path)))
cat(sprintf("  Saved: %s\n", basename(tbl3b_path)))

# =============================================================================
# 17. FINAL SUMMARY REPORT
# =============================================================================
cat("\n=== Analysis Complete ===\n\n")

cat("OUTPUT FILES\n")
cat("Figures (", length(list.files(FIG_DIR)), "):\n", sep = "")
cat(paste0("  - ", list.files(FIG_DIR)), sep = "\n")
cat("\nTables  (", length(list.files(TBL_DIR)), "):\n", sep = "")
cat(paste0("  - ", list.files(TBL_DIR)), sep = "\n")

cat("\nKEY FINDINGS\n")
cat(sprintf("Basketball athletes: n = %d | Monitoring: weeks %d–%d\n",
            n_distinct(master$athlete_id[master$sport == "Basketball"]),
            min(master$week[master$sport == "Basketball"], na.rm = TRUE),
            max(master$week[master$sport == "Basketball"], na.rm = TRUE)))
cat(sprintf("Rugby athletes:      n = %d | Monitoring: weeks %d–%d\n",
            n_distinct(master$athlete_id[master$sport == "Rugby"]),
            min(master$week[master$sport == "Rugby"], na.rm = TRUE),
            max(master$week[master$sport == "Rugby"], na.rm = TRUE)))
cat(sprintf("\nDSI Basketball: %.3f ± %.3f | DSI < 1.0: %.1f%%\n",
            mean(dsi_data$dsi[dsi_data$sport == "Basketball"], na.rm = TRUE),
            sd(dsi_data$dsi[dsi_data$sport == "Basketball"],   na.rm = TRUE),
            mean(dsi_data$dsi[dsi_data$sport == "Basketball"] < 1.0) * 100))
cat(sprintf("DSI Rugby:      %.3f ± %.3f | DSI < 1.0: %.1f%%\n",
            mean(dsi_data$dsi[dsi_data$sport == "Rugby"], na.rm = TRUE),
            sd(dsi_data$dsi[dsi_data$sport == "Rugby"],   na.rm = TRUE),
            mean(dsi_data$dsi[dsi_data$sport == "Rugby"] < 1.0) * 100))
