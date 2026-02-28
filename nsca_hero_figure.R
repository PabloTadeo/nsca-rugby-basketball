# =============================================================================
# Hero Figure — NSCA Poster
# Temporal Neuromuscular Performance Trajectories: Basketball vs. Rugby
#
# Authors: Pablo Tadeo Ríos-Gallardo, PhD & Samuel Montalvo, PhD
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(lubridate)
  library(scales)
})

BASE_DIR <- "/Volumes/ADATA HD710 PRO/fod_rugby"
RAW_DIR  <- file.path(BASE_DIR, "Raw_data")
FIG_DIR  <- file.path(BASE_DIR, "Figures")

cat("=== Hero Figure — NSCA Poster ===\n\n")

# =============================================================================
# 1. DATA PIPELINE (identical to nsca_poster_analysis.R)
# =============================================================================
cat("--- Loading & building master dataset ---\n")

read_utf8 <- function(f)
  read_csv(file.path(RAW_DIR, f), locale = locale(encoding = "UTF-8"),
           show_col_types = FALSE, name_repair = "minimal")

num <- function(x) suppressWarnings(as.numeric(x))

# --- Basketball CMJ ---
bas_cmj <- read_utf8("CMJ_Basket - CMJ.csv") |>
  rename(athlete_name = Name, athlete_id = ExternalId, week = Week,
         date_chr = Date, bw_kg = `BW [KG]`,
         jump_height_cm = `Jump Height (Imp-Mom) [cm]`,
         peak_power_w = `Peak Power [W]`, rsi_mod = `RSI-modified [m/s]`,
         cmj_peak_force_n = `Concentric Peak Force [N]`) |>
  select(athlete_id, athlete_name, week, date_chr, bw_kg,
         jump_height_cm, peak_power_w, rsi_mod, cmj_peak_force_n) |>
  mutate(sport = "Basketball", week = as.integer(num(week)),
         date = suppressWarnings(dmy(date_chr)),
         across(c(bw_kg, jump_height_cm, peak_power_w,
                  rsi_mod, cmj_peak_force_n), num),
         jump_height_cm   = ifelse(jump_height_cm < 0 | jump_height_cm > 120, NA, jump_height_cm),
         peak_power_w     = ifelse(peak_power_w < 0, NA, peak_power_w),
         rsi_mod          = ifelse(rsi_mod < 0 | rsi_mod > 5, NA, rsi_mod),
         cmj_peak_force_n = ifelse(cmj_peak_force_n <= 0, NA, cmj_peak_force_n)) |>
  filter(!is.na(athlete_id), !is.na(week))

# --- Basketball IMTP ---
bas_imtp <- read_utf8("IMPT_Basket - Hoja 1.csv") |>
  rename(athlete_id = ExternalId, week = Week,
         imtp_peak_force_n = `Peak Vertical Force [N]`) |>
  select(athlete_id, week, imtp_peak_force_n) |>
  mutate(week = as.integer(num(week)),
         imtp_peak_force_n = num(imtp_peak_force_n),
         imtp_peak_force_n = ifelse(imtp_peak_force_n <= 0, NA, imtp_peak_force_n)) |>
  filter(!is.na(athlete_id), !is.na(week))

# --- Rugby Demographics ---
rug_demo <- read_utf8("CMJ_IMTP_RUGBY - Hoja 1.csv") |>
  rename_with(~ trimws(.x)) |>
  rename(athlete_id = externalId, weight_kg = `peso (kg)`,
         height_cm = `altura(cm)`, birth_date = fecha_nacimiento) |>
  mutate(weight_kg = num(weight_kg), height_cm = num(height_cm)) |>
  select(athlete_id, weight_kg, height_cm)

# --- Rugby CMJ ---
rug_cmj <- read_utf8("CMJ_Rugby - Hoja 1.csv") |>
  rename(athlete_name = Name, athlete_id = ExternalId,
         date_chr = Date, bw_kg = `BW [KG]`,
         jump_height_cm = `Jump Height (Imp-Mom) [cm]`,
         peak_power_w = `Peak Power [W]`, rsi_mod = `RSI-modified [m/s]`,
         cmj_peak_force_n = `Concentric Peak Force [N]`) |>
  select(athlete_id, athlete_name, week, date_chr, bw_kg,
         jump_height_cm, peak_power_w, rsi_mod, cmj_peak_force_n) |>
  mutate(sport = "Rugby", week = as.integer(num(week)),
         date = suppressWarnings(dmy(date_chr)),
         across(c(bw_kg, jump_height_cm, peak_power_w,
                  rsi_mod, cmj_peak_force_n), num),
         jump_height_cm   = ifelse(jump_height_cm < 0 | jump_height_cm > 120, NA, jump_height_cm),
         peak_power_w     = ifelse(peak_power_w < 0, NA, peak_power_w),
         rsi_mod          = ifelse(rsi_mod < 0 | rsi_mod > 5, NA, rsi_mod),
         cmj_peak_force_n = ifelse(cmj_peak_force_n <= 0, NA, cmj_peak_force_n)) |>
  filter(!is.na(athlete_id), !is.na(week))

# --- Rugby IMTP ---
rug_imtp <- read_utf8("IMPT_Rugby - Hoja 1.csv") |>
  rename(athlete_id = ExternalId,
         imtp_peak_force_n = `Peak Vertical Force [N]`) |>
  select(athlete_id, week, imtp_peak_force_n) |>
  mutate(week = as.integer(num(week)),
         imtp_peak_force_n = num(imtp_peak_force_n),
         imtp_peak_force_n = ifelse(imtp_peak_force_n <= 0, NA, imtp_peak_force_n)) |>
  filter(!is.na(athlete_id), !is.na(week))

# --- Aggregate: CMJ (mean) | IMTP (max = best trial) ---
agg_cmj <- function(df)
  df |> group_by(athlete_id, athlete_name, sport, week) |>
    summarise(across(c(bw_kg, jump_height_cm, peak_power_w,
                       rsi_mod, cmj_peak_force_n),
                     ~ mean(.x, na.rm = TRUE)), .groups = "drop") |>
    mutate(across(where(is.numeric), ~ ifelse(is.nan(.x), NA_real_, .x)))

agg_imtp <- function(df)
  df |> group_by(athlete_id, week) |>
    summarise(imtp_peak_force_n = max(imtp_peak_force_n, na.rm = TRUE),
              .groups = "drop") |>
    mutate(imtp_peak_force_n = ifelse(is.infinite(imtp_peak_force_n),
                                      NA_real_, imtp_peak_force_n))

# --- Merge within sports ---
bas_perf <- left_join(agg_cmj(bas_cmj), agg_imtp(bas_imtp),
                      by = c("athlete_id", "week")) |>
  mutate(dsi = cmj_peak_force_n / imtp_peak_force_n,
         dsi = ifelse(dsi < 0.3 | dsi > 2.0, NA_real_, dsi))

rug_perf <- left_join(agg_cmj(rug_cmj), agg_imtp(rug_imtp),
                      by = c("athlete_id", "week")) |>
  left_join(rug_demo, by = "athlete_id") |>
  mutate(bw_kg = coalesce(weight_kg, bw_kg),
         dsi   = cmj_peak_force_n / imtp_peak_force_n,
         dsi   = ifelse(dsi < 0.3 | dsi > 2.0, NA_real_, dsi)) |>
  select(-weight_kg)

# --- Master dataset ---
master <- bind_rows(bas_perf, rug_perf) |>
  mutate(sport = factor(sport, levels = c("Basketball", "Rugby")))

cat(sprintf("Master: %d rows | Basketball %d obs | Rugby %d obs\n\n",
            nrow(master),
            sum(master$sport == "Basketball"),
            sum(master$sport == "Rugby")))

# =============================================================================
# 2. RESHAPE TO LONG FORMAT (4 key variables)
# =============================================================================
cat("--- Reshaping to long format ---\n")

# Ordered so the 2×2 grid reads: CMJ PF | IMTP PF
#                                  RSI-mod | DSI
PANEL_LEVELS <- c("cmj_peak_force_n",  "imtp_peak_force_n",
                  "rsi_mod",           "dsi")
PANEL_LABELS <- c("CMJ Peak Force (N)", "IMTP Peak Force (N)",
                  "RSI-modified (m/s)", "Dynamic Strength Index (DSI)")

hero_long <- master |>
  select(sport, week, athlete_id,
         cmj_peak_force_n, imtp_peak_force_n, rsi_mod, dsi) |>
  pivot_longer(
    cols      = all_of(PANEL_LEVELS),
    names_to  = "variable_key",
    values_to = "value"
  ) |>
  filter(!is.na(value), !is.na(week)) |>
  mutate(
    Variable = factor(variable_key,
                      levels = PANEL_LEVELS,
                      labels = PANEL_LABELS)
  )

cat(sprintf("Long format: %d rows across %d panels\n",
            nrow(hero_long), n_distinct(hero_long$Variable)))
hero_long |>
  count(Variable, sport) |>
  print()
cat("\n")

# =============================================================================
# 3. DESIGN TOKENS
# =============================================================================

# Okabe-Ito — consistent with the rest of the poster
SPORT_COLS  <- c(Basketball = "#0072B2", Rugby = "#E69F00")
SPORT_FILLS <- c(Basketball = "#0072B2", Rugby = "#E69F00")

# =============================================================================
# 4. BUILD HERO FIGURE
# =============================================================================
cat("--- Building hero figure ---\n")

hero <- ggplot(hero_long,
               aes(x = week, y = value,
                   color = sport, fill = sport)) +

  # ── Layer 1: raw athlete-week observations (ghost points) ──────────────────
  geom_point(
    alpha  = 0.25,
    size   = 1.4,
    stroke = 0,
    shape  = 16
  ) +

  # ── Layer 2: LOESS trend + shaded 95% CI ───────────────────────────────────
  geom_smooth(
    method  = "loess",
    formula = y ~ x,
    se      = TRUE,
    span    = 0.8,
    alpha   = 0.18,           # CI ribbon transparency
    linewidth = 1.6
  ) +

  # ── Facets ─────────────────────────────────────────────────────────────────
  facet_wrap(~ Variable,
             scales = "free_y",
             ncol   = 2) +

  # ── Scales ─────────────────────────────────────────────────────────────────
  scale_color_manual(values = SPORT_COLS,  name = "Sport") +
  scale_fill_manual(values  = SPORT_FILLS, name = "Sport") +
  scale_x_continuous(
    name   = "Monitoring Week",
    breaks = seq(0, 35, by = 5),
    expand = expansion(add = c(0.3, 0.8))
  ) +
  scale_y_continuous(expand = expansion(mult = 0.07)) +

  # ── Labels ─────────────────────────────────────────────────────────────────
  labs(
    title   = NULL,       # title goes on the physical poster
    y       = NULL        # facet strip labels carry variable identity
  ) +

  # ── Theme ──────────────────────────────────────────────────────────────────
  theme_classic(base_size = 18) +
  theme(
    # Facet strips — grey fill, bold text, no border
    strip.background  = element_rect(fill = "grey91", color = NA),
    strip.text        = element_text(face = "bold", size = 17,
                                     margin = margin(t = 5, b = 5)),

    # Axes
    axis.title.x      = element_text(face = "bold", size = 18,
                                     margin = margin(t = 10)),
    axis.text         = element_text(size = 15, color = "grey20"),
    axis.line         = element_line(color = "grey30", linewidth = 0.5),
    axis.ticks        = element_line(color = "grey30", linewidth = 0.4),

    # Subtle horizontal grid to guide the eye
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.35),
    panel.spacing      = unit(1.2, "lines"),

    # Legend — bottom, clean
    legend.position    = "bottom",
    legend.direction   = "horizontal",
    legend.title       = element_text(face = "bold", size = 17),
    legend.text        = element_text(size = 16),
    legend.key.width   = unit(2.0, "cm"),
    legend.key.height  = unit(0.6, "cm"),
    legend.margin      = margin(t = 8),

    # Overall margins
    plot.margin        = margin(14, 18, 10, 18)
  ) +

  # ── Legend guide tweaks ────────────────────────────────────────────────────
  guides(
    color = guide_legend(override.aes = list(
      alpha     = 1,
      size      = 4,
      linewidth = 2
    )),
    fill = "none"
  )

# =============================================================================
# 5. EXPORT AT 600 DPI
# =============================================================================
out_path <- file.path(FIG_DIR, "Hero_Temporal_Facet.png")

ggsave(
  filename = out_path,
  plot     = hero,
  width    = 12,
  height   = 10,
  dpi      = 600,
  bg       = "white",
  units    = "in"
)

file_size_mb <- round(file.info(out_path)$size / 1024^2, 1)
cat(sprintf("\n  Saved: Hero_Temporal_Facet.png\n"))
cat(sprintf("  Dimensions : 12 × 10 in\n"))
cat(sprintf("  Resolution : 600 dpi  (pixel size: 7200 × 6000 px)\n"))
cat(sprintf("  File size  : %s MB\n", file_size_mb))
cat("\n=== Hero Figure Complete ===\n")
