# =============================================================================
# NSCA Poster — Hero Composite Figure
# Neuromuscular Performance: Basketball vs. Rugby
# Variables: DSI, Jump Height, CMJ & IMTP Peak Force, Propulsive Impulse
#
# Authors: Pablo Tadeo Ríos-Gallardo, PhD & Samuel Montalvo, PhD
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(lubridate)
  library(patchwork)
  library(scales)
  library(gt)
  library(gtExtras)
  if (requireNamespace("cowplot", quietly = TRUE)) library(cowplot)
})

BASE_DIR <- "/Volumes/ADATA HD710 PRO/fod_rugby"
RAW_DIR  <- file.path(BASE_DIR, "Raw_data")
FIG_DIR  <- file.path(BASE_DIR, "Figures")
TBL_DIR  <- file.path(BASE_DIR, "Tables")

cat("=== NSCA Hero Composite — Basketball & Rugby ===\n\n")

# =============================================================================
# 1. DATA PIPELINE (extended: + Jump Height, CMJ Impulse, IMTP Impulse 100ms)
# =============================================================================
cat("--- Loading & processing data ---\n")

read_utf8 <- function(f)
  read_csv(file.path(RAW_DIR, f), locale = locale(encoding = "UTF-8"),
           show_col_types = FALSE, name_repair = "minimal")
num <- function(x) suppressWarnings(as.numeric(x))

# ---- Basketball CMJ (+ Concentric Impulse) ----
bas_cmj <- read_utf8("CMJ_Basket - CMJ.csv") |>
  rename(
    athlete_name     = Name,
    athlete_id       = ExternalId,
    week             = Week,
    date_chr         = Date,
    bw_kg            = `BW [KG]`,
    jump_height_cm   = `Jump Height (Imp-Mom) [cm]`,
    peak_power_w     = `Peak Power [W]`,
    rsi_mod          = `RSI-modified [m/s]`,
    cmj_peak_force_n = `Concentric Peak Force [N]`,
    cmj_impulse_ns   = `Concentric Impulse [N s]`   # NEW
  ) |>
  select(athlete_id, athlete_name, week, date_chr, bw_kg,
         jump_height_cm, peak_power_w, rsi_mod,
         cmj_peak_force_n, cmj_impulse_ns) |>
  mutate(
    sport = "Basketball", week = as.integer(num(week)),
    date  = suppressWarnings(dmy(date_chr)),
    across(c(bw_kg, jump_height_cm, peak_power_w, rsi_mod,
             cmj_peak_force_n, cmj_impulse_ns), num),
    jump_height_cm   = ifelse(jump_height_cm < 0 | jump_height_cm > 120, NA, jump_height_cm),
    peak_power_w     = ifelse(peak_power_w < 0, NA, peak_power_w),
    rsi_mod          = ifelse(rsi_mod < 0 | rsi_mod > 5, NA, rsi_mod),
    cmj_peak_force_n = ifelse(cmj_peak_force_n <= 0, NA, cmj_peak_force_n),
    cmj_impulse_ns   = ifelse(cmj_impulse_ns   <= 0, NA, cmj_impulse_ns)
  ) |>
  filter(!is.na(athlete_id), !is.na(week))

# ---- Basketball IMTP (+ Absolute Impulse 100ms) ----
bas_imtp <- read_utf8("IMPT_Basket - Hoja 1.csv") |>
  rename(
    athlete_id             = ExternalId,
    week                   = Week,
    imtp_peak_force_n      = `Peak Vertical Force [N]`,
    imtp_impulse_100ms_ns  = `Absolute Impulse - 100ms [N s]`  # NEW
  ) |>
  select(athlete_id, week, imtp_peak_force_n, imtp_impulse_100ms_ns) |>
  mutate(
    week = as.integer(num(week)),
    across(c(imtp_peak_force_n, imtp_impulse_100ms_ns), num),
    imtp_peak_force_n     = ifelse(imtp_peak_force_n     <= 0, NA, imtp_peak_force_n),
    imtp_impulse_100ms_ns = ifelse(imtp_impulse_100ms_ns <= 0, NA, imtp_impulse_100ms_ns)
  ) |>
  filter(!is.na(athlete_id), !is.na(week))

# ---- Rugby Demographics ----
rug_demo <- read_utf8("CMJ_IMTP_RUGBY - Hoja 1.csv") |>
  rename_with(~ trimws(.x)) |>
  rename(athlete_id = externalId, weight_kg = `peso (kg)`,
         height_cm  = `altura(cm)`, birth_date = fecha_nacimiento) |>
  mutate(weight_kg = num(weight_kg), height_cm = num(height_cm)) |>
  select(athlete_id, weight_kg, height_cm)

# ---- Rugby CMJ (+ Concentric Impulse) ----
rug_cmj <- read_utf8("CMJ_Rugby - Hoja 1.csv") |>
  rename(
    athlete_name     = Name,
    athlete_id       = ExternalId,
    date_chr         = Date,
    bw_kg            = `BW [KG]`,
    jump_height_cm   = `Jump Height (Imp-Mom) [cm]`,
    peak_power_w     = `Peak Power [W]`,
    rsi_mod          = `RSI-modified [m/s]`,
    cmj_peak_force_n = `Concentric Peak Force [N]`,
    cmj_impulse_ns   = `Concentric Impulse [N s]`   # NEW
  ) |>
  select(athlete_id, athlete_name, week, date_chr, bw_kg,
         jump_height_cm, peak_power_w, rsi_mod,
         cmj_peak_force_n, cmj_impulse_ns) |>
  mutate(
    sport = "Rugby", week = as.integer(num(week)),
    date  = suppressWarnings(dmy(date_chr)),
    across(c(bw_kg, jump_height_cm, peak_power_w, rsi_mod,
             cmj_peak_force_n, cmj_impulse_ns), num),
    jump_height_cm   = ifelse(jump_height_cm < 0 | jump_height_cm > 120, NA, jump_height_cm),
    peak_power_w     = ifelse(peak_power_w < 0, NA, peak_power_w),
    rsi_mod          = ifelse(rsi_mod < 0 | rsi_mod > 5, NA, rsi_mod),
    cmj_peak_force_n = ifelse(cmj_peak_force_n <= 0, NA, cmj_peak_force_n),
    cmj_impulse_ns   = ifelse(cmj_impulse_ns   <= 0, NA, cmj_impulse_ns)
  ) |>
  filter(!is.na(athlete_id), !is.na(week))

# ---- Rugby IMTP (+ Absolute Impulse 100ms) ----
rug_imtp <- read_utf8("IMPT_Rugby - Hoja 1.csv") |>
  rename(
    athlete_id            = ExternalId,
    imtp_peak_force_n     = `Peak Vertical Force [N]`,
    imtp_impulse_100ms_ns = `Absolute Impulse - 100ms [N s]`   # NEW
  ) |>
  select(athlete_id, week, imtp_peak_force_n, imtp_impulse_100ms_ns) |>
  mutate(
    week = as.integer(num(week)),
    across(c(imtp_peak_force_n, imtp_impulse_100ms_ns), num),
    imtp_peak_force_n     = ifelse(imtp_peak_force_n     <= 0, NA, imtp_peak_force_n),
    imtp_impulse_100ms_ns = ifelse(imtp_impulse_100ms_ns <= 0, NA, imtp_impulse_100ms_ns)
  ) |>
  filter(!is.na(athlete_id), !is.na(week))

# ---- Aggregate per athlete × week ----
agg_cmj <- function(df)
  df |> group_by(athlete_id, athlete_name, sport, week) |>
    summarise(across(c(bw_kg, jump_height_cm, peak_power_w, rsi_mod,
                       cmj_peak_force_n, cmj_impulse_ns),
                     ~ mean(.x, na.rm = TRUE)), .groups = "drop") |>
    mutate(across(where(is.numeric), ~ ifelse(is.nan(.x), NA_real_, .x)))

agg_imtp <- function(df)
  df |> group_by(athlete_id, week) |>
    summarise(
      imtp_peak_force_n     = max(imtp_peak_force_n,     na.rm = TRUE),
      imtp_impulse_100ms_ns = max(imtp_impulse_100ms_ns, na.rm = TRUE),
      .groups = "drop"
    ) |>
    mutate(across(c(imtp_peak_force_n, imtp_impulse_100ms_ns),
                  ~ ifelse(is.infinite(.x), NA_real_, .x)))

# ---- Merge within sports + DSI ----
bas_perf <- left_join(agg_cmj(bas_cmj), agg_imtp(bas_imtp),
                      by = c("athlete_id", "week")) |>
  mutate(dsi = cmj_peak_force_n / imtp_peak_force_n,
         dsi = ifelse(dsi < 0.3 | dsi > 2.0, NA_real_, dsi),
         height_cm = NA_real_)

rug_perf <- left_join(agg_cmj(rug_cmj), agg_imtp(rug_imtp),
                      by = c("athlete_id", "week")) |>
  left_join(rug_demo, by = "athlete_id") |>
  mutate(bw_kg = coalesce(weight_kg, bw_kg),
         dsi   = cmj_peak_force_n / imtp_peak_force_n,
         dsi   = ifelse(dsi < 0.3 | dsi > 2.0, NA_real_, dsi)) |>
  select(-weight_kg)

master <- bind_rows(bas_perf, rug_perf) |>
  mutate(sport = factor(sport, levels = c("Basketball", "Rugby")))

# Quick data audit
cat(sprintf("Master: %d rows | BAS %d obs | RUG %d obs\n",
            nrow(master), sum(master$sport=="Basketball"), sum(master$sport=="Rugby")))
cat(sprintf("CMJ Impulse non-NA: BAS=%d | RUG=%d\n",
            sum(!is.na(master$cmj_impulse_ns[master$sport=="Basketball"])),
            sum(!is.na(master$cmj_impulse_ns[master$sport=="Rugby"]))))
cat(sprintf("IMTP Impulse 100ms non-NA: BAS=%d | RUG=%d\n\n",
            sum(!is.na(master$imtp_impulse_100ms_ns[master$sport=="Basketball"])),
            sum(!is.na(master$imtp_impulse_100ms_ns[master$sport=="Rugby"]))))

# =============================================================================
# 2. DESIGN TOKENS
# =============================================================================
SPORT_COLS  <- c(Basketball = "#0072B2", Rugby = "#E69F00")
SPORT_FILLS <- c(Basketball = "#0072B2", Rugby = "#E69F00")
COL_REF     <- "#D55E00"

# Base poster theme — used by all panels
theme_P <- function(bs = 18) {
  theme_classic(base_size = bs) %+replace%
    theme(
      axis.title         = element_text(face = "bold", size = bs),
      axis.text          = element_text(size = bs - 2, color = "grey15"),
      axis.line          = element_line(color = "grey30", linewidth = 0.5),
      axis.ticks         = element_line(color = "grey30", linewidth = 0.4),
      legend.position    = "none",        # collected at composite level
      strip.background   = element_rect(fill = "grey92", color = NA),
      strip.text         = element_text(face = "bold", size = bs - 1,
                                        margin = margin(4, 0, 4, 0)),
      panel.grid.major.y = element_line(color = "grey90", linewidth = 0.35),
      plot.subtitle      = element_text(size = bs - 3, color = "grey30",
                                        margin = margin(b = 6)),
      plot.caption       = element_text(size = bs - 6, color = "grey45", hjust = 0),
      plot.margin        = margin(10, 14, 6, 14)
    )
}

# =============================================================================
# 3. PANEL A — DSI by Sport (Box + Jitter + Force-Velocity zones)
# =============================================================================
cat("--- Panel A: DSI by Sport ---\n")

dsi_data <- master |> filter(!is.na(dsi))

dsi_pct <- dsi_data |>
  group_by(sport) |>
  summarise(pct = round(mean(dsi < 1.0) * 100, 1),
            y   = quantile(dsi, 0.98, na.rm = TRUE) + 0.03,
            .groups = "drop") |>
  mutate(lbl = paste0(pct, "% FV\nDeficit"))

p_A <- ggplot(dsi_data, aes(x = sport, y = dsi, fill = sport, color = sport)) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 1.0,
           fill = "#FFF5E0", alpha = 0.65) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 1.0, ymax = Inf,
           fill = "#E5F0F8", alpha = 0.65) +
  geom_boxplot(alpha = 0.60, outlier.shape = NA, width = 0.44,
               color = "grey20", linewidth = 0.65) +
  geom_jitter(width = 0.12, alpha = 0.50, size = 2.4, stroke = 0) +
  geom_hline(yintercept = 1.0, linetype = "dashed",
             color = COL_REF, linewidth = 1.1) +
  geom_text(data = dsi_pct,
            aes(x = sport, y = y, label = lbl, color = sport),
            inherit.aes = FALSE, size = 5.5, fontface = "bold", lineheight = 0.9) +
  annotate("text", x = 0.52, y = 0.80, label = "Force-Velocity\nDeficit Zone",
           hjust = 0, size = 5.2, color = "#9E6000", fontface = "italic", lineheight = 0.9) +
  annotate("text", x = 0.52, y = 1.19, label = "Adequate\nStrength Zone",
           hjust = 0, size = 5.2, color = "#004F8E", fontface = "italic", lineheight = 0.9) +
  scale_fill_manual(values = SPORT_FILLS, name = "Sport") +
  scale_color_manual(values = SPORT_COLS,  name = "Sport") +
  scale_y_continuous(breaks = seq(0.4, 1.6, 0.1),
                     expand = expansion(mult = 0.09)) +
  labs(subtitle = "Dynamic Strength Index by Sport  |  Dashed line = DSI 1.0",
       x = NULL, y = "Dynamic Strength Index (DSI)") +
  theme_P(18) +
  theme(axis.text.x = element_text(size = 19, face = "bold"))

# =============================================================================
# 4. PANEL B — Jump Height by Sport (Box + Jitter)
# =============================================================================
cat("--- Panel B: Jump Height by Sport ---\n")

jh_stats <- master |>
  filter(!is.na(jump_height_cm)) |>
  group_by(sport) |>
  summarise(m  = mean(jump_height_cm, na.rm = TRUE),
            sd = sd(jump_height_cm,   na.rm = TRUE),
            y  = quantile(jump_height_cm, 0.97, na.rm = TRUE) + 0.8,
            .groups = "drop") |>
  mutate(lbl = sprintf("%.1f ± %.1f cm", m, sd))

p_B <- ggplot(master |> filter(!is.na(jump_height_cm)),
              aes(x = sport, y = jump_height_cm, fill = sport, color = sport)) +
  geom_boxplot(alpha = 0.60, outlier.shape = NA, width = 0.44,
               color = "grey20", linewidth = 0.65) +
  geom_jitter(width = 0.12, alpha = 0.50, size = 2.4, stroke = 0) +
  geom_text(data = jh_stats,
            aes(x = sport, y = y, label = lbl, color = sport),
            inherit.aes = FALSE, size = 5.2, fontface = "bold") +
  scale_fill_manual(values = SPORT_FILLS, name = "Sport") +
  scale_color_manual(values = SPORT_COLS,  name = "Sport") +
  scale_y_continuous(expand = expansion(mult = 0.09)) +
  labs(subtitle = "Countermovement Jump Height by Sport  |  Mean ± SD annotated",
       x = NULL, y = "Jump Height (cm)") +
  theme_P(18) +
  theme(axis.text.x = element_text(size = 19, face = "bold"))

# =============================================================================
# 5. PANEL C — Longitudinal LOESS (4-panel, full width)
#    Variables: CMJ Peak Force | IMTP Peak Force |
#               CMJ Propulsive Impulse | IMTP Impulse 100ms
# =============================================================================
cat("--- Panel C: Longitudinal LOESS trends ---\n")

TEMPORAL_VARS <- c(
  "cmj_peak_force_n"     = "CMJ Peak Force (N)",
  "imtp_peak_force_n"    = "IMTP Peak Force (N)",
  "cmj_impulse_ns"       = "CMJ Propulsive Impulse (N·s)",
  "imtp_impulse_100ms_ns"= "IMTP Impulse 100ms (N·s)"
)

temporal_long <- master |>
  select(sport, week, all_of(names(TEMPORAL_VARS))) |>
  pivot_longer(
    cols      = all_of(names(TEMPORAL_VARS)),
    names_to  = "var_key",
    values_to = "value"
  ) |>
  filter(!is.na(value), !is.na(week)) |>
  mutate(
    Variable = factor(var_key,
                      levels = names(TEMPORAL_VARS),
                      labels = unname(TEMPORAL_VARS))
  )

# Observation counts per facet for caption
obs_counts <- temporal_long |>
  count(Variable, sport) |>
  pivot_wider(names_from = sport, values_from = n)

p_C <- ggplot(temporal_long,
              aes(x = week, y = value, color = sport, fill = sport)) +
  # Background raw points
  geom_point(alpha = 0.18, size = 1.3, stroke = 0, shape = 16) +
  # LOESS trend + 95% CI ribbon
  geom_smooth(
    method    = "loess",
    formula   = y ~ x,
    se        = TRUE,
    span      = 0.8,
    alpha     = 0.17,
    linewidth = 1.6
  ) +
  facet_wrap(~ Variable, scales = "free_y", ncol = 4) +
  scale_color_manual(values = SPORT_COLS,  name = "Sport") +
  scale_fill_manual(values  = SPORT_FILLS, name = "Sport") +
  scale_x_continuous(
    name   = "Monitoring Week",
    breaks = seq(0, 35, by = 5),
    expand = expansion(add = c(0.2, 0.8))
  ) +
  scale_y_continuous(expand = expansion(mult = 0.07)) +
  labs(
    subtitle = "Longitudinal Neuromuscular Trends: LOESS ± 95% CI  |  Ghost points = individual athlete-sessions",
    y        = NULL
  ) +
  theme_P(16) +
  theme(
    strip.text    = element_text(face = "bold", size = 15.5,
                                 margin = margin(5, 0, 5, 0)),
    panel.spacing = unit(1.1, "lines"),
    axis.title.x  = element_text(face = "bold", size = 16,
                                 margin = margin(t = 8))
  )

# =============================================================================
# 6. PANEL D — CMJ Propulsive Impulse vs. Jump Height (scatter)
# =============================================================================
cat("--- Panel D: CMJ Propulsive Impulse vs. Jump Height ---\n")

imp_jh_data <- master |>
  filter(!is.na(cmj_impulse_ns), !is.na(jump_height_cm))

r_imp_jh <- imp_jh_data |>
  group_by(sport) |>
  summarise(
    r   = cor(cmj_impulse_ns, jump_height_cm, use = "complete.obs"),
    x   = min(cmj_impulse_ns, na.rm = TRUE) * 1.02,
    y   = max(jump_height_cm, na.rm = TRUE) * 0.97,
    lbl = paste0("r = ", round(r, 2)),
    .groups = "drop"
  )

p_D <- ggplot(imp_jh_data,
              aes(x = cmj_impulse_ns, y = jump_height_cm,
                  color = sport, fill = sport)) +
  geom_point(alpha = 0.45, size = 2.2, stroke = 0) +
  geom_smooth(method = "lm", formula = y ~ x,
              se = TRUE, alpha = 0.15, linewidth = 1.3) +
  geom_text(data = r_imp_jh,
            aes(x = x, y = y, label = lbl, color = sport),
            inherit.aes = FALSE, size = 5.5, fontface = "bold", hjust = 0) +
  scale_color_manual(values = SPORT_COLS,  name = "Sport") +
  scale_fill_manual(values  = SPORT_FILLS, name = "Sport") +
  scale_x_continuous(labels = comma, expand = expansion(mult = 0.05)) +
  scale_y_continuous(expand = expansion(mult = 0.07)) +
  labs(
    subtitle = "CMJ Propulsive Impulse → Jump Height  |  r = Pearson correlation",
    x        = "CMJ Propulsive Impulse (N·s)",
    y        = "Jump Height (cm)"
  ) +
  theme_P(18) +
  theme(legend.position = "none")

# =============================================================================
# 7. PANEL E — CMJ Peak Force vs. IMTP Peak Force (scatter)
# =============================================================================
cat("--- Panel E: CMJ Peak Force vs. IMTP Peak Force ---\n")

force_data <- master |>
  filter(!is.na(cmj_peak_force_n), !is.na(imtp_peak_force_n))

r_forces <- force_data |>
  group_by(sport) |>
  summarise(
    r   = cor(imtp_peak_force_n, cmj_peak_force_n, use = "complete.obs"),
    x   = min(imtp_peak_force_n, na.rm = TRUE) * 1.01,
    y   = max(cmj_peak_force_n,  na.rm = TRUE) * 0.97,
    lbl = paste0("r = ", round(r, 2)),
    .groups = "drop"
  )

p_E <- ggplot(force_data,
              aes(x = imtp_peak_force_n, y = cmj_peak_force_n,
                  color = sport, fill = sport)) +
  geom_point(alpha = 0.45, size = 2.2, stroke = 0) +
  geom_smooth(method = "lm", formula = y ~ x,
              se = TRUE, alpha = 0.15, linewidth = 1.3) +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted",
              color = "grey55", linewidth = 0.9) +
  geom_text(data = r_forces,
            aes(x = x, y = y, label = lbl, color = sport),
            inherit.aes = FALSE, size = 5.5, fontface = "bold", hjust = 0) +
  scale_color_manual(values = SPORT_COLS,  name = "Sport") +
  scale_fill_manual(values  = SPORT_FILLS, name = "Sport") +
  scale_x_continuous(labels = comma, expand = expansion(mult = 0.05)) +
  scale_y_continuous(labels = comma, expand = expansion(mult = 0.05)) +
  labs(
    subtitle = "CMJ vs. IMTP Peak Force  |  Dotted line = identity (slope = 1)",
    x        = "IMTP Peak Force (N)",
    y        = "CMJ Peak Force (N)"
  ) +
  theme_P(18) +
  theme(legend.position = "none")

# =============================================================================
# 8. PATCHWORK COMPOSITE (3 rows)
#    Row 1 (h=5.5): [p_A | p_B]
#    Row 2 (h=7.5): [p_C — full width, 4 facets]
#    Row 3 (h=5.0): [p_D | p_E]
# =============================================================================
cat("--- Composing patchwork layout ---\n")

composite <- (p_A | p_B) /
             p_C          /
             (p_D | p_E)  +
  plot_layout(
    guides  = "collect",
    heights = c(5.5, 7.5, 5.0)
  ) +
  plot_annotation(
    tag_levels = "A",
    theme = theme(
      plot.tag = element_text(size = 30, face = "bold", color = "grey10",
                              margin = margin(r = 4, b = 4))
    )
  ) &
  theme(
    legend.position    = "bottom",
    legend.direction   = "horizontal",
    legend.title       = element_text(face = "bold", size = 19),
    legend.text        = element_text(size = 18),
    legend.key.width   = unit(2.2, "cm"),
    legend.key.height  = unit(0.65, "cm"),
    legend.margin      = margin(t = 6, b = 4),
    legend.box.spacing = unit(0, "pt")
  )

# =============================================================================
# 9. EXPORT HERO COMPOSITE (24 × 18 in, 600 dpi)
# =============================================================================
out_path <- file.path(FIG_DIR, "Hero_NSCA_Poster_Snapshot.png")

cat(sprintf(
  "\n--- Rendering composite: 24 × 18 in @ 600 dpi (14400 × 10800 px) ---\n"))
cat("This may take 3-7 minutes. Please wait...\n\n")
t0 <- Sys.time()

ggsave(
  filename = out_path,
  plot     = composite,
  width    = 24,
  height   = 18,
  dpi      = 600,
  bg       = "white",
  units    = "in"
)

elapsed    <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
file_mb    <- round(file.info(out_path)$size / 1024^2, 1)
cat(sprintf("Saved: %s\n", basename(out_path)))
cat(sprintf("Size : %.1f MB | Time: %.0f s\n\n", file_mb, elapsed))

# =============================================================================
# 10. UPDATED DESCRIPTIVE TABLES (now includes Jump Height + Impulses)
# =============================================================================
cat("--- Regenerating descriptive tables with new variables ---\n")

PERF_VARS <- c(
  "Jump Height (cm)"              = "jump_height_cm",
  "CMJ Peak Force (N)"            = "cmj_peak_force_n",
  "CMJ Propulsive Impulse (N·s)"  = "cmj_impulse_ns",
  "IMTP Peak Force (N)"           = "imtp_peak_force_n",
  "IMTP Impulse 100ms (N·s)"      = "imtp_impulse_100ms_ns",
  "RSI-modified (m/s)"            = "rsi_mod",
  "Dynamic Strength Index"        = "dsi"
)

# Helper: apply shared gt styling in one call
apply_gt_style <- function(gt_obj) {
  gt_obj |>
    tab_options(
      table.font.size                   = px(12),
      heading.title.font.size           = px(15),
      heading.subtitle.font.size        = px(12),
      column_labels.font.weight         = "bold",
      row_group.font.weight             = "bold",
      table.border.top.color            = "black",
      table.border.top.width            = px(2),
      column_labels.border.bottom.color = "black",
      column_labels.border.bottom.width = px(1),
      table.border.bottom.color         = "black",
      table.border.bottom.width         = px(2),
      source_notes.font.size            = px(10)
    )
}

# ---- Table 1: Descriptive Statistics by Sport ----
desc_data <- master |>
  select(sport, all_of(unname(PERF_VARS))) |>
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
    Variable    = names(PERF_VARS)[match(var_col, PERF_VARS)],
    Variable    = factor(Variable, levels = names(PERF_VARS))
  ) |>
  select(Variable, Sport = sport, n, `Mean ± SD`, Median, Min, Max) |>
  arrange(Variable, Sport)

tbl1_gt <- desc_data |>
  gt(groupname_col = "Variable") |>
  tab_header(
    title    = md("**Descriptive Statistics by Sport**"),
    subtitle = "Female university athletes — Basketball and Rugby"
  ) |>
  cols_label(Sport = "Sport", n = md("*n*"), Median = "Median",
             Min = "Min", Max = "Max") |>
  tab_style(style = cell_text(weight = "bold"),
            locations = cells_column_labels()) |>
  tab_style(style = cell_text(weight = "bold", color = "grey30"),
            locations = cells_row_groups()) |>
  tab_style(style = cell_fill(color = "#EEF4FB"),
            locations = cells_body(rows = Sport == "Basketball")) |>
  tab_source_note(md(paste0(
    "CMJ = Countermovement Jump; IMTP = Isometric Mid-Thigh Pull; ",
    "RSI-mod = RSI-modified; DSI = CMJ Peak Force / IMTP Peak Force.  \n",
    "Impulses: CMJ Propulsive = Concentric Impulse; IMTP Impulse = Absolute Impulse at 100ms."
  ))) |>
  opt_table_font(font = list(google_font("Roboto"), default_fonts())) |>
  apply_gt_style()

tbl1_path <- file.path(TBL_DIR, "Table1_Descriptive_Statistics.html")
gtsave(tbl1_gt, tbl1_path)
cat(sprintf("  Updated: %s\n", basename(tbl1_path)))

# ---- Table 2: Demographics ----
bas_demo_tbl <- master |>
  filter(sport == "Basketball") |>
  distinct(athlete_id, .keep_all = TRUE) |>
  summarise(
    Sport            = "Basketball",
    n                = n(),
    `BM (kg)`        = sprintf("%.1f ± %.1f", mean(bw_kg, na.rm=T), sd(bw_kg, na.rm=T)),
    `Height (cm)`    = "Not available",
    `Age (yr)`       = "Not available"
  )

rug_demo_tbl <- rug_demo |>
  left_join(
    read_utf8("CMJ_IMTP_RUGBY - Hoja 1.csv") |>
      rename_with(~ trimws(.x)) |>
      rename(athlete_id = externalId, birth_date = fecha_nacimiento) |>
      mutate(birth_date = suppressWarnings(dmy(birth_date)),
             age_years = as.numeric(difftime(as.Date("2025-06-01"),
                                             birth_date, units="days"))/365.25) |>
      select(athlete_id, age_years),
    by = "athlete_id"
  ) |>
  summarise(
    Sport         = "Rugby",
    n             = n(),
    `BM (kg)`     = sprintf("%.1f ± %.1f", mean(weight_kg, na.rm=T), sd(weight_kg, na.rm=T)),
    `Height (cm)` = sprintf("%.1f ± %.1f", mean(height_cm, na.rm=T), sd(height_cm, na.rm=T)),
    `Age (yr)`    = sprintf("%.1f ± %.1f", mean(age_years, na.rm=T), sd(age_years, na.rm=T))
  )

tbl2_gt <- bind_rows(bas_demo_tbl, rug_demo_tbl) |>
  gt() |>
  tab_header(title    = md("**Athlete Demographic Characteristics**"),
             subtitle = "Female university athletes — Basketball and Rugby") |>
  cols_label(Sport = "Sport", n = md("*n*")) |>
  tab_style(style = cell_text(weight = "bold"), locations = cells_column_labels()) |>
  tab_style(style = cell_fill(color = "#EEF4FB"), locations = cells_body(rows = 1)) |>
  tab_source_note(md("Values = mean ± SD. BM = Body Mass. Basketball height/age not in provided files.")) |>
  opt_table_font(font = list(google_font("Roboto"), default_fonts())) |>
  apply_gt_style()

tbl2_path <- file.path(TBL_DIR, "Table2_Demographics.html")
gtsave(tbl2_gt, tbl2_path)
cat(sprintf("  Updated: %s\n", basename(tbl2_path)))

# =============================================================================
# 11. FINAL SUMMARY
# =============================================================================
cat("\n=== Complete ===\n")
cat(sprintf("Hero composite : %s (%.1f MB)\n",
            basename(out_path), file_mb))
cat(sprintf("Dimensions     : 24 × 18 in @ 600 dpi (14400 × 10800 px)\n"))
cat(sprintf("Tables updated : %s\n", paste(basename(c(tbl1_path, tbl2_path)),
                                            collapse = ", ")))
cat(sprintf("\nKey statistics:\n"))
cat(sprintf("  DSI Basketball: %.3f ± %.3f | DSI < 1.0: %.1f%%\n",
            mean(dsi_data$dsi[dsi_data$sport=="Basketball"], na.rm=TRUE),
            sd(dsi_data$dsi[dsi_data$sport=="Basketball"],   na.rm=TRUE),
            mean(dsi_data$dsi[dsi_data$sport=="Basketball"] < 1.0)*100))
cat(sprintf("  DSI Rugby:      %.3f ± %.3f | DSI < 1.0: %.1f%%\n",
            mean(dsi_data$dsi[dsi_data$sport=="Rugby"], na.rm=TRUE),
            sd(dsi_data$dsi[dsi_data$sport=="Rugby"],   na.rm=TRUE),
            mean(dsi_data$dsi[dsi_data$sport=="Rugby"] < 1.0)*100))
cat(sprintf("  CMJ Impulse Basketball: %.1f ± %.1f N·s\n",
            mean(master$cmj_impulse_ns[master$sport=="Basketball"], na.rm=TRUE),
            sd(master$cmj_impulse_ns[master$sport=="Basketball"],   na.rm=TRUE)))
cat(sprintf("  CMJ Impulse Rugby:      %.1f ± %.1f N·s\n",
            mean(master$cmj_impulse_ns[master$sport=="Rugby"], na.rm=TRUE),
            sd(master$cmj_impulse_ns[master$sport=="Rugby"],   na.rm=TRUE)))
