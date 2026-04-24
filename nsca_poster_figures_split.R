# =============================================================================
# NSCA Poster — Split Hero Figures
#
# Produces TWO publication-quality figures from the same LMM pipeline:
#
#  Hero_Fig1_Baseline_Mechanisms.png  (24 × 11 in @ 600 dpi)
#    Row 1: DSI Boxplot  |  Jump Height Boxplot          (Baseline)
#    Row 2: CMJ Impulse × JH Scatter  |  CMJ × IMTP Scatter  (Mechanisms)
#
#  Hero_Fig2_LMM_Trajectories.png    (24 × 7 in @ 600 dpi)
#    Single row: DSI | Jump Height | CMJ PF | IMTP PF | Peak Power
#    Lines = LMM predicted marginal means  |  Ribbon = 95% CI
#
# All LMMs, predictions, coord_cartesian Y-clips, and themes are identical
# to nsca_poster_hero_composite_FINAL.R — only patchwork assembly changes.
#
# Authors: Pablo Tadeo Ríos-Gallardo, PhD & Samuel Montalvo, PhD
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(lubridate)
  library(lme4)
  library(lmerTest)
  library(ggeffects)
  library(patchwork)
  library(scales)
  library(rmcorr)   # repeated-measures correlation (Bakdash & Marusich, 2017)
})

BASE_DIR <- "/Volumes/ADATA HD710 PRO/fod_rugby"
RAW_DIR  <- file.path(BASE_DIR, "Raw_data")
FIG_DIR  <- file.path(BASE_DIR, "Figures")

cat("=== NSCA Poster — Split Figures (Fig1 + Fig2) ===\n\n")

# =============================================================================
# 1. DATA PIPELINE
# =============================================================================
cat("--- 1. Loading data ---\n")

read_utf8 <- function(f)
  read_csv(file.path(RAW_DIR, f),
           locale        = locale(encoding = "UTF-8"),
           show_col_types = FALSE,
           name_repair   = "minimal")

num <- function(x) suppressWarnings(as.numeric(x))

bas_cmj <- read_utf8("CMJ_Basket - CMJ.csv") |>
  rename(
    athlete_name     = Name,       athlete_id       = ExternalId,
    week             = Week,       bw_kg            = `BW [KG]`,
    jump_height_cm   = `Jump Height (Imp-Mom) [cm]`,
    peak_power_w     = `Peak Power [W]`,
    rsi_mod          = `RSI-modified [m/s]`,
    cmj_peak_force_n = `Concentric Peak Force [N]`,
    cmj_impulse_ns   = `Concentric Impulse [N s]`
  ) |>
  select(athlete_id, athlete_name, week, bw_kg,
         jump_height_cm, peak_power_w, rsi_mod,
         cmj_peak_force_n, cmj_impulse_ns) |>
  mutate(
    sport            = "Basketball",
    week             = as.integer(num(week)),
    across(c(bw_kg, jump_height_cm, peak_power_w, rsi_mod,
             cmj_peak_force_n, cmj_impulse_ns), num),
    jump_height_cm   = ifelse(jump_height_cm < 0 | jump_height_cm > 120, NA, jump_height_cm),
    peak_power_w     = ifelse(peak_power_w   <= 0, NA, peak_power_w),
    rsi_mod          = ifelse(rsi_mod < 0 | rsi_mod > 5, NA, rsi_mod),
    cmj_peak_force_n = ifelse(cmj_peak_force_n <= 0, NA, cmj_peak_force_n),
    cmj_impulse_ns   = ifelse(cmj_impulse_ns   <= 0, NA, cmj_impulse_ns)
  ) |>
  filter(!is.na(athlete_id), !is.na(week))

bas_imtp <- read_utf8("IMPT_Basket - Hoja 1.csv") |>
  rename(athlete_id = ExternalId, week = Week,
         imtp_peak_force_n     = `Peak Vertical Force [N]`,
         imtp_impulse_100ms_ns = `Absolute Impulse - 100ms [N s]`) |>
  select(athlete_id, week, imtp_peak_force_n, imtp_impulse_100ms_ns) |>
  mutate(week = as.integer(num(week)),
         across(c(imtp_peak_force_n, imtp_impulse_100ms_ns), num),
         imtp_peak_force_n     = ifelse(imtp_peak_force_n     <= 0, NA, imtp_peak_force_n),
         imtp_impulse_100ms_ns = ifelse(imtp_impulse_100ms_ns <= 0, NA, imtp_impulse_100ms_ns)) |>
  filter(!is.na(athlete_id), !is.na(week))

rug_demo <- read_utf8("CMJ_IMTP_RUGBY - Hoja 1.csv") |>
  rename_with(~ trimws(.x)) |>
  rename(athlete_id = externalId,
         weight_kg  = `peso (kg)`,
         height_cm  = `altura(cm)`) |>
  mutate(weight_kg = num(weight_kg), height_cm = num(height_cm)) |>
  select(athlete_id, weight_kg, height_cm)

rug_cmj <- read_utf8("CMJ_Rugby - Hoja 1.csv") |>
  rename(
    athlete_name     = Name,       athlete_id       = ExternalId,
    bw_kg            = `BW [KG]`,
    jump_height_cm   = `Jump Height (Imp-Mom) [cm]`,
    peak_power_w     = `Peak Power [W]`,
    rsi_mod          = `RSI-modified [m/s]`,
    cmj_peak_force_n = `Concentric Peak Force [N]`,
    cmj_impulse_ns   = `Concentric Impulse [N s]`
  ) |>
  select(athlete_id, athlete_name, week, bw_kg,
         jump_height_cm, peak_power_w, rsi_mod,
         cmj_peak_force_n, cmj_impulse_ns) |>
  mutate(
    sport            = "Rugby",
    week             = as.integer(num(week)),
    across(c(bw_kg, jump_height_cm, peak_power_w, rsi_mod,
             cmj_peak_force_n, cmj_impulse_ns), num),
    jump_height_cm   = ifelse(jump_height_cm < 0 | jump_height_cm > 120, NA, jump_height_cm),
    peak_power_w     = ifelse(peak_power_w   <= 0, NA, peak_power_w),
    rsi_mod          = ifelse(rsi_mod < 0 | rsi_mod > 5, NA, rsi_mod),
    cmj_peak_force_n = ifelse(cmj_peak_force_n <= 0, NA, cmj_peak_force_n),
    cmj_impulse_ns   = ifelse(cmj_impulse_ns   <= 0, NA, cmj_impulse_ns)
  ) |>
  filter(!is.na(athlete_id), !is.na(week))

rug_imtp <- read_utf8("IMPT_Rugby - Hoja 1.csv") |>
  rename(athlete_id = ExternalId,
         imtp_peak_force_n     = `Peak Vertical Force [N]`,
         imtp_impulse_100ms_ns = `Absolute Impulse - 100ms [N s]`) |>
  select(athlete_id, week, imtp_peak_force_n, imtp_impulse_100ms_ns) |>
  mutate(week = as.integer(num(week)),
         across(c(imtp_peak_force_n, imtp_impulse_100ms_ns), num),
         imtp_peak_force_n     = ifelse(imtp_peak_force_n     <= 0, NA, imtp_peak_force_n),
         imtp_impulse_100ms_ns = ifelse(imtp_impulse_100ms_ns <= 0, NA, imtp_impulse_100ms_ns)) |>
  filter(!is.na(athlete_id), !is.na(week))

agg_cmj <- function(df)
  df |> group_by(athlete_id, athlete_name, sport, week) |>
    summarise(
      across(c(bw_kg, jump_height_cm, peak_power_w, rsi_mod,
               cmj_peak_force_n, cmj_impulse_ns),
             ~ mean(.x, na.rm = TRUE)), .groups = "drop"
    ) |>
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

cat(sprintf("Master: %d rows | Basketball n=%d (%d obs) | Rugby n=%d (%d obs)\n\n",
            nrow(master),
            n_distinct(master$athlete_id[master$sport == "Basketball"]),
            sum(master$sport == "Basketball"),
            n_distinct(master$athlete_id[master$sport == "Rugby"]),
            sum(master$sport == "Rugby")))

# =============================================================================
# 2. DESIGN TOKENS & SHARED THEME
# =============================================================================
SPORT_COLS  <- c(Basketball = "#0072B2", Rugby = "#E69F00")
SPORT_FILLS <- c(Basketball = "#0072B2", Rugby = "#E69F00")
COL_REF     <- "#D55E00"

# Shared poster theme ─ used by ALL panels across both figures
theme_P <- function(bs = 18) {
  theme_classic(base_size = bs) %+replace%
    theme(
      axis.title         = element_text(face = "bold", size = bs),
      axis.text          = element_text(size = bs - 2, color = "grey15"),
      axis.line          = element_line(color = "grey30", linewidth = 0.5),
      axis.ticks         = element_line(color = "grey30", linewidth = 0.4),
      legend.position    = "none",
      strip.background   = element_rect(fill = "grey92", color = NA),
      strip.text         = element_text(face = "bold", size = bs - 1,
                                        margin = margin(4, 0, 4, 0)),
      panel.grid.major.y = element_line(color = "grey90", linewidth = 0.35),
      plot.subtitle      = element_text(size = bs - 3, color = "grey30",
                                        margin = margin(b = 6)),
      plot.margin        = margin(8, 12, 6, 12)   # tight margins throughout
    )
}

# =============================================================================
# 3. LMM FITTING (5 outcomes; BOBYQA; scaled variables back-transformed)
# =============================================================================
cat("--- 3. Fitting LMMs ---\n")

LMM_CTRL <- lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))

fit_lmm_silent <- function(outcome_col, label, data, scale_factor = 1) {
  df <- data |>
    filter(!is.na(.data[[outcome_col]]), !is.na(week), !is.na(sport)) |>
    mutate(y_fit = .data[[outcome_col]] / scale_factor)

  fit <- lmer(y_fit ~ sport * week + (1 | athlete_id),
              data = df, REML = TRUE, control = LMM_CTRL)

  conv <- tryCatch(fit@optinfo$conv$lme4$messages, error = function(e) NULL)
  if (!is.null(conv) && length(conv) > 0) {
    fit <- lmer(y_fit ~ sport * week + (1 | athlete_id),
                data = df, REML = TRUE,
                control = lmerControl(optimizer = "Nelder_Mead",
                                      optCtrl   = list(maxfun = 2e5)))
    conv <- tryCatch(fit@optinfo$conv$lme4$messages, error = function(e) NULL)
  }
  ok  <- is.null(conv) || length(conv) == 0
  vc  <- as.data.frame(VarCorr(fit))
  icc <- vc$sdcor[vc$grp == "athlete_id"]^2 /
         (vc$sdcor[vc$grp == "athlete_id"]^2 + vc$sdcor[vc$grp == "Residual"]^2)

  cat(sprintf("  %s %-26s  n=%d  ICC=%.3f%s\n",
              ifelse(ok, "\u2713", "\u26A0"), label, nrow(df), icc,
              ifelse(scale_factor != 1,
                     sprintf("  [\u00f7%g back-transformed]", scale_factor), "")))
  fit
}

m_dsi  <- fit_lmm_silent("dsi",              "DSI",              master, 1)
m_jh   <- fit_lmm_silent("jump_height_cm",   "Jump Height (cm)", master, 1)
m_cmj  <- fit_lmm_silent("cmj_peak_force_n", "CMJ PF (N)",       master, 1000)
m_imtp <- fit_lmm_silent("imtp_peak_force_n","IMTP PF (N)",      master, 1000)
m_pp   <- fit_lmm_silent("peak_power_w",     "Peak Power (W)",   master, 1000)
cat("\n")

# =============================================================================
# 4. LMM PREDICTION EXTRACTION
#    Population-average fixed effects (ggpredict default for merMod)
#    Clipped to observed week range per sport — no extrapolation
# =============================================================================
cat("--- 4. Extracting LMM predictions ---\n")

get_predictions <- function(fit, var_col, scale_factor = 1) {
  ranges <- master |>
    filter(!is.na(.data[[var_col]]), !is.na(week)) |>
    group_by(sport) |>
    summarise(w_min = min(week), w_max = max(week), .groups = "drop")

  pred_raw <- ggpredict(
    fit,
    terms = list(week  = seq(min(ranges$w_min), max(ranges$w_max)),
                 sport = c("Basketball", "Rugby"))
  )

  as.data.frame(pred_raw) |>
    rename(week = x) |>
    mutate(sport     = factor(group, levels = c("Basketball", "Rugby")),
           predicted = predicted * scale_factor,
           conf.low  = conf.low  * scale_factor,
           conf.high = conf.high * scale_factor) |>
    left_join(ranges, by = "sport") |>
    filter(week >= w_min, week <= w_max)
}

pred_dsi  <- get_predictions(m_dsi,  "dsi",              1)
pred_jh   <- get_predictions(m_jh,   "jump_height_cm",   1)
pred_cmj  <- get_predictions(m_cmj,  "cmj_peak_force_n", 1000)
pred_imtp <- get_predictions(m_imtp, "imtp_peak_force_n",1000)
pred_pp   <- get_predictions(m_pp,   "peak_power_w",     1000)
cat("  Done.\n\n")

# =============================================================================
# 5. PANEL BUILDERS
# =============================================================================

# ---- 5a. LMM trajectory panel (used in Fig 2) ----
# Ghost points + LMM ribbon + LMM line + coord_cartesian Y clip
make_lmm_panel <- function(pred_df, var_col, var_label, bs = 16) {
  df_raw <- master |>
    filter(!is.na(.data[[var_col]]), !is.na(week)) |>
    mutate(y_val = .data[[var_col]])

  q_lo  <- quantile(df_raw$y_val, 0.01, na.rm = TRUE)
  q_hi  <- quantile(df_raw$y_val, 0.99, na.rm = TRUE)
  pad   <- (q_hi - q_lo) * 0.06

  ggplot() +
    geom_point(data = df_raw,
               aes(x = week, y = y_val, color = sport),
               alpha = 0.15, size = 1.4, stroke = 0, shape = 16) +
    geom_ribbon(data = pred_df,
                aes(x = week, ymin = conf.low, ymax = conf.high, fill = sport),
                alpha = 0.25, color = NA) +
    geom_line(data      = pred_df,
              aes(x = week, y = predicted, color = sport),
              linewidth = 1.75) +
    coord_cartesian(ylim = c(q_lo - pad, q_hi + pad)) +
    scale_color_manual(values = SPORT_COLS,  name = "Sport") +
    scale_fill_manual(values  = SPORT_FILLS, name = "Sport") +
    scale_x_continuous(name   = "Monitoring Week",
                       breaks = seq(0, 35, by = 5),
                       expand = expansion(add = c(0.2, 0.8))) +
    scale_y_continuous(labels = scales::comma,
                       expand = expansion(mult = c(0, 0))) +
    labs(title = var_label, y = NULL) +
    theme_P(bs) +
    theme(plot.title   = element_text(face = "bold", size = bs - 0.5,
                                      hjust = 0.5, margin = margin(b = 4)),
          axis.title.x = element_text(face = "bold", size = bs - 1,
                                      margin = margin(t = 8)),
          legend.position = "none")
}

# ---- 5b. DSI boxplot ----
cat("--- 5. Building panels ---\n")
cat("  Panel A: DSI\n")

dsi_data <- master |> filter(!is.na(dsi))

# FV% labels — fixed y ABOVE the data range (scale max ~1.6 + expansion)
# avoids overlap with top whisker / jitter regardless of sport
dsi_pct <- dsi_data |>
  group_by(sport) |>
  summarise(pct = round(mean(dsi < 1.0) * 100, 1),
            .groups = "drop") |>
  mutate(
    lbl = paste0(pct, "% FV Deficit"),
    y   = 1.51     # fixed: above the data range, below the upper y-limit
  )

p_A <- ggplot(dsi_data, aes(x = sport, y = dsi, fill = sport, color = sport)) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 1.0,
           fill = "#FFF5E0", alpha = 0.65) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 1.0,  ymax = Inf,
           fill = "#E5F0F8", alpha = 0.65) +
  geom_boxplot(alpha = 0.60, outlier.shape = NA, width = 0.44,
               color = "grey20", linewidth = 0.65) +
  geom_jitter(width = 0.12, alpha = 0.50, size = 2.4, stroke = 0) +
  geom_hline(yintercept = 1.0, linetype = "dashed",
             color = COL_REF, linewidth = 1.1) +
  # FV% label: white-fill label box, anchored above the data range
  geom_label(data = dsi_pct,
             aes(x = sport, y = y, label = lbl, color = sport),
             inherit.aes    = FALSE,
             size           = 5.2,
             fontface       = "bold",
             fill           = "white",
             label.size     = 0.3,
             label.padding  = unit(0.25, "lines")) +
  # Zone annotations — right side of panel (x > 2), away from both boxes
  annotate("text", x = 2.54, y = 0.53,
           label    = "Force-Velocity Deficit Zone",
           hjust    = 1, size = 4.5, color = "#9E6000", fontface = "italic") +
  annotate("text", x = 2.54, y = 1.46,
           label    = "Adequate Strength Zone",
           hjust    = 1, size = 4.5, color = "#004F8E", fontface = "italic") +
  scale_fill_manual(values  = SPORT_FILLS, name = "Sport") +
  scale_color_manual(values = SPORT_COLS,  name = "Sport") +
  scale_y_continuous(breaks = seq(0.4, 1.6, 0.1),
                     expand = expansion(mult = 0.09)) +
  labs(subtitle = "Dynamic Strength Index  |  Dashed line = DSI 1.0",
       x = NULL, y = "Dynamic Strength Index (DSI)") +
  theme_P(18) +
  theme(axis.text.x = element_text(size = 19, face = "bold"))

# ---- 5c. Jump Height boxplot ----
cat("  Panel B: Jump Height\n")

jh_y_top <- quantile(master$jump_height_cm, 0.99, na.rm = TRUE) + 2.5

jh_stats <- master |>
  filter(!is.na(jump_height_cm)) |>
  group_by(sport) |>
  summarise(m  = mean(jump_height_cm, na.rm = TRUE),
            sd = sd(jump_height_cm,   na.rm = TRUE),
            .groups = "drop") |>
  mutate(lbl = sprintf("%.1f \u00b1 %.1f cm", m, sd),
         y   = jh_y_top)       # fixed above 99th pctile — clear of all jitter

p_B <- ggplot(master |> filter(!is.na(jump_height_cm)),
              aes(x = sport, y = jump_height_cm, fill = sport, color = sport)) +
  geom_boxplot(alpha = 0.60, outlier.shape = NA, width = 0.44,
               color = "grey20", linewidth = 0.65) +
  geom_jitter(width = 0.12, alpha = 0.50, size = 2.4, stroke = 0) +
  # Mean±SD label: white-fill label box anchored above all data points
  geom_label(data = jh_stats,
             aes(x = sport, y = y, label = lbl, color = sport),
             inherit.aes   = FALSE,
             size          = 5.2,
             fontface      = "bold",
             fill          = "white",
             label.size    = 0.3,
             label.padding = unit(0.25, "lines")) +
  scale_fill_manual(values  = SPORT_FILLS, name = "Sport") +
  scale_color_manual(values = SPORT_COLS,  name = "Sport") +
  scale_y_continuous(expand = expansion(mult = 0.09)) +
  labs(subtitle = "Countermovement Jump Height  |  Mean \u00b1 SD annotated",
       x = NULL, y = "Jump Height (cm)") +
  theme_P(18) +
  theme(axis.text.x = element_text(size = 19, face = "bold"))

# ---- 5d. Scatter: CMJ Impulse × Jump Height — rmcorr ----
cat("  Panel C: CMJ Impulse vs Jump Height (rmcorr)\n")

imp_jh <- master |>
  filter(!is.na(cmj_impulse_ns), !is.na(jump_height_cm)) |>
  mutate(athlete_fct = factor(athlete_id))   # rmcorr requires factor/character

# rmcorr per sport (≥2 obs per athlete enforced by the grouping structure)
rmc_C <- lapply(c("Basketball", "Rugby"), function(s) {
  sub <- imp_jh |> filter(sport == s) |>
    group_by(athlete_fct) |> filter(n() >= 2) |> ungroup()
  rmcorr::rmcorr(participant = athlete_fct,
                 measure1    = cmj_impulse_ns,
                 measure2    = jump_height_cm,
                 dataset     = sub)
})
names(rmc_C) <- c("Basketball", "Rugby")

cat(sprintf("    Basketball: r_rm = %.3f  (df=%d, p=%s)\n",
            rmc_C$Basketball$r, rmc_C$Basketball$df,
            ifelse(rmc_C$Basketball$p < .001, "<.001",
                   sprintf("%.3f", rmc_C$Basketball$p))))
cat(sprintf("    Rugby:      r_rm = %.3f  (df=%d, p=%s)\n",
            rmc_C$Rugby$r, rmc_C$Rugby$df,
            ifelse(rmc_C$Rugby$p < .001, "<.001",
                   sprintf("%.3f", rmc_C$Rugby$p))))

# Label positions: top-right sparse area (high impulse, high JH)
r_imp_jh <- tibble(
  sport = factor(c("Basketball", "Rugby"), levels = c("Basketball", "Rugby")),
  r_rm  = c(rmc_C$Basketball$r, rmc_C$Rugby$r),
  p_rm  = c(rmc_C$Basketball$p, rmc_C$Rugby$p)
) |>
  mutate(
    lbl = paste0("r_rm = ", sprintf("%.2f", r_rm),
                 ifelse(p_rm < .001, "***",
                        ifelse(p_rm < .01, "**",
                               ifelse(p_rm < .05, "*", "")))),
    # Top-right corner: high impulse × high JH (sparse for Rugby athletes)
    x   = quantile(imp_jh$cmj_impulse_ns, 0.93, na.rm = TRUE),
    y   = quantile(imp_jh$jump_height_cm, 0.96, na.rm = TRUE) -
            (row_number() - 1L) * 3.5   # stagger vertically
  )

p_C <- ggplot(imp_jh, aes(x = cmj_impulse_ns, y = jump_height_cm,
                           color = sport, fill = sport)) +
  geom_point(alpha = 0.45, size = 2.2, stroke = 0) +
  geom_smooth(method = "lm", formula = y ~ x,
              se = TRUE, alpha = 0.15, linewidth = 1.3) +
  geom_label(data = r_imp_jh,
             aes(x = x, y = y, label = lbl, color = sport),
             inherit.aes   = FALSE,
             size          = 5.2,
             fontface      = "bold",
             fill          = "white",
             label.size    = 0.3,
             label.padding = unit(0.25, "lines"),
             hjust         = 1) +
  scale_color_manual(values = SPORT_COLS,  name = "Sport") +
  scale_fill_manual(values  = SPORT_FILLS, name = "Sport") +
  scale_x_continuous(labels = comma, expand = expansion(mult = 0.05)) +
  scale_y_continuous(expand = expansion(mult = 0.07)) +
  labs(subtitle = paste0("CMJ Propulsive Impulse \u2192 Jump Height  |  ",
                         "r_rm = Repeated Measures Correlation (Bakdash & Marusich, 2017)"),
       x = "CMJ Propulsive Impulse (N\u00b7s)",
       y = "Jump Height (cm)") +
  theme_P(18)

# ---- 5e. Scatter: CMJ Peak Force × IMTP Peak Force — rmcorr ----
cat("  Panel D: CMJ Force vs IMTP Force (rmcorr)\n")

force_scatter <- master |>
  filter(!is.na(cmj_peak_force_n), !is.na(imtp_peak_force_n)) |>
  mutate(athlete_fct = factor(athlete_id))

# rmcorr per sport
rmc_D <- lapply(c("Basketball", "Rugby"), function(s) {
  sub <- force_scatter |> filter(sport == s) |>
    group_by(athlete_fct) |> filter(n() >= 2) |> ungroup()
  rmcorr::rmcorr(participant = athlete_fct,
                 measure1    = imtp_peak_force_n,
                 measure2    = cmj_peak_force_n,
                 dataset     = sub)
})
names(rmc_D) <- c("Basketball", "Rugby")

cat(sprintf("    Basketball: r_rm = %.3f  (df=%d, p=%s)\n",
            rmc_D$Basketball$r, rmc_D$Basketball$df,
            ifelse(rmc_D$Basketball$p < .001, "<.001",
                   sprintf("%.3f", rmc_D$Basketball$p))))
cat(sprintf("    Rugby:      r_rm = %.3f  (df=%d, p=%s)\n",
            rmc_D$Rugby$r, rmc_D$Rugby$df,
            ifelse(rmc_D$Rugby$p < .001, "<.001",
                   sprintf("%.3f", rmc_D$Rugby$p))))

r_forces <- tibble(
  sport = factor(c("Basketball", "Rugby"), levels = c("Basketball", "Rugby")),
  r_rm  = c(rmc_D$Basketball$r, rmc_D$Rugby$r),
  p_rm  = c(rmc_D$Basketball$p, rmc_D$Rugby$p)
) |>
  mutate(
    lbl = paste0("r_rm = ", sprintf("%.2f", r_rm),
                 ifelse(p_rm < .001, "***",
                        ifelse(p_rm < .01, "**",
                               ifelse(p_rm < .05, "*", "")))),
    # Bottom-right corner: high IMTP × low CMJ — sparse area in scatter
    x = quantile(force_scatter$imtp_peak_force_n, 0.93, na.rm = TRUE),
    y = quantile(force_scatter$cmj_peak_force_n,  0.06, na.rm = TRUE) +
          (row_number() - 1L) * 130
  )

p_D <- ggplot(force_scatter,
              aes(x = imtp_peak_force_n, y = cmj_peak_force_n,
                  color = sport, fill = sport)) +
  geom_point(alpha = 0.45, size = 2.2, stroke = 0) +
  geom_smooth(method = "lm", formula = y ~ x,
              se = TRUE, alpha = 0.15, linewidth = 1.3) +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted",
              color = "grey55", linewidth = 0.9) +
  # r_rm labels: bottom-right sparse corner, white-fill box
  geom_label(data = r_forces,
             aes(x = x, y = y, label = lbl, color = sport),
             inherit.aes   = FALSE,
             size          = 5.5,
             fontface      = "bold",
             fill          = "white",
             label.size    = 0.3,
             label.padding = unit(0.25, "lines"),
             hjust         = 1) +
  scale_color_manual(values = SPORT_COLS,  name = "Sport") +
  scale_fill_manual(values  = SPORT_FILLS, name = "Sport") +
  scale_x_continuous(labels = comma, expand = expansion(mult = 0.05)) +
  scale_y_continuous(labels = comma, expand = expansion(mult = 0.05)) +
  labs(subtitle = paste0("CMJ vs. IMTP Peak Force  |  Dotted line = identity  |  ",
                         "r_rm = Repeated Measures Correlation"),
       x = "IMTP Peak Force (N)",
       y = "CMJ Peak Force (N)") +
  theme_P(18)

# ---- 5f. LMM trajectory panels ----
cat("  Panels E1–E5: LMM trajectories\n\n")

p_e1 <- make_lmm_panel(pred_dsi,  "dsi",              "Dynamic Strength Index (DSI)")
p_e2 <- make_lmm_panel(pred_jh,   "jump_height_cm",   "Jump Height (cm)")
p_e3 <- make_lmm_panel(pred_cmj,  "cmj_peak_force_n", "CMJ Peak Force (N)")
p_e4 <- make_lmm_panel(pred_imtp, "imtp_peak_force_n","IMTP Peak Force (N)")
p_e5 <- make_lmm_panel(pred_pp,   "peak_power_w",     "Peak Power (W)")

# =============================================================================
# 6. FIGURE 1 — Baseline & Mechanisms (24 × 11 in)
#
#    Row 1 (5.5 in): [A: DSI Box | B: Jump Height Box]
#    Row 2 (5.5 in): [C: Impulse × JH | D: CMJ × IMTP]
#    Tags: A B C D  |  Single bottom legend
# =============================================================================
cat("--- 6. Assembling Figure 1 (Baseline & Mechanisms) ---\n")

# Shared legend theme for Fig 1 — applied with & to all panels
FIG1_LEGEND <- theme(
  legend.position    = "bottom",
  legend.direction   = "horizontal",
  legend.title       = element_text(face = "bold", size = 20),
  legend.text        = element_text(size = 19),
  legend.key.width   = unit(2.0, "cm"),
  legend.key.height  = unit(0.60, "cm"),
  legend.margin      = margin(t = 5, b = 3),
  legend.box.spacing = unit(0, "pt")
)

fig1 <- (p_A | p_B) /
         (p_C | p_D)  +
  plot_layout(
    guides  = "collect",
    heights = c(5.5, 5.5)   # 11 in total
  ) +
  plot_annotation(
    tag_levels = "A",
    theme = theme(
      plot.tag = element_text(size = 28, face = "bold", color = "grey10",
                              margin = margin(r = 3, b = 3))
    )
  ) &
  FIG1_LEGEND

path_fig1 <- file.path(FIG_DIR, "Hero_Fig1_Baseline_Mechanisms.png")

cat("  Rendering Hero_Fig1_Baseline_Mechanisms.png (24 × 11 in @ 600 dpi)...\n")
t1 <- Sys.time()
ggsave(path_fig1, plot = fig1,
       width = 24, height = 11, dpi = 600, bg = "white", units = "in")
e1   <- round(as.numeric(difftime(Sys.time(), t1, units = "secs")))
mb1  <- round(file.info(path_fig1)$size / 1024^2, 1)
cat(sprintf("  \u2713 Saved: %s  (%.1f MB, %d s)\n\n", basename(path_fig1), mb1, e1))

# =============================================================================
# 7. FIGURE 2 — LMM Predicted Trajectories (24 × 7 in)
#
#    Single row: E1 DSI | E2 Jump Height | E3 CMJ PF | E4 IMTP PF | E5 Peak Power
#    Lines = LMM fixed-effect predicted marginal means (Sport × Week)
#    Ribbon = 95% CI on fixed-effect predictions
#    Y window = coord_cartesian on 1st–99th pctile of raw data ± 6% pad
#    Single bottom legend  |  No panel tags (variable titles serve that role)
# =============================================================================
cat("--- 7. Assembling Figure 2 (LMM Trajectories) ---\n")

# Shared legend theme for Fig 2 — slightly larger for the banner format
FIG2_LEGEND <- theme(
  legend.position    = "bottom",
  legend.direction   = "horizontal",
  legend.title       = element_text(face = "bold", size = 19),
  legend.text        = element_text(size = 18),
  legend.key.width   = unit(2.2, "cm"),
  legend.key.height  = unit(0.60, "cm"),
  legend.margin      = margin(t = 4, b = 2),
  legend.box.spacing = unit(0, "pt")
)

fig2 <- wrap_plots(list(p_e1, p_e2, p_e3, p_e4, p_e5), nrow = 1) +
  plot_layout(guides = "collect") +
  plot_annotation(
    subtitle = paste0(
      "Linear Mixed Model Predicted Marginal Means \u00b1 95% CI  \u2502  ",
      "Points = individual athlete-session observations  \u2502  ",
      "Model: outcome ~ Sport \u00d7 Week + (1 | Athlete)  \u2502  ",
      "Basketball: weeks 1\u201333  \u2502  Rugby: weeks 1\u201331"
    ),
    theme = theme(
      plot.subtitle = element_text(
        size   = 13, color  = "grey20", face   = "italic",
        hjust  = 0.5, margin = margin(b = 3, t = 4)
      )
    )
  ) &
  FIG2_LEGEND

path_fig2 <- file.path(FIG_DIR, "Hero_Fig2_LMM_Trajectories.png")

cat("  Rendering Hero_Fig2_LMM_Trajectories.png (24 × 7 in @ 600 dpi)...\n")
t2 <- Sys.time()
ggsave(path_fig2, plot = fig2,
       width = 24, height = 7, dpi = 600, bg = "white", units = "in")
e2  <- round(as.numeric(difftime(Sys.time(), t2, units = "secs")))
mb2 <- round(file.info(path_fig2)$size / 1024^2, 1)
cat(sprintf("  \u2713 Saved: %s  (%.1f MB, %d s)\n\n", basename(path_fig2), mb2, e2))

# =============================================================================
# 8. FINAL VERIFICATION REPORT
# =============================================================================
cat("\u2554", strrep("\u2550", 62), "\u2557\n", sep = "")
cat("\u2551  Output Verification                                          \u2551\n")
cat("\u255a", strrep("\u2550", 62), "\u255d\n\n", sep = "")

check_file <- function(path, w, h) {
  info <- file.info(path)
  px_w <- round(w * 600)
  px_h <- round(h * 600)
  cat(sprintf("  \u2713 %-44s\n",    basename(path)))
  cat(sprintf("    Size     : %.1f MB\n",  info$size / 1024^2))
  cat(sprintf("    Canvas   : %d \u00d7 %d in  \u2192 %s \u00d7 %s px @ 600 dpi\n",
              w, h, format(px_w, big.mark = ","),
              format(px_h, big.mark = ",")))
  cat(sprintf("    Modified : %s\n\n",    format(info$mtime, "%Y-%m-%d %H:%M:%S")))
}

check_file(path_fig1, 24, 11)
check_file(path_fig2, 24,  7)

cat("LMM ICCs (between-athlete variance proportion):\n")
report_icc <- function(label, fit) {
  vc  <- as.data.frame(VarCorr(fit))
  icc <- vc$sdcor[vc$grp == "athlete_id"]^2 /
         (vc$sdcor[vc$grp == "athlete_id"]^2 + vc$sdcor[vc$grp == "Residual"]^2)
  cat(sprintf("  %-26s  ICC = %.3f\n", label, icc))
}
report_icc("DSI",              m_dsi)
report_icc("Jump Height (cm)", m_jh)
report_icc("CMJ Peak Force",   m_cmj)
report_icc("IMTP Peak Force",  m_imtp)
report_icc("Peak Power",       m_pp)

cat("\n=== Done ===\n")
