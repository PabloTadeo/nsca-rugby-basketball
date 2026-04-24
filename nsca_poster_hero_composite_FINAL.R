# =============================================================================
# NSCA Poster — Hero Composite FINAL
# Neuromuscular Performance: Basketball vs. Rugby
#
# LAYOUT (biomechanical story, 3 rows):
#   Row 1 — Baseline Distributions : DSI boxplot | Jump Height boxplot
#   Row 2 — Mechanical Relationships: CMJ Impulse × JH scatter | CMJ × IMTP scatter
#   Row 3 — LMM Predicted Trajectories: 5 panels (Sport × Week fixed-effect predictions)
#
# Row 3 uses LINEAR MIXED MODEL predicted marginal means (NOT LOESS).
# Models: outcome ~ Sport * Week + (1|Athlete)
# Predictions via ggeffects::ggpredict() — population-average fixed effects.
#
# Authors: Pablo Tadeo Ríos-Gallardo, PhD & Samuel Montalvo, PhD
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(lubridate)
  library(lme4)
  library(lmerTest)   # Satterthwaite p-values; also makes ggeffects use them
  library(ggeffects)  # ggpredict() for LMM marginal predictions
  library(patchwork)
  library(scales)
  library(gt)
  library(gtExtras)
})

BASE_DIR <- "/Volumes/ADATA HD710 PRO/fod_rugby"
RAW_DIR  <- file.path(BASE_DIR, "Raw_data")
FIG_DIR  <- file.path(BASE_DIR, "Figures")
TBL_DIR  <- file.path(BASE_DIR, "Tables")

cat("=== NSCA Hero Composite FINAL (LMM Trajectories) ===\n\n")

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

# ---- Basketball CMJ ----
bas_cmj <- read_utf8("CMJ_Basket - CMJ.csv") |>
  rename(
    athlete_name     = Name,
    athlete_id       = ExternalId,
    week             = Week,
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
    sport = "Basketball",
    week  = as.integer(num(week)),
    across(c(bw_kg, jump_height_cm, peak_power_w, rsi_mod,
             cmj_peak_force_n, cmj_impulse_ns), num),
    jump_height_cm   = ifelse(jump_height_cm < 0 | jump_height_cm > 120, NA, jump_height_cm),
    peak_power_w     = ifelse(peak_power_w   <= 0, NA, peak_power_w),
    rsi_mod          = ifelse(rsi_mod < 0 | rsi_mod > 5, NA, rsi_mod),
    cmj_peak_force_n = ifelse(cmj_peak_force_n <= 0, NA, cmj_peak_force_n),
    cmj_impulse_ns   = ifelse(cmj_impulse_ns   <= 0, NA, cmj_impulse_ns)
  ) |>
  filter(!is.na(athlete_id), !is.na(week))

# ---- Basketball IMTP ----
bas_imtp <- read_utf8("IMPT_Basket - Hoja 1.csv") |>
  rename(
    athlete_id             = ExternalId,
    week                   = Week,
    imtp_peak_force_n      = `Peak Vertical Force [N]`,
    imtp_impulse_100ms_ns  = `Absolute Impulse - 100ms [N s]`
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
  rename(athlete_id = externalId,
         weight_kg  = `peso (kg)`,
         height_cm  = `altura(cm)`) |>
  mutate(weight_kg = num(weight_kg), height_cm = num(height_cm)) |>
  select(athlete_id, weight_kg, height_cm)

# ---- Rugby CMJ ----
rug_cmj <- read_utf8("CMJ_Rugby - Hoja 1.csv") |>
  rename(
    athlete_name     = Name,
    athlete_id       = ExternalId,
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
    sport = "Rugby",
    week  = as.integer(num(week)),
    across(c(bw_kg, jump_height_cm, peak_power_w, rsi_mod,
             cmj_peak_force_n, cmj_impulse_ns), num),
    jump_height_cm   = ifelse(jump_height_cm < 0 | jump_height_cm > 120, NA, jump_height_cm),
    peak_power_w     = ifelse(peak_power_w   <= 0, NA, peak_power_w),
    rsi_mod          = ifelse(rsi_mod < 0 | rsi_mod > 5, NA, rsi_mod),
    cmj_peak_force_n = ifelse(cmj_peak_force_n <= 0, NA, cmj_peak_force_n),
    cmj_impulse_ns   = ifelse(cmj_impulse_ns   <= 0, NA, cmj_impulse_ns)
  ) |>
  filter(!is.na(athlete_id), !is.na(week))

# ---- Rugby IMTP ----
rug_imtp <- read_utf8("IMPT_Rugby - Hoja 1.csv") |>
  rename(
    athlete_id            = ExternalId,
    imtp_peak_force_n     = `Peak Vertical Force [N]`,
    imtp_impulse_100ms_ns = `Absolute Impulse - 100ms [N s]`
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
  df |>
    group_by(athlete_id, athlete_name, sport, week) |>
    summarise(
      across(c(bw_kg, jump_height_cm, peak_power_w, rsi_mod,
               cmj_peak_force_n, cmj_impulse_ns),
             ~ mean(.x, na.rm = TRUE)),
      .groups = "drop"
    ) |>
    mutate(across(where(is.numeric), ~ ifelse(is.nan(.x), NA_real_, .x)))

agg_imtp <- function(df)
  df |>
    group_by(athlete_id, week) |>
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
  mutate(
    dsi       = cmj_peak_force_n / imtp_peak_force_n,
    dsi       = ifelse(dsi < 0.3 | dsi > 2.0, NA_real_, dsi),
    height_cm = NA_real_
  )

rug_perf <- left_join(agg_cmj(rug_cmj), agg_imtp(rug_imtp),
                      by = c("athlete_id", "week")) |>
  left_join(rug_demo, by = "athlete_id") |>
  mutate(
    bw_kg = coalesce(weight_kg, bw_kg),
    dsi   = cmj_peak_force_n / imtp_peak_force_n,
    dsi   = ifelse(dsi < 0.3 | dsi > 2.0, NA_real_, dsi)
  ) |>
  select(-weight_kg)

master <- bind_rows(bas_perf, rug_perf) |>
  mutate(sport = factor(sport, levels = c("Basketball", "Rugby")))

cat(sprintf(
  "Master: %d rows | Basketball n=%d (%d obs) | Rugby n=%d (%d obs)\n\n",
  nrow(master),
  n_distinct(master$athlete_id[master$sport == "Basketball"]),
  sum(master$sport == "Basketball"),
  n_distinct(master$athlete_id[master$sport == "Rugby"]),
  sum(master$sport == "Rugby")
))

# =============================================================================
# 2. DESIGN TOKENS
# =============================================================================
SPORT_COLS  <- c(Basketball = "#0072B2", Rugby = "#E69F00")
SPORT_FILLS <- c(Basketball = "#0072B2", Rugby = "#E69F00")
COL_REF     <- "#D55E00"

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
      plot.caption       = element_text(size = bs - 6, color = "grey45", hjust = 0),
      plot.margin        = margin(10, 14, 6, 14)
    )
}

# =============================================================================
# 3. FIT LINEAR MIXED MODELS (5 outcomes)
# =============================================================================
cat("--- 3. Fitting LMMs (outcome ~ Sport * Week + (1|Athlete)) ---\n")

LMM_CTRL <- lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))

# Helper: fit with convergence check and console report
fit_lmm_silent <- function(outcome_col, label, data, scale_factor = 1) {
  df <- data |>
    filter(!is.na(.data[[outcome_col]]), !is.na(week), !is.na(sport)) |>
    mutate(y_fit = .data[[outcome_col]] / scale_factor)

  fit <- lmer(y_fit ~ sport * week + (1 | athlete_id),
              data = df, REML = TRUE, control = LMM_CTRL)

  conv <- tryCatch(fit@optinfo$conv$lme4$messages, error = function(e) NULL)
  ok   <- is.null(conv) || length(conv) == 0

  # Fallback to Nelder_Mead on convergence warning
  if (!ok) {
    fit <- lmer(y_fit ~ sport * week + (1 | athlete_id),
                data = df, REML = TRUE,
                control = lmerControl(optimizer  = "Nelder_Mead",
                                      optCtrl    = list(maxfun = 2e5)))
    conv <- tryCatch(fit@optinfo$conv$lme4$messages, error = function(e) NULL)
    ok   <- is.null(conv) || length(conv) == 0
  }

  vc      <- as.data.frame(VarCorr(fit))
  icc_val <- vc$sdcor[vc$grp == "athlete_id"]^2 /
             (vc$sdcor[vc$grp == "athlete_id"]^2 + vc$sdcor[vc$grp == "Residual"]^2)

  cat(sprintf("  %s %-28s  n=%d  ICC=%.3f%s\n",
              ifelse(ok, "\u2713", "\u26A0"),
              label,
              nrow(df),
              icc_val,
              ifelse(scale_factor != 1,
                     sprintf("  [fitted in scaled units \u00f7%g]", scale_factor), "")))
  fit
}

m_dsi  <- fit_lmm_silent("dsi",              "DSI",             master, 1)
m_jh   <- fit_lmm_silent("jump_height_cm",   "Jump Height (cm)",master, 1)
m_cmj  <- fit_lmm_silent("cmj_peak_force_n", "CMJ PF (N)",      master, 1000)
m_imtp <- fit_lmm_silent("imtp_peak_force_n","IMTP PF (N)",     master, 1000)
m_pp   <- fit_lmm_silent("peak_power_w",     "Peak Power (W)",  master, 1000)

cat("\n")

# =============================================================================
# 4. LMM PREDICTION HELPER
#    Returns a data frame ready for geom_ribbon + geom_line
#    Predictions are CLIPPED to the observed week range per sport
#    (no extrapolation beyond each sport's monitoring period)
# =============================================================================
get_predictions <- function(fit, var_col, scale_factor = 1) {
  # Observed week range per sport for the given variable
  sport_ranges <- master |>
    filter(!is.na(.data[[var_col]]), !is.na(week)) |>
    group_by(sport) |>
    summarise(w_min = min(week), w_max = max(week), .groups = "drop")

  max_wk <- max(sport_ranges$w_max)

  # Population-average predictions from fixed effects only
  # ggpredict(..., type = "fixed") is the default for merMod
  pred_raw <- ggpredict(
    fit,
    terms = list(
      week  = seq(min(sport_ranges$w_min), max_wk),
      sport = c("Basketball", "Rugby")
    )
  )

  as.data.frame(pred_raw) |>
    rename(week = x) |>
    mutate(
      sport     = factor(group, levels = c("Basketball", "Rugby")),
      predicted = predicted * scale_factor,
      conf.low  = conf.low  * scale_factor,
      conf.high = conf.high * scale_factor
    ) |>
    left_join(sport_ranges, by = "sport") |>
    filter(week >= w_min, week <= w_max)    # clip to observed range
}

cat("--- 4. Extracting LMM predictions ---\n")
pred_dsi  <- get_predictions(m_dsi,  "dsi",              1)
pred_jh   <- get_predictions(m_jh,   "jump_height_cm",   1)
pred_cmj  <- get_predictions(m_cmj,  "cmj_peak_force_n", 1000)
pred_imtp <- get_predictions(m_imtp, "imtp_peak_force_n",1000)
pred_pp   <- get_predictions(m_pp,   "peak_power_w",     1000)
cat("  Done.\n\n")

# =============================================================================
# 5. LMM PANEL BUILDER
#    Ghost points (raw obs) + LMM ribbon (95% CI) + LMM line (fixed effects)
#    Y window: coord_cartesian on 1st–99th pctile of raw data ± 6% pad
#    (coord_cartesian clips the VIEW without removing data from the model)
# =============================================================================
make_lmm_panel <- function(pred_df, var_col, var_label) {
  df_raw <- master |>
    filter(!is.na(.data[[var_col]]), !is.na(week)) |>
    mutate(y_val = .data[[var_col]])

  q_lo  <- quantile(df_raw$y_val, 0.01, na.rm = TRUE)
  q_hi  <- quantile(df_raw$y_val, 0.99, na.rm = TRUE)
  pad   <- (q_hi - q_lo) * 0.06
  ylims <- c(q_lo - pad, q_hi + pad)

  ggplot() +
    # ── Layer 1: individual athlete-session observations (ghost) ──────────────
    geom_point(
      data  = df_raw,
      aes(x = week, y = y_val, color = sport),
      alpha = 0.15, size = 1.4, stroke = 0, shape = 16
    ) +
    # ── Layer 2: LMM 95% CI ribbon ────────────────────────────────────────────
    geom_ribbon(
      data  = pred_df,
      aes(x = week, ymin = conf.low, ymax = conf.high, fill = sport),
      alpha = 0.25, color = NA
    ) +
    # ── Layer 3: LMM predicted marginal mean line ─────────────────────────────
    geom_line(
      data      = pred_df,
      aes(x = week, y = predicted, color = sport),
      linewidth = 1.75
    ) +
    # ── Y window ─────────────────────────────────────────────────────────────
    coord_cartesian(ylim = ylims) +
    # ── Scales ───────────────────────────────────────────────────────────────
    scale_color_manual(values = SPORT_COLS,  name = "Sport") +
    scale_fill_manual(values  = SPORT_FILLS, name = "Sport") +
    scale_x_continuous(
      name   = "Monitoring Week",
      breaks = seq(0, 35, by = 5),
      expand = expansion(add = c(0.2, 0.8))
    ) +
    scale_y_continuous(
      labels = scales::comma,
      expand = expansion(mult = c(0, 0))
    ) +
    # ── Labels ───────────────────────────────────────────────────────────────
    labs(title = var_label, y = NULL) +
    theme_P(16) +
    theme(
      plot.title      = element_text(face = "bold", size = 15.5,
                                     hjust = 0.5, margin = margin(b = 4)),
      axis.title.x    = element_text(face = "bold", size = 15,
                                     margin = margin(t = 8)),
      legend.position = "none"        # collected at composite level
    )
}

# =============================================================================
# 6. ROW 1 — PANEL A: DSI Boxplot + Jitter
# =============================================================================
cat("--- 6. Panel A: DSI by Sport ---\n")

dsi_data <- master |> filter(!is.na(dsi))

dsi_pct <- dsi_data |>
  group_by(sport) |>
  summarise(
    pct = round(mean(dsi < 1.0) * 100, 1),
    y   = quantile(dsi, 0.98, na.rm = TRUE) + 0.03,
    .groups = "drop"
  ) |>
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
            inherit.aes = FALSE, size = 5.5, fontface = "bold",
            lineheight  = 0.9) +
  annotate("text", x = 0.52, y = 0.80,
           label = "Force-Velocity\nDeficit Zone",
           hjust = 0, size = 5.0, color = "#9E6000",
           fontface = "italic", lineheight = 0.9) +
  annotate("text", x = 0.52, y = 1.19,
           label = "Adequate\nStrength Zone",
           hjust = 0, size = 5.0, color = "#004F8E",
           fontface = "italic", lineheight = 0.9) +
  scale_fill_manual(values  = SPORT_FILLS, name = "Sport") +
  scale_color_manual(values = SPORT_COLS,  name = "Sport") +
  scale_y_continuous(breaks  = seq(0.4, 1.6, 0.1),
                     expand  = expansion(mult = 0.09)) +
  labs(
    subtitle = "Dynamic Strength Index by Sport  |  Dashed line = DSI 1.0",
    x = NULL,
    y = "Dynamic Strength Index (DSI)"
  ) +
  theme_P(18) +
  theme(axis.text.x = element_text(size = 19, face = "bold"))

# =============================================================================
# 7. ROW 1 — PANEL B: Jump Height Boxplot + Jitter
# =============================================================================
cat("--- 7. Panel B: Jump Height by Sport ---\n")

jh_stats <- master |>
  filter(!is.na(jump_height_cm)) |>
  group_by(sport) |>
  summarise(
    m  = mean(jump_height_cm, na.rm = TRUE),
    sd = sd(jump_height_cm,   na.rm = TRUE),
    y  = quantile(jump_height_cm, 0.97, na.rm = TRUE) + 0.8,
    .groups = "drop"
  ) |>
  mutate(lbl = sprintf("%.1f \u00b1 %.1f cm", m, sd))

p_B <- ggplot(master |> filter(!is.na(jump_height_cm)),
              aes(x = sport, y = jump_height_cm,
                  fill = sport, color = sport)) +
  geom_boxplot(alpha = 0.60, outlier.shape = NA, width = 0.44,
               color = "grey20", linewidth = 0.65) +
  geom_jitter(width = 0.12, alpha = 0.50, size = 2.4, stroke = 0) +
  geom_text(data = jh_stats,
            aes(x = sport, y = y, label = lbl, color = sport),
            inherit.aes = FALSE, size = 5.2, fontface = "bold") +
  scale_fill_manual(values  = SPORT_FILLS, name = "Sport") +
  scale_color_manual(values = SPORT_COLS,  name = "Sport") +
  scale_y_continuous(expand = expansion(mult = 0.09)) +
  labs(
    subtitle = "Countermovement Jump Height by Sport  |  Mean \u00b1 SD annotated",
    x = NULL,
    y = "Jump Height (cm)"
  ) +
  theme_P(18) +
  theme(axis.text.x = element_text(size = 19, face = "bold"))

# =============================================================================
# 8. ROW 2 — PANEL C: CMJ Propulsive Impulse × Jump Height (scatter)
# =============================================================================
cat("--- 8. Panel C: CMJ Impulse vs Jump Height ---\n")

imp_jh <- master |> filter(!is.na(cmj_impulse_ns), !is.na(jump_height_cm))

r_imp_jh <- imp_jh |>
  group_by(sport) |>
  summarise(
    r   = cor(cmj_impulse_ns, jump_height_cm, use = "complete.obs"),
    x   = min(cmj_impulse_ns,  na.rm = TRUE) * 1.02,
    y   = max(jump_height_cm,  na.rm = TRUE) * 0.97,
    lbl = paste0("r = ", round(r, 2)),
    .groups = "drop"
  )

p_C <- ggplot(imp_jh, aes(x = cmj_impulse_ns, y = jump_height_cm,
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
    subtitle = "CMJ Propulsive Impulse \u2192 Jump Height  |  r = Pearson correlation",
    x = "CMJ Propulsive Impulse (N\u00b7s)",
    y = "Jump Height (cm)"
  ) +
  theme_P(18)

# =============================================================================
# 9. ROW 2 — PANEL D: CMJ Peak Force × IMTP Peak Force (scatter)
# =============================================================================
cat("--- 9. Panel D: CMJ Force vs IMTP Force ---\n")

force_scatter <- master |>
  filter(!is.na(cmj_peak_force_n), !is.na(imtp_peak_force_n))

r_forces <- force_scatter |>
  group_by(sport) |>
  summarise(
    r   = cor(imtp_peak_force_n, cmj_peak_force_n, use = "complete.obs"),
    x   = min(imtp_peak_force_n, na.rm = TRUE) * 1.01,
    y   = max(cmj_peak_force_n,  na.rm = TRUE) * 0.97,
    lbl = paste0("r = ", round(r, 2)),
    .groups = "drop"
  )

p_D <- ggplot(force_scatter,
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
    x = "IMTP Peak Force (N)",
    y = "CMJ Peak Force (N)"
  ) +
  theme_P(18)

# =============================================================================
# 10. ROW 3 — LMM TRAJECTORY PANELS (5 panels)
#     Uses LMM predicted marginal means (NOT LOESS)
#     Lines = fixed-effect Sport × Week predictions
#     Ribbon = 95% CI on fixed-effect predictions
# =============================================================================
cat("--- 10. Building LMM trajectory panels ---\n")

p_e1 <- make_lmm_panel(pred_dsi,  "dsi",              "Dynamic Strength Index (DSI)")
p_e2 <- make_lmm_panel(pred_jh,   "jump_height_cm",   "Jump Height (cm)")
p_e3 <- make_lmm_panel(pred_cmj,  "cmj_peak_force_n", "CMJ Peak Force (N)")
p_e4 <- make_lmm_panel(pred_imtp, "imtp_peak_force_n","IMTP Peak Force (N)")
p_e5 <- make_lmm_panel(pred_pp,   "peak_power_w",     "Peak Power (W)")

cat("  LMM panels built.\n\n")

# Assemble Row 3 as a sub-patchwork, then wrap as a single opaque cell
# so the outer composite preserves A/B/C/D/[E] tag sequence.
p_E_inner <- wrap_plots(list(p_e1, p_e2, p_e3, p_e4, p_e5), nrow = 1) +
  plot_annotation(
    subtitle = paste0(
      "Linear Mixed Model Predicted Marginal Means \u00b1 95% CI  |  ",
      "Points = individual athlete-session observations  |  ",
      "Model: outcome ~ Sport \u00d7 Week + (1 | Athlete)"
    ),
    theme = theme(
      plot.subtitle = element_text(
        size   = 13.5, color  = "grey20", face   = "italic",
        hjust  = 0,    margin = margin(b = 4, l = 2, t = 2)
      )
    )
  )

p_E <- wrap_elements(full = p_E_inner)

# =============================================================================
# 11. PATCHWORK COMPOSITE — BIOMECHANICAL STORY
#
#   Row 1 (h=5.5) : [A: DSI Box | B: Jump Height Box]       Baseline
#   Row 2 (h=5.5) : [C: Impulse×JH scatter | D: CMJ×IMTP]  Mechanisms
#   Row 3 (h=7.0) : [E: 5 LMM Trajectory panels]            Longitudinal
#   Total height   : 18.0 in
# =============================================================================
cat("--- 11. Composing patchwork layout ---\n")

composite <- (p_A | p_B) /
             (p_C | p_D) /
             p_E          +
  plot_layout(
    guides  = "collect",
    heights = c(5.5, 5.5, 7.0)   # total = 18.0 in
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
# 12. EXPORT (24 × 18 in @ 600 dpi — 14,400 × 10,800 px)
# =============================================================================
out_path <- file.path(FIG_DIR, "Hero_NSCA_Poster_Snapshot.png")

cat(sprintf("--- 12. Rendering: 24 \u00d7 18 in @ 600 dpi ---\n"))
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

elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
size_mb <- round(file.info(out_path)$size / 1024^2, 1)
cat(sprintf("  \u2713 Saved : %s\n",   basename(out_path)))
cat(sprintf("  Size  : %.1f MB\n",     size_mb))
cat(sprintf("  Time  : %.0f s\n\n",    elapsed))

# =============================================================================
# 13. FINAL SUMMARY
# =============================================================================
cat("\u2554", strrep("\u2550", 64), "\u2557\n", sep = "")
cat("\u2551  Hero Composite FINAL \u2014 Complete                          \u2551\n")
cat("\u255a", strrep("\u2550", 64), "\u255d\n\n", sep = "")

cat("LAYOUT:\n")
cat("  Row 1 (5.5 in) : A = DSI Boxplot  |  B = Jump Height Boxplot\n")
cat("  Row 2 (5.5 in) : C = CMJ Impulse \u00d7 JH Scatter  |  D = CMJ \u00d7 IMTP Scatter\n")
cat("  Row 3 (7.0 in) : E = LMM Trajectories (5 panels)\n\n")

cat("LMM PANELS IN ROW 3:\n")
cat("  e1 Dynamic Strength Index (DSI)\n")
cat("  e2 Jump Height (cm)\n")
cat("  e3 CMJ Peak Force (N)\n")
cat("  e4 IMTP Peak Force (N)\n")
cat("  e5 Peak Power (W)\n\n")

cat("KEY STATISTICS:\n")
cat(sprintf(
  "  DSI Basketball : %.3f \u00b1 %.3f | FV-Deficit: %.1f%%\n",
  mean(dsi_data$dsi[dsi_data$sport == "Basketball"], na.rm = TRUE),
  sd(dsi_data$dsi[dsi_data$sport == "Basketball"],   na.rm = TRUE),
  mean(dsi_data$dsi[dsi_data$sport == "Basketball"] < 1.0) * 100
))
cat(sprintf(
  "  DSI Rugby      : %.3f \u00b1 %.3f | FV-Deficit: %.1f%%\n",
  mean(dsi_data$dsi[dsi_data$sport == "Rugby"], na.rm = TRUE),
  sd(dsi_data$dsi[dsi_data$sport == "Rugby"],   na.rm = TRUE),
  mean(dsi_data$dsi[dsi_data$sport == "Rugby"] < 1.0) * 100
))

cat(sprintf("\nOutput: %s (%.1f MB, 24\u00d718 in @ 600 dpi)\n",
            basename(out_path), size_mb))
cat("=== Done ===\n")
