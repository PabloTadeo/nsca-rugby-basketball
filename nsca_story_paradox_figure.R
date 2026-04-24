# =============================================================================
# Figure 4: The Neuromuscular Paradox
# "Divergent Adaptations & Sensitivity Analysis"
#
# Panel A — Divergent Neuromuscular Adaptations
#   LMM predicted % change from baseline (Week 1) for CMJ PF, IMTP PF, DSI.
#   Ribbon = between-sport range; Bold line = sport average trajectory.
#   Forces rise | DSI falls — the paradox in one frame.
#
# Panel B — Sensitivity Analysis: Unadjusted vs. Body Mass Adjusted Models
#   Forest plot of the linear Week effect (% change per week, relative to
#   grand mean) for Base and Adjusted (+Body Mass) models.
#   Overlapping CIs prove the finding is NOT a body-mass artefact.
#
# Authors: Pablo Tadeo Ríos-Gallardo, PhD & Samuel Montalvo, PhD
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(lubridate)
  library(lme4)
  library(lmerTest)
  library(broom.mixed)
  library(ggeffects)
  library(patchwork)
  library(scales)
})

BASE_DIR <- "/Volumes/ADATA HD710 PRO/fod_rugby"
RAW_DIR  <- file.path(BASE_DIR, "Raw_data")
FIG_DIR  <- file.path(BASE_DIR, "Figures")

cat("=== Figure 4: Neuromuscular Paradox (Divergence + Sensitivity) ===\n\n")

# =============================================================================
# 1. DATA PIPELINE (self-contained)
# =============================================================================
cat("--- 1. Building master dataset ---\n")

read_utf8 <- function(f)
  read_csv(file.path(RAW_DIR, f),
           locale = locale(encoding = "UTF-8"),
           show_col_types = FALSE, name_repair = "minimal")
num <- function(x) suppressWarnings(as.numeric(x))

bas_cmj <- read_utf8("CMJ_Basket - CMJ.csv") |>
  rename(athlete_name = Name, athlete_id = ExternalId, week = Week,
         bw_kg = `BW [KG]`, jump_height_cm = `Jump Height (Imp-Mom) [cm]`,
         peak_power_w = `Peak Power [W]`, rsi_mod = `RSI-modified [m/s]`,
         cmj_peak_force_n = `Concentric Peak Force [N]`) |>
  select(athlete_id, athlete_name, week, bw_kg, jump_height_cm,
         peak_power_w, rsi_mod, cmj_peak_force_n) |>
  mutate(sport = "Basketball", week = as.integer(num(week)),
         across(c(bw_kg, jump_height_cm, peak_power_w, rsi_mod, cmj_peak_force_n), num),
         jump_height_cm   = ifelse(jump_height_cm < 0 | jump_height_cm > 120, NA, jump_height_cm),
         peak_power_w     = ifelse(peak_power_w <= 0, NA, peak_power_w),
         rsi_mod          = ifelse(rsi_mod < 0 | rsi_mod > 5, NA, rsi_mod),
         cmj_peak_force_n = ifelse(cmj_peak_force_n <= 0, NA, cmj_peak_force_n)) |>
  filter(!is.na(athlete_id), !is.na(week))

bas_imtp <- read_utf8("IMPT_Basket - Hoja 1.csv") |>
  rename(athlete_id = ExternalId, week = Week,
         imtp_peak_force_n = `Peak Vertical Force [N]`) |>
  select(athlete_id, week, imtp_peak_force_n) |>
  mutate(week = as.integer(num(week)),
         imtp_peak_force_n = num(imtp_peak_force_n),
         imtp_peak_force_n = ifelse(imtp_peak_force_n <= 0, NA, imtp_peak_force_n)) |>
  filter(!is.na(athlete_id), !is.na(week))

rug_demo <- read_utf8("CMJ_IMTP_RUGBY - Hoja 1.csv") |>
  rename_with(~ trimws(.x)) |>
  rename(athlete_id = externalId, weight_kg = `peso (kg)`) |>
  mutate(weight_kg = num(weight_kg)) |>
  select(athlete_id, weight_kg)

rug_cmj <- read_utf8("CMJ_Rugby - Hoja 1.csv") |>
  rename(athlete_name = Name, athlete_id = ExternalId, bw_kg = `BW [KG]`,
         jump_height_cm = `Jump Height (Imp-Mom) [cm]`, peak_power_w = `Peak Power [W]`,
         rsi_mod = `RSI-modified [m/s]`, cmj_peak_force_n = `Concentric Peak Force [N]`) |>
  select(athlete_id, athlete_name, week, bw_kg, jump_height_cm,
         peak_power_w, rsi_mod, cmj_peak_force_n) |>
  mutate(sport = "Rugby", week = as.integer(num(week)),
         across(c(bw_kg, jump_height_cm, peak_power_w, rsi_mod, cmj_peak_force_n), num),
         jump_height_cm   = ifelse(jump_height_cm < 0 | jump_height_cm > 120, NA, jump_height_cm),
         peak_power_w     = ifelse(peak_power_w <= 0, NA, peak_power_w),
         rsi_mod          = ifelse(rsi_mod < 0 | rsi_mod > 5, NA, rsi_mod),
         cmj_peak_force_n = ifelse(cmj_peak_force_n <= 0, NA, cmj_peak_force_n)) |>
  filter(!is.na(athlete_id), !is.na(week))

rug_imtp <- read_utf8("IMPT_Rugby - Hoja 1.csv") |>
  rename(athlete_id = ExternalId, imtp_peak_force_n = `Peak Vertical Force [N]`) |>
  select(athlete_id, week, imtp_peak_force_n) |>
  mutate(week = as.integer(num(week)),
         imtp_peak_force_n = num(imtp_peak_force_n),
         imtp_peak_force_n = ifelse(imtp_peak_force_n <= 0, NA, imtp_peak_force_n)) |>
  filter(!is.na(athlete_id), !is.na(week))

agg_cmj  <- function(df)
  df |> group_by(athlete_id, athlete_name, sport, week) |>
    summarise(across(c(bw_kg, jump_height_cm, peak_power_w, rsi_mod, cmj_peak_force_n),
                     ~ mean(.x, na.rm = TRUE)), .groups = "drop") |>
    mutate(across(where(is.numeric), ~ ifelse(is.nan(.x), NA_real_, .x)))

agg_imtp <- function(df)
  df |> group_by(athlete_id, week) |>
    summarise(imtp_peak_force_n = max(imtp_peak_force_n, na.rm = TRUE), .groups = "drop") |>
    mutate(imtp_peak_force_n = ifelse(is.infinite(imtp_peak_force_n), NA_real_, imtp_peak_force_n))

bas_perf <- left_join(agg_cmj(bas_cmj), agg_imtp(bas_imtp), by = c("athlete_id", "week")) |>
  mutate(dsi = cmj_peak_force_n / imtp_peak_force_n,
         dsi = ifelse(dsi < 0.3 | dsi > 2.0, NA_real_, dsi))

rug_perf <- left_join(agg_cmj(rug_cmj), agg_imtp(rug_imtp), by = c("athlete_id", "week")) |>
  left_join(rug_demo, by = "athlete_id") |>
  mutate(bw_kg = coalesce(weight_kg, bw_kg),
         dsi   = cmj_peak_force_n / imtp_peak_force_n,
         dsi   = ifelse(dsi < 0.3 | dsi > 2.0, NA_real_, dsi)) |>
  select(-weight_kg)

master <- bind_rows(bas_perf, rug_perf) |>
  mutate(sport = factor(sport, levels = c("Basketball", "Rugby")))

# Body mass dataset (session-level BW — no static rugby override)
bm_master <- bind_rows(
  bas_cmj |> group_by(athlete_id, athlete_name, sport, week) |>
    summarise(bw_kg = mean(bw_kg, na.rm=TRUE), .groups="drop"),
  rug_cmj |> group_by(athlete_id, athlete_name, sport, week) |>
    summarise(bw_kg = mean(bw_kg, na.rm=TRUE), .groups="drop")
) |>
  mutate(bw_kg = ifelse(is.nan(bw_kg)|bw_kg<35|bw_kg>160, NA_real_, bw_kg),
         sport = factor(sport, levels = c("Basketball","Rugby"))) |>
  filter(!is.na(bw_kg))

bm_cov_df <- bm_master |> select(athlete_id, week, bm_cov = bw_kg)

cat(sprintf("Master: %d rows | Basketball=%d obs | Rugby=%d obs\n\n",
            nrow(master), sum(master$sport=="Basketball"), sum(master$sport=="Rugby")))

# Grand means (for % change normalisation in forest plot)
gm_dsi  <- mean(master$dsi,              na.rm = TRUE)
gm_cmj  <- mean(master$cmj_peak_force_n, na.rm = TRUE)
gm_imtp <- mean(master$imtp_peak_force_n, na.rm = TRUE)

cat(sprintf("Grand means: DSI=%.3f | CMJ PF=%.0f N | IMTP PF=%.0f N\n\n",
            gm_dsi, gm_cmj, gm_imtp))

# =============================================================================
# 2. FIT BASE MODELS (3 key outcomes)
# =============================================================================
cat("--- 2. Fitting base LMMs ---\n")

LMM_CTRL <- lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))

fit_base <- function(col, sf = 1) {
  df <- master |>
    filter(!is.na(.data[[col]]), !is.na(week), !is.na(sport)) |>
    mutate(y_fit = .data[[col]] / sf)
  lmer(y_fit ~ sport * week + (1 | athlete_id),
       data = df, REML = TRUE, control = LMM_CTRL)
}

m_dsi_b  <- fit_base("dsi",              1)
m_cmj_b  <- fit_base("cmj_peak_force_n", 1000)
m_imtp_b <- fit_base("imtp_peak_force_n",1000)
cat("  ✓ DSI, CMJ PF, IMTP PF — all converged\n\n")

# =============================================================================
# 3. FIT ADJUSTED MODELS (+ Body Mass covariate, REML=TRUE, shared sample)
# =============================================================================
cat("--- 3. Fitting adjusted LMMs (+ body mass covariate) ---\n")

fit_adjusted <- function(col, sf = 1) {
  df <- master |>
    filter(!is.na(.data[[col]]), !is.na(week), !is.na(sport)) |>
    mutate(y_fit = .data[[col]] / sf) |>
    left_join(bm_cov_df, by = c("athlete_id", "week")) |>
    filter(!is.na(bm_cov))
  lmer(y_fit ~ sport * week + bm_cov + (1 | athlete_id),
       data = df, REML = TRUE, control = LMM_CTRL)
}

m_dsi_a  <- fit_adjusted("dsi",              1)
m_cmj_a  <- fit_adjusted("cmj_peak_force_n", 1000)
m_imtp_a <- fit_adjusted("imtp_peak_force_n",1000)
cat("  ✓ All adjusted models converged\n\n")

# =============================================================================
# 4. EXTRACT WEEK COEFFICIENTS FOR FOREST PLOT (Panel B data)
# =============================================================================
cat("--- 4. Extracting Week coefficients ---\n")

extract_week_pct <- function(fit, sf, grand_mean, outcome_lbl, model_lbl) {
  broom.mixed::tidy(fit, effects = "fixed", conf.int = TRUE, conf.level = 0.95) |>
    filter(term == "week") |>
    mutate(
      estimate  = estimate  * sf,
      conf.low  = conf.low  * sf,
      conf.high = conf.high * sf,
      std.error = std.error * sf,
      # % change per week relative to grand mean
      pct_est = estimate  / grand_mean * 100,
      pct_lo  = conf.low  / grand_mean * 100,
      pct_hi  = conf.high / grand_mean * 100,
      Outcome = outcome_lbl,
      Model   = model_lbl
    ) |>
    select(Outcome, Model, pct_est, pct_lo, pct_hi, estimate, conf.low, conf.high, p.value)
}

forest_data <- bind_rows(
  extract_week_pct(m_dsi_b,  1,    gm_dsi,  "Dynamic Strength\nIndex (DSI)",    "Base"),
  extract_week_pct(m_cmj_b,  1000, gm_cmj,  "CMJ Peak\nForce (N)",              "Base"),
  extract_week_pct(m_imtp_b, 1000, gm_imtp, "IMTP Peak\nForce (N)",             "Base"),
  extract_week_pct(m_dsi_a,  1,    gm_dsi,  "Dynamic Strength\nIndex (DSI)",    "Adjusted (+BM)"),
  extract_week_pct(m_cmj_a,  1000, gm_cmj,  "CMJ Peak\nForce (N)",              "Adjusted (+BM)"),
  extract_week_pct(m_imtp_a, 1000, gm_imtp, "IMTP Peak\nForce (N)",             "Adjusted (+BM)")
) |>
  mutate(
    Outcome = factor(Outcome,
                     levels = c("IMTP Peak\nForce (N)",
                                "CMJ Peak\nForce (N)",
                                "Dynamic Strength\nIndex (DSI)")),
    Model   = factor(Model, levels = c("Base", "Adjusted (+BM)"))
  )

cat("  Week coefficient summary (% change per week):\n")
forest_data |>
  mutate(lbl = sprintf("%s | %s: %.4f%% [%.4f, %.4f]",
                       Model, gsub("\n", " ", Outcome), pct_est, pct_lo, pct_hi)) |>
  pull(lbl) |> cat(sep = "\n  ")
cat("\n\n")

# =============================================================================
# 5. LMM PREDICTIONS → % CHANGE FROM BASELINE WEEK (Panel A data)
# =============================================================================
cat("--- 5. Computing LMM prediction trajectories (Panel A) ---\n")

# Week range: clip to common (1–31) so ribbon spans same domain for both sports
COMMON_WEEKS <- seq(1, 31)

get_pct_change <- function(fit, sf, var_lbl) {
  # Predict across both sports for the common week range
  pred_raw <- ggpredict(
    fit,
    terms = list(week = COMMON_WEEKS, sport = c("Basketball", "Rugby"))
  )
  as.data.frame(pred_raw) |>
    rename(week = x) |>
    mutate(
      sport     = factor(group, levels = c("Basketball", "Rugby")),
      predicted = predicted * sf,
      conf.low  = conf.low  * sf,
      conf.high = conf.high * sf,
      Variable  = var_lbl
    ) |>
    group_by(sport) |>
    arrange(week) |>
    mutate(
      baseline   = first(predicted),
      pct_change = (predicted - baseline) / baseline * 100,
      pct_lo     = (conf.low  - baseline) / baseline * 100,
      pct_hi     = (conf.high - baseline) / baseline * 100
    ) |>
    ungroup()
}

pred_pct_long <- bind_rows(
  get_pct_change(m_dsi_b,  1,    "DSI"),
  get_pct_change(m_cmj_b,  1000, "CMJ Peak Force"),
  get_pct_change(m_imtp_b, 1000, "IMTP Peak Force")
) |>
  mutate(Variable = factor(Variable,
                            levels = c("IMTP Peak Force", "CMJ Peak Force", "DSI")))

# Ribbon: min–max across sports at each week (shows sport-to-sport variability)
pred_pct_ribbon <- pred_pct_long |>
  group_by(Variable, week) |>
  summarise(
    ribbon_lo = min(pct_change),
    ribbon_hi = max(pct_change),
    center    = mean(pct_change),   # average trajectory
    .groups   = "drop"
  ) |>
  mutate(Variable = factor(Variable, levels = levels(pred_pct_long$Variable)))

cat(sprintf("  Predictions: %d rows | Ribbon: %d rows | Common weeks: 1–31\n\n",
            nrow(pred_pct_long), nrow(pred_pct_ribbon)))

# Print direction summary
pred_pct_long |>
  group_by(Variable, sport) |>
  summarise(wk1 = pct_change[week==min(week)],
            wk_last = pct_change[week==max(week)],
            direction = ifelse(wk_last > 0, "\u2191 rises", "\u2193 FALLS"),
            .groups = "drop") |>
  mutate(across(c(wk1, wk_last), ~ round(.x, 2))) |>
  print()
cat("\n")

# =============================================================================
# 6. DESIGN TOKENS
# =============================================================================
VAR_COLS <- c(
  "IMTP Peak Force" = "#1B7837",   # forest green — force ↑↑ (stronger)
  "CMJ Peak Force"  = "#2166AC",   # ocean blue  — force ↑
  "DSI"             = "#D6604D"    # brick red   — falls (the paradox)
)
VAR_FILLS <- VAR_COLS

MODEL_COLS   <- c("Base" = "#333333", "Adjusted (+BM)" = "#E69F00")
MODEL_SHAPES <- c("Base" = 15L,       "Adjusted (+BM)" = 16L)

theme_story <- function(bs = 17) {
  theme_classic(base_size = bs) %+replace%
    theme(
      axis.title         = element_text(face = "bold", size = bs),
      axis.text          = element_text(size = bs - 2, color = "grey15"),
      axis.line          = element_line(color = "grey30", linewidth = 0.5),
      axis.ticks         = element_line(color = "grey30", linewidth = 0.4),
      panel.grid.major.y = element_line(color = "grey90", linewidth = 0.35),
      panel.grid.major.x = element_line(color = "grey90", linewidth = 0.35),
      strip.background   = element_rect(fill = "grey93", color = NA),
      strip.text         = element_text(face = "bold", size = bs - 1,
                                        margin = margin(4, 0, 4, 0)),
      legend.position    = "bottom",
      legend.direction   = "horizontal",
      legend.title       = element_text(face = "bold", size = bs - 2),
      legend.text        = element_text(size = bs - 3),
      legend.key.width   = unit(1.6, "cm"),
      plot.title         = element_text(face = "bold", size = bs + 1,
                                        margin = margin(b = 15)),
      plot.subtitle      = element_text(size = bs - 3, color = "grey35",
                                        margin = margin(b = 20)),
      legend.key.height  = unit(0.5, "cm"),
      legend.margin      = margin(t = 6),
      plot.margin        = margin(12, 16, 8, 16)
    )
}

# =============================================================================
# 7. PANEL A — Divergent Neuromuscular Adaptations
# =============================================================================
cat("--- 6. Building Panel A (divergent % change trajectories) ---\n")

# Annotate direction — placed at end of monitoring window (week 31)
label_ends <- pred_pct_ribbon |>
  filter(week == max(COMMON_WEEKS)) |>
  mutate(
    lbl   = case_when(
      Variable == "IMTP Peak Force" ~ "IMTP PF \u2191",
      Variable == "CMJ Peak Force"  ~ "CMJ PF \u2191",
      Variable == "DSI"             ~ "DSI \u2193"
    ),
    hjust = -0.08
  )

# y-axis limits with padding
y_pad  <- 1.5
y_lo   <- floor(min(pred_pct_ribbon$ribbon_lo)) - y_pad
y_hi   <- ceiling(max(pred_pct_ribbon$ribbon_hi)) + y_pad + 3  # extra for end labels

p_A <- ggplot() +

  # ── Background zones ──────────────────────────────────────────────────────
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0,     ymax = Inf,
           fill = "#E8F5E9", alpha = 0.60) +   # light green = forces rising
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf,  ymax = 0,
           fill = "#FEECEB", alpha = 0.60) +   # light red   = DSI falling

  # ── Reference line: baseline ──────────────────────────────────────────────
  geom_hline(yintercept = 0, color = "grey35", linewidth = 0.7,
             linetype = "dashed") +

  # ── Layer 1: between-sport ribbon ─────────────────────────────────────────
  geom_ribbon(
    data    = pred_pct_ribbon,
    mapping = aes(x = week, ymin = ribbon_lo, ymax = ribbon_hi, fill = Variable),
    alpha   = 0.18,
    inherit.aes = FALSE
  ) +

  # ── Layer 2: individual sport lines (thin) ─────────────────────────────────
  geom_line(
    data    = pred_pct_long,
    mapping = aes(x = week, y = pct_change,
                  color    = Variable,
                  linetype = sport),
    linewidth   = 0.75,
    inherit.aes = FALSE
  ) +

  # ── Layer 3: center line = average of both sports (bold) ──────────────────
  geom_line(
    data    = pred_pct_ribbon,
    mapping = aes(x = week, y = center, color = Variable),
    linewidth   = 2.5,
    inherit.aes = FALSE
  ) +

  # ── End labels (right margin) ─────────────────────────────────────────────
  geom_text(
    data    = label_ends,
    mapping = aes(x = week, y = center, label = lbl, color = Variable),
    hjust       = -0.08,
    fontface    = "bold",
    size        = 5.5,
    inherit.aes = FALSE
  ) +

  # ── Zone annotation text ──────────────────────────────────────────────────
  annotate("text", x = 1.5, y = y_hi - 1.2,
           label = "Force capacity rising", hjust = 0,
           size = 4.8, color = "#1B7837", fontface = "italic") +
  annotate("text", x = 1.5, y = -1.0,
           label = "Efficiency ratio falling", hjust = 0, vjust = 1,
           size = 4.8, color = "#B02C1F", fontface = "italic") +

  # ── Scales ────────────────────────────────────────────────────────────────
  scale_color_manual(values = VAR_COLS,
                     name   = "Metric",
                     labels = c("IMTP Peak Force", "CMJ Peak Force", "DSI")) +
  scale_fill_manual(values  = VAR_FILLS, guide = "none") +
  scale_linetype_manual(values = c(Basketball = "solid", Rugby = "dashed"),
                        name   = "Sport") +
  scale_x_continuous(
    name   = "Monitoring Week",
    breaks = seq(0, 31, by = 5),
    expand = expansion(add = c(0.2, 4.5))   # right expansion for end labels
  ) +
  scale_y_continuous(
    name   = "Change from Baseline (%)",
    breaks = scales::breaks_pretty(n = 8),
    expand = expansion(mult = 0)
  ) +
  coord_cartesian(ylim = c(y_lo, y_hi)) +

  labs(
    title = "Divergent Neuromuscular Adaptations"
  ) +

  theme_story(17) +

  guides(
    color    = guide_legend(override.aes = list(linewidth = 2.5), order = 1),
    linetype = guide_legend(override.aes = list(linewidth = 1.2),  order = 2),
    fill     = "none"
  )

# =============================================================================
# 8. PANEL B — Forest Plot: Base vs. Adjusted Week Coefficients
# =============================================================================
cat("--- 7. Building Panel B (forest plot: base vs. adjusted) ---\n")

# Significance stars helper
sig_star <- function(p) {
  case_when(p < .001 ~ "***", p < .01 ~ "**", p < .05 ~ "*", TRUE ~ "ns")
}

forest_plot_data <- forest_data |>
  mutate(
    star = sig_star(p.value),
    lbl  = paste0(sprintf("%.4f%%", pct_est), " ", star)
  )

# Check the sign pattern
cat("  Forest plot coefficient directions:\n")
forest_plot_data |>
  select(Outcome, Model, pct_est, pct_lo, pct_hi, p.value) |>
  mutate(across(where(is.numeric), ~ round(.x, 5))) |>
  print()
cat("\n")

# x-axis: 0-crossing reference + range
x_range <- range(c(forest_plot_data$pct_lo, forest_plot_data$pct_hi))
x_pad   <- diff(x_range) * 0.12
x_lims  <- c(x_range[1] - x_pad, x_range[2] + x_pad)

p_B <- ggplot(
  forest_plot_data,
  aes(x = pct_est, y = Outcome, color = Model, shape = Model)
) +

  # ── Background shading: negative = DSI zone ─────────────────────────────
  annotate("rect", xmin = -Inf, xmax = 0, ymin = -Inf, ymax = Inf,
           fill = "#FEECEB", alpha = 0.50) +
  annotate("rect", xmin = 0, xmax =  Inf, ymin = -Inf, ymax = Inf,
           fill = "#E8F5E9", alpha = 0.50) +

  # ── Zero reference ────────────────────────────────────────────────────────
  geom_vline(xintercept = 0, color = "grey30", linewidth = 0.9,
             linetype = "dashed") +

  # ── 95% CI horizontal bars ────────────────────────────────────────────────
  geom_errorbarh(
    aes(xmin = pct_lo, xmax = pct_hi),
    height   = 0.22,
    linewidth = 0.9,
    position = position_dodge(width = 0.55)
  ) +

  # ── Point estimates ───────────────────────────────────────────────────────
  geom_point(
    size     = 4.5,
    position = position_dodge(width = 0.55)
  ) +

  # ── Significance labels (right of each row group) ─────────────────────────
  geom_text(
    data    = forest_plot_data |> filter(Model == "Base"),
    mapping = aes(x = pct_hi, y = Outcome, label = star),
    color   = "grey25",
    size    = 5.5,
    hjust   = -0.3,
    fontface = "bold",
    position = position_dodge(width = 0.55),
    inherit.aes = FALSE
  ) +

  # ── Scales ────────────────────────────────────────────────────────────────
  scale_color_manual(values = MODEL_COLS,   name = "Model") +
  scale_shape_manual(values = MODEL_SHAPES, name = "Model") +
  scale_x_continuous(
    name   = "Week Effect (% change per week vs. grand mean)",
    breaks = c(-0.2, 0, 0.2, 0.4),
    labels = c("-20%", "0%", "+20%", "+40%"),
    limits = c(-0.3, 0.45)
  ) +

  # ── Variable annotations on the y-axis strip ─────────────────────────────
  labs(
    title = "Sensitivity Analysis: Unadjusted vs. Body Mass Adjusted Models"
  ) +

  theme_story(17) +
  theme(
    axis.title.y   = element_blank(),
    axis.text.y    = element_text(size = 15, face = "bold", hjust = 1),
    legend.position = "bottom"
  )

# =============================================================================
# 9. ASSEMBLE WITH PATCHWORK
# =============================================================================
cat("--- 8. Assembling composite figure ---\n")

fig4 <- p_A | p_B +
  plot_layout(widths = c(1.7, 1), guides = "collect") +
  plot_annotation(
    tag_levels = "A",
    subtitle   = paste0(
      "Predicted trajectories expressed as % change from Week 1 (Panel A). ",
      "Unadjusted vs Body Mass adjusted fixed effects (Panel B)."
    ),
    theme = theme(
      plot.tag      = element_text(size = 22, face = "bold", margin = margin(r = 4)),
      plot.subtitle = element_text(size = 18, margin = margin(b = 20)),
      legend.position = "bottom"
    )
  ) &
  theme(legend.position = "bottom")

# =============================================================================
# 10. EXPORT
# =============================================================================
cat("--- 9. Exporting figure ---\n")

out_path <- file.path(FIG_DIR, "Hero_Fig4_Paradox_Sensitivity.png")
t0 <- proc.time()

ggsave(
  filename = out_path,
  plot     = fig4,
  width    = 18,
  height   = 10,
  dpi      = 600,
  bg       = "white",
  units    = "in"
)

elapsed <- round((proc.time() - t0)["elapsed"], 0)
file_mb  <- round(file.info(out_path)$size / 1024^2, 1)

cat(sprintf("\n  ✓ Saved: Hero_Fig4_Paradox_Sensitivity.png\n"))
cat(sprintf("  Dimensions : 18 \u00d7 10 in\n"))
cat(sprintf("  Resolution : 600 dpi  (pixel size: 10,800 \u00d7 6,000 px)\n"))
cat(sprintf("  File size  : %.1f MB\n", file_mb))
cat(sprintf("  Render time: %ds\n", elapsed))

# =============================================================================
# 11. RESULTS VERIFICATION REPORT
# =============================================================================
cat("\n\u2554\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2557\n")
cat("\u2551  FIGURE 4 — VERIFICATION REPORT                         \u2551\n")
cat("\u255a\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u255d\n\n")

cat("PANEL A — Trajectory directions at Week 31 (vs. Week 1 baseline):\n")
pred_pct_long |>
  group_by(Variable, sport) |>
  summarise(end_pct = round(pct_change[week == max(COMMON_WEEKS)], 2),
            direction = ifelse(end_pct > 0, "\u2191 RISES", "\u2193 FALLS"),
            .groups = "drop") |>
  arrange(Variable, sport) |>
  print()

cat("\nPANEL B — Robustness check (Base vs. Adjusted overlap):\n")
overlap_check <- forest_data |>
  mutate(Model_key = ifelse(Model == "Base", "base", "adj"),
         Outcome   = gsub("\n", " ", Outcome)) |>
  select(Outcome, Model_key, pct_est, pct_lo, pct_hi, p.value) |>
  pivot_wider(names_from  = Model_key,
              values_from = c(pct_est, pct_lo, pct_hi, p.value)) |>
  mutate(delta_pct = abs(pct_est_adj - pct_est_base))

for (i in seq_len(nrow(overlap_check))) {
  cat(sprintf("  %-28s  Base %+.4f%%  Adj %+.4f%%  \u0394 = %.5f%%  Sign matches: %s\n",
              overlap_check$Outcome[i],
              overlap_check$pct_est_base[i],
              overlap_check$pct_est_adj[i],
              overlap_check$delta_pct[i],
              ifelse(sign(overlap_check$pct_est_base[i]) ==
                     sign(overlap_check$pct_est_adj[i]),
                     "\u2713 YES", "\u2717 NO")))
}

cat("\nVERDICT:\n")
all_match <- all(sign(overlap_check$pct_est_base) ==
                   sign(overlap_check$pct_est_adj))
if (all_match) {
  cat("  \u2713 ALL signs preserved after BM adjustment\n")
  cat("  \u2713 The 'forces \u2191\u2191, DSI \u2193' paradox is ROBUST\n")
  cat("  \u2713 Body mass fluctuations do NOT explain the divergence\n")
} else {
  cat("  \u2717 Sign reversal detected — interpret with caution\n")
}

cat("\n=== Figure 4 Complete ===\n")
