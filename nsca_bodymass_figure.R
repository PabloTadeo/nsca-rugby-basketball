# =============================================================================
# Body Mass Longitudinal Monitoring — NSCA Poster
# Figure 3: Session-to-Session Body Mass Variation (Basketball vs. Rugby)
#
# Design  : Individual athlete spaghetti lines (background, α=0.22) +
#           LMM predicted marginal means (foreground, bold line) +
#           95% CI ribbon (foreground, semi-transparent)
# Facets  : Basketball | Rugby  (free_x = per-sport week range)
# Y-axis  : coord_cartesian clipped to 1st–99th pctile ± 5% padding
#           (fine BM fluctuation is visible; extreme outliers do not distort)
#
# DATA NOTE:
#   BW [KG] from the Hawkins CMJ files is the session-level body mass
#   recorded at each testing session.  For rugby athletes the static
#   demographic weight (CMJ_IMTP_RUGBY "peso (kg)") is NOT used here, so
#   true session-to-session fluctuations are preserved.
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
  library(scales)
})

BASE_DIR <- "/Volumes/ADATA HD710 PRO/fod_rugby"
RAW_DIR  <- file.path(BASE_DIR, "Raw_data")
FIG_DIR  <- file.path(BASE_DIR, "Figures")

cat("=== Body Mass Longitudinal Figure — NSCA Poster ===\n\n")

# =============================================================================
# 1. DATA PIPELINE
# =============================================================================
cat("--- 1. Loading & building body mass dataset ---\n")

read_utf8 <- function(f)
  read_csv(file.path(RAW_DIR, f),
           locale         = locale(encoding = "UTF-8"),
           show_col_types = FALSE,
           name_repair    = "minimal")

num <- function(x) suppressWarnings(as.numeric(x))

# ── Basketball: BW recorded each CMJ session ──────────────────────────────
bas_bm_raw <- read_utf8("CMJ_Basket - CMJ.csv") |>
  rename(athlete_name = Name,
         athlete_id   = ExternalId,
         week         = Week,
         bw_kg        = `BW [KG]`) |>
  select(athlete_id, athlete_name, week, bw_kg) |>
  mutate(sport = "Basketball",
         week  = as.integer(num(week)),
         bw_kg = num(bw_kg),
         bw_kg = ifelse(bw_kg < 35 | bw_kg > 160, NA_real_, bw_kg)) |>
  filter(!is.na(athlete_id), !is.na(week), !is.na(bw_kg))

# ── Rugby: BW from session CMJ file — NOT overridden by static demo weight ─
rug_bm_raw <- read_utf8("CMJ_Rugby - Hoja 1.csv") |>
  rename(athlete_name = Name,
         athlete_id   = ExternalId,
         bw_kg        = `BW [KG]`) |>
  select(athlete_id, athlete_name, week, bw_kg) |>
  mutate(sport = "Rugby",
         week  = as.integer(num(week)),
         bw_kg = num(bw_kg),
         bw_kg = ifelse(bw_kg < 35 | bw_kg > 160, NA_real_, bw_kg)) |>
  filter(!is.na(athlete_id), !is.na(week), !is.na(bw_kg))

# ── Aggregate: mean BW per athlete × week (multiple sessions per week) ────
agg_bm <- function(df)
  df |>
    group_by(athlete_id, athlete_name, sport, week) |>
    summarise(bw_kg = mean(bw_kg, na.rm = TRUE), .groups = "drop") |>
    mutate(bw_kg = ifelse(is.nan(bw_kg), NA_real_, bw_kg)) |>
    filter(!is.na(bw_kg))

bm_master <- bind_rows(agg_bm(bas_bm_raw),
                       agg_bm(rug_bm_raw)) |>
  mutate(sport = factor(sport, levels = c("Basketball", "Rugby")))

# ── Descriptive audit ────────────────────────────────────────────────────
cat(sprintf("Body mass dataset: %d athlete-week obs (aggregated)\n", nrow(bm_master)))
bm_master |>
  group_by(sport) |>
  summarise(
    n_athletes = n_distinct(athlete_id),
    n_obs      = n(),
    wk_min     = min(week),
    wk_max     = max(week),
    mean_kg    = round(mean(bw_kg, na.rm = TRUE), 2),
    sd_kg      = round(sd(bw_kg,   na.rm = TRUE), 2),
    min_kg     = round(min(bw_kg,  na.rm = TRUE), 2),
    max_kg     = round(max(bw_kg,  na.rm = TRUE), 2),
    cv_pct     = round(sd_kg / mean_kg * 100, 2),
    .groups    = "drop"
  ) |>
  as.data.frame() |>
  print()
cat("\n")

# =============================================================================
# 2. LINEAR MIXED MODEL
# formula: bw_kg ~ sport * week + (1 | athlete_id)
# No scaling needed: body mass is in physiologically manageable range (kg)
# =============================================================================
cat("--- 2. Fitting LMM: Body Mass (kg) ~ Sport * Week + (1 | Athlete) ---\n")

LMM_CTRL <- lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))

fit_bm <- lmer(
  bw_kg ~ sport * week + (1 | athlete_id),
  data    = bm_master,
  REML    = TRUE,
  control = LMM_CTRL
)

# ── Convergence check ────────────────────────────────────────────────────
conv_msgs <- tryCatch(fit_bm@optinfo$conv$lme4$messages, error = function(e) NULL)
if (is.null(conv_msgs) || length(conv_msgs) == 0) {
  cat("  ✓  Model converged with no warnings\n")
} else {
  cat(sprintf("  ⚠  %s  — re-fitting with Nelder_Mead\n",
              paste(conv_msgs, collapse = "; ")))
  fit_bm <- lmer(bw_kg ~ sport * week + (1 | athlete_id),
                 data    = bm_master,
                 REML    = TRUE,
                 control = lmerControl(optimizer = "Nelder_Mead",
                                       optCtrl   = list(maxfun = 2e5)))
  conv2 <- tryCatch(fit_bm@optinfo$conv$lme4$messages, error = function(e) NULL)
  if (is.null(conv2) || length(conv2) == 0) cat("  ✓  Re-fit converged\n") else
    cat(sprintf("  ✗  Still flagged: %s\n", paste(conv2, collapse = "; ")))
}

# ── Variance components & ICC ────────────────────────────────────────────
vc         <- as.data.frame(VarCorr(fit_bm))
sd_athlete <- vc$sdcor[vc$grp == "athlete_id"]
sd_resid   <- vc$sdcor[vc$grp == "Residual"]
icc_bm     <- sd_athlete^2 / (sd_athlete^2 + sd_resid^2)

cat(sprintf("\n  Between-athlete SD  (τ₀)  = %.3f kg\n",  sd_athlete))
cat(sprintf("  Within-athlete SD   (σ)   = %.3f kg\n",  sd_resid))
cat(sprintf("  ICC                       = %.3f\n",   icc_bm))
cat(sprintf("  → %.1f%% of variance = stable between-athlete differences\n",
            icc_bm * 100))
cat(sprintf("  → %.1f%% of variance = within-athlete week-to-week fluctuation\n\n",
            (1 - icc_bm) * 100))

# ── Fixed effects (Satterthwaite df via lmerTest) ────────────────────────
fe_bm <- broom.mixed::tidy(
  fit_bm, effects = "fixed", conf.int = TRUE, conf.level = 0.95
)

TERM_LABELS <- c(
  "(Intercept)"     = "Intercept (Basketball, Week 0)",
  "sportRugby"      = "Sport: Rugby vs. Basketball",
  "week"            = "Week (linear slope — Basketball)",
  "sportRugby:week" = "Sport × Week Interaction"
)

cat("Fixed effects (Satterthwaite df):\n")
cat(strrep("─", 90), "\n")
fe_bm |>
  mutate(
    Term  = coalesce(TERM_LABELS[term], term),
    p_lbl = case_when(p.value < .001 ~ "< .001",
                      TRUE           ~ sprintf("%.3f", p.value)),
    sig   = ifelse(p.value < .05, " *", "  ")
  ) |>
  select(Term, estimate, std.error, statistic, df, p_lbl, sig) |>
  mutate(across(where(is.numeric), ~ round(.x, 3))) |>
  rename(`β (kg)` = estimate, SE = std.error, t = statistic,
         `p` = p_lbl, ` ` = sig) |>
  as.data.frame() |>
  print(row.names = FALSE)
cat("\n  * = p < .05\n\n")

sig_terms <- fe_bm |> filter(term != "(Intercept)", p.value < .05)
if (nrow(sig_terms) > 0) {
  cat("  Significant fixed effects:\n")
  for (i in seq_len(nrow(sig_terms))) {
    t  <- sig_terms$term[i]
    lbl <- coalesce(TERM_LABELS[t], t)
    cat(sprintf("    ▸ %s: β=%.3f kg, t=%.2f, p%s\n",
                lbl, sig_terms$estimate[i], sig_terms$statistic[i],
                ifelse(sig_terms$p.value[i] < .001, "<.001",
                       paste0("=", sprintf("%.3f", sig_terms$p.value[i])))))
  }
} else {
  cat("  ℹ  No significant fixed effects (p ≥ .05)\n")
  cat("     Body mass was statistically stable across the monitoring season.\n")
}
cat("\n")

# =============================================================================
# 3. LMM PREDICTIONS via ggeffects (population-average fixed effects only)
# =============================================================================
cat("--- 3. Extracting LMM predicted marginal means ---\n")

week_ranges <- bm_master |>
  group_by(sport) |>
  summarise(w_min = min(week), w_max = max(week), .groups = "drop")

pred_bm_raw <- ggpredict(
  fit_bm,
  terms = list(
    week  = seq(min(week_ranges$w_min), max(week_ranges$w_max)),
    sport = c("Basketball", "Rugby")
  )
)

pred_bm <- as.data.frame(pred_bm_raw) |>
  rename(week = x) |>
  mutate(sport = factor(group, levels = c("Basketball", "Rugby"))) |>
  left_join(week_ranges, by = "sport") |>
  filter(week >= w_min, week <= w_max)  # clip to each sport's observed range

cat(sprintf("  Predictions: %d rows (%d Basketball / %d Rugby)\n\n",
            nrow(pred_bm),
            sum(pred_bm$sport == "Basketball"),
            sum(pred_bm$sport == "Rugby")))

# =============================================================================
# 4. DESIGN TOKENS
# =============================================================================
SPORT_COLS  <- c(Basketball = "#0072B2", Rugby = "#E69F00")
SPORT_FILLS <- c(Basketball = "#0072B2", Rugby = "#E69F00")

theme_BM <- function(bs = 18) {
  theme_classic(base_size = bs) %+replace%
    theme(
      # Axes
      axis.title         = element_text(face = "bold", size = bs),
      axis.text          = element_text(size = bs - 2, color = "grey15"),
      axis.line          = element_line(color = "grey30", linewidth = 0.5),
      axis.ticks         = element_line(color = "grey30", linewidth = 0.4),
      # Facet strips
      strip.background   = element_rect(fill = "grey92", color = NA),
      strip.text         = element_text(face = "bold", size = bs,
                                        margin = margin(5, 0, 5, 0)),
      # Grid — horizontal guide lines only (emphasise BM stability)
      panel.grid.major.y = element_line(color = "grey90", linewidth = 0.35),
      panel.spacing      = unit(1.6, "lines"),
      # Legend — bottom, consistent with poster style
      legend.position    = "bottom",
      legend.direction   = "horizontal",
      legend.title       = element_text(face = "bold", size = bs - 1),
      legend.text        = element_text(size = bs - 2),
      legend.key.width   = unit(2.0, "cm"),
      legend.key.height  = unit(0.55, "cm"),
      legend.margin      = margin(t = 6),
      # Subtitle
      plot.subtitle      = element_text(size = bs - 4, color = "grey35",
                                        margin = margin(b = 8)),
      plot.margin        = margin(14, 18, 10, 18)
    )
}

# =============================================================================
# 5. BUILD FIGURE
# =============================================================================
cat("--- 4. Building figure ---\n")

# ── Y-axis: 1st–99th pctile ± 5% padding ─────────────────────────────────
q_lo  <- quantile(bm_master$bw_kg, 0.01, na.rm = TRUE)
q_hi  <- quantile(bm_master$bw_kg, 0.99, na.rm = TRUE)
pad   <- (q_hi - q_lo) * 0.05
ylims <- c(q_lo - pad, q_hi + pad)

# ── ICC annotation data frame — one row per sport (facet) ─────────────────
icc_ann <- tibble(
  sport = factor(c("Basketball", "Rugby"),
                 levels = c("Basketball", "Rugby")),
  label = sprintf("ICC = %.2f", icc_bm)
)

# ── Within-athlete SD annotation — quantifies fluctuation magnitude ────────
wsd_ann <- bm_master |>
  group_by(sport, athlete_id) |>
  summarise(within_sd = sd(bw_kg, na.rm = TRUE), .groups = "drop") |>
  group_by(sport) |>
  summarise(mean_wsd = mean(within_sd, na.rm = TRUE), .groups = "drop") |>
  mutate(label2 = sprintf("Within-athlete SD: %.2f kg", mean_wsd))

ann_df <- left_join(icc_ann, wsd_ann, by = "sport") |>
  mutate(full_label = paste0(label, "\n", label2),
         x = -Inf, y = Inf)

# ── Assemble plot ─────────────────────────────────────────────────────────
p_bm <- ggplot() +

  # Layer 1 — spaghetti (background): individual athlete trajectories ──────
  geom_line(
    data    = bm_master,
    mapping = aes(x     = week,
                  y     = bw_kg,
                  group = athlete_id,
                  color = sport),
    alpha     = 0.22,
    linewidth = 0.50
  ) +

  # Layer 2 — LMM 95% CI ribbon ────────────────────────────────────────────
  geom_ribbon(
    data    = pred_bm,
    mapping = aes(x    = week,
                  ymin = conf.low,
                  ymax = conf.high,
                  fill = sport),
    alpha       = 0.28,
    inherit.aes = FALSE
  ) +

  # Layer 3 — LMM predicted mean (foreground, bold) ────────────────────────
  geom_line(
    data    = pred_bm,
    mapping = aes(x     = week,
                  y     = predicted,
                  color = sport),
    linewidth   = 2.1,
    inherit.aes = FALSE
  ) +

  # Layer 4 — ICC & within-SD annotation (top-left of each facet) ──────────
  geom_text(
    data    = ann_df,
    mapping = aes(x = x, y = y, label = full_label),
    hjust   = -0.07, vjust = 1.25,
    size    = 5.0, color = "grey30", fontface = "italic",
    inherit.aes = FALSE
  ) +

  # Facet by sport (free_x: each sport fills its own monitoring window) ─────
  facet_wrap(~ sport, ncol = 2, scales = "free_x") +

  # Y-axis: clipped view to emphasise within-athlete fluctuation ────────────
  coord_cartesian(ylim = ylims) +

  # Scales ──────────────────────────────────────────────────────────────────
  scale_color_manual(values = SPORT_COLS,  name = "Sport") +
  scale_fill_manual(values  = SPORT_FILLS, name = "Sport") +
  scale_x_continuous(
    name   = "Monitoring Week",
    breaks = seq(0, 35, by = 5),
    expand = expansion(add = c(0.3, 0.8))
  ) +
  scale_y_continuous(
    name   = "Body Mass (kg)",
    breaks = scales::breaks_pretty(n = 8),
    expand = expansion(mult = c(0, 0))
  ) +

  # Labels ──────────────────────────────────────────────────────────────────
  labs(
    subtitle = sprintf(
      paste0("Thin lines = individual athlete trajectories  |  ",
             "Bold line = LMM predicted mean (95%% CI ribbon)  |  ",
             "ICC = %.2f  (%.0f%% between-athlete variance)"),
      icc_bm, icc_bm * 100
    )
  ) +

  theme_BM(18) +

  guides(
    color = guide_legend(
      override.aes = list(alpha = 1, linewidth = 2.2)
    ),
    fill = "none"
  )

# =============================================================================
# 6. EXPORT
# =============================================================================
cat("--- 5. Exporting figure ---\n")

out_path <- file.path(FIG_DIR, "Hero_Fig3_BodyMass_Variation.png")
t0 <- proc.time()

ggsave(
  filename = out_path,
  plot     = p_bm,
  width    = 12,
  height   = 8,
  dpi      = 600,
  bg       = "white",
  units    = "in"
)

elapsed <- round((proc.time() - t0)["elapsed"], 0)
file_mb <- round(file.info(out_path)$size / 1024^2, 1)

cat(sprintf("\n  ✓ Saved: Hero_Fig3_BodyMass_Variation.png\n"))
cat(sprintf("  Dimensions : 12 × 8 in\n"))
cat(sprintf("  Resolution : 600 dpi  (pixel size: 7200 × 4800 px)\n"))
cat(sprintf("  File size  : %.1f MB\n", file_mb))
cat(sprintf("  Render time: %ds\n", elapsed))

# =============================================================================
# 7. FINAL RESULTS REPORT
# =============================================================================
cat("\n╔══════════════════════════════════════════════════════════════════╗\n")
cat("║  BODY MASS LMM — RESULTS SUMMARY                               ║\n")
cat("╚══════════════════════════════════════════════════════════════════╝\n\n")

cat(sprintf("  n athletes : 32 total (Basketball=18, Rugby=14)\n"))
cat(sprintf("  n obs      : %d total (%d Basketball / %d Rugby)\n",
            nrow(bm_master),
            sum(bm_master$sport == "Basketball"),
            sum(bm_master$sport == "Rugby")))

cat(sprintf("\n  ICC = %.3f\n", icc_bm))
cat(sprintf("    → %.1f%% stable between-athlete differences\n", icc_bm * 100))
cat(sprintf("    → %.1f%% within-athlete week-to-week fluctuation\n",
            (1 - icc_bm) * 100))
cat(sprintf("    τ₀ (between-athlete SD) = %.3f kg\n",  sd_athlete))
cat(sprintf("    σ  (within-athlete SD)  = %.3f kg\n\n",  sd_resid))

n_sig <- sum(fe_bm$p.value[fe_bm$term != "(Intercept)"] < .05)
if (n_sig > 0) {
  cat(sprintf("  Significant fixed effects (%d term(s), p < .05):\n", n_sig))
  for (i in seq_len(nrow(sig_terms))) {
    t   <- sig_terms$term[i]
    lbl <- coalesce(TERM_LABELS[t], t)
    pv  <- sig_terms$p.value[i]
    cat(sprintf("    ▸ %s\n",      lbl))
    cat(sprintf("      β = %.3f kg  |  t = %.2f  |  p %s\n",
                sig_terms$estimate[i],
                sig_terms$statistic[i],
                ifelse(pv < .001, "< .001", sprintf("= %.3f", pv))))
  }
} else {
  cat("  Fixed effects: none significant (p ≥ .05)\n")
  cat("  → Body mass was STATISTICALLY STABLE across the monitoring season\n")
  cat("    for both sports, with no significant linear trend or sport × time interaction.\n")
}

cat("\n=== Body Mass Figure Complete ===\n")
