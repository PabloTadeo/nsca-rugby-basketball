# =============================================================================
# NSCA Poster — Linear Mixed Models (LMM) Analysis
# Neuromuscular Performance Temporal Trends: Basketball vs. Rugby
#
# Models : DSI | Jump Height (cm) | CMJ Peak Force (N) | Peak Power (W)
# Formula: outcome ~ sport * week + (1 | athlete_id)
# RE     : Random intercept per athlete (accounts for repeated measures)
# p-vals : Satterthwaite approximation via lmerTest
# Optim  : BOBYQA (robust for small-N repeated-measures designs)
#
# Authors: Pablo Tadeo Ríos-Gallardo, PhD & Samuel Montalvo, PhD
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(lubridate)
  library(lme4)
  library(lmerTest)    # Satterthwaite df + p-values for lmer objects
  library(broom.mixed) # tidy() extractor for lmer
  library(gt)
  library(gtExtras)
  library(scales)
})

BASE_DIR <- "/Volumes/ADATA HD710 PRO/fod_rugby"
RAW_DIR  <- file.path(BASE_DIR, "Raw_data")
TBL_DIR  <- file.path(BASE_DIR, "Tables")

cat("=== NSCA Poster — LMM Analysis (Basketball & Rugby) ===\n\n")

# =============================================================================
# 1. DATA PIPELINE  (self-contained; includes peak_power_w)
# =============================================================================
cat("--- 1. Loading and preparing master dataset ---\n")

read_utf8 <- function(f)
  read_csv(file.path(RAW_DIR, f),
           locale        = locale(encoding = "UTF-8"),
           show_col_types = FALSE,
           name_repair   = "minimal")

num <- function(x) suppressWarnings(as.numeric(x))

# ---- Basketball CMJ (includes Peak Power) ----
bas_cmj <- read_utf8("CMJ_Basket - CMJ.csv") |>
  rename(
    athlete_name     = Name,
    athlete_id       = ExternalId,
    week             = Week,
    bw_kg            = `BW [KG]`,
    jump_height_cm   = `Jump Height (Imp-Mom) [cm]`,
    peak_power_w     = `Peak Power [W]`,
    rsi_mod          = `RSI-modified [m/s]`,
    cmj_peak_force_n = `Concentric Peak Force [N]`
  ) |>
  select(athlete_id, athlete_name, week, bw_kg,
         jump_height_cm, peak_power_w, rsi_mod, cmj_peak_force_n) |>
  mutate(
    sport            = "Basketball",
    week             = as.integer(num(week)),
    across(c(bw_kg, jump_height_cm, peak_power_w,
             rsi_mod, cmj_peak_force_n), num),
    jump_height_cm   = ifelse(jump_height_cm < 0 | jump_height_cm > 120,
                               NA, jump_height_cm),
    peak_power_w     = ifelse(peak_power_w   <= 0, NA, peak_power_w),
    rsi_mod          = ifelse(rsi_mod         < 0 | rsi_mod > 5,
                               NA, rsi_mod),
    cmj_peak_force_n = ifelse(cmj_peak_force_n <= 0, NA, cmj_peak_force_n)
  ) |>
  filter(!is.na(athlete_id), !is.na(week))

# ---- Basketball IMTP ----
bas_imtp <- read_utf8("IMPT_Basket - Hoja 1.csv") |>
  rename(
    athlete_id        = ExternalId,
    week              = Week,
    imtp_peak_force_n = `Peak Vertical Force [N]`
  ) |>
  select(athlete_id, week, imtp_peak_force_n) |>
  mutate(
    week              = as.integer(num(week)),
    imtp_peak_force_n = num(imtp_peak_force_n),
    imtp_peak_force_n = ifelse(imtp_peak_force_n <= 0, NA, imtp_peak_force_n)
  ) |>
  filter(!is.na(athlete_id), !is.na(week))

# ---- Rugby Demographics ----
rug_demo <- read_utf8("CMJ_IMTP_RUGBY - Hoja 1.csv") |>
  rename_with(~ trimws(.x)) |>
  rename(
    athlete_id = externalId,
    weight_kg  = `peso (kg)`,
    height_cm  = `altura(cm)`
  ) |>
  mutate(weight_kg = num(weight_kg), height_cm = num(height_cm)) |>
  select(athlete_id, weight_kg, height_cm)

# ---- Rugby CMJ (includes Peak Power) ----
rug_cmj <- read_utf8("CMJ_Rugby - Hoja 1.csv") |>
  rename(
    athlete_name     = Name,
    athlete_id       = ExternalId,
    bw_kg            = `BW [KG]`,
    jump_height_cm   = `Jump Height (Imp-Mom) [cm]`,
    peak_power_w     = `Peak Power [W]`,
    rsi_mod          = `RSI-modified [m/s]`,
    cmj_peak_force_n = `Concentric Peak Force [N]`
  ) |>
  select(athlete_id, athlete_name, week, bw_kg,
         jump_height_cm, peak_power_w, rsi_mod, cmj_peak_force_n) |>
  mutate(
    sport            = "Rugby",
    week             = as.integer(num(week)),
    across(c(bw_kg, jump_height_cm, peak_power_w,
             rsi_mod, cmj_peak_force_n), num),
    jump_height_cm   = ifelse(jump_height_cm < 0 | jump_height_cm > 120,
                               NA, jump_height_cm),
    peak_power_w     = ifelse(peak_power_w   <= 0, NA, peak_power_w),
    rsi_mod          = ifelse(rsi_mod         < 0 | rsi_mod > 5,
                               NA, rsi_mod),
    cmj_peak_force_n = ifelse(cmj_peak_force_n <= 0, NA, cmj_peak_force_n)
  ) |>
  filter(!is.na(athlete_id), !is.na(week))

# ---- Rugby IMTP ----
rug_imtp <- read_utf8("IMPT_Rugby - Hoja 1.csv") |>
  rename(
    athlete_id        = ExternalId,
    imtp_peak_force_n = `Peak Vertical Force [N]`
  ) |>
  select(athlete_id, week, imtp_peak_force_n) |>
  mutate(
    week              = as.integer(num(week)),
    imtp_peak_force_n = num(imtp_peak_force_n),
    imtp_peak_force_n = ifelse(imtp_peak_force_n <= 0, NA, imtp_peak_force_n)
  ) |>
  filter(!is.na(athlete_id), !is.na(week))

# ---- Aggregate per athlete × week ----
agg_cmj <- function(df)
  df |>
    group_by(athlete_id, athlete_name, sport, week) |>
    summarise(
      across(c(bw_kg, jump_height_cm, peak_power_w,
               rsi_mod, cmj_peak_force_n),
             ~ mean(.x, na.rm = TRUE)),
      .groups = "drop"
    ) |>
    mutate(across(where(is.numeric), ~ ifelse(is.nan(.x), NA_real_, .x)))

agg_imtp <- function(df)
  df |>
    group_by(athlete_id, week) |>
    summarise(
      imtp_peak_force_n = max(imtp_peak_force_n, na.rm = TRUE),
      .groups = "drop"
    ) |>
    mutate(imtp_peak_force_n = ifelse(is.infinite(imtp_peak_force_n),
                                      NA_real_, imtp_peak_force_n))

# ---- Merge within sports + DSI ----
bas_perf <- left_join(agg_cmj(bas_cmj), agg_imtp(bas_imtp),
                      by = c("athlete_id", "week")) |>
  mutate(
    dsi = cmj_peak_force_n / imtp_peak_force_n,
    dsi = ifelse(dsi < 0.3 | dsi > 2.0, NA_real_, dsi)
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

# ── Body mass dataset (session-level BW from CMJ files, NO demo override) ──
# Rugby athletes' static "peso (kg)" is intentionally excluded so that
# genuine week-to-week fluctuations in measured body mass are preserved.
bm_master <- bind_rows(
  bas_cmj |>
    group_by(athlete_id, athlete_name, sport, week) |>
    summarise(bw_kg = mean(bw_kg, na.rm = TRUE), .groups = "drop"),
  rug_cmj |>
    group_by(athlete_id, athlete_name, sport, week) |>
    summarise(bw_kg = mean(bw_kg, na.rm = TRUE), .groups = "drop")
) |>
  mutate(
    bw_kg = ifelse(is.nan(bw_kg) | bw_kg < 35 | bw_kg > 160, NA_real_, bw_kg),
    sport = factor(sport, levels = c("Basketball", "Rugby"))
  ) |>
  filter(!is.na(bw_kg))

# Audit
cat(sprintf(
  "Master: %d athlete-week rows | Basketball n_athletes=%d (%d obs) | Rugby n_athletes=%d (%d obs)\n\n",
  nrow(master),
  n_distinct(master$athlete_id[master$sport == "Basketball"]),
  sum(master$sport == "Basketball"),
  n_distinct(master$athlete_id[master$sport == "Rugby"]),
  sum(master$sport == "Rugby")
))

# Data completeness per LMM variable
lmm_vars <- c(
  "Dynamic Strength Index (DSI)" = "dsi",
  "Jump Height (cm)"             = "jump_height_cm",
  "CMJ Peak Force (N)"           = "cmj_peak_force_n",
  "Peak Power (W)"               = "peak_power_w",
  "Body Mass (kg)"               = "bw_kg"
)

cat("Non-missing observations per variable:\n")
for (lbl in names(lmm_vars)) {
  col <- lmm_vars[[lbl]]
  n_total <- sum(!is.na(master[[col]]))
  n_bas   <- sum(!is.na(master[[col]][master$sport == "Basketball"]))
  n_rug   <- sum(!is.na(master[[col]][master$sport == "Rugby"]))
  cat(sprintf("  %-32s  total=%d  Basketball=%d  Rugby=%d\n",
              lbl, n_total, n_bas, n_rug))
}
cat("\n")

# =============================================================================
# 2. LMM CONFIGURATION
# =============================================================================
cat("--- 2. LMM configuration ---\n")
cat("  Formula  : outcome ~ sport * week + (1 | athlete_id)\n")
cat("  RE       : Random intercept per athlete (within-person correlation)\n")
cat("  p-values : Satterthwaite approximation (lmerTest)\n")
cat("  Optimizer: BOBYQA (maxfun = 200,000) — robust for small N\n")
cat("  Scaling  : CMJ Peak Force and Peak Power fitted in scaled units (/1000)\n")
cat("             to prevent ill-conditioning; estimates back-transformed.\n\n")

LMM_CTRL <- lmerControl(
  optimizer = "bobyqa",
  optCtrl   = list(maxfun = 2e5)
)

# =============================================================================
# 3. FIT MODELS
# =============================================================================

#' Fit LMM, check convergence, extract fixed effects
#'
#' @param outcome_col  Column name in `master`
#' @param label        Display label for the table
#' @param data         Data frame (master)
#' @param scale_factor Divisor applied before fitting; estimates back-transformed
fit_lmm <- function(outcome_col, label, data, scale_factor = 1) {

  df_fit <- data |>
    filter(!is.na(.data[[outcome_col]]),
           !is.na(week), !is.na(sport)) |>
    mutate(y_fit = .data[[outcome_col]] / scale_factor)

  n_obs      <- nrow(df_fit)
  n_athletes <- n_distinct(df_fit$athlete_id)

  cat(sprintf("  %-35s  n_obs=%d  n_athletes=%d",
              label, n_obs, n_athletes))
  if (scale_factor != 1)
    cat(sprintf("  [scaled ÷%g]", scale_factor))
  cat("\n")

  fit <- lmer(
    y_fit ~ sport * week + (1 | athlete_id),
    data    = df_fit,
    REML    = TRUE,
    control = LMM_CTRL
  )

  # ---- Convergence check ----
  conv_msgs <- tryCatch(fit@optinfo$conv$lme4$messages,
                        error = function(e) NULL)
  if (is.null(conv_msgs) || length(conv_msgs) == 0) {
    cat("    ✓  Model converged with no warnings\n")
  } else {
    cat(sprintf("    ⚠  %s\n", paste(conv_msgs, collapse = "; ")))
    # Attempt fallback with Nelder_Mead if bobyqa warns
    cat("    ↻  Re-fitting with Nelder_Mead optimizer...\n")
    fit <- lmer(
      y_fit ~ sport * week + (1 | athlete_id),
      data    = df_fit,
      REML    = TRUE,
      control = lmerControl(optimizer    = "Nelder_Mead",
                            optCtrl      = list(maxfun = 2e5))
    )
    conv2 <- tryCatch(fit@optinfo$conv$lme4$messages, error = function(e) NULL)
    if (is.null(conv2) || length(conv2) == 0) {
      cat("    ✓  Re-fit converged with Nelder_Mead\n")
    } else {
      cat(sprintf("    ✗  Still flagged: %s\n", paste(conv2, collapse = "; ")))
    }
  }

  # ---- Variance components ----
  vc             <- as.data.frame(VarCorr(fit))
  sd_athlete_raw <- vc$sdcor[vc$grp == "athlete_id"]
  sd_resid_raw   <- vc$sdcor[vc$grp == "Residual"]
  sd_athlete     <- sd_athlete_raw * scale_factor
  sd_resid       <- sd_resid_raw   * scale_factor
  icc            <- sd_athlete^2 / (sd_athlete^2 + sd_resid^2)
  cat(sprintf("    RE SD = %.4f | Residual SD = %.4f | ICC = %.3f\n\n",
              sd_athlete, sd_resid, icc))

  # ---- Fixed effects (Satterthwaite df) ----
  res <- broom.mixed::tidy(
    fit,
    effects    = "fixed",
    conf.int   = TRUE,
    conf.level = 0.95
  ) |>
    mutate(
      Outcome   = label,
      # Back-transform to original units
      estimate  = estimate  * scale_factor,
      conf.low  = conf.low  * scale_factor,
      conf.high = conf.high * scale_factor,
      std.error = std.error * scale_factor
    )

  list(fit       = fit,
       results   = res,
       icc       = icc,
       n_obs     = n_obs,
       n_athletes = n_athletes)
}

cat("--- 3. Fitting models ---\n\n")

m_dsi  <- fit_lmm("dsi",              "Dynamic Strength Index (DSI)", master,
                   scale_factor = 1)
m_jh   <- fit_lmm("jump_height_cm",   "Jump Height (cm)",             master,
                   scale_factor = 1)
m_cmj  <- fit_lmm("cmj_peak_force_n",  "CMJ Peak Force (N)",  master,
                   scale_factor = 1000)   # fit in kN → estimates in N
m_imtp <- fit_lmm("imtp_peak_force_n", "IMTP Peak Force (N)", master,
                   scale_factor = 1000)   # fit in kN → estimates in N
m_pp   <- fit_lmm("peak_power_w",      "Peak Power (W)",      master,
                   scale_factor = 1000)   # fit in kW → estimates in W
# Body mass fitted on bm_master (session-level BW, no rugby demo override)
m_bm   <- fit_lmm("bw_kg",            "Body Mass (kg)",       bm_master,
                   scale_factor = 1)

# =============================================================================
# 4. COMPILE RESULTS
# =============================================================================
cat("--- 4. Compiling fixed effects results ---\n\n")

TERM_LABELS <- c(
  "(Intercept)"     = "Intercept (Basketball, Week 0)",
  "sportRugby"      = "Sport: Rugby vs. Basketball",
  "week"            = "Week (linear trend — Basketball)",
  "sportRugby:week" = "Sport × Week Interaction"
)

OUTCOME_LEVELS <- c(
  "Dynamic Strength Index (DSI)",
  "Jump Height (cm)",
  "CMJ Peak Force (N)",
  "IMTP Peak Force (N)",
  "Peak Power (W)",
  "Body Mass (kg)"
)

all_results <- bind_rows(
  m_dsi$results,
  m_jh$results,
  m_cmj$results,
  m_imtp$results,
  m_pp$results,
  m_bm$results
)

tbl_data <- all_results |>
  # Drop intercept — not interpretable for poster audience
  filter(term != "(Intercept)") |>
  mutate(
    Term    = TERM_LABELS[term],
    Term    = ifelse(is.na(Term), term, Term),
    Outcome = factor(Outcome, levels = OUTCOME_LEVELS),
    # Smart decimal precision per outcome
    decimals = case_when(
      Outcome %in% c("Dynamic Strength Index (DSI)")       ~ 3L,
      Outcome %in% c("Jump Height (cm)", "Body Mass (kg)") ~ 2L,
      TRUE                                                  ~ 1L   # N and W
    ),
    # Format Estimate (95% CI)
    est_ci = pmap_chr(
      list(estimate, conf.low, conf.high, decimals),
      function(e, lo, hi, d) {
        fmt <- paste0("%.", d, "f")
        sprintf(paste0(fmt, " [", fmt, ", ", fmt, "]"), e, lo, hi)
      }
    ),
    se_fmt  = pmap_chr(list(std.error, decimals),
                       function(s, d) sprintf(paste0("%.", d, "f"), s)),
    t_fmt   = sprintf("%.2f",   statistic),
    df_fmt  = sprintf("%.1f",   df),
    p_fmt   = case_when(
      p.value < .001 ~ "<.001",
      TRUE           ~ sprintf("%.3f", p.value)
    ),
    p_sig   = p.value < .05
  ) |>
  arrange(Outcome, match(term, c("sportRugby", "week", "sportRugby:week"))) |>
  select(Outcome, Term, est_ci, se_fmt, t_fmt, df_fmt, p_fmt, p_sig, p.value)

# Print to terminal
cat("Fixed effects summary (excluding intercept):\n")
cat(strrep("─", 100), "\n")
tbl_data |>
  select(Outcome, Term, est_ci, p_fmt) |>
  mutate(sig = ifelse(p_fmt != "<.001" & as.numeric(p_fmt) >= .05, " ", "*")) |>
  print(n = Inf)
cat("\n")

# =============================================================================
# 5. ICC SUMMARY
# =============================================================================
icc_tbl <- tibble(
  Outcome  = OUTCOME_LEVELS,
  ICC      = c(m_dsi$icc, m_jh$icc, m_cmj$icc, m_imtp$icc, m_pp$icc, m_bm$icc),
  n_obs    = c(m_dsi$n_obs, m_jh$n_obs, m_cmj$n_obs, m_imtp$n_obs,
               m_pp$n_obs, m_bm$n_obs),
  n_ath    = c(m_dsi$n_athletes, m_jh$n_athletes, m_cmj$n_athletes,
               m_imtp$n_athletes, m_pp$n_athletes, m_bm$n_athletes)
)

cat("Intraclass Correlation Coefficients (ICC — between-athlete variance):\n")
icc_tbl |> print()
cat("\n")

icc_footnote <- sprintf(
  "ICC (athlete random intercept): DSI = %.2f, Jump Height = %.2f, CMJ PF = %.2f, IMTP PF = %.2f, Peak Power = %.2f, Body Mass = %.2f. ICC reflects the proportion of total variance attributable to stable between-athlete differences.",
  m_dsi$icc, m_jh$icc, m_cmj$icc, m_imtp$icc, m_pp$icc, m_bm$icc
)

# =============================================================================
# 5b. BODY MASS COVARIATE TEST (Likelihood Ratio Tests)
# =============================================================================
# For each performance outcome, we test whether adding weekly body mass as a
# time-varying covariate significantly improves model fit over the base model.
# LRT requires ML (not REML) and identical samples for both models.
# Body mass covariate is session-level BW from bm_master (no static demo override).
# =============================================================================
cat("--- 5b. Body Mass Covariate Testing (LRT) ---\n")
cat("  Formula (base)    : outcome ~ sport * week + (1 | athlete_id)\n")
cat("  Formula (adjusted): outcome ~ sport * week + bm_cov + (1 | athlete_id)\n")
cat("  Note: fitted on shared sample (bm non-NA); REML=FALSE for fair LRT.\n\n")

bm_cov_df <- bm_master |> select(athlete_id, week, bm_cov = bw_kg)

perf_specs <- list(
  list(label = "Dynamic Strength Index (DSI)",  col = "dsi",               sf = 1,    data = master),
  list(label = "Jump Height (cm)",              col = "jump_height_cm",    sf = 1,    data = master),
  list(label = "CMJ Peak Force (N)",            col = "cmj_peak_force_n",  sf = 1000, data = master),
  list(label = "IMTP Peak Force (N)",           col = "imtp_peak_force_n", sf = 1000, data = master),
  list(label = "Peak Power (W)",                col = "peak_power_w",      sf = 1000, data = master)
)

lrt_results <- purrr::map_dfr(perf_specs, function(ps) {
  df_joint <- ps$data |>
    filter(!is.na(.data[[ps$col]]), !is.na(week), !is.na(sport)) |>
    mutate(y_fit = .data[[ps$col]] / ps$sf) |>
    left_join(bm_cov_df, by = c("athlete_id", "week")) |>
    filter(!is.na(bm_cov))

  fit_base_ml <- lmer(y_fit ~ sport * week + (1 | athlete_id),
                      data = df_joint, REML = FALSE, control = LMM_CTRL)
  fit_adj_ml  <- lmer(y_fit ~ sport * week + bm_cov + (1 | athlete_id),
                      data = df_joint, REML = FALSE, control = LMM_CTRL)

  lrt_out  <- anova(fit_base_ml, fit_adj_ml, refit = FALSE)
  chi2     <- lrt_out[2, "Chisq"]
  pval     <- lrt_out[2, "Pr(>Chisq)"]
  delta_aic <- AIC(fit_adj_ml) - AIC(fit_base_ml)

  # BM coefficient (back-transformed)
  bm_b   <- fixef(fit_adj_ml)["bm_cov"] * ps$sf
  bm_se  <- sqrt(vcov(fit_adj_ml)["bm_cov", "bm_cov"]) * ps$sf
  bm_t   <- fixef(fit_adj_ml)["bm_cov"] /
              sqrt(vcov(fit_adj_ml)["bm_cov", "bm_cov"])

  cat(sprintf("  %-35s n=%d  χ²=%5.2f  p=%s  ΔAIC=%+5.1f  β_BM=%+.4f\n",
              ps$label, nrow(df_joint), chi2,
              ifelse(pval < .001, "< .001", sprintf("%.3f", pval)),
              delta_aic, bm_b))

  tibble(
    Outcome   = ps$label,
    n_joint   = nrow(df_joint),
    chi2      = chi2,
    p_lrt     = pval,
    delta_aic = delta_aic,
    bm_b      = bm_b,
    bm_se     = bm_se,
    bm_t      = bm_t,
    bm_sig    = pval < .05
  )
})

any_sig  <- any(lrt_results$bm_sig)
n_sig_bm <- sum(lrt_results$bm_sig)

cat(sprintf("\n  Significant BM covariate effects: %d / %d outcomes (p < .05)\n",
            n_sig_bm, nrow(lrt_results)))

bm_footnote_text <- if (any_sig) {
  # Build per-outcome BM coefficient summary string
  bm_detail <- lrt_results |>
    filter(bm_sig) |>
    mutate(row_txt = sprintf("%s: \u03b2\u209a\u209c = %+.2f (p%s)",
                             Outcome,
                             bm_b,
                             ifelse(p_lrt < .001, " < .001",
                                    paste0(" = ", sprintf("%.3f", p_lrt))))) |>
    pull(row_txt) |>
    paste(collapse = "; ")
  sprintf(paste0(
    "Body mass (kg/session, time-varying) was tested as a covariate via Likelihood Ratio Test ",
    "(formula: outcome \u223c Sport \u00d7 Week + Body Mass + (1 | Athlete); REML = FALSE). ",
    "Significant improvement in model fit (p < .05) for %d/%d outcomes. ",
    "Between-athlete BM variance is largely captured by the random intercept (ICC\u209a\u209c = 0.97), ",
    "so only the within-athlete BM signal (session-level fluctuation, \u03c3 = 1.47 kg) contributes. ",
    "Time-varying BM coefficients from adjusted models: %s. ",
    "Base-model estimates (Sport, Week, Sport\u00d7Week) are reported in all table rows above; ",
    "these are nearly identical to adjusted-model estimates due to the high BM ICC."),
    n_sig_bm, nrow(lrt_results), bm_detail)
} else {
  paste0(
    "Body mass (kg/session) was tested as a time-varying covariate via Likelihood Ratio Test ",
    "(formula: outcome \u223c Sport \u00d7 Week + Body Mass + (1 | Athlete); REML = FALSE). ",
    "Addition of body mass did not significantly improve model fit for any outcome ",
    "(all LRT \u03c7\u00b2 p > .05). ",
    "This is consistent with ICC\u209a\u209c = 0.97: stable between-athlete differences in body mass ",
    "are already captured by the random intercept. Base models are reported."
  )
}

cat("\n╔═══════════════════════════════════════════════════════════════════════╗\n")
cat("║  BODY MASS COVARIATE — METHODOLOGICAL EVALUATION                     ║\n")
cat("╚═══════════════════════════════════════════════════════════════════════╝\n\n")
cat("Q1: Was it statistically necessary to adjust for body mass?\n")
if (!any_sig) {
  cat("  ✗  NO — Body mass did not significantly improve fit for any outcome.\n")
  cat("  Reason: ICC_BM = 0.972 → 97.2% of BM variance is between-athlete.\n")
  cat("  The random intercept (1|athlete_id) already captures that structure.\n")
  cat("  The time-varying BM signal (σ = 1.47 kg) adds negligible predictive power.\n")
} else {
  cat(sprintf("  ✓  YES for %d outcome(s):\n", n_sig_bm))
  for (i in seq_len(nrow(lrt_results))) {
    if (lrt_results$bm_sig[i]) {
      cat(sprintf("    ▸ %s: χ²=%.2f, ΔAIC=%.1f, β_BM=%+.3f\n",
                  lrt_results$Outcome[i],
                  lrt_results$chi2[i],
                  lrt_results$delta_aic[i],
                  lrt_results$bm_b[i]))
    }
  }
}

cat("\nQ2: Is DSI sensitive to body mass fluctuations or time?\n")
dsi_lrt <- lrt_results |> filter(Outcome == "Dynamic Strength Index (DSI)")
dsi_fe  <- tbl_data |>   filter(Outcome == "Dynamic Strength Index (DSI)")
cat(sprintf("  BM covariate: χ²=%.2f, p=%s, ΔAIC=%+.1f → %s\n",
            dsi_lrt$chi2,
            ifelse(dsi_lrt$p_lrt < .001, "<.001", sprintf("%.3f", dsi_lrt$p_lrt)),
            dsi_lrt$delta_aic,
            ifelse(dsi_lrt$bm_sig, "SENSITIVE to BM ✓", "NOT sensitive to BM")))
dsi_week_p <- tbl_data |>
  filter(Outcome == "Dynamic Strength Index (DSI)",
         grepl("Week", Term)) |>
  pull(p.value)
cat(sprintf("  Week trend:   p = %s → %s\n\n",
            ifelse(dsi_week_p[1] < .001, "<.001", sprintf("%.3f", dsi_week_p[1])),
            ifelse(dsi_week_p[1] < .05, "SIGNIFICANT temporal trend", "no significant trend")))

# =============================================================================
# 6. BUILD GT TABLE (NSCA STYLE)
# =============================================================================
cat("--- 5. Building Table 4 (gt HTML — NSCA style) ---\n")

# Header row colour (NSCA dark slate)
HDR_COL  <- "#1A2B3C"   # dark navy
ALT_COL  <- "#F4F7FA"   # very light blue-grey for alternating rows
SIG_COL  <- "#B22222"   # firebrick red for significant p-values

tbl4_gt <- tbl_data |>
  select(Outcome, Term, est_ci, se_fmt, t_fmt, df_fmt, p_fmt, p_sig, p.value) |>
  gt(groupname_col = "Outcome") |>

  # ── Header ──────────────────────────────────────────────────────────────────
  tab_header(
    title    = md("**Table 4. Linear Mixed Model Results: Neuromuscular Performance**"),
    subtitle = md(paste0(
      "*outcome* ~ Sport × Week + (1 | Athlete)  |  ",
      "REML estimation  |  Satterthwaite *df*  |  ",
      "Reference level: Basketball"
    ))
  ) |>

  # ── Column labels ───────────────────────────────────────────────────────────
  cols_label(
    Term    = "Fixed Effect",
    est_ci  = md("Estimate (95% CI)"),
    se_fmt  = md("*SE*"),
    t_fmt   = md("*t*"),
    df_fmt  = md("*df*"),
    p_fmt   = md("*p*")
  ) |>

  # ── Column alignment ────────────────────────────────────────────────────────
  cols_align(align = "left",
             columns = c(Term, est_ci)) |>
  cols_align(align = "center",
             columns = c(se_fmt, t_fmt, df_fmt, p_fmt)) |>

  # ── Row group styling — dark navy header ────────────────────────────────────
  tab_style(
    style     = list(
      cell_fill(color = HDR_COL),
      cell_text(weight = "bold", color = "white", size = px(12.5))
    ),
    locations = cells_row_groups()
  ) |>

  # ── Alternating row shading ─────────────────────────────────────────────────
  tab_style(
    style     = cell_fill(color = ALT_COL),
    locations = cells_body(rows = seq(2, nrow(tbl_data), 2))
  ) |>

  # ── Column label bold ───────────────────────────────────────────────────────
  tab_style(
    style     = cell_text(weight = "bold"),
    locations = cells_column_labels()
  ) |>

  # ── Significant p-values: bold + red ────────────────────────────────────────
  tab_style(
    style     = cell_text(weight = "bold", color = SIG_COL),
    locations = cells_body(
      columns = p_fmt,
      rows    = p_sig == TRUE
    )
  ) |>

  # ── Significant rows: subtle highlight ──────────────────────────────────────
  tab_style(
    style     = cell_text(style = "italic"),
    locations = cells_body(
      columns = Term,
      rows    = p_sig == TRUE
    )
  ) |>

  # ── Hide helper columns ──────────────────────────────────────────────────────
  cols_hide(columns = c(p_sig, p.value)) |>

  # ── Footnotes ───────────────────────────────────────────────────────────────
  tab_footnote(
    footnote  = md(paste0(
      "**Bold red** *p*-values indicate statistical significance (*p* < .05). ",
      "*SE* = Standard Error; *df* = degrees of freedom (Satterthwaite); CI = 95% Confidence Interval."
    )),
    locations = cells_column_labels(columns = p_fmt)
  ) |>
  tab_footnote(
    footnote  = paste0(
      "Sport coefficient = mean difference at Week 0 (Rugby minus Basketball). ",
      "Week coefficient = linear slope per monitoring week (Basketball). ",
      "Sport × Week = differential slope (Rugby slope minus Basketball slope)."
    ),
    locations = cells_column_labels(columns = Term)
  ) |>
  tab_footnote(
    footnote  = paste0(
      "CMJ Peak Force, IMTP Peak Force, and Peak Power were fitted in scaled units (÷1,000) ",
      "to prevent numerical ill-conditioning; all estimates are back-transformed and reported ",
      "in original units (N and W, respectively)."
    ),
    locations = cells_row_groups(
      groups = c("CMJ Peak Force (N)", "IMTP Peak Force (N)", "Peak Power (W)")
    )
  ) |>
  tab_footnote(
    footnote  = paste0(
      "Body Mass was modelled using the session-level BW recorded by the Hawkins force ",
      "platform at each CMJ testing session (BW [KG] column). For rugby athletes, the ",
      "static demographic weight was intentionally excluded so that genuine ",
      "session-to-session fluctuations are preserved in the random-effects structure."
    ),
    locations = cells_row_groups(groups = "Body Mass (kg)")
  ) |>
  tab_footnote(
    footnote  = md(bm_footnote_text),
    locations = cells_title(groups = "subtitle")
  ) |>
  tab_footnote(
    footnote  = icc_footnote,
    locations = cells_title(groups = "title")
  ) |>

  # ── Source note ─────────────────────────────────────────────────────────────
  tab_source_note(md(paste0(
    "CMJ = Countermovement Jump; IMTP = Isometric Mid-Thigh Pull; ",
    "DSI = Dynamic Strength Index (CMJ Peak Force ÷ IMTP Peak Force).  \n",
    "Models fitted with **lme4** v", packageVersion("lme4"), " + ",
    "**lmerTest** v", packageVersion("lmerTest"), " in R ",
    paste(R.version$major, R.version$minor, sep = "."),
    ". Optimizer: BOBYQA (maxfun = 200,000)."
  ))) |>

  # ── Table options (NSCA style) ───────────────────────────────────────────────
  opt_table_font(font = list(google_font("Roboto"), default_fonts())) |>
  tab_options(
    table.font.size                    = px(12),
    heading.title.font.size            = px(14),
    heading.subtitle.font.size         = px(11),
    column_labels.font.weight          = "bold",
    row_group.font.weight              = "bold",
    table.border.top.color             = "black",
    table.border.top.width             = px(2),
    column_labels.border.bottom.color  = "black",
    column_labels.border.bottom.width  = px(1.5),
    table.border.bottom.color          = "black",
    table.border.bottom.width          = px(2),
    row_group.border.top.color         = "grey20",
    row_group.border.top.width         = px(1),
    row_group.border.bottom.color      = "grey20",
    row_group.border.bottom.width      = px(0.5),
    footnotes.font.size                = px(10),
    source_notes.font.size             = px(10),
    data_row.padding                   = px(5),
    column_labels.padding              = px(7)
  )

tbl4_path <- file.path(TBL_DIR, "Table4_Mixed_Models_Results.html")
gtsave(tbl4_gt, tbl4_path)
cat(sprintf("  ✓ Saved: %s\n\n", basename(tbl4_path)))

# =============================================================================
# 7. FINAL SUMMARY REPORT
# =============================================================================
cat("╔══════════════════════════════════════════════════════════════════╗\n")
cat("║              LMM ANALYSIS — COMPLETE                            ║\n")
cat("╚══════════════════════════════════════════════════════════════════╝\n\n")

cat("CONVERGENCE & MODEL FIT:\n")
models_list <- list(m_dsi, m_jh, m_cmj, m_imtp, m_pp, m_bm)
labels_list <- OUTCOME_LEVELS
for (i in seq_along(models_list)) {
  m   <- models_list[[i]]
  lbl <- labels_list[i]
  cat(sprintf("  %-35s  n=%d athletes  %d obs  ICC=%.3f\n",
              lbl, m$n_athletes, m$n_obs, m$icc))
}

cat("\nSIGNIFICANT FIXED EFFECTS (p < .05):\n")
tbl_data |>
  filter(p_sig) |>
  mutate(p_display = ifelse(p_fmt == "<.001", "p<.001",
                            paste0("p=", p_fmt))) |>
  select(Outcome, Term, est_ci, p_display) |>
  as.data.frame() |>
  print(row.names = FALSE)

cat(sprintf("\nTable 4 saved → %s\n", basename(tbl4_path)))

# =============================================================================
# 5c. ADJUSTED MODELS (REML=TRUE) — BM COEFFICIENT EXTRACTION
# =============================================================================
# Refit adjusted models with REML=TRUE for proper SE/CI; extract β_BM per outcome.
# =============================================================================
cat("\n--- 5c. Adjusted REML models — BM coefficient rows ---\n")

adj_bm_rows <- purrr::map_dfr(perf_specs, function(ps) {
  df_joint <- ps$data |>
    filter(!is.na(.data[[ps$col]]), !is.na(week), !is.na(sport)) |>
    mutate(y_fit = .data[[ps$col]] / ps$sf) |>
    left_join(bm_cov_df, by = c("athlete_id", "week")) |>
    filter(!is.na(bm_cov))

  fit_adj_reml <- lmer(y_fit ~ sport * week + bm_cov + (1 | athlete_id),
                       data = df_joint, REML = TRUE, control = LMM_CTRL)

  broom.mixed::tidy(fit_adj_reml, effects = "fixed", conf.int = TRUE) |>
    filter(term == "bm_cov") |>
    mutate(
      across(c(estimate, conf.low, conf.high, std.error), ~ .x * ps$sf),
      Outcome = ps$label
    )
})

cat("  BM coefficients (adjusted REML models):\n")
for (i in seq_len(nrow(adj_bm_rows))) {
  r <- adj_bm_rows[i, ]
  cat(sprintf("  %-35s  beta=%+.4f [%.4f, %.4f]  p=%s\n",
              r$Outcome, r$estimate, r$conf.low, r$conf.high,
              ifelse(r$p.value < .001, "<.001", sprintf("%.3f", r$p.value))))
}
cat("\n")

# =============================================================================
# 5d. RANDOM EFFECTS SUMMARY — tau00, sigma2, ICC FROM BASE MODELS
# =============================================================================
cat("--- 5d. Random effects (VarCorr) from base models ---\n")

base_models_list_re <- list(
  list(fit = m_dsi$fit,  label = "Dynamic Strength Index (DSI)",  sf = 1),
  list(fit = m_jh$fit,   label = "Jump Height (cm)",              sf = 1),
  list(fit = m_cmj$fit,  label = "CMJ Peak Force (N)",            sf = 1000),
  list(fit = m_imtp$fit, label = "IMTP Peak Force (N)",           sf = 1000),
  list(fit = m_pp$fit,   label = "Peak Power (W)",                sf = 1000),
  list(fit = m_bm$fit,   label = "Body Mass (kg)",                sf = 1)
)

re_summary <- purrr::map_dfr(base_models_list_re, function(m) {
  vc    <- as.data.frame(VarCorr(m$fit))
  tau00 <- vc$vcov[vc$grp == "athlete_id"] * m$sf^2
  sig2  <- vc$vcov[vc$grp == "Residual"]   * m$sf^2
  tibble(
    Outcome = m$label,
    tau00   = tau00,
    tau0    = sqrt(tau00),
    sigma2  = sig2,
    sigma   = sqrt(sig2),
    icc     = tau00 / (tau00 + sig2)
  )
})

cat("  Random effects summary:\n")
print(re_summary, digits = 4)
cat("\n")

# =============================================================================
# 5e. CATEGORICAL WEEK MODELS — FATIGUE VALLEY DETECTION
# =============================================================================
# Fit: outcome ~ Sport * as.factor(week) + (1 | athlete_id)
# Purpose: identify which specific weeks show the largest Basketball deviations
# from Week 1 for DSI, CMJ PF, and IMTP PF.
# =============================================================================
cat("--- 5e. Categorical week models (fatigue valley detection) ---\n")
cat("  Formula: outcome ~ Sport * as.factor(week) + (1 | athlete_id)\n\n")

cat_outcomes <- list(
  list(col = "dsi",               label = "Dynamic Strength Index (DSI)",  sf = 1),
  list(col = "cmj_peak_force_n",  label = "CMJ Peak Force (N)",            sf = 1000),
  list(col = "imtp_peak_force_n", label = "IMTP Peak Force (N)",           sf = 1000)
)

cat_week_results <- purrr::map(cat_outcomes, function(ps) {
  df <- master |>
    filter(!is.na(.data[[ps$col]]), !is.na(week), !is.na(sport)) |>
    mutate(
      y_fit  = .data[[ps$col]] / ps$sf,
      week_f = factor(week, levels = sort(unique(week)))
    )

  cat(sprintf("  Fitting: %s  (n=%d obs, weeks %d-%d)\n",
              ps$label, nrow(df),
              min(df$week, na.rm = TRUE), max(df$week, na.rm = TRUE)))

  fit <- tryCatch(
    lmer(y_fit ~ sport * week_f + (1 | athlete_id),
         data = df, REML = TRUE, control = LMM_CTRL),
    error = function(e) {
      cat("    BOBYQA failed, retrying with Nelder_Mead...\n")
      lmer(y_fit ~ sport * week_f + (1 | athlete_id),
           data = df, REML = TRUE,
           control = lmerControl(optimizer = "Nelder_Mead",
                                 optCtrl   = list(maxfun = 2e5)))
    }
  )

  conv_msgs <- tryCatch(fit@optinfo$conv$lme4$messages, error = function(e) NULL)
  if (is.null(conv_msgs) || length(conv_msgs) == 0) {
    cat("    Converged\n")
  } else {
    cat(sprintf("    Warning: %s\n", paste(conv_msgs, collapse = "; ")))
  }

  # Main week effects = Basketball deviations from Week 1 (reference)
  fe <- broom.mixed::tidy(fit, effects = "fixed", conf.int = TRUE) |>
    filter(grepl("^week_f", term), !grepl("^sportRugby:week_f", term)) |>
    mutate(
      week_num   = as.integer(gsub("week_f", "", term)),
      est_orig   = estimate  * ps$sf,
      ci_lo_orig = conf.low  * ps$sf,
      ci_hi_orig = conf.high * ps$sf,
      sig        = p.value < .05
    ) |>
    arrange(est_orig)

  valleys <- head(fe, 3)   # top-3 most negative Basketball deviations

  unit_lbl <- if (ps$sf == 1000) "N" else "DSI units"
  cat(sprintf("    Biggest Basketball drops vs. Week 1 (%s):\n", unit_lbl))
  for (j in seq_len(nrow(valleys))) {
    r <- valleys[j, ]
    p_disp <- ifelse(r$p.value < .001, "<.001", sprintf("%.3f", r$p.value))
    cat(sprintf("      Week %2d: delta = %+.4f  p = %s  sig = %s\n",
                r$week_num, r$est_orig, p_disp,
                ifelse(r$sig, "YES", "no")))
  }
  cat("\n")

  list(fit = fit, fe = fe, valleys = valleys, label = ps$label, sf = ps$sf)
})
names(cat_week_results) <- sapply(cat_outcomes, `[[`, "label")

# =============================================================================
# 6b. BUILD COMPREHENSIVE TABLE (Table4_Mixed_Models_Complete.html)
# =============================================================================
cat("--- 6b. Building Table4_Mixed_Models_Complete.html ---\n")

# ── 1. Fixed effects rows (base model, excluding intercept) ─────────────────
fe_rows_c <- tbl_data |>
  mutate(row_type = "fe") |>
  select(Outcome, Term, est_ci, se_fmt, t_fmt, df_fmt, p_fmt, p_sig, p.value, row_type)

# ── 2. BM coefficient rows — adjusted REML, significant outcomes only ────────
sig_bm_outcomes <- lrt_results |> filter(bm_sig) |> pull(Outcome)

bm_fmt_rows <- adj_bm_rows |>
  filter(Outcome %in% sig_bm_outcomes) |>
  mutate(
    decimals = case_when(
      Outcome == "Dynamic Strength Index (DSI)"             ~ 3L,
      Outcome %in% c("Jump Height (cm)", "Body Mass (kg)") ~ 2L,
      TRUE                                                   ~ 1L
    ),
    Term    = "Body Mass covariate (\u03b2_BM)\u2020",
    est_ci  = pmap_chr(
      list(estimate, conf.low, conf.high, decimals),
      function(e, lo, hi, d) {
        fmt <- paste0("%.", d, "f")
        sprintf(paste0(fmt, " [", fmt, ", ", fmt, "]"), e, lo, hi)
      }
    ),
    se_fmt  = pmap_chr(list(std.error, decimals),
                       function(s, d) sprintf(paste0("%.", d, "f"), s)),
    t_fmt   = sprintf("%.2f",  statistic),
    df_fmt  = sprintf("%.1f",  df),
    p_fmt   = case_when(p.value < .001 ~ "<.001", TRUE ~ sprintf("%.3f", p.value)),
    p_sig   = p.value < .05,
    row_type = "bm"
  ) |>
  select(Outcome, Term, est_ci, se_fmt, t_fmt, df_fmt, p_fmt, p_sig, p.value, row_type)

# ── 3. Random effects rows (tau0, sigma, ICC) ────────────────────────────────
format_re_val <- function(val, outcome_str, re_type_str) {
  if (re_type_str == "icc") return(sprintf("%.3f", val))
  d <- dplyr::case_when(
    grepl("DSI", outcome_str)                                           ~ 4L,
    outcome_str %in% c("Jump Height (cm)", "Body Mass (kg)")           ~ 2L,
    TRUE                                                                 ~ 1L
  )
  sprintf(paste0("%.", d, "f"), val)
}

re_fmt_rows <- re_summary |>
  mutate(Outcome = factor(as.character(Outcome), levels = OUTCOME_LEVELS)) |>
  arrange(Outcome) |>
  tidyr::pivot_longer(
    cols      = c(tau0, sigma, icc),
    names_to  = "re_type",
    values_to = "re_val"
  ) |>
  mutate(
    Term     = case_when(
      re_type == "tau0"  ~ "Between-athlete SD (\u03c4\u2080)",
      re_type == "sigma" ~ "Within-athlete residual SD (\u03c3)",
      re_type == "icc"   ~ "ICC (athlete random intercept)"
    ),
    est_ci   = mapply(format_re_val, re_val, as.character(Outcome), re_type),
    se_fmt   = "\u2014",
    t_fmt    = "\u2014",
    df_fmt   = "\u2014",
    p_fmt    = "\u2014",
    p_sig    = FALSE,
    p.value  = NA_real_,
    row_type = "re"
  ) |>
  select(Outcome, Term, est_ci, se_fmt, t_fmt, df_fmt, p_fmt, p_sig, p.value, row_type)

# ── 4. Combine in correct display order ──────────────────────────────────────
row_order <- c(
  "Sport: Rugby vs. Basketball",
  "Week (linear trend \u2014 Basketball)",
  "Sport \u00d7 Week Interaction",
  "Body Mass covariate (\u03b2_BM)\u2020",
  "Between-athlete SD (\u03c4\u2080)",
  "Within-athlete residual SD (\u03c3)",
  "ICC (athlete random intercept)"
)

tbl_complete <- bind_rows(fe_rows_c, bm_fmt_rows, re_fmt_rows) |>
  mutate(
    Outcome  = factor(as.character(Outcome), levels = OUTCOME_LEVELS),
    Term_ord = match(Term, row_order)
  ) |>
  arrange(Outcome, Term_ord) |>
  select(-Term_ord)

cat(sprintf("  Combined table: %d rows  (fe=%d  bm=%d  re=%d)\n",
            nrow(tbl_complete),
            sum(tbl_complete$row_type == "fe"),
            sum(tbl_complete$row_type == "bm"),
            sum(tbl_complete$row_type == "re")))

# ── Categorical week footnote strings ────────────────────────────────────────
make_cat_fn <- function(cwr, unit) {
  v    <- cwr$fe |> arrange(est_orig)   # all week deviations, most negative first
  top3 <- head(v, 3)

  # Choose wording based on whether actual drops occurred
  has_drops <- min(top3$est_orig, na.rm = TRUE) < 0
  any_sig   <- any(top3$sig, na.rm = TRUE)

  p_str  <- ifelse(top3$p.value < .001, "<.001",
                   paste0("=", sprintf("%.3f", top3$p.value)))
  wk_str <- paste(
    sprintf("Week %d (\u0394=%+.4g %s, p%s%s)",
            top3$week_num, top3$est_orig, unit, p_str,
            ifelse(top3$sig, "*", "")),
    collapse = "; "
  )

  phrase <- if (has_drops) {
    "largest Basketball decreases from Week 1 baseline"
  } else {
    "no Basketball drops below Week 1; weeks with smallest gains"
  }

  sprintf(
    "Categorical week model (%s ~ Sport \u00d7 Week_factor + (1|Athlete)): %s: %s.",
    cwr$label, phrase, wk_str
  )
}

cat_fn_dsi  <- make_cat_fn(cat_week_results[["Dynamic Strength Index (DSI)"]], "DSI")
cat_fn_cmj  <- make_cat_fn(cat_week_results[["CMJ Peak Force (N)"]],           "N")
cat_fn_imtp <- make_cat_fn(cat_week_results[["IMTP Peak Force (N)"]],          "N")

# ── Colour constants (re-use or define if not in scope) ──────────────────────
AMBER_COL <- "#FFF3CD"
RE_COL    <- "#EBEBEB"

# ── Build gt table ────────────────────────────────────────────────────────────
tbl_complete_gt <- tbl_complete |>
  select(Outcome, Term, est_ci, se_fmt, t_fmt, df_fmt, p_fmt,
         p_sig, p.value, row_type) |>
  gt(groupname_col = "Outcome") |>

  tab_header(
    title    = md("**Table 4. Comprehensive Linear Mixed Model Results: Neuromuscular Performance**"),
    subtitle = md(paste0(
      "*outcome* ~ Sport \u00d7 Week + (1 | Athlete) \u00b7 REML \u00b7 ",
      "Satterthwaite *df* \u00b7 Reference: Basketball  \n",
      "\u2020 Body mass (BM) covariate from adjusted REML model \u2014 ",
      "shown only for outcomes with significant LRT improvement (*p* < .05). ",
      "Grey rows = random effects. Amber rows = BM covariate."
    ))
  ) |>

  cols_label(
    Term   = "Effect / Parameter",
    est_ci = md("Estimate (95% CI)"),
    se_fmt = md("*SE*"),
    t_fmt  = md("*t*"),
    df_fmt = md("*df*"),
    p_fmt  = md("*p*")
  ) |>

  cols_align(align = "left",   columns = c(Term, est_ci)) |>
  cols_align(align = "center", columns = c(se_fmt, t_fmt, df_fmt, p_fmt)) |>

  # Row group headers — dark navy
  tab_style(
    style     = list(
      cell_fill(color = HDR_COL),
      cell_text(weight = "bold", color = "white", size = px(12.5))
    ),
    locations = cells_row_groups()
  ) |>

  # BM rows — amber
  tab_style(
    style     = list(
      cell_fill(color = AMBER_COL),
      cell_text(style = "italic")
    ),
    locations = cells_body(rows = row_type == "bm")
  ) |>

  # RE rows — light grey background
  tab_style(
    style     = cell_fill(color = RE_COL),
    locations = cells_body(rows = row_type == "re")
  ) |>

  # RE rows — italic grey Term label
  tab_style(
    style     = cell_text(style = "italic", color = "grey40"),
    locations = cells_body(columns = Term, rows = row_type == "re")
  ) |>

  # Column labels bold
  tab_style(
    style     = cell_text(weight = "bold"),
    locations = cells_column_labels()
  ) |>

  # Significant p-values: bold red
  tab_style(
    style     = cell_text(weight = "bold", color = SIG_COL),
    locations = cells_body(columns = p_fmt, rows = p_sig == TRUE)
  ) |>

  # Significant FE rows: italic term
  tab_style(
    style     = cell_text(style = "italic"),
    locations = cells_body(
      columns = Term,
      rows    = p_sig == TRUE & row_type == "fe"
    )
  ) |>

  cols_hide(columns = c(p_sig, p.value, row_type)) |>

  # Footnotes
  tab_footnote(
    footnote  = md(paste0(
      "**Bold red** *p*: statistically significant (*p* < .05).  ",
      "**Amber rows**: body mass covariate (\u03b2_BM) from adjusted REML model ",
      "(*outcome* ~ Sport \u00d7 Week + BM + (1|Athlete)); ",
      "included only where BM improved fit (LRT *p* < .05, 4/5 outcomes).  ",
      "**Grey rows**: random effects \u2014 \u03c4\u2080 = between-athlete SD; ",
      "\u03c3 = within-athlete residual SD; ICC = proportion of variance from ",
      "stable between-athlete differences."
    )),
    locations = cells_column_labels(columns = Term)
  ) |>
  tab_footnote(
    footnote  = md("*SE* = Standard Error; *df* = Satterthwaite degrees of freedom."),
    locations = cells_column_labels(columns = p_fmt)
  ) |>
  tab_footnote(
    footnote  = paste0(
      "CMJ Peak Force, IMTP Peak Force, and Peak Power fitted in scaled units (\u00f71,000); ",
      "estimates back-transformed to original units (N and W)."
    ),
    locations = cells_row_groups(
      groups = c("CMJ Peak Force (N)", "IMTP Peak Force (N)", "Peak Power (W)")
    )
  ) |>
  tab_footnote(
    footnote  = cat_fn_dsi,
    locations = cells_row_groups(groups = "Dynamic Strength Index (DSI)")
  ) |>
  tab_footnote(
    footnote  = cat_fn_cmj,
    locations = cells_row_groups(groups = "CMJ Peak Force (N)")
  ) |>
  tab_footnote(
    footnote  = cat_fn_imtp,
    locations = cells_row_groups(groups = "IMTP Peak Force (N)")
  ) |>

  tab_source_note(md(paste0(
    "CMJ = Countermovement Jump; IMTP = Isometric Mid-Thigh Pull; ",
    "DSI = Dynamic Strength Index (CMJ PF \u00f7 IMTP PF); BM = Body Mass.  \n",
    "Models: **lme4** v", packageVersion("lme4"), " + ",
    "**lmerTest** v", packageVersion("lmerTest"), " in R ",
    paste(R.version$major, R.version$minor, sep = "."),
    ". Optimizer: BOBYQA (maxfun = 200,000). ",
    "Categorical week model *p*-values are unadjusted for multiple comparisons."
  ))) |>

  opt_table_font(font = list(google_font("Roboto"), default_fonts())) |>
  tab_options(
    table.font.size                    = px(12),
    heading.title.font.size            = px(14),
    heading.subtitle.font.size         = px(11),
    column_labels.font.weight          = "bold",
    row_group.font.weight              = "bold",
    table.border.top.color             = "black",
    table.border.top.width             = px(2),
    column_labels.border.bottom.color  = "black",
    column_labels.border.bottom.width  = px(1.5),
    table.border.bottom.color          = "black",
    table.border.bottom.width          = px(2),
    row_group.border.top.color         = "grey20",
    row_group.border.top.width         = px(1),
    row_group.border.bottom.color      = "grey20",
    row_group.border.bottom.width      = px(0.5),
    footnotes.font.size                = px(10),
    source_notes.font.size             = px(10),
    data_row.padding                   = px(5),
    column_labels.padding              = px(7)
  )

tbl_complete_path <- file.path(TBL_DIR, "Table4_Mixed_Models_Complete.html")
gtsave(tbl_complete_gt, tbl_complete_path)
cat(sprintf("  Saved: %s\n\n", basename(tbl_complete_path)))

cat("=== Done ===\n")
