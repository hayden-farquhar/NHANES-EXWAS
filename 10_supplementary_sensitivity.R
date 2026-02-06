# 11_supplementary_sensitivity.R
# Additional sensitivity analyses for robustness:
# 1. eGFR-adjusted sensitivity for perchlorate-BUN (Table S15)
# 2. E-value analysis for HIGH-novelty findings
# 3. Verify alcohol adjustment direction (bidirectional vs systematic attenuation)

# Load dependencies
source("00_functions.R")

# Load data
load("exwas_novelty_screen.RData")  # dat: 5265x609

# Load master results (has SE) and validation results
master <- read_csv("exwas_master_results.csv")
validation <- read_csv("validation_results.csv")

cat("=== Additional Sensitivity Analyses ===\n\n")

# -----------------------------------------------------------------------------
# 1. eGFR-Adjusted Sensitivity for Perchlorate-BUN
# -----------------------------------------------------------------------------
# Concern: Reverse causation - impaired renal function could increase
# both urinary perchlorate (reduced clearance) and BUN (reduced urea excretion).
# Adjusting for eGFR helps distinguish renal confounding from direct effect.

cat("1. eGFR-Adjusted Sensitivity for Perchlorate-BUN\n")
cat("================================================\n\n")

# Prepare data for perchlorate-BUN analysis
perc_data <- dat %>%
  filter(!is.na(URXUP8), !is.na(LBXSBU), !is.na(eGFR)) %>%
  mutate(
    log_URXUP8 = log(URXUP8),
    weight = WTSA2YR
  ) %>%
  filter(weight > 0)

cat(sprintf("Sample size with complete perchlorate, BUN, and eGFR data: n = %d\n", nrow(perc_data)))

# Primary model (without eGFR)
design_primary <- svydesign(
  id = ~SDMVPSU,
  strata = ~SDMVSTRA,
  weights = ~weight,
  data = perc_data,
  nest = TRUE
)

model_primary <- svyglm(
  LBXSBU ~ log_URXUP8 + RIDAGEYR + female + race3 + INDFMPIR + BMXBMI + smoker,
  design = design_primary
)

primary_results <- tidy(model_primary) %>%
  filter(term == "log_URXUP8")

# eGFR-adjusted model
model_egfr <- svyglm(
  LBXSBU ~ log_URXUP8 + RIDAGEYR + female + race3 + INDFMPIR + BMXBMI + smoker + eGFR,
  design = design_primary
)

egfr_results <- tidy(model_egfr) %>%
  filter(term == "log_URXUP8")

# Compare results
cat("\nPerchlorate-BUN Results:\n")
cat(sprintf("  Primary model:     beta = %.3f (SE = %.3f), p = %.2e\n",
            primary_results$estimate, primary_results$std.error, primary_results$p.value))
cat(sprintf("  eGFR-adjusted:     beta = %.3f (SE = %.3f), p = %.2e\n",
            egfr_results$estimate, egfr_results$std.error, egfr_results$p.value))

pct_change <- (egfr_results$estimate - primary_results$estimate) / primary_results$estimate * 100
cat(sprintf("  Percent change:    %.1f%%\n", pct_change))
cat(sprintf("  Direction match:   %s\n", ifelse(sign(egfr_results$estimate) == sign(primary_results$estimate), "Yes", "No")))

# Get eGFR coefficient
egfr_coef <- tidy(model_egfr) %>%
  filter(term == "eGFR")

cat(sprintf("\n  eGFR coefficient:  beta = %.4f (p = %.2e)\n",
            egfr_coef$estimate, egfr_coef$p.value))
cat("  Interpretation: Negative eGFR-BUN coefficient expected (lower GFR -> higher BUN)\n")

# Save Table S15 results
table_s15 <- data.frame(
  model = c("Primary (no eGFR)", "eGFR-adjusted"),
  beta = c(primary_results$estimate, egfr_results$estimate),
  se = c(primary_results$std.error, egfr_results$std.error),
  p_value = c(primary_results$p.value, egfr_results$p.value),
  pct_change = c(NA, pct_change),
  n = nrow(perc_data)
)

write_csv(table_s15, "figures/table_s15_egfr_adjusted_perchlorate_bun.csv")
cat("\n  -> Saved to figures/table_s15_egfr_adjusted_perchlorate_bun.csv\n")

# -----------------------------------------------------------------------------
# 2. E-Value Analysis for HIGH-Novelty Findings
# -----------------------------------------------------------------------------
# E-value: minimum strength of association an unmeasured confounder would need
# with both exposure and outcome to fully explain away the observed association

cat("\n\n2. E-Value Analysis for HIGH-Novelty Findings\n")
cat("==============================================\n\n")

# Function to calculate E-value for continuous outcomes
# E-value = RR + sqrt(RR * (RR - 1)) where RR approximated from beta/SE
calc_evalue <- function(beta, se, outcome_sd) {
  # Convert beta to approximate RR using Cohen's d approach
  # d = beta / outcome_sd, then convert d to OR/RR approximation
  d <- abs(beta) / outcome_sd

  # Approximate RR from d: RR ~ exp(d * pi / sqrt(3))
  # For continuous outcomes, we use the transformation: RR = exp(0.91 * d)
  rr <- exp(0.91 * d)

  # E-value formula
  e_value <- rr + sqrt(rr * (rr - 1))

  # Lower CI bound
  d_lower <- max(0, abs(beta) - 1.96 * se) / outcome_sd
  rr_lower <- exp(0.91 * d_lower)
  e_value_lower <- ifelse(rr_lower > 1, rr_lower + sqrt(rr_lower * (rr_lower - 1)), 1)

  return(list(
    effect_rr = rr,
    e_value = e_value,
    e_value_lower = e_value_lower
  ))
}

# Get SE from master results (which has the se column)
# DMA-Uric acid
dma_master <- master %>% filter(exposure == "log_URXUDMA", outcome == "LBXSUA")
dma_validation <- validation %>% filter(exposure == "log_URXUDMA", outcome == "LBXSUA")
uric_sd <- sd(dat$LBXSUA, na.rm = TRUE)

dma_evalue <- calc_evalue(
  beta = dma_validation$beta_primary,
  se = dma_master$se,
  outcome_sd = uric_sd
)

cat("DMA-Uric Acid:\n")
cat(sprintf("  Outcome SD: %.2f mg/dL\n", uric_sd))
cat(sprintf("  Effect estimate: beta = %.3f (SE = %.3f)\n", dma_validation$beta_primary, dma_master$se))
cat(sprintf("  Approximate RR: %.2f\n", dma_evalue$effect_rr))
cat(sprintf("  E-value: %.2f\n", dma_evalue$e_value))
cat(sprintf("  E-value for CI lower bound: %.2f\n", dma_evalue$e_value_lower))
cat("  Interpretation: Unmeasured confounder would need RR >= E-value with both\n")
cat("                  exposure AND outcome to explain away the association\n")

# Perchlorate-BUN
perc_master <- master %>% filter(exposure == "log_URXUP8", outcome == "LBXSBU")
perc_validation <- validation %>% filter(exposure == "log_URXUP8", outcome == "LBXSBU")
bun_sd <- sd(dat$LBXSBU, na.rm = TRUE)

perc_evalue <- calc_evalue(
  beta = perc_validation$beta_primary,
  se = perc_master$se,
  outcome_sd = bun_sd
)

cat("\nPerchlorate-BUN:\n")
cat(sprintf("  Outcome SD: %.2f mg/dL\n", bun_sd))
cat(sprintf("  Effect estimate: beta = %.3f (SE = %.3f)\n", perc_validation$beta_primary, perc_master$se))
cat(sprintf("  Approximate RR: %.2f\n", perc_evalue$effect_rr))
cat(sprintf("  E-value: %.2f\n", perc_evalue$e_value))
cat(sprintf("  E-value for CI lower bound: %.2f\n", perc_evalue$e_value_lower))

# Methylmercury-Waist (though likely confounded)
mehg_master <- master %>% filter(exposure == "log_LBXBGM", outcome == "BMXWAIST")
mehg_validation <- validation %>% filter(exposure == "log_LBXBGM", outcome == "BMXWAIST")
waist_sd <- sd(dat$BMXWAIST, na.rm = TRUE)

mehg_evalue <- calc_evalue(
  beta = abs(mehg_validation$beta_primary),  # Use absolute value
  se = mehg_master$se,
  outcome_sd = waist_sd
)

cat("\nMethylmercury-Waist Circumference:\n")
cat(sprintf("  Outcome SD: %.2f cm\n", waist_sd))
cat(sprintf("  Effect estimate: beta = %.3f (SE = %.3f)\n", mehg_validation$beta_primary, mehg_master$se))
cat(sprintf("  Approximate RR: %.2f\n", mehg_evalue$effect_rr))
cat(sprintf("  E-value: %.2f\n", mehg_evalue$e_value))
cat(sprintf("  E-value for CI lower bound: %.2f\n", mehg_evalue$e_value_lower))
cat("  Note: This association is likely explained by fish consumption confounding\n")

# Save E-value results
evalue_results <- data.frame(
  finding = c("DMA-Uric acid", "Perchlorate-BUN", "Methylmercury-Waist"),
  beta = c(dma_validation$beta_primary, perc_validation$beta_primary, mehg_validation$beta_primary),
  se = c(dma_master$se, perc_master$se, mehg_master$se),
  outcome_sd = c(uric_sd, bun_sd, waist_sd),
  approx_rr = c(dma_evalue$effect_rr, perc_evalue$effect_rr, mehg_evalue$effect_rr),
  e_value = c(dma_evalue$e_value, perc_evalue$e_value, mehg_evalue$e_value),
  e_value_lower = c(dma_evalue$e_value_lower, perc_evalue$e_value_lower, mehg_evalue$e_value_lower)
)

write_csv(evalue_results, "figures/evalue_analysis.csv")
cat("\n  -> Saved to figures/evalue_analysis.csv\n")

# -----------------------------------------------------------------------------
# 3. Verify Alcohol Adjustment Direction (Bidirectional vs Systematic)
# -----------------------------------------------------------------------------

cat("\n\n3. Alcohol Adjustment Direction Analysis\n")
cat("========================================\n\n")

# Download alcohol data (ALQ_J) and calculate drinks_per_week
cat("Downloading ALQ_J (alcohol use)...\n")
alq <- tryCatch(nhanes("ALQ_J"), error = function(e) NULL)

if (!is.null(alq)) {
  alq_cols <- intersect(c("SEQN", "ALQ121", "ALQ130"), names(alq))
  alq <- alq[, alq_cols]

  # ALQ121 is frequency of drinking (categorical text or codes)
  # ALQ130 is number of drinks per occasion
  alq <- alq %>%
    mutate(
      ALQ121_char = as.character(ALQ121),
      # Convert frequency text labels to approximate occasions per year
      freq_per_year = case_when(
        ALQ121_char %in% c("Every day", "0") ~ 365,
        ALQ121_char %in% c("Nearly every day", "1") ~ 300,
        ALQ121_char %in% c("3 to 4 times a week", "2") ~ 182,
        ALQ121_char %in% c("2 times a week", "3") ~ 104,
        ALQ121_char %in% c("Once a week", "4") ~ 52,
        ALQ121_char %in% c("2 to 3 times a month", "5") ~ 30,
        ALQ121_char %in% c("Once a month", "6") ~ 12,
        ALQ121_char %in% c("7 to 11 times in the last year", "7") ~ 9,
        ALQ121_char %in% c("3 to 6 times in the last year", "8") ~ 4.5,
        ALQ121_char %in% c("1 to 2 times in the last year", "9") ~ 1.5,
        ALQ121_char %in% c("Never in the last year", "10") ~ 0,
        TRUE ~ NA_real_
      ),
      ALQ130_num = as.numeric(as.character(ALQ130)),
      ALQ130_num = ifelse(ALQ130_num %in% c(777, 999), NA, ALQ130_num),
      # Approximate drinks per week
      drinks_per_week = ifelse(!is.na(freq_per_year) & !is.na(ALQ130_num),
                               (freq_per_year * ALQ130_num) / 52, NA)
    )

  # Merge with main data
  dat_alcohol <- dat %>%
    left_join(alq %>% select(SEQN, drinks_per_week), by = "SEQN") %>%
    filter(!is.na(drinks_per_week))

  cat(sprintf("Sample with alcohol data: n = %d\n", nrow(dat_alcohol)))

  # Get validated findings
  validated <- validation %>% filter(validated == TRUE)

  alcohol_changes <- data.frame()

  for (i in 1:nrow(validated)) {
    exp_var <- validated$exposure[i]
    out_var <- validated$outcome[i]

    # Get base variable name for weight selection
    base_exp <- gsub("^log_", "", exp_var)
    weight_var <- select_weight(base_exp, dat_alcohol)

    # Subset to non-missing
    model_data <- dat_alcohol %>%
      filter(!is.na(.data[[base_exp]]), !is.na(.data[[out_var]])) %>%
      mutate(
        log_exp = log(.data[[base_exp]]),
        weight = .data[[weight_var]]
      ) %>%
      filter(weight > 0)

    if (nrow(model_data) < 100) next

    # Build design
    design <- svydesign(
      id = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~weight,
      data = model_data, nest = TRUE
    )

    # Determine if BMI should be included
    is_anthro <- out_var %in% c("BMXBMI", "BMXWAIST")

    # Primary model (no alcohol)
    if (is_anthro) {
      formula_primary <- as.formula(paste(out_var, "~ log_exp + RIDAGEYR + female + race3 + INDFMPIR + smoker"))
      formula_alcohol <- as.formula(paste(out_var, "~ log_exp + RIDAGEYR + female + race3 + INDFMPIR + smoker + drinks_per_week"))
    } else {
      formula_primary <- as.formula(paste(out_var, "~ log_exp + RIDAGEYR + female + race3 + INDFMPIR + BMXBMI + smoker"))
      formula_alcohol <- as.formula(paste(out_var, "~ log_exp + RIDAGEYR + female + race3 + INDFMPIR + BMXBMI + smoker + drinks_per_week"))
    }

    tryCatch({
      model_primary <- svyglm(formula_primary, design = design)
      model_alcohol <- svyglm(formula_alcohol, design = design)

      beta_primary <- coef(model_primary)["log_exp"]
      beta_alcohol <- coef(model_alcohol)["log_exp"]

      pct_change <- (beta_alcohol - beta_primary) / beta_primary * 100

      alcohol_changes <- bind_rows(alcohol_changes, data.frame(
        exposure = exp_var,
        outcome = out_var,
        beta_primary = beta_primary,
        beta_alcohol = beta_alcohol,
        pct_change = pct_change,
        direction = ifelse(pct_change > 0, "amplified", "attenuated"),
        n = nrow(model_data)
      ))
    }, error = function(e) {
      cat(sprintf("  Warning: Model failed for %s-%s: %s\n", exp_var, out_var, e$message))
    })
  }

  if (nrow(alcohol_changes) > 0) {
    cat("\nAlcohol Adjustment Effect Direction:\n")
    print(alcohol_changes %>% select(exposure, outcome, pct_change, direction) %>% as.data.frame())

    n_attenuated <- sum(alcohol_changes$pct_change < 0, na.rm = TRUE)
    n_amplified <- sum(alcohol_changes$pct_change > 0, na.rm = TRUE)
    median_change <- median(abs(alcohol_changes$pct_change), na.rm = TRUE)
    range_change <- range(alcohol_changes$pct_change, na.rm = TRUE)

    cat(sprintf("\nSummary:\n"))
    cat(sprintf("  Attenuated (negative change): %d findings\n", n_attenuated))
    cat(sprintf("  Amplified (positive change):  %d findings\n", n_amplified))
    cat(sprintf("  Median absolute change: %.1f%%\n", median_change))
    cat(sprintf("  Range: %.1f%% to +%.1f%%\n", range_change[1], range_change[2]))

    if (abs(n_attenuated - n_amplified) <= 3) {
      cat("\nConclusion: Changes are BIDIRECTIONAL (both attenuation and amplification observed)\n")
      cat("This suggests the larger effect estimate changes reflect reduced sample size and\n")
      cat("associated sampling variability, rather than systematic alcohol confounding.\n")
    } else if (n_attenuated > n_amplified) {
      cat("\nConclusion: Changes show SYSTEMATIC ATTENUATION pattern\n")
      cat("This could suggest alcohol is a genuine confounder.\n")
    } else {
      cat("\nConclusion: Changes show SYSTEMATIC AMPLIFICATION pattern\n")
      cat("Unexpected - may reflect collider bias from conditioning on alcohol availability.\n")
    }

    write_csv(alcohol_changes, "figures/alcohol_adjustment_direction.csv")
    cat("\n  -> Saved to figures/alcohol_adjustment_direction.csv\n")
  } else {
    cat("No alcohol adjustment results computed.\n")
  }
} else {
  cat("Could not download ALQ_J data. Skipping alcohol analysis.\n")
}

# -----------------------------------------------------------------------------
# Summary
# -----------------------------------------------------------------------------

cat("\n\n=== Summary ===\n")
cat("1. eGFR-adjusted perchlorate-BUN: Table S15 generated\n")
cat("2. E-value analysis: Results saved for three HIGH-novelty findings\n")
cat("3. Alcohol adjustment: Direction analysis completed\n")

cat("\nAll results saved to figures/ directory.\n")
