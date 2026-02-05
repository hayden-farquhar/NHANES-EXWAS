# ============================================================================
# 05_sensitivity_analyses.R
# Robustness checks for validated findings
# ============================================================================

source("00_functions.R")

# --- 1. Load data ------------------------------------------------------------

load("validation_results.RData")

# Load primary analysis data (2017-2018)
if (file.exists("exwas_novelty_screen.RData")) {
  load("exwas_novelty_screen.RData")  # gives `dat`
  dat_full <- dat; rm(dat)
} else if (file.exists("exwas_expanded_results.RData")) {
  load("exwas_expanded_results.RData")  # gives `expanded`
  dat_full <- expanded; rm(expanded)
} else {
  stop("Cannot find primary analysis data.")
}
cat(sprintf("Primary data: %d obs x %d cols\n", nrow(dat_full), ncol(dat_full)))

# Select findings to test: validated or top FDR-significant
sens_findings <- validated %>%
  filter(validated | suggestive_val) %>%
  select(exposure, outcome, exp_label, out_label, type, beta_primary)

if (nrow(sens_findings) == 0) {
  cat("No validated findings; using top 15 FDR-significant.\n")
  sens_findings <- validated %>%
    arrange(fdr_primary) %>%
    head(15) %>%
    select(exposure, outcome, exp_label, out_label, type, beta_primary)
}

cat(sprintf("Sensitivity analyses for %d findings\n\n", nrow(sens_findings)))

# --- 2. Define sensitivity analysis specifications ----------------------------

run_sensitivity <- function(exposure, outcome, data, label,
                            covariates = COVARS_FULL,
                            weight_var = NULL,
                            subset_expr = NULL,
                            binary = FALSE) {
  # Apply subsetting if specified
  if (!is.null(subset_expr)) {
    sub_data <- data %>% filter(!!rlang::parse_expr(subset_expr))
  } else {
    sub_data <- data
  }

  if (is.null(weight_var)) {
    weight_var <- select_weight(exposure, sub_data)
  }

  res <- tryCatch(
    run_exwas_model(exposure, outcome, sub_data, covariates, weight_var, binary),
    error = function(e) NULL
  )

  if (is.null(res)) {
    return(data.frame(
      exposure = exposure, outcome = outcome,
      analysis = label, beta = NA, se = NA, p_value = NA, n = NA,
      stringsAsFactors = FALSE
    ))
  }

  data.frame(
    exposure = exposure, outcome = outcome,
    analysis = label,
    beta     = res$beta,
    se       = res$se,
    p_value  = res$p_value,
    n        = res$n,
    stringsAsFactors = FALSE
  )
}

# --- 3. Run all sensitivity analyses -----------------------------------------

all_sens <- list()

for (i in 1:nrow(sens_findings)) {
  exp_v  <- sens_findings$exposure[i]
  out_v  <- sens_findings$outcome[i]
  binary <- sens_findings$type[i] == "binary"

  if (!(exp_v %in% names(dat_full)) || !(out_v %in% names(dat_full))) next

  wt <- select_weight(exp_v, dat_full)

  cat(sprintf("\n--- %s -> %s ---\n", exp_v, out_v))

  # a) Primary model (reference)
  all_sens[[length(all_sens) + 1]] <- run_sensitivity(
    exp_v, out_v, dat_full, "Primary model",
    covariates = COVARS_FULL, weight_var = wt, binary = binary
  )

  # b) Stratify by sex: females only
  all_sens[[length(all_sens) + 1]] <- run_sensitivity(
    exp_v, out_v, dat_full, "Females only",
    covariates = gsub("\\+ female|female \\+", "", COVARS_FULL) %>%
      trimws() %>% gsub("\\s+\\+\\s+\\+", " + ", .),
    weight_var = wt, subset_expr = "female == 1", binary = binary
  )

  # c) Stratify by sex: males only
  all_sens[[length(all_sens) + 1]] <- run_sensitivity(
    exp_v, out_v, dat_full, "Males only",
    covariates = gsub("\\+ female|female \\+", "", COVARS_FULL) %>%
      trimws() %>% gsub("\\s+\\+\\s+\\+", " + ", .),
    weight_var = wt, subset_expr = "female == 0", binary = binary
  )

  # d) Stratify by age: <50
  all_sens[[length(all_sens) + 1]] <- run_sensitivity(
    exp_v, out_v, dat_full, "Age < 50",
    weight_var = wt, subset_expr = "RIDAGEYR < 50", binary = binary
  )

  # e) Stratify by age: >=50
  all_sens[[length(all_sens) + 1]] <- run_sensitivity(
    exp_v, out_v, dat_full, "Age >= 50",
    weight_var = wt, subset_expr = "RIDAGEYR >= 50", binary = binary
  )

  # f) Exclude extreme outliers (> 99th percentile of exposure)
  p99 <- quantile(dat_full[[exp_v]], 0.99, na.rm = TRUE)
  outlier_expr <- sprintf("%s <= %f", exp_v, p99)
  all_sens[[length(all_sens) + 1]] <- run_sensitivity(
    exp_v, out_v, dat_full, "Excl. outliers (>P99)",
    weight_var = wt, subset_expr = outlier_expr, binary = binary
  )

  # g) Additional covariate: education (if available)
  if ("DMDEDUC2" %in% names(dat_full)) {
    covars_edu <- paste(COVARS_FULL, "+ DMDEDUC2")
    all_sens[[length(all_sens) + 1]] <- run_sensitivity(
      exp_v, out_v, dat_full, "Adjust for education",
      covariates = covars_edu, weight_var = wt, binary = binary
    )
  }

  # h) Additional covariate: cotinine (continuous smoking biomarker)
  if ("LBXCOT" %in% names(dat_full)) {
    # Log-transform cotinine
    if (!"log_LBXCOT" %in% names(dat_full)) {
      cot_vals <- as.numeric(dat_full$LBXCOT)
      min_pos <- min(cot_vals[cot_vals > 0], na.rm = TRUE)
      dat_full$log_LBXCOT <- log(ifelse(cot_vals <= 0, min_pos / 2, cot_vals))
    }
    covars_cot <- paste(COVARS_FULL, "+ log_LBXCOT")
    # Replace 'smoker' with cotinine
    covars_cot <- gsub("smoker", "log_LBXCOT", COVARS_FULL)
    all_sens[[length(all_sens) + 1]] <- run_sensitivity(
      exp_v, out_v, dat_full, "Cotinine instead of smoking",
      covariates = covars_cot, weight_var = wt, binary = binary
    )
  }

  # i) Restrict to adults 20+ (exclude adolescents)
  all_sens[[length(all_sens) + 1]] <- run_sensitivity(
    exp_v, out_v, dat_full, "Adults 20+ only",
    weight_var = wt, subset_expr = "RIDAGEYR >= 20", binary = binary
  )
}

# --- 4. Compile results ------------------------------------------------------

sens_df <- bind_rows(all_sens) %>%
  left_join(
    sens_findings %>% select(exposure, outcome, beta_primary, exp_label, out_label),
    by = c("exposure", "outcome")
  ) %>%
  mutate(
    direction_match = sign(beta) == sign(beta_primary),
    pct_change      = ifelse(!is.na(beta_primary) & beta_primary != 0,
                             (beta - beta_primary) / abs(beta_primary) * 100, NA)
  )

# --- 5. Summary --------------------------------------------------------------

cat("\n=== SENSITIVITY ANALYSIS SUMMARY ===\n\n")

# For each finding, show whether it's robust across analyses
sens_summary <- sens_df %>%
  filter(!is.na(beta)) %>%
  group_by(exposure, outcome, exp_label, out_label) %>%
  summarise(
    n_analyses     = n(),
    n_dir_match    = sum(direction_match, na.rm = TRUE),
    n_nominal_sig  = sum(p_value < 0.05, na.rm = TRUE),
    n_suggestive   = sum(p_value < 0.10, na.rm = TRUE),
    median_pct_change = median(abs(pct_change), na.rm = TRUE),
    robust = n_dir_match == n_analyses & n_nominal_sig >= n_analyses * 0.5,
    .groups = "drop"
  ) %>%
  arrange(desc(robust), desc(n_nominal_sig))

cat("Robustness overview:\n")
as.data.frame(sens_summary) %>% print()

# Show full detail for each finding
for (exp_v in unique(sens_df$exposure)) {
  for (out_v in unique(sens_df$outcome[sens_df$exposure == exp_v])) {
    sub <- sens_df %>% filter(exposure == exp_v, outcome == out_v)
    if (nrow(sub) == 0) next
    cat(sprintf("\n--- %s -> %s ---\n",
                sub$exp_label[1], sub$out_label[1]))
    sub %>%
      select(analysis, beta, se, p_value, n, direction_match, pct_change) %>%
      print()
  }
}

# --- 6. Save -----------------------------------------------------------------

save(sens_df, sens_summary, file = "sensitivity_results.RData")
write.csv(sens_df, "sensitivity_results.csv", row.names = FALSE)
write.csv(sens_summary, "sensitivity_summary.csv", row.names = FALSE)

cat("\nSaved: sensitivity_results.RData, sensitivity_results.csv, sensitivity_summary.csv\n")
