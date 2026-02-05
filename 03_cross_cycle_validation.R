# ============================================================================
# 03_cross_cycle_validation.R
# Replicate FDR-significant findings from 2017-2018 in NHANES 2015-2016
# ============================================================================

source("00_functions.R")

# --- 1. Load data ------------------------------------------------------------

load("exwas_master_results.RData")        # master, priority_findings
load("nhanes_2015_2016_validation.RData") # val_data

cat(sprintf("Validation data: %d obs x %d cols\n", nrow(val_data), ncol(val_data)))

# --- 2. Identify findings to validate ---------------------------------------

# All FDR < 0.05 findings from primary analysis
to_validate <- master %>%
  filter(sig_fdr05) %>%
  select(exposure, outcome, beta, p_fdr_global, n, direction, exp_label, out_label,
         chem_class, out_domain, type) %>%
  rename(
    beta_primary   = beta,
    fdr_primary    = p_fdr_global,
    n_primary      = n,
    dir_primary    = direction
  )

cat(sprintf("\nFindings to validate: %d\n", nrow(to_validate)))

# --- 3. Run validation models ------------------------------------------------

cat("\n=== RUNNING VALIDATION MODELS ===\n\n")

val_results <- list()

for (i in 1:nrow(to_validate)) {
  exp_v  <- to_validate$exposure[i]
  out_v  <- to_validate$outcome[i]
  binary <- to_validate$type[i] == "binary"

  # Check variables exist in validation data
  if (!(exp_v %in% names(val_data))) {
    cat(sprintf("  SKIP %s -> %s: exposure not available in 2015-2016\n", exp_v, out_v))
    val_results[[i]] <- data.frame(
      exposure = exp_v, outcome = out_v,
      beta_val = NA, se_val = NA, p_val = NA, n_val = NA, df_val = NA,
      status = "exposure_missing"
    )
    next
  }
  if (!(out_v %in% names(val_data))) {
    cat(sprintf("  SKIP %s -> %s: outcome not available in 2015-2016\n", exp_v, out_v))
    val_results[[i]] <- data.frame(
      exposure = exp_v, outcome = out_v,
      beta_val = NA, se_val = NA, p_val = NA, n_val = NA, df_val = NA,
      status = "outcome_missing"
    )
    next
  }

  wt <- select_weight(exp_v, val_data)

  res <- tryCatch(
    run_exwas_model(
      exposure   = exp_v,
      outcome    = out_v,
      data       = val_data,
      covariates = COVARS_FULL,
      weight_var = wt,
      binary     = binary
    ),
    error = function(e) {
      cat(sprintf("  ERROR %s -> %s: %s\n", exp_v, out_v, e$message))
      NULL
    }
  )

  if (is.null(res)) {
    val_results[[i]] <- data.frame(
      exposure = exp_v, outcome = out_v,
      beta_val = NA, se_val = NA, p_val = NA, n_val = NA, df_val = NA,
      status = "model_failed"
    )
  } else {
    cat(sprintf("  %s -> %s: beta=%.3f, p=%.4f, n=%d\n",
                exp_v, out_v, res$beta, res$p_value, res$n))
    val_results[[i]] <- data.frame(
      exposure = exp_v, outcome = out_v,
      beta_val = res$beta, se_val = res$se,
      p_val = res$p_value, n_val = res$n, df_val = res$df,
      status = "ok"
    )
  }
}

val_results_df <- bind_rows(val_results)

# --- 4. Merge primary + validation results -----------------------------------

validated <- to_validate %>%
  left_join(val_results_df, by = c("exposure", "outcome")) %>%
  mutate(
    dir_val          = case_when(
      is.na(beta_val) ~ NA_character_,
      beta_val > 0    ~ "Positive",
      TRUE            ~ "Negative"
    ),
    direction_match  = dir_primary == dir_val,
    val_nominal      = !is.na(p_val) & p_val < 0.05,
    val_suggestive   = !is.na(p_val) & p_val < 0.10,
    validated        = direction_match & val_nominal,
    suggestive_val   = direction_match & val_suggestive & !val_nominal
  )

# --- 5. Summary --------------------------------------------------------------

cat("\n=== VALIDATION SUMMARY ===\n\n")

n_attempted <- sum(validated$status == "ok")
n_missing   <- sum(validated$status != "ok")
n_validated <- sum(validated$validated, na.rm = TRUE)
n_suggest   <- sum(validated$suggestive_val, na.rm = TRUE)
n_dir_match <- sum(validated$direction_match, na.rm = TRUE)

cat(sprintf("Attempted:                  %d / %d\n", n_attempted, nrow(validated)))
cat(sprintf("Exposure/outcome missing:   %d\n", n_missing))
cat(sprintf("Direction consistent:        %d / %d (%.0f%%)\n",
            n_dir_match, n_attempted, n_dir_match / n_attempted * 100))
cat(sprintf("Validated (p<0.05 + dir):   %d\n", n_validated))
cat(sprintf("Suggestive (p<0.10 + dir):  %d\n", n_suggest))

cat("\n--- VALIDATED FINDINGS (replicated, p<0.05, same direction) ---\n")
validated %>%
  filter(validated) %>%
  select(exp_label, out_label, beta_primary, beta_val, p_val, n_primary, n_val) %>%
  as.data.frame() %>%
  print()

cat("\n--- SUGGESTIVE (p<0.10, same direction) ---\n")
validated %>%
  filter(suggestive_val) %>%
  select(exp_label, out_label, beta_primary, beta_val, p_val) %>%
  as.data.frame() %>%
  print()

cat("\n--- NOT REPLICATED ---\n")
validated %>%
  filter(status == "ok", !validated, !suggestive_val) %>%
  select(exp_label, out_label, beta_primary, beta_val, p_val, direction_match) %>%
  as.data.frame() %>%
  print()

# --- 6. Save -----------------------------------------------------------------

save(validated, file = "validation_results.RData")
write.csv(validated, "validation_results.csv", row.names = FALSE)

cat("\nSaved: validation_results.RData, validation_results.csv\n")
