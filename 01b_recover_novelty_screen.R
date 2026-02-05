# ============================================================================
# 01b_recover_novelty_screen.R
# Recover the novelty-candidate ExWAS results that were lost when .RDataTmp
# became corrupted. Downloads additional NHANES 2017-2018 tables, merges
# with existing data, runs ExWAS, and saves results for consolidation.
# ============================================================================

source("00_functions.R")

# --- 1. Load existing dataset ------------------------------------------------

load("exwas_expanded_results.RData")  # gives us `expanded` (5265 x 502)
dat <- expanded
cat(sprintf("Base dataset: %d obs x %d vars\n", nrow(dat), ncol(dat)))

# --- 2. Download additional NHANES 2017-2018 tables --------------------------

additional_tables <- list(
  VOCWB  = "VOCWB_J",    # Whole blood VOCs (benzene = LBXVBZ)
  IHGEM  = "IHGEM_J",    # Inorganic/methylmercury (LBXBGM)
  UIO    = "UIO_J",       # Urinary iodine (URXUIO)
  PERNT  = "PERNT_J",     # Perchlorate/nitrate/thiocyanate (URXUP8)
  UAS    = "UAS_J",       # Urinary arsenic species (URXUAS, URXUDMA)
  SSGLYP = "SSGLYP_J",   # Glyphosate (surplus serum)
  UCOBALT = "UM_J"        # Urinary metals (might contain cobalt)
)

# Also try some alternative table names for cobalt
alt_tables <- list(
  COBALT2 = "PBCD_J"  # Blood metals panel (may include cobalt)
)

downloaded <- list()
for (nm in names(additional_tables)) {
  tbl_id <- additional_tables[[nm]]
  cat(sprintf("Downloading %s (%s)... ", nm, tbl_id))
  tryCatch({
    df <- nhanes(tbl_id)
    downloaded[[nm]] <- df
    cat(sprintf("OK (%d x %d)\n", nrow(df), ncol(df)))
  }, error = function(e) {
    cat(sprintf("FAILED: %s\n", e$message))
  })
}

# Try alternative tables
for (nm in names(alt_tables)) {
  tbl_id <- alt_tables[[nm]]
  cat(sprintf("Downloading %s (%s)... ", nm, tbl_id))
  tryCatch({
    df <- nhanes(tbl_id)
    downloaded[[nm]] <- df
    cat(sprintf("OK (%d x %d)\n", nrow(df), ncol(df)))
  }, error = function(e) {
    cat(sprintf("FAILED: %s\n", e$message))
  })
}

# --- 3. Merge new tables into existing data ----------------------------------

cat("\nMerging new tables...\n")
for (nm in names(downloaded)) {
  df <- downloaded[[nm]]
  if (!"SEQN" %in% names(df)) next

  # Only add NEW columns (avoid duplicates)
  new_cols <- setdiff(names(df), names(dat))
  if (length(new_cols) == 0) {
    cat(sprintf("  %s: no new columns\n", nm))
    next
  }
  df_sub <- df[, c("SEQN", new_cols), drop = FALSE]
  dat <- left_join(dat, df_sub, by = "SEQN")
  cat(sprintf("  + %s: added %d columns (now %d total)\n", nm, length(new_cols), ncol(dat)))
}

# --- 4. Identify target novelty chemicals ------------------------------------

# These are the chemicals from the R history's novelty assessment
target_chems <- c(
  "LBXVBZ",    # Blood benzene
  "LBXBGM",    # Methylmercury (blood)
  "LBXBCO",    # Blood cobalt
  "URXUIO",    # Urinary iodine
  "URXUP8",    # Urinary perchlorate
  "URXUAS",    # Total arsenic (urinary)
  "URXUDMA",   # Dimethylarsonic acid (urinary)
  "SSGLYP"     # Glyphosate (surplus serum)
)

cat("\n=== NOVELTY CHEMICAL AVAILABILITY ===\n")
available_chems <- c()
for (v in target_chems) {
  if (v %in% names(dat)) {
    vals <- as.numeric(dat[[v]])
    n_valid <- sum(!is.na(vals))
    n_pos <- sum(vals > 0, na.rm = TRUE)
    cat(sprintf("  %s: %d non-missing, %d positive (%.0f%%)\n",
                v, n_valid, n_pos, n_pos / max(n_valid, 1) * 100))
    if (n_valid >= 50 && n_pos / n_valid >= 0.5) {
      available_chems <- c(available_chems, v)
    }
  } else {
    cat(sprintf("  %s: NOT FOUND\n", v))
  }
}

# --- 5. Log-transform available chemicals ------------------------------------

cat("\nLog-transforming...\n")
new_log_vars <- c()
for (v in available_chems) {
  log_name <- paste0("log_", v)
  if (log_name %in% names(dat)) {
    cat(sprintf("  %s already exists\n", log_name))
    new_log_vars <- c(new_log_vars, log_name)
    next
  }
  vals <- as.numeric(dat[[v]])
  min_pos <- min(vals[vals > 0], na.rm = TRUE)
  vals_clean <- ifelse(!is.na(vals) & vals <= 0, min_pos / 2, vals)
  dat[[log_name]] <- log(vals_clean)
  new_log_vars <- c(new_log_vars, log_name)
  cat(sprintf("  Created %s\n", log_name))
}

# --- 6. Ensure covariates exist ----------------------------------------------

# Recreate race3 and smoker if missing
if (!"race3" %in% names(dat)) {
  dat$race3 <- case_when(
    dat$RIDRETH3 == "Non-Hispanic White" | dat$RIDRETH3 == 3 ~ "White",
    dat$RIDRETH3 == "Non-Hispanic Black" | dat$RIDRETH3 == 4 ~ "Black",
    TRUE ~ "Other"
  )
  dat$race3 <- factor(dat$race3, levels = c("White", "Black", "Other"))
}
if (!"smoker" %in% names(dat)) {
  if ("smoking_status" %in% names(dat)) {
    dat$smoker <- as.numeric(dat$smoking_status == "Current")
  } else if ("LBXCOT" %in% names(dat)) {
    dat$smoker <- as.numeric(as.numeric(dat$LBXCOT) > 10)
  }
}
if (!"female" %in% names(dat)) {
  dat$female <- as.numeric(dat$RIAGENDR == 2)
}

# --- 7. Define outcomes (matching original screen) ---------------------------

continuous_outcomes <- c("LBXSATSI", "LBXSASSI", "LBXSGB", "LBXSAPSI", "LBXSC3SI",
                         "LBXSBU", "LBXSKSI", "LBXSNASI", "LBXSCLSI", "LBXSPH",
                         "LBXSTB", "LBXSTP", "LBXSUA", "LBXSTR", "LBXSCH",
                         "LBXTC", "LBDHDD", "LBXGH", "LBXGLU", "LBXHSCRP",
                         "LBXWBCSI", "LBXRBCSI", "LBXHGB", "LBXPLTSI",
                         "mean_sys", "mean_dia", "BMXBMI", "BMXWAIST",
                         "PHQ9_total", "SLD012", "eGFR")
binary_outcomes <- c("depression_binary")
all_outcomes <- c(continuous_outcomes, binary_outcomes)
avail_outcomes <- all_outcomes[all_outcomes %in% names(dat)]
cat(sprintf("\nAvailable outcomes: %d\n", length(avail_outcomes)))

# --- 8. Run ExWAS for novelty chemicals --------------------------------------

cat(sprintf("\n=== RUNNING ExWAS: %d NOVELTY EXPOSURES x %d OUTCOMES ===\n\n",
            length(new_log_vars), length(avail_outcomes)))

test_grid <- expand.grid(
  exposure = new_log_vars,
  outcome  = avail_outcomes,
  stringsAsFactors = FALSE
)
cat(sprintf("Total tests: %d\n", nrow(test_grid)))

t0 <- Sys.time()
results_novelty <- list()
n_done <- 0
n_fail <- 0

for (i in 1:nrow(test_grid)) {
  exp_v <- test_grid$exposure[i]
  out_v <- test_grid$outcome[i]
  binary <- (out_v == "depression_binary")
  wt <- select_weight(exp_v, dat)
  covars <- if (out_v %in% c("BMXBMI", "BMXWAIST")) COVARS_NO_BMI else COVARS_FULL

  res <- tryCatch(
    run_exwas_model(exp_v, out_v, dat, covars, wt, binary),
    error = function(e) NULL
  )

  if (!is.null(res)) {
    results_novelty[[length(results_novelty) + 1]] <- res
    n_done <- n_done + 1
  } else {
    n_fail <- n_fail + 1
  }

  if (i %% 100 == 0) {
    cat(sprintf("  %d/%d done (%.0fs)\n", i, nrow(test_grid),
                as.numeric(Sys.time() - t0, units = "secs")))
  }
}

elapsed <- as.numeric(Sys.time() - t0, units = "secs")
cat(sprintf("\nCompleted: %d successful, %d failed (%.0fs)\n", n_done, n_fail, elapsed))

# --- 9. Compile and apply FDR ------------------------------------------------

novelty_df <- bind_rows(results_novelty)

# Check for LOD pile-up: flag chemicals with >70% at single value
cat("\n=== LOD CHECK ===\n")
lod_problem <- c()
for (v in available_chems) {
  vals <- as.numeric(dat[[v]])
  vals <- vals[!is.na(vals)]
  if (length(vals) < 50) next
  pct_at_min <- mean(vals == min(vals)) * 100
  cat(sprintf("  %s: %.0f%% at minimum value\n", v, pct_at_min))
  if (pct_at_min > 70) lod_problem <- c(lod_problem, v)
}

if (length(lod_problem) > 0) {
  cat(sprintf("\nRemoving LOD-problematic chemicals: %s\n", paste(lod_problem, collapse = ", ")))
  lod_log <- paste0("log_", lod_problem)
  novelty_df <- novelty_df %>% filter(!(exposure %in% lod_log))
}

# Filter implausibly large effects
novelty_df <- novelty_df %>%
  filter(
    (outcome != "depression_binary" & abs(beta) < 50) |
    (outcome == "depression_binary" & abs(beta) < 5)
  )

# FDR correction (within this novelty screen)
novelty_df <- novelty_df %>%
  mutate(
    p_fdr     = p.adjust(p_value, method = "BH"),
    sig_fdr05 = p_fdr < 0.05,
    sig_fdr01 = p_fdr < 0.01
  ) %>%
  arrange(p_value)

cat(sprintf("\n=== NOVELTY SCREEN RESULTS ===\n"))
cat(sprintf("Total tests (after filtering): %d\n", nrow(novelty_df)))
cat(sprintf("Nominal p<0.05: %d (%.1f%%)\n",
            sum(novelty_df$p_value < 0.05), mean(novelty_df$p_value < 0.05) * 100))
cat(sprintf("FDR<0.05: %d\n", sum(novelty_df$sig_fdr05)))
cat(sprintf("FDR<0.01: %d\n", sum(novelty_df$sig_fdr01)))

# --- 10. Display top findings ------------------------------------------------

cat("\nTop 30 hits:\n")
novelty_df %>%
  head(30) %>%
  mutate(
    exp_name = coalesce(CHEM_LABELS[exposure], exposure),
    out_name = coalesce(OUTCOME_LABELS[outcome], outcome)
  ) %>%
  select(exp_name, out_name, beta, p_value, p_fdr, n, df) %>%
  as.data.frame() %>%
  print()

if (sum(novelty_df$sig_fdr05) > 0) {
  cat("\n=== FDR-SIGNIFICANT NOVELTY FINDINGS ===\n")
  novelty_df %>%
    filter(sig_fdr05) %>%
    mutate(
      exp_name = coalesce(CHEM_LABELS[exposure], exposure),
      out_name = coalesce(OUTCOME_LABELS[outcome], outcome),
      direction = ifelse(beta > 0, "Positive", "Negative")
    ) %>%
    select(exp_name, out_name, beta, p_value, p_fdr, n, direction) %>%
    as.data.frame() %>%
    print()
}

# --- 11. Save ----------------------------------------------------------------

# Save as CSV for consolidation script
write.csv(novelty_df, "exwas_novelty_screen_results.csv", row.names = FALSE)

# Save full objects
save(novelty_df, dat, available_chems, new_log_vars,
     file = "exwas_novelty_screen.RData")

cat(sprintf("\nSaved: exwas_novelty_screen_results.csv (%d rows)\n", nrow(novelty_df)))
cat("Saved: exwas_novelty_screen.RData\n")
cat("\nRun 01_consolidate_results.R next to incorporate these into the master table.\n")
