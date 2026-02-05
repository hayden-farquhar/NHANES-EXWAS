# ============================================================================
# 10_reviewer_analyses.R
# Analyses addressing specific reviewer comments
# Run after scripts 00-09. Requires existing .RData files.
# ============================================================================

source("00_functions.R")

dir.create("figures", showWarnings = FALSE)

# --- 0. Load existing results and data ----------------------------------------

load("validation_results.RData")        # validated
load("exwas_master_results.RData")      # master, priority_findings
load("additional_analyses_results.RData")  # dat_full (has fish, PA, alcohol, creatinine)

cat(sprintf("Primary data: %d obs x %d cols\n", nrow(dat_full), ncol(dat_full)))

# Validated findings for additional sensitivity
val_findings <- validated %>%
  filter(validated) %>%
  select(exposure, outcome, exp_label, out_label, type, beta_primary)

cat(sprintf("Validated findings: %d\n\n", nrow(val_findings)))


# ============================================================================
# PART 1: Download 24-hour dietary recall data for protein intake
# ============================================================================

cat("=== PART 1: Downloading dietary data for protein intake ===\n\n")

# Helper: download with retry
nhanes_retry <- function(table_name, max_tries = 3) {
  for (attempt in 1:max_tries) {
    result <- tryCatch(nhanes(table_name), error = function(e) {
      cat(sprintf("  Attempt %d/%d failed: %s\n", attempt, max_tries, e$message))
      NULL
    })
    if (!is.null(result)) return(result)
    if (attempt < max_tries) Sys.sleep(2)
  }
  return(NULL)
}

# Download 24-hour dietary recall (Day 1 totals)
cat("Downloading DR1TOT_J (24-hour dietary recall totals)...\n")
diet <- nhanes_retry("DR1TOT_J")

if (!is.null(diet)) {
  # DR1TPROT = Total protein (gm)
  # DR1TKCAL = Total energy (kcal)
  if ("DR1TPROT" %in% names(diet)) {
    diet <- diet %>%
      transmute(
        SEQN,
        DR1TPROT = as.numeric(as.character(DR1TPROT)),
        DR1TKCAL = as.numeric(as.character(DR1TKCAL))
      ) %>%
      mutate(
        protein_g = ifelse(DR1TPROT < 0, NA, DR1TPROT),  # -1 = missing
        energy_kcal = ifelse(DR1TKCAL < 0, NA, DR1TKCAL),
        # Protein density (g per 1000 kcal) - accounts for total intake
        protein_density = ifelse(!is.na(protein_g) & !is.na(energy_kcal) & energy_kcal > 0,
                                 protein_g / energy_kcal * 1000, NA)
      )

    dat_full <- dat_full %>%
      left_join(diet %>% select(SEQN, protein_g, energy_kcal, protein_density),
                by = "SEQN")

    n_prot <- sum(!is.na(dat_full$protein_g))
    cat(sprintf("  Protein intake: %d non-missing (mean = %.1f g/day)\n",
                n_prot, mean(dat_full$protein_g, na.rm = TRUE)))
    cat(sprintf("  Protein density: mean = %.1f g/1000 kcal\n",
                mean(dat_full$protein_density, na.rm = TRUE)))
  } else {
    cat("  Warning: DR1TPROT not found in DR1TOT_J\n")
  }
} else {
  cat("  Skipping dietary data (download failed)\n")
}


# ============================================================================
# PART 2: Protein intake sensitivity for DMA-uric acid and perchlorate-BUN
# ============================================================================

cat("\n=== PART 2: Protein intake sensitivity analysis ===\n\n")

# Re-define run_sensitivity
run_sensitivity <- function(exposure, outcome, data, label,
                            covariates = COVARS_FULL,
                            weight_var = NULL,
                            subset_expr = NULL,
                            binary = FALSE) {
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

protein_sens <- list()

# DMA -> uric acid with protein adjustment
if ("protein_g" %in% names(dat_full)) {
  cat("Running protein-adjusted models for HIGH-novelty findings...\n")

  # a) DMA-uric acid
  covars_prot <- paste(COVARS_FULL, "+ protein_g")
  protein_sens[[1]] <- run_sensitivity(
    "log_URXUDMA", "LBXSUA", dat_full, "Adjust for protein intake",
    covariates = covars_prot, weight_var = "WTSA2YR", binary = FALSE
  )

  # b) Perchlorate-BUN
  protein_sens[[2]] <- run_sensitivity(
    "log_URXUP8", "LBXSBU", dat_full, "Adjust for protein intake",
    covariates = covars_prot, weight_var = "WTSA2YR", binary = FALSE
  )

  # c) Also with protein density (g/1000 kcal) to account for total intake
  covars_dens <- paste(COVARS_FULL, "+ protein_density")
  protein_sens[[3]] <- run_sensitivity(
    "log_URXUDMA", "LBXSUA", dat_full, "Adjust for protein density",
    covariates = covars_dens, weight_var = "WTSA2YR", binary = FALSE
  )

  protein_sens[[4]] <- run_sensitivity(
    "log_URXUP8", "LBXSBU", dat_full, "Adjust for protein density",
    covariates = covars_dens, weight_var = "WTSA2YR", binary = FALSE
  )
}

protein_sens_df <- bind_rows(protein_sens) %>%
  mutate(
    exp_label = case_when(
      exposure == "log_URXUDMA" ~ "Dimethylarsonic acid",
      exposure == "log_URXUP8" ~ "Urinary perchlorate",
      TRUE ~ exposure
    ),
    out_label = case_when(
      outcome == "LBXSUA" ~ "Uric acid",
      outcome == "LBXSBU" ~ "BUN",
      TRUE ~ outcome
    )
  )

# Get primary betas for comparison
primary_dma_ua <- master %>%
  filter(exposure == "log_URXUDMA", outcome == "LBXSUA") %>%
  pull(beta)
primary_perc_bun <- master %>%
  filter(exposure == "log_URXUP8", outcome == "LBXSBU") %>%
  pull(beta)

protein_sens_df <- protein_sens_df %>%
  mutate(
    beta_primary = case_when(
      exposure == "log_URXUDMA" ~ primary_dma_ua,
      exposure == "log_URXUP8" ~ primary_perc_bun
    ),
    pct_change = (beta - beta_primary) / abs(beta_primary) * 100,
    direction_match = sign(beta) == sign(beta_primary)
  )

cat("\n=== PROTEIN INTAKE SENSITIVITY RESULTS ===\n")
protein_sens_df %>%
  transmute(
    Chemical = exp_label,
    Outcome = out_label,
    Analysis = analysis,
    Beta_primary = round(beta_primary, 3),
    Beta_adj = round(beta, 3),
    `%Change` = round(pct_change, 1),
    P = signif(p_value, 3),
    N = n,
    Robust = ifelse(direction_match & p_value < 0.05, "Yes", "No")
  ) %>%
  as.data.frame() %>%
  print()


# ============================================================================
# PART 3: Full 6-level race/ethnicity sensitivity analysis
# ============================================================================

cat("\n=== PART 3: 6-level race/ethnicity sensitivity ===\n\n")

# RIDRETH3 codes:
# 1 = Mexican American
# 2 = Other Hispanic
# 3 = Non-Hispanic White
# 4 = Non-Hispanic Black
# 6 = Non-Hispanic Asian
# 7 = Other Race - Including Multi-Racial

# Create 6-level race factor from RIDRETH3 if not already done
if (!"race6" %in% names(dat_full)) {
  # First check current coding
  if ("RIDRETH3" %in% names(dat_full)) {
    dat_full$RIDRETH3_num <- as.numeric(as.character(dat_full$RIDRETH3))

    dat_full <- dat_full %>%
      mutate(
        race6 = case_when(
          RIDRETH3_num == 1 ~ "Mexican American",
          RIDRETH3_num == 2 ~ "Other Hispanic",
          RIDRETH3_num == 3 ~ "Non-Hispanic White",
          RIDRETH3_num == 4 ~ "Non-Hispanic Black",
          RIDRETH3_num == 6 ~ "Non-Hispanic Asian",
          RIDRETH3_num == 7 ~ "Other/Multiracial",
          TRUE ~ NA_character_
        ),
        race6 = factor(race6, levels = c("Non-Hispanic White", "Non-Hispanic Black",
                                          "Mexican American", "Other Hispanic",
                                          "Non-Hispanic Asian", "Other/Multiracial"))
      )

    cat(sprintf("6-level race distribution:\n"))
    print(table(dat_full$race6, useNA = "ifany"))
  }
}

# Define covariate string with 6-level race
COVARS_RACE6 <- "RIDAGEYR + female + race6 + INDFMPIR + BMXBMI + smoker"

# Run for all blood biomarker validated findings (not urinary - those have separate weights)
blood_findings <- val_findings %>%
  filter(grepl("^log_LBX", exposure))

race6_sens <- list()

for (i in 1:nrow(blood_findings)) {
  exp_v <- blood_findings$exposure[i]
  out_v <- blood_findings$outcome[i]
  binary <- isTRUE(blood_findings$type[i] == "binary")

  if (!(exp_v %in% names(dat_full)) || !(out_v %in% names(dat_full))) next

  cat(sprintf("  6-level race: %s -> %s\n",
              blood_findings$exp_label[i], blood_findings$out_label[i]))

  # Adjust covariates if outcome is BMI/waist
  covars_use <- COVARS_RACE6
  if (out_v %in% c("BMXBMI", "BMXWAIST")) {
    covars_use <- gsub("\\+ BMXBMI|BMXBMI \\+", "", covars_use)
    covars_use <- trimws(gsub("\\s+\\+\\s+\\+", " + ", covars_use))
  }

  race6_sens[[length(race6_sens) + 1]] <- run_sensitivity(
    exp_v, out_v, dat_full, "6-level race/ethnicity",
    covariates = covars_use, weight_var = "WTMEC2YR", binary = binary
  )
}

race6_sens_df <- bind_rows(race6_sens) %>%
  left_join(
    blood_findings %>% select(exposure, outcome, beta_primary, exp_label, out_label),
    by = c("exposure", "outcome")
  ) %>%
  mutate(
    direction_match = sign(beta) == sign(beta_primary),
    pct_change = (beta - beta_primary) / abs(beta_primary) * 100
  )

cat("\n=== 6-LEVEL RACE/ETHNICITY SENSITIVITY RESULTS ===\n")
race6_sens_df %>%
  filter(!is.na(beta)) %>%
  transmute(
    Chemical = exp_label,
    Outcome = out_label,
    Beta_3level = round(beta_primary, 3),
    Beta_6level = round(beta, 3),
    `%Change` = round(pct_change, 1),
    P = signif(p_value, 3),
    N = n,
    Robust = ifelse(direction_match & p_value < 0.05, "Yes", "No")
  ) %>%
  as.data.frame() %>%
  print()


# ============================================================================
# PART 4: Within-round FDR correction
# ============================================================================

cat("\n=== PART 4: Within-round FDR correction ===\n\n")

# Apply FDR separately within each screening round
within_fdr <- master %>%
  group_by(screen) %>%
  mutate(
    p_fdr_within = p.adjust(p_value, method = "BH"),
    sig_fdr_within = p_fdr_within < 0.05
  ) %>%
  ungroup()

# Summary by round
fdr_summary <- within_fdr %>%
  group_by(screen) %>%
  summarise(
    n_tests = n(),
    n_sig_global = sum(sig_fdr05, na.rm = TRUE),
    n_sig_within = sum(sig_fdr_within, na.rm = TRUE),
    .groups = "drop"
  )

cat("FDR correction comparison (global vs within-round):\n\n")
as.data.frame(fdr_summary) %>% print()

# Which findings lose significance under within-round FDR?
global_only <- within_fdr %>%
  filter(sig_fdr05 & !sig_fdr_within) %>%
  select(exposure, outcome, exp_label, out_label, screen,
         p_value, p_fdr_global, p_fdr_within)

cat("\nFindings significant under global FDR but NOT within-round FDR:\n")
if (nrow(global_only) > 0) {
  global_only %>%
    mutate(
      p_global = signif(p_fdr_global, 3),
      p_within = signif(p_fdr_within, 3)
    ) %>%
    select(exp_label, out_label, screen, p_global, p_within) %>%
    as.data.frame() %>%
    print()
} else {
  cat("  None - all globally-significant findings remain significant under within-round FDR\n")
}

# Export within-round FDR results
write.csv(
  within_fdr %>%
    filter(sig_fdr_within) %>%
    select(screen, exposure, outcome, exp_label, out_label, beta, p_value,
           p_fdr_global, p_fdr_within),
  "within_round_fdr_results.csv",
  row.names = FALSE
)


# ============================================================================
# PART 5: LOD 40-70% chemical analysis
# ============================================================================

cat("\n=== PART 5: Chemicals with 40-70% detection (LOD analysis) ===\n\n")

# Identify chemicals by their detection frequency
# Need to check LOD comment codes (LC suffix variables)

# Get list of unique chemical exposures
chem_vars <- unique(gsub("^log_", "", master$exposure))

lod_analysis <- data.frame(
  variable = character(),
  n_total = integer(),
  n_above_lod = integer(),
  pct_detect = numeric(),
  stringsAsFactors = FALSE
)

for (chem in chem_vars) {
  if (!(chem %in% names(dat_full))) next

  # Check for corresponding LOD comment code variable
  lc_var <- sub("LBX|URX|SS", "LC", chem)  # e.g., LBXBPB -> LCBPB

  vals <- dat_full[[chem]]
  n_total <- sum(!is.na(vals))

  if (lc_var %in% names(dat_full)) {
    # LC codes: typically 0 = at/below LOD, 1 = above LOD
    lc <- as.numeric(as.character(dat_full[[lc_var]]))
    n_above <- sum(lc == 1, na.rm = TRUE)
  } else {
    # If no LC variable, estimate based on value clustering at minimum
    min_val <- min(vals, na.rm = TRUE)
    n_at_min <- sum(abs(vals - min_val) < 1e-6, na.rm = TRUE)
    n_above <- n_total - n_at_min
  }

  if (n_total > 0) {
    lod_analysis <- rbind(lod_analysis, data.frame(
      variable = chem,
      n_total = n_total,
      n_above_lod = n_above,
      pct_detect = 100 * n_above / n_total,
      stringsAsFactors = FALSE
    ))
  }
}

# Filter to 40-70% detection range
lod_40_70 <- lod_analysis %>%
  filter(pct_detect >= 40, pct_detect <= 70) %>%
  arrange(pct_detect)

cat(sprintf("Chemicals with 40-70%% detection frequency: %d\n", nrow(lod_40_70)))

if (nrow(lod_40_70) > 0) {
  cat("\nChemicals in 40-70% detection range:\n")
  lod_40_70 %>%
    left_join(
      data.frame(
        variable = gsub("^log_", "", master$exposure),
        exp_label = master$exp_label
      ) %>% distinct(),
      by = "variable"
    ) %>%
    transmute(
      Variable = variable,
      Label = coalesce(exp_label, variable),
      `% Detected` = round(pct_detect, 1),
      N = n_total
    ) %>%
    as.data.frame() %>%
    print()
}

# How many chemicals below 40% (our threshold was 30% for inclusion)
n_below_40 <- sum(lod_analysis$pct_detect < 40)
n_below_30 <- sum(lod_analysis$pct_detect < 30)
cat(sprintf("\nChemicals with <40%% detection: %d\n", n_below_40))
cat(sprintf("Chemicals with <30%% detection (excluded): %d\n", n_below_30))


# ============================================================================
# PART 6: Generate supplementary tables S6-S11
# ============================================================================

cat("\n=== PART 6: Generating supplementary tables ===\n\n")

# Table S6: Within-round FDR correction comparison
table_s6 <- fdr_summary %>%
  rename(
    `Screening Round` = screen,
    `N Tests` = n_tests,
    `N Sig (Global FDR)` = n_sig_global,
    `N Sig (Within-Round FDR)` = n_sig_within
  )
write.csv(table_s6, "figures/table_s6_within_round_fdr.csv", row.names = FALSE)
cat("  Saved: figures/table_s6_within_round_fdr.csv\n")

# Table S7: 6-level race/ethnicity sensitivity
table_s7 <- race6_sens_df %>%
  filter(!is.na(beta)) %>%
  transmute(
    Chemical = exp_label,
    Outcome = out_label,
    `β (3-level race)` = round(beta_primary, 3),
    `β (6-level race)` = round(beta, 3),
    `% Change` = round(pct_change, 1),
    `P-value` = signif(p_value, 3),
    N = n
  )
write.csv(table_s7, "figures/table_s7_race6_sensitivity.csv", row.names = FALSE)
cat("  Saved: figures/table_s7_race6_sensitivity.csv\n")

# Table S8: Protein intake adjustment
table_s8 <- protein_sens_df %>%
  filter(!is.na(beta)) %>%
  transmute(
    Chemical = exp_label,
    Outcome = out_label,
    Adjustment = analysis,
    `β (Primary)` = round(beta_primary, 3),
    `β (Adjusted)` = round(beta, 3),
    `% Change` = round(pct_change, 1),
    `P-value` = signif(p_value, 3),
    N = n,
    Robust = ifelse(direction_match & p_value < 0.05, "Yes", "No")
  )
write.csv(table_s8, "figures/table_s8_protein_sensitivity.csv", row.names = FALSE)
cat("  Saved: figures/table_s8_protein_sensitivity.csv\n")

# Table S9: LOD 40-70% chemicals
table_s9 <- lod_40_70 %>%
  left_join(
    data.frame(
      variable = gsub("^log_", "", master$exposure),
      exp_label = master$exp_label
    ) %>% distinct(),
    by = "variable"
  ) %>%
  transmute(
    Chemical = coalesce(exp_label, variable),
    Variable = variable,
    `% Detected` = round(pct_detect, 1),
    `N Samples` = n_total
  )
write.csv(table_s9, "figures/table_s9_lod_40_70.csv", row.names = FALSE)
cat("  Saved: figures/table_s9_lod_40_70.csv\n")

# Table S10: Full power analysis (from script 09)
# Copy existing power_analysis.csv as table_s10
if (file.exists("power_analysis.csv")) {
  file.copy("power_analysis.csv", "figures/table_s10_power_analysis.csv", overwrite = TRUE)
  cat("  Saved: figures/table_s10_power_analysis.csv\n")
}

# Table S11: STROBE checklist
strobe_checklist <- tribble(
  ~Item, ~`Checklist Item`, ~`Manuscript Section`,

  "Title and abstract", "1a. Indicate the study's design with a commonly used term in the title or the abstract", "Title, Abstract",
  "", "1b. Provide in the abstract an informative and balanced summary of what was done and what was found", "Abstract",

  "Introduction", "2. Explain the scientific background and rationale for the investigation", "Introduction ¶1-2",
  "", "3. State specific objectives, including any prespecified hypotheses", "Introduction ¶3",

  "Methods", "4. Present key elements of study design early in the paper", "Methods 2.1",
  "", "5. Describe the setting, locations, and relevant dates", "Methods 2.1",
  "", "6. Give the eligibility criteria, and the sources and methods of selection of participants", "Methods 2.1",
  "", "7. Clearly define all outcomes, exposures, predictors, potential confounders, and effect modifiers", "Methods 2.2-2.3",
  "", "8. For each variable of interest, give sources of data and details of methods of assessment", "Methods 2.2-2.3",
  "", "9. Describe any efforts to address potential sources of bias", "Methods 2.3, 2.4.1",
  "", "10. Explain how the study size was arrived at", "Methods 2.1",
  "", "11. Explain how quantitative variables were handled in the analyses", "Methods 2.3",
  "", "12. Describe all statistical methods", "Methods 2.4",

  "Results", "13. Report numbers of individuals at each stage of study", "Table 1, Results 3.1",
  "", "14. Give characteristics of study participants", "Table 1",
  "", "15. Report numbers of outcome events or summary measures", "Results 3.2-3.5",
  "", "16. Give unadjusted estimates and, if applicable, confounder-adjusted estimates", "Table 2, Figures",
  "", "17. Report other analyses performed", "Results 3.4-3.5",

  "Discussion", "18. Summarise key results with reference to study objectives", "Discussion ¶1",
  "", "19. Discuss limitations of the study, taking into account sources of potential bias", "Discussion - Limitations",
  "", "20. Give a cautious overall interpretation of results considering objectives, limitations, multiplicity of analyses", "Discussion ¶6-7",
  "", "21. Discuss the generalizability of the results", "Discussion ¶2-3",
  "", "22. Give the source of funding and the role of the funders", "Not applicable (unfunded study)"
)

write.csv(strobe_checklist, "figures/table_s11_strobe_checklist.csv", row.names = FALSE)
cat("  Saved: figures/table_s11_strobe_checklist.csv\n")


# ============================================================================
# PART 7: Summary and save
# ============================================================================

cat("\n=== SUMMARY OF REVIEWER-REQUESTED ANALYSES ===\n\n")

cat("1. PROTEIN INTAKE SENSITIVITY:\n")
cat(sprintf("   DMA-uric acid: β changed by %.1f%% with protein adjustment (p=%.3f)\n",
            protein_sens_df$pct_change[protein_sens_df$exposure == "log_URXUDMA" &
                                         protein_sens_df$analysis == "Adjust for protein intake"],
            protein_sens_df$p_value[protein_sens_df$exposure == "log_URXUDMA" &
                                      protein_sens_df$analysis == "Adjust for protein intake"]))
cat(sprintf("   Perchlorate-BUN: β changed by %.1f%% with protein adjustment (p=%.3f)\n",
            protein_sens_df$pct_change[protein_sens_df$exposure == "log_URXUP8" &
                                         protein_sens_df$analysis == "Adjust for protein intake"],
            protein_sens_df$p_value[protein_sens_df$exposure == "log_URXUP8" &
                                      protein_sens_df$analysis == "Adjust for protein intake"]))

cat("\n2. 6-LEVEL RACE/ETHNICITY:\n")
n_robust_race6 <- sum(race6_sens_df$direction_match & race6_sens_df$p_value < 0.05, na.rm = TRUE)
cat(sprintf("   %d/%d blood biomarker findings robust with 6-level race\n",
            n_robust_race6, nrow(race6_sens_df)))

cat("\n3. WITHIN-ROUND FDR:\n")
cat(sprintf("   Global FDR significant: %d findings\n", sum(fdr_summary$n_sig_global)))
cat(sprintf("   Within-round FDR significant: %d findings\n", sum(fdr_summary$n_sig_within)))
cat(sprintf("   Findings losing significance: %d\n", nrow(global_only)))

cat("\n4. LOD ANALYSIS:\n")
cat(sprintf("   Chemicals with 40-70%% detection: %d\n", nrow(lod_40_70)))
cat(sprintf("   Chemicals excluded (<30%% detection): %d\n", n_below_30))

cat("\n")

# Save all results
save(protein_sens_df, race6_sens_df, within_fdr, fdr_summary, global_only,
     lod_analysis, lod_40_70, strobe_checklist, dat_full,
     file = "reviewer_analyses_results.RData")

cat("Saved: reviewer_analyses_results.RData\n")
cat("Saved: within_round_fdr_results.csv\n")
cat("\nSupplementary tables saved to figures/:\n")
cat("  table_s6_within_round_fdr.csv\n")
cat("  table_s7_race6_sensitivity.csv\n")
cat("  table_s8_protein_sensitivity.csv\n")
cat("  table_s9_lod_40_70.csv\n")
cat("  table_s10_power_analysis.csv\n")
cat("  table_s11_strobe_checklist.csv\n")

cat("\n=== 10_reviewer_analyses.R complete ===\n")
