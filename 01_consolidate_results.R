# ============================================================================
# 01_consolidate_results.R
# Merge all three ExWAS screening rounds into a single master results table
# ============================================================================

source("00_functions.R")

# --- 1. Load all result files ------------------------------------------------

# Round 1: PFAS-thyroid (60 tests, 0 FDR-sig)
r1 <- read.csv("exwas_results_pfas_thyroid.csv", stringsAsFactors = FALSE) %>%
  mutate(screen = "PFAS-Thyroid", .before = 1)

# Round 2: Broad screen (648 tests, 0 FDR-sig)
r2 <- read.csv("exwas_broad_screen_all_results.csv", stringsAsFactors = FALSE) %>%
  mutate(screen = "Broad", .before = 1)

# Round 3: Expanded screen (1920 tests, 27 FDR-sig)
r3 <- read.csv("exwas_expanded_all_results.csv", stringsAsFactors = FALSE) %>%
  mutate(screen = "Expanded", .before = 1)

cat(sprintf("Round 1 (PFAS-Thyroid): %d tests\n", nrow(r1)))
cat(sprintf("Round 2 (Broad):        %d tests\n", nrow(r2)))
cat(sprintf("Round 3 (Expanded):     %d tests\n", nrow(r3)))

# Round 4: Novelty screen (from 01b_recover_novelty_screen.R)
has_novelty <- file.exists("exwas_novelty_screen_results.csv")
if (has_novelty) {
  r4 <- read.csv("exwas_novelty_screen_results.csv", stringsAsFactors = FALSE) %>%
    mutate(screen = "Novelty", .before = 1)
  cat(sprintf("Round 4 (Novelty):      %d tests\n", nrow(r4)))
} else {
  cat("Round 4 (Novelty):      NOT FOUND - run 01b_recover_novelty_screen.R first\n")
}

# --- 2. Standardise columns -------------------------------------------------

# Round 1 columns: exposure, outcome, outcome_type, beta, se, t_value, p_value,
#   n, p_fdr, p_bonferroni, sig_fdr05, sig_fdr01, sig_bonf,
#   exposure_label, outcome_label
r1_std <- r1 %>%
  transmute(
    screen,
    exposure,
    outcome,
    beta, se, t_value, p_value, n,
    df       = NA_integer_,
    weight   = "WTSSBJ2Y",   # PFAS were surplus samples
    exp_label  = exposure_label,
    out_label  = outcome_label,
    chem_class = "PFAS",
    out_domain = "Thyroid",
    type       = outcome_type
  )

# Round 2 columns: exposure, outcome, beta, se, t_value, p_value, n,
#   exp_label, chem_class, is_urinary, out_label, out_domain, type,
#   p_fdr, p_bonf, sig_nom, sig_fdr05, sig_fdr01, sig_bonf05, direction
r2_std <- r2 %>%
  transmute(
    screen,
    exposure, outcome,
    beta, se, t_value, p_value, n,
    df       = NA_integer_,
    weight   = ifelse(is_urinary, "WTSA2YR", "WTMEC2YR"),
    exp_label, out_label, chem_class, out_domain, type
  )

# Round 3 columns: exposure, outcome, beta, se, t_value, p_value, n, df,
#   p_fdr, sig_fdr05, sig_fdr01, sig_bonf
r3_std <- r3 %>%
  transmute(
    screen,
    exposure, outcome,
    beta, se, t_value, p_value, n, df,
    weight   = NA_character_,
    exp_label  = CHEM_LABELS[exposure],
    out_label  = OUTCOME_LABELS[outcome],
    chem_class = case_when(
      grepl("^log_LBXPF|^log_LBXMPAH", exposure)  ~ "PFAS",
      grepl("^log_LBX(BPB|BCD|THG|BSE|BMN|BGM|BCO)", exposure) ~ "Heavy metals",
      grepl("^log_URX(M|EC)", exposure) ~ "Phthalates",
      grepl("^log_URXP", exposure)  ~ "PAHs",
      grepl("^log_URXU", exposure)  ~ "Urinary metals/elements",
      grepl("^log_URX(HPM|OXY|AAM|AMC|BMA|CYH|DHB|CEM|ATC|BPM)", exposure) ~ "VOC metabolites",
      grepl("^log_LBXV", exposure)  ~ "Blood VOCs",
      grepl("^log_SS", exposure)    ~ "Surplus serum",
      TRUE                          ~ "Other"
    ),
    out_domain = case_when(
      outcome %in% c("BMXBMI", "BMXWAIST")   ~ "Anthropometric",
      outcome %in% c("eGFR", "LBXSBU")       ~ "Kidney",
      outcome %in% c("LBXSATSI", "LBXSASSI", "LBXSGB", "LBXSAPSI",
                      "log_LBXSAPSI", "LBXSTB") ~ "Liver",
      outcome %in% c("LBXGH", "LBXGLU", "LBXTC", "LBDHDD", "LBXSTR",
                      "LBXSUA", "log_HOMA_IR", "log_LBXIN") ~ "Metabolic",
      outcome %in% c("LBXWBCSI", "LBXRBCSI", "LBXHGB", "LBXPLTSI") ~ "Hematology",
      outcome %in% c("mean_sys", "mean_dia", "avg_SBP", "avg_DBP") ~ "Cardiovascular",
      outcome %in% c("PHQ9_total", "depression_binary") ~ "Mental health",
      outcome %in% c("LBXHSCRP", "log_LBXHSCRP") ~ "Inflammation",
      outcome == "SLD012" ~ "Sleep",
      TRUE ~ "Other"
    ),
    type = ifelse(outcome == "depression_binary", "binary", "continuous")
  )

# Round 4 standardisation (novelty screen)
if (has_novelty) {
  r4_std <- r4 %>%
    transmute(
      screen,
      exposure, outcome,
      beta, se, t_value, p_value, n, df,
      weight,
      exp_label  = CHEM_LABELS[exposure],
      out_label  = OUTCOME_LABELS[outcome],
      chem_class = case_when(
        grepl("^log_SS", exposure)    ~ "Surplus serum",
        grepl("^log_LBXV", exposure)  ~ "Blood VOCs",
        grepl("^log_LBXBGM", exposure) ~ "Heavy metals",
        grepl("^log_LBXBCO", exposure) ~ "Heavy metals",
        grepl("^log_URXU(IO|P8)", exposure) ~ "Urinary elements",
        grepl("^log_URXU(AS|DMA)", exposure) ~ "Urinary metals/elements",
        TRUE ~ "Other"
      ),
      out_domain = case_when(
        outcome %in% c("BMXBMI", "BMXWAIST")   ~ "Anthropometric",
        outcome %in% c("eGFR", "LBXSBU")       ~ "Kidney",
        outcome %in% c("LBXSATSI", "LBXSASSI", "LBXSGB", "LBXSAPSI",
                        "log_LBXSAPSI", "LBXSTB") ~ "Liver",
        outcome %in% c("LBXGH", "LBXGLU", "LBXTC", "LBDHDD", "LBXSTR",
                        "LBXSUA", "log_HOMA_IR", "log_LBXIN") ~ "Metabolic",
        outcome %in% c("LBXWBCSI", "LBXRBCSI", "LBXHGB", "LBXPLTSI") ~ "Hematology",
        outcome %in% c("mean_sys", "mean_dia") ~ "Cardiovascular",
        outcome %in% c("PHQ9_total", "depression_binary") ~ "Mental health",
        outcome %in% c("LBXHSCRP", "log_LBXHSCRP") ~ "Inflammation",
        outcome == "SLD012" ~ "Sleep",
        TRUE ~ "Other"
      ),
      type = ifelse(outcome == "depression_binary", "binary", "continuous")
    )
}

# --- 3. Combine all rounds ---------------------------------------------------

if (has_novelty) {
  master <- bind_rows(r1_std, r2_std, r3_std, r4_std)
} else {
  master <- bind_rows(r1_std, r2_std, r3_std)
}

cat(sprintf("\nMaster table: %d total tests\n", nrow(master)))

# --- 4. Apply global FDR correction across all tests -------------------------

master <- master %>%
  mutate(
    p_fdr_global    = p.adjust(p_value, method = "BH"),
    p_bonf_global   = p.adjust(p_value, method = "bonferroni"),
    sig_nominal     = p_value < 0.05,
    sig_fdr05       = p_fdr_global < 0.05,
    sig_fdr01       = p_fdr_global < 0.01,
    sig_bonf        = p_bonf_global < 0.05,
    direction       = ifelse(beta > 0, "Positive", "Negative")
  ) %>%
  arrange(p_value)

# Fill in missing labels
master <- master %>%
  mutate(
    exp_label = coalesce(exp_label, CHEM_LABELS[exposure], exposure),
    out_label = coalesce(out_label, OUTCOME_LABELS[outcome], outcome)
  )

# --- 5. Summary --------------------------------------------------------------

cat(sprintf("\n=== MASTER RESULTS SUMMARY ===\n"))
cat(sprintf("Total tests:           %d\n", nrow(master)))
cat(sprintf("Nominal p < 0.05:      %d (%.1f%%)\n",
            sum(master$sig_nominal), mean(master$sig_nominal) * 100))
cat(sprintf("FDR < 0.05 (global):   %d\n", sum(master$sig_fdr05)))
cat(sprintf("FDR < 0.01 (global):   %d\n", sum(master$sig_fdr01)))
cat(sprintf("Bonferroni < 0.05:     %d\n", sum(master$sig_bonf)))

cat("\n--- FDR-significant findings ---\n")
master %>%
  filter(sig_fdr05) %>%
  select(screen, exp_label, out_label, beta, p_value, p_fdr_global, n, direction) %>%
  as.data.frame() %>%
  print()

cat("\n--- By chemical class ---\n")
master %>%
  filter(sig_fdr05) %>%
  count(chem_class, sort = TRUE) %>%
  as.data.frame() %>%
  print()

cat("\n--- By outcome domain ---\n")
master %>%
  filter(sig_fdr05) %>%
  count(out_domain, sort = TRUE) %>%
  as.data.frame() %>%
  print()

# --- 6. Create priority table for novelty assessment -------------------------

priority_findings <- master %>%
  filter(sig_fdr05) %>%
  mutate(
    novelty_notes = case_when(
      # HIGH novelty
      grepl("SSGLYP", exposure) ~
        "HIGH: Glyphosate NHANES surplus data is new; very few publications",
      grepl("LBXVBZ", exposure) & grepl("depression|PHQ", outcome, ignore.case = TRUE) ~
        "HIGH: Blood benzene-depression link largely unstudied",
      grepl("LBXVBZ", exposure) & grepl("WBC", outcome) ~
        "LOW: Benzene hematotoxicity is textbook",
      grepl("LBXVBZ", exposure) ~
        "MODERATE: Blood-level benzene (vs urinary metabolite) less studied in ExWAS",

      # MODERATE novelty
      grepl("URXUIO", exposure) & grepl("BMI|WAIST", outcome) ~
        "MODERATE: Iodine-obesity emerging but limited ExWAS",
      grepl("URXUIO", exposure) ~
        "MODERATE: Iodine-health links understudied at population level",
      grepl("URXUP8", exposure) ~
        "MODERATE: Perchlorate known thyroid disruptor; kidney/other links less studied",
      grepl("LBXBCO", exposure) ~
        "MODERATE: Cobalt population-level associations less studied",
      grepl("LBXBGM", exposure) ~
        "MODERATE: Methylmercury (vs total mercury) specific associations less common",
      grepl("URXM", exposure) & grepl("LBXSTB", outcome) ~
        "MODERATE: Phthalate-bilirubin link less common in literature",
      grepl("URXOXY", exposure) ~
        "MODERATE: Oxychlordane-health associations sparsely studied in NHANES",

      # LOW-MODERATE novelty
      grepl("LBXBMN", exposure) & grepl("RBC|eGFR|HGB", outcome) ~
        "LOW-MODERATE: Manganese-hematology partially known",
      grepl("LBXBMN", exposure) & grepl("BMI|WAIST", outcome) ~
        "LOW-MODERATE: Manganese-adiposity emerging evidence",
      grepl("URXU(CD|CO|CS|PB|TL)", exposure) & grepl("eGFR", outcome) ~
        "LOW-MODERATE: Urinary metals-kidney partially established",
      grepl("URXUCD", exposure) & grepl("BUN|SBU", outcome) ~
        "LOW-MODERATE: Urinary cadmium-kidney function partially established",
      grepl("URXU(AS|DMA)", exposure) ~
        "LOW-MODERATE: Arsenic species-health links partially documented",
      grepl("LBXTHG", exposure) & grepl("APSI|ALP", outcome, ignore.case = TRUE) ~
        "LOW-MODERATE: Mercury-liver enzyme associations emerging",

      # LOW novelty (well-established)
      grepl("LBXBPB", exposure) ~
        "LOW: Lead-health associations well established",
      grepl("LBXBSE", exposure) ~
        "LOW: Selenium-health associations well established",
      grepl("LBXBCD", exposure) ~
        "LOW: Cadmium-health associations well established",

      TRUE ~ "NEEDS ASSESSMENT"
    )
  ) %>%
  arrange(novelty_notes, p_value) %>%
  select(exposure, exp_label, outcome, out_label, beta, se, p_value, p_fdr_global,
         n, direction, chem_class, out_domain, novelty_notes)

cat("\n=== PRIORITY FINDINGS FOR LITERATURE CHECK ===\n")
print(as.data.frame(priority_findings))

# --- 7. Save -----------------------------------------------------------------

write.csv(master, "exwas_master_results.csv", row.names = FALSE)
write.csv(priority_findings, "exwas_priority_findings.csv", row.names = FALSE)
save(master, priority_findings, file = "exwas_master_results.RData")

cat("\nSaved: exwas_master_results.csv, exwas_priority_findings.csv, exwas_master_results.RData\n")
