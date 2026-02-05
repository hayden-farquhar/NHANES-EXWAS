# ============================================================================
# 07_novelty_assessment.R
# Literature-based novelty assessment for validated ExWAS findings
# Based on systematic PubMed searches conducted 2026-02-05
# ============================================================================

source("00_functions.R")

# --- 1. Load validated findings -----------------------------------------------

load("validation_results.RData")  # validated

findings <- validated %>%
  filter(validated) %>%
  select(exposure, outcome, exp_label, out_label, beta_primary, beta_val,
         p_val, n_primary, n_val, chem_class, out_domain, dir_primary) %>%
  arrange(chem_class, exposure, outcome)

cat(sprintf("Validated findings to assess: %d\n\n", nrow(findings)))

# --- 2. Novelty assessment based on PubMed literature search ------------------
#
# PubMed searches conducted 2026-02-05 for each exposure-outcome pair.
# Novelty levels:
#   HIGH       = 0-1 directly relevant publications
#   MODERATE   = 2-5 publications, or only studied in different populations
#   LOW-MOD    = 5-15 publications, partially established
#   LOW        = >15 publications, well-established in NHANES/epi literature
#
# Each entry includes:
#   - PubMed hit count (approximate, from targeted searches)
#   - Key references (PMIDs)
#   - Direction consistency with our findings
#   - Assessment rationale

novelty <- tribble(
  ~exposure, ~outcome, ~novelty_level, ~pubmed_hits, ~key_pmids, ~lit_direction, ~rationale,

  # === BLOOD LEAD (4 findings) ===
  "log_LBXBPB", "LBXTC", "LOW", 53,
  "24553633,39511585,27555608",
  "Consistent",
  "Lead-cholesterol association well-established in NHANES literature. Multiple studies show positive lead-TC association, likely via oxidative stress and lipid metabolism disruption.",

  "log_LBXBPB", "BMXBMI", "LOW", 29,
  "37878298,33453048,36455248",
  "Consistent",
  "Lead-BMI inverse association documented in multiple NHANES analyses and meta-analyses. Reverse causation (higher BMI dilutes blood lead) is a known concern.",

  "log_LBXBPB", "BMXWAIST", "LOW", 29,
  "37878298,33453048,36455248",
  "Consistent",
  "Lead-waist circumference follows same pattern as lead-BMI. Same literature base. Well-established.",

  "log_LBXBPB", "LBXGH", "LOW", 262,
  "35420698,30198134,35002965",
  "Mixed",
  "Lead-diabetes/HbA1c extensively studied. Our inverse association (lead linked to lower HbA1c) contrasts with some studies showing positive association at higher exposures. May reflect healthy-worker/survivor bias or non-linear relationship at low-level exposure.",

  # === BLOOD SELENIUM (3 findings) ===
  "log_LBXBSE", "LBXHGB", "LOW-MODERATE", 4,
  "27885969,40965713",
  "Consistent",
  "Selenium-hemoglobin association has limited but consistent literature. Selenium is essential for erythropoiesis. Few large population studies specifically examine this pair.",

  "log_LBXBSE", "LBXRBCSI", "LOW-MODERATE", 4,
  "27885969,40965713",
  "Consistent",
  "Same literature base as selenium-hemoglobin. Selenium's role in erythropoiesis makes this biologically plausible but sparsely studied at population level.",

  "log_LBXBSE", "LBXTC", "LOW", 22,
  "20102763,2537285,39267744,18689378",
  "Consistent",
  "Selenium-cholesterol positive association well-documented in NHANES and other population studies. High selenium linked to dyslipidemia. Established finding.",

  # === MERCURY / METHYLMERCURY (3 findings) ===
  "log_LBXBGM", "LBXSAPSI", "MODERATE", 1,
  "36649900",
  "Partially consistent",
  "Li et al. 2023 (NHANES 2011-2018) studied blood metals and liver function but focused on Hg-bilirubin, not methylmercury-ALP specifically. Our finding of methylmercury inversely associated with ALP is novel in the ExWAS/NHANES context. Biological plausibility via hepatotoxicity.",

  "log_LBXTHG", "log_LBXSAPSI", "MODERATE", 1,
  "36649900",
  "Partially consistent",
  "Same single study (Li et al. 2023). Total mercury-ALP association largely unstudied at population level. Our cross-cycle validated finding adds new evidence.",

  "log_LBXBGM", "BMXWAIST", "HIGH", 0,
  "",
  "Novel",
  "No publications found for methylmercury-adiposity association. Inverse association (higher MeHg linked to lower waist) likely confounded by fish consumption (fish = MeHg source but also associated with lower obesity). Important to discuss confounding.",

  # === BLOOD MANGANESE (2 findings) ===
  "log_LBXBMN", "BMXBMI", "MODERATE", 2,
  "34993913,39887166",
  "Mixed",
  "Very limited literature. Only 2 studies found linking manganese to obesity/metabolic syndrome, both recent. Our positive Mn-BMI association is partially novel. Prenatal studies exist (Smith et al. 2022) but adult population data scarce.",

  "log_LBXBMN", "BMXWAIST", "MODERATE", 2,
  "34993913,39887166",
  "Mixed",
  "Same sparse literature as Mn-BMI. Positive Mn-waist association adds to emerging but limited evidence base.",

  # === URINARY IODINE (1 finding) ===
  "log_URXUIO", "BMXBMI", "MODERATE", 8,
  "33256547,21813915,32295049",
  "Consistent",
  "Small literature on iodine-obesity in population studies. De Angelis et al. 2021 found lower urinary iodine in obese children. Our positive association (higher iodine linked to higher BMI) in adults may reflect reverse causation (obese individuals have higher caloric/iodine intake) or thyroid-mediated pathway. Limited NHANES-specific evidence in adults.",

  # === URINARY DMA ARSENIC (1 finding) ===
  "log_URXUDMA", "LBXSUA", "HIGH", 0,
  "",
  "Novel",
  "No publications found linking dimethylarsonic acid (DMA) specifically to uric acid. Arsenic exposure has been linked to metabolic disruption broadly, but the DMA-uric acid association appears unstudied. Biologically plausible via purine metabolism disruption or renal tubular effects.",

  # === URINARY PERCHLORATE (1 finding) ===
  "log_URXUP8", "LBXSBU", "HIGH", 2,
  "37154820,40441702",
  "Partially consistent",
  "Only 2 publications on perchlorate-kidney in NHANES. Li et al. 2023 found perchlorate positively associated with eGFR (protective); Xue et al. 2025 studied adolescents. Our finding of perchlorate positively associated with BUN in adults is novel and not directly addressed in existing literature. Perchlorate primarily studied as thyroid disruptor, not nephrotoxin."
)

# --- 3. Merge with finding details -------------------------------------------

novelty_full <- findings %>%
  left_join(novelty, by = c("exposure", "outcome")) %>%
  mutate(
    novelty_level = factor(novelty_level,
                           levels = c("HIGH", "MODERATE", "LOW-MODERATE", "LOW")),
    publishability = case_when(
      novelty_level == "HIGH"         ~ "Lead story candidate",
      novelty_level == "MODERATE"     ~ "Strong supporting finding",
      novelty_level == "LOW-MODERATE" ~ "Confirmatory / context",
      novelty_level == "LOW"          ~ "Replication / background",
      TRUE                            ~ "Unassessed"
    )
  ) %>%
  arrange(novelty_level, exposure)

# --- 4. Display results -------------------------------------------------------

cat("=== NOVELTY ASSESSMENT: VALIDATED ExWAS FINDINGS ===\n\n")

for (lvl in c("HIGH", "MODERATE", "LOW-MODERATE", "LOW")) {
  sub <- novelty_full %>% filter(novelty_level == lvl)
  if (nrow(sub) == 0) next

  cat(sprintf("--- %s NOVELTY (%d findings) ---\n\n", lvl, nrow(sub)))

  for (i in 1:nrow(sub)) {
    cat(sprintf("  %s -> %s\n", sub$exp_label[i], sub$out_label[i]))
    cat(sprintf("    Beta (2017-18): %.3f | Beta (2015-16): %.3f | p_val: %.4f\n",
                sub$beta_primary[i], sub$beta_val[i], sub$p_val[i]))
    cat(sprintf("    Direction: %s | PubMed hits: %d\n",
                sub$dir_primary[i], sub$pubmed_hits[i]))
    cat(sprintf("    Literature: %s\n", sub$lit_direction[i]))
    cat(sprintf("    Rationale: %s\n\n", sub$rationale[i]))
  }
}

# --- 5. Publication strategy summary ------------------------------------------

cat("\n=== PUBLICATION STRATEGY ===\n\n")

cat("LEAD FINDINGS (highest novelty, most publishable):\n")
novelty_full %>%
  filter(novelty_level == "HIGH") %>%
  select(exp_label, out_label, dir_primary, pubmed_hits, publishability) %>%
  as.data.frame() %>%
  print()

cat("\nSTRONG SUPPORTING (moderate novelty, robust):\n")
novelty_full %>%
  filter(novelty_level == "MODERATE") %>%
  select(exp_label, out_label, dir_primary, pubmed_hits, publishability) %>%
  as.data.frame() %>%
  print()

cat("\nCONFIRMATORY (established associations replicated):\n")
novelty_full %>%
  filter(novelty_level %in% c("LOW-MODERATE", "LOW")) %>%
  select(exp_label, out_label, dir_primary, pubmed_hits, publishability) %>%
  as.data.frame() %>%
  print()

cat("\n=== RECOMMENDED MANUSCRIPT FRAMING ===\n\n")
cat("1. Frame as comprehensive ExWAS with cross-cycle validation (methodological strength)\n")
cat("2. Lead with the 3 HIGH-novelty findings:\n")
cat("   a) DMA arsenic -> Uric acid (entirely novel, 0 prior publications)\n")
cat("   b) Perchlorate -> BUN (novel renal association, only 2 tangential papers)\n")
cat("   c) Methylmercury -> Waist circumference (novel but discuss fish confounding)\n")
cat("3. Highlight MODERATE-novelty findings as supporting evidence:\n")
cat("   - Mercury/methylmercury-ALP (hepatotoxicity, minimal prior NHANES evidence)\n")
cat("   - Manganese-adiposity (emerging, very limited adult data)\n")
cat("   - Iodine-BMI (limited adult NHANES evidence)\n")
cat("4. Include LOW-novelty findings as validation of ExWAS methodology\n")
cat("   - Successfully replicated well-known lead, selenium associations\n")
cat("   - Demonstrates analytical pipeline is sound\n")
cat("5. Key limitations to address:\n")
cat("   - Cross-sectional design (no causal inference)\n")
cat("   - MeHg-waist: fish consumption confounding\n")
cat("   - Lead-BMI: reverse causation possibility\n")
cat("   - Iodine-BMI: dietary confounding\n")
cat("   - Low degrees of freedom for subsample analyses\n")

# --- 6. Save -----------------------------------------------------------------

write.csv(novelty_full, "novelty_assessment.csv", row.names = FALSE)
save(novelty_full, file = "novelty_assessment.RData")

cat("\nSaved: novelty_assessment.csv, novelty_assessment.RData\n")
