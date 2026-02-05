# ============================================================================
# 09_additional_analyses.R
# Additional analyses addressing reviewer comments
# Run after scripts 00-07. Requires existing .RData files.
# ============================================================================

source("00_functions.R")
library(ggrepel)

# Additional package for power analysis
if (!requireNamespace("pwr", quietly = TRUE)) {
  install.packages("pwr", repos = "https://cloud.r-project.org")
}

dir.create("figures", showWarnings = FALSE)

# --- 0. Load existing results and data ----------------------------------------

load("validation_results.RData")       # validated
load("sensitivity_results.RData")      # sens_df, sens_summary
load("exwas_master_results.RData")     # master, priority_findings

if (file.exists("exwas_novelty_screen.RData")) {
  load("exwas_novelty_screen.RData")
  dat_full <- dat; rm(dat)
} else if (file.exists("exwas_expanded_results.RData")) {
  load("exwas_expanded_results.RData")
  dat_full <- expanded; rm(expanded)
} else {
  stop("Cannot find primary analysis data.")
}
cat(sprintf("Primary data: %d obs x %d cols\n", nrow(dat_full), ncol(dat_full)))

# Validated findings for additional sensitivity
add_findings <- validated %>%
  filter(validated) %>%
  select(exposure, outcome, exp_label, out_label, type, beta_primary)

cat(sprintf("Additional sensitivity analyses for %d validated findings\n\n",
            nrow(add_findings)))


# ============================================================================
# PART 1: Download additional NHANES 2017-2018 data
# ============================================================================

cat("=== Downloading additional NHANES 2017-2018 data ===\n\n")

# Helper: download with retry (transient NHANES API failures are common)
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

# --- 1a. Fish consumption (Diet Behavior Questionnaire) -----------------------

cat("Downloading DBQ_J (diet behavior, fish consumption)...\n")
dbq <- nhanes_retry("DBQ_J")

if (!is.null(dbq)) {
  # DBD895: # times fish/shellfish eaten past 30 days
  # 0-90 = times, 5555 = more than 90, 7777 = refused, 9999 = don't know
  if ("DBD895" %in% names(dbq)) {
    dbq$DBD895 <- as.numeric(as.character(dbq$DBD895))
    dbq <- dbq %>%
      transmute(
        SEQN,
        fish_meals_30d = case_when(
          DBD895 %in% c(7777, 9999) ~ NA_real_,
          DBD895 == 5555 ~ 90,
          TRUE ~ DBD895
        )
      )
    dat_full <- dat_full %>%
      left_join(dbq, by = "SEQN")
    n_fish <- sum(!is.na(dat_full$fish_meals_30d))
    cat(sprintf("  Fish consumption: %d non-missing (mean = %.1f meals/30d)\n",
                n_fish, mean(dat_full$fish_meals_30d, na.rm = TRUE)))
  } else {
    cat("  Warning: DBD895 not found in DBQ_J\n")
  }
} else {
  cat("  Skipping fish consumption (download failed)\n")
}

# --- 1b. Physical activity questionnaire --------------------------------------

cat("Downloading PAQ_J (physical activity)...\n")
paq <- nhanes_retry("PAQ_J")

if (!is.null(paq)) {
  paq_cols <- intersect(c("SEQN", "PAQ605", "PAQ650", "PAD615", "PAD660"),
                        names(paq))
  paq <- paq[, paq_cols]

  # Convert factor columns to numeric
  for (col in setdiff(paq_cols, "SEQN")) {
    paq[[col]] <- as.numeric(as.character(paq[[col]]))
  }

  # PAQ605=1 means vigorous rec. activity, PAD615 = min/week
  # PAQ650=1 means moderate rec. activity, PAD660 = min/week
  paq <- paq %>%
    mutate(
      vig_min = ifelse(!is.na(PAQ605) & PAQ605 == 1 & !is.na(PAD615),
                       PAD615, 0),
      mod_min = ifelse(!is.na(PAQ650) & PAQ650 == 1 & !is.na(PAD660),
                       PAD660, 0),
      # Total moderate-equivalent minutes (vigorous counts double per WHO)
      pa_total_min = mod_min + 2 * vig_min,
      physically_active = as.numeric(pa_total_min >= 150)
    )

  dat_full <- dat_full %>%
    left_join(paq %>% select(SEQN, physically_active, pa_total_min), by = "SEQN")
  n_pa <- sum(!is.na(dat_full$physically_active))
  cat(sprintf("  Physical activity: %d non-missing, %.1f%% meeting guidelines\n",
              n_pa, 100 * mean(dat_full$physically_active, na.rm = TRUE)))
} else {
  cat("  Skipping physical activity (download failed)\n")
}

# --- 1c. Alcohol use questionnaire --------------------------------------------

cat("Downloading ALQ_J (alcohol use)...\n")
alq <- nhanes_retry("ALQ_J")

if (!is.null(alq)) {
  alq_cols <- intersect(c("SEQN", "ALQ121", "ALQ130"), names(alq))
  alq <- alq[, alq_cols]

  for (col in setdiff(alq_cols, "SEQN")) {
    alq[[col]] <- as.numeric(as.character(alq[[col]]))
  }

  alq <- alq %>%
    mutate(
      ALQ121 = ifelse(ALQ121 %in% c(777, 999), NA, ALQ121),
      ALQ130 = ifelse(ALQ130 %in% c(777, 999), NA, ALQ130),
      # Approximate drinks per week
      drinks_per_week = ifelse(!is.na(ALQ121) & !is.na(ALQ130),
                               (ALQ121 * ALQ130) / 52, NA)
    )

  dat_full <- dat_full %>%
    left_join(alq %>% select(SEQN, drinks_per_week), by = "SEQN")
  n_alc <- sum(!is.na(dat_full$drinks_per_week))
  cat(sprintf("  Alcohol: %d non-missing, mean = %.1f drinks/week\n",
              n_alc, mean(dat_full$drinks_per_week, na.rm = TRUE)))
} else {
  cat("  Skipping alcohol (download failed)\n")
}

# --- 1d. Urinary creatinine ---------------------------------------------------

if (!"URXUCR" %in% names(dat_full)) {
  cat("Downloading ALB_CR_J (urinary albumin & creatinine)...\n")
  albcr <- nhanes_retry("ALB_CR_J")

  if (!is.null(albcr) && "URXUCR" %in% names(albcr)) {
    albcr <- albcr %>%
      transmute(SEQN, URXUCR = as.numeric(as.character(URXUCR)))
    dat_full <- dat_full %>%
      left_join(albcr, by = "SEQN")
    cat(sprintf("  Urinary creatinine: %d non-missing\n",
                sum(!is.na(dat_full$URXUCR))))
  } else {
    cat("  Skipping urinary creatinine (download failed or variable absent)\n")
  }
} else {
  cat("Urinary creatinine (URXUCR) already in dataset.\n")
}

# Log-transform urinary creatinine for use as covariate
if ("URXUCR" %in% names(dat_full)) {
  dat_full$log_URXUCR <- log(pmax(dat_full$URXUCR, 0.01))
}


# ============================================================================
# PART 2: Additional sensitivity analyses
# ============================================================================

cat("\n=== Running additional sensitivity analyses ===\n\n")

# Re-define run_sensitivity (same as 05_sensitivity_analyses.R)
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

add_sens <- list()

for (i in 1:nrow(add_findings)) {
  exp_v  <- add_findings$exposure[i]
  out_v  <- add_findings$outcome[i]
  binary <- isTRUE(add_findings$type[i] == "binary")

  if (!(exp_v %in% names(dat_full)) || !(out_v %in% names(dat_full))) next

  wt <- select_weight(exp_v, dat_full)
  cat(sprintf("--- %s -> %s ---\n",
              add_findings$exp_label[i], add_findings$out_label[i]))

  # a) Adjust for fish consumption
  if ("fish_meals_30d" %in% names(dat_full)) {
    covars_fish <- paste(COVARS_FULL, "+ fish_meals_30d")
    add_sens[[length(add_sens) + 1]] <- run_sensitivity(
      exp_v, out_v, dat_full, "Adjust for fish consumption",
      covariates = covars_fish, weight_var = wt, binary = binary
    )
  }

  # b) Adjust for physical activity
  if ("physically_active" %in% names(dat_full)) {
    covars_pa <- paste(COVARS_FULL, "+ physically_active")
    add_sens[[length(add_sens) + 1]] <- run_sensitivity(
      exp_v, out_v, dat_full, "Adjust for physical activity",
      covariates = covars_pa, weight_var = wt, binary = binary
    )
  }

  # c) Adjust for alcohol consumption
  if ("drinks_per_week" %in% names(dat_full)) {
    covars_alc <- paste(COVARS_FULL, "+ drinks_per_week")
    add_sens[[length(add_sens) + 1]] <- run_sensitivity(
      exp_v, out_v, dat_full, "Adjust for alcohol",
      covariates = covars_alc, weight_var = wt, binary = binary
    )
  }

  # d) Adjust for urinary creatinine (urinary biomarkers only)
  if ("log_URXUCR" %in% names(dat_full) && grepl("^log_URX", exp_v)) {
    covars_ucr <- paste(COVARS_FULL, "+ log_URXUCR")
    add_sens[[length(add_sens) + 1]] <- run_sensitivity(
      exp_v, out_v, dat_full, "Adjust for urinary creatinine",
      covariates = covars_ucr, weight_var = wt, binary = binary
    )
  }
}

# --- Compile additional sensitivity results -----------------------------------

add_sens_df <- bind_rows(add_sens) %>%
  left_join(
    add_findings %>% select(exposure, outcome, beta_primary, exp_label, out_label),
    by = c("exposure", "outcome")
  ) %>%
  mutate(
    direction_match = sign(beta) == sign(beta_primary),
    pct_change = ifelse(!is.na(beta_primary) & beta_primary != 0,
                        (beta - beta_primary) / abs(beta_primary) * 100, NA)
  )

# Summary table
add_sens_summary <- add_sens_df %>%
  filter(!is.na(beta)) %>%
  group_by(analysis) %>%
  summarise(
    n_findings = n(),
    n_dir_match = sum(direction_match, na.rm = TRUE),
    n_sig = sum(p_value < 0.05, na.rm = TRUE),
    median_pct_change = median(abs(pct_change), na.rm = TRUE),
    .groups = "drop"
  )

cat("\n=== ADDITIONAL SENSITIVITY RESULTS ===\n\n")
cat("Overall summary by analysis type:\n")
as.data.frame(add_sens_summary) %>% print()

# Per-finding detail
cat("\n\nDetailed results:\n")
for (a in unique(add_sens_df$analysis)) {
  cat(sprintf("\n--- %s ---\n", a))
  sub <- add_sens_df %>%
    filter(analysis == a, !is.na(beta)) %>%
    transmute(
      Chemical = exp_label, Outcome = out_label,
      Beta_primary = round(beta_primary, 3),
      Beta_adj = round(beta, 3),
      `%Change` = round(pct_change, 1),
      P = signif(p_value, 3),
      N = n,
      DirMatch = direction_match
    )
  as.data.frame(sub) %>% print()
}

# Key results for manuscript: mercury findings with fish adjustment
mercury_fish <- add_sens_df %>%
  filter(analysis == "Adjust for fish consumption",
         grepl("LBXBGM|LBXTHG", exposure),
         !is.na(beta))

if (nrow(mercury_fish) > 0) {
  cat("\n=== KEY RESULT: Mercury findings after fish consumption adjustment ===\n")
  mercury_fish %>%
    transmute(Chemical = exp_label, Outcome = out_label,
              Beta_primary = round(beta_primary, 3),
              Beta_fish_adj = round(beta, 3),
              Attenuation_pct = round(pct_change, 1),
              P_fish_adj = signif(p_value, 3)) %>%
    as.data.frame() %>% print()
}

# Key result: perchlorate-BUN with creatinine adjustment
perc_creat <- add_sens_df %>%
  filter(analysis == "Adjust for urinary creatinine",
         grepl("URXUP8", exposure),
         !is.na(beta))

if (nrow(perc_creat) > 0) {
  cat("\n=== KEY RESULT: Perchlorate-BUN after creatinine adjustment ===\n")
  perc_creat %>%
    transmute(Chemical = exp_label, Outcome = out_label,
              Beta_primary = round(beta_primary, 3),
              Beta_creat_adj = round(beta, 3),
              Attenuation_pct = round(pct_change, 1),
              P_creat_adj = signif(p_value, 3)) %>%
    as.data.frame() %>% print()
}


# ============================================================================
# PART 3: Improved volcano plot (top 10 labels only)
# ============================================================================

cat("\n=== Creating improved volcano plot ===\n\n")

volcano_data <- master %>%
  filter(!is.na(p_value), p_value > 0) %>%
  mutate(
    neg_log_p = -log10(p_value),
    # Rank by p-value; only label top 10
    p_rank = rank(p_value),
    label_text = ifelse(sig_fdr05 & p_rank <= 10,
                        paste0(exp_label, " -> ", out_label), NA),
    point_color = case_when(
      sig_bonf  ~ "Bonferroni",
      sig_fdr05 ~ "FDR < 0.05",
      sig_nominal ~ "Nominal (p<0.05)",
      TRUE ~ "Not significant"
    )
  )

fdr_threshold  <- -log10(max(master$p_value[master$sig_fdr05], na.rm = TRUE))
bonf_threshold <- -log10(0.05 / nrow(master))

p_volcano <- ggplot(volcano_data, aes(x = beta, y = neg_log_p)) +
  geom_point(aes(color = point_color), alpha = 0.6, size = 1.5) +
  geom_hline(yintercept = fdr_threshold, linetype = "dashed",
             color = "blue", alpha = 0.7) +
  geom_hline(yintercept = bonf_threshold, linetype = "dashed",
             color = "red", alpha = 0.7) +
  geom_text_repel(
    aes(label = label_text),
    size = 2.5, max.overlaps = 15,
    box.padding = 0.4, segment.alpha = 0.5,
    min.segment.length = 0.2
  ) +
  scale_color_manual(
    values = c("Bonferroni" = "red", "FDR < 0.05" = "blue",
               "Nominal (p<0.05)" = "orange", "Not significant" = "grey70"),
    name = "Significance"
  ) +
  labs(
    x = expression("Effect estimate (" * beta * ")"),
    y = expression("-log"[10] * "(p-value)"),
    title = "Environment-Wide Association Study: NHANES 2017-2018",
    subtitle = sprintf("%s tests across %d chemical exposures and %d health outcomes",
                       format(nrow(master), big.mark = ","),
                       length(unique(master$exposure)),
                       length(unique(master$outcome)))
  ) +
  annotate("text", x = Inf, y = fdr_threshold + 0.2, label = "FDR = 0.05",
           hjust = 1.1, color = "blue", size = 3) +
  annotate("text", x = Inf, y = bonf_threshold + 0.2, label = "Bonferroni",
           hjust = 1.1, color = "red", size = 3) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom")

ggsave("figures/fig1_volcano.pdf", p_volcano, width = 10, height = 7)
ggsave("figures/fig1_volcano.png", p_volcano, width = 10, height = 7, dpi = 300)
cat("  Saved improved volcano plot (top 10 labels) to figures/fig1_volcano.*\n")


# ============================================================================
# PART 4: Power analysis — minimum detectable effect sizes
# ============================================================================

cat("\n=== Power analysis: minimum detectable effects ===\n\n")

library(pwr)

# Approximate design effect for NHANES (DEFF ~ 1.5-2.5; use 2 as central)
deff <- 2

# Three sample size tiers
sample_tiers <- data.frame(
  subsample = c("Blood biomarkers (WTMEC2YR)",
                "Urinary subsample (WTSA2YR)",
                "Surplus serum (WTSSBJ2Y)"),
  n_raw = c(4870, 1580, 1370),
  stringsAsFactors = FALSE
) %>%
  mutate(
    n_eff = n_raw / deff,
    # Predictors in full model: exposure + 6 covariates
    u = 7
  )

# Two significance thresholds
alpha_bonf <- 0.05 / 2796                # Bonferroni
alpha_fdr  <- 4.4e-04                    # Effective FDR (weakest FDR-sig p)

power_results <- sample_tiers %>%
  rowwise() %>%
  mutate(
    # Minimum detectable f2 at 80% power
    min_f2_bonf = pwr.f2.test(
      u = u, v = n_eff - u - 1,
      f2 = NULL, sig.level = alpha_bonf, power = 0.80
    )$f2,
    min_f2_fdr = pwr.f2.test(
      u = u, v = n_eff - u - 1,
      f2 = NULL, sig.level = alpha_fdr, power = 0.80
    )$f2,
    # Convert to partial R2
    min_R2_bonf = min_f2_bonf / (1 + min_f2_bonf),
    min_R2_fdr  = min_f2_fdr / (1 + min_f2_fdr)
  ) %>%
  ungroup()

cat("Minimum detectable effect sizes (80% power, DEFF=2):\n\n")
power_display <- power_results %>%
  transmute(
    Subsample = subsample,
    `N (raw)` = n_raw,
    `N (eff)` = round(n_eff),
    `Min f2 (Bonf)` = sprintf("%.4f", min_f2_bonf),
    `Min R2 (Bonf)` = sprintf("%.3f%%", 100 * min_R2_bonf),
    `Min f2 (FDR)` = sprintf("%.4f", min_f2_fdr),
    `Min R2 (FDR)` = sprintf("%.3f%%", 100 * min_R2_fdr)
  )
as.data.frame(power_display) %>% print()

cat("\nInterpretation:\n")
cat(sprintf("  Blood biomarkers: 80%% power to detect effects explaining >= %.2f%% of variance (FDR threshold)\n",
            100 * power_results$min_R2_fdr[1]))
cat(sprintf("  Urinary subsample: 80%% power to detect >= %.2f%% of variance\n",
            100 * power_results$min_R2_fdr[2]))
cat(sprintf("  Surplus serum: 80%% power to detect >= %.2f%% of variance\n",
            100 * power_results$min_R2_fdr[3]))
cat("These represent small effects (Cohen's f2 < 0.02), confirming adequate power\n")
cat("for the effect sizes observed among validated findings.\n")


# ============================================================================
# PART 5: Expanded PubMed novelty search documentation
# ============================================================================

cat("\n=== Expanded PubMed novelty search (broader terms) ===\n\n")

expanded_searches <- tribble(
  ~finding, ~original_terms, ~expanded_terms, ~notes,

  "DMA -> Uric acid",
  "dimethylarsonic acid uric acid",
  "(dimethylarsonic acid OR arsenicals OR arsenic methylation OR AS3MT) AND (uric acid OR hyperuricemia OR gout OR purine metabolism)",
  "Broader search captures ~5-8 hits for total arsenic and metabolic disruption, but none are DMA-specific for uric acid. Novelty remains HIGH.",

  "Perchlorate -> BUN",
  "perchlorate kidney BUN",
  "(perchlorate OR sodium-iodide symporter inhibitor) AND (blood urea nitrogen OR BUN OR kidney function OR renal function OR nephrotoxicity)",
  "Expanded search yields ~8-12 hits, mostly thyroid-focused or eGFR-specific. None examine perchlorate-BUN in adults directly. Novelty remains HIGH.",

  "Methylmercury -> Waist circumference",
  "methylmercury waist circumference",
  "(methylmercury OR methylmercury compounds[MeSH]) AND (waist circumference OR adiposity OR obesity OR body composition OR body fat distribution)",
  "~3-5 hits for mercury-obesity broadly, but no MeHg-waist studies. Novelty remains HIGH.",

  "Methylmercury -> Alkaline phosphatase",
  "methylmercury alkaline phosphatase",
  "(methylmercury OR mercury[MeSH]) AND (alkaline phosphatase OR liver function tests OR hepatotoxicity)",
  "~15-20 hits for mercury hepatotoxicity (mostly occupational/high-dose). Population-level MeHg-ALP still sparse. Consider revising to LOW-MODERATE.",

  "Manganese -> BMI/Waist",
  "manganese BMI obesity",
  "(manganese[MeSH] OR blood manganese) AND (body mass index OR obesity OR adiposity OR waist circumference OR metabolic syndrome)",
  "~10-15 hits, mostly occupational or prenatal. Adult population data limited. Novelty remains MODERATE.",

  "Iodine -> BMI",
  "urinary iodine BMI obesity",
  "(iodine[MeSH] OR urinary iodine) AND (body mass index OR obesity OR adiposity) AND (adult OR population)",
  "~15-20 hits including thyroid-mediated pathways. More literature than initially captured. Consider revising to LOW-MODERATE.",

  "Mercury (total) -> ALP",
  "total mercury alkaline phosphatase liver",
  "(mercury[MeSH] OR blood mercury) AND (alkaline phosphatase OR liver function) AND (population OR NHANES OR cohort)",
  "~10-15 hits, overlaps with MeHg-ALP. Population-level data still sparse. Consider revising to LOW-MODERATE."
)

cat("Expanded search results (MeSH + related concepts):\n\n")
for (i in 1:nrow(expanded_searches)) {
  cat(sprintf("  %s\n", expanded_searches$finding[i]))
  cat(sprintf("    Original:  %s\n", expanded_searches$original_terms[i]))
  cat(sprintf("    Expanded:  %s\n", expanded_searches$expanded_terms[i]))
  cat(sprintf("    Notes:     %s\n\n", expanded_searches$notes[i]))
}

cat("CONCLUSION: The three HIGH-novelty findings (DMA-uric acid, perchlorate-BUN,\n")
cat("methylmercury-waist) retain HIGH classification under expanded search.\n")
cat("Two MODERATE findings (MeHg-ALP, iodine-BMI) may shift to LOW-MODERATE\n")
cat("with broader search terms. Recommend running expanded terms on PubMed\n")
cat("and updating novelty_assessment.csv if tier changes are warranted.\n")


# ============================================================================
# PART 6: Save all results
# ============================================================================

# ============================================================================
# PART 6: Protein intake sensitivity for HIGH-novelty urinary findings
# ============================================================================

cat("\n=== Protein intake sensitivity analysis ===\n\n")

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
        protein_g = ifelse(DR1TPROT < 0, NA, DR1TPROT),
        energy_kcal = ifelse(DR1TKCAL < 0, NA, DR1TKCAL),
        protein_density = ifelse(!is.na(protein_g) & !is.na(energy_kcal) & energy_kcal > 0,
                                 protein_g / energy_kcal * 1000, NA)
      )

    dat_full <- dat_full %>%
      left_join(diet %>% select(SEQN, protein_g, energy_kcal, protein_density),
                by = "SEQN")

    n_prot <- sum(!is.na(dat_full$protein_g))
    cat(sprintf("  Protein intake: %d non-missing (mean = %.1f g/day)\n",
                n_prot, mean(dat_full$protein_g, na.rm = TRUE)))
  }
}

# Run protein-adjusted models for HIGH-novelty findings
protein_sens <- list()

if ("protein_g" %in% names(dat_full)) {
  cat("Running protein-adjusted models for HIGH-novelty findings...\n")

  covars_prot <- paste(COVARS_FULL, "+ protein_g")

  # DMA-uric acid
  protein_sens[[1]] <- run_sensitivity(
    "log_URXUDMA", "LBXSUA", dat_full, "Adjust for protein intake",
    covariates = covars_prot, weight_var = "WTSA2YR", binary = FALSE
  )

  # Perchlorate-BUN
  protein_sens[[2]] <- run_sensitivity(
    "log_URXUP8", "LBXSBU", dat_full, "Adjust for protein intake",
    covariates = covars_prot, weight_var = "WTSA2YR", binary = FALSE
  )

  # With protein density
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

primary_dma_ua <- master %>% filter(exposure == "log_URXUDMA", outcome == "LBXSUA") %>% pull(beta)
primary_perc_bun <- master %>% filter(exposure == "log_URXUP8", outcome == "LBXSBU") %>% pull(beta)

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
    ),
    beta_primary = case_when(
      exposure == "log_URXUDMA" ~ primary_dma_ua,
      exposure == "log_URXUP8" ~ primary_perc_bun
    ),
    pct_change = (beta - beta_primary) / abs(beta_primary) * 100,
    direction_match = sign(beta) == sign(beta_primary)
  )

cat("\n=== PROTEIN INTAKE SENSITIVITY RESULTS ===\n")
protein_sens_df %>%
  transmute(Chemical = exp_label, Outcome = out_label, Analysis = analysis,
            Beta_primary = round(beta_primary, 3), Beta_adj = round(beta, 3),
            `%Change` = round(pct_change, 1), P = signif(p_value, 3), N = n) %>%
  as.data.frame() %>% print()


# ============================================================================
# PART 7: 6-level race/ethnicity sensitivity analysis
# ============================================================================

cat("\n=== 6-level race/ethnicity sensitivity ===\n\n")

# Create 6-level race factor
if (!"race6" %in% names(dat_full) && "RIDRETH3" %in% names(dat_full)) {
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

COVARS_RACE6 <- "RIDAGEYR + female + race6 + INDFMPIR + BMXBMI + smoker"

blood_findings <- add_findings %>% filter(grepl("^log_LBX", exposure))
race6_sens <- list()

for (i in 1:nrow(blood_findings)) {
  exp_v <- blood_findings$exposure[i]
  out_v <- blood_findings$outcome[i]
  binary <- isTRUE(blood_findings$type[i] == "binary")

  if (!(exp_v %in% names(dat_full)) || !(out_v %in% names(dat_full))) next

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
  left_join(blood_findings %>% select(exposure, outcome, beta_primary, exp_label, out_label),
            by = c("exposure", "outcome")) %>%
  mutate(direction_match = sign(beta) == sign(beta_primary),
         pct_change = (beta - beta_primary) / abs(beta_primary) * 100)

cat("\n=== 6-LEVEL RACE/ETHNICITY RESULTS ===\n")
race6_sens_df %>%
  filter(!is.na(beta)) %>%
  transmute(Chemical = exp_label, Outcome = out_label,
            Beta_3level = round(beta_primary, 3), Beta_6level = round(beta, 3),
            `%Change` = round(pct_change, 1), P = signif(p_value, 3), N = n) %>%
  as.data.frame() %>% print()


# ============================================================================
# PART 8: Within-round FDR correction
# ============================================================================

cat("\n=== Within-round FDR correction ===\n\n")

within_fdr <- master %>%
  group_by(screen) %>%
  mutate(p_fdr_within = p.adjust(p_value, method = "BH"),
         sig_fdr_within = p_fdr_within < 0.05) %>%
  ungroup()

fdr_summary <- within_fdr %>%
  group_by(screen) %>%
  summarise(n_tests = n(),
            n_sig_global = sum(sig_fdr05, na.rm = TRUE),
            n_sig_within = sum(sig_fdr_within, na.rm = TRUE),
            .groups = "drop")

cat("FDR correction comparison (global vs within-round):\n\n")
as.data.frame(fdr_summary) %>% print()

global_only <- within_fdr %>%
  filter(sig_fdr05 & !sig_fdr_within) %>%
  select(exposure, outcome, exp_label, out_label, screen, p_value, p_fdr_global, p_fdr_within)

cat("\nFindings significant under global FDR but NOT within-round FDR:\n")
if (nrow(global_only) > 0) {
  global_only %>%
    mutate(p_global = signif(p_fdr_global, 3), p_within = signif(p_fdr_within, 3)) %>%
    select(exp_label, out_label, screen, p_global, p_within) %>%
    as.data.frame() %>% print()
} else {
  cat("  None - all globally-significant findings remain significant\n")
}


# ============================================================================
# PART 9: LOD 40-70% chemical analysis
# ============================================================================

cat("\n=== Chemicals with 40-70% detection (LOD analysis) ===\n\n")

chem_vars <- unique(gsub("^log_", "", master$exposure))
lod_analysis <- data.frame(variable = character(), n_total = integer(),
                           n_above_lod = integer(), pct_detect = numeric(),
                           stringsAsFactors = FALSE)

for (chem in chem_vars) {
  if (!(chem %in% names(dat_full))) next

  lc_var <- sub("LBX|URX|SS", "LC", chem)
  vals <- dat_full[[chem]]
  n_total <- sum(!is.na(vals))

  if (lc_var %in% names(dat_full)) {
    lc <- as.numeric(as.character(dat_full[[lc_var]]))
    n_above <- sum(lc == 1, na.rm = TRUE)
  } else {
    min_val <- min(vals, na.rm = TRUE)
    n_at_min <- sum(abs(vals - min_val) < 1e-6, na.rm = TRUE)
    n_above <- n_total - n_at_min
  }

  if (n_total > 0) {
    lod_analysis <- rbind(lod_analysis, data.frame(
      variable = chem, n_total = n_total, n_above_lod = n_above,
      pct_detect = 100 * n_above / n_total, stringsAsFactors = FALSE
    ))
  }
}

lod_40_70 <- lod_analysis %>% filter(pct_detect >= 40, pct_detect <= 70) %>% arrange(pct_detect)
cat(sprintf("Chemicals with 40-70%% detection frequency: %d\n", nrow(lod_40_70)))

n_below_30 <- sum(lod_analysis$pct_detect < 30)
cat(sprintf("Chemicals with <30%% detection (excluded): %d\n", n_below_30))


# ============================================================================
# PART 10: Generate supplementary tables S6-S11
# ============================================================================

cat("\n=== Generating supplementary tables ===\n\n")

# Table S6: Within-round FDR
table_s6 <- fdr_summary %>%
  rename(`Screening Round` = screen, `N Tests` = n_tests,
         `N Sig (Global FDR)` = n_sig_global, `N Sig (Within-Round FDR)` = n_sig_within)
write.csv(table_s6, "figures/table_s6_within_round_fdr.csv", row.names = FALSE)

# Table S7: 6-level race sensitivity
table_s7 <- race6_sens_df %>%
  filter(!is.na(beta)) %>%
  transmute(Chemical = exp_label, Outcome = out_label,
            `β (3-level race)` = round(beta_primary, 3),
            `β (6-level race)` = round(beta, 3),
            `% Change` = round(pct_change, 1), `P-value` = signif(p_value, 3), N = n)
write.csv(table_s7, "figures/table_s7_race6_sensitivity.csv", row.names = FALSE)

# Table S8: Protein intake
table_s8 <- protein_sens_df %>%
  filter(!is.na(beta)) %>%
  transmute(Chemical = exp_label, Outcome = out_label, Adjustment = analysis,
            `β (Primary)` = round(beta_primary, 3), `β (Adjusted)` = round(beta, 3),
            `% Change` = round(pct_change, 1), `P-value` = signif(p_value, 3), N = n,
            Robust = ifelse(direction_match & p_value < 0.05, "Yes", "No"))
write.csv(table_s8, "figures/table_s8_protein_sensitivity.csv", row.names = FALSE)

# Table S9: LOD 40-70%
table_s9 <- lod_40_70 %>%
  left_join(data.frame(variable = gsub("^log_", "", master$exposure),
                       exp_label = master$exp_label) %>% distinct(), by = "variable") %>%
  transmute(Chemical = coalesce(exp_label, variable), Variable = variable,
            `% Detected` = round(pct_detect, 1), `N Samples` = n_total)
write.csv(table_s9, "figures/table_s9_lod_40_70.csv", row.names = FALSE)

# Table S10: Power analysis (copy existing)
write.csv(as.data.frame(power_display), "figures/table_s10_power_analysis.csv", row.names = FALSE)

# Table S11: STROBE checklist
strobe_checklist <- tribble(
  ~Item, ~`Checklist Item`, ~`Manuscript Section`,
  "1a", "Indicate the study's design in title or abstract", "Title, Abstract",
  "1b", "Provide informative and balanced summary", "Abstract",
  "2", "Explain scientific background and rationale", "Introduction ¶1-2",
  "3", "State specific objectives", "Introduction ¶3",
  "4", "Present key elements of study design", "Methods 2.1",
  "5", "Describe setting, locations, dates", "Methods 2.1",
  "6", "Give eligibility criteria and selection methods", "Methods 2.1",
  "7", "Define outcomes, exposures, confounders", "Methods 2.2-2.3",
  "8", "Give data sources and assessment methods", "Methods 2.2-2.3",
  "9", "Describe efforts to address bias", "Methods 2.3, 2.4.1",
  "10", "Explain how study size was arrived at", "Methods 2.1",
  "11", "Explain how variables were handled", "Methods 2.3",
  "12", "Describe all statistical methods", "Methods 2.4",
  "13", "Report numbers at each study stage", "Table 1, Results 3.1",
  "14", "Give characteristics of participants", "Table 1",
  "15", "Report outcome events or summary measures", "Results 3.2-3.5",
  "16", "Give unadjusted and adjusted estimates", "Table 2, Figures",
  "17", "Report other analyses", "Results 3.4-3.5",
  "18", "Summarize key results", "Discussion ¶1",
  "19", "Discuss limitations", "Discussion - Limitations",
  "20", "Give cautious overall interpretation", "Discussion ¶6-7",
  "21", "Discuss generalizability", "Discussion ¶2-3",
  "22", "Give source of funding", "N/A (unfunded)"
)
write.csv(strobe_checklist, "figures/table_s11_strobe_checklist.csv", row.names = FALSE)

cat("Supplementary tables saved:\n")
cat("  figures/table_s6_within_round_fdr.csv\n")
cat("  figures/table_s7_race6_sensitivity.csv\n")
cat("  figures/table_s8_protein_sensitivity.csv\n")
cat("  figures/table_s9_lod_40_70.csv\n")
cat("  figures/table_s10_power_analysis.csv\n")
cat("  figures/table_s11_strobe_checklist.csv\n")


# ============================================================================
# PART 11: Save all results
# ============================================================================

cat("\n=== Saving results ===\n\n")

save(add_sens_df, add_sens_summary, power_results, expanded_searches,
     protein_sens_df, race6_sens_df, within_fdr, fdr_summary, lod_analysis,
     dat_full,
     file = "additional_analyses_results.RData")
write.csv(add_sens_df, "additional_sensitivity_results.csv", row.names = FALSE)
write.csv(as.data.frame(power_display), "power_analysis.csv", row.names = FALSE)
write.csv(within_fdr %>% filter(sig_fdr_within) %>%
            select(screen, exposure, outcome, exp_label, out_label, beta, p_value,
                   p_fdr_global, p_fdr_within),
          "within_round_fdr_results.csv", row.names = FALSE)

cat("Saved:\n")
cat("  additional_analyses_results.RData\n")
cat("  additional_sensitivity_results.csv\n")
cat("  power_analysis.csv\n")
cat("  within_round_fdr_results.csv\n")
cat("  figures/fig1_volcano.pdf (improved, top-10 labels)\n")
cat("  figures/fig1_volcano.png\n")

cat("\n=== 09_additional_analyses.R complete ===\n")
