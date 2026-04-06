# ============================================================================
# 12_perchlorate_bun_focused.R
# Focused analysis: Perchlorate--BUN association
#
# Generates all new analyses, figures, and tables for the standalone
# manuscript on the perchlorate--BUN association.
#
# Run from the repository directory (same level as 00_functions.R).
# Requires parent .RData files in the parent directory (../):
#   exwas_novelty_screen.RData, nhanes_2015_2016_validation.RData,
#   validation_results.RData, sensitivity_results.RData,
#   additional_analyses_results.RData, dose_response_results.RData
# ============================================================================

# --- Setup -------------------------------------------------------------------

library(tidyverse)
library(survey)
library(broom)
library(splines)

for (pkg in c("tableone", "patchwork")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
}
library(tableone)
library(patchwork)

# Paths -- relative to the repository directory
P08_DIR <- ".."          # Parent directory contains .RData files
OUT_DIR <- "perchlorate_bun_outputs"

dir.create(file.path(OUT_DIR, "figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUT_DIR, "tables"), recursive = TRUE, showWarnings = FALSE)

# Source shared functions (same directory)
source("00_functions.R")

# --- Load P08 data -----------------------------------------------------------

cat("=== Loading P08 data ===\n\n")

load(file.path(P08_DIR, "exwas_novelty_screen.RData"))       # dat
dat_discovery <- dat; rm(dat)
cat(sprintf("Discovery data (2017-2018): %d obs x %d cols\n",
            nrow(dat_discovery), ncol(dat_discovery)))

load(file.path(P08_DIR, "nhanes_2015_2016_validation.RData")) # val_data
cat(sprintf("Validation data (2015-2016): %d obs x %d cols\n",
            nrow(val_data), ncol(val_data)))

load(file.path(P08_DIR, "validation_results.RData"))           # validated
load(file.path(P08_DIR, "sensitivity_results.RData"))          # sens_df, sens_summary
load(file.path(P08_DIR, "dose_response_results.RData"))        # dr_results, dr_summary

# Load additional analyses (includes dat_full with fish/PA/alcohol/creatinine/protein)
load(file.path(P08_DIR, "additional_analyses_results.RData"))
if (exists("dat_full")) {
  dat_discovery <- dat_full
  rm(dat_full)
  cat("Using dat_full from additional_analyses (includes fish/PA/alcohol/creatinine/protein)\n")
}

cat("\n")

# --- Helper: survey design builder ------------------------------------------

make_design <- function(data, weight_var = "WTSA2YR") {
  svydesign(
    id      = ~SDMVPSU,
    strata  = ~SDMVSTRA,
    weights = as.formula(paste0("~", weight_var)),
    nest    = TRUE,
    data    = data
  )
}

# --- Helper: complete-case subset for perchlorate--BUN ----------------------

perc_subset <- function(data, extra_vars = NULL, weight_var = "WTSA2YR") {
  core_vars <- c("log_URXUP8", "URXUP8", "LBXSBU", "RIDAGEYR", "female", "race3",
                  "INDFMPIR", "BMXBMI", "smoker",
                  "SDMVPSU", "SDMVSTRA", weight_var)
  all_vars <- unique(c(core_vars, extra_vars))
  all_vars <- all_vars[all_vars %in% names(data)]
  sub <- data[complete.cases(data[, all_vars]), ]
  sub <- sub[sub[[weight_var]] > 0, ]
  return(sub)
}

# Ensure log transforms exist
if (!"log_URXUP8" %in% names(dat_discovery)) {
  dat_discovery$log_URXUP8 <- log(pmax(dat_discovery$URXUP8, 0.001))
}
if (!"log_URXNO3" %in% names(dat_discovery) && "URXNO3" %in% names(dat_discovery)) {
  dat_discovery$log_URXNO3 <- log(pmax(dat_discovery$URXNO3, 0.001))
}
if (!"log_URXSCN" %in% names(dat_discovery) && "URXSCN" %in% names(dat_discovery)) {
  dat_discovery$log_URXSCN <- log(pmax(dat_discovery$URXSCN, 0.001))
}

if (!"log_URXUP8" %in% names(val_data)) {
  val_data$log_URXUP8 <- log(pmax(val_data$URXUP8, 0.001))
}


# ============================================================================
# ANALYSIS 1: Table 1 (descriptive statistics by perchlorate quartile)
# ============================================================================

cat("=== ANALYSIS 1: Table 1 ===\n\n")

disc_sub <- perc_subset(dat_discovery)
disc_sub$perc_q <- ntile(disc_sub$log_URXUP8, 4)
disc_sub$perc_q_factor <- factor(disc_sub$perc_q)

cat(sprintf("Analytic sample (discovery): n = %d\n", nrow(disc_sub)))

# Quartile boundaries
q_cuts <- quantile(disc_sub$URXUP8, probs = c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE)
cat("Perchlorate quartile boundaries (ug/L):\n")
print(round(q_cuts, 2))

# Table 1 using tableone
tab1_vars <- intersect(
  c("RIDAGEYR", "female", "BMXBMI", "INDFMPIR", "smoker",
    "LBXSBU", "eGFR", "LBXSCR", "URXUIO", "race3"),
  names(disc_sub)
)
cat_vars <- intersect(c("female", "smoker", "race3"), names(disc_sub))

tab1 <- CreateTableOne(vars = tab1_vars, strata = "perc_q", data = disc_sub,
                        factorVars = cat_vars, test = TRUE)
tab1_csv <- print(tab1, showAllLevels = TRUE, printToggle = FALSE,
                   noSpaces = TRUE, smd = FALSE)
write.csv(tab1_csv, file.path(OUT_DIR, "tables", "table1_discovery.csv"))
cat("Saved Table 1\n")

# Validation cohort
val_sub <- val_data[complete.cases(val_data[, intersect(
  c("log_URXUP8", "LBXSBU", "RIDAGEYR", "female", "race3",
    "INDFMPIR", "BMXBMI", "smoker", "SDMVPSU", "SDMVSTRA", "WTSA2YR"),
  names(val_data))]), ]
val_sub <- val_sub[val_sub$WTSA2YR > 0, ]
val_sub$perc_q <- ntile(val_sub$log_URXUP8, 4)
val_sub$perc_q_factor <- factor(val_sub$perc_q)

tab1_val_vars <- intersect(tab1_vars, names(val_sub))
cat_vars_val <- intersect(cat_vars, names(val_sub))
tab1_val <- CreateTableOne(vars = tab1_val_vars, strata = "perc_q", data = val_sub,
                            factorVars = cat_vars_val, test = TRUE)
tab1_val_csv <- print(tab1_val, showAllLevels = TRUE, printToggle = FALSE,
                       noSpaces = TRUE, smd = FALSE)
write.csv(tab1_val_csv, file.path(OUT_DIR, "tables", "table_s1_validation.csv"))
cat("Saved Table S1 (validation)\n\n")


# ============================================================================
# ANALYSIS 2: Natural cubic splines (non-linearity assessment)
# ============================================================================

cat("=== ANALYSIS 2: Natural cubic splines ===\n\n")

# Compute knots and create basis columns
knots_perc <- quantile(disc_sub$log_URXUP8, probs = c(0.25, 0.50, 0.75))
boundary_knots <- range(disc_sub$log_URXUP8)

ns_basis <- ns(disc_sub$log_URXUP8, knots = knots_perc,
               Boundary.knots = boundary_knots)
disc_sub$ns1 <- ns_basis[, 1]
disc_sub$ns2 <- ns_basis[, 2]
disc_sub$ns3 <- ns_basis[, 3]
disc_sub$ns4 <- ns_basis[, 4]

disc_sub$race3 <- factor(disc_sub$race3)
disc_sub$smoker <- factor(disc_sub$smoker)
disc_sub$female <- factor(disc_sub$female)

des <- make_design(disc_sub)

mod_spline <- svyglm(
  LBXSBU ~ ns1 + ns2 + ns3 + ns4 + RIDAGEYR + female + race3 +
    INDFMPIR + BMXBMI + smoker,
  design = des
)

# Non-linearity test
nonlin_test <- regTermTest(mod_spline, ~ ns2 + ns3 + ns4)
cat(sprintf("Non-linearity test: F = %.3f (df1 = %d, df2 = %d), p = %.4f\n",
            nonlin_test$Ftest, nonlin_test$df, nonlin_test$ddf, nonlin_test$p))

# Predictions for spline curve
pred_x <- seq(boundary_knots[1], boundary_knots[2], length.out = 200)
ns_pred <- ns(pred_x, knots = knots_perc, Boundary.knots = boundary_knots)

coefs <- coef(mod_spline)
V <- vcov(mod_spline)
ns_names <- c("ns1", "ns2", "ns3", "ns4")
ns_coefs <- coefs[ns_names]
ns_vcov <- V[ns_names, ns_names]

ns_pred_mat <- as.matrix(ns_pred)
colnames(ns_pred_mat) <- ns_names

partial_fit <- ns_pred_mat %*% ns_coefs
partial_se <- sqrt(diag(ns_pred_mat %*% ns_vcov %*% t(ns_pred_mat)))

median_x <- median(disc_sub$log_URXUP8)
median_idx <- which.min(abs(pred_x - median_x))
ref_fit <- partial_fit[median_idx]

pred_grid <- data.frame(
  log_URXUP8 = pred_x,
  fit_centred = as.numeric(partial_fit - ref_fit),
  se = as.numeric(partial_se),
  ci_lower_centred = as.numeric(partial_fit - ref_fit - 1.96 * partial_se),
  ci_upper_centred = as.numeric(partial_fit - ref_fit + 1.96 * partial_se)
)

write.csv(pred_grid, file.path(OUT_DIR, "tables", "spline_predictions.csv"),
          row.names = FALSE)
cat("Saved spline predictions\n\n")


# ============================================================================
# ANALYSIS 3: Interaction analyses (sex, age)
# ============================================================================

cat("=== ANALYSIS 3: Interaction analyses ===\n\n")

disc_sub$female_num <- as.numeric(as.character(disc_sub$female))
des <- make_design(disc_sub)

# Sex interaction
mod_sex_int <- svyglm(
  LBXSBU ~ log_URXUP8 * female_num + RIDAGEYR + race3 + INDFMPIR + BMXBMI + smoker,
  design = des
)
sex_ct <- summary(mod_sex_int)$coefficients
cat(sprintf("Sex interaction: beta = %.3f, p = %.4f\n",
            sex_ct["log_URXUP8:female_num", "Estimate"],
            sex_ct["log_URXUP8:female_num", "Pr(>|t|)"]))

# Stratum-specific
mod_male <- svyglm(LBXSBU ~ log_URXUP8 + RIDAGEYR + race3 + INDFMPIR + BMXBMI + smoker,
                    design = subset(des, female_num == 0))
mod_female <- svyglm(LBXSBU ~ log_URXUP8 + RIDAGEYR + race3 + INDFMPIR + BMXBMI + smoker,
                      design = subset(des, female_num == 1))

male_c <- summary(mod_male)$coefficients["log_URXUP8", ]
female_c <- summary(mod_female)$coefficients["log_URXUP8", ]
cat(sprintf("  Males:   beta = %.3f (SE = %.3f), p = %.4f\n",
            male_c["Estimate"], male_c["Std. Error"], male_c["Pr(>|t|)"]))
cat(sprintf("  Females: beta = %.3f (SE = %.3f), p = %.4f\n",
            female_c["Estimate"], female_c["Std. Error"], female_c["Pr(>|t|)"]))

# Age interaction
mod_age_int <- svyglm(
  LBXSBU ~ log_URXUP8 * RIDAGEYR + female + race3 + INDFMPIR + BMXBMI + smoker,
  design = des
)
age_ct <- summary(mod_age_int)$coefficients
cat(sprintf("\nAge interaction: beta = %.4f, p = %.4f\n",
            age_ct["log_URXUP8:RIDAGEYR", "Estimate"],
            age_ct["log_URXUP8:RIDAGEYR", "Pr(>|t|)"]))

# Age strata
mod_young <- svyglm(LBXSBU ~ log_URXUP8 + female + race3 + INDFMPIR + BMXBMI + smoker,
                     design = subset(des, RIDAGEYR < 50))
mod_old <- svyglm(LBXSBU ~ log_URXUP8 + female + race3 + INDFMPIR + BMXBMI + smoker,
                   design = subset(des, RIDAGEYR >= 50))

young_c <- summary(mod_young)$coefficients["log_URXUP8", ]
old_c <- summary(mod_old)$coefficients["log_URXUP8", ]

interaction_results <- data.frame(
  interaction = c("Sex (female)", "Age (continuous)"),
  interaction_beta = c(sex_ct["log_URXUP8:female_num", "Estimate"],
                       age_ct["log_URXUP8:RIDAGEYR", "Estimate"]),
  interaction_p = c(sex_ct["log_URXUP8:female_num", "Pr(>|t|)"],
                    age_ct["log_URXUP8:RIDAGEYR", "Pr(>|t|)"]),
  stratum1_label = c("Males", "Age < 50"),
  stratum1_beta = c(male_c["Estimate"], young_c["Estimate"]),
  stratum1_se = c(male_c["Std. Error"], young_c["Std. Error"]),
  stratum1_p = c(male_c["Pr(>|t|)"], young_c["Pr(>|t|)"]),
  stratum2_label = c("Females", "Age >= 50"),
  stratum2_beta = c(female_c["Estimate"], old_c["Estimate"]),
  stratum2_se = c(female_c["Std. Error"], old_c["Std. Error"]),
  stratum2_p = c(female_c["Pr(>|t|)"], old_c["Pr(>|t|)"])
)

write.csv(interaction_results, file.path(OUT_DIR, "tables", "interaction_results.csv"),
          row.names = FALSE)
cat("Saved interaction results\n\n")


# ============================================================================
# ANALYSIS 4: Iodine status interaction
# ============================================================================

cat("=== ANALYSIS 4: Iodine status interaction ===\n\n")

if ("URXUIO" %in% names(disc_sub)) {
  disc_sub$iodine_cat <- cut(disc_sub$URXUIO,
    breaks = c(0, 100, 300, Inf),
    labels = c("Deficient (<100)", "Sufficient (100-299)", "Excessive (>=300)"),
    right = FALSE)

  cat("Iodine status distribution:\n")
  print(table(disc_sub$iodine_cat, useNA = "ifany"))

  des_iod <- make_design(disc_sub[!is.na(disc_sub$iodine_cat), ])

  # Interaction test
  mod_iod <- svyglm(
    LBXSBU ~ log_URXUP8 * iodine_cat + RIDAGEYR + female + race3 +
      INDFMPIR + BMXBMI + smoker,
    design = des_iod)
  iod_ct <- summary(mod_iod)$coefficients

  # Stratified
  iodine_strata <- list()
  for (lvl in levels(disc_sub$iodine_cat)) {
    mod_sub <- tryCatch(
      svyglm(LBXSBU ~ log_URXUP8 + RIDAGEYR + female + race3 + INDFMPIR + BMXBMI + smoker,
             design = subset(des_iod, iodine_cat == lvl)),
      error = function(e) NULL)
    if (!is.null(mod_sub)) {
      cs <- summary(mod_sub)$coefficients["log_URXUP8", ]
      iodine_strata[[lvl]] <- data.frame(
        iodine_status = lvl, beta = cs["Estimate"], se = cs["Std. Error"],
        p_value = cs["Pr(>|t|)"], n = nobs(mod_sub), stringsAsFactors = FALSE)
      cat(sprintf("  %s: beta = %.3f, p = %.4f, n = %d\n",
                  lvl, cs["Estimate"], cs["Pr(>|t|)"], nobs(mod_sub)))
    }
  }

  iodine_results <- bind_rows(iodine_strata)
  write.csv(iodine_results, file.path(OUT_DIR, "tables", "table_s3_iodine.csv"),
            row.names = FALSE)
  cat("Saved iodine results\n\n")
} else {
  cat("URXUIO not available. Skipping.\n\n")
}


# ============================================================================
# ANALYSIS 5: Co-NIS-inhibitor adjustment (nitrate, thiocyanate)
# ============================================================================

cat("=== ANALYSIS 5: Co-NIS-inhibitor adjustment ===\n\n")

if (all(c("log_URXNO3", "log_URXSCN") %in% names(disc_sub))) {
  nis_sub <- disc_sub[complete.cases(disc_sub[, c("log_URXNO3", "log_URXSCN")]), ]
  des_nis <- make_design(nis_sub)
  cat(sprintf("NIS-complete sample: n = %d\n", nrow(nis_sub)))

  mod_base <- svyglm(LBXSBU ~ log_URXUP8 + RIDAGEYR + female + race3 + INDFMPIR + BMXBMI + smoker,
                      design = des_nis)
  mod_no3 <- svyglm(LBXSBU ~ log_URXUP8 + log_URXNO3 + RIDAGEYR + female + race3 + INDFMPIR + BMXBMI + smoker,
                     design = des_nis)
  mod_scn <- svyglm(LBXSBU ~ log_URXUP8 + log_URXSCN + RIDAGEYR + female + race3 + INDFMPIR + BMXBMI + smoker,
                     design = des_nis)
  mod_both <- svyglm(LBXSBU ~ log_URXUP8 + log_URXNO3 + log_URXSCN + RIDAGEYR + female + race3 + INDFMPIR + BMXBMI + smoker,
                      design = des_nis)

  get_b <- function(m) summary(m)$coefficients["log_URXUP8", ]

  nis_results <- data.frame(
    model = c("Primary (NIS-complete)", "Adjusted for nitrate",
              "Adjusted for thiocyanate", "Adjusted for both"),
    beta = c(get_b(mod_base)["Estimate"], get_b(mod_no3)["Estimate"],
             get_b(mod_scn)["Estimate"], get_b(mod_both)["Estimate"]),
    se = c(get_b(mod_base)["Std. Error"], get_b(mod_no3)["Std. Error"],
           get_b(mod_scn)["Std. Error"], get_b(mod_both)["Std. Error"]),
    p_value = c(get_b(mod_base)["Pr(>|t|)"], get_b(mod_no3)["Pr(>|t|)"],
                get_b(mod_scn)["Pr(>|t|)"], get_b(mod_both)["Pr(>|t|)"]),
    n = nrow(nis_sub)
  )
  nis_results$pct_change <- (nis_results$beta - nis_results$beta[1]) / abs(nis_results$beta[1]) * 100

  cat("\nCo-NIS-inhibitor results:\n")
  print(as.data.frame(nis_results))

  write.csv(nis_results, file.path(OUT_DIR, "tables", "table_s4_nis.csv"),
            row.names = FALSE)
  cat("Saved NIS results\n\n")
} else {
  cat("Nitrate/thiocyanate not available. Skipping.\n\n")
}


# ============================================================================
# ANALYSIS 6: Other kidney markers (specificity)
# ============================================================================

cat("=== ANALYSIS 6: Kidney marker specificity ===\n\n")

des <- make_design(disc_sub)

kidney_markers <- list(
  list(var = "LBXSCR", label = "Serum creatinine", transform = FALSE),
  list(var = "eGFR", label = "eGFR (CKD-EPI 2021)", transform = FALSE),
  list(var = "URDACT", label = "Albumin-creatinine ratio", transform = TRUE)
)

kidney_results <- list()
for (km in kidney_markers) {
  if (!(km$var %in% names(disc_sub))) {
    cat(sprintf("  %s: not available\n", km$label)); next
  }
  out_var <- if (km$transform) {
    disc_sub[[paste0("log_", km$var)]] <- log(pmax(disc_sub[[km$var]], 0.01))
    paste0("log_", km$var)
  } else km$var

  sub <- disc_sub[!is.na(disc_sub[[out_var]]), ]
  des_k <- make_design(sub)
  mod_k <- svyglm(as.formula(paste(out_var, "~ log_URXUP8 + RIDAGEYR + female + race3 + INDFMPIR + BMXBMI + smoker")),
                   design = des_k)
  kc <- summary(mod_k)$coefficients["log_URXUP8", ]
  kidney_results[[km$label]] <- data.frame(
    outcome = km$label, beta = kc["Estimate"], se = kc["Std. Error"],
    p_value = kc["Pr(>|t|)"], n = nobs(mod_k), stringsAsFactors = FALSE)
  cat(sprintf("  %s: beta = %.4f, p = %.4f\n", km$label, kc["Estimate"], kc["Pr(>|t|)"]))
}

# Add BUN for comparison
bun_mod <- svyglm(LBXSBU ~ log_URXUP8 + RIDAGEYR + female + race3 + INDFMPIR + BMXBMI + smoker,
                   design = des)
bun_c <- summary(bun_mod)$coefficients["log_URXUP8", ]
kidney_results[["BUN (primary)"]] <- data.frame(
  outcome = "BUN (primary)", beta = bun_c["Estimate"], se = bun_c["Std. Error"],
  p_value = bun_c["Pr(>|t|)"], n = nobs(bun_mod), stringsAsFactors = FALSE)

kidney_df <- bind_rows(kidney_results)
write.csv(kidney_df, file.path(OUT_DIR, "tables", "table_s2_kidney_markers.csv"),
          row.names = FALSE)
cat("Saved kidney marker results\n\n")


# ============================================================================
# ANALYSIS 7: Population-level considerations
# ============================================================================

cat("=== ANALYSIS 7: Population-level considerations ===\n\n")

q_bun <- svyby(~LBXSBU, ~perc_q, design = des, svymean, na.rm = TRUE)
cat("Survey-weighted mean BUN by perchlorate quartile:\n")
print(q_bun)

diff_bun <- q_bun$LBXSBU[4] - q_bun$LBXSBU[1]
cat(sprintf("\nQ4 - Q1 difference: %.2f mg/dL (%.1f%% of clinical range 7-20)\n",
            diff_bun, diff_bun / 13 * 100))

disc_sub$bun_elevated <- as.numeric(disc_sub$LBXSBU > 20)
des_e <- make_design(disc_sub)
q_elev <- svyby(~bun_elevated, ~perc_q, design = des_e, svymean, na.rm = TRUE)
cat("\nProportion with BUN > 20 mg/dL by quartile:\n")
print(q_elev)

pop_attrib <- data.frame(
  quartile = 1:4, mean_bun = q_bun$LBXSBU, se_bun = q_bun$se,
  pct_elevated = q_elev$bun_elevated * 100, se_elevated = q_elev$se * 100)
pop_attrib$diff_from_q1 <- pop_attrib$mean_bun - pop_attrib$mean_bun[1]

write.csv(pop_attrib, file.path(OUT_DIR, "tables", "population_attributable.csv"),
          row.names = FALSE)
cat("Saved population attributable results\n\n")


# ============================================================================
# ANALYSIS 8: Compile P08 results for perchlorate--BUN
# ============================================================================

cat("=== ANALYSIS 8: Extracting P08 perchlorate-BUN results ===\n\n")

master <- read.csv(file.path(P08_DIR, "exwas_master_results.csv"))
primary <- master %>% filter(exposure == "log_URXUP8", outcome == "LBXSBU")

validation <- read.csv(file.path(P08_DIR, "validation_results.csv"))
val_result <- validation %>% filter(exposure == "log_URXUP8", outcome == "LBXSBU")

sens <- read.csv(file.path(P08_DIR, "sensitivity_results.csv"))
perc_sens <- sens %>% filter(exposure == "log_URXUP8", outcome == "LBXSBU")

add_sens <- read.csv(file.path(P08_DIR, "additional_sensitivity_results.csv"))
perc_add <- add_sens %>% filter(exposure == "log_URXUP8", outcome == "LBXSBU")

egfr_sens <- read.csv(file.path(P08_DIR, "figures", "table_s15_egfr_adjusted_perchlorate_bun.csv"))
evalue <- read.csv(file.path(P08_DIR, "figures", "evalue_analysis.csv"))

table2 <- data.frame(
  cohort = c("Discovery (2017-2018)", "Validation (2015-2016)"),
  beta = c(primary$beta, val_result$beta_val),
  se = c(primary$se, val_result$se_val),
  ci_lower = c(primary$beta - 1.96 * primary$se,
               val_result$beta_val - 1.96 * val_result$se_val),
  ci_upper = c(primary$beta + 1.96 * primary$se,
               val_result$beta_val + 1.96 * val_result$se_val),
  p_value = c(primary$p_value, val_result$p_val),
  n = c(primary$n, val_result$n_val)
)
write.csv(table2, file.path(OUT_DIR, "tables", "table2_primary_validation.csv"),
          row.names = FALSE)

all_sens_compiled <- bind_rows(
  perc_sens %>% transmute(category = "Demographic/Model", analysis, beta, se, p_value, n,
                           pct_change = (beta - primary$beta) / abs(primary$beta) * 100),
  perc_add %>% transmute(category = "Dietary/Lifestyle", analysis, beta, se, p_value, n, pct_change),
  egfr_sens %>% transmute(category = "Clinical", analysis = model, beta, se, p_value, n, pct_change)
)
write.csv(all_sens_compiled, file.path(OUT_DIR, "tables", "table4_sensitivity.csv"),
          row.names = FALSE)

cat(sprintf("Primary: beta = %.3f (95%% CI: %.2f, %.2f), p = %.2e\n",
            primary$beta, primary$beta - 1.96 * primary$se,
            primary$beta + 1.96 * primary$se, primary$p_value))
cat(sprintf("Validation: beta = %.3f (95%% CI: %.2f, %.2f), p = %.4f\n",
            val_result$beta_val, val_result$beta_val - 1.96 * val_result$se_val,
            val_result$beta_val + 1.96 * val_result$se_val, val_result$p_val))
cat("Saved compiled results\n\n")


# ============================================================================
# FIGURES
# ============================================================================

cat("=== Generating figures ===\n\n")

# --- Figure 1: Dose-response panel ---

# Panel A: Quartile plot
perc_dr_detail <- NULL
for (r in dr_results) {
  if (r$exposure == "log_URXUP8" && r$outcome == "LBXSBU") {
    perc_dr_detail <- r; break
  }
}

# Discovery quartiles
if (!is.null(perc_dr_detail)) {
  dr_disc <- perc_dr_detail$q_estimates %>% mutate(cohort = "Discovery (2017-2018)")
} else {
  dr_disc <- data.frame(quartile = 1:4, beta = c(0, NA, NA, NA),
                         se = 0, ci_lower = 0, ci_upper = 0,
                         cohort = "Discovery (2017-2018)")
}

# Validation quartiles
des_val <- make_design(val_sub)
mod_vq <- svyglm(LBXSBU ~ perc_q_factor + RIDAGEYR + female + race3 + INDFMPIR + BMXBMI + smoker,
                  design = des_val)
vq_ct <- summary(mod_vq)$coefficients
vq_rows <- grep("^perc_q_factor", rownames(vq_ct))

dr_val <- data.frame(
  quartile = 1:4,
  beta = c(0, vq_ct[vq_rows, "Estimate"]),
  se = c(0, vq_ct[vq_rows, "Std. Error"]),
  ci_lower = c(0, vq_ct[vq_rows, "Estimate"] - 1.96 * vq_ct[vq_rows, "Std. Error"]),
  ci_upper = c(0, vq_ct[vq_rows, "Estimate"] + 1.96 * vq_ct[vq_rows, "Std. Error"]),
  cohort = "Validation (2015-2016)")

dr_combined <- bind_rows(
  dr_disc %>% select(quartile, beta, se, ci_lower, ci_upper, cohort),
  dr_val %>% select(quartile, beta, se, ci_lower, ci_upper, cohort))

p_quartile <- ggplot(dr_combined, aes(x = factor(quartile), y = beta,
                                       color = cohort, group = cohort)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_pointrange(aes(ymin = ci_lower, ymax = ci_upper),
                  position = position_dodge(width = 0.3), size = 0.6) +
  scale_color_manual(values = c("Discovery (2017-2018)" = "#2166AC",
                                 "Validation (2015-2016)" = "#B2182B"), name = NULL) +
  labs(x = "Perchlorate exposure quartile",
       y = expression(Delta * " BUN (mg/dL) vs Q1"),
       title = "A. Dose-response across exposure quartiles") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank(),
        plot.title = element_text(size = 13, face = "bold"))

# Panel B: Spline
q05 <- quantile(disc_sub$log_URXUP8, 0.05)
q95 <- quantile(disc_sub$log_URXUP8, 0.95)
spline_trim <- pred_grid %>% filter(log_URXUP8 >= q05 & log_URXUP8 <= q95)

p_spline <- ggplot(spline_trim, aes(x = log_URXUP8)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_ribbon(aes(ymin = ci_lower_centred, ymax = ci_upper_centred),
              fill = "#2166AC", alpha = 0.2) +
  geom_line(aes(y = fit_centred), color = "#2166AC", linewidth = 1) +
  geom_rug(data = disc_sub %>% filter(log_URXUP8 >= q05 & log_URXUP8 <= q95),
           aes(x = log_URXUP8), sides = "b", alpha = 0.05, inherit.aes = FALSE) +
  labs(x = expression("Log urinary perchlorate (" * mu * "g/L)"),
       y = expression(Delta * " BUN (mg/dL) vs median exposure"),
       title = "B. Continuous dose-response (natural cubic spline)") +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(size = 13, face = "bold"))

fig1 <- p_quartile / p_spline +
  plot_annotation(title = "Figure 1. Perchlorate and blood urea nitrogen",
                  theme = theme(plot.title = element_text(size = 14, face = "bold")))

ggsave(file.path(OUT_DIR, "figures", "fig1_dose_response.pdf"), fig1, width = 8, height = 10, dpi = 300)
ggsave(file.path(OUT_DIR, "figures", "fig1_dose_response.png"), fig1, width = 8, height = 10, dpi = 300)
cat("Saved Figure 1\n")

# --- Figure 2: Forest plot ---

primary_row <- data.frame(category = "Reference", analysis = "Primary model",
                           beta = primary$beta, se = primary$se,
                           p_value = primary$p_value, n = primary$n, pct_change = 0)

forest_data <- bind_rows(primary_row, all_sens_compiled %>% filter(analysis != "Primary (no eGFR)"))
if (exists("nis_results")) {
  nis_rows <- nis_results %>% filter(model != "Primary (NIS-complete)") %>%
    transmute(category = "Co-exposure", analysis = model, beta, se, p_value, n, pct_change)
  forest_data <- bind_rows(forest_data, nis_rows)
}

forest_data <- forest_data %>%
  filter(!is.na(beta)) %>%
  mutate(ci_lower = beta - 1.96 * se, ci_upper = beta + 1.96 * se,
         label = ifelse(!is.na(n), paste0(analysis, " (n=", n, ")"), analysis),
         label = factor(label, levels = rev(unique(label))))

cat_colours <- c("Reference" = "black", "Demographic/Model" = "#2166AC",
                 "Dietary/Lifestyle" = "#4DAF4A", "Clinical" = "#984EA3",
                 "Co-exposure" = "#FF7F00")

p_forest <- ggplot(forest_data, aes(x = beta, y = label, color = category)) +
  geom_vline(xintercept = primary$beta, linetype = "dashed", color = "grey40") +
  geom_vline(xintercept = 0, linetype = "solid", color = "grey70") +
  geom_pointrange(aes(xmin = ci_lower, xmax = ci_upper), size = 0.5) +
  scale_color_manual(values = cat_colours, name = "Category") +
  labs(x = expression(beta * " (mg/dL per log-unit perchlorate)"), y = NULL,
       title = "Figure 2. Sensitivity analyses") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank(),
        plot.title = element_text(size = 13, face = "bold"))

ggsave(file.path(OUT_DIR, "figures", "fig2_forest.pdf"), p_forest, width = 10, height = 8, dpi = 300)
ggsave(file.path(OUT_DIR, "figures", "fig2_forest.png"), p_forest, width = 10, height = 8, dpi = 300)
cat("Saved Figure 2\n")

# --- Figure 3: Kidney marker specificity ---

kidney_plot <- kidney_df %>%
  mutate(ci_lower = beta - 1.96 * se, ci_upper = beta + 1.96 * se,
         sig = ifelse(p_value < 0.05, "p < 0.05", "p >= 0.05"),
         outcome = factor(outcome, levels = rev(c("BUN (primary)", "Serum creatinine",
                                                   "eGFR (CKD-EPI 2021)",
                                                   "Albumin-creatinine ratio"))))

p_kidney <- ggplot(kidney_plot, aes(x = beta, y = outcome, color = sig)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_pointrange(aes(xmin = ci_lower, xmax = ci_upper), size = 0.7) +
  scale_color_manual(values = c("p < 0.05" = "#2166AC", "p >= 0.05" = "grey50"), name = NULL) +
  labs(x = expression(beta * " per log-unit perchlorate"), y = NULL,
       title = "Figure 3. Kidney marker specificity") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank(),
        plot.title = element_text(size = 13, face = "bold"))

ggsave(file.path(OUT_DIR, "figures", "fig3_kidney_specificity.pdf"), p_kidney, width = 8, height = 5, dpi = 300)
ggsave(file.path(OUT_DIR, "figures", "fig3_kidney_specificity.png"), p_kidney, width = 8, height = 5, dpi = 300)
cat("Saved Figure 3\n")


# ============================================================================
# SUMMARY
# ============================================================================

cat("\n=== All outputs saved to", OUT_DIR, "===\n\n")
cat("Tables:\n")
list.files(file.path(OUT_DIR, "tables")) %>% paste("  ", .) %>% cat(sep = "\n")
cat("\n\nFigures:\n")
list.files(file.path(OUT_DIR, "figures")) %>% paste("  ", .) %>% cat(sep = "\n")
cat("\n\n=== 12_perchlorate_bun_focused.R complete ===\n")
