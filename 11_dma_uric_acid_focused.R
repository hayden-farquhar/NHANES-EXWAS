# ============================================================================
# 11_dma_uric_acid_focused.R
# Focused analysis: DMA--Uric Acid association
#
# Generates all new analyses, figures, and tables for the standalone
# manuscript on the DMA--uric acid association.
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
library(nhanesA)

# Paths -- relative to the repository directory
P08_DIR <- ".."          # Parent directory contains .RData files
OUT_DIR <- "dma_uric_acid_outputs"

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
load(file.path(P08_DIR, "additional_analyses_results.RData"))  # add_sens_df, etc.
# Note: additional_analyses_results.RData also contains dat_full with extra
# variables (fish, PA, alcohol, creatinine, protein). Use that version.
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

# --- Helper: complete-case subset for DMA--uric acid ------------------------

dma_subset <- function(data, extra_vars = NULL, weight_var = "WTSA2YR") {
  core_vars <- c("log_URXUDMA", "LBXSUA", "RIDAGEYR", "female", "race3",
                  "INDFMPIR", "BMXBMI", "smoker",
                  "SDMVPSU", "SDMVSTRA", weight_var)
  all_vars <- unique(c(core_vars, extra_vars))
  all_vars <- all_vars[all_vars %in% names(data)]
  sub <- data[complete.cases(data[, all_vars]), ]
  sub <- sub[sub[[weight_var]] > 0, ]
  return(sub)
}


# ============================================================================
# ANALYSIS 1: Sample-specific Table 1 (stratified by DMA quartiles)
# ============================================================================

cat("=== ANALYSIS 1: Table 1 (weighted demographics by DMA quartile) ===\n\n")

disc_sub <- dma_subset(dat_discovery)
disc_sub$dma_quartile <- ntile(disc_sub$log_URXUDMA, 4)
disc_sub$dma_q_label <- factor(paste0("Q", disc_sub$dma_quartile),
                                levels = paste0("Q", 1:4))

cat(sprintf("Analytic sample (discovery): n = %d\n", nrow(disc_sub)))

# Compute quartile cut points (on original scale for interpretability)
if ("URXUDMA" %in% names(disc_sub)) {
  q_cuts <- quantile(disc_sub$URXUDMA, probs = c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE)
  cat("DMA quartile boundaries (original scale, ug/L):\n")
  print(round(q_cuts, 2))
}

des <- make_design(disc_sub)

# Weighted means by quartile
# For continuous variables, report mean; for binary indicators, create numeric versions
disc_sub$female_num <- as.numeric(disc_sub$female == 1 |
  as.character(disc_sub$female) == "1")
disc_sub$smoker_num <- as.numeric(as.character(disc_sub$smoker) == "current")

table1_vars <- c("RIDAGEYR", "female_num", "BMXBMI", "INDFMPIR", "smoker_num", "LBXSUA")
table1_labels <- c("Age (years)", "Female (%)", "BMI (kg/m2)",
                    "Poverty-income ratio", "Current smoker (%)", "Uric acid (mg/dL)")
table1_is_pct <- c(FALSE, TRUE, FALSE, FALSE, TRUE, FALSE)

des <- make_design(disc_sub)  # re-create with new columns

table1_rows <- list()

for (j in seq_along(table1_vars)) {
  v <- table1_vars[j]
  if (!(v %in% names(disc_sub))) next

  row_data <- list(Variable = table1_labels[j])
  mult <- if (table1_is_pct[j]) 100 else 1

  # Overall
  overall_mean <- svymean(as.formula(paste0("~", v)), des, na.rm = TRUE)
  row_data[["Overall"]] <- sprintf("%.1f", mult * coef(overall_mean))

  # By quartile
  for (q in 1:4) {
    q_des <- subset(des, dma_quartile == q)
    q_mean <- svymean(as.formula(paste0("~", v)), q_des, na.rm = TRUE)
    row_data[[paste0("Q", q)]] <- sprintf("%.1f", mult * coef(q_mean))
  }

  table1_rows[[j]] <- as.data.frame(row_data, stringsAsFactors = FALSE)
}

# Sample sizes per quartile
n_row <- list(Variable = "N (unweighted)")
n_row[["Overall"]] <- as.character(nrow(disc_sub))
for (q in 1:4) {
  n_row[[paste0("Q", q)]] <- as.character(sum(disc_sub$dma_quartile == q))
}
table1_rows <- c(list(as.data.frame(n_row, stringsAsFactors = FALSE)), table1_rows)

# DMA concentration per quartile (log-scale)
dma_row <- list(Variable = "log DMA (mean)")
dma_overall <- svymean(~log_URXUDMA, des, na.rm = TRUE)
dma_row[["Overall"]] <- sprintf("%.2f", coef(dma_overall))
for (q in 1:4) {
  q_des <- subset(des, dma_quartile == q)
  q_mean <- svymean(~log_URXUDMA, q_des, na.rm = TRUE)
  dma_row[[paste0("Q", q)]] <- sprintf("%.2f", coef(q_mean))
}
table1_rows <- c(table1_rows, list(as.data.frame(dma_row, stringsAsFactors = FALSE)))

# Hyperuricemia prevalence (primary: >6.8 mg/dL for both sexes, MSU saturation threshold)
disc_sub$hyperuricemia <- as.numeric(disc_sub$LBXSUA > 6.8)
des_hu <- make_design(disc_sub)
hu_row <- list(Variable = "Hyperuricemia (%)")
hu_overall <- svymean(~hyperuricemia, des_hu, na.rm = TRUE)
hu_row[["Overall"]] <- sprintf("%.1f", 100 * coef(hu_overall))
for (q in 1:4) {
  q_des <- subset(des_hu, dma_quartile == q)
  q_mean <- svymean(~hyperuricemia, q_des, na.rm = TRUE)
  hu_row[[paste0("Q", q)]] <- sprintf("%.1f", 100 * coef(q_mean))
}
table1_rows <- c(table1_rows, list(as.data.frame(hu_row, stringsAsFactors = FALSE)))

# Race/ethnicity distribution
race_row_nh_white <- list(Variable = "Non-Hispanic White (%)")
race_row_nh_black <- list(Variable = "Non-Hispanic Black (%)")
race_row_other <- list(Variable = "Other race/ethnicity (%)")

# Detect actual race3 levels
race3_levels <- unique(as.character(disc_sub$race3))
cat("Race3 levels in data:", paste(race3_levels, collapse = ", "), "\n")

# Create binary indicators for each level
disc_sub$nh_white <- as.numeric(grepl("White|white", as.character(disc_sub$race3)))
disc_sub$nh_black <- as.numeric(grepl("Black|black", as.character(disc_sub$race3)))
disc_sub$other_race <- as.numeric(grepl("Other|other", as.character(disc_sub$race3)))
des_race <- make_design(disc_sub)

for (rv in c("nh_white", "nh_black", "other_race")) {
  target_row <- switch(rv,
    nh_white = race_row_nh_white,
    nh_black = race_row_nh_black,
    other_race = race_row_other
  )
  ov <- svymean(as.formula(paste0("~", rv)), des_race, na.rm = TRUE)
  target_row[["Overall"]] <- sprintf("%.1f", 100 * coef(ov))
  for (q in 1:4) {
    q_des <- subset(des_race, dma_quartile == q)
    q_mean <- svymean(as.formula(paste0("~", rv)), q_des, na.rm = TRUE)
    target_row[[paste0("Q", q)]] <- sprintf("%.1f", 100 * coef(q_mean))
  }
  if (rv == "nh_white") race_row_nh_white <- target_row
  if (rv == "nh_black") race_row_nh_black <- target_row
  if (rv == "other_race") race_row_other <- target_row
}

table1_rows <- c(table1_rows,
  list(as.data.frame(race_row_nh_white, stringsAsFactors = FALSE)),
  list(as.data.frame(race_row_nh_black, stringsAsFactors = FALSE)),
  list(as.data.frame(race_row_other, stringsAsFactors = FALSE))
)

table1 <- bind_rows(table1_rows)
write.csv(table1, file.path(OUT_DIR, "tables", "table1_demographics.csv"),
          row.names = FALSE)
cat("\nTable 1 saved.\n")
print(table1)


# ============================================================================
# ANALYSIS 2: Hyperuricemia (binary outcome)
# ============================================================================

cat("\n=== ANALYSIS 2: Hyperuricemia (binary outcome) ===\n\n")

# --- 2a. Discovery cycle (2017-2018) ---

# Primary threshold: >6.8 mg/dL for both sexes (MSU saturation point)
disc_sub$hyperuricemia <- as.numeric(disc_sub$LBXSUA > 6.8)

des_hu <- make_design(disc_sub)

# Continuous log-DMA model
hu_model_cont <- svyglm(
  hyperuricemia ~ log_URXUDMA + RIDAGEYR + female + race3 + INDFMPIR + BMXBMI + smoker,
  design = des_hu, family = quasibinomial()
)
hu_cont_tidy <- tidy(hu_model_cont, conf.int = TRUE, exponentiate = TRUE)
hu_cont_row <- hu_cont_tidy %>% filter(term == "log_URXUDMA")
cat("Discovery: OR per log-unit DMA =",
    sprintf("%.2f (95%% CI: %.2f-%.2f), p = %.4f\n",
            hu_cont_row$estimate, hu_cont_row$conf.low,
            hu_cont_row$conf.high, hu_cont_row$p.value))

# Quartile model
disc_sub$dma_q_factor <- factor(disc_sub$dma_quartile)
des_hu_q <- make_design(disc_sub)

hu_model_q <- svyglm(
  hyperuricemia ~ dma_q_factor + RIDAGEYR + female + race3 + INDFMPIR + BMXBMI + smoker,
  design = des_hu_q, family = quasibinomial()
)
hu_q_tidy <- tidy(hu_model_q, conf.int = TRUE, exponentiate = TRUE)
hu_q_rows <- hu_q_tidy %>% filter(grepl("dma_q_factor", term))

cat("\nDiscovery quartile ORs (vs Q1):\n")
hu_q_results <- bind_rows(
  tibble(Quartile = "Q1 (ref)", OR = 1.00, CI_low = NA, CI_high = NA, p = NA),
  hu_q_rows %>% transmute(
    Quartile = gsub("dma_q_factor", "Q", term),
    OR = round(estimate, 2),
    CI_low = round(conf.low, 2),
    CI_high = round(conf.high, 2),
    p = round(p.value, 4)
  )
)
print(hu_q_results)

# Trend test
disc_sub$dma_q_numeric <- disc_sub$dma_quartile
des_hu_trend <- make_design(disc_sub)
hu_trend <- svyglm(
  hyperuricemia ~ dma_q_numeric + RIDAGEYR + female + race3 + INDFMPIR + BMXBMI + smoker,
  design = des_hu_trend, family = quasibinomial()
)
hu_trend_p <- summary(hu_trend)$coefficients["dma_q_numeric", "Pr(>|t|)"]
cat(sprintf("Trend test p = %.4f\n", hu_trend_p))

# --- 2b. Validation cycle (2015-2016) ---

val_sub <- dma_subset(val_data)
val_sub$hyperuricemia <- as.numeric(val_sub$LBXSUA > 6.8)

des_val_hu <- make_design(val_sub)

hu_val_model <- svyglm(
  hyperuricemia ~ log_URXUDMA + RIDAGEYR + female + race3 + INDFMPIR + BMXBMI + smoker,
  design = des_val_hu, family = quasibinomial()
)
hu_val_tidy <- tidy(hu_val_model, conf.int = TRUE, exponentiate = TRUE)
hu_val_row <- hu_val_tidy %>% filter(term == "log_URXUDMA")
cat(sprintf("\nValidation: OR per log-unit DMA = %.2f (95%% CI: %.2f-%.2f), p = %.4f\n",
            hu_val_row$estimate, hu_val_row$conf.low,
            hu_val_row$conf.high, hu_val_row$p.value))

# Save hyperuricemia results
hu_results <- tibble(
  Cycle = c("2017-2018", "2015-2016"),
  OR = c(hu_cont_row$estimate, hu_val_row$estimate),
  CI_low = c(hu_cont_row$conf.low, hu_val_row$conf.low),
  CI_high = c(hu_cont_row$conf.high, hu_val_row$conf.high),
  p_value = c(hu_cont_row$p.value, hu_val_row$p.value),
  N = c(nrow(disc_sub), nrow(val_sub))
)

write.csv(hu_results,
          file.path(OUT_DIR, "tables", "table_s1_hyperuricemia.csv"),
          row.names = FALSE)
write.csv(hu_q_results,
          file.path(OUT_DIR, "tables", "table_s1_hyperuricemia_quartiles.csv"),
          row.names = FALSE)
cat("Hyperuricemia results saved.\n")


# ============================================================================
# ANALYSIS 3: Restricted cubic spline
# ============================================================================

cat("\n=== ANALYSIS 3: Restricted cubic spline ===\n\n")

if (!requireNamespace("rms", quietly = TRUE)) {
  install.packages("rms", repos = "https://cloud.r-project.org")
}
library(rms)

# Fit RCS with 4 knots using ns() from splines (more compatible with svyglm/predict)
disc_rcs <- disc_sub
disc_rcs$wt <- disc_rcs$WTSA2YR

des_rcs <- make_design(disc_rcs)

library(splines)

# Use ns() (natural splines) with 3 df (similar flexibility to rcs with 4 knots)
rcs_model <- svyglm(
  LBXSUA ~ ns(log_URXUDMA, df = 3) + RIDAGEYR + female + race3 + INDFMPIR + BMXBMI + smoker,
  design = des_rcs
)

# Linear model for comparison
linear_model <- svyglm(
  LBXSUA ~ log_URXUDMA + RIDAGEYR + female + race3 + INDFMPIR + BMXBMI + smoker,
  design = des_rcs
)

# Test for non-linearity via Wald test on spline terms
rcs_anova <- tryCatch({
  regTermTest(rcs_model, ~ ns(log_URXUDMA, df = 3), df = NULL)
}, error = function(e) NULL)

if (!is.null(rcs_anova)) {
  cat(sprintf("Overall association test: F = %.2f, p = %.4f\n",
              rcs_anova$Ftest[1], rcs_anova$p[1]))
}

# Test non-linearity: compare linear vs spline
nonlin_test <- tryCatch({
  anova(linear_model, rcs_model)
}, error = function(e) {
  cat("Note: anova comparison not available; using AIC\n")
  cat(sprintf("Linear AIC: %.1f, Spline AIC: %.1f\n", AIC(linear_model), AIC(rcs_model)))
  NULL
})

# Generate predictions across DMA range
# Key: use the ns() basis constructed from the TRAINING data
dma_range <- range(disc_rcs$log_URXUDMA, na.rm = TRUE)
dma_seq <- seq(dma_range[1], dma_range[2], length.out = 200)

# Build the same ns basis using the knots from the training fit
ns_attr <- attributes(ns(disc_rcs$log_URXUDMA, df = 3))
knots_internal <- ns_attr$knots
boundary_knots <- ns_attr$Boundary.knots

# Create prediction frame with proper factor levels
pred_data <- data.frame(
  log_URXUDMA = dma_seq,
  RIDAGEYR = mean(disc_rcs$RIDAGEYR, na.rm = TRUE),
  female = 0,
  race3 = factor(levels(disc_rcs$race3)[1], levels = levels(disc_rcs$race3)),
  INDFMPIR = mean(disc_rcs$INDFMPIR, na.rm = TRUE),
  BMXBMI = mean(disc_rcs$BMXBMI, na.rm = TRUE),
  smoker = factor(levels(disc_rcs$smoker)[1], levels = levels(disc_rcs$smoker))
)

# Predict using model matrix + vcov approach (avoids predict.svyglm issues)
mm_pred <- model.matrix(delete.response(terms(rcs_model)), data = pred_data)
V <- vcov(rcs_model)
common_cols <- intersect(colnames(mm_pred), colnames(V))
mm_pred <- mm_pred[, common_cols, drop = FALSE]
V <- V[common_cols, common_cols, drop = FALSE]

pred_data$fit <- as.numeric(mm_pred %*% coef(rcs_model)[common_cols])
pred_data$se <- sqrt(rowSums((mm_pred %*% V) * mm_pred))
pred_data$ci_low <- pred_data$fit - 1.96 * pred_data$se
pred_data$ci_high <- pred_data$fit + 1.96 * pred_data$se

# Centre at median DMA
median_dma <- median(disc_rcs$log_URXUDMA, na.rm = TRUE)
median_idx <- which.min(abs(pred_data$log_URXUDMA - median_dma))
pred_data$fit_centered <- pred_data$fit - pred_data$fit[median_idx]
pred_data$ci_low_centered <- pred_data$ci_low - pred_data$fit[median_idx]
pred_data$ci_high_centered <- pred_data$ci_high - pred_data$fit[median_idx]

# Figure 3: RCS plot
fig3_rcs <- ggplot(pred_data, aes(x = log_URXUDMA, y = fit_centered)) +
  geom_ribbon(aes(ymin = ci_low_centered, ymax = ci_high_centered),
              fill = "steelblue", alpha = 0.2) +
  geom_line(color = "steelblue", linewidth = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_rug(data = disc_rcs, aes(x = log_URXUDMA, y = NULL),
           sides = "b", alpha = 0.05) +
  labs(
    x = "Log urinary DMA (ug/L)",
    y = "Change in serum uric acid (mg/dL)",
    title = "Dose-response relationship: DMA and serum uric acid",
    subtitle = "Natural cubic splines (3 df), centred at median DMA"
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))

ggsave(file.path(OUT_DIR, "figures", "fig3_rcs.pdf"),
       fig3_rcs, width = 8, height = 6)
ggsave(file.path(OUT_DIR, "figures", "fig3_rcs.png"),
       fig3_rcs, width = 8, height = 6, dpi = 300)
cat("RCS figure saved.\n")


# ============================================================================
# ANALYSIS 4: E-value calculation
# ============================================================================

cat("\n=== ANALYSIS 4: E-values ===\n\n")

# E-value formula for continuous outcomes:
# E-value = RR + sqrt(RR * (RR - 1)) where RR is approximated from the
# standardized effect size
# For continuous outcomes, convert beta/SD(outcome) to approximate RR

sd_ua <- sd(disc_sub$LBXSUA, na.rm = TRUE)
beta_discovery <- 0.201
se_discovery <- 0.027
beta_ci_low <- beta_discovery - 1.96 * se_discovery

# Approximate RR using VanderWeele & Ding (2017) formula:
# RR_approx = exp(0.91 * beta / SD_outcome) for continuous outcomes
rr_point <- exp(0.91 * beta_discovery / sd_ua)
rr_ci_low <- exp(0.91 * beta_ci_low / sd_ua)

# E-value
evalue_point <- rr_point + sqrt(rr_point * (rr_point - 1))
evalue_ci <- rr_ci_low + sqrt(rr_ci_low * (rr_ci_low - 1))

cat(sprintf("SD(uric acid) = %.2f mg/dL\n", sd_ua))
cat(sprintf("Approximate RR (point): %.2f\n", rr_point))
cat(sprintf("Approximate RR (lower CI): %.2f\n", rr_ci_low))
cat(sprintf("E-value (point estimate): %.2f\n", evalue_point))
cat(sprintf("E-value (lower CI bound): %.2f\n", evalue_ci))

# Also for creatinine-adjusted estimate
beta_creat <- 0.135
se_creat <- 0.135 / qnorm(1 - 0.012/2)  # approximate from p=0.012
rr_creat <- exp(0.91 * beta_creat / sd_ua)
rr_creat_ci <- exp(0.91 * (beta_creat - 1.96 * se_creat) / sd_ua)
evalue_creat <- rr_creat + sqrt(rr_creat * (rr_creat - 1))
evalue_creat_ci <- rr_creat_ci + sqrt(rr_creat_ci * (rr_creat_ci - 1))

cat(sprintf("\nCreatinine-adjusted:\n"))
cat(sprintf("E-value (point): %.2f\n", evalue_creat))
cat(sprintf("E-value (lower CI): %.2f\n", evalue_creat_ci))

evalue_table <- tibble(
  Model = c("Primary (discovery)", "Primary (lower CI)",
            "Creatinine-adjusted", "Creatinine-adjusted (lower CI)"),
  Beta = c(beta_discovery, beta_ci_low, beta_creat, beta_creat - 1.96 * se_creat),
  Approx_RR = c(rr_point, rr_ci_low, rr_creat, rr_creat_ci),
  E_value = c(evalue_point, evalue_ci, evalue_creat, evalue_creat_ci)
)

write.csv(evalue_table,
          file.path(OUT_DIR, "tables", "table_s3_evalues.csv"),
          row.names = FALSE)
cat("E-value table saved.\n")


# ============================================================================
# ANALYSIS 5: eGFR adjustment (renal mediation)
# ============================================================================

cat("\n=== ANALYSIS 5: eGFR adjustment ===\n\n")

# Check if eGFR is available in discovery data
if ("eGFR" %in% names(dat_discovery)) {
  disc_egfr <- dma_subset(dat_discovery, extra_vars = "eGFR")
  des_egfr <- make_design(disc_egfr)

  # Model with eGFR
  egfr_model <- svyglm(
    LBXSUA ~ log_URXUDMA + RIDAGEYR + female + race3 + INDFMPIR + BMXBMI + smoker + eGFR,
    design = des_egfr
  )
  egfr_tidy <- tidy(egfr_model, conf.int = TRUE)
  egfr_row <- egfr_tidy %>% filter(term == "log_URXUDMA")

  # Primary model in same subset (for fair comparison)
  primary_egfr <- svyglm(
    LBXSUA ~ log_URXUDMA + RIDAGEYR + female + race3 + INDFMPIR + BMXBMI + smoker,
    design = des_egfr
  )
  primary_tidy <- tidy(primary_egfr, conf.int = TRUE)
  primary_row <- primary_tidy %>% filter(term == "log_URXUDMA")

  pct_change <- (egfr_row$estimate - primary_row$estimate) / abs(primary_row$estimate) * 100

  cat(sprintf("Primary (in eGFR subset, n=%d): beta = %.3f (p = %.4f)\n",
              nrow(disc_egfr), primary_row$estimate, primary_row$p.value))
  cat(sprintf("eGFR-adjusted: beta = %.3f (p = %.4f)\n",
              egfr_row$estimate, egfr_row$p.value))
  cat(sprintf("Change: %.1f%%\n", pct_change))

  egfr_result <- tibble(
    Model = c("Primary (eGFR subset)", "eGFR-adjusted"),
    Beta = c(primary_row$estimate, egfr_row$estimate),
    SE = c(primary_row$std.error, egfr_row$std.error),
    CI_low = c(primary_row$conf.low, egfr_row$conf.low),
    CI_high = c(primary_row$conf.high, egfr_row$conf.high),
    P_value = c(primary_row$p.value, egfr_row$p.value),
    N = nrow(disc_egfr),
    Pct_change = c(NA, pct_change)
  )

  write.csv(egfr_result,
            file.path(OUT_DIR, "tables", "egfr_adjustment.csv"),
            row.names = FALSE)
} else {
  cat("eGFR not in discovery data; computing from serum creatinine.\n")

  # Try to compute eGFR if serum creatinine is available
  if ("LBXSCR" %in% names(dat_discovery)) {
    dat_discovery$eGFR <- calc_eGFR(
      creatinine = dat_discovery$LBXSCR,
      age = dat_discovery$RIDAGEYR,
      female = dat_discovery$female
    )
    cat("Computed eGFR from serum creatinine. Re-running eGFR analysis...\n")

    disc_egfr <- dma_subset(dat_discovery, extra_vars = "eGFR")
    des_egfr <- make_design(disc_egfr)

    egfr_model <- svyglm(
      LBXSUA ~ log_URXUDMA + RIDAGEYR + female + race3 + INDFMPIR + BMXBMI + smoker + eGFR,
      design = des_egfr
    )
    egfr_tidy <- tidy(egfr_model, conf.int = TRUE)
    egfr_row <- egfr_tidy %>% filter(term == "log_URXUDMA")

    primary_egfr <- svyglm(
      LBXSUA ~ log_URXUDMA + RIDAGEYR + female + race3 + INDFMPIR + BMXBMI + smoker,
      design = des_egfr
    )
    primary_tidy <- tidy(primary_egfr, conf.int = TRUE)
    primary_row <- primary_tidy %>% filter(term == "log_URXUDMA")

    pct_change <- (egfr_row$estimate - primary_row$estimate) / abs(primary_row$estimate) * 100

    cat(sprintf("Primary (in eGFR subset, n=%d): beta = %.3f (p = %.4f)\n",
                nrow(disc_egfr), primary_row$estimate, primary_row$p.value))
    cat(sprintf("eGFR-adjusted: beta = %.3f (p = %.4f)\n",
                egfr_row$estimate, egfr_row$p.value))
    cat(sprintf("Change: %.1f%%\n", pct_change))

    egfr_result <- tibble(
      Model = c("Primary (eGFR subset)", "eGFR-adjusted"),
      Beta = c(primary_row$estimate, egfr_row$estimate),
      SE = c(primary_row$std.error, egfr_row$std.error),
      CI_low = c(primary_row$conf.low, egfr_row$conf.low),
      CI_high = c(primary_row$conf.high, egfr_row$conf.high),
      P_value = c(primary_row$p.value, egfr_row$p.value),
      N = nrow(disc_egfr),
      Pct_change = c(NA, pct_change)
    )

    write.csv(egfr_result,
              file.path(OUT_DIR, "tables", "egfr_adjustment.csv"),
              row.names = FALSE)
  } else {
    cat("LBXSCR not available; downloading serum biochemistry...\n")
    biopro <- tryCatch(nhanes("BIOPRO_J"), error = function(e) NULL)
    if (!is.null(biopro) && "LBXSCR" %in% names(biopro)) {
      biopro <- biopro %>% transmute(SEQN, LBXSCR = as.numeric(as.character(LBXSCR)))
      dat_discovery <- dat_discovery %>% left_join(biopro, by = "SEQN")
      dat_discovery$eGFR <- calc_eGFR(dat_discovery$LBXSCR, dat_discovery$RIDAGEYR, dat_discovery$female)

      disc_egfr <- dma_subset(dat_discovery, extra_vars = "eGFR")
      des_egfr <- make_design(disc_egfr)

      egfr_model <- svyglm(
        LBXSUA ~ log_URXUDMA + RIDAGEYR + female + race3 + INDFMPIR + BMXBMI + smoker + eGFR,
        design = des_egfr
      )
      egfr_tidy <- tidy(egfr_model, conf.int = TRUE)
      egfr_row <- egfr_tidy %>% filter(term == "log_URXUDMA")

      primary_egfr <- svyglm(
        LBXSUA ~ log_URXUDMA + RIDAGEYR + female + race3 + INDFMPIR + BMXBMI + smoker,
        design = des_egfr
      )
      primary_tidy <- tidy(primary_egfr, conf.int = TRUE)
      primary_row <- primary_tidy %>% filter(term == "log_URXUDMA")

      pct_change <- (egfr_row$estimate - primary_row$estimate) / abs(primary_row$estimate) * 100

      cat(sprintf("Primary (in eGFR subset, n=%d): beta = %.3f (p = %.4f)\n",
                  nrow(disc_egfr), primary_row$estimate, primary_row$p.value))
      cat(sprintf("eGFR-adjusted: beta = %.3f (p = %.4f)\n",
                  egfr_row$estimate, egfr_row$p.value))
      cat(sprintf("Change: %.1f%%\n", pct_change))

      egfr_result <- tibble(
        Model = c("Primary (eGFR subset)", "eGFR-adjusted"),
        Beta = c(primary_row$estimate, egfr_row$estimate),
        SE = c(primary_row$std.error, egfr_row$std.error),
        CI_low = c(primary_row$conf.low, egfr_row$conf.low),
        CI_high = c(primary_row$conf.high, egfr_row$conf.high),
        P_value = c(primary_row$p.value, egfr_row$p.value),
        N = nrow(disc_egfr),
        Pct_change = c(NA, pct_change)
      )

      write.csv(egfr_result,
                file.path(OUT_DIR, "tables", "egfr_adjustment.csv"),
                row.names = FALSE)
    } else {
      cat("WARNING: Could not obtain serum creatinine for eGFR calculation.\n")
    }
  }
}


# ============================================================================
# ANALYSIS 6: Total arsenic correlation
# ============================================================================

cat("\n=== ANALYSIS 6: Total arsenic--DMA correlation ===\n\n")

if ("log_URXUAS" %in% names(disc_sub) && "log_URXUDMA" %in% names(disc_sub)) {
  arsenic_sub <- disc_sub %>%
    filter(!is.na(log_URXUAS), !is.na(log_URXUDMA))

  r_pearson <- cor(arsenic_sub$log_URXUAS, arsenic_sub$log_URXUDMA)
  r_spearman <- cor(arsenic_sub$log_URXUAS, arsenic_sub$log_URXUDMA, method = "spearman")

  cat(sprintf("n = %d participants with both total As and DMA\n", nrow(arsenic_sub)))
  cat(sprintf("Pearson r (log scale): %.3f\n", r_pearson))
  cat(sprintf("Spearman rho: %.3f\n", r_spearman))

  # Summary statistics
  if ("URXUDMA" %in% names(arsenic_sub) && "URXUAS" %in% names(arsenic_sub)) {
    cat(sprintf("\nDMA: median = %.1f ug/L, IQR = %.1f-%.1f\n",
                median(arsenic_sub$URXUDMA, na.rm = TRUE),
                quantile(arsenic_sub$URXUDMA, 0.25, na.rm = TRUE),
                quantile(arsenic_sub$URXUDMA, 0.75, na.rm = TRUE)))
    cat(sprintf("Total As: median = %.1f ug/L, IQR = %.1f-%.1f\n",
                median(arsenic_sub$URXUAS, na.rm = TRUE),
                quantile(arsenic_sub$URXUAS, 0.25, na.rm = TRUE),
                quantile(arsenic_sub$URXUAS, 0.75, na.rm = TRUE)))

    # DMA as fraction of total arsenic
    arsenic_sub$dma_fraction <- arsenic_sub$URXUDMA / arsenic_sub$URXUAS
    cat(sprintf("DMA as %% of total As: median = %.1f%%, IQR = %.1f%%-%.1f%%\n",
                100 * median(arsenic_sub$dma_fraction, na.rm = TRUE),
                100 * quantile(arsenic_sub$dma_fraction, 0.25, na.rm = TRUE),
                100 * quantile(arsenic_sub$dma_fraction, 0.75, na.rm = TRUE)))
  }

  # Scatter plot
  fig_s1_scatter <- ggplot(arsenic_sub, aes(x = log_URXUAS, y = log_URXUDMA)) +
    geom_point(alpha = 0.3, size = 1) +
    geom_smooth(method = "lm", color = "steelblue", se = TRUE) +
    annotate("text", x = Inf, y = -Inf,
             label = sprintf("r = %.2f\nrho = %.2f\nn = %d",
                             r_pearson, r_spearman, nrow(arsenic_sub)),
             hjust = 1.1, vjust = -0.5, size = 4) +
    labs(
      x = "Log total urinary arsenic (ug/L)",
      y = "Log urinary DMA (ug/L)",
      title = "Correlation between total arsenic and DMA"
    ) +
    theme_minimal(base_size = 12)

  ggsave(file.path(OUT_DIR, "figures", "fig_s1_arsenic_correlation.pdf"),
         fig_s1_scatter, width = 7, height = 6)
  ggsave(file.path(OUT_DIR, "figures", "fig_s1_arsenic_correlation.png"),
         fig_s1_scatter, width = 7, height = 6, dpi = 300)
  cat("Arsenic correlation figure saved.\n")

} else {
  cat("Total arsenic (URXUAS) not available in dataset.\n")
}


# ============================================================================
# ANALYSIS 7: Validation cycle dose-response
# ============================================================================

cat("\n=== ANALYSIS 7: Validation cycle dose-response ===\n\n")

val_dr <- dma_subset(val_data)
val_dr$dma_quartile <- ntile(val_dr$log_URXUDMA, 4)
val_dr$dma_q_factor <- factor(val_dr$dma_quartile)
val_dr$dma_q_numeric <- val_dr$dma_quartile

des_val_dr <- make_design(val_dr)

# Quartile model
val_q_model <- svyglm(
  LBXSUA ~ dma_q_factor + RIDAGEYR + female + race3 + INDFMPIR + BMXBMI + smoker,
  design = des_val_dr
)
val_q_ct <- summary(val_q_model)$coefficients
q_rows <- grep("^dma_q_factor", rownames(val_q_ct))

val_q_estimates <- data.frame(
  quartile = c(1, 2, 3, 4),
  beta = c(0, val_q_ct[q_rows, "Estimate"]),
  se = c(0, val_q_ct[q_rows, "Std. Error"]),
  p_value = c(NA, val_q_ct[q_rows, "Pr(>|t|)"]),
  stringsAsFactors = FALSE
) %>%
  mutate(ci_low = beta - 1.96 * se, ci_high = beta + 1.96 * se)

# Trend test
val_trend <- svyglm(
  LBXSUA ~ dma_q_numeric + RIDAGEYR + female + race3 + INDFMPIR + BMXBMI + smoker,
  design = des_val_dr
)
val_trend_p <- summary(val_trend)$coefficients["dma_q_numeric", "Pr(>|t|)"]
val_monotonic <- all(diff(val_q_estimates$beta) >= 0) | all(diff(val_q_estimates$beta) <= 0)

cat("Validation cycle (2015-2016) dose-response:\n")
cat(sprintf("  Q2: beta = %.2f\n  Q3: beta = %.2f\n  Q4: beta = %.2f\n",
            val_q_estimates$beta[2], val_q_estimates$beta[3], val_q_estimates$beta[4]))
cat(sprintf("  p_trend = %.4f, monotonic = %s\n", val_trend_p, val_monotonic))

# Extract discovery dose-response from P08 results
disc_dr <- NULL
for (r in dr_results) {
  if (r$exposure == "log_URXUDMA" && r$outcome == "LBXSUA") {
    disc_dr <- r
    break
  }
}

# Figure 1: Combined dose-response (both cycles)
if (!is.null(disc_dr)) {
  disc_q_df <- disc_dr$q_estimates %>%
    mutate(cycle = "2017-2018 (Discovery)")

  val_q_df <- val_q_estimates %>%
    mutate(cycle = "2015-2016 (Validation)")

  combined_dr <- bind_rows(disc_q_df, val_q_df)

  fig1_dr <- ggplot(combined_dr, aes(x = quartile, y = beta, color = cycle)) +
    geom_pointrange(aes(ymin = ci_low, ymax = ci_high),
                    position = position_dodge(width = 0.3), size = 0.8) +
    # No connecting lines: quartiles are ordinal, not continuous
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    scale_x_continuous(breaks = 1:4, labels = paste0("Q", 1:4)) +
    scale_color_manual(values = c("2017-2018 (Discovery)" = "steelblue",
                                   "2015-2016 (Validation)" = "darkorange")) +
    annotate("text", x = 3.5, y = max(combined_dr$ci_high, na.rm = TRUE) * 0.9,
             label = sprintf("Discovery p-trend = %.4f\nValidation p-trend = %.4f",
                             disc_dr$p_trend, val_trend_p),
             size = 3.5, hjust = 0.5) +
    labs(
      x = "DMA quartile",
      y = "Difference in serum uric acid vs Q1 (mg/dL)",
      title = "Dose-response: urinary DMA and serum uric acid",
      subtitle = "Quartile analysis in two independent NHANES cycles",
      color = "NHANES cycle"
    ) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "bottom",
          plot.title = element_text(face = "bold"))

  ggsave(file.path(OUT_DIR, "figures", "fig1_dose_response.pdf"),
         fig1_dr, width = 8, height = 6)
  ggsave(file.path(OUT_DIR, "figures", "fig1_dose_response.png"),
         fig1_dr, width = 8, height = 6, dpi = 300)
  cat("Figure 1 (combined dose-response) saved.\n")
}

# Save validation DR results
write.csv(val_q_estimates,
          file.path(OUT_DIR, "tables", "validation_dose_response.csv"),
          row.names = FALSE)


# ============================================================================
# ANALYSIS 8: Gout medication sensitivity
# ============================================================================

cat("\n=== ANALYSIS 8: Gout medication sensitivity ===\n\n")

# Try to download prescription data
cat("Downloading RXQ_RX_J (prescription medication data)...\n")
rxq <- tryCatch(nhanes("RXQ_RX_J"), error = function(e) {
  cat("  Note: Could not download prescription data:", e$message, "\n")
  NULL
})

if (!is.null(rxq)) {
  # Look for urate-lowering therapy
  # RXDDRUG contains drug names; look for allopurinol, febuxostat, probenecid
  rxq_char <- rxq %>%
    mutate(drug_name = toupper(as.character(RXDDRUG)))

  ult_users <- rxq_char %>%
    filter(grepl("ALLOPURINOL|FEBUXOSTAT|PROBENECID|COLCHICINE", drug_name)) %>%
    pull(SEQN) %>%
    unique()

  cat(sprintf("  Found %d participants on urate-lowering/gout therapy\n", length(ult_users)))

  # Exclude ULT users and re-run
  disc_no_ult <- disc_sub %>% filter(!(SEQN %in% ult_users))
  cat(sprintf("  Sample after exclusion: n = %d (excluded %d)\n",
              nrow(disc_no_ult), nrow(disc_sub) - nrow(disc_no_ult)))

  des_no_ult <- make_design(disc_no_ult)

  ult_model <- svyglm(
    LBXSUA ~ log_URXUDMA + RIDAGEYR + female + race3 + INDFMPIR + BMXBMI + smoker,
    design = des_no_ult
  )
  ult_tidy <- tidy(ult_model, conf.int = TRUE)
  ult_row <- ult_tidy %>% filter(term == "log_URXUDMA")

  pct_change <- (ult_row$estimate - 0.201) / 0.201 * 100
  cat(sprintf("  Beta (excl ULT): %.3f (95%% CI: %.3f-%.3f), p = %.4f\n",
              ult_row$estimate, ult_row$conf.low, ult_row$conf.high, ult_row$p.value))
  cat(sprintf("  Change from primary: %.1f%%\n", pct_change))

  ult_result <- tibble(
    Model = c("Primary", "Excluding ULT users"),
    Beta = c(0.201, ult_row$estimate),
    SE = c(0.027, ult_row$std.error),
    CI_low = c(0.201 - 1.96*0.027, ult_row$conf.low),
    CI_high = c(0.201 + 1.96*0.027, ult_row$conf.high),
    P_value = c(1.5e-04, ult_row$p.value),
    N = c(nrow(disc_sub), nrow(disc_no_ult)),
    ULT_excluded = c(0, nrow(disc_sub) - nrow(disc_no_ult))
  )

  write.csv(ult_result,
            file.path(OUT_DIR, "tables", "gout_medication_sensitivity.csv"),
            row.names = FALSE)

} else {
  cat("  Skipping gout medication sensitivity (data unavailable).\n")
}


# ============================================================================
# ANALYSIS 9: 6-level race/ethnicity
# ============================================================================

cat("\n=== ANALYSIS 9: 6-level race/ethnicity ===\n\n")

if ("RIDRETH3" %in% names(dat_discovery)) {
  # Create 6-level race factor
  ridreth3_char <- as.character(dat_discovery$RIDRETH3)

  dat_discovery$race6 <- case_when(
    ridreth3_char == "Mexican American" ~ "Mexican American",
    ridreth3_char == "Other Hispanic" ~ "Other Hispanic",
    ridreth3_char == "Non-Hispanic White" ~ "Non-Hispanic White",
    ridreth3_char == "Non-Hispanic Black" ~ "Non-Hispanic Black",
    ridreth3_char == "Non-Hispanic Asian" ~ "Non-Hispanic Asian",
    ridreth3_char == "Other Race - Including Multi-Racial" ~ "Other/Multiracial",
    ridreth3_char == "1" ~ "Mexican American",
    ridreth3_char == "2" ~ "Other Hispanic",
    ridreth3_char == "3" ~ "Non-Hispanic White",
    ridreth3_char == "4" ~ "Non-Hispanic Black",
    ridreth3_char == "6" ~ "Non-Hispanic Asian",
    ridreth3_char == "7" ~ "Other/Multiracial",
    TRUE ~ NA_character_
  )
  dat_discovery$race6 <- factor(dat_discovery$race6,
    levels = c("Non-Hispanic White", "Non-Hispanic Black",
               "Mexican American", "Other Hispanic",
               "Non-Hispanic Asian", "Other/Multiracial"))

  disc_race6 <- dma_subset(dat_discovery, extra_vars = "race6")

  cat("6-level race distribution:\n")
  print(table(disc_race6$race6, useNA = "ifany"))

  des_race6 <- make_design(disc_race6)

  # 3-level model (primary, for comparison in same subset)
  model_3level <- svyglm(
    LBXSUA ~ log_URXUDMA + RIDAGEYR + female + race3 + INDFMPIR + BMXBMI + smoker,
    design = des_race6
  )
  tidy_3 <- tidy(model_3level, conf.int = TRUE) %>% filter(term == "log_URXUDMA")

  # 6-level model
  model_6level <- svyglm(
    LBXSUA ~ log_URXUDMA + RIDAGEYR + female + race6 + INDFMPIR + BMXBMI + smoker,
    design = des_race6
  )
  tidy_6 <- tidy(model_6level, conf.int = TRUE) %>% filter(term == "log_URXUDMA")

  pct_change <- (tidy_6$estimate - tidy_3$estimate) / abs(tidy_3$estimate) * 100

  cat(sprintf("\n3-level race: beta = %.3f (p = %.4f)\n", tidy_3$estimate, tidy_3$p.value))
  cat(sprintf("6-level race: beta = %.3f (p = %.4f)\n", tidy_6$estimate, tidy_6$p.value))
  cat(sprintf("Change: %.1f%%\n", pct_change))

  race6_result <- tibble(
    Model = c("3-level race", "6-level race"),
    Beta = c(tidy_3$estimate, tidy_6$estimate),
    SE = c(tidy_3$std.error, tidy_6$std.error),
    CI_low = c(tidy_3$conf.low, tidy_6$conf.low),
    CI_high = c(tidy_3$conf.high, tidy_6$conf.high),
    P_value = c(tidy_3$p.value, tidy_6$p.value),
    N = nrow(disc_race6),
    Pct_change = c(NA, pct_change)
  )

  write.csv(race6_result,
            file.path(OUT_DIR, "tables", "race6_sensitivity.csv"),
            row.names = FALSE)
} else {
  cat("RIDRETH3 not available in dataset.\n")
}


# ============================================================================
# ANALYSIS 10: Combined sensitivity forest plot (Figure 2)
# ============================================================================

cat("\n=== ANALYSIS 10: Combined sensitivity forest plot ===\n\n")

# Collect all results into a single data frame for the forest plot

# Start with results we know from P08
forest_data <- tibble(
  label = character(),
  beta = double(),
  ci_low = double(),
  ci_high = double(),
  p_value = double(),
  n = integer(),
  category = character()
)

# 1. Primary model (discovery)
forest_data <- bind_rows(forest_data, tibble(
  label = "Primary model (2017-2018)",
  beta = 0.201, ci_low = 0.201 - 1.96*0.027, ci_high = 0.201 + 1.96*0.027,
  p_value = 1.5e-04, n = 1593L, category = "Primary"
))

# 2. Cross-cycle validation
forest_data <- bind_rows(forest_data, tibble(
  label = "Validation (2015-2016)",
  beta = 0.206, ci_low = 0.206 - 1.96*0.068, ci_high = 0.206 + 1.96*0.068,
  p_value = 0.003, n = 1575L, category = "Primary"
))

# 3. Extract sensitivity results from P08 sens_df
dma_sens <- sens_df %>%
  filter(exposure == "log_URXUDMA", outcome == "LBXSUA", !is.na(beta))

for (i in seq_len(nrow(dma_sens))) {
  if (dma_sens$analysis[i] == "Primary model") next  # already included
  forest_data <- bind_rows(forest_data, tibble(
    label = dma_sens$analysis[i],
    beta = dma_sens$beta[i],
    ci_low = dma_sens$beta[i] - 1.96 * dma_sens$se[i],
    ci_high = dma_sens$beta[i] + 1.96 * dma_sens$se[i],
    p_value = dma_sens$p_value[i],
    n = as.integer(dma_sens$n[i]),
    category = "Stratification/Subgroup"
  ))
}

# 4. Extract additional sensitivity results from P08 add_sens_df
dma_add_sens <- add_sens_df %>%
  filter(exposure == "log_URXUDMA", outcome == "LBXSUA", !is.na(beta))

for (i in seq_len(nrow(dma_add_sens))) {
  forest_data <- bind_rows(forest_data, tibble(
    label = dma_add_sens$analysis[i],
    beta = dma_add_sens$beta[i],
    ci_low = dma_add_sens$beta[i] - 1.96 * dma_add_sens$se[i],
    ci_high = dma_add_sens$beta[i] + 1.96 * dma_add_sens$se[i],
    p_value = dma_add_sens$p_value[i],
    n = as.integer(dma_add_sens$n[i]),
    category = "Confounder adjustment"
  ))
}

# 5. Protein adjustment from P08
if (exists("protein_sens_df")) {
  prot_dma <- protein_sens_df %>%
    filter(exposure == "log_URXUDMA", outcome == "LBXSUA", !is.na(beta))
  for (i in seq_len(nrow(prot_dma))) {
    forest_data <- bind_rows(forest_data, tibble(
      label = prot_dma$analysis[i],
      beta = prot_dma$beta[i],
      ci_low = prot_dma$beta[i] - 1.96 * prot_dma$se[i],
      ci_high = prot_dma$beta[i] + 1.96 * prot_dma$se[i],
      p_value = prot_dma$p_value[i],
      n = as.integer(prot_dma$n[i]),
      category = "Confounder adjustment"
    ))
  }
}

# 6. Add new results from this script (if they exist)
if (exists("egfr_result")) {
  egfr_adj <- egfr_result %>% filter(Model == "eGFR-adjusted")
  forest_data <- bind_rows(forest_data, tibble(
    label = "Adjust for eGFR",
    beta = egfr_adj$Beta,
    ci_low = egfr_adj$CI_low,
    ci_high = egfr_adj$CI_high,
    p_value = egfr_adj$P_value,
    n = as.integer(egfr_adj$N),
    category = "Confounder adjustment"
  ))
}

if (exists("ult_result")) {
  ult_adj <- ult_result %>% filter(Model == "Excluding ULT users")
  forest_data <- bind_rows(forest_data, tibble(
    label = "Exclude gout medication",
    beta = ult_adj$Beta,
    ci_low = ult_adj$CI_low,
    ci_high = ult_adj$CI_high,
    p_value = ult_adj$P_value,
    n = as.integer(ult_adj$N),
    category = "Sample restriction"
  ))
}

if (exists("race6_result")) {
  r6 <- race6_result %>% filter(Model == "6-level race")
  forest_data <- bind_rows(forest_data, tibble(
    label = "6-level race/ethnicity",
    beta = r6$Beta,
    ci_low = r6$CI_low,
    ci_high = r6$CI_high,
    p_value = r6$P_value,
    n = as.integer(r6$N),
    category = "Confounder adjustment"
  ))
}

# Remove duplicates
forest_data <- forest_data %>% distinct(label, .keep_all = TRUE)

# Remove protein density (redundant with protein intake)
forest_data <- forest_data %>% filter(label != "Adjust for protein density")

# Order for display
category_order <- c("Primary", "Stratification/Subgroup",
                     "Confounder adjustment", "Sample restriction")
forest_data$category <- factor(forest_data$category, levels = category_order)
forest_data <- forest_data %>% arrange(category, desc(abs(beta)))

# Assign row positions
forest_data$row <- rev(seq_len(nrow(forest_data)))

# Figure 2: Forest plot
fig2_forest <- ggplot(forest_data, aes(x = beta, y = row)) +
  geom_vline(xintercept = 0, linetype = "solid", color = "grey40", linewidth = 0.5) +
  geom_vline(xintercept = 0.201, linetype = "dotted", color = "steelblue", alpha = 0.5) +
  geom_pointrange(aes(xmin = ci_low, xmax = ci_high, color = category),
                  size = 0.5) +
  scale_y_continuous(
    breaks = forest_data$row,
    labels = paste0(forest_data$label, " (n=", format(forest_data$n, big.mark = ","), ")")
  ) +
  scale_color_manual(values = c(
    "Primary" = "steelblue",
    "Stratification/Subgroup" = "darkorange",
    "Confounder adjustment" = "forestgreen",
    "Sample restriction" = "purple"
  )) +
  labs(
    x = expression("Effect on serum uric acid (mg/dL per log-unit DMA, " * beta * " [95% CI])"),
    y = NULL,
    title = "Sensitivity analysis: DMA and serum uric acid",
    subtitle = "All model specifications",
    color = "Analysis type"
  ) +
  coord_cartesian(xlim = c(-0.05, NA)) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold"),
    axis.text.y = element_text(size = 9)
  )

ggsave(file.path(OUT_DIR, "figures", "fig2_forest.pdf"),
       fig2_forest, width = 10, height = max(6, nrow(forest_data) * 0.4))
ggsave(file.path(OUT_DIR, "figures", "fig2_forest.png"),
       fig2_forest, width = 10, height = max(6, nrow(forest_data) * 0.4), dpi = 300)
cat("Figure 2 (forest plot) saved.\n")

# Save forest data as Table 2
forest_table <- forest_data %>%
  arrange(category, desc(abs(beta))) %>%
  transmute(
    Category = as.character(category),
    Model = label,
    Beta = round(beta, 3),
    `CI_low` = round(ci_low, 3),
    `CI_high` = round(ci_high, 3),
    P_value = signif(p_value, 3),
    N = n
  )

write.csv(forest_table,
          file.path(OUT_DIR, "tables", "table2_sensitivity_summary.csv"),
          row.names = FALSE)


# ============================================================================
# SUPPLEMENTARY: Distribution plots
# ============================================================================

cat("\n=== Supplementary: Distribution plots ===\n\n")

# Figure S2: Distribution of DMA and uric acid in both cycles
disc_plot <- disc_sub %>%
  select(log_URXUDMA, LBXSUA) %>%
  mutate(cycle = "2017-2018 (Discovery)")

val_plot <- val_sub %>%
  select(log_URXUDMA, LBXSUA) %>%
  mutate(cycle = "2015-2016 (Validation)")

combined_plot <- bind_rows(disc_plot, val_plot)

fig_s2a <- ggplot(combined_plot, aes(x = log_URXUDMA, fill = cycle)) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = c("2017-2018 (Discovery)" = "steelblue",
                                "2015-2016 (Validation)" = "darkorange")) +
  labs(x = "Log urinary DMA (ug/L)", y = "Density",
       title = "Distribution of urinary DMA", fill = "Cycle") +
  theme_minimal(base_size = 11)

fig_s2b <- ggplot(combined_plot, aes(x = LBXSUA, fill = cycle)) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = c("2017-2018 (Discovery)" = "steelblue",
                                "2015-2016 (Validation)" = "darkorange")) +
  labs(x = "Serum uric acid (mg/dL)", y = "Density",
       title = "Distribution of serum uric acid", fill = "Cycle") +
  theme_minimal(base_size = 11)

# Combine into single figure
if (!requireNamespace("patchwork", quietly = TRUE)) {
  install.packages("patchwork", repos = "https://cloud.r-project.org")
}
library(patchwork)

fig_s2 <- fig_s2a / fig_s2b + plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

ggsave(file.path(OUT_DIR, "figures", "fig_s2_distributions.pdf"),
       fig_s2, width = 8, height = 8)
ggsave(file.path(OUT_DIR, "figures", "fig_s2_distributions.png"),
       fig_s2, width = 8, height = 8, dpi = 300)
cat("Distribution figure saved.\n")


# ============================================================================
# SAVE SESSION
# ============================================================================

cat("\n=== Saving session ===\n\n")

# Save key results for manuscript drafting
save(
  table1, forest_data, forest_table,
  hu_results, hu_q_results,
  val_q_estimates, val_trend_p, val_monotonic,
  evalue_table,
  file = file.path(OUT_DIR, "analysis_results.RData")
)

cat("All analyses complete.\n")
cat("Outputs saved to:\n")
cat(sprintf("  %s/figures/\n", OUT_DIR))
cat(sprintf("  %s/tables/\n", OUT_DIR))
cat(sprintf("  %s/analysis_results.RData\n", OUT_DIR))

# List all output files
cat("\nGenerated files:\n")
for (d in c("figures", "tables")) {
  fpath <- file.path(OUT_DIR, "outputs", d)
  files <- list.files(fpath, full.names = FALSE)
  for (f in files) cat(sprintf("  outputs/%s/%s\n", d, f))
}

cat("\n=== 01_focused_analysis.R complete ===\n")
