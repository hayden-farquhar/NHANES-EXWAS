# ============================================================================
# 04_dose_response.R
# Quartile-based dose-response analysis for validated findings
# ============================================================================

source("00_functions.R")

# --- 1. Load data ------------------------------------------------------------

load("validation_results.RData")  # validated

# Load primary analysis data (2017-2018)
# .RDataTmp is corrupted; use novelty screen data (most complete) or expanded results
if (file.exists("exwas_novelty_screen.RData")) {
  load("exwas_novelty_screen.RData")  # gives `dat` (expanded + novelty chemicals)
  dat_full <- dat
  rm(dat)
  cat(sprintf("Loaded novelty screen data: %d obs x %d cols\n",
              nrow(dat_full), ncol(dat_full)))
} else if (file.exists("exwas_expanded_results.RData")) {
  load("exwas_expanded_results.RData")  # gives `expanded`
  dat_full <- expanded
  rm(expanded)
  cat(sprintf("Loaded expanded results data: %d obs x %d cols\n",
              nrow(dat_full), ncol(dat_full)))
} else {
  stop("Cannot find primary analysis data. Run 01b_recover_novelty_screen.R first.")
}

# --- 2. Select findings for dose-response ------------------------------------

# Use validated findings (replicated in 2015-2016) if available,
# otherwise fall back to all FDR-significant primary findings
dr_findings <- validated %>%
  filter(validated | suggestive_val) %>%
  select(exposure, outcome, exp_label, out_label, type, beta_primary, dir_primary)

if (nrow(dr_findings) == 0) {
  cat("No validated findings. Using all FDR-significant primary findings.\n")
  dr_findings <- validated %>%
    filter(status == "ok" | status == "exposure_missing") %>%
    head(20) %>%
    select(exposure, outcome, exp_label, out_label, type, beta_primary, dir_primary)
}

cat(sprintf("Dose-response analysis for %d findings\n\n", nrow(dr_findings)))

# --- 3. Run quartile analysis ------------------------------------------------

dr_results <- list()

for (i in 1:nrow(dr_findings)) {
  exp_v  <- dr_findings$exposure[i]
  out_v  <- dr_findings$outcome[i]
  binary <- dr_findings$type[i] == "binary"

  if (!(exp_v %in% names(dat_full)) || !(out_v %in% names(dat_full))) {
    cat(sprintf("  SKIP %s -> %s: variable missing\n", exp_v, out_v))
    next
  }

  wt <- select_weight(exp_v, dat_full)

  # Determine covariates
  covars <- if (out_v %in% c("BMXBMI", "BMXWAIST")) COVARS_NO_BMI else COVARS_FULL

  # Identify complete cases
  all_vars_needed <- c(out_v, exp_v, all.vars(as.formula(paste("~", covars))),
                       "SDMVPSU", "SDMVSTRA", wt)
  all_vars_needed <- all_vars_needed[all_vars_needed %in% names(dat_full)]
  sub_dat <- dat_full[complete.cases(dat_full[, all_vars_needed]), ]

  if (nrow(sub_dat) < 100) {
    cat(sprintf("  SKIP %s -> %s: n=%d too small\n", exp_v, out_v, nrow(sub_dat)))
    next
  }

  # Create quartiles of exposure (both factor and numeric for trend test)
  sub_dat$exp_quartile <- ntile(sub_dat[[exp_v]], 4)
  sub_dat$exp_q_factor <- factor(sub_dat$exp_quartile)
  sub_dat$exp_q_numeric <- sub_dat$exp_quartile

  # Compute quartile medians for plotting
  q_medians <- sub_dat %>%
    group_by(exp_quartile) %>%
    summarise(
      exp_median = median(.data[[exp_v]], na.rm = TRUE),
      n = n(),
      .groups = "drop"
    )

  # a) Quartile model (Q1 = reference)
  formula_q <- as.formula(paste(out_v, "~ exp_q_factor +", covars))

  des <- svydesign(
    id = ~SDMVPSU, strata = ~SDMVSTRA,
    weights = as.formula(paste0("~", wt)),
    nest = TRUE, data = sub_dat
  )

  if (binary) {
    mod_q <- svyglm(formula_q, design = des, family = quasibinomial())
  } else {
    mod_q <- svyglm(formula_q, design = des)
  }

  ct_q <- summary(mod_q)$coefficients

  # Extract Q2, Q3, Q4 vs Q1 (reference)
  q_rows <- grep("^exp_q_factor", rownames(ct_q))
  q_estimates <- data.frame(
    quartile = c(1, 2, 3, 4),
    beta     = c(0, ct_q[q_rows, "Estimate"]),
    se       = c(0, ct_q[q_rows, "Std. Error"]),
    p_value  = c(NA, ct_q[q_rows, "Pr(>|t|)"]),
    stringsAsFactors = FALSE
  ) %>%
    left_join(q_medians, by = c("quartile" = "exp_quartile")) %>%
    mutate(
      ci_lower = beta - 1.96 * se,
      ci_upper = beta + 1.96 * se
    )

  # b) Linear trend test (quartile as ordinal)
  formula_trend <- as.formula(paste(out_v, "~ exp_q_numeric +", covars))
  if (binary) {
    mod_trend <- svyglm(formula_trend, design = des, family = quasibinomial())
  } else {
    mod_trend <- svyglm(formula_trend, design = des)
  }
  ct_trend <- summary(mod_trend)$coefficients
  p_trend <- ct_trend["exp_q_numeric", "Pr(>|t|)"]
  beta_trend <- ct_trend["exp_q_numeric", "Estimate"]

  cat(sprintf("  %s -> %s: p_trend=%.4f, monotonic=%s\n",
              exp_v, out_v, p_trend,
              ifelse(all(diff(q_estimates$beta) >= 0) | all(diff(q_estimates$beta) <= 0),
                     "YES", "NO")))

  dr_results[[length(dr_results) + 1]] <- list(
    exposure     = exp_v,
    outcome      = out_v,
    exp_label    = dr_findings$exp_label[i],
    out_label    = dr_findings$out_label[i],
    q_estimates  = q_estimates,
    p_trend      = p_trend,
    beta_trend   = beta_trend,
    monotonic    = all(diff(q_estimates$beta) >= 0) | all(diff(q_estimates$beta) <= 0),
    n_total      = nrow(sub_dat)
  )
}

# --- 4. Summary table --------------------------------------------------------

cat("\n=== DOSE-RESPONSE SUMMARY ===\n\n")

dr_summary <- map_dfr(dr_results, function(r) {
  data.frame(
    exposure  = r$exp_label,
    outcome   = r$out_label,
    Q2_beta   = r$q_estimates$beta[2],
    Q3_beta   = r$q_estimates$beta[3],
    Q4_beta   = r$q_estimates$beta[4],
    p_trend   = r$p_trend,
    monotonic = r$monotonic,
    n         = r$n_total,
    stringsAsFactors = FALSE
  )
})

as.data.frame(dr_summary) %>% print()

# --- 5. Save -----------------------------------------------------------------

save(dr_results, dr_summary, file = "dose_response_results.RData")
write.csv(dr_summary, "dose_response_summary.csv", row.names = FALSE)

cat("\nSaved: dose_response_results.RData, dose_response_summary.csv\n")
