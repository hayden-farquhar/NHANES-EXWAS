# ============================================================================
# 14_dma_enhanced_analyses.R
# Enhanced Multi-Cycle Analyses for DMA--Uric Acid
#
# Analyses: cycle-specific meta-analysis, mediation (eGFR + GGT),
# PAF, interaction tests, arsenobetaine adjustment, temporal trend,
# metabolic syndrome covariates.
#
# Requires: dma_uric_acid_outputs/dma_multicycle_harmonised.rds
# Output:   dma_uric_acid_outputs/*.csv
# ============================================================================

library(tidyverse)
library(survey)
library(broom)
library(meta)

REPO_DIR <- getwd()
DATA_PATH <- file.path(REPO_DIR, "dma_uric_acid_outputs", "dma_multicycle_harmonised.rds")
OUT_DIR <- file.path(REPO_DIR, "dma_uric_acid_outputs")
stopifnot(file.exists(DATA_PATH))

dat <- readRDS(DATA_PATH)
analytic <- dat[dat$analytic == TRUE, ]
cat(sprintf("Loaded: %d analytic rows, %d cycles\n", nrow(analytic), length(unique(analytic$cycle))))

make_design <- function(data) {
  svydesign(id = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~WTSA2YR, nest = TRUE, data = data)
}

core_formula <- LBXSUA ~ log_URXUDMA + RIDAGEYR + female + race3 + INDFMPIR + BMXBMI + smoker

# --- A1: Cycle-specific + meta-analysis ---
cat("\n=== Cycle-specific meta-analysis ===\n")
cycle_res <- analytic %>% group_by(cycle) %>%
  group_modify(~ { tryCatch({
    des <- make_design(.x)
    mod <- svyglm(core_formula, design = des)
    tidy(mod, conf.int = TRUE) %>% filter(term == "log_URXUDMA") %>% mutate(n = nrow(.x))
  }, error = function(e) tibble()) }) %>% ungroup() %>% filter(!is.na(estimate))

ma <- metagen(TE = cycle_res$estimate, seTE = cycle_res$std.error,
              studlab = cycle_res$cycle, sm = "MD", method.tau = "REML",
              common = FALSE, random = TRUE)
cat(sprintf("Pooled: %.3f (%.3f, %.3f), I2=%.1f%%\n", ma$TE.random, ma$lower.random, ma$upper.random, ma$I2*100))
write.csv(cycle_res %>% select(cycle, estimate, std.error, conf.low, conf.high, p.value, n),
          file.path(OUT_DIR, "cycle_specific_results.csv"), row.names = FALSE)

# --- A3/A4: Mediation ---
cat("\n=== Mediation ===\n")
for (mediator in c("eGFR", "log_GGT")) {
  md <- analytic[!is.na(analytic[[mediator]]), ]
  des <- make_design(md)
  te <- coef(svyglm(core_formula, design = des))["log_URXUDMA"]
  nde_f <- update(core_formula, paste0(". ~ . + ", mediator))
  nde <- coef(svyglm(nde_f, design = des))["log_URXUDMA"]
  nie <- te - nde
  cat(sprintf("  %s: TE=%.4f, NDE=%.4f, NIE=%.4f, prop=%.1f%%\n", mediator, te, nde, nie, 100*nie/te))
}

# --- A5: PAF ---
cat("\n=== PAF ===\n")
analytic$high_dma <- as.integer(analytic$dma_quartile_pooled >= 3)
paf_des <- make_design(analytic)
paf_mod <- svyglm(hyperuricemia ~ high_dma + RIDAGEYR + female + race3 + INDFMPIR + BMXBMI + smoker,
                   design = paf_des, family = quasibinomial())
or_h <- exp(coef(paf_mod)["high_dma"])
p_exp <- svymean(~high_dma, paf_des)[1]
paf <- (p_exp * (or_h - 1)) / (p_exp * (or_h - 1) + 1)
cat(sprintf("PAF: %.1f%%\n", 100 * paf))

# --- A6: Interactions ---
cat("\n=== Interactions ===\n")
int_des <- make_design(analytic)
for (term in c("female", "BMXBMI")) {
  mod <- svyglm(as.formula(paste0("LBXSUA ~ log_URXUDMA * ", term,
    " + RIDAGEYR + race3 + INDFMPIR + smoker", if (term != "BMXBMI") " + BMXBMI" else "")),
    design = int_des)
  int_row <- tidy(mod) %>% filter(grepl(paste0("log_URXUDMA:", term), term))
  cat(sprintf("  DMA x %s: beta=%.4f, p=%.4f\n", term, int_row$estimate, int_row$p.value))
}

# --- A7: Arsenobetaine ---
cat("\n=== Arsenobetaine ===\n")
ab <- analytic[!is.na(analytic$log_URXUAB), ]
ab_des <- make_design(ab)
b_no <- coef(svyglm(core_formula, design = ab_des))["log_URXUDMA"]
b_ab <- coef(svyglm(update(core_formula, . ~ . + log_URXUAB), design = ab_des))["log_URXUDMA"]
cat(sprintf("Without AsB: %.3f, With: %.3f, Attenuation: %.0f%%\n", b_no, b_ab, 100*(1 - b_ab/b_no)))

# --- A8: Temporal trend ---
cat("\n=== Temporal trend ===\n")
analytic$cycle_midyear <- as.numeric(substr(analytic$cycle, 1, 4)) + 0.5
trend_des <- make_design(analytic)
mod_t <- svyglm(LBXSUA ~ log_URXUDMA * cycle_midyear + RIDAGEYR + female + race3 + INDFMPIR + BMXBMI + smoker,
                design = trend_des)
t_int <- tidy(mod_t) %>% filter(grepl("cycle_midyear", term))
cat(sprintf("DMA x year: beta=%.5f, p=%.4f\n", t_int$estimate, t_int$p.value))

# --- A9: Sensitivity ---
cat("\n=== Sensitivity analyses ===\n")
run_sens <- function(label, data, formula) {
  tryCatch({
    d <- data[complete.cases(data[, all.vars(formula)]), ]
    if (nrow(d) < 500) return(NULL)
    des <- make_design(d)
    res <- tidy(svyglm(formula, design = des), conf.int = TRUE) %>% filter(term == "log_URXUDMA")
    tibble(label = label, beta = res$estimate, lo = res$conf.low, hi = res$conf.high, p = res$p.value, n = nrow(d))
  }, error = function(e) NULL)
}

core_c <- LBXSUA ~ log_URXUDMA + RIDAGEYR + female + race3 + INDFMPIR + BMXBMI + smoker + cycle
sens <- compact(list(
  run_sens("Primary", analytic, core_c),
  run_sens("Females", analytic[analytic$female==1,], LBXSUA ~ log_URXUDMA + RIDAGEYR + race3 + INDFMPIR + BMXBMI + smoker + cycle),
  run_sens("Males", analytic[analytic$female==0,], LBXSUA ~ log_URXUDMA + RIDAGEYR + race3 + INDFMPIR + BMXBMI + smoker + cycle),
  run_sens("+AsB", analytic[!is.na(analytic$log_URXUAB),], update(core_c, . ~ . + log_URXUAB)),
  run_sens("+Creatinine", analytic[!is.na(analytic$log_URXUCR),], update(core_c, . ~ . + log_URXUCR)),
  run_sens("+eGFR", analytic[!is.na(analytic$eGFR),], update(core_c, . ~ . + eGFR)),
  run_sens("+HbA1c", analytic[!is.na(analytic$LBXGH),], update(core_c, . ~ . + LBXGH)),
  run_sens("+GGT", analytic[!is.na(analytic$log_GGT),], update(core_c, . ~ . + log_GGT))
))
sens_df <- bind_rows(sens)
print(sens_df)
write.csv(sens_df, file.path(OUT_DIR, "sensitivity_results.csv"), row.names = FALSE)
cat("\nDone.\n")
