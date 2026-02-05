# ============================================================================
# 00_functions.R - Shared utilities for NHANES ExWAS project
# ============================================================================

# --- Packages ----------------------------------------------------------------
required_pkgs <- c("nhanesA", "survey", "tidyverse", "broom", "ggrepel",
                   "forestplot", "pheatmap", "kableExtra", "tableone")

for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
}

library(tidyverse)
library(survey)
library(broom)
library(nhanesA)

# --- CKD-EPI 2021 eGFR (race-free equation) ---------------------------------
calc_eGFR <- function(creatinine, age, female) {
  # creatinine in mg/dL, age in years, female = 1/0
  kappa  <- ifelse(female == 1, 0.7, 0.9)
  alpha  <- ifelse(female == 1, -0.241, -0.302)
  female_mult <- ifelse(female == 1, 1.012, 1.0)

  scr_ratio <- creatinine / kappa
  eGFR <- 142 *
    pmin(scr_ratio, 1)^alpha *
    pmax(scr_ratio, 1)^(-1.200) *
    0.9938^age *
    female_mult
  return(eGFR)
}

# --- Survey-weighted ExWAS model runner --------------------------------------
run_exwas_model <- function(exposure, outcome, data,
                            covariates = "RIDAGEYR + female + race3 + INDFMPIR + BMXBMI + smoker",
                            weight_var = "WTMEC2YR",
                            binary = FALSE) {
  # Drop BMI from covariates if outcome is BMI or waist
  if (outcome %in% c("BMXBMI", "BMXWAIST")) {
    covariates <- gsub("\\+ BMXBMI|BMXBMI \\+", "", covariates)
    covariates <- trimws(gsub("\\s+\\+\\s+\\+", " + ", covariates))
  }

  formula_str <- paste(outcome, "~", exposure, "+", covariates)
  formula_obj <- as.formula(formula_str)

  # Identify all variables needed
  all_vars <- all.vars(formula_obj)
  design_vars <- c("SDMVPSU", "SDMVSTRA", weight_var)
  needed <- unique(c(all_vars, design_vars))
  needed <- needed[needed %in% names(data)]

  sub_dat <- data[complete.cases(data[, needed]), ]
  if (nrow(sub_dat) < 50) return(NULL)

  des <- svydesign(
    id      = ~SDMVPSU,
    strata  = ~SDMVSTRA,
    weights = as.formula(paste0("~", weight_var)),
    nest    = TRUE,
    data    = sub_dat
  )

  if (binary) {
    mod <- svyglm(formula_obj, design = des, family = quasibinomial())
  } else {
    mod <- svyglm(formula_obj, design = des)
  }

  ct <- summary(mod)$coefficients
  if (!(exposure %in% rownames(ct))) return(NULL)

  data.frame(
    exposure = exposure,
    outcome  = outcome,
    beta     = ct[exposure, "Estimate"],
    se       = ct[exposure, "Std. Error"],
    t_value  = ct[exposure, "t value"],
    p_value  = ct[exposure, "Pr(>|t|)"],
    n        = nrow(sub_dat),
    df       = mod$df.residual,
    weight   = weight_var,
    stringsAsFactors = FALSE
  )
}

# --- Weight selector based on variable prefix --------------------------------
select_weight <- function(exposure_var, data) {
  raw <- sub("^log_", "", exposure_var)
  if (grepl("^SS", raw))       wt <- "WTSSBJ2Y"
  else if (grepl("^LBX", raw)) wt <- "WTMEC2YR"
  else if (grepl("^URX", raw)) wt <- "WTSA2YR"
  else                          wt <- "WTMEC2YR"
  # Fall back to WTMEC2YR if chosen weight doesn't exist
  if (!(wt %in% names(data))) wt <- "WTMEC2YR"
  return(wt)
}

# --- Standard covariate strings ----------------------------------------------
COVARS_FULL    <- "RIDAGEYR + female + race3 + INDFMPIR + BMXBMI + smoker"
COVARS_NO_BMI  <- "RIDAGEYR + female + race3 + INDFMPIR + smoker"

# --- Chemical and outcome label lookups --------------------------------------
CHEM_LABELS <- c(
  "log_LBXBPB"  = "Blood lead",
  "log_LBXBCD"  = "Blood cadmium",
  "log_LBXTHG"  = "Blood mercury (total)",
  "log_LBXBSE"  = "Blood selenium",
  "log_LBXBMN"  = "Blood manganese",
  "log_LBXBGM"  = "Methylmercury (blood)",
  "log_LBXBCO"  = "Blood cobalt",
  "log_LBXVBZ"  = "Blood benzene",
  "log_LBXPFNA" = "PFNA",
  "log_LBXPFHS" = "PFHxS",
  "log_LBXPFDE" = "PFDeA",
  "log_LBXPFUA" = "PFUA",
  "log_LBXMPAH" = "Methylparaben (blood)",
  "log_URXMHH"  = "MEHHP (phthalate)",
  "log_URXMOH"  = "MEOHP (phthalate)",
  "log_URXECP"  = "MECPP (phthalate)",
  "log_URXMC1"  = "MCNP (phthalate)",
  "log_URXMEP"  = "MEP (phthalate)",
  "log_URXMIB"  = "MiBP (phthalate)",
  "log_URXMHBP" = "MHBP (phthalate)",
  "log_URXMZP"  = "MBzP (phthalate)",
  "log_URXMNP"  = "MnBP (phthalate)",
  "log_URXMCOH" = "MCOCH (phthalate)",
  "log_URXMHNC" = "MHNCH (phthalate)",
  "log_URXMHP"  = "MEHP (phthalate)",
  "log_URXP01"  = "1-Hydroxynaphthalene (PAH)",
  "log_URXP03"  = "3-Hydroxyfluorene (PAH)",
  "log_URXP04"  = "2-Hydroxyfluorene (PAH)",
  "log_URXP06"  = "1-Hydroxyphenanthrene (PAH)",
  "log_URXP10"  = "1-Hydroxypyrene (PAH)",
  "log_URXUCD"  = "Urinary cadmium",
  "log_URXUCO"  = "Urinary cobalt",
  "log_URXUCS"  = "Urinary cesium",
  "log_URXUPB"  = "Urinary lead",
  "log_URXUTL"  = "Urinary thallium",
  "log_URXUIO"  = "Urinary iodine",
  "log_URXUP8"  = "Urinary perchlorate",
  "log_URXUAS"  = "Total arsenic (urinary)",
  "log_URXUDMA" = "Dimethylarsonic acid (urinary)",
  "log_URXHPM"  = "HPMA (acrolein metabolite)",
  "log_URXOXY"  = "Oxychlordane (urinary)",
  "log_SSGLYP"  = "Glyphosate (surplus serum)",
  "log_SSGENX"  = "GenX (surplus serum)"
)

OUTCOME_LABELS <- c(
  "LBXSATSI"     = "ALT",
  "LBXSASSI"     = "AST",
  "LBXSGB"       = "GGT",
  "LBXSAPSI"     = "Alk Phosphatase",
  "log_LBXSAPSI" = "Alk Phosphatase (log)",
  "LBXSC3SI"     = "Bicarbonate",
  "LBXSBU"       = "BUN",
  "LBXSKSI"      = "Potassium",
  "LBXSNASI"     = "Sodium",
  "LBXSCLSI"     = "Chloride",
  "LBXSPH"       = "Phosphorus",
  "LBXSTB"       = "Total bilirubin",
  "LBXSTP"       = "Total protein",
  "LBXSUA"       = "Uric acid",
  "LBXSTR"       = "Triglycerides",
  "LBXSCH"       = "Total cholesterol (serum)",
  "LBXTC"        = "Total cholesterol",
  "LBDHDD"       = "HDL cholesterol",
  "LBXGH"        = "HbA1c",
  "LBXGLU"       = "Fasting glucose",
  "LBXHSCRP"     = "CRP",
  "log_LBXHSCRP" = "CRP (log)",
  "LBXWBCSI"     = "WBC",
  "LBXRBCSI"     = "RBC count",
  "LBXHGB"       = "Hemoglobin",
  "LBXPLTSI"     = "Platelets",
  "mean_sys"     = "Systolic BP",
  "mean_dia"     = "Diastolic BP",
  "BMXBMI"       = "BMI",
  "BMXWAIST"     = "Waist circumference",
  "PHQ9_total"   = "PHQ-9 depression score",
  "SLD012"       = "Sleep duration",
  "eGFR"         = "eGFR (kidney)",
  "depression_binary" = "Depression (PHQ-9>=10)"
)

cat("00_functions.R loaded.\n")
