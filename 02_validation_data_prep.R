# ============================================================================
# 02_validation_data_prep.R
# Download and prepare NHANES 2015-2016 data for cross-cycle validation
# ============================================================================

source("00_functions.R")

# --- 1. Define NHANES 2015-2016 tables to download ---------------------------
# Suffix _I = 2015-2016 cycle (vs _J = 2017-2018)

tables_to_download <- list(
  # Demographics & weights
  DEMO    = "DEMO_I",

  # Health outcomes
  BIOPRO  = "BIOPRO_I",   # Biochemistry (ALT, AST, ALP, BUN, creatinine, etc.)
  CBC     = "CBC_I",       # Complete blood count
  GHB     = "GHB_I",       # Glycohemoglobin (HbA1c)
  GLU     = "GLU_I",       # Fasting glucose
  TCHOL   = "TCHOL_I",     # Total cholesterol
  HDL     = "HDL_I",       # HDL cholesterol
  TRIGLY  = "TRIGLY_I",    # Triglycerides
  HSCRP   = "HSCRP_I",    # C-reactive protein
  BMX     = "BMX_I",       # Body measures (BMI, waist)
  BPX     = "BPX_I",       # Blood pressure
  DPQ     = "DPQ_I",       # Depression (PHQ-9)
  SLQ     = "SLQ_I",       # Sleep

  # Chemical exposures - blood metals
  PBCD    = "PBCD_I",      # Lead, cadmium, total mercury, selenium, manganese

  # Chemical exposures - urinary metals
  UHM     = "UM_I",        # Urinary heavy metals (Cd, Co, Cs, Pb, Tl, etc.)
  UIO     = "UIO_I",       # Urinary iodine

  # Chemical exposures - phthalates
  PHTHTE  = "PHTHTE_I",    # Phthalate metabolites

  # Chemical exposures - PAHs
  PAH     = "PAH_I",       # Urinary PAH metabolites

  # Chemical exposures - VOCs
  VOCWB   = "VOCWB_I",     # Volatile organic compounds (whole blood)
  UVOC    = "UVOC2_I",     # Urinary VOC metabolites

  # Chemical exposures - other
  PERNT   = "PERNT_I",     # Perchlorate, nitrate, thiocyanate
  UAS     = "UAS_I",       # Urinary arsenic species

  # Chemical exposures - methylmercury
  IHGEM   = "IHGEM_I",     # Inorganic/methylmercury (LBXBGM)

  # Smoking / cotinine
  COT     = "COT_I",       # Cotinine
  SMQ     = "SMQ_I"        # Smoking questionnaire
)

# --- 2. Download all tables --------------------------------------------------

downloaded <- list()
failed <- c()

for (tbl_name in names(tables_to_download)) {
  tbl_id <- tables_to_download[[tbl_name]]
  cat(sprintf("Downloading %s (%s)... ", tbl_name, tbl_id))
  tryCatch({
    df <- nhanes(tbl_id)
    downloaded[[tbl_name]] <- df
    cat(sprintf("OK (%d rows, %d cols)\n", nrow(df), ncol(df)))
  }, error = function(e) {
    cat(sprintf("FAILED: %s\n", e$message))
    failed <<- c(failed, tbl_id)
  })
}

# Remove empty/corrupt downloads (downloaded length 0 issue)
empty_tbls <- c()
for (nm in names(downloaded)) {
  if (is.null(downloaded[[nm]]) || nrow(downloaded[[nm]]) == 0) {
    empty_tbls <- c(empty_tbls, nm)
  }
}
if (length(empty_tbls) > 0) {
  cat(sprintf("Removing %d empty downloads: %s\n", length(empty_tbls),
              paste(empty_tbls, collapse = ", ")))
  downloaded <- downloaded[!names(downloaded) %in% empty_tbls]
}

# Try alternate table names for failed/empty tables
alt_names <- list(
  COT  = c("COT_I", "COTNAL_I"),
  SMQ  = c("SMQ_I", "SMQRTU_I"),
  UVOC = c("UVOC_I", "UVOC2_I")
)
for (nm in names(alt_names)) {
  if (nm %in% names(downloaded)) next  # already have it
  for (alt_id in alt_names[[nm]]) {
    if (alt_id %in% failed) next
    cat(sprintf("Retrying %s as %s... ", nm, alt_id))
    tryCatch({
      df <- nhanes(alt_id)
      if (!is.null(df) && nrow(df) > 0) {
        downloaded[[nm]] <- df
        cat(sprintf("OK (%d rows, %d cols)\n", nrow(df), ncol(df)))
        break
      } else {
        cat("empty\n")
      }
    }, error = function(e) cat(sprintf("FAILED\n")))
  }
}

cat(sprintf("\nDownloaded: %d tables, Failed: %d tables\n",
            length(downloaded), length(failed)))
if (length(failed) > 0) cat("Failed tables:", paste(failed, collapse = ", "), "\n")

# --- 3. Merge all tables by SEQN ---------------------------------------------

val_data <- downloaded[["DEMO"]]
cat(sprintf("\nStarting merge with DEMO: %d rows\n", nrow(val_data)))

for (tbl_name in setdiff(names(downloaded), "DEMO")) {
  df <- downloaded[[tbl_name]]
  if ("SEQN" %in% names(df)) {
    # Remove any columns already in val_data (except SEQN)
    dups <- setdiff(intersect(names(val_data), names(df)), "SEQN")
    if (length(dups) > 0) {
      df <- df[, !(names(df) %in% dups)]
    }
    val_data <- left_join(val_data, df, by = "SEQN")
    cat(sprintf("  + %s: now %d cols\n", tbl_name, ncol(val_data)))
  }
}

cat(sprintf("\nMerged dataset: %d rows x %d cols\n", nrow(val_data), ncol(val_data)))

# --- 4. Create derived variables (matching 2017-2018 pipeline) ---------------

cat("\nCreating derived variables...\n")

# Age
val_data$RIDAGEYR <- as.numeric(val_data$RIDAGEYR)

# Sex: female indicator
val_data$female <- as.numeric(val_data$RIAGENDR == 2)

# Race: collapse to 3 categories (handle numeric, character, or labelled)
ridreth3 <- as.character(val_data$RIDRETH3)
val_data$race3 <- case_when(
  ridreth3 %in% c("3", "Non-Hispanic White")  ~ "White",
  ridreth3 %in% c("4", "Non-Hispanic Black")  ~ "Black",
  TRUE ~ "Other"
)
val_data$race3 <- factor(val_data$race3, levels = c("White", "Black", "Other"))
cat(sprintf("  race3: %s\n", paste(names(table(val_data$race3)), table(val_data$race3),
            sep = "=", collapse = ", ")))

# Income-to-poverty ratio
val_data$INDFMPIR <- as.numeric(val_data$INDFMPIR)

# Smoking: derive from questionnaire + cotinine
if ("SMQ020" %in% names(val_data) && "SMQ040" %in% names(val_data)) {
  smq020 <- as.character(val_data$SMQ020)
  smq040 <- as.character(val_data$SMQ040)
  val_data$smoking_status <- case_when(
    smq020 %in% c("2", "No")  ~ "Never",
    smq040 %in% c("1", "2", "Every day", "Some days") ~ "Current",
    smq040 %in% c("3", "Not at all") ~ "Former",
    TRUE ~ NA_character_
  )
  val_data$smoker <- as.numeric(val_data$smoking_status == "Current")
  # If SMQ labels didn't match, fall back to cotinine
  if (sum(!is.na(val_data$smoker)) < 100 && "LBXCOT" %in% names(val_data)) {
    cat("  WARNING: SMQ labels may not match, falling back to cotinine\n")
    cot_vals <- as.numeric(as.character(val_data$LBXCOT))
    val_data$smoker <- as.numeric(!is.na(cot_vals) & cot_vals > 10)
    val_data$smoker[is.na(cot_vals)] <- NA
  }
  cat(sprintf("  Smoking: derived from SMQ (%d non-missing)\n",
              sum(!is.na(val_data$smoker))))
} else if ("LBXCOT" %in% names(val_data)) {
  # Fallback: use cotinine
  cot_vals <- as.numeric(val_data$LBXCOT)
  val_data$smoker <- as.numeric(!is.na(cot_vals) & cot_vals > 10)
  val_data$smoker[is.na(cot_vals)] <- NA
  cat(sprintf("  Smoking: derived from cotinine (%d non-missing)\n",
              sum(!is.na(val_data$smoker))))
} else {
  # Neither available — set to NA (will limit covariate adjustment)
  val_data$smoker <- NA_real_
  cat("  WARNING: No smoking data available (COT and SMQ both missing)\n")
  cat("  Smoker covariate will be NA — models will need to handle this\n")
}

# BMI
val_data$BMXBMI <- as.numeric(val_data$BMXBMI)

# eGFR from serum creatinine (CKD-EPI 2021)
if ("LBXSCR" %in% names(val_data)) {
  val_data$eGFR <- calc_eGFR(
    creatinine = as.numeric(val_data$LBXSCR),
    age        = val_data$RIDAGEYR,
    female     = val_data$female
  )
  cat(sprintf("  eGFR: %d non-missing\n", sum(!is.na(val_data$eGFR))))
}

# Blood pressure: average of available readings
bp_sys_cols <- grep("^BPXSY[0-9]", names(val_data), value = TRUE)
bp_dia_cols <- grep("^BPXDI[0-9]", names(val_data), value = TRUE)
if (length(bp_sys_cols) > 0) {
  bp_sys_mat <- sapply(val_data[, bp_sys_cols], function(x) as.numeric(as.character(x)))
  bp_dia_mat <- sapply(val_data[, bp_dia_cols], function(x) as.numeric(as.character(x)))
  # Replace 0 diastolic with NA (artifact)
  bp_dia_mat[bp_dia_mat == 0] <- NA
  val_data$mean_sys <- rowMeans(bp_sys_mat, na.rm = TRUE)
  val_data$mean_dia <- rowMeans(bp_dia_mat, na.rm = TRUE)
  val_data$mean_sys[is.nan(val_data$mean_sys)] <- NA
  val_data$mean_dia[is.nan(val_data$mean_dia)] <- NA
}

# PHQ-9 depression score
dpq_cols <- grep("^DPQ0[1-9]0$", names(val_data), value = TRUE)
if (length(dpq_cols) >= 9) {
  # nhanesA may return text labels; recode to numeric
  recode_dpq <- function(x) {
    x <- as.character(x)
    dplyr::case_when(
      x %in% c("0", "Not at all")               ~ 0L,
      x %in% c("1", "Several days")              ~ 1L,
      x %in% c("2", "More than half the days")   ~ 2L,
      x %in% c("3", "Nearly every day")          ~ 3L,
      TRUE ~ NA_integer_  # 7=Refused, 9=Don't know, other
    )
  }
  dpq_mat <- sapply(val_data[, dpq_cols], recode_dpq)
  val_data$PHQ9_total <- rowSums(dpq_mat, na.rm = FALSE)
  val_data$depression_binary <- as.numeric(val_data$PHQ9_total >= 10)
  cat(sprintf("  PHQ-9: %d non-missing (range %d-%d)\n",
              sum(!is.na(val_data$PHQ9_total)),
              min(val_data$PHQ9_total, na.rm = TRUE),
              max(val_data$PHQ9_total, na.rm = TRUE)))
}

# Log-transform outcome variables that were logged in primary analysis
if ("LBXSAPSI" %in% names(val_data)) {
  val_data$log_LBXSAPSI <- log(as.numeric(as.character(val_data$LBXSAPSI)))
}
if ("LBXHSCRP" %in% names(val_data)) {
  val_data$log_LBXHSCRP <- log(as.numeric(as.character(val_data$LBXHSCRP)))
}

# --- 5. Log-transform chemical exposures -------------------------------------

cat("\nLog-transforming chemical exposures...\n")

# Identify all potential chemical exposure variables
chem_prefixes <- c("LBXBPB", "LBXBCD", "LBXTHG", "LBXBSE", "LBXBMN",
                   "LBXBGM", "LBXBCO", "LBXVBZ",
                   "URXMHH", "URXMOH", "URXECP", "URXMC1", "URXMEP",
                   "URXMIB", "URXMHBP", "URXMZP", "URXMNP",
                   "URXP01", "URXP03", "URXP04", "URXP06", "URXP10",
                   "URXUCD", "URXUCO", "URXUCS", "URXUPB", "URXUTL",
                   "URXUIO", "URXUP8", "URXUAS", "URXUDMA",
                   "URXHPM", "URXOXY")

n_logged <- 0
for (v in chem_prefixes) {
  if (v %in% names(val_data)) {
    vals <- as.numeric(as.character(val_data[[v]]))
    n_valid <- sum(!is.na(vals))
    if (n_valid < 50) next

    n_pos <- sum(vals > 0, na.rm = TRUE)
    if (n_pos / n_valid < 0.6) next

    # Handle non-positive values (below LOD)
    min_pos <- min(vals[vals > 0], na.rm = TRUE)
    vals_clean <- ifelse(!is.na(vals) & vals <= 0, min_pos / 2, vals)
    log_name <- paste0("log_", v)
    val_data[[log_name]] <- log(vals_clean)
    n_logged <- n_logged + 1
  }
}
cat(sprintf("  Log-transformed: %d exposures\n", n_logged))

# --- 6. Check variable availability ------------------------------------------

cat("\n=== VARIABLE AVAILABILITY CHECK (2015-2016) ===\n")

# Exposures that were FDR-significant in 2017-2018
sig_exposures <- c("log_LBXTHG", "log_LBXBPB", "log_LBXBSE", "log_LBXBMN",
                   "log_LBXBCD",
                   "log_URXMHH", "log_URXMOH", "log_URXECP",
                   "log_URXUCD", "log_URXUCO", "log_URXUCS",
                   "log_URXUPB", "log_URXUTL",
                   "log_URXHPM", "log_URXOXY")

sig_outcomes <- c("log_LBXSAPSI", "LBXSTB", "eGFR", "LBXSBU",
                  "LBXGH", "LBXTC", "LBXRBCSI", "LBXHGB",
                  "BMXBMI", "BMXWAIST")

cat("\nExposures:\n")
for (v in sig_exposures) {
  avail <- v %in% names(val_data)
  n_val <- if (avail) sum(!is.na(val_data[[v]])) else 0
  cat(sprintf("  %s: %s (n=%d)\n", v, ifelse(avail, "YES", "NO"), n_val))
}

cat("\nOutcomes:\n")
for (v in sig_outcomes) {
  avail <- v %in% names(val_data)
  n_val <- if (avail) sum(!is.na(val_data[[v]])) else 0
  cat(sprintf("  %s: %s (n=%d)\n", v, ifelse(avail, "YES", "NO"), n_val))
}

cat("\nCovariates:\n")
for (v in c("RIDAGEYR", "female", "race3", "INDFMPIR", "BMXBMI", "smoker",
            "SDMVPSU", "SDMVSTRA", "WTMEC2YR", "WTSA2YR")) {
  avail <- v %in% names(val_data)
  n_val <- if (avail) sum(!is.na(val_data[[v]])) else 0
  cat(sprintf("  %s: %s (n=%d)\n", v, ifelse(avail, "YES", "NO"), n_val))
}

# --- 7. Save -----------------------------------------------------------------

save(val_data, file = "nhanes_2015_2016_validation.RData")
cat(sprintf("\nSaved: nhanes_2015_2016_validation.RData (%d obs x %d vars)\n",
            nrow(val_data), ncol(val_data)))
