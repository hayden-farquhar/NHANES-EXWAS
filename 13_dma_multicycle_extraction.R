# ============================================================================
# 13_dma_multicycle_extraction.R
# Multi-Cycle NHANES Data Extraction for DMA--Uric Acid Analysis
#
# Downloads and harmonises urinary DMA, serum uric acid, covariates,
# metals panel, arsenobetaine, GGT, and metabolic syndrome variables
# across 7 NHANES cycles (2005--2018).
#
# Output: dma_multicycle_harmonised.rds (~12,400 rows x 57 cols)
#
# Run from the repository directory. Requires internet for first run
# (downloads from CDC via nhanesA). Subsequent runs use cached data.
# ============================================================================

library(nhanesA)
library(tidyverse)

REPO_DIR <- getwd()
OUT_DIR <- file.path(REPO_DIR, "dma_uric_acid_outputs")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
CACHE_DIR <- file.path(OUT_DIR, "cycle_cache")
dir.create(CACHE_DIR, recursive = TRUE, showWarnings = FALSE)

cycles <- list(
  list(label = "2005-2006", suffix = "_D"),
  list(label = "2007-2008", suffix = "_E"),
  list(label = "2009-2010", suffix = "_F"),
  list(label = "2011-2012", suffix = "_G"),
  list(label = "2013-2014", suffix = "_H"),
  list(label = "2015-2016", suffix = "_I"),
  list(label = "2017-2018", suffix = "_J")
)

safe_nhanes <- function(tbl_name, max_retries = 3) {
  for (i in seq_len(max_retries)) {
    result <- tryCatch(nhanes(tbl_name), error = function(e) NULL)
    if (!is.null(result) && nrow(result) > 0) return(result)
    if (i < max_retries) Sys.sleep(2)
  }
  NULL
}

as_num <- function(x) suppressWarnings(as.numeric(as.character(x)))

extract_cycle <- function(cyc) {
  s <- cyc$suffix
  label <- cyc$label
  cache_file <- file.path(CACHE_DIR, paste0("cycle_", gsub("-", "_", label), ".rds"))
  if (file.exists(cache_file)) { cat(sprintf("[%s] Cached\n", label)); return(readRDS(cache_file)) }

  cat(sprintf("\n=== Extracting %s ===\n", label))

  tbl_names <- list(
    demo = paste0("DEMO", s), biopro = paste0("BIOPRO", s),
    uas = paste0("UAS", s), bmx = paste0("BMX", s),
    smq = paste0("SMQ", s), cot = paste0("COT", s),
    ghb = paste0("GHB", s), glu = paste0("GLU", s),
    trigly = paste0("TRIGLY", s), bpx = paste0("BPX", s),
    uhm = paste0("UM", s), pbcd = paste0("PBCD", s)
  )
  uhm_alts <- c(paste0("UM", s), paste0("UHM", s))

  tbls <- list()
  for (nm in names(tbl_names)) {
    cat(sprintf("  %s (%s)... ", nm, tbl_names[[nm]]))
    df <- safe_nhanes(tbl_names[[nm]])
    if (is.null(df) && nm == "uhm") for (alt in uhm_alts) { df <- safe_nhanes(alt); if (!is.null(df)) break }
    if (!is.null(df)) { cat(sprintf("OK (%d)\n", nrow(df))); tbls[[nm]] <- df } else cat("FAILED\n")
  }

  if (is.null(tbls$biopro)) {
    cat(sprintf("  BIOPRO missing, trying L40%s... ", s))
    tbls$biopro <- safe_nhanes(paste0("L40", s))
    if (!is.null(tbls$biopro)) cat("OK\n") else cat("FAILED\n")
  }

  if (is.null(tbls$demo) || is.null(tbls$uas) || is.null(tbls$biopro)) {
    cat(sprintf("  [%s] SKIPPING\n", label)); return(NULL)
  }

  merged <- tbls$demo
  for (nm in setdiff(names(tbls), "demo")) {
    if (!is.null(tbls[[nm]])) {
      shared <- setdiff(intersect(names(merged), names(tbls[[nm]])), "SEQN")
      join_df <- tbls[[nm]][, !(names(tbls[[nm]]) %in% shared) | names(tbls[[nm]]) == "SEQN"]
      merged <- merge(merged, join_df, by = "SEQN", all.x = TRUE)
    }
  }

  n_m <- nrow(merged)
  sc <- function(col) if (col %in% names(merged)) as_num(merged[[col]]) else rep(NA_real_, n_m)

  d <- data.frame(SEQN = merged$SEQN, cycle = label, stringsAsFactors = FALSE)
  d$SDMVPSU <- sc("SDMVPSU"); d$SDMVSTRA <- sc("SDMVSTRA")
  d$WTSA2YR <- sc("WTSA2YR"); if (all(is.na(d$WTSA2YR))) d$WTSA2YR <- sc("WTSAS2YR")
  d$RIDAGEYR <- sc("RIDAGEYR")

  gender_raw <- if ("RIAGENDR" %in% names(merged)) as.character(merged$RIAGENDR) else rep(NA_character_, n_m)
  gender_num <- suppressWarnings(as.numeric(gender_raw))
  d$female <- as.integer(grepl("Female", gender_raw, ignore.case = TRUE) | (!is.na(gender_num) & gender_num == 2))

  rcn <- if ("RIDRETH3" %in% names(merged)) "RIDRETH3" else if ("RIDRETH1" %in% names(merged)) "RIDRETH1" else NULL
  if (!is.null(rcn)) {
    rc <- as.character(merged[[rcn]]); rn <- suppressWarnings(as.numeric(rc))
    d$race3 <- case_when(rn %in% 3 | grepl("Non-Hispanic White", rc, TRUE) ~ "NHW",
                         rn %in% 4 | grepl("Non-Hispanic Black", rc, TRUE) ~ "NHB", TRUE ~ "Other")
    d$race6 <- rn
  } else { d$race3 <- "Other"; d$race6 <- NA_real_ }

  d$INDFMPIR <- sc("INDFMPIR"); d$BMXBMI <- sc("BMXBMI")

  s040 <- if ("SMQ040" %in% names(merged)) as.character(merged$SMQ040) else rep(NA_character_, n_m)
  s020 <- if ("SMQ020" %in% names(merged)) as.character(merged$SMQ020) else rep(NA_character_, n_m)
  s040n <- suppressWarnings(as.numeric(s040)); s020n <- suppressWarnings(as.numeric(s020))
  cur <- (!is.na(s040n) & s040n %in% 1:2) | grepl("Every day|Some days", s040, TRUE)
  nev <- (!is.na(s020n) & s020n == 2) | grepl("^No$", s020, TRUE)
  d$smoker <- ifelse(cur, "current", ifelse(nev, "never", "former"))

  d$URXUDMA <- sc("URXUDMA"); d$log_URXUDMA <- log(pmax(d$URXUDMA, 0.01, na.rm = TRUE))
  d$URXUAB <- sc("URXUAB"); d$log_URXUAB <- ifelse(!is.na(d$URXUAB) & d$URXUAB > 0, log(d$URXUAB), NA)
  d$LBXSUA <- sc("LBXSUA"); d$hyperuricemia <- as.integer(d$LBXSUA > 6.8)
  d$LBXGGT <- sc("LBXSGTSI"); if (all(is.na(d$LBXGGT))) d$LBXGGT <- sc("LBXGGT")
  d$log_GGT <- ifelse(!is.na(d$LBXGGT) & d$LBXGGT > 0, log(d$LBXGGT), NA)

  scr <- sc("LBXSCR"); age <- d$RIDAGEYR; fem <- d$female
  kappa <- ifelse(fem == 1, 0.7, 0.9); alpha <- ifelse(fem == 1, -0.241, -0.302)
  d$eGFR <- 142 * pmin(scr/kappa, 1)^alpha * pmax(scr/kappa, 1)^(-1.200) * 0.9938^age * ifelse(fem == 1, 1.012, 1)

  d$URXUCR <- sc("URXUCR"); d$log_URXUCR <- ifelse(!is.na(d$URXUCR) & d$URXUCR > 0, log(d$URXUCR), NA)
  d$LBXCOT <- sc("LBXCOT"); d$LBXGH <- sc("LBXGH"); d$LBXGLU <- sc("LBXGLU"); d$LBXTR <- sc("LBXTR")

  bp_sys <- grep("^BPXSY[0-9]", names(merged), value = TRUE)
  if (length(bp_sys) > 0) {
    m <- sapply(bp_sys, function(x) as_num(merged[[x]])); if (is.matrix(m)) d$SBP_mean <- rowMeans(m, na.rm=TRUE) else d$SBP_mean <- m
  } else d$SBP_mean <- NA_real_
  bp_dia <- grep("^BPXDI[0-9]", names(merged), value = TRUE)
  if (length(bp_dia) > 0) {
    m <- sapply(bp_dia, function(x) as_num(merged[[x]])); if (is.matrix(m)) d$DBP_mean <- rowMeans(m, na.rm=TRUE) else d$DBP_mean <- m
  } else d$DBP_mean <- NA_real_
  d$hypertension <- as.integer(d$SBP_mean >= 130 | d$DBP_mean >= 80)

  for (mv in c("URXUCD","URXUPB","URXUBA","URXUCO","URXUMN","URXUTL","URXUTU")) {
    d[[mv]] <- sc(mv); d[[paste0("log_", mv)]] <- ifelse(!is.na(d[[mv]]) & d[[mv]] > 0, log(d[[mv]]), NA)
  }
  d$URXUAS <- sc("URXUAS3"); if (all(is.na(d$URXUAS))) d$URXUAS <- sc("URXUAS")
  d$log_URXUAS <- ifelse(!is.na(d$URXUAS) & d$URXUAS > 0, log(d$URXUAS), NA)

  d$LBXBCD <- sc("LBXBCD"); d$LBXBPB <- sc("LBXBPB"); d$LBXTHG <- sc("LBXTHG")

  d <- d[!is.na(d$RIDAGEYR) & d$RIDAGEYR >= 18 & d$RIDAGEYR <= 79, ]
  d <- d[!is.na(d$URXUDMA) & !is.na(d$WTSA2YR) & d$WTSA2YR > 0, ]
  cat(sprintf("  [%s] %d rows\n", label, nrow(d)))
  saveRDS(d, cache_file)
  d
}

all_cycles <- compact(map(cycles, ~ tryCatch(extract_cycle(.x), error = function(e) NULL)))
common <- Reduce(intersect, map(all_cycles, names))
stacked <- bind_rows(map(all_cycles, ~ .x[, common]))

stacked$analytic <- complete.cases(stacked[, c("log_URXUDMA","LBXSUA","RIDAGEYR","female","race3","INDFMPIR","BMXBMI","smoker")])
stacked$dma_quartile_pooled <- ntile(stacked$log_URXUDMA, 4)

saveRDS(stacked, file.path(OUT_DIR, "dma_multicycle_harmonised.rds"))
cat(sprintf("\nSaved: %d rows, %d analytic, %d cycles\n", nrow(stacked), sum(stacked$analytic), length(unique(stacked$cycle))))
