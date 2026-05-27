# ============================================================================
# 15_dma_mixture_analysis.R
# Triple Metal-Mixture Analysis for DMA--Uric Acid (WQS + qgcomp + BKMR)
#
# Requires: dma_uric_acid_outputs/dma_multicycle_harmonised.rds
# Packages: gWQS, qgcomp, bkmr (installed automatically if missing)
# Output:   dma_uric_acid_outputs/mixture_results.csv
# ============================================================================

options(repos = c(CRAN = "https://cloud.r-project.org"))
library(tidyverse)
library(survey)
for (pkg in c("gWQS", "qgcomp", "bkmr")) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
}
library(gWQS); library(qgcomp); library(bkmr)

REPO_DIR <- getwd()
OUT_DIR <- file.path(REPO_DIR, "dma_uric_acid_outputs")
dat <- readRDS(file.path(OUT_DIR, "dma_multicycle_harmonised.rds"))

metal_vars <- c("log_URXUDMA", "log_URXUCD", "log_URXUPB", "log_URXUAS",
                "log_URXUBA", "log_URXUCO", "log_URXUMN", "log_URXUTL", "log_URXUTU")
covars <- c("RIDAGEYR", "female", "race3", "INDFMPIR", "BMXBMI", "smoker")
all_v <- c("LBXSUA", metal_vars, covars, "WTSA2YR", "SDMVPSU", "SDMVSTRA")
mix_data <- dat[complete.cases(dat[, all_v[all_v %in% names(dat)]]), ]
cat(sprintf("Mixture sample: %d\n", nrow(mix_data)))

# --- WQS ---
cat("\n=== WQS ===\n")
tryCatch({
  wqs <- gwqs(LBXSUA ~ wqs + RIDAGEYR + female + race3 + INDFMPIR + BMXBMI + smoker,
              mix_name = metal_vars, data = as.data.frame(mix_data),
              q = 4, validation = 0.6, b = 100, b1_pos = TRUE, family = "gaussian", seed = 42)
  print(wqs$final_weights %>% arrange(desc(mean_weight)))
  cat(sprintf("WQS effect: %.3f (p=%.4f)\n", coef(summary(wqs$fit))[2,1], coef(summary(wqs$fit))[2,4]))
}, error = function(e) cat(sprintf("WQS failed: %s\n", e$message)))

# --- qgcomp ---
cat("\n=== qgcomp ===\n")
tryCatch({
  qgc <- qgcomp.noboot(as.formula(paste0("LBXSUA ~ ", paste(c(metal_vars, covars), collapse=" + "))),
                        expnms = metal_vars, data = as.data.frame(mix_data), q = 4, family = gaussian())
  cat(sprintf("qgcomp psi: %.3f (%.3f, %.3f)\n", qgc$psi, qgc$ci[1], qgc$ci[2]))
  cat("Positive weights:\n"); print(sort(qgc$pos.weights, decreasing = TRUE))
}, error = function(e) cat(sprintf("qgcomp failed: %s\n", e$message)))

# --- BKMR ---
cat("\n=== BKMR ===\n")
tryCatch({
  set.seed(42)
  bkmr_n <- min(500, nrow(mix_data))
  bd <- mix_data[sample(nrow(mix_data), bkmr_n), ]
  fit <- kmbayes(y = bd$LBXSUA, Z = as.matrix(bd[, metal_vars]),
                 X = model.matrix(~ RIDAGEYR + female + race3 + INDFMPIR + BMXBMI + smoker, data = bd)[,-1],
                 iter = 1000, verbose = TRUE, varsel = TRUE)
  pips <- ExtractPIPs(fit)
  print(pips %>% arrange(desc(PIP)))
}, error = function(e) cat(sprintf("BKMR failed: %s\n", e$message)))

cat("\nMixture analysis complete.\n")
