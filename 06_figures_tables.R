# ============================================================================
# 06_figures_tables.R
# Publication-quality figures and tables
# ============================================================================

source("00_functions.R")
library(ggrepel)

# Create output directory
dir.create("figures", showWarnings = FALSE)

# --- 1. Load all results -----------------------------------------------------

load("exwas_master_results.RData")   # master, priority_findings
load("validation_results.RData")     # validated
if (file.exists("dose_response_results.RData")) load("dose_response_results.RData")
if (file.exists("sensitivity_results.RData"))   load("sensitivity_results.RData")

# Load primary data for Table 1
if (file.exists("exwas_novelty_screen.RData")) {
  load("exwas_novelty_screen.RData")  # gives `dat`
  dat_full <- dat; rm(dat)
} else if (file.exists("exwas_expanded_results.RData")) {
  load("exwas_expanded_results.RData")  # gives `expanded`
  dat_full <- expanded; rm(expanded)
} else {
  stop("Cannot find primary analysis data.")
}

# =============================================================================
# FIGURE 1: Volcano Plot (main ExWAS result)
# =============================================================================

cat("Creating Figure 1: Volcano plot...\n")

volcano_data <- master %>%
  filter(!is.na(p_value), p_value > 0) %>%
  mutate(
    neg_log_p  = -log10(p_value),
    label_text = ifelse(sig_fdr05,
                        paste0(exp_label, " -> ", out_label), NA),
    point_color = case_when(
      sig_bonf  ~ "Bonferroni",
      sig_fdr05 ~ "FDR < 0.05",
      sig_nominal ~ "Nominal (p<0.05)",
      TRUE ~ "Not significant"
    )
  )

# Threshold lines
fdr_threshold  <- -log10(max(master$p_value[master$sig_fdr05], na.rm = TRUE))
bonf_threshold <- -log10(0.05 / nrow(master))

p_volcano <- ggplot(volcano_data, aes(x = beta, y = neg_log_p)) +
  geom_point(aes(color = point_color), alpha = 0.6, size = 1.5) +
  geom_hline(yintercept = fdr_threshold, linetype = "dashed", color = "blue", alpha = 0.7) +
  geom_hline(yintercept = bonf_threshold, linetype = "dashed", color = "red", alpha = 0.7) +
  geom_text_repel(
    aes(label = label_text),
    size = 2.2, max.overlaps = 20,
    box.padding = 0.3, segment.alpha = 0.5
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

# =============================================================================
# FIGURE 2: Forest Plot of Validated Findings
# =============================================================================

cat("Creating Figure 2: Forest plot...\n")

forest_data <- validated %>%
  filter(status == "ok") %>%
  mutate(
    finding = paste0(exp_label, " -> ", out_label),
    ci_lower_primary = beta_primary - 1.96 * se_val,  # approximation
    ci_upper_primary = beta_primary + 1.96 * se_val
  ) %>%
  arrange(beta_primary)

# Restructure for paired forest plot (primary + validation)
fp_primary <- forest_data %>%
  transmute(
    finding,
    beta = beta_primary,
    se = NA_real_,  # SE not stored in primary summary
    cycle = "2017-2018 (primary)",
    validated = validated
  )

fp_validation <- forest_data %>%
  transmute(
    finding,
    beta = beta_val,
    se = se_val,
    cycle = "2015-2016 (validation)",
    validated = validated
  )

fp_combined <- bind_rows(fp_primary, fp_validation) %>%
  mutate(
    ci_lower = beta - 1.96 * se,
    ci_upper = beta + 1.96 * se,
    finding = factor(finding, levels = unique(forest_data$finding))
  )

p_forest <- ggplot(fp_combined, aes(x = beta, y = finding, color = cycle, shape = cycle)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(position = position_dodge(width = 0.5), size = 2.5) +
  geom_errorbarh(
    aes(xmin = ci_lower, xmax = ci_upper),
    position = position_dodge(width = 0.5),
    height = 0.2, na.rm = TRUE
  ) +
  scale_color_manual(values = c("2017-2018 (primary)" = "#2166AC",
                                "2015-2016 (validation)" = "#B2182B")) +
  labs(
    x = expression("Effect estimate (" * beta * ")"),
    y = NULL,
    title = "Cross-Cycle Validation of ExWAS Findings",
    color = "NHANES Cycle",
    shape = "NHANES Cycle"
  ) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom")

ggsave("figures/fig2_forest.pdf", p_forest, width = 9, height = max(6, nrow(forest_data) * 0.3 + 2))
ggsave("figures/fig2_forest.png", p_forest, width = 9, height = max(6, nrow(forest_data) * 0.3 + 2), dpi = 300)

# =============================================================================
# FIGURE 3: Dose-Response Curves
# =============================================================================

if (exists("dr_results") && length(dr_results) > 0) {
  cat("Creating Figure 3: Dose-response curves...\n")

  # Combine all quartile estimates
  dr_plot_data <- map_dfr(dr_results, function(r) {
    r$q_estimates %>%
      mutate(
        finding = paste0(r$exp_label, " -> ", r$out_label),
        p_trend = r$p_trend
      )
  })

  p_dr <- ggplot(dr_plot_data, aes(x = quartile, y = beta)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    geom_point(size = 2.5) +
    geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.15) +
    geom_line(alpha = 0.5) +
    geom_text(
      data = dr_plot_data %>% filter(quartile == 4),
      aes(label = sprintf("p-trend=%.3f", p_trend)),
      hjust = -0.1, size = 2.8, nudge_y = 0
    ) +
    facet_wrap(~ finding, scales = "free_y", ncol = 3) +
    scale_x_continuous(breaks = 1:4, labels = paste0("Q", 1:4)) +
    labs(
      x = "Exposure Quartile",
      y = expression("Change in outcome vs. Q1 (" * beta * ")"),
      title = "Dose-Response Relationships",
      subtitle = "Quartile analysis with Q1 as reference"
    ) +
    theme_minimal(base_size = 10) +
    theme(strip.text = element_text(size = 8))

  n_facets <- length(unique(dr_plot_data$finding))
  fig_height <- ceiling(n_facets / 3) * 3 + 1

  ggsave("figures/fig3_dose_response.pdf", p_dr, width = 10, height = fig_height)
  ggsave("figures/fig3_dose_response.png", p_dr, width = 10, height = fig_height, dpi = 300)
}

# =============================================================================
# FIGURE 4: Chemical-Outcome Heatmap
# =============================================================================

cat("Creating Figure 4: Heatmap...\n")

# Create matrix of -log10(p) * sign(beta) for nominally significant results
# Deduplicate: keep most significant result per exposure-outcome pair (across rounds)
heatmap_data <- master %>%
  filter(sig_nominal) %>%
  mutate(
    signed_log_p = -log10(p_value) * sign(beta),
    exp_short = coalesce(exp_label, exposure),
    out_short = coalesce(out_label, outcome)
  ) %>%
  group_by(exp_short, out_short) %>%
  slice_min(p_value, n = 1, with_ties = FALSE) %>%
  ungroup()

# Pivot to matrix
heat_matrix <- heatmap_data %>%
  select(exp_short, out_short, signed_log_p) %>%
  pivot_wider(names_from = out_short, values_from = signed_log_p) %>%
  column_to_rownames("exp_short") %>%
  as.matrix()

# Only plot if matrix is reasonable size
if (nrow(heat_matrix) > 2 && ncol(heat_matrix) > 2) {
  # Replace NA with 0 for heatmap
  heat_matrix[is.na(heat_matrix)] <- 0

  # Limit to exposures/outcomes with at least 2 associations
  row_keep <- rowSums(heat_matrix != 0) >= 2
  col_keep <- colSums(heat_matrix != 0) >= 2
  heat_sub <- heat_matrix[row_keep, col_keep, drop = FALSE]

  if (nrow(heat_sub) > 2 && ncol(heat_sub) > 2) {
    pdf("figures/fig4_heatmap.pdf", width = 12, height = max(8, nrow(heat_sub) * 0.3 + 2))
    pheatmap::pheatmap(
      heat_sub,
      color = colorRampPalette(c("#B2182B", "white", "#2166AC"))(100),
      breaks = seq(-max(abs(heat_sub)), max(abs(heat_sub)), length.out = 101),
      main = "Chemical-Outcome Association Heatmap\n(Signed -log10 p-value, nominally significant only)",
      fontsize_row = 7,
      fontsize_col = 8,
      clustering_method = "ward.D2"
    )
    dev.off()
    cat("  Heatmap saved.\n")
  }
}

# =============================================================================
# FIGURE 5: Sensitivity Analysis Forest Plot
# =============================================================================

if (exists("sens_df")) {
  cat("Creating Figure 5: Sensitivity analysis...\n")

  sens_plot_data <- sens_df %>%
    filter(!is.na(beta)) %>%
    mutate(
      finding = paste0(exp_label, " -> ", out_label),
      ci_lower = beta - 1.96 * se,
      ci_upper = beta + 1.96 * se
    )

  # One panel per finding
  findings_to_plot <- unique(sens_plot_data$finding)

  for (f in findings_to_plot) {
    sub <- sens_plot_data %>% filter(finding == f)

    p_sens <- ggplot(sub, aes(x = beta, y = reorder(analysis, beta))) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
      geom_point(aes(color = p_value < 0.05), size = 2.5) +
      geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper), height = 0.2) +
      scale_color_manual(values = c("TRUE" = "blue", "FALSE" = "grey50"),
                         name = "p < 0.05") +
      labs(x = expression(beta), y = NULL, title = f) +
      theme_minimal(base_size = 10)

    fname <- gsub("[^a-zA-Z0-9]", "_", f)
    ggsave(sprintf("figures/fig5_sensitivity_%s.png", fname),
           p_sens, width = 8, height = 4, dpi = 300)
  }
}

# =============================================================================
# TABLE 1: Weighted Demographics
# =============================================================================

cat("Creating Table 1: Weighted demographics...\n")

if (requireNamespace("tableone", quietly = TRUE)) {
  library(tableone)

  # Build survey design with full sample
  tbl1_vars <- c("RIDAGEYR", "female", "race3", "INDFMPIR", "BMXBMI",
                  "smoker", "BMXWAIST")
  tbl1_vars <- tbl1_vars[tbl1_vars %in% names(dat_full)]

  tbl1_data <- dat_full[complete.cases(dat_full[, c(tbl1_vars, "SDMVPSU", "SDMVSTRA", "WTMEC2YR")]), ]

  des_tbl1 <- svydesign(
    id = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~WTMEC2YR,
    nest = TRUE, data = tbl1_data
  )

  cat_vars <- c("female", "race3", "smoker")
  cat_vars <- cat_vars[cat_vars %in% tbl1_vars]

  tbl1 <- svyCreateTableOne(
    vars    = tbl1_vars,
    factorVars = cat_vars,
    data    = des_tbl1
  )

  tbl1_print <- print(tbl1, showAllLevels = TRUE, smd = FALSE, printToggle = FALSE)
  write.csv(tbl1_print, "figures/table1_demographics.csv")
  cat("  Table 1 saved to figures/table1_demographics.csv\n")
}

# =============================================================================
# TABLE 2: Summary of FDR-Significant Findings
# =============================================================================

cat("Creating Table 2: Significant findings...\n")

tbl2 <- master %>%
  filter(sig_fdr05) %>%
  select(exp_label, out_label, chem_class, out_domain,
         beta, se, p_value, p_fdr_global, n, direction) %>%
  arrange(p_fdr_global) %>%
  rename(
    `Chemical`        = exp_label,
    `Outcome`         = out_label,
    `Chemical Class`  = chem_class,
    `Outcome Domain`  = out_domain,
    `Beta`            = beta,
    `SE`              = se,
    `P-value`         = p_value,
    `FDR`             = p_fdr_global,
    `N`               = n,
    `Direction`       = direction
  )

write.csv(tbl2, "figures/table2_significant_findings.csv", row.names = FALSE)

# =============================================================================
# TABLE 3: Validation Results
# =============================================================================

cat("Creating Table 3: Validation results...\n")

tbl3 <- validated %>%
  filter(status == "ok") %>%
  select(exp_label, out_label, beta_primary, beta_val, p_val,
         n_primary, n_val, direction_match, validated) %>%
  arrange(desc(validated), p_val) %>%
  rename(
    `Chemical`           = exp_label,
    `Outcome`            = out_label,
    `Beta (2017-18)`     = beta_primary,
    `Beta (2015-16)`     = beta_val,
    `P (2015-16)`        = p_val,
    `N (2017-18)`        = n_primary,
    `N (2015-16)`        = n_val,
    `Direction Match`    = direction_match,
    `Validated`          = validated
  )

write.csv(tbl3, "figures/table3_validation.csv", row.names = FALSE)

cat("\nAll figures and tables saved to figures/ directory.\n")
