# NHANES ExWAS 2017-2018: Environment-Wide Association Study

Environment-wide association study (ExWAS) of 92 chemical biomarkers and 48 health outcomes using NHANES 2017-2018, with cross-cycle validation in NHANES 2015-2016.

## Abstract

We screened 2,796 chemical-health associations in a nationally representative US sample (NHANES 2017-2018) using survey-weighted regression with FDR correction. Of 26 globally significant findings (FDR < 0.05), 15 validated in an independent cycle (2015-2016) with concordant direction and nominal significance. All 15 showed significant dose-response trends (14/15 monotonic) and were robust across nine sensitivity analyses. Three associations were classified as HIGH novelty: dimethylarsonic acid with uric acid, perchlorate with blood urea nitrogen, and methylmercury with waist circumference. Five additional associations were MODERATE novelty.

## Repository Contents

### Analysis Scripts (run in order)

| Script | Purpose | Internet required |
|--------|---------|:-:|
| `00_functions.R` | Shared utilities: eGFR calculation, model runner, weight selector, label lookups | No |
| `01_consolidate_results.R` | Merge 3 screening rounds into master table, global FDR, novelty flags | No |
| `01b_recover_novelty_screen.R` | Recovery script for novelty round screening results | No |
| `02_validation_data_prep.R` | Download and prepare NHANES 2015-2016 for cross-cycle validation | **Yes** |
| `03_cross_cycle_validation.R` | Replicate FDR-significant findings in 2015-2016 | No |
| `04_dose_response.R` | Quartile analysis and linear trend tests for validated findings | No |
| `05_sensitivity_analyses.R` | Robustness checks: sex/age strata, outlier exclusion, education, cotinine, adults-only | No |
| `06_figures_tables.R` | Volcano plot, forest plot, dose-response curves, heatmap, Table 1-3 | No |
| `07_novelty_assessment.R` | Structured PubMed literature search and novelty classification | **Yes** |
| `09_additional_analyses.R` | Fish/PA/alcohol/creatinine sensitivity, improved volcano plot, power analysis, expanded PubMed search | **Yes** |
| `10_extended_sensitivity.R` | 24h dietary fish adjustment, quadratic age, standardized volcano plots, DAG, systematic novelty search | **Yes** |

### DMA--Uric Acid Enhanced Analysis (Scripts 13--15)

These scripts implement the multi-cycle analysis for the DMA--uric acid focused paper, extending the original 2-cycle finding (Script 11) to 7 independent NHANES cycles (2005--2018) with metal-mixture analysis, formal mediation, and comprehensive sensitivity testing.

| Script | Purpose | Internet required |
|--------|---------|:-:|
| `13_dma_multicycle_extraction.R` | Download and harmonise 7 NHANES cycles (2005--2018): DMA, uric acid, covariates, metals panel, GGT, arsenobetaine, metabolic syndrome variables | **Yes** (first run) |
| `14_dma_enhanced_analyses.R` | Cycle-specific meta-analysis, eGFR + GGT mediation, PAF, interaction tests, arsenobetaine adjustment, temporal trend, sensitivity analyses | No |
| `15_dma_mixture_analysis.R` | Triple metal-mixture analysis (WQS + qgcomp + BKMR) on 9 urinary metals | No |

Outputs are saved to `dma_uric_acid_outputs/`.

### Input Data (CSV)

These three files are outputs of the original interactive screening phase (run in R console, not captured in scripts). They are included as inputs required by `01_consolidate_results.R`:

| File | Description |
|------|-------------|
| `exwas_results_pfas_thyroid.csv` | Round 1: PFAS-thyroid screening (60 tests, 0 FDR-sig) |
| `exwas_broad_screen_all_results.csv` | Round 2: Broad screen (648 tests, 0 FDR-sig) |
| `exwas_expanded_all_results.csv` | Round 3: Expanded screen (1,920 tests, 27 FDR-sig) |

### Output Data (CSV, regenerable from scripts)

| File | Created by |
|------|------------|
| `exwas_master_results.csv` | `01_consolidate_results.R` |
| `exwas_priority_findings.csv` | `01_consolidate_results.R` |
| `exwas_novelty_screen_results.csv` | `01b_recover_novelty_screen.R` |
| `validation_results.csv` | `03_cross_cycle_validation.R` |
| `dose_response_summary.csv` | `04_dose_response.R` |
| `sensitivity_results.csv` | `05_sensitivity_analyses.R` |
| `sensitivity_summary.csv` | `05_sensitivity_analyses.R` |
| `novelty_assessment.csv` | `07_novelty_assessment.R` |
| `additional_sensitivity_results.csv` | `09_additional_analyses.R` |
| `power_analysis.csv` | `09_additional_analyses.R` |
| `figures/table_s12_fish_24h_sensitivity.csv` | `10_extended_sensitivity.R` |
| `figures/table_s13_quadratic_age.csv` | `10_extended_sensitivity.R` |
| `figures/table_s14_systematic_novelty_search.csv` | `10_extended_sensitivity.R` |

### Figures and Tables (`figures/`)

| File | Description |
|------|-------------|
| `fig1_volcano.*` | Volcano plot of all 2,796 tests (top 10 labeled) |
| `fig2_forest.*` | Forest plot of 15 validated findings |
| `fig3_dose_response.*` | Dose-response curves (quartile means) |
| `fig4_heatmap.pdf` | Sensitivity analysis heatmap |
| `fig5_sensitivity_*.png` | Individual sensitivity forest plots (15 panels) |
| `table1_demographics.csv` | Weighted sample demographics |
| `table2_significant_findings.csv` | All FDR-significant findings |
| `table3_validation.csv` | Cross-cycle validation results |
| `fig_s16_volcano_standardized.png` | Volcano plot with standardized effect sizes |
| `fig_s17_volcano_partial_r2.png` | Volcano plot with partial R² |
| `supplementary_dag.txt` | DAG for covariate selection rationale |

## Reproducibility

### Prerequisites

- **R** >= 4.0
- **R packages**: `nhanesA`, `survey`, `tidyverse`, `broom`, `knitr`, `rmarkdown`
- **Internet access** for scripts that download NHANES data (02, 07)

Install packages:

```r
install.packages(c("nhanesA", "survey", "tidyverse", "broom", "knitr", "rmarkdown"))
```

### Execution

Run scripts sequentially in R from the project root directory:

```r
# Step 1: Load shared functions
source("00_functions.R")

# Step 2: Consolidate screening results into master table
source("01_consolidate_results.R")     # ~1 min, offline
source("01b_recover_novelty_screen.R") # ~1 min, offline

# Step 3: Download NHANES 2015-2016 validation data
source("02_validation_data_prep.R")    # ~10-20 min, requires internet

# Step 4: Cross-cycle validation
source("03_cross_cycle_validation.R")  # ~2 min

# Step 5: Dose-response and sensitivity analyses
source("04_dose_response.R")           # ~5 min
source("05_sensitivity_analyses.R")    # ~10 min

# Step 6: Generate figures and tables
source("06_figures_tables.R")          # ~2 min

# Step 7: Novelty assessment (PubMed search)
source("07_novelty_assessment.R")      # ~5 min, requires internet

# Step 8: Additional sensitivity analyses
source("09_additional_analyses.R")     # ~10 min, requires internet

# Step 9: Extended sensitivity (24h fish, quadratic age, standardized plots, DAG)
source("10_extended_sensitivity.R")    # ~10 min, requires internet
```

### Note on the initial screening phase

The original ExWAS screening (Rounds 1-3) was conducted interactively in the R console using NHANES 2017-2018 data downloaded via the `nhanesA` package. The console history (`.Rhistory`, not tracked) documents this work, but the three input CSVs are provided so all downstream analyses are fully reproducible from scripts.

## Data Availability

All data come from the [National Health and Nutrition Examination Survey (NHANES)](https://www.cdc.gov/nchs/nhanes/index.htm), which is publicly available. Scripts 02 and 07 download data directly from the NHANES API via the `nhanesA` R package. No restricted-use data are required.

## Key Results

| Stage | Count |
|-------|-------|
| Total exposure-outcome tests | 2,796 |
| Global FDR < 0.05 | 26 |
| Cross-cycle validated (2015-2016) | 15 / 21 testable (71%) |
| Dose-response significant | 15 / 15 (100%) |
| Robust across sensitivity analyses | 15 / 15 (100%) |
| HIGH novelty | 3 |
| MODERATE novelty | 5 |

## Statistical Methods

- **Design**: Survey-weighted linear regression (`survey::svyglm`) respecting NHANES complex sampling
- **Covariates**: Age, sex, race/ethnicity (3-level), poverty-income ratio, BMI, smoking status
- **Multiple testing**: Benjamini-Hochberg FDR applied globally across all 2,796 tests
- **Weights**: `WTMEC2YR` (blood), `WTSA2YR` (urinary subsample), `WTSSBJ2Y` (surplus serum)

## Related Projects

Two HIGH-novelty findings from this ExWAS were developed into standalone focused papers with deeper mechanistic investigation:

- **Project 29 — DMA and Uric Acid**: Multi-cycle replication (7 NHANES cycles, 2005--2018), metal-mixture analysis, and GGT mediation. Enhanced analysis scripts 13--15 in this repository.
- **Project 30 — Perchlorate and BUN**: Perchlorate association with blood urea nitrogen. Repository: [hayden-farquhar/perchlorate-BUN](https://github.com/hayden-farquhar/perchlorate-BUN) (if available)

Both spin-off projects use the same NHANES 2017-2018 and 2015-2016 data sources as this parent ExWAS analysis.

## Citation

> Farquhar H. Environment-Wide Association Study of Chemical Biomarkers and Health Outcomes in NHANES 2017-2018: Discovery, Validation, and Dose-Response Analysis. *Preprint*. 2026. https://doi.org/10.64898/2026.02.07.26345792. Manuscript under consideration at a peer-reviewed journal.

## License

This project is licensed under the MIT License. See [LICENSE](LICENSE) for details.
