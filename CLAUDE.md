# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project

Environment-wide association study (ExWAS) screening 92 chemical biomarkers against 48 health outcomes in NHANES 2017-2018, with cross-cycle validation in 2015-2016. R-based, using `survey` for complex sampling design and `nhanesA` for data access. GitHub: https://github.com/hayden-farquhar/NHANES-EXWAS

## Running the Pipeline

Scripts run sequentially in R from the project root. All use relative paths and depend on the `.Rproj` setting the working directory.

```r
source("00_functions.R")              # Must run first — loads packages, defines utilities
source("01_consolidate_results.R")    # Merges 3 screening CSVs → master table
source("01b_recover_novelty_screen.R")# Re-runs novelty chemical screening (needs internet first time)
source("02_validation_data_prep.R")   # Downloads NHANES 2015-2016 (internet, ~15 min)
source("03_cross_cycle_validation.R") # Replicates FDR-sig findings in 2015-2016
source("04_dose_response.R")          # Quartile dose-response for validated findings
source("05_sensitivity_analyses.R")   # 9 robustness specifications per finding
source("06_figures_tables.R")         # All figures + Tables 1-3
source("07_novelty_assessment.R")     # PubMed literature search + novelty classification
source("09_additional_analyses.R")    # All sensitivity: fish/PA/alcohol/creatinine/protein, 6-level race, within-round FDR, LOD, STROBE
```

Required packages: `nhanesA`, `survey`, `tidyverse`, `broom`, `ggrepel`, `forestplot`, `pheatmap`, `kableExtra`, `tableone`, `knitr`, `rmarkdown`. Script `00_functions.R` auto-installs missing packages.

## Architecture

### Data Flow

```
3 input CSVs (interactive screening results, committed)
    ↓
01 + 01b → exwas_master_results.RData (master, priority_findings)
    ↓
02 → nhanes_2015_2016_validation.RData (val_data, downloaded from CDC)
    ↓
03 → validation_results.RData (validated)
    ↓
04, 05, 07 (parallel) → dose_response, sensitivity, novelty .RData
    ↓
06 → figures/ (PDFs, PNGs, table CSVs)
    ↓
09 → all sensitivity analyses, volcano, power, within-round FDR, LOD analysis, STROBE
```

### Key Components in `00_functions.R`

- `run_exwas_model()` — Core model runner. Takes exposure/outcome names, builds survey design, runs `svyglm()`, returns tidy one-row result. Auto-drops BMI from covariates when outcome is BMI/waist.
- `select_weight()` — Picks survey weight by variable prefix: `LBX` → WTMEC2YR, `URX` → WTSA2YR, `SS` → WTSSBJ2Y.
- `calc_eGFR()` — Race-free CKD-EPI 2021 equation.
- `CHEM_LABELS` / `OUTCOME_LABELS` — Named character vectors mapping variable names to readable labels.
- `COVARS_FULL` / `COVARS_NO_BMI` — Covariate formula strings.

### Statistical Model

```
outcome ~ log_exposure + RIDAGEYR + female + race3 + INDFMPIR + BMXBMI + smoker
```

Survey design: `svydesign(id=~SDMVPSU, strata=~SDMVSTRA, weights=~weight_var, nest=TRUE)`

FDR: Benjamini-Hochberg applied **globally** across all 2,796 tests in `01_consolidate_results.R`.

### Input CSVs (not regenerable from scripts)

The three screening CSVs (`exwas_results_pfas_thyroid.csv`, `exwas_broad_screen_all_results.csv`, `exwas_expanded_all_results.csv`) were produced interactively and are committed as pipeline inputs. Everything downstream is regenerable.

## Critical Technical Notes

**Data files:**
- `.RDataTmp` is **CORRUPTED** — never load it. Use `exwas_novelty_screen.RData` (has `dat`: 5265x609) or `exwas_expanded_results.RData` (has `expanded`: 5265x502).
- `.RData` files are gitignored. CSVs are tracked and are the canonical interchange format.

**nhanesA quirks (especially 2015-2016 cycle):**
- Columns often returned as **text factors**, not numeric codes. Always convert: `as.numeric(as.character(x))`.
- RIDRETH3 race codes: match on **both** `"3"` and `"Non-Hispanic White"` (nhanesA inconsistently returns codes vs labels).
- PHQ-9 depression items come as text: `"Not at all"`=0, `"Several days"`=1, `"More than half the days"`=2, `"Nearly every day"`=3.
- Table names differ across cycles (e.g., `UHM_I` → `UM_I`, `UVOC2_I` may need `UVOC_I` fallback). Scripts handle this with try/fallback logic.
- `COT_I`/`SMQ_I` can have transient download failures — retry logic is built in.

**Survey design:**
- Variables referenced in `svydesign()` (SDMVPSU, SDMVSTRA, weight) must exist in the data **before** the call.
- Low residual df (7-9) with surplus/subsample weights is **expected** for NHANES, not a bug.
- Binary outcomes use `quasibinomial()` family to avoid separation issues.

**Filtering:**
- Chemicals with >70% values at LOD are excluded (insufficient variability).
- Implausibly large effects filtered: |beta| > 50 for continuous outcomes, >5 for binary.
- Below-LOD values replaced with LOD/sqrt(2) before log-transformation.

**R quirks:**
- `print(n=50)` only works on tibbles, not data.frames — wrap with `as_tibble()` first.

## Results Summary

2,796 tests → 26 FDR-significant → 15/21 cross-cycle validated (71%) → 14 robust across all sensitivity analyses (iodine-BMI eliminated as dilution artifact; methylmercury-waist suspected dietary confounder) → 2 confirmed HIGH novelty (DMA-uric acid, perchlorate-BUN), 5 MODERATE novelty.
