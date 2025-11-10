# Drug-Specific EEG Biomarkers for Anti-Seizure Medication Selection in Focal Epilepsy  
### A Bayesian Precision Medicine Approach


---

## 🧩 Overview
This repository contains R code used for the **Bayesian statistical analyses** and **EEG-based biomarker discovery** described in  
**_Ricci et al., "Drug-Specific EEG Biomarkers for Anti-Seizure Medication Selection in Focal Epilepsy: A Bayesian Precision Medicine Approach"_**.

It provides a fully reproducible workflow for developing and validating predictive models of anti-seizure medication (ASM) response using EEG and clinical features.

---

## 📖 Citation
> **Ricci L.**, *et al.* (2025).  
> *Drug-Specific EEG Biomarkers for Anti-Seizure Medication Selection in Focal Epilepsy: A Bayesian Precision Medicine Approach.*  
> _[Journal name / DOI to be updated]_  

---

## 📂 Input Data Format

Prepare a single Excel file named, for example:

`database.xlsx`

Each **row** = one patient  
Each **column** = one variable  

| OUTCOME_bin | SEX | AGE | AETIOLOGY | EEG | PSYCHIATRIC | PSD_DELTA_LOG_z | PSD_THETA_LOG_z | AEC_BETA_z | Farmaco |
|--------------|-----|-----|------------|-----|--------------|-----------------|-----------------|-------------|----------|
| 1 | F | 42 | STR | yes | no | 0.53 | 0.21 | -0.34 | LEV |
| 0 | M | 55 | UNK | no | yes | -0.27 | 0.10 | 0.14 | NACB |

### Required columns

- **OUTCOME_bin**: binary outcome (`0 = responder`, `1 = non-responder`)
- **EEG variables**: all numerical EEG-derived metrics
- **Clinical variables**: categorical/demographic predictors
- **Optional:** `Farmaco` (drug classification: `LEV`, `NACB`, etc.)

---

## 🧩 Execution Order

Follow this order to reproduce the complete analysis pipeline.

### 1️⃣ cross_validation.R
**Goal:** Perform 5-fold cross-validation using Bayesian logistic regression (`stan_glm`) for EEG-only and Clinical-only models.

**Input:**
```r
eeg_features_to_use <- c("PSD_DELTA_LOG_z", "PSD_THETA_LOG_z", "AEC_BETA_z", ...)
clinical_features_to_use <- c("SEX", "AGE", "AETIOLOGY", "EEG", "PSYCHIATRIC")

**Output:**
`EEG_vs_Clinical_Comparisons.xlsx` — tables of AUC differences and p-values  
Console output includes significant results (*, **, ***)

**Description:**  
Performs all relevant statistical tests:

- EEG vs Clinical models  
- LEV vs NACB EEG models  
- ALL vs Drug-specific EEG models  

Also generates structured Excel outputs and clinical interpretations.

---

### 3️⃣ feature_importance_extraction.R
**Goal:** Identify and rank the most significant EEG and clinical predictors using Bayesian posterior distributions.

**Input:**  
List of `stan_glm` model objects or `cv_results_list`.

**Output:**  
`Feature_Importance_Summary.xlsx`  
Includes:  
- Variable names  
- Posterior mean ± SD  
- 95% credible intervals  
- Rank order (importance)

**Description:**  
Extracts posterior summaries from each fitted model, computes absolute effect sizes, and ranks variables by contribution to model discrimination (AUC).

---

### 4️⃣ sensitivity_analysis.R
**Goal:** Evaluate model robustness across parameter perturbations or data subsets.

**Input:**  
Outputs from `cross_validation.R`.

**Output:**  
Summary tables showing variation in AUC and metrics under simulated noise or subgroup conditions.

---

### 5️⃣ power_analysis.R (optional)
**Goal:** Estimate minimum required sample size for detecting significant AUC differences.

**Input:**  
Mean AUCs and standard deviations from CV results.

**Output:**  
Power estimation tables and suggested N for replication studies.

---

## 🧪 Recommended Execution Flow

1. Prepare your Excel data (`database.xlsx`) following the specified format.  
2. Run `cross_validation.R` to fit and validate models.  
3. Run `statistical_comparison.R` to compare model AUCs.  
4. Run `feature_importance_extraction.R` to find top predictive features.  
5. (Optional) Run `sensitivity_analysis.R` or `power_analysis.R` for robustness checks.

---

## 📊 Output Summary

| File | Description |
|------|--------------|
| **PERSEFO_5Fold_CV_Results.xlsx** | Performance metrics for all models |
| **EEG_vs_Clinical_Comparisons.xlsx** | Statistical significance (DeLong test results) |
| **Feature_Importance_Summary.xlsx** | Ranked variable importance across models |
| **PERSEFO_Complete_Results.xlsx** | Aggregated metrics + confidence intervals (final summary) |

---

## 🚀 Quick Start Example

```r
# Load required packages
library(dplyr)
library(rstanarm)
library(pROC)
library(caret)
library(writexl)

# 1. Load your dataset
database <- readxl::read_excel("database.xlsx")

# 2. Run scripts in sequence
source("R/cross_validation.R")
source("R/statistical_comparison.R")
source("R/feature_importance_extraction.R")

# 3. Execute main functions
cv_results_list <- run_separate_models_cv(database, eeg_features_to_use, clinical_features_to_use)
comparison_results <- run_eeg_vs_clinical_comparisons(cv_results_list)
importance_summary <- extract_feature_importance(cv_results_list)

