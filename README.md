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

## 🧩 Script Execution Order

Follow this exact order to reproduce the full analysis workflow.

1️⃣ **cross_validation.R**  
Goal: Perform 5-Fold Cross-Validation using Bayesian logistic regression (`stan_glm`) for EEG-only, Clinical-only, and Mixed (EEG + Clinical) models.  

2️⃣ **statistical_comparison.R**  
Goal: Statistically compare model AUCs using the DeLong test (paired/unpaired).  

3️⃣ **feature_importance_extraction.R**  
Goal: Identify and rank the most significant EEG and clinical predictors using Bayesian posterior distributions.  

4️⃣ **sensitivity_analysis.R**  
Goal: Evaluate model robustness across parameter perturbations or data subsets.  

5️⃣ **power_analysis.R (optional)**  
Goal: Estimate minimum required sample size for detecting significant AUC differences.  

