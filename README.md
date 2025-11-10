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

## ⚙️ Requirements
- **R version:** 4.3.2 or higher  
- **Required packages:**  
  `rstanarm`, `bayestestR`, `car`, `pROC`, `caret`, `dplyr`, `knitr`, `writexl`

You can install all dependencies with:
```r
install.packages(c(
  "rstanarm",      # Bayesian modeling
  "bayestestR",    # Bayesian inference
  "car",           # Multicollinearity diagnostics
  "pROC",          # ROC and AUC analysis
  "caret",         # Cross-validation
  "dplyr",         # Data manipulation
  "knitr",         # Table formatting
  "writexl"        # Excel export
))



