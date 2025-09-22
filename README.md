# Drug-Specific EEG Biomarkers for Anti-Seizure Medication Selection in Focal Epilepsy: A Bayesian Precision Medicine Approach.

R code for Bayesian Analysis

## Citation
[to be updated]

## Requirements
- R version 4.3.2+
- Packages: rstanarm, bayestestR, car, pROC, caret, dplyr

## Usage
1. Install dependencies
2. Load your data (format described in example_data/)
3. Run scripts in order: 01 → 02 → 03 → 04 → 05

## Contact
lorenzo.ricci@ucalgary.ca

# Note on Data Preprocessing
Raw EEG and clinical data preprocessing scripts are not included due to:
- Dataset-specific file paths and variable extraction procedures
- Privacy considerations

This repository focuses on the core statistical and machine learning 
methodology that can be adapted to similar datasets.

# Repository Structure

R/power_analysis.R - Statistical power calculations - 01
R/multicollinearity_vif.R - VIF analysis for variable selection - 02
R/cross_validation.R - 5-fold CV Bayesian modeling - 03
R/statistical_comparisons.R - DeLong tests for AUC comparisons - 04
R/sensitivity_analysis.R - Bayesian prior sensitivity testing - 05


# Usage

Install required packages
Load your data 
Adapt variable names in scripts to match your dataset
Run analysis scripts as needed


### Required Packages
```r
install.packages(c(
  "rstanarm",      # Bayesian modeling
  "bayestestR",    # Bayesian inference
  "car",           # VIF calculations
  "pROC",          # ROC analysis
  "caret",         # Cross-validation
  "dplyr",         # Data manipulation
  "knitr",         # Table formatting
  "writexl"        # Excel export
))


