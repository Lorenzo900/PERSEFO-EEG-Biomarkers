# =================================================================
# PERSEFO: MULTICOLLINEARITY ASSESSMENT AND VIF ANALYSIS
# Robust VIF checking for EEG and clinical variables separately
# =================================================================

# Load required libraries
suppressPackageStartupMessages({
  library(dplyr)
  library(car)  # For VIF calculation
  library(writexl)
  library(knitr)
})

# Set seed for reproducibility
set.seed(123)

# =================================================================
# FUNCTION: VIF-BASED MULTICOLLINEARITY CHECK
# =================================================================

check_and_remove_vif <- function(data, variables, outcome_var = "OUTCOME_bin", vif_threshold = 10) {
  
  cat("\n=== VIF MULTICOLLINEARITY CHECK ===\n")
  cat(sprintf("VIF Threshold: %.1f\n", vif_threshold))
  cat(sprintf("Initial variables: %d\n", length(variables)))
  
  # Remove cases with missing data
  data_complete <- data[complete.cases(data[c(outcome_var, variables)]), ]
  
  if(nrow(data_complete) < 50) {
    cat("WARNING: Too few cases for VIF check\n")
    return(list(remaining_variables = variables, removed_variables = c(), n_removed = 0))
  }
  
  # STEP 1: Remove zero/near-zero variance variables
  cat("\n--- STEP 1: Zero variance check ---\n")
  
  # Identify numeric and factor variables
  numeric_vars <- variables[sapply(variables, function(v) is.numeric(data_complete[[v]]))]
  factor_vars <- setdiff(variables, numeric_vars)
  
  # Check variance for numeric variables
  low_variance_vars <- c()
  if(length(numeric_vars) > 0) {
    variances <- sapply(numeric_vars, function(v) var(data_complete[[v]], na.rm = TRUE))
    low_variance_vars <- numeric_vars[variances < 1e-8]
    
    if(length(low_variance_vars) > 0) {
      cat(sprintf("Removed %d near-zero variance variables:\n", length(low_variance_vars)))
      for(v in low_variance_vars) {
        cat(sprintf("  - %s (var = %e)\n", v, variances[v]))
      }
    }
  }
  
  # Check categorical variables with single level
  single_level_vars <- c()
  if(length(factor_vars) > 0) {
    for(v in factor_vars) {
      n_levels <- length(unique(data_complete[[v]]))
      if(n_levels <= 1) {
        single_level_vars <- c(single_level_vars, v)
      }
    }
    
    if(length(single_level_vars) > 0) {
      cat(sprintf("Removed %d single-level categorical variables:\n", length(single_level_vars)))
      for(v in single_level_vars) {
        cat(sprintf("  - %s\n", v))
      }
    }
  }
  
  # Remove problematic variables
  problematic_vars <- c(low_variance_vars, single_level_vars)
  remaining_vars <- setdiff(variables, problematic_vars)
  
  cat(sprintf("Variables after initial cleaning: %d\n", length(remaining_vars)))
  
  # STEP 2: Check perfect correlations between numeric variables
  cat("\n--- STEP 2: Perfect correlation check ---\n")
  
  current_numeric_vars <- remaining_vars[sapply(remaining_vars, function(v) is.numeric(data_complete[[v]]))]
  perfect_corr_vars <- c()
  
  if(length(current_numeric_vars) > 1) {
    # Calculate correlation matrix
    cor_matrix <- cor(data_complete[current_numeric_vars], use = "complete.obs")
    
    # Find perfect correlations (|r| > 0.999) excluding diagonal
    for(i in 1:(ncol(cor_matrix)-1)) {
      for(j in (i+1):ncol(cor_matrix)) {
        if(abs(cor_matrix[i,j]) > 0.999) {
          var_to_remove <- colnames(cor_matrix)[j]  # Remove second variable
          if(!var_to_remove %in% perfect_corr_vars) {
            perfect_corr_vars <- c(perfect_corr_vars, var_to_remove)
            cat(sprintf("Perfect correlation: %s <-> %s (r = %.3f) - removing %s\n", 
                        colnames(cor_matrix)[i], colnames(cor_matrix)[j], 
                        cor_matrix[i,j], var_to_remove))
          }
        }
      }
    }
  }
  
  # Remove perfectly correlated variables
  remaining_vars <- setdiff(remaining_vars, perfect_corr_vars)
  cat(sprintf("Variables after perfect correlation removal: %d\n", length(remaining_vars)))
  
  # STEP 3: Iterative VIF check
  cat("\n--- STEP 3: Iterative VIF check ---\n")
  
  removed_vif_vars <- c()
  iteration <- 1
  
  while(length(remaining_vars) >= 2 && iteration <= 50) {
    
    cat(sprintf("\nVIF iteration %d: %d variables\n", iteration, length(remaining_vars)))
    
    # Create formula
    formula_str <- paste(outcome_var, "~", paste(remaining_vars, collapse = " + "))
    
    tryCatch({
      # Fit model
      temp_model <- glm(as.formula(formula_str), 
                        data = data_complete, 
                        family = binomial,
                        control = glm.control(maxit = 100))
      
      # Check for aliased coefficients
      if(any(is.na(coef(temp_model)))) {
        # Identify aliased coefficients
        coefs <- coef(temp_model)
        na_coefs <- names(coefs)[is.na(coefs)]
        
        cat("Aliased coefficients found:\n")
        for(coef_name in na_coefs) {
          cat(sprintf("  - %s\n", coef_name))
        }
        
        # Extract variable names from aliased coefficients
        aliased_vars <- c()
        for(var in remaining_vars) {
          if(any(grepl(var, na_coefs))) {
            aliased_vars <- c(aliased_vars, var)
          }
        }
        
        if(length(aliased_vars) > 0) {
          # Remove first aliased variable
          var_to_remove <- aliased_vars[1]
          remaining_vars <- setdiff(remaining_vars, var_to_remove)
          removed_vif_vars <- c(removed_vif_vars, var_to_remove)
          cat(sprintf("Removed aliased variable: %s\n", var_to_remove))
          next
        }
      }
      
      # If model converged without aliased coefficients, calculate VIF
      if(temp_model$converged && !any(is.na(coef(temp_model)))) {
        
        # Calculate VIF
        vif_values <- car::vif(temp_model)
        
        # Handle VIF output for categorical variables
        if(is.matrix(vif_values)) {
          # For categorical variables, use GVIF^(1/(2*Df))
          vif_values <- vif_values[, "GVIF^(1/(2*Df))"]
        }
        
        # Find maximum VIF
        max_vif <- max(vif_values, na.rm = TRUE)
        max_vif_var <- names(vif_values)[which.max(vif_values)]
        
        cat(sprintf("Maximum VIF: %.2f (%s)\n", max_vif, max_vif_var))
        
        # If all VIFs are below threshold, exit
        if(max_vif < vif_threshold) {
          cat("All VIFs below threshold - COMPLETED\n")
          break
        }
        
        # Remove variable with highest VIF
        remaining_vars <- setdiff(remaining_vars, max_vif_var)
        removed_vif_vars <- c(removed_vif_vars, max_vif_var)
        cat(sprintf("Removed high VIF: %s (VIF = %.2f)\n", max_vif_var, max_vif))
        
      } else {
        cat("Model did not converge - stopping VIF check\n")
        break
      }
      
    }, error = function(e) {
      cat(sprintf("Error in VIF calculation: %s\n", e$message))
      # In case of error, remove last variable added
      if(length(remaining_vars) > 0) {
        var_to_remove <- remaining_vars[length(remaining_vars)]
        remaining_vars <<- setdiff(remaining_vars, var_to_remove)
        removed_vif_vars <<- c(removed_vif_vars, var_to_remove)
        cat(sprintf("Removed due to error: %s\n", var_to_remove))
      } else {
        break
      }
    })
    
    iteration <- iteration + 1
  }
  
  # FINAL SUMMARY
  total_removed <- c(problematic_vars, perfect_corr_vars, removed_vif_vars)
  
  cat(sprintf("\n=== FINAL VIF CHECK RESULTS ===\n"))
  cat(sprintf("Initial variables: %d\n", length(variables)))
  cat(sprintf("Final variables: %d\n", length(remaining_vars)))
  cat(sprintf("Total removed variables: %d\n", length(total_removed)))
  
  if(length(total_removed) > 0) {
    cat("\nRemoval breakdown:\n")
    if(length(problematic_vars) > 0) {
      cat(sprintf("- Zero variance/single level: %d\n", length(problematic_vars)))
    }
    if(length(perfect_corr_vars) > 0) {
      cat(sprintf("- Perfect correlations: %d\n", length(perfect_corr_vars)))
    }
    if(length(removed_vif_vars) > 0) {
      cat(sprintf("- High VIF: %d\n", length(removed_vif_vars)))
    }
  }
  
  return(list(
    remaining_variables = remaining_vars,
    removed_variables = total_removed,
    n_removed = length(total_removed),
    removed_zero_var = problematic_vars,
    removed_perfect_corr = perfect_corr_vars,
    removed_high_vif = removed_vif_vars
  ))
}

# =================================================================
# FUNCTION: VIF ANALYSIS FOR EEG AND CLINICAL SEPARATELY
# =================================================================

run_vif_analysis_separate <- function(database, eeg_variables, clinical_variables) {
  
  cat("=== VIF ANALYSIS: EEG AND CLINICAL SEPARATELY ===\n")
  
  # Results storage
  vif_results <- list()
  
  # 1. Clinical variables only (ALL dataset)
  cat("\n--- Dataset ALL - Clinical variables only ---\n")
  vif_clinical <- check_and_remove_vif(database, clinical_variables, "OUTCOME_bin")
  vif_results[["clinical_all"]] <- vif_clinical
  
  # 2. EEG variables only (ALL dataset)
  cat("\n--- Dataset ALL - EEG variables only ---\n")
  vif_all_eeg <- check_and_remove_vif(database, eeg_variables, "OUTCOME_bin")
  vif_results[["eeg_all"]] <- vif_all_eeg
  
  # 3. Drug-specific datasets (if applicable)
  if("Farmaco" %in% names(database)) {
    
    # LEV dataset - EEG only
    cat("\n--- Dataset LEV - EEG variables only ---\n")
    lev_data <- database %>% filter(Farmaco == "LEV")
    if(nrow(lev_data) > 50) {  # Minimum sample size check
      vif_lev_eeg <- check_and_remove_vif(lev_data, eeg_variables, "OUTCOME_bin")
      vif_results[["eeg_lev"]] <- vif_lev_eeg
    } else {
      cat("LEV dataset too small for VIF analysis\n")
    }
    
    # NACB dataset - EEG only
    cat("\n--- Dataset NACB - EEG variables only ---\n")
    nacb_data <- database %>% filter(Farmaco == "NACB")
    if(nrow(nacb_data) > 50) {  # Minimum sample size check
      vif_nacb_eeg <- check_and_remove_vif(nacb_data, eeg_variables, "OUTCOME_bin")
      vif_results[["eeg_nacb"]] <- vif_nacb_eeg
    } else {
      cat("NACB dataset too small for VIF analysis\n")
    }
  }
  
  # Create summary table
  vif_summary <- data.frame(
    Dataset = names(vif_results),
    Variable_Type = c("Clinical", "EEG", 
                      if("eeg_lev" %in% names(vif_results)) "EEG" else NULL,
                      if("eeg_nacb" %in% names(vif_results)) "EEG" else NULL),
    Original_Vars = sapply(vif_results, function(x) length(x$removed_variables) + length(x$remaining_variables)),
    Final_Vars = sapply(vif_results, function(x) length(x$remaining_variables)),
    Removed_Vars = sapply(vif_results, function(x) x$n_removed),
    Retention_Rate = round(sapply(vif_results, function(x) 
      length(x$remaining_variables) / (length(x$removed_variables) + length(x$remaining_variables))), 3),
    stringsAsFactors = FALSE
  )
  
  # Clean up dataset names for display
  vif_summary$Dataset <- case_when(
    vif_summary$Dataset == "clinical_all" ~ "ALL - Clinical",
    vif_summary$Dataset == "eeg_all" ~ "ALL - EEG", 
    vif_summary$Dataset == "eeg_lev" ~ "LEV - EEG",
    vif_summary$Dataset == "eeg_nacb" ~ "NACB - EEG",
    TRUE ~ vif_summary$Dataset
  )
  
  print(knitr::kable(vif_summary, caption = "VIF Analysis Summary: EEG and Clinical Variables"))
  
  # Save results
  write_xlsx(list(
    VIF_Summary = vif_summary,
    Clinical_Results = if("clinical_all" %in% names(vif_results)) vif_results[["clinical_all"]] else NULL,
    EEG_ALL_Results = if("eeg_all" %in% names(vif_results)) vif_results[["eeg_all"]] else NULL,
    EEG_LEV_Results = if("eeg_lev" %in% names(vif_results)) vif_results[["eeg_lev"]] else NULL,
    EEG_NACB_Results = if("eeg_nacb" %in% names(vif_results)) vif_results[["eeg_nacb"]] else NULL
  ), "VIF_Analysis_Separate_Variables.xlsx")
  
  cat("\n✅ VIF analysis results saved to 'VIF_Analysis_Separate_Variables.xlsx'\n")
  
  return(vif_results)
}

# =================================================================
# EXAMPLE USAGE
# =================================================================

# Example variable definitions - replace with your actual variables
example_eeg_variables <- c(
  "PSD_DELTA_LOG_z", "PSD_THETA_LOG_z", "PSD_ALFA_LOG_z", "PSD_BETA_LOG_z", "PSD_GAMMA_LOG_z",
  "AEC_DELTA_z", "AEC_THETA_z", "AEC_ALFA_z", "AEC_BETA_z", "AEC_GAMMA_z"
  # Add your actual EEG variables here
)

example_clinical_variables <- c("SEX", "SEIZURE_TYPE", "AETIOLOGY", "EEG", "PSYCHIATRIC", "AGE")

# To run the analysis:
# vif_results <- run_vif_analysis_separate(your_database, your_eeg_variables, your_clinical_variables)

cat("VIF analysis functions loaded successfully!\n")
cat("Use run_vif_analysis_separate() to analyze EEG and clinical variables separately.\n")
