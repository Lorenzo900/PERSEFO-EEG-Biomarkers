# =================================================================
# PERSEFO: 5-FOLD CROSS-VALIDATION FRAMEWORK
# Robust cross-validation for EEG and clinical models separately
# =================================================================

# Load required libraries
suppressPackageStartupMessages({
  library(dplyr)
  library(rstanarm)
  library(pROC)
  library(caret)
  library(knitr)
  library(writexl)
})

# Set seed for reproducibility
set.seed(123)

# =================================================================
# 5-FOLD CROSS-VALIDATION FUNCTION
# =================================================================

run_5fold_cv <- function(data, variables, model_name, outcome_var = "OUTCOME_bin", n_folds = 5) {
  
  cat("\n", paste(rep("=", 60), collapse = ""), "\n")
  cat("Running 5-Fold Cross-Validation for:", model_name, "\n")
  
  # Remove cases with missing data
  data_complete <- data[complete.cases(data[c(outcome_var, variables)]), ]
  
  # Check minimum sample size
  if(nrow(data_complete) < 50) {
    cat("WARNING: Insufficient patients for 5-fold CV. Found:", nrow(data_complete), "\n")
    return(NULL)
  }
  
  # Check outcome variability
  if(length(unique(data_complete[[outcome_var]])) < 2) {
    cat("WARNING: Outcome has only one class. Cannot perform CV.\n")
    return(NULL)
  }
  
  cat(sprintf("  Using %d patients and %d variables for 5-fold validation\n", 
              nrow(data_complete), length(variables)))
  
  # Create stratified folds
  set.seed(123)
  folds <- createFolds(data_complete[[outcome_var]], k = n_folds, list = TRUE, returnTrain = FALSE)
  
  all_predictions <- data.frame()
  formula_str <- paste(outcome_var, "~", paste(variables, collapse = " + "))
  
  # 5-Fold CV loop
  for(fold_i in 1:n_folds) {
    cat(sprintf("  - Fold %d/%d\n", fold_i, n_folds))
    
    # Train/test split
    test_indices <- folds[[fold_i]]
    train_data <- data_complete[-test_indices, ]
    test_data <- data_complete[test_indices, ]
    
    cat(sprintf("    Training: %d patients, Test: %d patients\n", nrow(train_data), nrow(test_data)))
    
    # Check minimum training data
    if(nrow(train_data) < 20 || length(unique(train_data[[outcome_var]])) < 2) {
      cat("    -> Insufficient training data or single outcome class. Skipping fold.\n")
      next
    }
    
    # Train model
    tryCatch({
      model <- stan_glm(
        as.formula(formula_str),
        data = train_data,
        family = binomial(link = "logit"),
        prior = normal(0, 1),
        prior_intercept = normal(0, 2),
        chains = 2, 
        iter = 1000, 
        cores = 2, 
        refresh = 0,
        seed = 123 + fold_i
      )
      
      # Predict probabilities on test fold
      pred_probs <- posterior_predict(model, newdata = test_data, draws = 200)
      pred_mean <- apply(pred_probs, 2, mean)
      
      # Store fold predictions
      fold_predictions <- data.frame(
        prob = pred_mean,
        observed = test_data[[outcome_var]],
        fold = fold_i,
        patient_id = rownames(test_data)
      )
      
      all_predictions <- rbind(all_predictions, fold_predictions)
      
      cat(sprintf("    -> Completed. Predictions: %d patients\n", nrow(fold_predictions)))
      
    }, error = function(e) {
      cat(sprintf("    -> ERROR in fold %d: %s\n", fold_i, e$message))
    })
  }
  
  # Check if any predictions were generated
  if(nrow(all_predictions) == 0) {
    cat("ERROR: No predictions generated from 5-fold CV.\n")
    return(NULL)
  }
  
  # Calculate aggregate performance
  roc_obj <- roc(all_predictions$observed, all_predictions$prob, quiet = TRUE)
  auc_val <- auc(roc_obj)
  
  # Calculate CI with bootstrap
  tryCatch({
    ci_val <- ci.auc(roc_obj, method = "bootstrap", n.boot = 500)
    auc_final <- ci_val[2]
    ci_lower <- ci_val[1]
    ci_upper <- ci_val[3]
  }, error = function(e) {
    auc_final <- auc_val
    ci_lower <- auc_val
    ci_upper <- auc_val
  })
  
  # Calculate per-fold metrics
  fold_aucs <- sapply(1:n_folds, function(f) {
    fold_data <- all_predictions[all_predictions$fold == f, ]
    if(nrow(fold_data) > 5 && length(unique(fold_data$observed)) > 1) {
      auc(roc(fold_data$observed, fold_data$prob, quiet = TRUE))
    } else {
      NA
    }
  })
  
  mean_fold_auc <- mean(fold_aucs, na.rm = TRUE)
  sd_fold_auc <- sd(fold_aucs, na.rm = TRUE)
  
  # Calculate optimal threshold using Youden's J statistic
  optimal_coords <- coords(roc_obj, "best", ret = c("threshold", "sensitivity", "specificity", "accuracy"))
  optimal_threshold <- optimal_coords$threshold
  
  # Calculate additional metrics at optimal threshold
  predicted_class <- ifelse(all_predictions$prob > optimal_threshold, 1, 0)
  confusion_matrix <- table(Predicted = predicted_class, Observed = all_predictions$observed)
  
  if(nrow(confusion_matrix) == 2 && ncol(confusion_matrix) == 2) {
    TP <- confusion_matrix[2,2]
    FP <- confusion_matrix[2,1] 
    FN <- confusion_matrix[1,2]
    TN <- confusion_matrix[1,1]
    
    sensitivity <- TP / (TP + FN)
    specificity <- TN / (TN + FP)
    precision <- TP / (TP + FP)
    balanced_accuracy <- (sensitivity + specificity) / 2
    f1_score <- 2 * (precision * sensitivity) / (precision + sensitivity)
  } else {
    sensitivity <- specificity <- precision <- balanced_accuracy <- f1_score <- NA
  }
  
  cat(sprintf("  ✅ 5-Fold CV for '%s' completed. AUC Pooled: %.3f [95%% CI: %.3f–%.3f]\n",
              model_name, auc_final, ci_lower, ci_upper))
  cat(sprintf("      AUC per fold: %.3f ± %.3f\n", mean_fold_auc, sd_fold_auc))
  
  return(list(
    model_name = model_name,
    n_predictions = nrow(all_predictions),
    n_folds_completed = length(unique(all_predictions$fold)),
    pooled_auc = auc_final,
    ci_lower = ci_lower,
    ci_upper = ci_upper,
    mean_fold_auc = mean_fold_auc,
    sd_fold_auc = sd_fold_auc,
    fold_aucs = fold_aucs,
    optimal_threshold = optimal_threshold,
    sensitivity = sensitivity,
    specificity = specificity,
    precision = precision,
    balanced_accuracy = balanced_accuracy,
    f1_score = f1_score,
    predictions = all_predictions,
    variables_used = variables,
    n_variables = length(variables)
  ))
}

# =================================================================
# SEPARATE MODEL COMPARISON FUNCTION (EEG vs CLINICAL)
# =================================================================

run_separate_models_cv <- function(database, eeg_variables, clinical_variables, drug_stratification = TRUE) {
  
  cat("=== 5-FOLD CV: EEG AND CLINICAL MODELS SEPARATELY ===\n")
  
  # Initialize results storage
  cv_results <- list()
  
  # 1. Clinical-only model (ALL dataset)
  cat("\n=== CLINICAL-ONLY MODEL ===\n")
  
  cat("Running Clinical model...\n")
  result_clinical <- run_5fold_cv(database, clinical_variables, "Clinical_Only")
  cv_results[["Clinical"]] <- result_clinical
  
  # 2. EEG-only models
  cat("\n=== EEG-ONLY MODELS ===\n")
  
  # All patients - EEG only
  cat("Running ALL_EEG model...\n")
  result_all_eeg <- run_5fold_cv(database, eeg_variables, "ALL_EEG_Only")
  cv_results[["ALL_EEG"]] <- result_all_eeg
  
  # Drug-specific EEG models if stratification is enabled
  if(drug_stratification && "Farmaco" %in% names(database)) {
    
    # LEV patients - EEG only
    cat("Running LEV_EEG model...\n")
    lev_data <- database %>% filter(Farmaco == "LEV")
    result_lev_eeg <- run_5fold_cv(lev_data, eeg_variables, "LEV_EEG_Only")
    cv_results[["LEV_EEG"]] <- result_lev_eeg
    
    # NACB patients - EEG only
    cat("Running NACB_EEG model...\n")
    nacb_data <- database %>% filter(Farmaco == "NACB")
    result_nacb_eeg <- run_5fold_cv(nacb_data, eeg_variables, "NACB_EEG_Only")
    cv_results[["NACB_EEG"]] <- result_nacb_eeg
  }
  
  # =================================================================
  # RESULTS SUMMARY
  # =================================================================
  
  cat("\n=== CROSS-VALIDATION RESULTS SUMMARY ===\n")
  
  # Create summary table
  summary_df <- do.call(rbind, lapply(cv_results, function(res) {
    if (!is.null(res)) {
      data.frame(
        Model = res$model_name,
        Model_Type = case_when(
          grepl("Clinical", res$model_name) ~ "Clinical",
          grepl("EEG", res$model_name) ~ "EEG",
          TRUE ~ "Other"
        ),
        Dataset = case_when(
          grepl("ALL", res$model_name) ~ "ALL",
          grepl("LEV", res$model_name) ~ "LEV",
          grepl("NACB", res$model_name) ~ "NACB",
          grepl("Clinical", res$model_name) ~ "ALL",
          TRUE ~ "Unknown"
        ),
        N_Variables = res$n_variables,
        N_Predictions = res$n_predictions,
        N_Folds = res$n_folds_completed,
        Pooled_AUC = round(res$pooled_auc, 3),
        CI_Lower = round(res$ci_lower, 3),
        CI_Upper = round(res$ci_upper, 3),
        AUC_CI = paste0(round(res$pooled_auc, 3), " [", round(res$ci_lower, 3), "-", round(res$ci_upper, 3), "]"),
        Sensitivity = round(res$sensitivity, 3),
        Specificity = round(res$specificity, 3),
        Balanced_Accuracy = round(res$balanced_accuracy, 3),
        Clinical_Category = case_when(
          res$pooled_auc >= 0.8 ~ "Excellent",
          res$pooled_auc >= 0.7 ~ "Good",
          res$pooled_auc >= 0.6 ~ "Fair",
          TRUE ~ "Poor"
        )
      )
    }
  }))
  
  if (nrow(summary_df) > 0) {
    summary_df <- summary_df %>% arrange(desc(Pooled_AUC))
    
    print(knitr::kable(summary_df, caption = "5-Fold CV Results: EEG vs Clinical Models"))
    
    # Save results
    write_xlsx(list(
      Performance_Summary = summary_df,
      Clinical_Results = if(!is.null(cv_results[["Clinical"]])) cv_results[["Clinical"]] else NULL,
      EEG_ALL_Results = if(!is.null(cv_results[["ALL_EEG"]])) cv_results[["ALL_EEG"]] else NULL,
      EEG_LEV_Results = if(!is.null(cv_results[["LEV_EEG"]])) cv_results[["LEV_EEG"]] else NULL,
      EEG_NACB_Results = if(!is.null(cv_results[["NACB_EEG"]])) cv_results[["NACB_EEG"]] else NULL
    ), "Separate_Models_CV_Results.xlsx")
    
    cat("\n✅ Results saved to 'Separate_Models_CV_Results.xlsx'\n")
    
    # Clinical interpretation
    cat("\n=== CLINICAL INTERPRETATION ===\n")
    
    # Best overall model
    best_model <- summary_df[1, ]
    cat(sprintf("🏆 Best performing model: %s (AUC = %.3f)\n", best_model$Model, best_model$Pooled_AUC))
    
    # Compare EEG vs Clinical
    clinical_model <- summary_df[summary_df$Model_Type == "Clinical", ]
    eeg_models <- summary_df[summary_df$Model_Type == "EEG", ]
    
    if(nrow(clinical_model) > 0) {
      cat(sprintf("📊 Clinical model performance: AUC = %.3f (%s)\n", 
                  clinical_model$Pooled_AUC[1], clinical_model$Clinical_Category[1]))
    }
    
    if(nrow(eeg_models) > 0) {
      cat("🧠 EEG model performance:\n")
      for(i in 1:nrow(eeg_models)) {
        cat(sprintf("   - %s: AUC = %.3f (%s)\n", 
                    eeg_models$Dataset[i], eeg_models$Pooled_AUC[i], eeg_models$Clinical_Category[i]))
      }
      
      # Drug-specific comparison
      lev_performance <- eeg_models[eeg_models$Dataset == "LEV", ]
      nacb_performance <- eeg_models[eeg_models$Dataset == "NACB", ]
      
      if(nrow(lev_performance) > 0 && nrow(nacb_performance) > 0) {
        cat("\n💊 Drug-specific EEG biomarker evidence:\n")
        cat(sprintf("   - LEV: AUC = %.3f\n", lev_performance$Pooled_AUC[1]))
        cat(sprintf("   - NACB: AUC = %.3f\n", nacb_performance$Pooled_AUC[1]))
        
        if(nacb_performance$Pooled_AUC[1] > lev_performance$Pooled_AUC[1]) {
          cat("   → NACB shows superior EEG predictability\n")
        } else {
          cat("   → LEV shows superior EEG predictability\n")
        }
      }
    }
  }
  
  return(cv_results)
}

# =================================================================
# EXAMPLE USAGE
# =================================================================

# Example function to demonstrate usage
run_example_separate_cv <- function() {
  
  cat("Example: Running separate models CV on simulated data\n")
  
  # Create simulated dataset for demonstration
  set.seed(123)
  n <- 200
  example_data <- data.frame(
    OUTCOME_bin = rbinom(n, 1, 0.4),
    # EEG variables
    EEG_VAR1 = rnorm(n),
    EEG_VAR2 = rnorm(n), 
    EEG_VAR3 = rnorm(n),
    EEG_VAR4 = rnorm(n),
    # Clinical variables
    SEX = factor(sample(c("M", "F"), n, replace = TRUE)),
    AGE = rnorm(n, 45, 15),
    # Drug stratification
    Farmaco = factor(sample(c("LEV", "NACB"), n, replace = TRUE))
  )
  
  example_eeg_vars <- c("EEG_VAR1", "EEG_VAR2", "EEG_VAR3", "EEG_VAR4")
  example_clinical_vars <- c("SEX", "AGE")
  
  # Run separate models CV
  results <- run_separate_models_cv(
    database = example_data,
    eeg_variables = example_eeg_vars,
    clinical_variables = example_clinical_vars,
    drug_stratification = TRUE
  )
  
  return(results)
}

# Uncomment to run example:
# example_results <- run_example_separate_cv()

cat("Cross-validation functions loaded successfully!\n")
cat("Use run_5fold_cv() or run_separate_models_cv() for EEG and clinical models separately.\n")
