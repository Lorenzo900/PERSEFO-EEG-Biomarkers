# =================================================================
# PERSEFO: STATISTICAL COMPARISONS BETWEEN MODELS
# DeLong test for EEG vs Clinical model comparisons
# =================================================================

library(pROC)
library(dplyr)
library(knitr)
library(writexl)

# =================================================================
# DELONG TEST FUNCTION FOR AUC COMPARISONS
# =================================================================

compare_auc_models <- function(cv_results_list, comparison_name = "Model Comparison") {
  
  cat("\n=== ", comparison_name, " ===\n")
  
  # Filter valid models only
  valid_models <- cv_results_list[!sapply(cv_results_list, is.null)]
  
  if(length(valid_models) < 2) {
    cat("Less than 2 valid models for comparison\n")
    return(NULL)
  }
  
  # Extract model names and AUC values
  model_names <- names(valid_models)
  aucs <- sapply(valid_models, function(x) x$pooled_auc)
  
  cat("Available models:\n")
  for(i in 1:length(model_names)) {
    cat(sprintf("  %s: AUC = %.3f\n", model_names[i], aucs[i]))
  }
  
  # Initialize comparison matrices
  n_models <- length(valid_models)
  comparison_matrix <- matrix(NA, nrow = n_models, ncol = n_models,
                              dimnames = list(model_names, model_names))
  
  p_value_matrix <- matrix(NA, nrow = n_models, ncol = n_models,
                           dimnames = list(model_names, model_names))
  
  detailed_comparisons <- list()
  
  # Pairwise comparisons
  for(i in 1:(n_models-1)) {
    for(j in (i+1):n_models) {
      
      model1_name <- model_names[i]
      model2_name <- model_names[j]
      
      cat(sprintf("\n--- Comparison: %s vs %s ---\n", model1_name, model2_name))
      
      # Extract predictions
      pred1 <- valid_models[[i]]$predictions
      pred2 <- valid_models[[j]]$predictions
      
      # Check if same patients (for paired test)
      if(!all(pred1$patient_id %in% pred2$patient_id) || 
         !all(pred2$patient_id %in% pred1$patient_id)) {
        cat("WARNING: Different patients in models - using unpaired test\n")
        
        # Unpaired test (less powerful)
        roc1 <- roc(pred1$observed, pred1$prob, quiet = TRUE)
        roc2 <- roc(pred2$observed, pred2$prob, quiet = TRUE)
        
        tryCatch({
          test_result <- roc.test(roc1, roc2, method = "delong")
          p_value <- test_result$p.value
          test_type <- "Unpaired DeLong"
        }, error = function(e) {
          cat("Error in DeLong test:", e$message, "\n")
          p_value <- NA
          test_type <- "Failed"
        })
        
      } else {
        # Paired test (more powerful)
        # Order by patient_id to ensure matching
        pred1_ordered <- pred1[order(pred1$patient_id), ]
        pred2_ordered <- pred2[order(pred2$patient_id), ]
        
        # Verify that observed outcomes are identical
        if(!all(pred1_ordered$observed == pred2_ordered$observed)) {
          cat("WARNING: Different outcomes - possible data error\n")
        }
        
        tryCatch({
          # Paired DeLong test
          roc1 <- roc(pred1_ordered$observed, pred1_ordered$prob, quiet = TRUE)
          roc2 <- roc(pred2_ordered$observed, pred2_ordered$prob, quiet = TRUE)
          
          test_result <- roc.test(roc1, roc2, method = "delong", paired = TRUE)
          p_value <- test_result$p.value
          test_type <- "Paired DeLong"
        }, error = function(e) {
          cat("Error in paired DeLong test:", e$message, "\n")
          p_value <- NA
          test_type <- "Failed"
        })
      }
      
      # Store results
      auc_diff <- aucs[i] - aucs[j]
      comparison_matrix[i, j] <- auc_diff
      comparison_matrix[j, i] <- -auc_diff
      p_value_matrix[i, j] <- p_value
      p_value_matrix[j, i] <- p_value
      
      # Interpretation
      if(!is.na(p_value)) {
        significance <- if(p_value < 0.001) "***" else if(p_value < 0.01) "**" else if(p_value < 0.05) "*" else "ns"
        direction <- if(auc_diff > 0) paste(model1_name, ">", model2_name) else paste(model2_name, ">", model1_name)
        
        cat(sprintf("AUC difference: %.3f\n", abs(auc_diff)))
        cat(sprintf("P-value (%s): %.4f %s\n", test_type, p_value, significance))
        cat(sprintf("Result: %s\n", direction))
      }
      
      # Store detailed results
      detailed_comparisons[[paste(model1_name, "vs", model2_name)]] <- list(
        model1 = model1_name,
        model2 = model2_name,
        auc1 = aucs[i],
        auc2 = aucs[j],
        auc_difference = auc_diff,
        p_value = p_value,
        test_type = test_type,
        significant = if(!is.na(p_value)) p_value < 0.05 else FALSE
      )
    }
  }
  
  return(list(
    comparison_name = comparison_name,
    model_aucs = aucs,
    auc_differences = comparison_matrix,
    p_values = p_value_matrix,
    detailed_comparisons = detailed_comparisons
  ))
}

# =================================================================
# EEG VS CLINICAL COMPARISON FUNCTION
# =================================================================

run_eeg_vs_clinical_comparisons <- function(cv_results_list) {
  
  cat("=== EEG vs CLINICAL STATISTICAL COMPARISONS ===\n")
  
  # Check available models
  valid_models <- cv_results_list[!sapply(cv_results_list, is.null)]
  cat("Available models for comparison:", length(valid_models), "\n")
  cat("Model names:", paste(names(valid_models), collapse = ", "), "\n")
  
  if(length(valid_models) < 2) {
    cat("ERROR: Need at least 2 valid models for comparisons\n")
    return(NULL)
  }
  
  # Identify clinical models
  clinical_models <- valid_models[grepl("Clinical|CLINICAL", names(valid_models), ignore.case = TRUE)]
  
  # Identify EEG-only models
  eeg_only_models <- valid_models[grepl("EEG.*Only|EEG_Only", names(valid_models), ignore.case = TRUE)]
  
  # Remove any mixed/combined models if present
  eeg_only_models <- eeg_only_models[!grepl("Combined|Mixed|and", names(eeg_only_models), ignore.case = TRUE)]
  
  comparison_results <- list()
  
  # 1. Clinical vs EEG comparisons
  if(length(clinical_models) > 0 && length(eeg_only_models) > 0) {
    cat("\n--- Clinical vs EEG-Only Comparisons ---\n")
    
    for(eeg_name in names(eeg_only_models)) {
      for(clinical_name in names(clinical_models)) {
        comparison_models <- list()
        comparison_models[[clinical_name]] <- clinical_models[[clinical_name]]
        comparison_models[[eeg_name]] <- eeg_only_models[[eeg_name]]
        
        comparison_result <- compare_auc_models(comparison_models, 
                                               paste("Clinical vs", eeg_name))
        comparison_results[[paste(clinical_name, "vs", eeg_name)]] <- comparison_result
      }
    }
  }
  
  # 2. Drug-specific EEG comparisons
  if(any(grepl("LEV|NACB", names(eeg_only_models)))) {
    cat("\n--- Drug-Specific EEG Comparisons ---\n")
    
    # LEV vs NACB EEG comparisons
    lev_eeg_models <- eeg_only_models[grepl("LEV.*EEG", names(eeg_only_models))]
    nacb_eeg_models <- eeg_only_models[grepl("NACB.*EEG", names(eeg_only_models))]
    
    if(length(lev_eeg_models) > 0 && length(nacb_eeg_models) > 0) {
      for(lev_name in names(lev_eeg_models)) {
        for(nacb_name in names(nacb_eeg_models)) {
          comparison_models <- list()
          comparison_models[[lev_name]] <- lev_eeg_models[[lev_name]]
          comparison_models[[nacb_name]] <- nacb_eeg_models[[nacb_name]]
          
          comparison_result <- compare_auc_models(comparison_models, 
                                                 "LEV vs NACB EEG Biomarkers")
          comparison_results[[paste(lev_name, "vs", nacb_name)]] <- comparison_result
        }
      }
    }
  }
  
  # 3. ALL vs Drug-specific EEG comparisons
  all_eeg_models <- eeg_only_models[grepl("ALL.*EEG", names(eeg_only_models))]
  drug_specific_eeg <- eeg_only_models[grepl("LEV.*EEG|NACB.*EEG", names(eeg_only_models))]
  
  if(length(all_eeg_models) > 0 && length(drug_specific_eeg) > 0) {
    cat("\n--- ALL vs Drug-Specific EEG Comparisons ---\n")
    
    for(all_name in names(all_eeg_models)) {
      for(drug_name in names(drug_specific_eeg)) {
        comparison_models <- list()
        comparison_models[[all_name]] <- all_eeg_models[[all_name]]
        comparison_models[[drug_name]] <- drug_specific_eeg[[drug_name]]
        
        comparison_result <- compare_auc_models(comparison_models, 
                                               paste("ALL vs", drug_name))
        comparison_results[[paste(all_name, "vs", drug_name)]] <- comparison_result
      }
    }
  }
  
  # =================================================================
  # CREATE SUMMARY TABLES
  # =================================================================
  
  cat("\n=== CREATING SUMMARY TABLES ===\n")
  
  create_comparison_table <- function(comparison_results) {
    
    all_comparisons <- unlist(comparison_results, recursive = FALSE)
    detailed_comparisons <- all_comparisons[grepl("detailed_comparisons", names(all_comparisons))]
    
    if(length(detailed_comparisons) == 0) {
      return(data.frame())
    }
    
    table_data <- list()
    
    for(comp_list_name in names(detailed_comparisons)) {
      comp_list <- detailed_comparisons[[comp_list_name]]
      
      if(is.list(comp_list)) {
        for(comp_name in names(comp_list)) {
          comp <- comp_list[[comp_name]]
          
          if(is.list(comp) && "auc1" %in% names(comp)) {
            
            # Categorize comparison type
            comparison_type <- case_when(
              grepl("Clinical.*EEG", comp_name) ~ "Clinical vs EEG",
              grepl("LEV.*NACB", comp_name) ~ "Drug-Specific",
              grepl("ALL.*LEV|ALL.*NACB", comp_name) ~ "Population vs Drug-Specific",
              TRUE ~ "Other"
            )
            
            table_data[[length(table_data) + 1]] <- data.frame(
              Comparison_Type = comparison_type,
              Comparison = comp_name,
              Model1 = comp$model1,
              Model1_AUC = round(comp$auc1, 3),
              Model2 = comp$model2,
              Model2_AUC = round(comp$auc2, 3),
              AUC_Difference = round(abs(comp$auc_difference), 3),
              P_Value = if(!is.na(comp$p_value)) {
                if(comp$p_value < 0.001) "<0.001" else sprintf("%.3f", comp$p_value)
              } else "Failed",
              Significance = if(!is.na(comp$p_value)) {
                if(comp$p_value < 0.001) "***" else if(comp$p_value < 0.01) "**" else if(comp$p_value < 0.05) "*" else "ns"
              } else "-",
              Test_Type = comp$test_type,
              stringsAsFactors = FALSE
            )
          }
        }
      }
    }
    
    if(length(table_data) > 0) {
      return(do.call(rbind, table_data))
    } else {
      return(data.frame())
    }
  }
  
  summary_table <- create_comparison_table(comparison_results)
  
  if(nrow(summary_table) > 0) {
    # Sort by comparison type and significance
    summary_table <- summary_table %>%
      arrange(Comparison_Type, desc(AUC_Difference)) %>%
      mutate(
        Significant = Significance %in% c("*", "**", "***"),
        Better_Model = ifelse(Model1_AUC > Model2_AUC, Model1, Model2)
      )
    
    print(knitr::kable(summary_table, caption = "EEG vs Clinical Statistical Comparisons"))
    
    # Highlight significant findings
    significant_comparisons <- summary_table[summary_table$Significant, ]
    
    if(nrow(significant_comparisons) > 0) {
      cat("\n=== SIGNIFICANT DIFFERENCES FOUND ===\n")
      for(i in 1:nrow(significant_comparisons)) {
        row <- significant_comparisons[i, ]
        cat(sprintf("• %s: %s (AUC=%.3f) significantly outperforms %s (AUC=%.3f), p=%s %s\n",
                   row$Comparison_Type, row$Better_Model, 
                   max(row$Model1_AUC, row$Model2_AUC),
                   ifelse(row$Better_Model == row$Model1, row$Model2, row$Model1),
                   min(row$Model1_AUC, row$Model2_AUC),
                   row$P_Value, row$Significance))
      }
      
      # Clinical interpretation
      cat("\n=== CLINICAL INTERPRETATION ===\n")
      
      clinical_vs_eeg <- significant_comparisons[significant_comparisons$Comparison_Type == "Clinical vs EEG", ]
      if(nrow(clinical_vs_eeg) > 0) {
        cat("• EEG biomarkers demonstrate statistical superiority over clinical variables\n")
      }
      
      drug_specific <- significant_comparisons[significant_comparisons$Comparison_Type == "Drug-Specific", ]
      if(nrow(drug_specific) > 0) {
        cat("• Evidence for drug-specific EEG biomarker patterns\n")
      }
      
    } else {
      cat("\n=== NO SIGNIFICANT DIFFERENCES FOUND ===\n")
      cat("Consider emphasizing confidence intervals and effect sizes in the manuscript.\n")
    }
    
    # Save results
    write_xlsx(list(
      All_Comparisons = summary_table,
      Significant_Only = significant_comparisons,
      Clinical_vs_EEG = summary_table[summary_table$Comparison_Type == "Clinical vs EEG", ],
      Drug_Specific = summary_table[summary_table$Comparison_Type == "Drug-Specific", ]
    ), "EEG_vs_Clinical_Comparisons.xlsx")
    
    cat("\n✅ Results saved to 'EEG_vs_Clinical_Comparisons.xlsx'\n")
    
  } else {
    cat("WARNING: No valid comparisons could be generated\n")
  }
  
  return(list(
    summary_table = summary_table,
    comparison_results = comparison_results,
    significant_comparisons = if(exists("significant_comparisons")) significant_comparisons else data.frame()
  ))
}

# =================================================================
# CLINICAL UTILITY ASSESSMENT FOR EEG MODELS
# =================================================================

assess_eeg_clinical_utility <- function(cv_results_list) {
  
  cat("\n=== EEG BIOMARKER CLINICAL UTILITY ASSESSMENT ===\n")
  
  # Focus only on EEG and clinical models
  eeg_models <- cv_results_list[grepl("EEG", names(cv_results_list))]
  clinical_models <- cv_results_list[grepl("Clinical", names(cv_results_list))]
  
  utility_assessment <- data.frame(
    Model = character(),
    Model_Type = character(),
    AUC = numeric(),
    Clinical_Utility = character(),
    Recommendation = character(),
    stringsAsFactors = FALSE
  )
  
  # Assess clinical models
  for(model_name in names(clinical_models)) {
    result <- clinical_models[[model_name]]
    if(!is.null(result)) {
      auc <- result$pooled_auc
      
      utility_assessment <- rbind(utility_assessment, data.frame(
        Model = model_name,
        Model_Type = "Clinical",
        AUC = round(auc, 3),
        Clinical_Utility = "Baseline reference",
        Recommendation = "Current standard of care",
        stringsAsFactors = FALSE
      ))
    }
  }
  
  # Assess EEG models
  for(model_name in names(eeg_models)) {
    result <- eeg_models[[model_name]]
    if(!is.null(result)) {
      auc <- result$pooled_auc
      
      if(auc >= 0.8) {
        utility <- "Excellent - Ready for clinical validation"
        recommendation <- "Proceed to prospective validation study"
      } else if(auc >= 0.7) {
        utility <- "Good - Promising for clinical implementation"
        recommendation <- "Consider multicenter validation"
      } else if(auc >= 0.6) {
        utility <- "Fair - Research stage"
        recommendation <- "Investigate additional biomarkers"
      } else {
        utility <- "Poor - Limited predictive value"
        recommendation <- "Reconsider approach"
      }
      
      utility_assessment <- rbind(utility_assessment, data.frame(
        Model = model_name,
        Model_Type = "EEG",
        AUC = round(auc, 3),
        Clinical_Utility = utility,
        Recommendation = recommendation,
        stringsAsFactors = FALSE
      ))
    }
  }
  
  utility_assessment <- utility_assessment %>% arrange(desc(AUC))
  
  print(knitr::kable(utility_assessment, caption = "Clinical Utility Assessment: EEG vs Clinical"))
  
  return(utility_assessment)
}

# =================================================================
# EXAMPLE USAGE FUNCTION
# =================================================================

run_example_eeg_clinical_comparisons <- function() {
  
  cat("Example: EEG vs Clinical comparisons on simulated data\n")
  
  # Create simulated CV results for demonstration
  set.seed(123)
  
  # Simulate predictions for different models
  n_patients <- 100
  patient_ids <- paste0("P", 1:n_patients)
  true_outcomes <- rbinom(n_patients, 1, 0.4)
  
  # Clinical model (poor performance)
  clinical_probs <- runif(n_patients, 0.3, 0.7)
  clinical_predictions <- data.frame(
    patient_id = patient_ids,
    observed = true_outcomes,
    prob = clinical_probs,
    fold = sample(1:5, n_patients, replace = TRUE)
  )
  
  # NACB EEG model (good performance)
  nacb_probs <- ifelse(true_outcomes == 1, runif(n_patients, 0.5, 0.9), runif(n_patients, 0.1, 0.5))
  nacb_predictions <- data.frame(
    patient_id = patient_ids,
    observed = true_outcomes,
    prob = nacb_probs,
    fold = sample(1:5, n_patients, replace = TRUE)
  )
  
  # LEV EEG model (poor performance)
  lev_probs <- runif(n_patients, 0.3, 0.7)
  lev_predictions <- data.frame(
    patient_id = patient_ids,
    observed = true_outcomes,
    prob = lev_probs,
    fold = sample(1:5, n_patients, replace = TRUE)
  )
  
  # Create mock CV results
  example_cv_results <- list(
    "Clinical_Only" = list(
      model_name = "Clinical_Only",
      pooled_auc = auc(roc(clinical_predictions$observed, clinical_predictions$prob, quiet = TRUE)),
      predictions = clinical_predictions
    ),
    "NACB_EEG_Only" = list(
      model_name = "NACB_EEG_Only", 
      pooled_auc = auc(roc(nacb_predictions$observed, nacb_predictions$prob, quiet = TRUE)),
      predictions = nacb_predictions
    ),
    "LEV_EEG_Only" = list(
      model_name = "LEV_EEG_Only",
      pooled_auc = auc(roc(lev_predictions$observed, lev_predictions$prob, quiet = TRUE)),
      predictions = lev_predictions
    )
  )
  
  # Run comparisons
  comparison_results <- run_eeg_vs_clinical_comparisons(example_cv_results)
  
  # Assess clinical utility
  utility_results <- assess_eeg_clinical_utility(example_cv_results)
  
  return(list(
    comparisons = comparison_results,
    utility = utility_results
  ))
}

# Uncomment to run example:
# example_comparison_results <- run_example_eeg_clinical_comparisons()

cat("Statistical comparison functions loaded successfully!\n")
cat("Use run_eeg_vs_clinical_comparisons() with your CV results.\n")
