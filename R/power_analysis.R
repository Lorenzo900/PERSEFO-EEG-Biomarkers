# =================================================================
# PERSEFO: STATISTICAL POWER ANALYSIS FOR ROC CURVES
# Hanley & McNeil (1982) framework for sample size and detectable effect sizes
# =================================================================

library(pROC)
library(dplyr)
library(knitr)
library(writexl)

# =================================================================
# POWER ANALYSIS FUNCTIONS (HANLEY & MCNEIL 1982)
# =================================================================

calculate_roc_power <- function(n, prevalence, auc_null = 0.5, auc_alt = NULL, alpha = 0.05, power_target = 0.80) {
  
  # Calculate number of positive and negative cases
  n_pos <- round(n * prevalence)
  n_neg <- n - n_pos
  
  # Check minimum case requirements
  if(n_pos < 5 || n_neg < 5) {
    return(list(
      n_total = n,
      n_positive = n_pos,
      n_negative = n_neg,
      prevalence = prevalence,
      power = NA,
      detectable_auc = NA,
      se_auc = NA,
      error = "Insufficient cases"
    ))
  }
  
  # Hanley & McNeil (1982) simplified formula
  # Variance of AUC under H0 (AUC = 0.5)
  var_auc_h0 <- (1/12) * ((1/n_pos) + (1/n_neg))
  
  # Standard error
  se_auc <- sqrt(var_auc_h0)
  
  # Critical Z-score for alpha
  z_alpha <- qnorm(1 - alpha/2)
  
  # Z-score for desired power
  z_beta <- qnorm(power_target)
  
  # Minimum detectable AUC
  min_detectable_auc <- auc_null + (z_alpha + z_beta) * se_auc
  
  # Power for specific AUC
  if(!is.null(auc_alt)) {
    z_stat <- (auc_alt - auc_null) / se_auc
    power <- pnorm(z_stat - z_alpha) + pnorm(-z_stat - z_alpha)
  } else {
    power <- power_target
    auc_alt <- min_detectable_auc
  }
  
  return(list(
    n_total = n,
    n_positive = n_pos,
    n_negative = n_neg,
    prevalence = prevalence,
    power = power,
    detectable_auc = min_detectable_auc,
    se_auc = se_auc,
    alpha = alpha,
    power_target = power_target
  ))
}

# =================================================================
# COMPREHENSIVE POWER ANALYSIS FUNCTION
# =================================================================

run_comprehensive_power_analysis <- function(total_n = 232, drug_distribution = NULL, 
                                            prevalence_range = c(0.4, 0.5, 0.6),
                                            alpha = 0.05, power_target = 0.80) {
  
  cat("=== COMPREHENSIVE POWER ANALYSIS ===\n")
  cat("Based on Hanley & McNeil (1982) normal approximation for AUC\n")
  cat("Total sample size:", total_n, "\n")
  cat("Significance level (α):", alpha, "\n")
  cat("Target power:", power_target * 100, "%\n")
  
  # Default drug distribution if not provided
  if(is.null(drug_distribution)) {
    drug_distribution <- list(
      LEV = round(total_n * 0.51),   # ~51%
      NACB = round(total_n * 0.30),  # ~30%
      Other = round(total_n * 0.19)  # ~19%
    )
  }
  
  cat("Drug distribution:\n")
  for(drug in names(drug_distribution)) {
    cat(sprintf("  %s: %d patients\n", drug, drug_distribution[[drug]]))
  }
  
  # =================================================================
  # POWER ANALYSIS FOR COMPLETE DATASET
  # =================================================================
  
  cat("\n=== POWER ANALYSIS - COMPLETE DATASET ===\n")
  
  complete_dataset_power <- list()
  
  for(prev in prevalence_range) {
    result <- calculate_roc_power(n = total_n, prevalence = prev, 
                                 alpha = alpha, power_target = power_target)
    complete_dataset_power[[paste0("prev_", prev*100)]] <- result
    
    cat(sprintf("Prevalence %.0f%% (%d SF / %d NSF):\n", 
                prev*100, result$n_positive, result$n_negative))
    cat(sprintf("  Minimum detectable AUC (%.0f%% power): %.3f\n", 
                power_target*100, result$detectable_auc))
    cat(sprintf("  Standard Error AUC: %.3f\n", result$se_auc))
  }
  
  # =================================================================
  # POWER ANALYSIS FOR DRUG-SPECIFIC SUBGROUPS
  # =================================================================
  
  cat("\n=== POWER ANALYSIS - DRUG-SPECIFIC SUBGROUPS ===\n")
  
  drug_specific_power <- list()
  
  for(drug in names(drug_distribution)) {
    if(drug == "Other") next  # Focus on LEV and NACB
    
    n_drug <- drug_distribution[[drug]]
    cat(sprintf("\n--- %s Subgroup (n=%d) ---\n", drug, n_drug))
    
    drug_power <- list()
    
    for(prev in prevalence_range) {
      result <- calculate_roc_power(n = n_drug, prevalence = prev,
                                   alpha = alpha, power_target = power_target)
      drug_power[[paste0("prev_", prev*100)]] <- result
      
      cat(sprintf("Prevalence %.0f%% (%d SF / %d NSF):\n", 
                  prev*100, result$n_positive, result$n_negative))
      cat(sprintf("  Minimum detectable AUC: %.3f\n", result$detectable_auc))
    }
    
    drug_specific_power[[drug]] <- drug_power
  }
  
  # =================================================================
  # POWER FOR SPECIFIC AUC VALUES
  # =================================================================
  
  cat("\n=== POWER FOR SPECIFIC AUC VALUES ===\n")
  
  target_aucs <- c(0.60, 0.65, 0.70, 0.75, 0.80)
  prevalence_test <- 0.5  # Use 50% prevalence for this analysis
  
  power_for_aucs <- data.frame(
    Target_AUC = target_aucs,
    Power_ALL = sapply(target_aucs, function(auc) {
      calculate_roc_power(total_n, prevalence_test, auc_alt = auc, 
                         alpha = alpha, power_target = power_target)$power
    }),
    stringsAsFactors = FALSE
  )
  
  # Add drug-specific power calculations
  for(drug in c("LEV", "NACB")) {
    if(drug %in% names(drug_distribution)) {
      power_for_aucs[[paste0("Power_", drug)]] <- sapply(target_aucs, function(auc) {
        calculate_roc_power(drug_distribution[[drug]], prevalence_test, auc_alt = auc,
                           alpha = alpha, power_target = power_target)$power
      })
    }
  }
  
  # Round power values
  power_for_aucs[, -1] <- round(power_for_aucs[, -1], 3)
  
  print(knitr::kable(power_for_aucs, caption = "Statistical Power for Specific AUC Values (50% Prevalence)"))
  
  # =================================================================
  # CREATE SUMMARY TABLE
  # =================================================================
  
  cat("\n=== SUMMARY TABLE FOR MANUSCRIPT ===\n")
  
  power_summary <- data.frame(
    Dataset = character(),
    N_Total = numeric(),
    Prevalence = character(),
    N_SF = numeric(),
    N_NSF = numeric(),
    Min_Detectable_AUC_80pct = numeric(),
    stringsAsFactors = FALSE
  )
  
  # Add complete dataset results
  for(i in 1:length(prevalence_range)) {
    prev <- prevalence_range[i]
    result <- complete_dataset_power[[i]]
    
    power_summary <- rbind(power_summary, data.frame(
      Dataset = "ALL",
      N_Total = total_n,
      Prevalence = paste0(prev*100, "%"),
      N_SF = result$n_positive,
      N_NSF = result$n_negative,
      Min_Detectable_AUC_80pct = round(result$detectable_auc, 3),
      stringsAsFactors = FALSE
    ))
  }
  
  # Add drug-specific results
  for(drug in c("LEV", "NACB")) {
    if(drug %in% names(drug_specific_power)) {
      for(i in 1:length(prevalence_range)) {
        prev <- prevalence_range[i]
        result <- drug_specific_power[[drug]][[i]]
        
        power_summary <- rbind(power_summary, data.frame(
          Dataset = drug,
          N_Total = drug_distribution[[drug]],
          Prevalence = paste0(prev*100, "%"),
          N_SF = result$n_positive,
          N_NSF = result$n_negative,
          Min_Detectable_AUC_80pct = round(result$detectable_auc, 3),
          stringsAsFactors = FALSE
        ))
      }
    }
  }
  
  print(knitr::kable(power_summary, caption = "Power Analysis Summary (80% Power, α=0.05)"))
  
  # =================================================================
  # CLINICAL INTERPRETATION
  # =================================================================
  
  cat("\n=== CLINICAL INTERPRETATION ===\n")
  
  # Find best and worst case scenarios
  best_case_auc <- min(power_summary$Min_Detectable_AUC_80pct[power_summary$Dataset == "ALL"])
  worst_case_idx <- which.max(power_summary$Min_Detectable_AUC_80pct[power_summary$Dataset != "ALL"])
  worst_case_dataset <- power_summary$Dataset[power_summary$Dataset != "ALL"][worst_case_idx]
  worst_case_auc <- max(power_summary$Min_Detectable_AUC_80pct[power_summary$Dataset != "ALL"])
  
  cat(sprintf("• Complete dataset (n=%d): can detect AUC ≥ %.3f with 80%% power\n", 
              total_n, best_case_auc))
  cat(sprintf("• Drug-specific analyses: greatest limitation for %s subgroup (AUC ≥ %.3f)\n", 
              worst_case_dataset, worst_case_auc))
  
  # Recommendations based on power analysis
  if(worst_case_auc > 0.75) {
    cat("• Recommendation: Drug-specific analyses have limited power for detecting moderate effects\n")
    cat("• Consider emphasizing confidence intervals and effect sizes in interpretation\n")
  } else if(worst_case_auc > 0.70) {
    cat("• Recommendation: Drug-specific analyses have adequate power for clinically meaningful effects\n")
  } else {
    cat("• Recommendation: Drug-specific analyses have good power for detecting moderate to large effects\n")
  }
  
  # =================================================================
  # COMPARISON WITH PREVIOUS STUDIES
  # =================================================================
  
  cat("\n=== COMPARISON WITH LITERATURE ===\n")
  
  # Example comparison (adapt to your specific previous study)
  previous_n <- 202
  comparison_prev <- 0.5
  
  old_power <- calculate_roc_power(n = previous_n, prevalence = comparison_prev,
                                  alpha = alpha, power_target = power_target)
  new_power <- calculate_roc_power(n = total_n, prevalence = comparison_prev,
                                  alpha = alpha, power_target = power_target)
  
  cat("Comparison with previous study:\n")
  cat(sprintf("  Previous sample (n=%d): minimum detectable AUC = %.3f\n", 
              previous_n, old_power$detectable_auc))
  cat(sprintf("  Current study (n=%d): minimum detectable AUC = %.3f\n", 
              total_n, new_power$detectable_auc))
  cat(sprintf("  Improvement: %.3f AUC points\n", 
              old_power$detectable_auc - new_power$detectable_auc))
  
  # =================================================================
  # SAVE RESULTS
  # =================================================================
  
  write_xlsx(list(
    Power_Summary = power_summary,
    Power_for_AUCs = power_for_aucs,
    Complete_Dataset = do.call(rbind, lapply(complete_dataset_power, function(x) {
      data.frame(
        Prevalence = paste0(x$prevalence*100, "%"),
        N_Positive = x$n_positive,
        N_Negative = x$n_negative,
        Min_Detectable_AUC = round(x$detectable_auc, 3),
        SE_AUC = round(x$se_auc, 3)
      )
    })),
    Drug_Specific = if(length(drug_specific_power) > 0) {
      do.call(rbind, lapply(names(drug_specific_power), function(drug) {
        do.call(rbind, lapply(drug_specific_power[[drug]], function(x) {
          data.frame(
            Drug = drug,
            Prevalence = paste0(x$prevalence*100, "%"),
            N_Positive = x$n_positive,
            N_Negative = x$n_negative,
            Min_Detectable_AUC = round(x$detectable_auc, 3),
            SE_AUC = round(x$se_auc, 3)
          )
        }))
      }))
    } else {
      data.frame(Note = "No drug-specific analyses")
    }
  ), "Power_Analysis_Results.xlsx")
  
  cat("\n✅ Power analysis results saved to 'Power_Analysis_Results.xlsx'\n")
  
  return(list(
    power_summary = power_summary,
    power_for_aucs = power_for_aucs,
    complete_dataset_power = complete_dataset_power,
    drug_specific_power = drug_specific_power
  ))
}

# =================================================================
# SAMPLE SIZE CALCULATION FUNCTION
# =================================================================

calculate_required_sample_size <- function(target_auc = 0.70, prevalence = 0.5, 
                                          alpha = 0.05, power = 0.80, auc_null = 0.5) {
  
  cat("=== SAMPLE SIZE CALCULATION ===\n")
  cat(sprintf("Target AUC: %.2f\n", target_auc))
  cat(sprintf("Prevalence: %.0f%%\n", prevalence * 100))
  cat(sprintf("α: %.2f, Power: %.0f%%\n", alpha, power * 100))
  
  # Z-scores
  z_alpha <- qnorm(1 - alpha/2)
  z_beta <- qnorm(power)
  
  # Effect size (difference from null)
  effect_size <- target_auc - auc_null
  
  # Approximate sample size calculation
  # Based on simplified Hanley & McNeil formula
  # n = (z_α + z_β)² / [12 * (AUC - 0.5)² * prevalence * (1 - prevalence)]
  
  denominator <- 12 * effect_size^2 * prevalence * (1 - prevalence)
  n_required <- (z_alpha + z_beta)^2 / denominator
  
  n_pos_required <- round(n_required * prevalence)
  n_neg_required <- round(n_required * (1 - prevalence))
  
  cat(sprintf("Required total sample size: %.0f patients\n", ceiling(n_required)))
  cat(sprintf("Required positive cases: %.0f\n", n_pos_required))
  cat(sprintf("Required negative cases: %.0f\n", n_neg_required))
  
  return(list(
    total_n = ceiling(n_required),
    n_positive = n_pos_required,
    n_negative = n_neg_required,
    target_auc = target_auc,
    prevalence = prevalence,
    alpha = alpha,
    power = power
  ))
}

# =================================================================
# EXAMPLE USAGE
# =================================================================

run_example_power_analysis <- function() {
  
  cat("Example: Power analysis for epilepsy biomarker study\n")
  
  # Run comprehensive power analysis
  power_results <- run_comprehensive_power_analysis(
    total_n = 232,
    drug_distribution = list(LEV = 119, NACB = 69, Other = 44),
    prevalence_range = c(0.4, 0.5, 0.6),
    alpha = 0.05,
    power_target = 0.80
  )
  
  # Calculate required sample sizes for different targets
  cat("\n=== SAMPLE SIZE REQUIREMENTS FOR DIFFERENT TARGETS ===\n")
  
  targets <- c(0.65, 0.70, 0.75, 0.80)
  for(target in targets) {
    cat(sprintf("\nFor AUC = %.2f:\n", target))
    sample_size_result <- calculate_required_sample_size(target_auc = target)
  }
  
  return(power_results)
}

# Uncomment to run example:
# example_power_results <- run_example_power_analysis()

cat("Power analysis functions loaded successfully!\n")
cat("Use run_comprehensive_power_analysis() for your study parameters.\n")
