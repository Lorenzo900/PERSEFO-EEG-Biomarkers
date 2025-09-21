# =================================================================
# PERSEFO: BAYESIAN PRIOR SENSITIVITY ANALYSIS
# Testing robustness of results across different prior specifications
# =================================================================

library(rstanarm)
library(dplyr)
library(pROC)
library(bayestestR)
library(writexl)

# =================================================================
# PRIOR SCENARIOS DEFINITION
# =================================================================

prior_scenarios <- list(
  "Conservative" = list(
    name = "Conservative (σ=0.5)",
    coef_prior = normal(0, 0.5),      # Small effects expected
    intercept_prior = normal(0, 1),
    description = "Conservative priors expecting small effect sizes"
  ),
  "Current" = list(
    name = "Current Study (σ=1)", 
    coef_prior = normal(0, 1),        # Current study setting
    intercept_prior = normal(0, 2),
    description = "Main analysis priors used in the study"
  ),
  "Weakly_Informative" = list(
    name = "Weakly Informative (σ=2.5)",
    coef_prior = normal(0, 2.5),      # Medium effects expected
    intercept_prior = normal(0, 2.5),
    description = "Weakly informative priors allowing moderate effects"
  ),
  "Non_Informative" = list(
    name = "Non-Informative (σ=5)",
    coef_prior = normal(0, 5),        # No expectations
    intercept_prior = normal(0, 5),
    description = "Non-informative priors with minimal constraints"
  )
)

# =================================================================
# SENSITIVITY ANALYSIS FUNCTIONS
# =================================================================

run_sensitivity_scenario <- function(data, variables, scenario_info, outcome_var = "OUTCOME_bin") {
  
  cat("\n--- Testing:", scenario_info$name, "---\n")
  cat("Description:", scenario_info$description, "\n")
  
  formula_str <- paste(outcome_var, "~", paste(variables, collapse = " + "))
  
  tryCatch({
    # Fit model with specified priors
    model <- stan_glm(
      as.formula(formula_str),
      data = data,
      family = binomial(link = "logit"),
      prior = scenario_info$coef_prior,
      prior_intercept = scenario_info$intercept_prior,
      chains = 2,
      iter = 1000,
      cores = 2,
      refresh = 0,
      seed = 123
    )
    
    # Extract posterior summary
    posterior_summary <- describe_posterior(model, ci = 0.95)
    
    # Calculate AUC
    predictions <- posterior_predict(model, draws = 100)
    pred_mean <- apply(predictions, 2, mean)
    auc_value <- auc(roc(data[[outcome_var]], pred_mean, quiet = TRUE))
    
    # Identify "significant" variables (PD > 95%)
    significant_vars <- posterior_summary %>%
      filter(Parameter != "(Intercept)", pd > 0.95) %>%
      pull(Parameter)
    
    # Extract coefficient statistics
    coef_stats <- posterior_summary %>%
      filter(Parameter != "(Intercept)") %>%
      select(Parameter, Median, CI_low, CI_high, pd) %>%
      mutate(
        Effect_Size = case_when(
          abs(Median) >= 1.0 ~ "Large",
          abs(Median) >= 0.5 ~ "Medium", 
          abs(Median) >= 0.2 ~ "Small",
          TRUE ~ "Negligible"
        ),
        Evidence_Strength = case_when(
          pd >= 0.99 ~ "Very Strong",
          pd >= 0.95 ~ "Strong",
          pd >= 0.90 ~ "Moderate",
          pd >= 0.80 ~ "Weak", 
          TRUE ~ "Insufficient"
        )
      )
    
    return(list(
      scenario = scenario_info$name,
      auc = auc_value,
      n_significant = length(significant_vars),
      significant_vars = significant_vars,
      posterior_summary = posterior_summary,
      coefficient_stats = coef_stats,
      converged = TRUE,
      model_object = model
    ))
    
  }, error = function(e) {
    cat("ERROR:", e$message, "\n")
    return(list(
      scenario = scenario_info$name,
      auc = NA,
      n_significant = 0,
      significant_vars = character(0),
      converged = FALSE,
      error_message = e$message
    ))
  })
}

# =================================================================
# MAIN SENSITIVITY ANALYSIS FUNCTION
# =================================================================

run_sensitivity_analysis <- function(data, variables, model_name = "Test_Model", outcome_var = "OUTCOME_bin") {
  
  cat("=== BAYESIAN PRIOR SENSITIVITY ANALYSIS ===\n")
  cat("Model:", model_name, "\n")
  cat("Variables:", length(variables), "\n")
  cat("Sample size:", nrow(data), "\n")
  
  # Remove missing data
  data_complete <- data[complete.cases(data[c(outcome_var, variables)]), ]
  
  if(nrow(data_complete) < 30) {
    cat("WARNING: Sample size too small for reliable sensitivity analysis\n")
    return(NULL)
  }
  
  cat("Complete cases:", nrow(data_complete), "\n")
  
  # Run all scenarios
  sensitivity_results <- list()
  
  for(i in 1:length(prior_scenarios)) {
    scenario_name <- names(prior_scenarios)[i]
    scenario_info <- prior_scenarios[[i]]
    
    result <- run_sensitivity_scenario(data_complete, variables, scenario_info, outcome_var)
    sensitivity_results[[scenario_name]] <- result
  }
  
  # =================================================================
  # RESULTS COMPARISON
  # =================================================================
  
  cat("\n=== SENSITIVITY ANALYSIS RESULTS ===\n")
  
  # Create comparison table
  comparison_table <- do.call(rbind, lapply(sensitivity_results, function(res) {
    if(res$converged) {
      data.frame(
        Scenario = res$scenario,
        AUC = round(res$auc, 3),
        N_Significant_Vars = res$n_significant,
        Converged = "Yes",
        stringsAsFactors = FALSE
      )
    } else {
      data.frame(
        Scenario = res$scenario,
        AUC = "Failed",
        N_Significant_Vars = 0,
        Converged = "No",
        stringsAsFactors = FALSE
      )
    }
  }))
  
  print(knitr::kable(comparison_table, caption = "Sensitivity Analysis Summary"))
  
  # =================================================================
  # VARIABLE ROBUSTNESS ANALYSIS
  # =================================================================
  
  cat("\n=== VARIABLE ROBUSTNESS ANALYSIS ===\n")
  
  # Collect significant variables from all scenarios
  all_significant_vars <- list()
  for(i in 1:length(sensitivity_results)) {
    if(sensitivity_results[[i]]$converged) {
      scenario_name <- names(sensitivity_results)[i]
      all_significant_vars[[scenario_name]] <- sensitivity_results[[i]]$significant_vars
    }
  }
  
  # Identify robust variables (significant in all scenarios)
  robust_vars <- character(0)
  unstable_vars <- character(0)
  
  if(length(all_significant_vars) > 1) {
    robust_vars <- Reduce(intersect, all_significant_vars)
    all_mentioned_vars <- unique(unlist(all_significant_vars))
    unstable_vars <- setdiff(all_mentioned_vars, robust_vars)
    
    cat("ROBUST variables (significant in ALL scenarios):\n")
    if(length(robust_vars) > 0) {
      for(var in robust_vars) {
        cat("  ✓", var, "\n")
      }
    } else {
      cat("  No variables significant in all scenarios\n")
    }
    
    if(length(unstable_vars) > 0) {
      cat("\nUNSTABLE variables (significant only in SOME scenarios):\n")
      for(var in unstable_vars) {
        scenarios_with_var <- names(all_significant_vars)[sapply(all_significant_vars, function(x) var %in% x)]
        cat("  ⚠", var, "- significant in:", paste(scenarios_with_var, collapse = ", "), "\n")
      }
    }
  }
  
  # =================================================================
  # ROBUSTNESS ASSESSMENT
  # =================================================================
  
  cat("\n=== ROBUSTNESS INTERPRETATION ===\n")
  
  # Calculate AUC range
  successful_aucs <- sapply(sensitivity_results, function(x) if(x$converged) x$auc else NA)
  successful_aucs <- successful_aucs[!is.na(successful_aucs)]
  
  robustness_assessment <- list()
  
  if(length(successful_aucs) > 1) {
    auc_range <- max(successful_aucs) - min(successful_aucs)
    auc_mean <- mean(successful_aucs)
    auc_sd <- sd(successful_aucs)
    
    cat("AUC across scenarios:\n")
    cat("  Range:", round(auc_range, 3), "\n")
    cat("  Mean:", round(auc_mean, 3), "\n")
    cat("  SD:", round(auc_sd, 3), "\n")
    
    if(auc_range < 0.05) {
      robustness_verdict <- "ROBUST: AUC stable across priors"
      cat("  ✓ ROBUST: AUC stable across priors\n")
    } else if(auc_range < 0.10) {
      robustness_verdict <- "MODERATELY STABLE: Small AUC variations"
      cat("  ⚠ MODERATELY STABLE: Small AUC variations\n")
    } else {
      robustness_verdict <- "UNSTABLE: AUC varies significantly with priors"
      cat("  ⚠ UNSTABLE: AUC varies significantly with priors\n")
    }
    
    robustness_assessment <- list(
      auc_range = auc_range,
      auc_mean = auc_mean,
      auc_sd = auc_sd,
      verdict = robustness_verdict,
      robust_variables = robust_vars,
      unstable_variables = unstable_vars
    )
  }
  
  # =================================================================
  # SAVE RESULTS
  # =================================================================
  
  # Prepare detailed results for export
  detailed_results <- list()
  
  for(scenario_name in names(sensitivity_results)) {
    result <- sensitivity_results[[scenario_name]]
    if(result$converged) {
      detailed_results[[scenario_name]] <- result$coefficient_stats
    }
  }
  
  # Save to Excel
  excel_output <- list(
    Summary = comparison_table,
    Robustness_Assessment = if(length(robustness_assessment) > 0) {
      data.frame(
        Metric = c("AUC_Range", "AUC_Mean", "AUC_SD", "Verdict"),
        Value = c(round(robustness_assessment$auc_range, 3),
                 round(robustness_assessment$auc_mean, 3),
                 round(robustness_assessment$auc_sd, 3),
                 robustness_assessment$verdict)
      )
    } else {
      data.frame(Metric = "No_Assessment", Value = "Insufficient_Data")
    },
    Robust_Variables = if(length(robust_vars) > 0) {
      data.frame(Variable = robust_vars)
    } else {
      data.frame(Variable = "None")
    },
    Unstable_Variables = if(length(unstable_vars) > 0) {
      data.frame(Variable = unstable_vars)
    } else {
      data.frame(Variable = "None")
    }
  )
  
  # Add detailed coefficient results for each scenario
  for(scenario_name in names(detailed_results)) {
    excel_output[[paste0("Coefficients_", scenario_name)]] <- detailed_results[[scenario_name]]
  }
  
  filename <- paste0("Sensitivity_Analysis_", gsub("[^A-Za-z0-9]", "_", model_name), ".xlsx")
  write_xlsx(excel_output, filename)
  
  cat("\n✅ Results saved to:", filename, "\n")
  
  return(list(
    comparison_table = comparison_table,
    robustness_assessment = robustness_assessment,
    detailed_results = detailed_results,
    sensitivity_results = sensitivity_results
  ))
}

# =================================================================
# EXAMPLE USAGE
# =================================================================

# Example function to demonstrate usage
run_example_sensitivity <- function() {
  
  cat("Example: Running sensitivity analysis on simulated data\n")
  
  # Create simulated dataset for demonstration
  set.seed(123)
  n <- 100
  example_data <- data.frame(
    OUTCOME_bin = rbinom(n, 1, 0.4),
    EEG_VAR1 = rnorm(n),
    EEG_VAR2 = rnorm(n), 
    EEG_VAR3 = rnorm(n),
    CLINICAL_VAR = factor(sample(c("A", "B"), n, replace = TRUE))
  )
  
  example_variables <- c("EEG_VAR1", "EEG_VAR2", "EEG_VAR3", "CLINICAL_VAR")
  
  # Run sensitivity analysis
  results <- run_sensitivity_analysis(
    data = example_data,
    variables = example_variables,
    model_name = "Example_Model"
  )
  
  return(results)
}

# Uncomment to run example:
# example_results <- run_example_sensitivity()

cat("Sensitivity analysis functions loaded successfully!\n")
cat("Use run_sensitivity_analysis() with your data and variables.\n")
