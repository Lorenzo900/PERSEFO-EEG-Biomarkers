# =================================================================
# FEATURE IMPORTANCE AND SIGNIFICANT VARIABLE EXTRACTION
# Works with models trained in cross_validation.R
# =================================================================

library(rstanarm)
library(dplyr)
library(knitr)
library(writexl)
library(broom.mixed)

# =================================================================
# FUNCTION: EXTRACT FEATURE IMPORTANCE FROM STAN MODELS
# =================================================================

extract_feature_importance <- function(stan_model, model_name = "Unknown_Model", ci_level = 0.95) {
  
  cat("\n=== Extracting feature importance for:", model_name, "===\n")
  
  if (is.null(stan_model)) {
    warning("Model object is NULL, skipping...")
    return(NULL)
  }
  
  # Extract posterior summaries
  posterior_summary <- broom.mixed::tidy(stan_model, effects = "fixed", conf.int = TRUE, conf.level = ci_level)
  
  posterior_summary <- posterior_summary %>%
    filter(term != "(Intercept)") %>%
    mutate(
      significant = ifelse(conf.low * conf.high > 0, TRUE, FALSE),
      importance = abs(estimate),
      ci_width = conf.high - conf.low
    ) %>%
    arrange(desc(importance))
  
  cat(sprintf("Extracted %d coefficients (%.0f%% credible interval)\n",
              nrow(posterior_summary), ci_level * 100))
  
  return(posterior_summary)
}

# =================================================================
# FUNCTION: SUMMARIZE MULTIPLE MODEL IMPORTANCES
# =================================================================

summarize_model_importances <- function(model_list) {
  
  importance_results <- list()
  
  for (model_name in names(model_list)) {
    model_obj <- model_list[[model_name]]
    
    if (!inherits(model_obj, "stanreg")) {
      cat(sprintf("Skipping %s (not a stan_glm object)\n", model_name))
      next
    }
    
    model_importance <- extract_feature_importance(model_obj, model_name)
    if (!is.null(model_importance)) {
      importance_results[[model_name]] <- model_importance %>%
        mutate(Model = model_name)
    }
  }
  
  if (length(importance_results) == 0) {
    cat("⚠️ No valid models found for importance extraction\n")
    return(NULL)
  }
  
  combined_importance <- bind_rows(importance_results)
  
  # Save Excel summary
  write_xlsx(
    list(
      Feature_Importance_All = combined_importance,
      Top10_Per_Model = combined_importance %>%
        group_by(Model) %>%
        slice_max(order_by = importance, n = 10) %>%
        ungroup()
    ),
    "Feature_Importance_Summary.xlsx"
  )
  
  cat("\n✅ Feature importance summary saved to 'Feature_Importance_Summary.xlsx'\n")
  
  print(knitr::kable(
    combined_importance %>% 
      group_by(Model) %>%
      slice_max(order_by = importance, n = 5) %>%
      ungroup(),
    caption = "Top 5 Most Important Variables per Model"
  ))
  
  return(combined_importance)
}

# =================================================================
# OPTIONAL: EXAMPLE DEMONSTRATION
# =================================================================

run_example_feature_importance <- function() {
  cat("Example: Simulating stan_glm models for demonstration\n")
  
  # Example dataset
  set.seed(42)
  n <- 150
  example_data <- data.frame(
    OUTCOME_bin = rbinom(n, 1, 0.4),
    EEG_VAR1 = rnorm(n),
    EEG_VAR2 = rnorm(n),
    EEG_VAR3 = rnorm(n),
    AGE = rnorm(n, 45, 15),
    SEX = factor(sample(c("M", "F"), n, replace = TRUE))
  )
  
  # Fit quick example model
  example_model <- stan_glm(
    OUTCOME_bin ~ EEG_VAR1 + EEG_VAR2 + EEG_VAR3 + AGE + SEX,
    data = example_data,
    family = binomial(link = "logit"),
    chains = 2, iter = 1000, cores = 2, refresh = 0
  )
  
  # Extract importance
  importance_summary <- summarize_model_importances(
    list(EEG_Model = example_model)
  )
  
  return(importance_summary)
}

# Uncomment to test:
# example_importance <- run_example_feature_importance()

cat("Feature importance extraction functions loaded successfully!\n")
cat("Use summarize_model_importances(list_of_stan_models) to extract and save variable importances.\n")
