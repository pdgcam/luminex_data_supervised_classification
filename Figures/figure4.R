

# --- Source functions
source(here('Models/train_binary_models.R'))
source(here('Models/train_multinomial_models.R'))

    
# Functions to select targets / classification questions 
select_targets <- function(preprocessed_data,
                           targets = c("flavi", "dengue", 
                                       "zika", "dengue_zika", 
                                       "dengue_serotype", "dengue_serotype_neg",
                                       "dengue_chik"),
                           drop_original_target = TRUE, 
                           negative_label = "negative") {
  
  if (!"target" %in% names(preprocessed_data)) {
    stop('Column "Target" not found in preprocessed_data.')
  }
  
  # Validate all targets
  valid_targets <- c("flavi", "dengue", 
                     "zika", "dengue_zika", 
                     "dengue_serotype", "dengue_serotype_neg",
                     "dengue_chik")
  
  invalid <- setdiff(targets, valid_targets)
  if (length(invalid) > 0) {
    stop("Invalid target(s): ", paste(invalid, collapse = ", "))
  }
  
  # mappings
  mapping_list <- list(
    flavi               = c(DENV1 = 1, DENV2 = 1, DENV3 = 1, DENV4 = 1, ZIKV = 1),
    dengue              = c(DENV1 = 1, DENV2 = 1, DENV3 = 1, DENV4 = 1),
    zika                = c(ZIKV = 1),
    dengue_zika         = c(DENV1 = 1, DENV2 = 1, DENV3 = 1, DENV4 = 1, ZIKV = 2),
    dengue_serotype     = c(DENV1 = 1, DENV2 = 2, DENV3 = 3, DENV4 = 4),
    dengue_serotype_neg = c(DENV1 = 1, DENV2 = 2, DENV3 = 3, DENV4 = 4),
    dengue_chik         = c(DENV1 = 0, DENV2 = 0, DENV3 = 0, DENV4 = 0, CHIKV = 1))
  
  # factor levels/labels
  label_specs <- list(
    flavi               = list(levels = c(0, 1),           labels = c("negative", "positive")),
    dengue              = list(levels = c(0, 1),           labels = c("negative", "positive")),
    zika                = list(levels = c(0, 1),           labels = c("negative", "positive")),
    dengue_zika         = list(levels = c(0, 1, 2),        labels = c("negative", "dengue", "zika")),
    dengue_serotype     = list(levels = c(1, 2, 3, 4),     labels = c("DENV1", "DENV2", "DENV3", "DENV4")),
    dengue_serotype_neg = list(levels = c(0, 1, 2, 3, 4),  labels = c("negative", "DENV1", "DENV2", "DENV3", "DENV4")),
    dengue_chik         = list(levels = c(0, 1),           labels = c("dengue", "CHIKV"))
  )
  
  
  # List to store results
  data_with_target <- list()
  
  # Process each target
  for (target in targets) {
    df_copy <- preprocessed_data
    
    # Build the chosen target column
    mp <- mapping_list[[target]]
    tgt_chr <- as.character(df_copy$target)
    
    # Map values: negative label get 0, mapped values get their code, rest get NA
    df_copy[[target]] <- ifelse(
      tgt_chr %in% negative_label, 
      0,
      as.numeric(dplyr::recode(tgt_chr, !!!mp, .default = NA_real_))
    )
    
    # Remove negatives for classifications that only include infected samples (eg: given infection, classify dengue vs chik)
    if (target %in% c("dengue_chik", "dengue_serotype")) {
      df_copy <- df_copy[!df_copy$target %in% negative_label, ]
    }

    # convert to factor with labels
    spec <- label_specs[[target]]
    df_copy[[target]] <-
      factor(df_copy[[target]], 
             levels = spec$levels, 
             labels = spec$labels)
    
    # drop original Target column if requested
    if (drop_original_target) {
      if ("target" %in% names(df_copy)) {
        df_copy[["target"]] <- NULL
      }
    }
    
    data_with_target[[target]] <- df_copy
  }
  return(data_with_target)
}




# Function to run binomial and multinomial classification with multiple targets simulatensouly 
train_multiple_targets <- function(
    data_list,  
    variables = NULL,
    k_fold = 5,
    metrics = c("AUROC", "AUPRC", "Brier", "StratBrier")) {
  
  # Store all results
  all_results <- list()
  all_comparisons <- list()
  all_predictions <- list()
  
  # Loop through each target dataset
  for (target_name in names(data_list)) {
    cat("Processing target:", target_name, "\n")
    
    current_data <- data_list[[target_name]]
    target_col <- paste0(target_name)
    
    # Check if target column exists
    if (!target_col %in% names(current_data)) {
      warning(paste("Target column", target_col, "not found in", target_name, "data. Skipping."))
      next
    }
    
    # Determine class type based on target variable
    target_values <- unique(current_data[[target_col]])
    n_classes <- length(target_values)
    print(n_classes)
    
    if (n_classes == 2) {
      # Train models for this target
      result <- train_binary_models(
        data = current_data,
        target = target_col,
        variables = variables,
        k_fold = k_fold,
        metrics = metrics
      )
      
      # Store results
      all_results[[target_name]] <- result
      
      # Add target name to comparison dataframe
      comparison_with_target <- result$comparison
      comparison_with_target$target <- target_name
      all_comparisons[[target_name]] <- comparison_with_target
      
      # Add target name to predictions dataframe
      predictions_with_target <- result$predictions
      predictions_with_target$target <- target_name
      all_predictions[[target_name]] <- predictions_with_target
    }
    
    else if (n_classes > 2){
      # Train models for this target
      result <- train_multinomial_models(
        data = current_data,
        target = target_col,
        variables = variables,
        k_fold = "LOOCV", 
        metrics = metrics
      )
      
      # Store results
      all_results[[target_name]] <- result
      
      # Add target name to comparison dataframe
      comparison_with_target <- result$comparison
      comparison_with_target$target <- target_name
      all_comparisons[[target_name]] <- comparison_with_target
      
      # Add target name to predictions dataframe
      predictions_with_target <- result$predictions
      predictions_with_target$target <- target_name
      all_predictions[[target_name]] <- predictions_with_target
    }
    
    else {
      warning(paste("Target", target_col, "has only 1 class. Skipping."))
      next
    }
  }
  
  # Combine all comparison dataframes
  combined_comparison <- do.call(rbind, all_comparisons)
  rownames(combined_comparison) <- NULL
  # Reorder columns to put Target first
  col_order <- c("target", setdiff(names(combined_comparison), "target"))
  combined_comparison <- combined_comparison[, col_order]
  
  # Combine all predictions dataframes
  combined_predictions <- dplyr::bind_rows(all_predictions)
  rownames(combined_predictions) <- NULL
  # Reorder columns to put Target first
  pred_col_order <- c("target", setdiff(names(combined_predictions), "target"))
  combined_predictions <- combined_predictions[, pred_col_order]
  
  return(list(
    results_by_target = all_results,
    combined_comparison = combined_comparison,
    combined_predictions = combined_predictions,
    summary = list(
      n_targets = length(all_results),
      target_names = names(all_results),
      k_fold = k_fold,
      metrics = metrics
    )
  ))
}


# ---- Import prepossessed datasets ---- 
ratio_df <- readRDS('Results/ratio_df.rds')
table(ratio_df$target)


# ---- Binary Results 
# Define targets / classification question
data_with_binomial_targets <- select_targets(
  preprocessed_data = ratio_df,
  targets = c("flavi", "dengue"),
  drop_original_target = FALSE,
  negative_label =  c("negative", "CHIKV")  # both treated as class 0
)

table(data_with_binomial_targets$flavi$flavi)
sum(is.na(data_with_binomial_targets$dengue$zika))
# Drop Nas
data_with_binomial_targets$flavi <- na.omit(data_with_binomial_targets$flavi)
data_with_binomial_targets$dengue <- na.omit(data_with_binomial_targets$dengue)
data_with_binomial_targets$zika <- na.omit(data_with_binomial_targets$zika)

# View distribution of targets
table(data_with_binomial_targets$flavi$flavi)
table(data_with_binomial_targets$dengue$dengue)

# Fit binomial models
binomial_modeling_results <- train_multiple_targets(
  data_list = data_with_binomial_targets,
  variables = NULL,  # Uses all columns except target
  k_fold = 5,
  metrics = c("AUROC", "AUPRC", "Brier")
)
warnings()
traceback()
binomial_modeling_results$combined_comparison

# --- Multinomial Results 
data_with_multinomial_targets <- select_targets(
  preprocessed_data = ratio_df,
  targets = c("dengue_serotype", "dengue_serotype_neg"),
  drop_original_target = TRUE
)

# Drop Nas
data_with_multinomial_targets$dengue_serotype <- na.omit(data_with_multinomial_targets$dengue_serotype)
data_with_multinomial_targets$dengue_serotype_neg <- na.omit(data_with_multinomial_targets$dengue_serotype_neg)

# View distribution of targets
table(data_with_multinomial_targets$dengue_serotype$dengue_serotype)
table(data_with_multinomial_targets$dengue_serotype_neg$dengue_serotype_neg)


# Given infection, classify stereotype + serotype and neg 
dengue_serotype_results <- train_multiple_targets(
  data = data_with_multinomial_targets,
  variables = NULL,
  k_fold = "LOOCV", 
  metrics = c("AUROC", "AUPRC", "Brier", "StratBrier"))



# Model results table 