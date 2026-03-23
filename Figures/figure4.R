

# --- Source functions
source(here('Models/train_binary_models.R'))
source(here('Models/train_multinomial_models.R'))

    
# Functions to select targets / classification questions 
select_targets <- function(preprocessed_data,
                           targets = c("flavi", "dengue", 
                                       "zika", "chik",
                                       "dengue_zika",  "dengue_chik",
                                       "dengue_serotype", "dengue_serotype_neg"), 
                                       drop_original_target = TRUE,  min_samples = 1) {
  
  if (!"target" %in% names(preprocessed_data)) {
    stop('Column "Target" not found in preprocessed_data.')
  }
  
  # Validate all targets
  valid_targets <- c("flavi", "dengue", 
                     "zika", "chik",
                     "dengue_zika", 
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
    chik                = c(CHIKV = 1),
    dengue_zika         = c(DENV1 = 1, DENV2 = 1, DENV3 = 1, DENV4 = 1, ZIKV = 2),
    dengue_serotype     = c(DENV1 = 1, DENV2 = 2, DENV3 = 3, DENV4 = 4),
    dengue_serotype_neg = c(DENV1 = 1, DENV2 = 2, DENV3 = 3, DENV4 = 4),
    dengue_chik        = c(DENV1 = 1, DENV2 = 1, DENV3 = 1, DENV4 = 1, CHIKV = 2))
  
  # factor levels/labels
  label_specs <- list(
    flavi               = list(levels = c(0, 1),           labels = c("negative", "positive")),
    dengue              = list(levels = c(0, 1),           labels = c("negative", "positive")),
    zika                = list(levels = c(0, 1),           labels = c("negative", "positive")),
    chik                = list(levels = c(0, 1),           labels = c("negative", "positive")),
    dengue_zika         = list(levels = c(0, 1, 2),        labels = c("negative", "dengue", "zika")),
    dengue_serotype     = list(levels = c(1, 2, 3, 4),     labels = c("DENV1", "DENV2", "DENV3", "DENV4")),
    dengue_serotype_neg = list(levels = c(0, 1, 2, 3, 4),  labels = c("negative", "DENV1", "DENV2", "DENV3", "DENV4")),
    dengue_chik         = list(levels = c(1, 2),           labels = c("dengue", "CHIKV"))
  )
  
  
  # List to store results
  data_with_target <- list()
  
  # Process each target
  for (target in targets) {
    df_copy <- preprocessed_data
    
    # Build the chosen target column
    mp <- mapping_list[[target]]
    tgt_chr <- as.character(df_copy$target)
    
   # Map positives, everything else is negative (0)
    df_copy[[target]] <- ifelse(
      tgt_chr %in% names(mp),
      as.numeric(dplyr::recode(tgt_chr, !!!mp, .default = NA_real_)),
      0
    )
    
    # Remove negatives for classifications that only include infected samples (eg: given infection, classify dengue vs chik)
    if (target %in% c("dengue_chik", "dengue_serotype")) {
      df_copy <- df_copy[df_copy[[target]] %in% c(1, 2, 3, 4), ]  # keep only mapped positives
    }
    # convert to factor with labels
    spec <- label_specs[[target]]
    df_copy[[target]] <-
      factor(df_copy[[target]], 
             levels = spec$levels, 
             labels = spec$labels)

    # Drop classes with fewer than min_samples
    class_counts <- table(df_copy[[target]])
    valid_classes <- names(class_counts[class_counts >= min_samples])
    
    if (length(valid_classes) < length(class_counts)) {
      dropped <- names(class_counts[class_counts < min_samples])
      cat(sprintf("Target '%s': dropping class(es) with < %d samples: %s\n",
                  target, min_samples, paste(dropped, collapse = ", ")))
      df_copy <- df_copy[df_copy[[target]] %in% valid_classes, ]
      df_copy[[target]] <- droplevels(df_copy[[target]])  # remove unused factor levels
    }
    
    
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


train_multiple_targets_univariate <- function(
    data_list,  
    variables = NULL,  # Vector of variable names to test one by one
    k_fold = 5,
    metrics = c("AUROC", "AUPRC", "Brier", "StratBrier"),
    univariate = TRUE) {
  
  
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
    cat("Number of classes:", n_classes, "\n")
    
    current_data[[target_col]] <-  factor(current_data[[target_col]])
    
    # If variables not specified, use all columns except target
      if (is.null(variables)) {
      variables <- setdiff(names(current_data), target_col)
    }
    
    # Check for only one class
    if (n_classes < 2) {
      warning(paste("Target", target_col, "has only 1 class. Skipping."))
      next
    }
    
    # Initialize storage for this target
    target_comparisons <- list()
    target_predictions <- list()
    target_results <- list()
    
    # Loop through each variable for univariate analysis
    cat("\nTesting", length(variables), "variables individually \n\n")
    
    for (i in seq_along(variables)) {
      var <- variables[i]
      cat(sprintf("[%d/%d] Testing variable: %s\n", i, length(variables), var))
      
      # Check if variable exists in data
      if (!var %in% names(current_data)) {
        warning(paste("Variable", var, "not found in data. Skipping."))
        next
      }
      
      # Create a unique identifier for this variable-target combination
      result_name <- paste0(target_name, "_", var)
      
      tryCatch({
        # Train model with ONLY this variable
        if (n_classes == 2) {
          cat("  -> Calling train_binary_models\n")
          result <- train_binary_models(
            data = current_data,
            target = target_col,
            variables = var,  
            k_fold = k_fold,
            metrics = metrics
          )
        } else if (n_classes > 2) {
          result <- train_multinomial_models(
            data = current_data,
            target = target_col,
            variables = var, 
            n_repeats = 5
          )
        }
       
         # Store results with variable name
        target_results[[var]] <- result
        
        # Add metadata to comparison dataframe
        comparison_with_meta <- result$comparison
        comparison_with_meta$Target <- target_name
        comparison_with_meta$Variable <- var
        target_comparisons[[var]] <- comparison_with_meta
        
        # Add metadata to predictions dataframe
        predictions_with_meta <- result$predictions
        predictions_with_meta$Target <- target_name
        predictions_with_meta$Variable <- var
        target_predictions[[var]] <- predictions_with_meta
        
        cat(sprintf("Completed"))
        
      }, error = function(e) {
        cat("  -> ERROR:", e$message, "\n") 
        warning(paste("Error processing variable", var, "for target", target_name, ":", e$message))
      })
    }
    
    # Store results for this target
    all_results[[target_name]] <- target_results
    
    # Combine comparisons and predictions for this target
    if (length(target_comparisons) > 0) {
      all_comparisons[[target_name]] <- dplyr::bind_rows(target_comparisons)
    }
    if (length(target_predictions) > 0) {
      all_predictions[[target_name]] <- dplyr::bind_rows(target_predictions)
    }
  }
  
  # Combine all comparison dataframes across targets
  combined_comparison <- NULL
  if (length(all_comparisons) > 0) {
    combined_comparison <- dplyr::bind_rows(all_comparisons)
    rownames(combined_comparison) <- NULL
    # Reorder columns to put Target and Variable first
    col_order <- c("Target", "Variable", setdiff(names(combined_comparison), c("Target", "Variable")))
    combined_comparison <- combined_comparison[, col_order]
  }
  
  # Combine all predictions dataframes across targets
  combined_predictions <- NULL
  if (length(all_predictions) > 0) {
    combined_predictions <- dplyr::bind_rows(all_predictions)
    rownames(combined_predictions) <- NULL
    # Reorder columns to put Target and Variable first
    pred_col_order <- c("Target", "Variable", setdiff(names(combined_predictions), c("Target", "Variable")))
    combined_predictions <- combined_predictions[, pred_col_order]
  }
  
  return(list(
    results_by_target = all_results,
    combined_comparison = combined_comparison,
    combined_predictions = combined_predictions,
    summary = list(
      n_targets = length(all_results),
      target_names = names(all_results),
      n_variables_tested = length(variables),
      variables_tested = variables,
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
  min_samples = 2
)


data_with_binomial_targets_chik <- select_targets(
  preprocessed_data = ratio_df,
  targets = c("chik", "dengue_chik"),
  drop_original_target = FALSE,
  min_samples = 2
)


sum(is.na(data_with_binomial_targets_chik$chik$chik))
sum(is.na(data_with_binomial_targets_chik$dengue_chik$dengue_chik))
sum(is.na(data_with_binomial_targets$flavi$flavi))
sum(is.na(data_with_binomial_targets$flavi$dengue))
sum(is.na(data_with_binomial_targets$flavi$zika))


# View distribution of targets
table(data_with_binomial_targets$flavi$flavi)
table(data_with_binomial_targets$dengue$dengue)
table((data_with_binomial_targets_chik$dengue_chik$dengue_chik))
table((data_with_binomial_targets_chik$chik$chik))


# Fit binomial models
binomial_modeling_results <- train_multiple_targets(
  data_list = data_with_binomial_targets,
  variables = NULL,  # Uses all columns except target
  k_fold = 5,
  metrics = c("AUROC", "AUPRC", "Brier"))


# Fit binomial model - for dengue vs chik 
binomial_modeling_results_dengue_chik <- train_multiple_targets(
  data_list = data_with_binomial_targets_chik,
  variables = NULL,  # Uses all columns except target
  k_fold = 5,
  metrics = c("AUROC", "AUPRC", "Brier")
)


# results
binomial_modeling_results$combined_comparison
binomial_modeling_results_dengue_chik$combined_comparison

# save 
saveRDS(binomial_modeling_results, 'Results/binomial_modeling_results.rds')
saveRDS(binomial_modeling_results_dengue_chik, 'Results/binomial_modeling_results_dengue_chik.rds')


# --- Multinomial Results 
data_with_multinomial_targets <- select_targets(
  preprocessed_data = ratio_df,
  targets = c("dengue_serotype", "dengue_serotype_neg"),
  drop_original_target = TRUE, 
  min_samples = 2
)

sum(is.na(data_with_multinomial_targets$dengue_serotype$dengue_serotype))
sum(is.na(data_with_multinomial_targets$dengue_serotype_neg$dengue_serotype_neg))


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


# Micro-average "takes imbalance into account" in the sense that the resulting performance is based on the proportion of every class
# i.e.the performance of a large class has more impact on the result than of a small class.
# Macro-average "doesn't take imbalance into account" in the sense that the resulting performance 
# is a simple average over the classes, so every class is given equal weight independently from their proportion.

# Reporting Micro Average AUC 

# ---- Best from binomial results ----
best_binomial <- bind_rows(
  binomial_modeling_results$combined_comparison,
  binomial_modeling_results_dengue_chik$combined_comparison
) %>%
  group_by(target) %>%
  slice_max(AUROC, n = 1) %>%
  ungroup() %>%
  rename(AUC = AUROC) %>%
  mutate(type = "binary")


# ---- Best from multinomial results ----
best_multinomial <- dengue_serotype_results$combined_comparison %>%
  group_by(target) %>%
  slice_max(AUC_Micro, n = 1) %>%
  ungroup() %>%
  rename(AUC = AUC_Micro) %>%
  mutate(type = "multinomial")


# ---- Combine into summary table ----
best_models <- bind_rows(
  best_binomial %>% dplyr::select(type, target, Model, AUC, AUPRC, Brier),
  best_multinomial %>% dplyr::select(type, target, Model, AUC = AUC, 
                               AUPRC = AUPRC_Micro, Brier)
) %>%
  arrange(type, target)

print(best_models)



# Univariate Analysis - look at each antigen independently 
serotype_variables_flavi <- grep("DENV|ZIKV", names(ratio_df), value = TRUE)
serotype_variables_alpha <- grep("CHIKV|ONNV|MAYV", names(ratio_df), value = TRUE)


univariate_results_flavi <- train_multiple_targets_univariate(
  data_list  = data_with_binomial_targets,
  variables  = serotype_variables_flavi,
  k_fold     = 5,
  metrics    = c("AUROC", "AUPRC", "Brier"),
)

univariate_results_alpha <- train_multiple_targets_univariate(
  data_list  = data_with_binomial_targets_chik,
  variables  = serotype_variables_alpha,
  k_fold     = 5,
  metrics    = c("AUROC", "AUPRC", "Brier"),
)

# save 
saveRDS(univariate_results_flavi, 'Results/univariate_results_flavi.rds')
saveRDS(univariate_results_alpha, 'Results/univariate_results_alpha.rds')



model_colours <- c(
  "GLMnet" = "#012b48",
  "Decision Tree" = "#0396f8",
  "SVM" = "#de5a7b"
)

plot_df_dengue <- univariate_results_flavi$combined_comparison %>%
  filter(Target == "dengue") %>%
  dplyr::select(Target, Variable, Model, AUROC) %>%
  pivot_longer(cols = AUROC, names_to = "Metric", values_to = "Value") %>%
  group_by(Variable) %>%
  mutate(mean_auc = mean(Value, na.rm = TRUE)) %>%
  ungroup()

var_order_dengue <- plot_df_dengue %>%
  distinct(Variable, mean_auc) %>%
  arrange(desc(mean_auc)) %>%
  mutate(
    Variable_clean = Variable %>%
      str_replace("SHERPADES_([^_]+)_DIII", "SHERPADES_\\1") %>%  # removes DIII only after SHERPADES
      str_replace("SHERPADES_", "SHERPADES ") %>%
      str_replace_all("_", " ") %>%
      str_wrap(width = 14)
  ) %>%
  mutate(
    Variable_clean = factor(Variable_clean, levels = Variable_clean)
  )

# Join labels back and explicitly set factor order
plot_df_dengue <- plot_df_dengue %>%
  left_join(
    var_order_dengue %>% dplyr::select(Variable, Variable_clean),
    by = "Variable"
  ) %>%
  mutate(
    Variable_clean = factor(
      Variable_clean,
      levels = var_order_dengue$Variable_clean
    )
  )


univariate_plot_dengue <- ggplot(
  plot_df_dengue,
  aes(x = Value, y = Variable_clean, color = Model, shape = Model)
) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "#8f0000") +
  geom_point(size = 7, alpha = 0.9) +
  facet_grid(Metric ~ ., scales = "free_y") +  # flipped facet direction
  scale_x_continuous(limits = c(0, 1)) +
  coord_flip() +
  scale_color_manual(values = model_colours) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 18), 
    axis.text.x = element_text(size = 18, angle = 40, hjust = 1),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.background = element_rect(fill = "white", colour = NA),  # white background
    plot.background  = element_rect(fill = "white", colour = NA),  # white outer bg
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(colour = "grey85"),
    legend.title = element_blank(),
    legend.text = element_text(size = 20),
    strip.text = element_blank(),
    legend.position = "bottom",
    plot.margin = margin(t = 10, r = 22, b = 10, l = 10)
  ) 

ggsave("Results/univariate_plot_dengue.png",
       univariate_plot_dengue, width = 20, height = 10, dpi = 300)




# CHIK Univariate results 
plot_df_alpha <- univariate_results_alpha$combined_comparison %>%
  filter(Target == "chik") %>%
  dplyr::select(Target, Variable, Model, AUROC) %>%
  pivot_longer(cols = AUROC, names_to = "Metric", values_to = "Value") %>%
  group_by(Variable) %>%
  mutate(mean_auc = mean(Value, na.rm = TRUE)) %>%
  ungroup()

var_order_alpha <- plot_df_alpha %>%
  distinct(Variable, mean_auc) %>%
  arrange(desc(mean_auc)) %>%
  mutate(
    Variable_clean = Variable %>%
      str_replace("SHERPADES_([^_]+)_DIII", "SHERPADES_\\1") %>%  # removes DIII only after SHERPADES
      str_replace("SHERPADES_", "SHERPADES ") %>%
      str_replace_all("_", " ") %>%
      str_wrap(width = 14)
  ) %>%
  mutate(
    Variable_clean = factor(Variable_clean, levels = Variable_clean)
  )

# Join labels back and explicitly set factor order
plot_df_alpha <- plot_df_alpha %>%
  left_join(
    var_order_alpha %>% dplyr::select(Variable, Variable_clean),
    by = "Variable"
  ) %>%
  mutate(
    Variable_clean = factor(
      Variable_clean,
      levels = var_order_alpha$Variable_clean
    )
  )


univariate_plot_alpha <- ggplot(
  plot_df_alpha,
  aes(x = Value, y = Variable_clean, color = Model, shape = Model)
) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "#8f0000") +
  geom_point(size = 7, alpha = 0.9) +
  facet_grid(Metric ~ ., scales = "free_y") +  # flipped facet direction
  scale_x_continuous(limits = c(0, 1)) +
  coord_flip() +
  scale_color_manual(values = model_colours) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 18), 
    axis.text.x = element_text(size = 18, angle = 40, hjust = 1),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.background = element_rect(fill = "white", colour = NA),  # white background
    plot.background  = element_rect(fill = "white", colour = NA),  # white outer bg
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(colour = "grey85"),
    legend.title = element_blank(),
    legend.text = element_text(size = 20),
    strip.text = element_blank(),
    legend.position = "bottom",
    plot.margin = margin(t = 10, r = 22, b = 10, l = 10)
  ) 

ggsave("Results/univariate_plot_alpha.png",
       univariate_plot_alpha, width = 20, height = 10, dpi = 300)


