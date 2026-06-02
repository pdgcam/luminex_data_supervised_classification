# --- Source functions
source(here('Models/train_binary_models.R'))
source(here('Models/train_multinomial_models.R'))
source(here('Functions.R'))

# function for univariate analysis across multiple targets
train_multiple_targets_univariate <- function(
    data_list,  
    variables = NULL,
    positive_class_map = list(),
    univariate = FALSE, 
    metrics = c("ROC", "AUPRC", "Brier", "StratBrier")) {
  
  # Store all results
  all_results <-  list()
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
    
    # Initialize storage for this target
    all_results[[target_name]] <- list()
    all_comparisons[[target_name]] <- list()
    all_predictions[[target_name]] <- list()
    
    # Loop through each variable for univariate analysis
    cat("\nTesting", length(variables), "variables individually \n\n")
    print(head(variables))

    pos_class <- positive_class_map[[target_name]]
    
    for (i in seq_along(variables)) {
      var <- variables[i]
      cat(sprintf("[%d/%d] Testing variable: %s\n", i, length(variables), var))
      
      # Check if variable exists in data
      if (!var %in% names(current_data)) {
        warning(paste("Variable", var, "not found in data. Skipping."))
        next
      }
      
        if (n_classes == 2) {
          result <- train_binary_models(
            data = current_data,
            target = target_col,
            variables = var,
            positive_class = pos_class,
            metrics = metrics
          )

          # Store results
          all_results[[target_name]][[var]] <- result
          
          
          # Add target name to comparison dataframe
          comparison_with_target <- result$comparison
          comparison_with_target$target <- target_name
          comparison_with_target$variable <- var
          all_comparisons[[target_name]][[var]] <- comparison_with_target
          
          # Add target name to predictions dataframe
          predictions_with_target <- result$predictions
          predictions_with_target$target <- target_name
          predictions_with_target$variable <- var
          all_predictions[[target_name]][[var]] <- predictions_with_target

        } else if (n_classes > 2) {
          result <- train_multinomial_models(
            data = current_data,
            target = target_col,
            variables = var,
            univariate = TRUE,
            metrics = metrics
          )
          all_results[[target_name]][[var]] <- result
          
          
          # Add target name to comparison dataframe
          comparison_with_target <- result$comparison
          comparison_with_target$target <- target_name
          comparison_with_target$variable <- var
          all_comparisons[[target_name]][[var]] <- comparison_with_target
          
          # Add target name to predictions dataframe
          predictions_with_target <- result$predictions
          predictions_with_target$target <- target_name
          predictions_with_target$variable <- var
          all_predictions[[target_name]][[var]] <- predictions_with_target
        
      }
    } }  
      
      # Combine all comparison dataframes
      combined_comparison  <- dplyr::bind_rows(unlist(all_comparisons,  recursive = FALSE))
      rownames(combined_comparison) <- NULL
      col_order <- c("target", "variable", setdiff(names(combined_comparison), c("target", "variable")))
      combined_comparison <- combined_comparison[, col_order]
      
      # Combine all predictions dataframes
      combined_predictions <- dplyr::bind_rows(unlist(all_predictions, recursive = FALSE))
      rownames(combined_predictions) <- NULL
      pred_col_order <- c("target", "variable", setdiff(names(combined_predictions), c("target", "variable")))
      combined_predictions <- combined_predictions[, pred_col_order]
      
      return(list(
        results_by_target = all_results,
        combined_comparison = combined_comparison,
        combined_predictions = combined_predictions,
        summary = list(
          n_targets = length(all_results),
          target_names = names(all_results),
          n_variables_tested = length(variables),
          variables_tested = variables,
          metrics = metrics
        )
      ))
    
}



serotype_variables_flavi <- grep("DENV|ZIKV", names(ratio_df), value = TRUE)
serotype_variables_alpha <- grep("CHIKV|ONNV|MAYV", names(ratio_df), value = TRUE)


# Define targets / classification question
data_with_binomial_targets <- select_targets(
  preprocessed_data = ratio_df_logged,
  targets = c("flavi", "dengue"),
  drop_original_target = TRUE,
  min_samples = 2
)

univariate_results_flavi <- train_multiple_targets_univariate(
  data_list  = data_with_binomial_targets$data,
  variables  = serotype_variables_flavi,
  metrics    = c("ROC", "AUPRC", "Brier"),
  positive_class_map = list(flavi = "positive", 
                            dengue = "positive") 
)


univariate_results_serotype <- train_multiple_targets_univariate(
  data_list  = data_with_multinomial_targets$data,
  variables  = serotype_variables_flavi,
  metrics    = c("ROC", "AUPRC", "Brier"),
  univariate = TRUE
)



univariate_results_alpha <- train_multiple_targets_univariate(
  data_list  = data_with_binomial_targets_chik,
  variables  = serotype_variables_alpha,
  metrics    = c("ROC", "AUPRC", "Brier"),
  positive_class_map = list(chik = "positive", 
                            dengue_chik = "dengue")   
)


# save 
saveRDS(univariate_results_flavi, 'Results/univariate_results_flavi.rds')
saveRDS(univariate_results_alpha, 'Results/univariate_results_alpha.rds')

univariate_results_flavi <- readRDS('Results/univariate_results_flavi.rds')
univariate_results_alpha <- readRDS('Results/univariate_results_alpha.rds')


model_colours <- c(
  "GLM" = "#012b48",
  "Decision Tree" = "#0396f8",
  "SVM" = "#de5a7b", 
  "NaiveBayes" = "#ed7f01"
)


plot_df_dengue <- univariate_results_flavi$combined_comparison %>%
  filter(target == "dengue") %>%
  dplyr::select(target, variable, Model, ROC, ROC_low, ROC_high) %>%
  pivot_longer(cols = ROC, names_to = "Metric", values_to = "Value") %>%
  group_by(variable) %>%
  mutate(mean_auc = mean(Value, na.rm = TRUE)) %>%
  ungroup()

var_order_dengue <- plot_df_dengue %>%
  distinct(variable, mean_auc) %>%
  arrange(desc(mean_auc)) %>%
  mutate(
    variable_clean = variable %>%
      str_replace("SHERPADES_([^_]+)_DIII", "SHERPADES_\\1") %>%  # removes DIII only after SHERPADES
      str_replace("SHERPADES_", "SHERPADES ") %>%
      str_replace_all("_", " ") %>%
      str_wrap(width = 14)
  ) %>%
  mutate(
    variable_clean = factor(variable_clean, levels = variable_clean)
  )

# Join labels back and explicitly set factor order
plot_df_dengue <- plot_df_dengue %>%
  left_join(
    var_order_dengue %>% dplyr::select(variable, variable_clean),
    by = "variable"
  ) %>%
  mutate(
    variable_clean = factor(
      variable_clean,
      levels = var_order_dengue$variable_clean
    )
  )

univariate_plot_dengue <- ggplot(
  plot_df_dengue,
  aes(x = Value, y = variable_clean, color = Model)
) +
  geom_vline(xintercept = 0.99, linetype = "dashed", colour = "#8f0000") +
  geom_errorbarh(
    aes(xmin = ROC_low, xmax = ROC_high),
    height = 0.2,
    alpha = 0.6,
    linewidth = 0.8,
    position = position_dodge(width = 0.8)
  ) + 
  geom_point(size = 7, alpha = 0.9,
  position = position_dodge(width = 0.8)) +
  facet_grid(Metric ~ ., scales = "free_y") +
  scale_x_continuous(limits = c(0, 1)) +
  coord_flip() +
  scale_color_manual(values = model_colours) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 18), 
    axis.text.x = element_text(size = 18, angle = 90, hjust = 1),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.background = element_rect(fill = "white", colour = NA),  
    plot.background  = element_rect(fill = "white", colour = NA), 
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(colour = "grey85"),
    legend.title = element_blank(),
    legend.text = element_text(size = 20),
    strip.text = element_blank(),
    legend.position = "bottom",
    plot.margin = margin(t = 10, r = 10, b = 10, l = 10)
  )

print(univariate_plot_dengue) 

ggsave("Results/NEW_univariate_plot_dengue.png",
       univariate_plot_dengue, width = 20, height = 10, dpi = 300)




# CHIK Univariate results 
plot_df_alpha <- univariate_results_alpha$combined_comparison %>%
  filter(target == "chik") %>%
  dplyr::select(target, variable, Model, ROC) %>%
  pivot_longer(cols = ROC, names_to = "Metric", values_to = "Value") %>%
  group_by(variable) %>%
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
  aes(x = Value, y = Variable_clean, color = Model)
) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "#8f0000") +
  geom_errorbarh(
    aes(xmin = ROC_low, xmax = ROC_high),
    height = 0.2,
    alpha = 0.6,
    linewidth = 0.8
  ) +
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

unique(univariate_results_serotype$combined_comparison$Model)



View(univariate_results_serotype$combined_comparison)

# dengue serotype univariate plot 
plot_df_dengue_serotype <- univariate_results_serotype$combined_comparison %>%
  filter(target == "dengue_serotype") %>%
  dplyr::select(target, variable, Model, AUC_Micro) %>%
  pivot_longer(cols = AUC_Micro, names_to = "Metric", values_to = "Value") %>%
  group_by(variable) %>%
  mutate(mean_auc = mean(Value, na.rm = TRUE)) %>%
  ungroup()

var_order_dengue_serotype <- plot_df_dengue_serotype %>%
  distinct(variable, mean_auc) %>%
  arrange(desc(mean_auc)) %>%
  mutate(
    variable_clean = variable %>%
      str_replace("SHERPADES_([^_]+)_DIII", "SHERPADES_\\1") %>%  # removes DIII only after SHERPADES
      str_replace("SHERPADES_", "SHERPADES ") %>%
      str_replace_all("_", " ") %>%
      str_wrap(width = 14)
  ) %>%
  mutate(
    variable_clean = factor(variable_clean, levels = variable_clean)
  )

# Join labels back and explicitly set factor order
plot_df_dengue_serotype <- plot_df_dengue_serotype %>%
  left_join(
    var_order_dengue_serotype %>% dplyr::select(variable, variable_clean),
    by = "variable"
  ) %>%
  mutate(
    variable_clean = factor(
      variable_clean,
      levels = var_order_dengue_serotype$variable_clean
    )
  )

univariate_plot_dengue_serotype <- ggplot(
  plot_df_dengue_serotype,
  aes(x = Value, y = variable_clean, color = Model)
) +
  geom_vline(xintercept = 0.99, linetype = "dashed", colour = "#8f0000") +
  geom_point(size = 7, alpha = 0.9,
  position = position_dodge(width = 1.2)) +
  facet_grid(Metric ~ ., scales = "free_y") +
  scale_x_continuous(limits = c(0, 1)) +
  coord_flip() +
  scale_color_manual(values = model_colours) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 18), 
    axis.text.x = element_text(size = 18, angle = 90, hjust = 1),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.background = element_rect(fill = "white", colour = NA),  
    plot.background  = element_rect(fill = "white", colour = NA), 
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(colour = "grey85"),
    legend.title = element_blank(),
    legend.text = element_text(size = 20),
    strip.text = element_blank(),
    legend.position = "bottom",
    plot.margin = margin(t = 10, r = 22, b = 10, l = 10)
  )
print(univariate_plot_dengue_serotype)


ggsave("Results/NEW_univariate_plot_dengue_serotype.png",
       univariate_plot_dengue_serotype, width = 20, height = 10, dpi = 300)




# Forward Selection - univariate models
forward_stepwise_selection <- function(
    univariate_results,   # output from train_multiple_targets_univariate
    data_list,
    positive_class_map = list(),
    metrics = c("ROC")) {
  
  all_stepwise <- list()
  
  for (target_name in names(univariate_results$results_by_target)) {
    cat("Forward selection for target:", target_name, "\n")
    
    current_data <- data_list[[target_name]]
    pos_class    <- positive_class_map[[target_name]]
    
    # ---- Rank variables by univariate AUC (best first) ----
    target_comparisons <- univariate_results$combined_comparison %>%
      filter(target == target_name) %>%
      arrange(desc(ROC))
    
    ranked_vars <- target_comparisons$variable
    
    # ---- Greedy forward steps ----
    selected  <- c()
    step_aucs <- c()
    
    for (var in ranked_vars) {
      selected <- c(selected, var)
      cat(sprintf("  Step %d: %s\n", length(selected), paste(selected, collapse = " + ")))
      
      result <- train_binary_models(
        data          = current_data,
        target        = target_name,
        variables     = selected,
        positive_class = pos_class,
        metrics       = metrics
      )
      
      # Pull GLM AUC (most stable for stepwise)
      auc_key  <- if (length(selected) == 1) "glm" else "glmnet"
      step_auc <- result$pooled_aucs[[auc_key]]$auc
      step_aucs <- c(step_aucs, step_auc)
    }
    
    # ---- Compile results ----
    all_stepwise[[target_name]] <- data.frame(
      step     = seq_along(ranked_vars),
      variable = ranked_vars,
      AUC      = step_aucs,
      gain     = c(NA, diff(step_aucs))
    )
  }
  
  return(all_stepwise)
}

plot_stepwise_results <- function(stepwise_results) {
  
  library(ggplot2)
  
  # Combine all targets into one df
  combined <- dplyr::bind_rows(
    lapply(names(stepwise_results), function(t) {
      stepwise_results[[t]]$target <- t
      stepwise_results[[t]]
    })
  )
  
  ggplot(combined, aes(x = step, y = AUC)) +
    geom_line(colour = "steelblue", linewidth = 0.8) +
    geom_point(colour = "steelblue", size = 2.5) +
    geom_text(aes(label = variable), angle = 30, hjust = 0, vjust = -0.8, size = 2.8) +
    scale_x_continuous(breaks = seq_len(max(combined$step))) +
    scale_y_continuous(limits = c(NA, 1.05)) +
    facet_wrap(~ target, scales = "free_x") +
    labs(x = "Step", y = "AUC", title = "Forward stepwise selection") +
    theme_minimal()
}

stepwise_results <- forward_stepwise_selection(
  univariate_results = univariate_results_flavi,
  data_list          = data_with_binomial_targets$data,
   positive_class_map = list(flavi = "positive", 
                            dengue = "positive") 
)

plot_stepwise_results(stepwise_results)


names(univariate_results_flavi)
