

train_binary_models <- function(
    data,
    target,
    variables = NULL,
    metrics = c("AUROC", "AUPRC", "Brier"),
    univariate = FALSE) {
  
  # ---- Input Validation ----
  if (!target %in% names(data)) {
    stop(paste("Target column", target, "not found in data"))
  }
  

  if (is.null(variables)) {
    variables <- setdiff(names(data), target)
  } else {
    missing_vars <- setdiff(variables, names(data))
    if (length(missing_vars) > 0) {
      stop(paste("Variables not found in data:", paste(missing_vars, collapse = ", ")))
    }
  }
  
  # --- get data to input to model ----
  model_data <- data[, c(variables, target)]
  n_predictors <- ncol(model_data) - 1
  
  
  valid_metrics <- c("AUROC", "AUPRC", "Brier", "StratBrier")
  invalid_metrics <- setdiff(metrics, valid_metrics)
  if (length(invalid_metrics) > 0) {
    stop(paste("Invalid metrics:", paste(invalid_metrics, collapse = ", "),
               "\nValid options:", paste(valid_metrics, collapse = ", ")))
  }

  
  # Ensure target is a factor
  if (!is.factor(model_data[[target]])) {
    model_data[[target]] <- factor(model_data[[target]])
  }

  if (!"positive" %in% levels(model_data[[target]])) {
    stop(paste("'positive' level not found in target — levels are:",
               paste(levels(model_data[[target]]), collapse = ", ")))
  }

  # expects positive class to be first level 
  model_data[[target]] <- relevel(model_data[[target]], ref = "positive")
  target_levels <- levels(model_data[[target]])  
  cat(sprintf("Binary classification: %s (class 1) vs %s (class 2)\n",
            target_levels[1], target_levels[2]))

  
  # ---- Set Up Train Control ----
  combinedSummary <- function(data, lev = NULL, model = NULL) {
    stopifnot(length(lev) == 2)
    
    # AUROC, Sens, Spec from twoClassSummary
    base <- caret::twoClassSummary(data, lev = lev, model = model)
    
    # lev[1] == positive class
    pos_name <- lev[1]
    neg_name <- lev[2] 
    p_pos <- data[[pos_name]]
    y_true <- as.integer(data$obs == pos_name)
    
    # AUPRC
    auprc <- tryCatch(
      MLmetrics::PRAUC(p_pos, y_true),
      error = function(e) NA_real_
    )
    
    # Standard Brier score
    brier <- mean((p_pos - y_true)^2)

    c(base, AUPRC = auprc, Brier = brier)
  }
  
  binary_control <- caret::trainControl(
    method  = "LOOCV", 
    summaryFunction = combinedSummary,
    classProbs = TRUE,
    verboseIter = FALSE,
    savePredictions = "all")

  # ---- Helper: pooled AUC from held-out predictions ----
  get_pooled_auc <- function(model, model_data, target) {
    pos_name    <- levels(model_data[[target]])[1]  # "positive"
    preds       <- model$pred
    best_tune   <- model$bestTune
    
    # Filter to best tuning parameters (glm has no tuning params — just take all)
    filter_mask <- rep(TRUE, nrow(preds))
    for (param in names(best_tune)) {
      if (param %in% names(preds)) {
        filter_mask <- filter_mask & (preds[[param]] == best_tune[[param]])
      }
    }
    pooled_best <- preds[filter_mask, ]
    
    roc_obj <- pROC::roc(
      response  = pooled_best$obs,
      predictor = pooled_best[[pos_name]],
      quiet     = TRUE
    )
    ci <- pROC::ci.auc(roc_obj, conf.level = 0.95)
    list(auc = as.numeric(ci[2]), ci_low = as.numeric(ci[1]), ci_high = as.numeric(ci[3]))
  }
  
  # ---- Train Models ----

  # --- GLM ---
  if (univariate || n_predictors == 1) {
    cat("Training GLM\n")
    glm_model <- caret::train(
      as.formula(paste(target, "~ .")),
      data = model_data,
      metric = "ROC",
      method = "glm",
      trControl = binary_control
    )
  } else {
    cat("Training GLM with elastic net\n")
    glm_model <- caret::train(
    as.formula(paste(target, "~ .")),
    data = model_data,
    metric = "ROC",
    method = "glmnet",
    tuneGrid = expand.grid(
    alpha = seq(0, 1, length.out = 5),
    lambda = 10^seq(-3, 1, length.out = 20)
    ),
  trControl = binary_control
  )
}
  
  # ---- Random Forest or Decision Tree ----
  # mtry values tuned to number of predictors 
  # mtry cant be > predictors 
  n_features <- length(variables)
  mtry_values <- unique(pmax(1, floor(c(
    sqrt(n_features),           
    n_features / 3,
    n_features / 2,
    n_features
  ))))
  mtry_values <- mtry_values[mtry_values <= n_features]
  
  if (univariate || n_predictors == 1) {
    cat("Training Decision Tree\n")
    tree_model <- caret::train(
      as.formula(paste(target, "~ .")),
      data = model_data,
      metric = "ROC",
      method = "rpart",
      tuneGrid = expand.grid(
        cp = c(0, 10^seq(-4, -1, length.out = 20))
      ),
      trControl = binary_control
    )
  } else {
    cat("Training Random Forest\n")
    rf_model <- caret::train(
      as.formula(paste(target, "~ .")),
      data = model_data,
      metric = "ROC",
      method = "ranger",
      tuneGrid = expand.grid(
        mtry = mtry_values,
        splitrule = c("gini", "extratrees"),
        min.node.size = c(1, 5, 10)),
      trControl = binary_control
    )
  }
  
  cat("Training SVM\n")
  svm_model <- caret::train(
    as.formula(paste(target, "~ .")),
    data = model_data,
    metric = "ROC",
    method = "svmRadial",
    tuneLength = 10,
    trControl = binary_control
  )

  # ---- Compute Pooled AUCs ----
  pooled_aucs <- list()
  pooled_aucs$glmnet <- get_pooled_auc(glm_model, model_data, target)
  if (univariate || n_predictors == 1) {
    pooled_aucs$tree <- get_pooled_auc(tree_model, model_data, target)
  } else {
    pooled_aucs$rf <- get_pooled_auc(rf_model, model_data, target)
  }
  pooled_aucs$svm <- get_pooled_auc(svm_model, model_data, target)

  # ---- Extract Predictions ----
  glm_preds <- glm_model$pred %>%
  mutate(Model = if (univariate || n_predictors == 1) "GLM" else "GLMnet")

  if (univariate || n_predictors == 1) {
    tree_preds <- tree_model$pred %>%
      mutate(Model = "Decision Tree")
  } else {
    tree_preds <- rf_model$pred %>%
      mutate(Model = "Random Forest")
  }

  svm_preds <- svm_model$pred %>%
    mutate(Model = "SVM")


  filter_to_best <- function(preds, best_tune) {
  mask <- rep(TRUE, nrow(preds))
  for (param in names(best_tune)) {
    if (param %in% names(preds)) {
      mask <- mask & (preds[[param]] == best_tune[[param]])
    }
  }
  preds[mask, ]
}

  all_predictions <- bind_rows(
  filter_to_best(glm_preds,  glm_model$bestTune),
  filter_to_best(tree_preds, if (univariate || n_predictors == 1) tree_model$bestTune else rf_model$bestTune),
  filter_to_best(svm_preds,  svm_model$bestTune)
)

  # Rename columns
  colnames(all_predictions)[colnames(all_predictions) == "obs"] <- "obs_class"
  colnames(all_predictions)[colnames(all_predictions) == "pred"] <- "pred_class"
  
  for (level in target_levels) {
    new_col <- paste0("prob_", level)
    if (level %in% colnames(all_predictions)) {
    colnames(all_predictions)[colnames(all_predictions) == level] <- new_col
  }
  }
  # Reorder columns for consistency
  prob_cols <- paste0("prob_", target_levels)
  all_predictions <- all_predictions %>%
    dplyr::select(Model, rowIndex, obs_class, pred_class, all_of(prob_cols))

  
  # ---- Build Model List ----
  model_list <- list(glmnet = glm_model, svm = svm_model)

  if (univariate || n_predictors == 1) {
    model_list$tree  <- tree_model
    display_names    <- c(if (univariate || n_predictors == 1) "GLM" else "GLMnet",
                          "Decision Tree", "SVM")
    model_names      <- c("glmnet", "tree", "svm")
  } else {
    model_list$rf    <- rf_model
    display_names    <- c("GLMnet", "Random Forest", "SVM")
    model_names      <- c("glmnet", "rf", "svm")
  }

  # ---- Compile Metrics using Pooled AUC ----
  # pooled_aucs already computed by get_pooled_auc() earlier in the function
  comparison_df <- data.frame(
    Model     = display_names,
    AUROC     = sapply(model_names, function(m) pooled_aucs[[m]]$auc),
    AUROC_low = sapply(model_names, function(m) pooled_aucs[[m]]$ci_low),
    AUROC_high= sapply(model_names, function(m) pooled_aucs[[m]]$ci_high)
  )
  
  # Calculate other metrics from model results
  for (metric in setdiff(metrics, "AUROC")) {
    comparison_df[[metric]] <- sapply(model_names, function(model) {
      model_obj <- model_list[[model]]
      if (length(model_obj$bestTune) == 0 || all(is.na(model_obj$bestTune))) {
        # No tuning parameters (e.g., GLM)
        model_obj$results[1, metric]
      } else {
        # Find row matching bestTune
        best_idx <- which(apply(model_obj$results[, names(model_obj$bestTune), drop = FALSE], 1, 
                                function(row) all(row == model_obj$bestTune)))
        if (length(best_idx) == 0) best_idx <- 1  # fallback
        model_obj$results[best_idx[1], metric]
      }
    })
  }
  
  # ---- Return Results ----
  return(list(
    models         = model_list,
    predictions    = all_predictions,
    comparison     = comparison_df,
    pooled_aucs    = pooled_aucs,
    variables_used = variables,
    target_used    = target,
    target_levels  = target_levels, 
    method         = "LOOCV",
    metrics        = metrics,
    model_type = if(univariate || n_predictors == 1) "tree" else "rf"
  ))
}




# Function to run binomial and multinomial classification with multiple targets simulatensouly 
train_multiple_targets <- function(
    data_list,  
    variables = NULL,
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
      metrics = metrics
    )
  ))
}

