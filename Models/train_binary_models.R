
train_binary_models <- function(
    data,
    target,
    variables = NULL,
    k_fold = 5,
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
  
  if (!k_fold %in% c(5, 10)) {
    warning("k_fold should be 5 or 10. Using provided value.")
  }
  
  # Ensure target is a factor
  if (!is.factor(model_data[[target]])) {
    model_data[[target]] <- factor(model_data[[target]])
  }
  
  # Get level names for later use
  target_levels <- levels(model_data[[target]])
  cat(sprintf("Binary classification: %s (class 1) vs %s (class 2)\n", 
              target_levels[1], target_levels[2]))
  
 
  # ---- Create Folds ----
  cv_folds <- caret::createFolds(model_data[[target]], k = k_fold)
  train_indices <- lapply(cv_folds, function(test_idx) {
    setdiff(seq_len(nrow(model_data)), test_idx)
  })
  
  for (i in seq_along(cv_folds)) {
    test_idx <- cv_folds[[i]]
    tbl <- table(model_data[[target]][test_idx])
    cat(sprintf("Fold %d assessment class counts: %s\n", i,
                paste(sprintf("%s=%d", names(tbl), as.integer(tbl)), collapse = ", ")))
  }
  
  # ---- Set Up Train Control ----
  combinedSummary <- function(data, lev = NULL, model = NULL) {
    stopifnot(length(lev) == 2)
    
    # AUROC, Sens, Spec from twoClassSummary
    base <- caret::twoClassSummary(data, lev = lev, model = model)
    
    # lev[2] == positive class
    pos_name <- lev[2]
    neg_name <- lev[1]
    
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
    summaryFunction = combinedSummary,
    classProbs = TRUE,
    verboseIter = FALSE,
    savePredictions = "final",
    index = train_indices
  )
  
  # ---- Train Models ----
  if (univariate || n_predictors == 1) {
    cat("Training GLM\n")
    glm_model <- caret::train(
      as.formula(paste(target, "~ .")),
      data = model_data,
      metric = "AUROC",
      method = "glm",
      trControl = binary_control
    )
  } else {
    cat("Training GLM with elastic net\n")
    glm_model <- caret::train(
    as.formula(paste(target, "~ .")),
    data = model_data,
    metric = "AUROC",
    method = "glmnet",
    tuneGrid = expand.grid(
    alpha = seq(0, 1, length.out = 5),
   lambda = 10^seq(-3, 1, length.out = 20)
    ),
  trControl = binary_control
  )
}
  
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
  
  # ---- Train Random Forest or Decision Tree ----
  if (univariate || n_predictors == 1) {
    cat("Training Decision Tree\n")
    tree_model <- caret::train(
      as.formula(paste(target, "~ .")),
      data = model_data,
      metric = "AUROC",
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
      metric = "AUROC",
      method = "ranger",
      tuneGrid = expand.grid(
        mtry = mtry_values,
        splitrule = c("gini", "extratrees"),
        min.node.size = c(1, 5, 10)),
      trControl = binary_control
    )
  }
  
  cat("Training SVM...\n")
  svm_model <- caret::train(
    as.formula(paste(target, "~ .")),
    data = model_data,
    metric = "AUROC",
    method = "svmRadial",
    tuneLength = 10,
    trControl = binary_control
  )
  
  # ---- Extract Predictions ----
  glm_preds <- glm_model$pred %>%
    mutate(Model = "GLMnet")
  
  if (univariate || n_predictors == 1) {
    tree_preds <- tree_model$pred %>%
      mutate(Model = "Decision Tree")
  } else {
    tree_preds <- rf_model$pred %>%
      mutate(Model = "Random Forest")
  }
  svm_preds <- svm_model$pred %>%
    mutate(Model = "SVM")
  
  # Combine all predictions
  all_predictions <- bind_rows(glm_preds, tree_preds, svm_preds)
  na_count <- sum(is.na(all_predictions$pred))
  if (na_count > 0) {
    warning(sprintf("Found %d NA predictions in results", na_count))
  }

  # Rename columns
  colnames(all_predictions)[colnames(all_predictions) == "obs"] <- "obs_class"
  colnames(all_predictions)[colnames(all_predictions) == "pred"] <- "pred_class"
  for (level in target_levels) {
    old_col <- level
    new_col <- paste0("prob_", level)
    if (old_col %in% colnames(all_predictions)) {
      colnames(all_predictions)[colnames(all_predictions) == old_col] <- new_col
    }
  }
  # Reorder columns for consistency
  prob_cols <- paste0("prob_", target_levels)
  
  
  all_predictions <- all_predictions %>%
    dplyr::select(Model, rowIndex, obs_class, pred_class, all_of(prob_cols))

  
  # ---- Extract and Compile Results ----
  model_list <- list(glmnet = glm_model, svm = svm_model)
  if (univariate || n_predictors == 1) {
    model_list$tree <- tree_model
    display_names <- c("GLMnet", "Decision Tree", "SVM")
    model_names <- c("glmnet", "tree", "svm")
  } else {
    model_list$rf <- rf_model
    display_names <- c("GLMnet", "Random Forest", "SVM")
    model_names <- c("glmnet", "rf", "svm")
  }
  
  resamp <- caret::resamples(model_list)
  comparison_df <- data.frame(Model = display_names)
  
  # Calculate metrics 
    for (metric in metrics) {
    resamp_metric <- if (metric == "AUROC") "ROC" else metric
    comparison_df[[metric]] <- sapply(model_names, function(model) {
      col_name <- paste0(model, "~", resamp_metric)
      # mean out-of-fold performance across all folds 
      mean(resamp$values[[col_name]], na.rm = TRUE)
    })
  }
  
  # ---- Return Results ----
  return(list(
    models         = model_list,
    predictions    = all_predictions,
    comparison     = comparison_df,
    cv_folds       = cv_folds,
    variables_used = variables,
    target_used    = target,
    target_levels  = target_levels, 
    k_fold         = k_fold,
    metrics        = metrics,
    model_type = if(univariate || n_predictors == 1) "tree" else "rf"
  ))
}

