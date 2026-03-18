  
# ---- Multinomial analysis ----
train_multinomial_models <- function(
    data,
    target,
    variables = NULL,
    k_fold = "LOOCV", 
    metrics = c("AUROC", "AUPRC", "Brier", "StratBrier")) {
  
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
  
  valid_metrics <- c("AUROC", "AUPRC", "Brier", "StratBrier", "logLoss")
  invalid_metrics <- setdiff(metrics, valid_metrics)
  if (length(invalid_metrics) > 0) {
    stop(paste("Invalid metrics:", paste(invalid_metrics, collapse = ", "),
               "\nValid options:", paste(valid_metrics, collapse = ", ")))
  }
  
  if (!k_fold %in% c(5, 10, "LOOCV")) {
    warning("k_fold should be 5, 10 or LOOCV. Using provided value.")
  }
  
  # ---- Prepare Data ----
  model_data <- data[, c(variables, target)]
  
  # Ensure target is a factor
  if (!is.factor(model_data[[target]])) {
    model_data[[target]] <- as.factor(model_data[[target]])
  }
  
  class_levels <- levels(model_data[[target]])
  n_classes <- length(class_levels) 
  n_samples    <- nrow(model_data)
  
  cat(sprintf("Multinomial classification with %d classes: %s\n", 
              n_classes, paste(class_levels, collapse = ", ")))
  
  
  # ---- Function to calculate additional metrics ----
  calculate_multiclass_metrics <- function(data, lev = NULL, model = NULL) {
    
    # Extract components
    y_true <- data$obs
    y_pred <- data$pred
    
    # Extract probability matrix (columns named after class levels)
    y_proba <- as.matrix(data[, lev, drop = FALSE])
    
    # Remove NA predictions
    idx     <- !is.na(y_pred)
    y_true  <- y_true[idx]
    y_pred  <- y_pred[idx]
    y_proba <- y_proba[idx, , drop = FALSE]
    
    if (length(y_true) == 0) {
      return(c( 
        Accuracy = NA_real_,
        AUC_Macro = NA_real_,
        AUC_Micro = NA_real_,
        AUPRC_Macro = NA_real_,
        AUPRC_Micro = NA_real_,
        Brier = NA_real_,
        StratBrier = NA_real_
      ))
    }
    
    # Ensure factors for confusion matrix
    y_true <- factor(y_true, levels = lev)  
    y_pred <- factor(y_pred, levels = lev)  
    
    # Accuracy
    acc <- mean(y_pred == y_true, na.rm = TRUE)
    
    # One-vs-Rest metrics for each class
    auc_per_class   <- numeric(length(lev))  
    auprc_per_class <- numeric(length(lev)) 
    
    # micro-averaging AUC/AUPRC
    all_y_binary <- c()
    all_p_class  <- c()
    
    for (k in seq_along(lev)) {  
      class_k  <- lev[k]  
      y_binary <- as.integer(y_true == class_k)
      p_class  <- y_proba[, k]
      
      # Store for micro calculation
      all_y_binary <- c(all_y_binary, y_binary)
      all_p_class  <- c(all_p_class, p_class)
      
      if (length(unique(y_binary)) < 2) {
        auc_per_class[k]   <- NA_real_
        auprc_per_class[k] <- NA_real_
        next
      }
      
      # AUC
      auc_per_class[k] <- tryCatch({
        MLmetrics::AUC(p_class, y_binary)
      }, error = function(e) NA_real_)
      
      # AUPRC
      auprc_per_class[k] <- tryCatch({
        MLmetrics::PRAUC(p_class, y_binary)
      }, error = function(e) NA_real_)
    }
    
    auc_macro  <- mean(auc_per_class,   na.rm = TRUE)
    auprc_macro <- mean(auprc_per_class, na.rm = TRUE)
    
    # Micro averages (all OvR comparisons pooled)
    auc_micro <- tryCatch({
      if(length(unique(all_y_binary)) >= 2) {
        MLmetrics::AUC(all_p_class, all_y_binary) 
      } else {
        NA_real_
      }
    }, error = function(e) NA_real_)
    
    auprc_micro <- tryCatch({
      if (length(unique(all_y_binary)) >= 2) {
        MLmetrics::PRAUC(all_p_class, all_y_binary) 
      } else {
        NA_real_
      }
    }, error = function(e) NA_real_)
    
    # ---- Multi-class Brier Score ----
    y_true_matrix <- matrix(0, nrow = length(y_true), ncol = length(lev))  
    colnames(y_true_matrix) <- lev  
    for (k in seq_along(lev)) {  
      y_true_matrix[, k] <- as.integer(y_true == lev[k])  
    }
    # multi-class Brier:
    brier <- mean((y_proba - y_true_matrix)^2)
    
    # ---- Stratified Brier (class-weighted) ----
    brier_per_class <- numeric(length(lev))  
    names(brier_per_class) <- lev 
    for (k in seq_along(lev)) {  
      idx_k <- y_true == lev[k]  
      if (sum(idx_k) > 0) {
        # Brier score for samples in class k
        brier_per_class[k] <- mean(rowSums(
          (y_proba[idx_k, , drop = FALSE] - y_true_matrix[idx_k, , drop = FALSE])^2
        ))
      } else {
        brier_per_class[k] <- NA_real_
      }
    }
    # Equal weight to each class (balanced)
    strat_brier <- mean(brier_per_class, na.rm = TRUE)
    
    c(
      Accuracy       = acc,
      AUC_Macro      = auc_macro,
      AUC_Micro      = auc_micro,
      AUPRC_Macro    = auprc_macro,
      AUPRC_Micro    = auprc_micro,
      Brier          = brier,
      StratBrier     = strat_brier
    )
  }
  
  # ---- Create trainControl ----
  if (k_fold == "LOOCV") {
    multi_control <- caret::trainControl(
      method = "LOOCV",
      summaryFunction = calculate_multiclass_metrics,
      classProbs = TRUE,
      verboseIter = FALSE,
      savePredictions = "all",
      allowParallel = FALSE
    )
    cv_folds <- NULL
  } 
  else {
    cv_folds <- caret::createFolds(model_data[[target]], k = k_fold)
    train_indices <- lapply(cv_folds, function(test_idx) {
      setdiff(seq_len(nrow(model_data)), test_idx)
    })
    multi_control <- caret::trainControl(
      summaryFunction = calculate_multiclass_metrics,
      method = "cv",
      classProbs = TRUE,
      verboseIter = FALSE,
      savePredictions = "all",
      index = train_indices
    )
  }
  
  # ---- Train Models ----
  # Random Forest 
  cat("Training Random Forest (multiclass)\n")
  rf_model <- tryCatch({
    caret::train(
      as.formula(paste(target, "~ .")),
      data = model_data,
      metric = "AUC_Micro",
      method = "ranger",
      tuneLength = 3,
      num.trees = 500,
      trControl = multi_control
    )
  }, error = function(e) {
    cat("Random Forest training failed:", e$message, "\n")
    NULL
  })
  
  # Naive Bayes
  cat("Training Naive Bayes\n")
  nb_model <- tryCatch({
    caret::train(
      as.formula(paste(target, "~ .")),
      data = model_data,
      method = "nb",
      metric = "AUC_Micro",
      trControl = multi_control
    )
  }, error = function(e) {
    cat("NB training failed:", e$message, "\n")
    NULL
  })
  
  # ---- Extract Predictions ----
  all_predictions <- list()
  
  if (!is.null(rf_model)) {
    rf_preds <- rf_model$pred %>%
      dplyr::mutate(Model = "Random Forest")
    all_predictions$rf <- rf_preds
  }
  
  if (!is.null(nb_model)) {
    nb_preds <- nb_model$pred %>%
      dplyr::mutate(Model = "NaiveBayes")
    all_predictions$nb <- nb_preds
  }
  
  # Ensure each preds df has all class prob columns
  ensure_prob_cols <- function(df, classes) {
    missing <- setdiff(classes, names(df))
    for (cl in missing) df[[cl]] <- NA_real_
    df[, c(setdiff(names(df), classes), classes)]
  }
  
  all_predictions <- lapply(all_predictions, ensure_prob_cols, classes = class_levels)
  
  # Combine all predictions
  combined_preds <- dplyr::bind_rows(all_predictions) %>%
    dplyr::select(Model, rowIndex, obs, pred, dplyr::all_of(class_levels)) %>%
    dplyr::rename(
      obs_class  = obs,
      pred_class = pred
    )
  
  # Prepare data in format expected by calculate_multiclass_metrics
  compute_metrics_wrapper <- function(df) {
    data_for_calc <- df %>%
      dplyr::rename(obs = obs_class, pred = pred_class)
    
    metrics <- calculate_multiclass_metrics(
      data = data_for_calc,
      lev = class_levels,  
      model = NULL
    )
    tibble::as_tibble(t(metrics))
  }
  
  # Compute OOF metrics
  oof_metrics <- combined_preds %>%
    dplyr::group_by(Model) %>%
    dplyr::group_modify(~ compute_metrics_wrapper(.x)) %>%
    dplyr::ungroup()
  
  
  # ---- Extract and Compile Results ----
  trained_models <- Filter(Negate(is.null), list(
    rf = rf_model,
    nb = nb_model
  ))

  comparison_df <- oof_metrics
  
  # ---- Return Results ----
  return(list(
    models = trained_models,
    predictions = combined_preds,
    oof_metrics = oof_metrics,
    comparison = comparison_df,
    cv_folds = cv_folds,
    variables_used = variables,
    target_used = target,
    k_fold = k_fold,
    metrics = metrics
  ))
}
  
