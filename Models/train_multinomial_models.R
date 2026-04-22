  # ---- Multinomial analysis ----
train_multinomial_models <- function(
    data,
    target,
    variables = NULL,
    metrics = c("ROC", "AUPRC", "Brier", "StratBrier")) {
  
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
  
  valid_metrics <- c("ROC", "AUPRC", "Brier", "StratBrier", "logLoss")
  invalid_metrics <- setdiff(metrics, valid_metrics)
  if (length(invalid_metrics) > 0) {
    stop(paste("Invalid metrics:", paste(invalid_metrics, collapse = ", "),
               "\nValid options:", paste(valid_metrics, collapse = ", ")))
  }
  
  
  # ---- Prepare Data ----
  model_data <- data[, c(variables, target)]
  
  # Ensure target is a factor
  if (!is.factor(model_data[[target]])) {
    model_data[[target]] <- as.factor(model_data[[target]])
  }
  
  class_levels <- levels(model_data[[target]])
  n_classes    <- length(class_levels) 
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
        Accuracy    = NA_real_,
        AUC_Macro   = NA_real_,
        AUC_Micro   = NA_real_,
        AUPRC_Macro = NA_real_,
        AUPRC_Micro = NA_real_,
        Brier       = NA_real_,
        StratBrier  = NA_real_
      ))
    }
    
    # Ensure factors
    y_true <- factor(y_true, levels = lev)  
    y_pred <- factor(y_pred, levels = lev)  
    
    # Accuracy
    acc <- mean(y_pred == y_true, na.rm = TRUE)
    
    # One-vs-Rest AUC and AUPRC per class
    auc_per_class   <- numeric(length(lev))  
    auprc_per_class <- numeric(length(lev)) 
    
    # Containers for micro-averaging
    all_y_binary <- c()
    all_p_class  <- c()
    
    for (k in seq_along(lev)) {  
      class_k  <- lev[k]  
      y_binary <- as.integer(y_true == class_k)
      p_class  <- y_proba[, k]
      
      # Accumulate for micro averaging
      all_y_binary <- c(all_y_binary, y_binary)
      all_p_class  <- c(all_p_class,  p_class)
      
      # Skip if only one class present (e.g. during LOOCV folds)
      if (length(unique(y_binary)) < 2) {
        auc_per_class[k]   <- NA_real_
        auprc_per_class[k] <- NA_real_
        next
      }
      
      # AUC (pROC with explicit direction to avoid flipping)
    auc_per_class[k] <- tryCatch({
    roc_obj <- pROC::roc(
      response  = y_binary,
      predictor = p_class,
      levels    = c(0, 1),    # 0 = negative, 1 = positive
      direction = "<",        # predictor increases with positive class
      quiet     = TRUE
    )
    as.numeric(pROC::auc(roc_obj))
  }, error = function(e) NA_real_)
        
      # AUPRC
      auprc_per_class[k] <- tryCatch({
        MLmetrics::PRAUC(p_class, y_binary)
      }, error = function(e) NA_real_)
    }
    
    # Macro AUC - average of per-class AUCs (
    auc_macro   <- mean(auc_per_class,   na.rm = TRUE)
    auprc_macro <- mean(auprc_per_class, na.rm = TRUE)
    
    # Micro AUC - computed on pooled binary labels and probabilities (across all samples and classes)
    auc_micro <- tryCatch({
      if (length(unique(all_y_binary)) >= 2) {
        as.numeric(pROC::auc(
          response  = all_y_binary,
          predictor = all_p_class,
          direction = "<",
          quiet     = TRUE
        ))
      } else {
        NA_real_
      }
    }, error = function(e) NA_real_)
    
    # Micro AUPRC
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
    brier <- mean((y_proba - y_true_matrix)^2)
    
    # ---- Stratified Brier (equal weight per class) ----
    brier_per_class <- numeric(length(lev))  
    names(brier_per_class) <- lev 
    for (k in seq_along(lev)) {  
      idx_k <- y_true == lev[k]  
      if (sum(idx_k) > 0) {
        brier_per_class[k] <- mean(rowSums(
          (y_proba[idx_k, , drop = FALSE] - y_true_matrix[idx_k, , drop = FALSE])^2
        ))
      } else {
        brier_per_class[k] <- NA_real_
      }
    }
    strat_brier <- mean(brier_per_class, na.rm = TRUE)
    
    c(
      Accuracy    = acc,
      AUC_Macro   = auc_macro,
      AUC_Micro   = auc_micro,
      AUPRC_Macro = auprc_macro,
      AUPRC_Micro = auprc_micro,
      Brier       = brier,
      StratBrier  = strat_brier
    )
  }
  
  # ---- Create trainControl ----
  multi_control <- caret::trainControl(
      method          = "LOOCV",
      summaryFunction = calculate_multiclass_metrics,
      classProbs      = TRUE,
      verboseIter     = FALSE,
      savePredictions = "all",
      allowParallel   = FALSE
    )

  
  # ---- Train Models ----
  cat("Training Random Forest (multiclass)\n")
  rf_model <- tryCatch({
    caret::train(
      as.formula(paste(target, "~ .")),
      data       = model_data,
      metric     = "AUC_Micro",
      method     = "ranger",
      tuneLength = 3,
      num.trees  = 500,
      trControl  = multi_control
    )
  }, error = function(e) {
    cat("Random Forest training failed:", e$message, "\n")
    NULL
  })
  
  cat("Training Naive Bayes\n")
  nb_model <- tryCatch({
    caret::train(
      as.formula(paste(target, "~ .")),
      data      = model_data,
      method    = "nb",
      metric    = "AUC_Micro",
      trControl = multi_control
    )
  }, error = function(e) {
    cat("Naive Bayes training failed:", e$message, "\n")
    NULL
  })
  
  # ---- Extract Predictions ----
  filter_to_best <- function(preds, best_tune) {
    mask <- rep(TRUE, nrow(preds))
    for (param in names(best_tune)) {
      if (param %in% names(preds)) {
        mask <- mask & (preds[[param]] == best_tune[[param]])
      }
    }
    preds[mask, ]
  }

  all_predictions <- list()

  if (!is.null(rf_model)) {
    all_predictions$rf <- filter_to_best(rf_model$pred, rf_model$bestTune) %>%
      dplyr::mutate(Model = "Random Forest")
  }

  if (!is.null(nb_model)) {
    all_predictions$nb <- filter_to_best(nb_model$pred, nb_model$bestTune) %>%
      dplyr::mutate(Model = "NaiveBayes")
  }

  # Combine predictions
  combined_preds <- dplyr::bind_rows(all_predictions) %>%
    dplyr::select(Model, rowIndex, obs, pred, dplyr::all_of(class_levels)) %>%
    dplyr::rename(
      obs_class  = obs,
      pred_class = pred
    )
  
  # ---- Compute OOF Metrics ----
  compute_metrics_wrapper <- function(df) {
    data_for_calc <- df %>%
      dplyr::rename(obs = obs_class, pred = pred_class)
    
    result <- calculate_multiclass_metrics(
      data  = data_for_calc,
      lev   = class_levels,  
      model = NULL
    )
    tibble::as_tibble(t(result))
  }
  
  oof_metrics <- combined_preds %>%
    dplyr::group_by(Model) %>%
    dplyr::group_modify(~ compute_metrics_wrapper(.x)) %>%
    dplyr::ungroup()
  
  # ---- Compile and Return ----
  trained_models <- Filter(Negate(is.null), list(
    rf = rf_model,
    nb = nb_model
  ))
  
  return(list(
    models         = trained_models,
    predictions    = combined_preds,
    oof_metrics    = oof_metrics,
    comparison     = oof_metrics,
    variables_used = variables,
    target_used    = target,
    metrics        = metrics
  ))
}