
train_binary_models <- function(
    data,
    target,
    variables = NULL,
    positive_class = NULL,
    metrics = c("ROC", "AUPRC", "Brier")) {

  # ---- Input Validation ----
  if (!target %in% names(data)) {
    stop(paste("Target column", target, "not found in data"))
  }

  if (is.null(variables)) {
    variables <- setdiff(names(data), c(target, "id_patient"))
  } else {
    missing_vars <- setdiff(variables, names(data))
    if (length(missing_vars) > 0) {
      stop(paste("Variables not found in data:", paste(missing_vars, collapse = ", ")))
    }
  }

  # --- get data to input to model ----
  model_data <- data[, c(variables, target)]
  n_predictors <- ncol(model_data) - 1

  if ("id_patient" %in% names(data)) {
    rownames(model_data) <- data[["id_patient"]]
  }


  valid_metrics <- c("ROC", "AUPRC", "Brier", "StratBrier")
  invalid_metrics <- setdiff(metrics, valid_metrics)
  if (length(invalid_metrics) > 0) {
    stop(paste("Invalid metrics:", paste(invalid_metrics, collapse = ", "),
               "\nValid options:", paste(valid_metrics, collapse = ", ")))
  }


  # Ensure target is a factor
  if (!is.factor(model_data[[target]])) {
    model_data[[target]] <- factor(model_data[[target]])
  }

  # If positive_class is not specified, use the first level of the factor as the positive class
  if (is.null(positive_class)) {
    positive_class <- levels(model_data[[target]])[1]
    message(sprintf("No positive_class specified — using '%s' as positive class", positive_class))
  }

  # if positive class not in target levels, throw error
  if (!positive_class %in% levels(model_data[[target]])) {
    stop(sprintf("positive_class '%s' not found in target levels: %s",
                 positive_class, paste(levels(model_data[[target]]), collapse = ", ")))
  }


  # caret expects positive class to be first level
  model_data[[target]] <- relevel(model_data[[target]], ref = positive_class)

  target_levels <- levels(model_data[[target]])
  cat(sprintf("Binary classification: %s (class 1) vs %s (class 2)\n",
              target_levels[1], target_levels[2]))


  # ---- Set Up Train Control ----
  combinedBinary <- function(data, lev = NULL, model = NULL) {
    stopifnot(length(lev) == 2)

    # ROC, Sensitivity, Specificity from twoClassSummary
    # this is for AUPRC and Brier
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
    method = "cv",
    number = 5,
    summaryFunction = combinedBinary,
    classProbs = TRUE,
    verboseIter = FALSE,
    savePredictions = "final")

  # mean AUC across all folds
  average_auc <- function(model, model_data, target) {

    results <- model$results
    best_tune <- model$bestTune
    n_resamples <- nrow(model$resample)

    if (length(best_tune) == 0 || ncol(best_tune) == 0) {
      best_row <- results[1, ]
    } else {
      idx <- which(apply(results[, names(best_tune), drop = FALSE], 1,
                          function(row) all(row == unlist(best_tune))))
      best_row <- results[if (length(idx) == 0) 1 else idx[1], ]
    }

    auc <- best_row$ROC
    se <- best_row$ROCSD / sqrt(n_resamples)

    list(auc = auc, sd = best_row$ROCSD, ci_low = best_row$ROC - 1.96 * se, ci_high = best_row$ROC + 1.96 * se)
  }


  # ---- Train Models ----
  # model set: GLMnet, Random Forest, XGBoost, PLS-DA

  # --- GLMnet ---
  cat("Training GLMnet\n")
  glmnet_model <- caret::train(
    as.formula(paste(target, "~ .")),
    data = model_data,
    metric = "ROC",
    method = "glmnet",
    tuneGrid = expand.grid(
      alpha = seq(0, 1, length.out = 5),
      lambda = 10^seq(-3, 1, length.out = 20)
    ),
    trControl = binary_control,
    preProcess = c("center", "scale")
  )

  # --- Random Forest ---
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
    trControl = binary_control,
    preProcess = c("center", "scale"),
    importance = "permutation" 
  )

  # --- XGBoost ---
  cat("Training XGBoost\n")
  xgb_model <- caret::train(
    as.formula(paste(target, "~ .")),
    data = model_data,
    metric = "ROC",
    method = "xgbTree",
    tuneGrid = expand.grid(
      nrounds = c(50, 100, 150),
      max_depth = c(2, 4, 6),
      eta = c(0.01, 0.1, 0.3),
      gamma = 0,
      colsample_bytree = 0.8,
      min_child_weight = 1,
      subsample = 0.8
    ),
    trControl = binary_control,
    preProcess = c("center", "scale"),
    verbose = 0
  )

  # --- PLS-DA ---
  # ncomp cant exceed number of predictors
  max_ncomp <- min(10, n_predictors)
  cat("Training PLS-DA\n")
  pls_model <- caret::train(
    as.formula(paste(target, "~ .")),
    data = model_data,
    metric = "ROC",
    method = "pls",
    tuneGrid = expand.grid(ncomp = seq_len(max_ncomp)),
    trControl = binary_control,
    preProcess = c("center", "scale")
  )

  # ---- Compute Pooled AUCs ----
  aucs <- list(
    glmnet = average_auc(glmnet_model, model_data, target),
    rf     = average_auc(rf_model, model_data, target),
    xgb    = average_auc(xgb_model, model_data, target),
    pls    = average_auc(pls_model, model_data, target)
  )

  # ---- Build Model List + Names ----
  model_list    <- list(glmnet = glmnet_model, rf = rf_model, xgb = xgb_model, pls = pls_model)
  model_names   <- c("glmnet", "rf", "xgb", "pls")
  display_names <- c("GLMnet", "Random Forest", "XGBoost", "PLS-DA")


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
    mapply(function(name, model) {
      filter_to_best(model$pred, model$bestTune) %>%
        mutate(Model = display_names[match(name, model_names)])
    }, model_names, model_list, SIMPLIFY = FALSE)
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


  # ---- Compile Metrics using Pooled AUC ----
  comparison_df <- data.frame(
    Model    = display_names,
    ROC      = sapply(model_names, function(m) aucs[[m]]$auc),
    ROC_low  = sapply(model_names, function(m) aucs[[m]]$ci_low),
    ROC_high = sapply(model_names, function(m) aucs[[m]]$ci_high)
  )

  # Calculate other metrics from model results
  for (metric in setdiff(metrics, "ROC")) {
    comparison_df[[metric]] <- sapply(model_names, function(model) {
      model_obj <- model_list[[model]]
      if (length(model_obj$bestTune) == 0 || all(is.na(model_obj$bestTune))) {
        # No tuning parameters
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
    aucs           = aucs,
    variables_used = variables,
    target_used    = target,
    target_levels  = target_levels,
    method         = "cv",
    metrics        = metrics))
}
