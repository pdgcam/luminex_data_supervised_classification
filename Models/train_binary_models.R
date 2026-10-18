match_best_row <- function(results, best_tune, model_label = "model") {
  if (is.null(best_tune) || length(best_tune) == 0 || ncol(best_tune) == 0) {
    return(1L)
  }
  keep <- rep(TRUE, nrow(results))
  for (p in names(best_tune)) {
    if (p %in% names(results)) {
      keep <- keep & (results[[p]] == best_tune[[p]])
    }
  }
  idx <- which(keep)
  if (length(idx) == 0) {
    warning(sprintf(
      "Could not match bestTune to results for '%s' - falling back to row 1. Reported metrics are probably wrong.",
      model_label
    ))
    return(1L)
  }
  idx[1]
}


train_binary_models <- function(
    data,
    target,
    variables = NULL,
    positive_class = NULL,
    cv_repeats = 1,
    fold_index = NULL,
    metrics = c("ROC", "AUPRC", "Brier")) {

  # ---- Input Validation ----
  if (!target %in% names(data)) {
    stop(paste("Target column", target, "not found in data"))
  }
  valid_metrics <- c("ROC", "AUPRC", "Brier", "StratBrier")
  invalid_metrics <- setdiff(metrics, valid_metrics)
  if (length(invalid_metrics) > 0) {
    stop(paste("Invalid metrics:", paste(invalid_metrics, collapse = ", "),
               "\nValid options:", paste(valid_metrics, collapse = ", ")))
  }
  if (is.null(variables)) {
    variables <- setdiff(names(data), c(target, "id_patient"))
  } else {
    missing_vars <- setdiff(variables, names(data))
    if (length(missing_vars) > 0) {
      stop(paste("Variables not found in data:", paste(missing_vars, collapse = ", ")))
    }
  }


  # --- Build model data ----
  model_data   <- as.data.frame(data[, c(variables, target)])
  n_predictors <- length(variables)
  id_lookup <- if ("id_patient" %in% names(data)) as.character(data[["id_patient"]]) else NULL
  n_incomplete <- sum(!stats::complete.cases(model_data))
  if (n_incomplete > 0) {
    if (!is.null(fold_index)) {
      stop(sprintf(
        "%d row(s) have missing values in target/variables, but fold_index was supplied. Filter to complete cases before creating folds.",
        n_incomplete))
    }
    warning(sprintf("%d row(s) with missing values will be dropped by caret.", n_incomplete))
  }
  if (!is.null(fold_index) && max(unlist(fold_index)) > nrow(model_data)) {
    stop("fold_index refers to row positions beyond nrow(data) - folds were built on a different frame.")
  }


  # ---- Ensure target is a factor
  if (!is.factor(model_data[[target]])) {
    model_data[[target]] <- factor(model_data[[target]])
  }
  levels(model_data[[target]]) <- make.names(levels(model_data[[target]]))
  if (is.null(positive_class)) {
    positive_class <- levels(model_data[[target]])[1]
    message(sprintf("No positive_class specified - using '%s' as positive class", positive_class))
  } else {
    positive_class <- make.names(positive_class)
    if (!positive_class %in% levels(model_data[[target]])) {
      stop(sprintf("positive_class '%s' not found in target levels: %s",
                   positive_class, paste(levels(model_data[[target]]), collapse = ", ")))
    }
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
    method  = if (cv_repeats > 1) "repeatedcv" else "cv",
    number = 5,
    repeats = if (cv_repeats > 1) cv_repeats else NA,
    index = fold_index,    
    summaryFunction = combinedBinary,
    classProbs = TRUE,
    verboseIter = FALSE,
    savePredictions = "final")

  # mean AUC across all folds 
  # mean across held out folds 
  # 50 held-out AUCs — 5 folds × 10 repeats
  average_auc <- function(model, model_label = "model") {

    results <- model$results
    n_resamples <- nrow(model$resample)  
    best_row    <- results[match_best_row(results, model$bestTune, model_label), ]

    se <- best_row$ROCSD / sqrt(n_resamples)

    list(auc     = best_row$ROC,
         sd      = best_row$ROCSD,
         ci_low  = max(0, best_row$ROC - 1.96 * se),
         ci_high = min(1, best_row$ROC + 1.96 * se))
  }


  # ---- Train Models ----
  # model set: GLMnet, Random Forest, XGBoost, PLS-DA


 # --- GLMnet ---
  if (n_predictors < 2) {
    cat("Training Bayesian Logistic GLM (single predictor)\n")
    glmnet_model <- caret::train(
      as.formula(paste(target, "~ .")),
      data = model_data,
      metric = "ROC",
      method = "bayesglm",
      family = "binomial",
      trControl = binary_control,
      preProcess = c("center", "scale")
    )
  } else {
    cat("Training GLMnet\n")
      set.seed(42); glmnet_model <- caret::train(
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
  }

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
    set.seed(42); rf_model <- caret::train(
    as.formula(paste(target, "~ .")),
    data = model_data,
    metric = "ROC",
    method = "ranger",
    tuneGrid = expand.grid(
      mtry = mtry_values,
      splitrule = c("gini", "extratrees"),
      min.node.size = c(1, 5, 10), stringsAsFactors = FALSE),
    trControl = binary_control,
    preProcess = c("center", "scale"),
    importance = "permutation" 
  )

  # --- XGBoost ---
  cat("Training XGBoost\n")
    set.seed(42);xgb_model <- caret::train(
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
    set.seed(42);cat("Training PLS-DA\n")
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
    glmnet = average_auc(glmnet_model, "glmnet"),
    rf     = average_auc(rf_model,     "rf"),
    xgb    = average_auc(xgb_model,    "xgb"),
    pls    = average_auc(pls_model,    "pls")
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


  all_predictions <- dplyr::bind_rows(
    mapply(function(name, model) {
      filter_to_best(model$pred, model$bestTune) |>
        dplyr::mutate(Model = display_names[match(name, model_names)])
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

  all_predictions$id_patient <- id_lookup[all_predictions$rowIndex]

  # ---- Pooled out-of-fold AUC ----
       
  # One AUC over all out-of-fold predictions, rather than the mean of per-fold
  # AUCs. With repeated CV each patient appears once per repeat, so average
  # within patient first - concatenating would treat correlated repeats as
  # independent observations.

  pos_col <- paste0("prob_", target_levels[1])
  pooled_auc <- function(preds, display_label) {
    d <- preds[preds$Model == display_label, , drop = FALSE]
    if (nrow(d) == 0) return(list(auc = NA_real_, ci_low = NA_real_, ci_high = NA_real_))

    d <- stats::aggregate(d[[pos_col]],
                          by  = list(rowIndex = d$rowIndex, obs_class = d$obs_class),
                          FUN = mean)
    names(d)[3] <- "p"

    tryCatch({
      r  <- pROC::roc(d$obs_class, d$p,
                      levels = rev(target_levels),
                      direction = "<", quiet = TRUE)
      ci <- as.numeric(pROC::ci.auc(r))
      list(auc = ci[2], ci_low = ci[1], ci_high = ci[3])
    }, error = function(e) list(auc = NA_real_, ci_low = NA_real_, ci_high = NA_real_))
  }
  
  pooled <- setNames(lapply(display_names, function(dn) pooled_auc(all_predictions, dn)),
                     model_names)


  # ---- Compile Metrics using Pooled AUC ----
    comparison_df <- data.frame(
    model_key = model_names,
    Model     = display_names,
    ROC       = unname(sapply(model_names, function(m) pooled[[m]]$auc)),
    ROC_low   = unname(sapply(model_names, function(m) pooled[[m]]$ci_low)),
    ROC_high  = unname(sapply(model_names, function(m) pooled[[m]]$ci_high)),
    ROC_mean_fold = unname(sapply(model_names, function(m) aucs[[m]]$auc)), #mean AUC on each held out fold 
    ROC_sd_fold   = unname(sapply(model_names, function(m) aucs[[m]]$sd)),
    stringsAsFactors = FALSE
  )
  rownames(comparison_df) <- NULL

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
    method         = binary_control$method,
    n_resamples    = nrow(rf_model$resample),
    glmnet_method  = if (n_predictors < 2) "bayesglm" else "glmnet",
    metrics        = metrics))
}
