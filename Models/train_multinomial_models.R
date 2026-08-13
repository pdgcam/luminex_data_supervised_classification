  # ---- Multinomial analysis ----
train_multinomial_models <- function(
    data,
    target,
    variables = NULL,
    metrics = c("ROC", "AUPRC", "Brier", "StratBrier"),
    univariate = FALSE) {

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

  valid_metrics <- c("ROC", "AUPRC", "Brier", "StratBrier", "logLoss")
  invalid_metrics <- setdiff(metrics, valid_metrics)
  if (length(invalid_metrics) > 0) {
    stop(paste("Invalid metrics:", paste(invalid_metrics, collapse = ", "),
               "\nValid options:", paste(valid_metrics, collapse = ", ")))
  }

  # ---- Prepare Data ----
  model_data <- data[, c(variables, target)]
  n_predictors <- ncol(model_data) - 1

  if ("id_patient" %in% names(data)) {
    rownames(model_data) <- data[["id_patient"]]
  }

  # Ensure target is a factor
  if (!is.factor(model_data[[target]])) {
    model_data[[target]] <- as.factor(model_data[[target]])
  }

  class_levels <- levels(model_data[[target]])
  n_classes    <- length(class_levels)
  n_samples    <- nrow(model_data)

  cat(sprintf("Multinomial classification with %d classes: %s\n",
              n_classes, paste(class_levels, collapse = ", ")))

  # single-predictor case now checked consistently up front (fixes a mismatch
  # in the original code, where training branched on `univariate` alone but
  # prediction extraction branched on `univariate || n_predictors == 1`)
  is_univariate <- univariate || n_predictors == 1

  # ---- Function to calculate additional metrics ----
  calculate_multiclass_metrics <- function(data, lev = NULL, model = NULL) {

    y_true <- data$obs
    y_pred <- data$pred
    y_proba <- as.matrix(data[, lev, drop = FALSE])

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

    y_true <- factor(y_true, levels = lev)
    y_pred <- factor(y_pred, levels = lev)

    acc <- mean(y_pred == y_true, na.rm = TRUE)

    auc_per_class   <- numeric(length(lev))
    auprc_per_class <- numeric(length(lev))

    all_y_binary <- c()
    all_p_class  <- c()

    for (k in seq_along(lev)) {
      class_k  <- lev[k]
      y_binary <- as.integer(y_true == class_k)
      p_class  <- y_proba[, k]

      all_y_binary <- c(all_y_binary, y_binary)
      all_p_class  <- c(all_p_class,  p_class)

      if (length(unique(y_binary)) < 2) {
        auc_per_class[k]   <- NA_real_
        auprc_per_class[k] <- NA_real_
        next
      }

      auc_per_class[k] <- tryCatch({
        roc_obj <- pROC::roc(
          response  = y_binary,
          predictor = p_class,
          levels    = c(0, 1),
          direction = "<",
          quiet     = TRUE
        )
        as.numeric(pROC::auc(roc_obj))
      }, error = function(e) NA_real_)

      auprc_per_class[k] <- tryCatch({
        MLmetrics::PRAUC(p_class, y_binary)
      }, error = function(e) NA_real_)
    }

    auc_macro   <- mean(auc_per_class,   na.rm = TRUE)
    auprc_macro <- mean(auprc_per_class, na.rm = TRUE)

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

    auprc_micro <- tryCatch({
      if (length(unique(all_y_binary)) >= 2) {
        MLmetrics::PRAUC(all_p_class, all_y_binary)
      } else {
        NA_real_
      }
    }, error = function(e) NA_real_)

    y_true_matrix <- matrix(0, nrow = length(y_true), ncol = length(lev))
    colnames(y_true_matrix) <- lev
    for (k in seq_along(lev)) {
      y_true_matrix[, k] <- as.integer(y_true == lev[k])
    }
    brier <- mean((y_proba - y_true_matrix)^2)

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

  n_features <- length(variables)
  mtry_values <- unique(pmax(1, floor(c(
    sqrt(n_features),
    n_features / 3,
    n_features / 2,
    n_features
  ))))
  mtry_values <- mtry_values[mtry_values <= n_features]

  # ---- Train Models ----
  glmnet_model <- rf_model <- xgb_model <- pls_model <- tree_model <- NULL

  if (is_univariate) {

    cat("Training Decision Tree\n")
    tree_model <- tryCatch({
      caret::train(
        as.formula(paste(target, "~ .")),
        data = model_data,
        metric = "AUC_Micro",
        method = "rpart",
        tuneGrid = expand.grid(cp = c(0, 10^seq(-4, -1, length.out = 20))),
        trControl = multi_control,
        control   = rpart::rpart.control(minsplit = 2, minbucket = 1),
        preProcess = c("center", "scale")
      )
    }, error = function(e) {
      cat("Decision Tree training failed:", e$message, "\n")
      NULL
    })

  } else {

    # --- GLMnet (multinomial family auto-selected by caret for >2 classes) ---
    cat("Training GLMnet\n")
    glmnet_model <- tryCatch({
      caret::train(
        as.formula(paste(target, "~ .")),
        data = model_data,
        metric = "AUC_Micro",
        method = "glmnet",
        tuneGrid = expand.grid(
          alpha = seq(0, 1, length.out = 5),
          lambda = 10^seq(-3, 1, length.out = 20)
        ),
        trControl = multi_control,
        preProcess = c("center", "scale")
      )
    }, error = function(e) {
      cat("GLMnet training failed:", e$message, "\n")
      NULL
    })

    # --- Random Forest (with permutation importance enabled) ---
    cat("Training Random Forest\n")
    rf_model <- tryCatch({
      caret::train(
        as.formula(paste(target, "~ .")),
        data = model_data,
        metric = "AUC_Micro",
        method = "ranger",
        tuneGrid = expand.grid(
          mtry = mtry_values,
          splitrule = c("gini", "extratrees"),
          min.node.size = c(1, 5, 10)),
        trControl = multi_control,
        preProcess = c("center", "scale"),
        importance = "permutation"    # <-- enables varImp() later, same fix as binary
      )
    }, error = function(e) {
      cat("Random Forest training failed:", e$message, "\n")
      NULL
    })

    # --- XGBoost (multi:softprob objective auto-selected by caret for >2 classes) ---
    cat("Training XGBoost\n")
    xgb_model <- tryCatch({
      caret::train(
        as.formula(paste(target, "~ .")),
        data = model_data,
        metric = "AUC_Micro",
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
        trControl = multi_control,
        preProcess = c("center", "scale"),
        verbose = 0
      )
    }, error = function(e) {
      cat("XGBoost training failed:", e$message, "\n")
      NULL
    })

    # --- PLS-DA ---
    max_ncomp <- min(10, n_predictors)
    cat("Training PLS-DA\n")
    pls_model <- tryCatch({
      caret::train(
        as.formula(paste(target, "~ .")),
        data = model_data,
        metric = "AUC_Micro",
        method = "pls",
        tuneGrid = expand.grid(ncomp = seq_len(max_ncomp)),
        trControl = multi_control,
        preProcess = c("center", "scale")
      )
    }, error = function(e) {
      cat("PLS-DA training failed:", e$message, "\n")
      NULL
    })
  }

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

  if (is_univariate) {
    if (!is.null(tree_model)) {
      all_predictions$tree <- filter_to_best(tree_model$pred, tree_model$bestTune) %>%
        dplyr::mutate(Model = "Decision Tree")
    }
  } else {
    if (!is.null(glmnet_model)) {
      all_predictions$glmnet <- filter_to_best(glmnet_model$pred, glmnet_model$bestTune) %>%
        dplyr::mutate(Model = "GLMnet")
    }
    if (!is.null(rf_model)) {
      all_predictions$rf <- filter_to_best(rf_model$pred, rf_model$bestTune) %>%
        dplyr::mutate(Model = "Random Forest")
    }
    if (!is.null(xgb_model)) {
      all_predictions$xgb <- filter_to_best(xgb_model$pred, xgb_model$bestTune) %>%
        dplyr::mutate(Model = "XGBoost")
    }
    if (!is.null(pls_model)) {
      all_predictions$pls <- filter_to_best(pls_model$pred, pls_model$bestTune) %>%
        dplyr::mutate(Model = "PLS-DA")
    }
  }

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

  coverage <- combined_preds %>%
    dplyr::group_by(Model) %>%
    dplyr::summarise(
      n_predictions = dplyr::n(),
      coverage = dplyr::n() / nrow(model_data),
      .groups = "drop"
    )

  oof_metrics <- oof_metrics %>%
    dplyr::left_join(coverage, by = "Model")

  # ---- Compile and Return ----
  trained_models <- if (is_univariate) {
    Filter(Negate(is.null), list(tree = tree_model))
  } else {
    Filter(Negate(is.null), list(glmnet = glmnet_model, rf = rf_model, xgb = xgb_model, pls = pls_model))
  }

  return(list(
    models         = trained_models,
    predictions    = combined_preds,
    oof_metrics    = oof_metrics,
    comparison     = oof_metrics,
    variables_used = variables,
    target_used    = target,
    metrics        = metrics,
    method         = "LOOCV"
  ))
}