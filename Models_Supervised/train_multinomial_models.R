 # ---- Multinomial analysis ----
train_multinomial_models <- function(
    data,
    target,
    variables = NULL,
    metrics = c("AUC_Macro", "AUC_Micro", "AUPRC_Macro", "Brier", "StratBrier")) {

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

   valid_metrics <- c("Accuracy", "AUC_Macro", "AUC_Micro",
                     "AUPRC_Macro", "AUPRC_Micro", "Brier", "StratBrier")

  invalid_metrics <- setdiff(metrics, valid_metrics)
  if (length(invalid_metrics) > 0) {
    stop(paste("Invalid metrics:", paste(invalid_metrics, collapse = ", "),
               "\nValid options:", paste(valid_metrics, collapse = ", ")))
  }

  # ---- Prepare Data ----
  model_data <- data[, c(variables, target)]
  n_predictors <- length(variables)
  id_lookup <- if ("id_patient" %in% names(data)) as.character(data[["id_patient"]]) else NULL

  if ("id_patient" %in% names(data)) {
    rownames(model_data) <- data[["id_patient"]]
  }

  # Ensure target is a factor
  if (!is.factor(model_data[[target]])) {
    model_data[[target]] <- factor(model_data[[target]])
  }
  levels(model_data[[target]]) <- make.names(levels(model_data[[target]]))

  class_levels <- levels(model_data[[target]])
  n_classes    <- length(class_levels)
  n_samples    <- nrow(model_data)

  cat(sprintf("Multinomial classification with %d classes: %s\n",
              n_classes, paste(class_levels, collapse = ", ")))


  # ---- Function to calculate additional metrics ----
  calculate_multiclass_metrics <- function(data, lev = NULL, model = NULL) {
    
    missing_cols <- setdiff(lev, names(data))
    if (length(missing_cols) > 0) {
      for (m in missing_cols) data[[m]] <- 0
    }

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
        StratBrier  = NA_real_,
        setNames(rep(NA_real_, length(lev)), paste0("AUC_",   lev)),
        setNames(rep(NA_real_, length(lev)), paste0("AUPRC_", lev))
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

    
    names(auc_per_class)   <- paste0("AUC_",   lev)
    names(auprc_per_class) <- paste0("AUPRC_", lev)
 

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
    
    brier <- mean(rowSums((y_proba - y_true_matrix)^2))
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
      StratBrier  = strat_brier,
      auc_per_class,
      auprc_per_class
    )
  }

  set.seed(42)
  cv_index <- caret::createMultiFolds(model_data[[target]], k = 5, times = 10)
  # ---- Create trainControl ----
  multi_control <- caret::trainControl(
    method          = "repeatedcv",
    number          = 5,
    repeats         = 10,
    index           = cv_index,
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
  # --- GLMnet, or multinomial logistic GLM if only one predictor ---
  # glmnet needs a design matrix with >= 2 columns
  if (n_predictors < 2) {
    cat("Training Multinomial GLM (single predictor)\n")
     set.seed(42); glmnet_model <- tryCatch({
      caret::train(
        as.formula(paste(target, "~ .")),
        data = model_data,
        metric = "AUC_Micro",
        method = "multinom",
        tuneGrid = expand.grid(decay = c(0, 1e-4, 1e-2, 1e-1)),
        trControl = multi_control,
        preProcess = c("center", "scale"),
        trace = FALSE
      )
    }, error = function(e) {
      cat("Multinomial GLM training failed:", e$message, "\n")
      NULL
    })

  } else {
    cat("Training GLMnet\n")
     set.seed(42); glmnet_model <- tryCatch({
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
  }

  # --- Random Forest (with permutation importance enabled) ---
  cat("Training Random Forest\n")
   set.seed(42); rf_model <- tryCatch({
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
      importance = "permutation"
    )
  }, error = function(e) {
    cat("Random Forest training failed:", e$message, "\n")
    NULL
  })

  # --- XGBoost (multi:softprob objective auto-selected by caret for >2 classes) ---
  cat("Training XGBoost\n")
   set.seed(42); xgb_model <- tryCatch({
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
   set.seed(42); pls_model <- tryCatch({
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
  if (!is.null(glmnet_model)) {
    all_predictions$glmnet <- filter_to_best(glmnet_model$pred, glmnet_model$bestTune) %>%
      dplyr::mutate(Model = "GLMnet", model_key = "glmnet")
  }
  if (!is.null(rf_model)) {
    all_predictions$rf <- filter_to_best(rf_model$pred, rf_model$bestTune) %>%
      dplyr::mutate(Model = "Random Forest", model_key = "rf")
  }
  if (!is.null(xgb_model)) {
    all_predictions$xgb <- filter_to_best(xgb_model$pred, xgb_model$bestTune) %>%
      dplyr::mutate(Model = "XGBoost", model_key = "xgb")
  }
  if (!is.null(pls_model)) {
    all_predictions$pls <- filter_to_best(pls_model$pred, pls_model$bestTune) %>%
      dplyr::mutate(Model = "PLS-DA", model_key = "pls")
  }


  combined_preds <- dplyr::bind_rows(all_predictions) %>%
    dplyr::mutate(Rep = sub("^Fold[0-9]+\\.", "", Resample)) %>%
    dplyr::select(model_key, Model, Resample, Rep, rowIndex, obs, pred,
                  dplyr::all_of(class_levels)) %>%
    dplyr::rename(obs_class = obs, pred_class = pred)


  # ---- Compute OOF Metrics ----
    # ---- Compute OOF Metrics ----
  compute_metrics_wrapper <- function(df) {
    data_for_calc <- df %>% dplyr::rename(obs = obs_class, pred = pred_class)
    result <- calculate_multiclass_metrics(data = data_for_calc,
                                           lev = class_levels, model = NULL)
    tibble::as_tibble(t(result))
  }

  # --- per-repeat metrics, then mean/sd across the 10 repeats ---
  oof_by_rep <- combined_preds %>%
    dplyr::group_by(model_key, Model, Rep) %>%
    dplyr::group_modify(~ compute_metrics_wrapper(.x)) %>%
    dplyr::ungroup()

  oof_metrics <- oof_by_rep %>%
    dplyr::group_by(model_key, Model) %>%
    dplyr::summarise(dplyr::across(
      dplyr::where(is.numeric),
      list(mean = ~mean(.x, na.rm = TRUE), sd = ~stats::sd(.x, na.rm = TRUE)),
      .names = "{.col}_{.fn}"), .groups = "drop")

  # --- pooled: average each patient's class probabilities over the 10 repeats ---
  pooled_preds <- combined_preds %>%
    dplyr::group_by(model_key, Model, rowIndex, obs_class) %>%
    dplyr::summarise(dplyr::across(dplyr::all_of(class_levels),
                                   ~mean(.x, na.rm = TRUE)), .groups = "drop")

  prob_mat <- as.matrix(pooled_preds[, class_levels, drop = FALSE])
  pooled_preds$pred_class <- factor(
    class_levels[max.col(prob_mat, ties.method = "first")], levels = class_levels)

  oof_pooled <- pooled_preds %>%
    dplyr::group_by(model_key, Model) %>%
    dplyr::group_modify(~ compute_metrics_wrapper(.x)) %>%
    dplyr::ungroup() %>%
    dplyr::rename_with(~paste0(.x, "_pooled"), -c(model_key, Model))

  oof_metrics <- dplyr::left_join(oof_metrics, oof_pooled,
                                  by = c("model_key", "Model"))

  coverage <- combined_preds %>%
    dplyr::group_by(model_key, Model) %>%
    dplyr::summarise(n_predictions = dplyr::n(),
                     coverage = dplyr::n() / (nrow(model_data) * 10),
                     .groups = "drop")

  oof_metrics <- dplyr::left_join(oof_metrics, coverage, by = c("model_key", "Model"))

  # ---- Compile and Return ----
   trained_models <- Filter(Negate(is.null),
    list(glmnet = glmnet_model, rf = rf_model, xgb = xgb_model, pls = pls_model))

  return(list(
    models         = trained_models,
    predictions    = combined_preds,
    comparison     = oof_metrics,
    variables_used = variables,
    target_used    = target,
    metrics        = metrics,
    method         = "repeatedCV",
    oof_by_rep     = oof_by_rep,
    pooled_predictions = pooled_preds
  ))
}