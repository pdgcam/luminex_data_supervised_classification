# ---- Function for Supervised Classification ----

# ---- Preprocess data ----
process_luminex_data <- function(raw_data, patient_mapping, pre_threshold = -1) {
  
  # Define metadata columns
  META_COLS <- c('id_patient', 'date_sample', 'PCR', 'days_since_infection', 'Target',
                 "HAI_DENV1" , "HAI_DENV2", "HAI_DENV3","HAI_DENV4" )
  
  # Label timepoints based on threshold
  processing_df <- raw_data %>%
    mutate(timepoint = ifelse(days_since_infection <= pre_threshold, 'pre', 'post'))
  
  # Identify antigen columns, exclude metadata
  antigen_cols <- setdiff(names(processing_df), c(META_COLS, 'timepoint'))
  
  # Dataset 1: Post/pre  Ratio
  grouped_means <- processing_df %>%
    group_by(id_patient, timepoint) %>%
    summarise(across(all_of(antigen_cols), ~mean(.x, na.rm = TRUE)), .groups = 'drop')
  
  averages_pivot <- grouped_means %>%
    pivot_wider(
      names_from = timepoint,
      values_from = all_of(antigen_cols),
      names_glue = "{.value}_{timepoint}"
    )
  
  post_cols <- grep("_post$", names(averages_pivot), value = TRUE)
  pre_cols <- grep("_pre$", names(averages_pivot), value = TRUE)
  base_names <- sub("_post$", "", post_cols)
  
  post <- averages_pivot %>%
    dplyr::select(all_of(post_cols)) %>%
    as.data.frame()
  names(post) <- base_names
  rownames(post) <- averages_pivot$id_patient
  
  pre <- averages_pivot %>%
    dplyr::select(all_of(pre_cols)) %>%
    as.data.frame()
  names(pre) <- base_names
  rownames(pre) <- averages_pivot$id_patient
  ratio_df <- (post / pre) %>%
    na.omit()

  # add target using the patient mapping
  ratio_df$Target <- patient_mapping$mapped_target[match(rownames(ratio_df), patient_mapping$id_patient)]
  
  # Dataset 2: Last draw (ie draw at max days since infection)
  cross_sectional_data <- raw_data %>%
    group_by(id_patient) %>%
    filter(days_since_infection == max(days_since_infection, na.rm = TRUE)) %>%
    ungroup()
  cross_sectional_data  <- cross_sectional_data  |> as.data.frame()
  rownames(cross_sectional_data) <- cross_sectional_data$id_patient
  # Keep patient ID, Target, and antigen columns
  cross_sectional_data$Target <- patient_mapping$mapped_target[match(rownames(cross_sectional_data), patient_mapping$id_patient)]


  # Remove metadata columns not required for further analysis
  cross_sectional_data$days_since_infection  <- NULL
  cross_sectional_data$id_patient <- NULL
  cross_sectional_data$HAI_DENV1 <- NULL
  cross_sectional_data$HAI_DENV2 <- NULL
  cross_sectional_data$HAI_DENV3 <- NULL
  cross_sectional_data$HAI_DENV4 <- NULL
  
 
  # Return all three datasets as a list
  return(list(
    ratio = ratio_df,
    cross_sectional_data = cross_sectional_data
  ))
}



# ---- Select targets / classification questions ----
select_targets <- function(preprocessed_data,
                           targets = c("flavi", "dengue", 
                                       "zika", "dengue_zika", 
                                       "dengue_serotype", "dengue_serotype_neg",
                                       "dengue_chik"),
                           drop_original_target = TRUE, 
                           negative_label = "no_infection") {
  
  if (!"Target" %in% names(preprocessed_data)) {
    stop('Column "Target" not found in preprocessed_data.')
  }
  
  # Validate all targets
  valid_targets <- c("flavi", "dengue", 
                     "zika", "dengue_zika", 
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
    dengue_zika         = c(DENV1 = 1, DENV2 = 1, DENV3 = 1, DENV4 = 1, ZIKV = 2),
    dengue_serotype     = c(DENV1 = 1, DENV2 = 2, DENV3 = 3, DENV4 = 4),
    dengue_serotype_neg = c(DENV1 = 1, DENV2 = 2, DENV3 = 3, DENV4 = 4),
    dengue_chik         = c(DENV1 = 0, DENV2 = 0, DENV3 = 0, DENV4 = 0, CHIKV = 1))
  
  # factor levels/labels
  label_specs <- list(
    flavi               = list(levels = c(0, 1),           labels = c("negative", "positive")),
    dengue              = list(levels = c(0, 1),           labels = c("negative", "positive")),
    zika                = list(levels = c(0, 1),           labels = c("negative", "positive")),
    dengue_zika         = list(levels = c(0, 1, 2),        labels = c("negative", "dengue", "zika")),
    dengue_serotype     = list(levels = c(1, 2, 3, 4),     labels = c("DENV1", "DENV2", "DENV3", "DENV4")),
    dengue_serotype_neg = list(levels = c(0, 1, 2, 3, 4),  labels = c("negative", "DENV1", "DENV2", "DENV3", "DENV4")),
    dengue_chik         = list(levels = c(0, 1),           labels = c("dengue", "CHIKV"))
  )
  
  
  # List to store results
  data_with_target <- list()
  
  # Process each target
  for (target in targets) {
    df_copy <- preprocessed_data
    
    # Build the chosen target column
    mp <- mapping_list[[target]]
    tgt_chr <- as.character(df_copy$Target)
    
    # Map values: negative label get 0, mapped values get their code, rest get NA
    df_copy[[target]] <- ifelse(
      tgt_chr == negative_label, 
      0,
      as.numeric(dplyr::recode(tgt_chr, !!!mp, .default = NA_real_))
    )
    
    # Remove negatives for classifications that only include infected samples (eg: given infection, classify dengue vs chik)
    if (target %in% c("dengue_chik", "dengue_serotype")) {
      df_copy <- df_copy[df_copy$Target != "no_infection", ]
    }

    # convert to factor with labels
    spec <- label_specs[[target]]
    df_copy[[target]] <-
      factor(df_copy[[target]], 
             levels = spec$levels, 
             labels = spec$labels)
    
    # drop original Target column if requested
    if (drop_original_target) {
      if ("Target" %in% names(df_copy)) {
        df_copy[["Target"]] <- NULL
      }
    }
    
    data_with_target[[target]] <- df_copy
  }
  return(data_with_target)
}



# ---- Binary analysis ----
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
  
  

# ---- Train with multiple targets simultaneously ----
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
      comparison_with_target$Target <- target_name
      all_comparisons[[target_name]] <- comparison_with_target
      
      # Add target name to predictions dataframe
      predictions_with_target <- result$predictions
      predictions_with_target$Target <- target_name
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
      comparison_with_target$Target <- target_name
      all_comparisons[[target_name]] <- comparison_with_target
      
      # Add target name to predictions dataframe
      predictions_with_target <- result$predictions
      predictions_with_target$Target <- target_name
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
  col_order <- c("Target", setdiff(names(combined_comparison), "Target"))
  combined_comparison <- combined_comparison[, col_order]
  
  # Combine all predictions dataframes
  combined_predictions <- dplyr::bind_rows(all_predictions)
  rownames(combined_predictions) <- NULL
  # Reorder columns to put Target first
  pred_col_order <- c("Target", setdiff(names(combined_predictions), "Target"))
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

                            

# ---- Univariate nalysis ----
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
      all_comparisons[[target_name]] <- do.call(rbind, target_comparisons)
    }
    if (length(target_predictions) > 0) {
      all_predictions[[target_name]] <- do.call(rbind, target_predictions)
    }
  }
  
  # Combine all comparison dataframes across targets
  combined_comparison <- NULL
  if (length(all_comparisons) > 0) {
    combined_comparison <- do.call(rbind, all_comparisons)
    rownames(combined_comparison) <- NULL
    # Reorder columns to put Target and Variable first
    col_order <- c("Target", "Variable", setdiff(names(combined_comparison), c("Target", "Variable")))
    combined_comparison <- combined_comparison[, col_order]
  }
  
  # Combine all predictions dataframes across targets
  combined_predictions <- NULL
  if (length(all_predictions) > 0) {
    combined_predictions <- do.call(rbind, all_predictions)
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

