

# ---- Preprocess data ----
process_luminex_data <- function(raw_data, patient_mapping, pre_threshold = -1) {
  
  # Define metadata columns
  META_COLS <- c('id_patient', 'date_sample', 'PCR', 'days_since_infection', 'Target',
                 "HAI_DENV1" , "HAI_DENV2", "HAI_DENV3","HAI_DENV4" )
  
 
  # Label timepoints based on threshold
  processing_df <- raw_data %>%
    mutate(timepoint = ifelse(days_since_infection <= pre_threshold, 'pre', 'post'))
  
  # Identify antigen columns (exclude metadata)
  antigen_cols <- setdiff(names(processing_df), c(META_COLS, 'timepoint'))
  
  # Dataset 1: Pre/Post Ratio
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
  
  
  
  # Dataset 2: Last Infection (max days since infection)
  last_infection <- raw_data %>%
    group_by(id_patient) %>%
    filter(days_since_infection == max(days_since_infection, na.rm = TRUE)) %>%
    ungroup()
  last_infection  <- last_infection  |> as.data.frame()
  rownames(last_infection) <- last_infection$id_patient
  # Keep patient ID, Target, and antigen columns
  last_infection$Target <- patient_mapping$mapped_target[match(rownames(last_infection), patient_mapping$id_patient)]
  
  
  
  # Dataset 3: First Infection (min days > 0)
  first_infection <- raw_data %>%
    filter(days_since_infection > 0) %>%
    group_by(id_patient) %>%
    filter(days_since_infection == min(days_since_infection, na.rm = TRUE)) %>%
    ungroup()
  first_infection <- first_infection |> as.data.frame()
  rownames(first_infection) <- first_infection$id_patient
  first_infection$Target <- patient_mapping$mapped_target[match(rownames(first_infection), patient_mapping$id_patient)]
  

  # Remove metadata columns not required for further analysis
  last_infection$days_since_infection  <- NULL
  last_infection$id_patient <- NULL
  last_infection$HAI_DENV1 <- NULL
  last_infection$HAI_DENV2 <- NULL
  last_infection$HAI_DENV3 <- NULL
  last_infection$HAI_DENV4 <- NULL
  
  first_infection$days_since_infection <- NULL
  first_infection$id_patient <- NULL
  
  
  # Return all three datasets as a list
  return(list(
    ratio = ratio_df,
    last_infection = last_infection,
    first_infection = first_infection
  ))
}


# ---- Select targets / classification questions ----
select_targets <- function(preprocessed_data,
                           targets = c("flavi", "dengue", 
                                       "zika", "dengue_zika", 
                                       "dengue_serotype", "dengue_serotype_neg",
                                       "dengue_chik"),
                           drop_original_target = TRUE, 
                           negative_label = "true_negative") {
  
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
    dengue_serotype     = c(DENV1 = 1, DENV2 = 2, DENV4 = 3),
    dengue_serotype_neg = c(DENV1 = 1, DENV2 = 2, DENV4 = 3),
    dengue_chik         = c(DENV1 = 0, DENV2 = 0, DENV3 = 0, DENV4 = 0, CHIKV = 1))
  
  # factor levels/labels
  label_specs <- list(
    flavi               = list(levels = c(0, 1),        labels = c("negative", "positive")),
    dengue              = list(levels = c(0, 1),        labels = c("negative", "positive")),
    zika                = list(levels = c(0, 1),        labels = c("negative", "positive")),
    dengue_zika         = list(levels = c(0, 1, 2),     labels = c("negative", "dengue", "zika")),
    dengue_serotype     = list(levels = c(1, 2, 3),     labels = c("DENV1", "DENV2", "DENV4")),
    dengue_serotype_neg = list(levels = c(0, 1, 2, 3),  labels = c("negative", "DENV1", "DENV2", "DENV4")),
    dengue_chik         = list(levels = c(0, 1),        labels = c("dengue", "CHIKV"))
  )
  
  
  # Create a list to store results
  data_with_target <- list()
  
  # Process each target
  for (target in targets) {
    df_copy <- preprocessed_data
    
    # build the chosen target column
    mp <- mapping_list[[target]]
    tgt_chr <- as.character(df_copy$Target)
    
    # Map values: negatives get 0, mapped values get their code, rest get NA
    df_copy[[target]] <- ifelse(
      tgt_chr == negative_label, 
      0,
      as.numeric(dplyr::recode(tgt_chr, !!!mp, .default = NA_real_))
    )
    
    # Remove negatives for classifications that only include infected samples (eg: given infection, classify dengue vs chik)
    if (target %in% c("dengue_chik", "dengue_serotype")) {
      df_copy <- df_copy[df_copy$Target != "true_negative", ]
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
    metrics = c("ROC", "AUPRC", "Brier", "StratBrier"),
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
  
  
  valid_metrics <- c("ROC", "AUPRC", "Brier", "StratBrier")
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
  
  # Convert test indices to train indices for caret
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
    
    # ROC, Sens, Spec from twoClassSummary
    base <- caret::twoClassSummary(data, lev = lev, model = model)
    
    # lev[2] is the positive class, lev[1] is the negative class
    pos_name <- lev[2]
    neg_name <- lev[1]
    
    p_pos <- data[[pos_name]]
    y_true <- as.integer(data$obs == pos_name)
    
    # AUPRC (threshold-free)
    auprc <- tryCatch(
      MLmetrics::PRAUC(p_pos, y_true),
      error = function(e) NA_real_
    )
    
    # Standard Brier score
    brier <- mean((p_pos - y_true)^2)
    

    #all metrics
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
    cat("Training GLM")
    glm_model <- caret::train(
      as.formula(paste(target, "~ .")),
      data = model_data,
      metric = "ROC",
      method = "glm",
      trControl = binary_control
    )
  } else {
    cat("Training GLM with elastic net...\n")
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
  
  # mtry values tuned to number of predictors 
  # mtry cant be > predictors 
  n_features <- length(variables)
  mtry_values <- unique(pmax(1, floor(c(
    sqrt(n_features),           # Default for classification
    n_features / 3,
    n_features / 2,
    n_features
  ))))
  mtry_values <- mtry_values[mtry_values <= n_features]
  
  # ---- Train Random Forest or Decision Tree ----
  if (univariate || n_predictors == 1) {
    cat("Training Decision Tree...\n")
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
    cat("Training Random Forest...\n")
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
  
  cat("Training SVM...\n")
  svm_model <- caret::train(
    as.formula(paste(target, "~ .")),
    data = model_data,
    metric = "ROC",
    method = "svmRadial",
    tuneLength = 10,
    trControl = binary_control
  )
  
  # ---- Extract Predictions ----
  glm_preds <- glm_model$pred %>%
    mutate(Model = "GLMnet")
  
  # Extract RF or Tree predictions depending on what was trained
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
  # Check for NA predictions and warn if found
  na_count <- sum(is.na(all_predictions$pred))
  if (na_count > 0) {
    warning(sprintf("Found %d NA predictions in results", na_count))
  }

  # Rename columns
  colnames(all_predictions)[colnames(all_predictions) == "obs"] <- "obs_class"
  colnames(all_predictions)[colnames(all_predictions) == "pred"] <- "pred_class"
  # Rename probability columns dynamically
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
  # Build the model list conditionally
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
  
  
  # Calculate metrics dynamically
  for (metric in metrics) {
    # Calculate mean for each model
    comparison_df[[metric]] <- sapply(model_names, function(model) {
      col_name <- paste0(model, "~", metric)
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
    n_repeats = 5) {
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
  
  model_data   <- data[, c(variables, target), drop = FALSE]
  X            <- model_data[, variables, drop = FALSE]
  y            <- model_data[[target]]
  y_factor     <- factor(y)
  n_samples    <- nrow(X)
  class_levels <- levels(y_factor)
  n_classes    <- length(class_levels)
  
  cat(sprintf("Multinomial classification with %d classes: %s\n", 
              n_classes, paste(class_levels, collapse = ", ")))
  
  # ---- Storage for results across repeats ----
  rf_all_true <- c(); rf_all_pred <- c(); rf_all_proba <- list(); rf_all_rowindex <- list()
  nb_all_true <- c(); nb_all_pred <- c(); nb_all_proba <- list(); nb_all_rowindex <- list()
  
  # ---- LOOCV with Repeats ----
    for (r in 1:n_repeats) {
      set.seed(r)
      cat(sprintf("Repeat %d/%d...\n", r, n_repeats))
      
      rf_pred  <- rep(NA_character_, n_samples)
      rf_proba <- matrix(NA_real_, nrow = n_samples, ncol = n_classes)
      colnames(rf_proba) <- class_levels
      
      nb_pred  <- rep(NA_character_, n_samples)
      nb_proba <- matrix(NA_real_, nrow = n_samples, ncol = n_classes)
      colnames(nb_proba) <- class_levels
      
      for (i in 1:n_samples) {
        X_train <- X[-i, , drop = FALSE]
        X_test  <- X[i, , drop = FALSE]
        y_train <- y_factor[-i]
        y_test  <- y_factor[i]
       
        # ---------- Random Forest ----------
        tryCatch({
          rf_model <- randomForest::randomForest(
            x = X_train, 
            y = y_train, 
            maxnodes = 4, 
            ntree = 100)
          
          p_rf <- predict(rf_model, X_test, type = "prob")
          p_rf <- as.matrix(p_rf)
          
          # Pad any missing classes with zero-prob columns
          missing_cols <- setdiff(class_levels, colnames(p_rf))
          if (length(missing_cols) > 0) {
            p_rf <- cbind(p_rf, matrix(0, nrow = nrow(p_rf), ncol = length(missing_cols),
                                       dimnames = list(NULL, missing_cols)))
          }
          
          p_rf <- p_rf[, class_levels, drop = FALSE]
          rf_proba[i, ] <- as.numeric(p_rf[1, ])
          rf_pred[i]    <- class_levels[which.max(p_rf[1, ])]
        }, error = function(e) {
          rf_pred[i]    <- NA_character_
          rf_proba[i, ] <- rep(NA_real_, length(class_levels))
        })
        
        # ---------- Naive Bayes ----------
        tryCatch({
          nb_model <- e1071::naiveBayes(x = X_train, y = y_train)
          p_nb <- predict(nb_model, X_test, type = "raw")
          p_nb <- as.matrix(p_nb)
          
          # Pad any missing classes with zero-prob columns
          missing_cols <- setdiff(class_levels, colnames(p_nb))
          if (length(missing_cols) > 0) {
            p_nb <- cbind(p_nb, matrix(0, nrow = nrow(p_nb), ncol = length(missing_cols),
                                       dimnames = list(NULL, missing_cols)))
          }
          p_nb <- p_nb[, class_levels, drop = FALSE]
          nb_proba[i, ] <- as.numeric(p_nb[1, ])
          nb_pred[i]    <- class_levels[which.max(p_nb[1, ])]
        }, error = function(e) {
          nb_pred[i]    <- NA_character_
          nb_proba[i, ] <- rep(NA_real_, length(class_levels))
        })
      }
      
      # ---- Store results for this repeat ----
      rf_all_true[[r]]     <- as.character(y_factor)
      rf_all_pred[[r]]     <- rf_pred
      rf_all_proba[[r]]    <- rf_proba
      rf_all_rowindex[[r]] <- 1:n_samples
      
      nb_all_true[[r]]     <- as.character(y_factor)
      nb_all_pred[[r]]     <- nb_pred
      nb_all_proba[[r]]    <- nb_proba
      nb_all_rowindex[[r]] <- 1:n_samples
    }
    
  # ---- Function to calculate additional metrics ----
  calculate_multiclass_metrics <- function(y_true, y_pred, y_proba, class_levels) {
    # Remove NA predictions
    idx     <- !is.na(y_pred)
    y_true  <- y_true[idx]
    y_pred  <- y_pred[idx]
    y_proba <- y_proba[idx, , drop = FALSE]
    
    if (length(y_true) == 0) {
      return(list(
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
    y_true <- factor(y_true, levels = class_levels)
    y_pred <- factor(y_pred, levels = class_levels)
    
    # Accuracy
    acc <- mean(y_pred == y_true, na.rm = TRUE)
    
    # One-vs-Rest metrics for each class
    auc_per_class   <- numeric(length(class_levels))
    auprc_per_class <- numeric(length(class_levels))
    
    # For micro-averaging AUC/AUPRC
    all_y_binary <- c()
    all_p_class  <- c()
    
    for (k in seq_along(class_levels)) {
      class_k  <- class_levels[k]
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
    y_true_matrix <- matrix(0, nrow = length(y_true), ncol = length(class_levels))
    colnames(y_true_matrix) <- class_levels
    for (k in seq_along(class_levels)) {
      y_true_matrix[, k] <- as.integer(y_true == class_levels[k])
    }
    # Standard multi-class Brier:
    brier <- mean((y_proba - y_true_matrix)^2)
    
    # ---- Stratified Brier (class-weighted) ----
    brier_per_class <- numeric(length(class_levels))
    names(brier_per_class) <- class_levels
    
    for (k in seq_along(class_levels)) {
      idx_k <- y_true == class_levels[k]
      if (sum(idx_k) > 0) {
        # Brier score for samples in class k
        brier_per_class[k] <- mean(rowSums(
          (y_proba[idx_k, , drop = FALSE] - y_true_matrix[idx_k, , drop = FALSE])^2
        ))
      } else {
        brier_per_class[k] <- NA_real_
      }
    }
    
    # Class proportions
    class_counts <- table(y_true)
    class_props  <- numeric(length(class_levels))
    names(class_props) <- class_levels
    
    for (cl in class_levels) {
      if (cl %in% names(class_counts)) {
        class_props[cl] <- class_counts[cl] / length(y_true)
      } else {
        class_props[cl] <- 0
      }
    }
    # Weighted average of per-class Brier scores
    strat_brier <- sum(class_props * brier_per_class, na.rm = TRUE)
    
    list(
      Accuracy       = acc,
      AUC_Macro      = auc_macro,
      AUC_Micro      = auc_micro,
      AUPRC_Macro    = auprc_macro,
      AUPRC_Micro    = auprc_micro,
      Brier          = brier,
      StratBrier     = strat_brier
    )
  }
  
  # ---- Calculate metrics for each repeat ----
  rf_metrics_list <- list()
  nb_metrics_list <- list()
  
  for (r in 1:n_repeats) {
    rf_metrics_list[[r]] <- calculate_multiclass_metrics(
      rf_all_true[[r]],
      rf_all_pred[[r]],
      rf_all_proba[[r]],
      class_levels
    )
    
    nb_metrics_list[[r]] <- calculate_multiclass_metrics(
      nb_all_true[[r]],
      nb_all_pred[[r]],
      nb_all_proba[[r]],
      class_levels
    )
  }
  
  # ---- Aggregate metrics across repeats ----
  aggregate_metrics <- function(metrics_list) {
    metrics_df <- do.call(rbind, lapply(metrics_list, as.data.frame))
    means <- colMeans(metrics_df, na.rm = TRUE)
    sds   <- apply(metrics_df, 2, sd, na.rm = TRUE)
    out <- data.frame(Mean = means, SD = sds, check.names = FALSE)
    out <- out[order(rownames(out)), , drop = FALSE]
    out
  }
  
  rf_summary <- aggregate_metrics(rf_metrics_list)
  nb_summary <- aggregate_metrics(nb_metrics_list)
  
  # ---- Compile comparison table ----
  comparison_df <- data.frame(
    Model           = c("Random Forest", "Naive Bayes"),
    AUC_Macro       = c(rf_summary["AUC_Macro", "Mean"],  nb_summary["AUC_Macro", "Mean"]),
    AUC_Macro_SD    = c(rf_summary["AUC_Macro", "SD"],    nb_summary["AUC_Macro", "SD"]),
    AUC_Micro       = c(rf_summary["AUC_Micro", "Mean"],  nb_summary["AUC_Micro", "Mean"]),
    AUC_Micro_SD    = c(rf_summary["AUC_Micro", "SD"],    nb_summary["AUC_Micro", "SD"]),
    AUPRC_Macro     = c(rf_summary["AUPRC_Macro", "Mean"],nb_summary["AUPRC_Macro", "Mean"]),
    AUPRC_Macro_SD  = c(rf_summary["AUPRC_Macro", "SD"],  nb_summary["AUPRC_Macro", "SD"]),
    AUPRC_Micro     = c(rf_summary["AUPRC_Micro", "Mean"],nb_summary["AUPRC_Micro", "Mean"]),
    AUPRC_Micro_SD  = c(rf_summary["AUPRC_Micro", "SD"],  nb_summary["AUPRC_Micro", "SD"]),
    Brier           = c(rf_summary["Brier", "Mean"],       nb_summary["Brier", "Mean"]),
    Brier_SD        = c(rf_summary["Brier", "SD"],         nb_summary["Brier", "SD"]),
    StratBrier      = c(rf_summary["StratBrier", "Mean"],  nb_summary["StratBrier", "Mean"]),
    StratBrier_SD   = c(rf_summary["StratBrier", "SD"],    nb_summary["StratBrier", "SD"]),
    check.names = FALSE
  )
  
  # ---- Compile predictions with probabilities ----
  # Flatten all repeats
  rf_pred_flat <- unlist(rf_all_pred)
  rf_true_flat <- unlist(rf_all_true)
  rf_proba_flat <- do.call(rbind, rf_all_proba)
  rf_rowidx_flat <- unlist(rf_all_rowindex)
  
  nb_pred_flat <- unlist(nb_all_pred)
  nb_true_flat <- unlist(nb_all_true)
  nb_proba_flat <- do.call(rbind, nb_all_proba)
  nb_rowidx_flat <- unlist(nb_all_rowindex)
  
  
  predictions <- data.frame(
    Model      = c(rep("Random Forest", length(rf_pred_flat)), 
                   rep("Naive Bayes", length(nb_pred_flat))),
    rowIndex   = c(rf_rowidx_flat, nb_rowidx_flat),
    obs_class  = c(rf_true_flat, nb_true_flat),
    pred_class = c(rf_pred_flat, nb_pred_flat),
    stringsAsFactors = FALSE
  )
  
  all_proba_combined <- rbind(rf_proba_flat, nb_proba_flat)
  colnames(all_proba_combined) <- paste0("prob_", class_levels)
  predictions <- cbind(predictions, all_proba_combined)
  
  # ---- Return Results ----
  return(list(
    predictions   = predictions,
    comparison    = comparison_df,
    variables_used = variables,
    target_used    = target,
    n_repeats      = n_repeats,
    n_classes      = n_classes,
    class_levels   = class_levels,
    model_type     = "multinomial"
  ))
}

  
# ---- Train with multiple targets simultaneously ----
train_multiple_targets <- function(
    data_list,  
    variables = NULL,
    k_fold = 5,
    metrics = c("ROC", "AUPRC", "Brier", "StratBrier")) {
  
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
        n_repeats = n_repeats
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
  combined_predictions <- do.call(rbind, all_predictions)
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
    metrics = c("AUC", "AUPRC", "Brier", "StratBrier"),
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
        cat("  -> Result structure:", names(result), "\n")  # Debug: see what's returned
       
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

