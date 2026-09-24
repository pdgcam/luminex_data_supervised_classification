# --- Source functions
source(here('Models/train_binary_models.R'))
source(here('Models/train_multinomial_models.R'))
source(here('Functions.R'))

# =============================================================================
# Forward selection wrapper for train_binary_models() / train_multinomial_models()
#
# Reuses your existing model-training functions unchanged. Greedily builds up
# a variable set: start empty, at each step try every remaining candidate
# antigen added to the current set, keep whichever addition performs best,
# repeat. Produces a step-by-step performance table you can plot.
#
# Requires: dplyr, ggplot2 (only for plot_forward_selection)
# =============================================================================


# ---- Metric direction lookup: is bigger better, or smaller better? ----------
# Used to decide which candidate "wins" at each step.
.metric_direction <- c(
  ROC          = "max",
  AUPRC        = "max",
  Brier        = "min",
  StratBrier   = "min",
  AUC_Macro    = "max",
  AUC_Micro    = "max",
  AUPRC_Macro  = "max",
  AUPRC_Micro  = "max",
  Accuracy     = "max",
  logLoss      = "min"
)


# ---- Helper: extract the best value of `metric` from a comparison df -------
# comparison_df has one row per model (Model, ROC, AUPRC, Brier, ...).
# If `model` is NULL, takes the best-performing model that step (max/min per
# metric direction). If `model` is given (e.g. "rf"), pins selection to that
# model's row instead, so selection is comparable model-to-model across steps.
.get_best_metric_value <- function(comparison_df, metric, model = NULL) {

  if (!metric %in% names(comparison_df)) {
    stop(sprintf(
      "Metric '%s' not found in comparison data (columns available: %s)",
      metric, paste(names(comparison_df), collapse = ", ")
    ))
  }

  direction <- .metric_direction[[metric]]
  if (is.null(direction)) direction <- "max"  # sensible fallback for unknown metrics

  df <- comparison_df
  if (!is.null(model)) {
    df <- df[df$Model == model, , drop = FALSE]
    if (nrow(df) == 0) {
      stop(sprintf(
        "Model '%s' not found among trained models (available: %s)",
        model, paste(unique(comparison_df$Model), collapse = ", ")
      ))
    }
  }

  vals <- df[[metric]]
  if (all(is.na(vals))) return(list(value = NA_real_, model = NA_character_))

  best_idx <- if (direction == "max") which.max(vals) else which.min(vals)
  list(value = vals[best_idx], model = df$Model[best_idx])
}





forward_selection <- function(
    data,
    target,
    variables,
    positive_class = NULL,
    metrics = c("ROC", "AUPRC", "Brier"),
    selection_metric = "auto",
    selection_model = NULL,
    max_vars = NULL,
    min_improvement = 0) {
 
  # ---- Input validation (mirrors train_binary_models) ----
  if (!target %in% names(data)) {
    stop(paste("Target column", target, "not found in data"))
  }
 
  missing_vars <- setdiff(variables, names(data))
  if (length(missing_vars) > 0) {
    stop(paste("Variables not found in data:", paste(missing_vars, collapse = ", ")))
  }
 
  n_classes <- length(unique(data[[target]]))
 
  # Resolve which metric ranks candidates
  target_selection_metric <- selection_metric
  if (identical(selection_metric, "auto")) {
    target_selection_metric <- if (n_classes == 2) "ROC" else "AUC_Micro"
  }
  cat(sprintf("Ranking candidates by: %s\n", target_selection_metric))
 
  direction   <- .metric_direction[[target_selection_metric]]
  if (is.null(direction)) direction <- "max"
  best_so_far <- if (direction == "max") -Inf else Inf
 
  remaining <- variables
  selected  <- character(0)
 
  step_rows  <- list()   # full comparison (all models) for the WINNER, per step
  candidate_rows <- list()   # full comparison for EVERY candidate tried, per step
 
  n_steps <- length(remaining)
  if (!is.null(max_vars)) n_steps <- min(n_steps, max_vars)
  fold_index <- NULL
  if (n_classes == 2) {
    set.seed(42)
    fold_index <- caret::createMultiFolds(
      factor(data[[target]]),
      k = 5,
      times = 10
    )
  }

 
  for (step in seq_len(n_steps)) {
 
    cat(sprintf("\n-- Step %d: %d candidate(s) remaining --\n",
                step, length(remaining)))
 
    step_candidate_metrics <- list()
    for (cand in remaining) {
 
      trial_vars <- c(selected, cand)
      cat(sprintf("  Trying + %-20s (set size = %d)\n", cand, length(trial_vars)))
 
      result <- tryCatch({
        if (n_classes == 2) {
          train_binary_models(
            data           = data,
            target         = target,
            variables      = trial_vars,
            positive_class = positive_class,
            metrics        = metrics, 
            #cv_repeats      = 10, 
            fold_index      = fold_index
          )
        } else {
          train_multinomial_models(
            data       = data,
            target     = target,
            variables  = trial_vars,
            univariate = FALSE,   # n_predictors==1 still forces decision-tree-only
            metrics    = metrics
          )
        }
      }, error = function(e) {
        warning(sprintf("Training failed for %s + [%s]: %s",
                         target, paste(trial_vars, collapse = ", "), e$message))
        NULL
      })
 
      if (is.null(result)) next
 
      comp <- result$comparison
      comp$candidate_variable <- cand
      comp$step     <- step
      comp$set_size <- length(trial_vars)
      step_candidate_metrics[[cand]] <- comp
    }
 
    if (length(step_candidate_metrics) == 0) {
      warning(sprintf("No candidates could be trained at step %d. Stopping.", step))
      break
    }
 
    candidate_metrics_df   <- dplyr::bind_rows(step_candidate_metrics)
    candidate_rows[[step]] <- candidate_metrics_df
 
    # Score each candidate (best model's value, or the pinned model's value)
    candidate_scores <- sapply(names(step_candidate_metrics), function(cand) {
      .get_best_metric_value(step_candidate_metrics[[cand]],
                              target_selection_metric,
                              model = selection_model)$value
    })
 
    if (all(is.na(candidate_scores))) {
      warning(sprintf("All candidate scores were NA at step %d. Stopping.", step))
      break
    }
 
    best_cand  <- if (direction == "max") {
      names(candidate_scores)[which.max(candidate_scores)]
    } else {
      names(candidate_scores)[which.min(candidate_scores)]
    }
    best_score <- candidate_scores[[best_cand]]
 
    cat(sprintf("  --> Selected '%s'  (%s = %.4f)\n",
                best_cand, target_selection_metric, best_score))
 
    improved <- if (direction == "max") (best_score - best_so_far) else (best_so_far - best_score)
    if (step > 1 && improved < min_improvement) {
      cat(sprintf(
        "  Improvement (%.4f) below min_improvement (%.4f) -- stopping before adding '%s'.\n",
        improved, min_improvement, best_cand))
      break
    }
    best_so_far <- best_score
 
    chosen_comp <- step_candidate_metrics[[best_cand]]
    chosen_comp$variable_added <- best_cand
    chosen_comp$selected_set   <- paste(c(selected, best_cand), collapse = " + ")
    step_rows[[step]] <- chosen_comp
 
    selected  <- c(selected, best_cand)
    remaining <- setdiff(remaining, best_cand)
  }
 
  combined_steps <- dplyr::bind_rows(step_rows)
  if (nrow(combined_steps) > 0) {
    combined_steps$target <- target
    col_order <- c("target", "step", "set_size", "variable_added", "selected_set",
                   setdiff(names(combined_steps),
                           c("target", "step", "set_size", "variable_added", "selected_set")))
    combined_steps <- combined_steps[, col_order]
    rownames(combined_steps) <- NULL
  }
 
  combined_candidates <- dplyr::bind_rows(candidate_rows)
  if (nrow(combined_candidates) > 0) {
    combined_candidates$target <- target
    rownames(combined_candidates) <- NULL
  }
 
  list(
    combined_steps         = combined_steps,
    combined_candidates    = combined_candidates,
    selected_vars          = selected,
    selection_metric_used  = target_selection_metric
  )
}
 


# ------ OLD CODE -------


# function for univariate analysis across multiple targets
train_multiple_targets_univariate <- function(
    data_list,  
    variables = NULL,
    positive_class_map = list(),
    univariate = FALSE, 
    metrics = c("ROC", "AUPRC", "Brier", "StratBrier")) {
  
  # Store all results
  all_results <-  list()
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
    
    # Initialize storage for this target
    all_results[[target_name]] <- list()
    all_comparisons[[target_name]] <- list()
    all_predictions[[target_name]] <- list()
    
    # Loop through each variable for univariate analysis
    cat("\nTesting", length(variables), "variables individually \n\n")
    print(head(variables))

    pos_class <- positive_class_map[[target_name]]
    
    for (i in seq_along(variables)) {
      var <- variables[i]
      cat(sprintf("[%d/%d] Testing variable: %s\n", i, length(variables), var))
      
      # Check if variable exists in data
      if (!var %in% names(current_data)) {
        warning(paste("Variable", var, "not found in data. Skipping."))
        next
      }
      
        if (n_classes == 2) {
          result <- train_binary_models(
            data = current_data,
            target = target_col,
            variables = var,
            positive_class = pos_class,
            metrics = metrics
          )

          # Store results
          all_results[[target_name]][[var]] <- result
          
          
          # Add target name to comparison dataframe
          comparison_with_target <- result$comparison
          comparison_with_target$target <- target_name
          comparison_with_target$variable <- var
          all_comparisons[[target_name]][[var]] <- comparison_with_target
          
          # Add target name to predictions dataframe
          predictions_with_target <- result$predictions
          predictions_with_target$target <- target_name
          predictions_with_target$variable <- var
          all_predictions[[target_name]][[var]] <- predictions_with_target

        } else if (n_classes > 2) {
          result <- train_multinomial_models(
            data = current_data,
            target = target_col,
            variables = var,
            univariate = TRUE,
            metrics = metrics
          )
          all_results[[target_name]][[var]] <- result
          
          
          # Add target name to comparison dataframe
          comparison_with_target <- result$comparison
          comparison_with_target$target <- target_name
          comparison_with_target$variable <- var
          all_comparisons[[target_name]][[var]] <- comparison_with_target
          
          # Add target name to predictions dataframe
          predictions_with_target <- result$predictions
          predictions_with_target$target <- target_name
          predictions_with_target$variable <- var
          all_predictions[[target_name]][[var]] <- predictions_with_target
        
      }
    } }  
      
      # Combine all comparison dataframes
      combined_comparison  <- dplyr::bind_rows(unlist(all_comparisons,  recursive = FALSE))
      rownames(combined_comparison) <- NULL
      col_order <- c("target", "variable", setdiff(names(combined_comparison), c("target", "variable")))
      combined_comparison <- combined_comparison[, col_order]
      
      # Combine all predictions dataframes
      combined_predictions <- dplyr::bind_rows(unlist(all_predictions, recursive = FALSE))
      rownames(combined_predictions) <- NULL
      pred_col_order <- c("target", "variable", setdiff(names(combined_predictions), c("target", "variable")))
      combined_predictions <- combined_predictions[, pred_col_order]
      
      return(list(
        results_by_target = all_results,
        combined_comparison = combined_comparison,
        combined_predictions = combined_predictions,
        summary = list(
          n_targets = length(all_results),
          target_names = names(all_results),
          n_variables_tested = length(variables),
          variables_tested = variables,
          metrics = metrics
        )
      ))
    
}

