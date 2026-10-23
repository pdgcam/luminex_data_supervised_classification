

# ---- Forward selection with Random Forest only ----
# Same logic as the original version; the only changes are:
#   - each candidate is scored with train_binary_models(which_models = "rf")
#   - folds are patient-grouped (make_patient_folds) instead of row-based
#   - binary targets only
forward_selection_rf <- function(
    data,
    target,
    variables,
    positive_class   = NULL,
    selection_metric = "ROC",    # "ROC", "ROC_mean_fold", "AUPRC" or "Brier" ("auto" = "ROC")
    max_vars         = NULL,
    min_improvement  = -Inf,     # -Inf = never stop early, run through all variables
    tune_grid        = NULL,     # NULL = the full RF grid in train_binary_models
    importance       = "none",
    cv_repeats       = 10,
    seed             = 42) {
 
  if (!target %in% names(data)) stop(paste("Target column", target, "not found in data"))
  missing_vars <- setdiff(variables, names(data))
  if (length(missing_vars) > 0) {
    stop(paste("Variables not found in data:", paste(missing_vars, collapse = ", ")))
  }
  variables <- unique(variables)
 
  keep <- stats::complete.cases(data[, c(target, variables), drop = FALSE])
  if (sum(!keep) > 0) {
    cat(sprintf("Dropping %d row(s) with missing values.\n", sum(!keep)))
    data <- data[keep, , drop = FALSE]
  }
 
  if (identical(selection_metric, "auto")) selection_metric <- "ROC"
  valid <- c("ROC", "ROC_mean_fold", "AUPRC", "Brier")
  if (!selection_metric %in% valid) {
    stop("selection_metric must be one of: ", paste(valid, collapse = ", "))
  }
 
  cat(sprintf("Ranking candidates by: %s (Random Forest)\n", selection_metric))
  direction   <- if (selection_metric == "Brier") "min" else "max"   # lower Brier is better
  best_so_far <- if (direction == "max") -Inf else Inf
 
  remaining <- variables
  selected  <- character(0)
  step_rows <- candidate_rows <- list()
 
  n_steps <- length(remaining)
  if (!is.null(max_vars)) n_steps <- min(n_steps, max_vars)
 
  # Patient-grouped folds, built once and shared by every candidate at every step
  groups     <- if ("id_patient" %in% names(data)) data$id_patient else seq_len(nrow(data))
  fold_index <- make_patient_folds(data[[target]], groups, repeats = cv_repeats, seed = seed)
 
  for (step in seq_len(n_steps)) {
 
    cat(sprintf("\n-- Step %d: %d candidate(s) remaining --\n", step, length(remaining)))
 
    step_candidates <- list()
    for (cand in remaining) {
 
      trial_vars <- c(selected, cand)
      cat(sprintf("  Trying + %-20s (set size = %d)\n", cand, length(trial_vars)))
 
      result <- tryCatch(
        suppressMessages(train_binary_models(
          data = data, target = target, variables = trial_vars,
          positive_class = positive_class, fold_index = fold_index,
          which_models = "rf", rf_grid = tune_grid)),
        error = function(e) {
          warning(sprintf("Training failed for %s + [%s]: %s",
                          target, paste(trial_vars, collapse = ", "), e$message))
          NULL
        })
      if (is.null(result)) next
 
      comp <- result$comparison
      comp$candidate_variable <- cand
      comp$step               <- step
      comp$set_size           <- length(trial_vars)
      step_candidates[[cand]] <- comp
    }
 
    if (length(step_candidates) == 0) {
      warning(sprintf("No candidates could be trained at step %d. Stopping.", step))
      break
    }
    candidate_rows[[step]] <- dplyr::bind_rows(step_candidates)
 
    scores <- sapply(step_candidates, function(x) x[[selection_metric]])
    if (all(is.na(scores))) {
      warning(sprintf("All candidate scores were NA at step %d. Stopping.", step))
      break
    }
 
    best_cand  <- names(scores)[if (direction == "max") which.max(scores) else which.min(scores)]
    best_score <- scores[[best_cand]]
    cat(sprintf("  --> Selected '%s'  (%s = %.4f)\n", best_cand, selection_metric, best_score))
 
    improved <- if (direction == "max") best_score - best_so_far else best_so_far - best_score
    if (step > 1 && improved < min_improvement) {
      cat(sprintf("  Improvement (%.4f) below min_improvement (%.4f) -- stopping before adding '%s'.\n",
                  improved, min_improvement, best_cand))
      break
    }
    best_so_far <- best_score
 
    chosen <- step_candidates[[best_cand]]
    chosen$variable_added <- best_cand
    chosen$selected_set   <- paste(c(selected, best_cand), collapse = " + ")
    step_rows[[step]]     <- chosen
 
    selected  <- c(selected, best_cand)
    remaining <- setdiff(remaining, best_cand)
  }
 
  combined_steps <- dplyr::bind_rows(step_rows)
  if (nrow(combined_steps) > 0) {
    combined_steps$target <- target
    front <- c("target", "step", "set_size", "variable_added", "selected_set")
    combined_steps <- combined_steps[, c(front, setdiff(names(combined_steps), front))]
    rownames(combined_steps) <- NULL
  }
 
  combined_candidates <- dplyr::bind_rows(candidate_rows)
  if (nrow(combined_candidates) > 0) {
    combined_candidates$target <- target
    rownames(combined_candidates) <- NULL
  }
 
  list(combined_steps        = combined_steps,
       combined_candidates   = combined_candidates,
       selected_vars         = selected,
       selection_metric_used = selection_metric,
       fold_index            = fold_index)
}
 