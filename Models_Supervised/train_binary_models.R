# ---- Patient-grouped, stratified, repeated CV folds ----
# Patients (not rows) are assigned to folds, so no patient is ever in train and test
# at once. Patients are stratified by which classes their rows contain (only negative,
# only positive, or both), so this also works when a patient's outcome changes over time.
# Pass groups = row numbers for ordinary stratified CV.
make_patient_folds <- function(y, groups, k = 5, repeats = 1, seed = 42) {
  pats  <- unique(groups)
  strat <- sapply(split(as.character(y), factor(groups, levels = pats)),
                  function(v) paste(sort(unique(v)), collapse = "+"))
  set.seed(seed)
  unlist(lapply(seq_len(repeats), function(r) {
    pat_fold <- caret::createFolds(factor(strat), k = k, list = FALSE)
    row_fold <- pat_fold[match(groups, pats)]
    setNames(lapply(seq_len(k), function(f) which(row_fold != f)),
             sprintf("Fold%d.Rep%02d", seq_len(k), r))
  }), recursive = FALSE)
}


train_binary_models <- function(data, target, variables = NULL,
                                positive_class = NULL,
                                cv_repeats = 1, fold_index = NULL,
                                which_models = c("glmnet", "rf", "xgb", "pls"),
                                rf_grid = NULL) {

  which_models <- match.arg(which_models, several.ok = TRUE)

  # ---- Inputs ----
  if (is.null(variables)) variables <- setdiff(names(data), c(target, "id_patient"))
  missing_cols <- setdiff(c(target, variables), names(data))
  if (length(missing_cols) > 0) stop("Not found in data: ", paste(missing_cols, collapse = ", "))

  model_data <- as.data.frame(data[, c(variables, target)])
  if (anyNA(model_data)) stop("Missing values in target/variables - filter to complete cases first.")
  id_patient <- if ("id_patient" %in% names(data)) as.character(data$id_patient) else NULL
  if (anyNA(id_patient)) stop("id_patient has missing values - every row needs a patient ID for grouped CV.")

  # ---- Target: 2-level factor, valid names, positive class first (caret convention) ----
  y <- factor(model_data[[target]])
  if (nlevels(y) != 2) stop("Target must have exactly 2 classes, found: ", nlevels(y))
  levels(y) <- make.names(levels(y))

  positive_class <- if (is.null(positive_class)) levels(y)[1] else make.names(positive_class)
  if (!positive_class %in% levels(y)) {
    stop("positive_class '", positive_class, "' not in: ", paste(levels(y), collapse = ", "))
  }
  model_data[[target]] <- relevel(y, ref = positive_class)
  target_levels <- levels(model_data[[target]])
  message("Positive class: ", target_levels[1], " | Negative: ", target_levels[2])

  # ---- Folds: shared by every model ----
  groups <- if (is.null(id_patient)) seq_len(nrow(model_data)) else id_patient
  if (is.null(fold_index)) {
    fold_index <- make_patient_folds(model_data[[target]], groups, repeats = cv_repeats)
  }
  if (max(unlist(fold_index)) > nrow(model_data)) stop("fold_index doesn't match the rows of data.")

  # Every held-out fold needs both classes, or AUC/AUPRC can't be computed
  one_class <- sapply(fold_index, function(tr) length(unique(model_data[[target]][-tr])) < 2)
  if (any(one_class)) {
    yy    <- model_data[[target]]
    n_pat <- sapply(levels(yy), function(l) length(unique(groups[yy == l])))
    n_row <- table(yy)
    stop(sum(one_class), " of ", length(fold_index), " test folds contain only one class. ",
         "Patients with >= 1 row of each class: ",
         paste(names(n_pat), n_pat, sep = " = ", collapse = ", "),
         " (rows: ", paste(names(n_row), n_row, sep = " = ", collapse = ", "), "). ",
         "5-fold CV needs at least 5 patients in each class.")
  }

  # ---- Per-fold metrics: ROC (= AUC) / Sens / Spec + AUPRC + Brier ----
  binary_summary <- function(data, lev, model = NULL) {
    p <- data[[lev[1]]]                 # prob of positive class
    y <- as.integer(data$obs == lev[1]) # 1 = positive
    c(caret::twoClassSummary(data, lev),
      AUPRC = MLmetrics::PRAUC(p, y),
      Brier = mean((p - y)^2))
  }

  ctrl <- caret::trainControl(
    method          = "cv",
    index           = fold_index,
    summaryFunction = binary_summary,
    classProbs      = TRUE,
    savePredictions = "final"
  )

  # ---- Train models ----
  form   <- as.formula(paste(target, "~ ."))
  n_vars <- length(variables)

  fit <- function(method, ...) {
    message("Training ", method)
    set.seed(42)
    caret::train(form, data = model_data, method = method,
                 metric = "ROC", trControl = ctrl, ...)
  }

  builders <- list(
    glmnet = function() {
      if (n_vars < 2) {
        fit("bayesglm", preProcess = c("center", "scale"))   # glmnet needs >= 2 predictors
      } else {
        fit("glmnet", preProcess = c("center", "scale"),
            tuneGrid = expand.grid(alpha  = seq(0, 1, length.out = 5),
                                   lambda = 10^seq(-3, 1, length.out = 20)))
      }
    },
    rf = function() {
      if (is.null(rf_grid)) {
        rf_grid <- expand.grid(
          mtry          = unique(pmax(1, floor(c(sqrt(n_vars), n_vars / 3, n_vars / 2, n_vars)))),
          splitrule     = c("gini", "extratrees"),
          min.node.size = c(1, 5, 10),
          stringsAsFactors = FALSE)
      }
      fit("ranger", importance = "permutation", tuneGrid = rf_grid)
    },
    xgb = function() {
      fit("xgbTree", verbose = 0,
          tuneGrid = expand.grid(nrounds = c(50, 100, 150), max_depth = c(2, 4, 6),
                                 eta = c(0.01, 0.1, 0.3), gamma = 0,
                                 colsample_bytree = 0.8, min_child_weight = 1,
                                 subsample = 0.8))
    },
    pls = function() {
      fit("pls", preProcess = c("center", "scale"),
          tuneGrid = expand.grid(ncomp = seq_len(min(10, n_vars))))
    }
  )

  models <- lapply(setNames(nm = which_models), function(m) builders[[m]]())
  display_names <- c(glmnet = "GLMnet", rf = "Random Forest", xgb = "XGBoost", pls = "PLS-DA")

  # ---- Out-of-fold predictions ----
  # savePredictions = "final" already keeps only the best tune, so no filtering needed
  prob_cols <- paste0("prob_", target_levels)

  predictions <- do.call(rbind, lapply(names(models), function(m) {
    p <- models[[m]]$pred
    data.frame(Model      = display_names[[m]],
               rowIndex   = p$rowIndex,
               obs_class  = p$obs,
               pred_class = p$pred,
               setNames(p[target_levels], prob_cols),
               row.names  = NULL)
  }))
  predictions$id_patient <- id_patient[predictions$rowIndex]

  # ---- Metrics ----
  # Mean of the per-fold metrics (model$resample = one row per held-out fold, best tune only).
  # CI for the fold mean uses the Nadeau-Bengio correction: folds share training data,
  # so the naive SD/sqrt(n) is far too narrow. Var = SD^2 * (1/n_folds + n_test/n_train).
  test_train_ratio <- mean(sapply(fold_index, function(tr) (nrow(model_data) - length(tr)) / length(tr)))
  fold_metrics <- function(model) {
    r  <- model$resample
    m  <- mean(r$ROC)
    se <- sd(r$ROC) * sqrt(1 / nrow(r) + test_train_ratio)
    tq <- qt(0.975, df = nrow(r) - 1)
    c(ROC_mean_fold = m, ROC_sd_fold = sd(r$ROC),
      ROC_fold_low = max(0, m - tq * se), ROC_fold_high = min(1, m + tq * se),
      AUPRC = mean(r$AUPRC), Brier = mean(r$Brier))
  }

  # One AUC over all out-of-fold predictions. With repeated CV each row appears
  # once per repeat, so average within row first (repeats aren't independent).
  pooled_auc <- function(label) {
    d <- predictions[predictions$Model == label, ]
    d <- aggregate(d[prob_cols[1]], by = d[c("rowIndex", "obs_class")], FUN = mean)
    r <- pROC::roc(d$obs_class, d[[prob_cols[1]]],
                   levels = rev(target_levels), direction = "<", quiet = TRUE)
    ci <- as.numeric(pROC::ci.auc(r))
    c(ROC = ci[2], ROC_low = ci[1], ROC_high = ci[3])
  }

  comparison <- data.frame(
    model_key = names(models),
    Model     = unname(display_names[names(models)]),
    t(sapply(names(models), function(m) {
      c(pooled_auc(display_names[[m]]), fold_metrics(models[[m]]))
    })),
    row.names = NULL
  )

  # ---- Return ----
  list(models         = models,
       predictions    = predictions,
       comparison     = comparison,
       variables_used = variables,
       target_used    = target,
       target_levels  = target_levels,
       n_resamples    = length(fold_index),
       glmnet_method  = models$glmnet$method)
}


# ---- Forward selection with Random Forest only ----
forward_selection_rf <- function(data, target, variables,
                                 positive_class   = NULL,
                                 selection_metric = c("ROC", "ROC_mean_fold", "AUPRC", "Brier"),
                                 max_vars         = NULL,
                                 min_improvement  = 0,
                                 cv_repeats       = 10,
                                 tune             = TRUE) {

  selection_metric <- match.arg(selection_metric)
  variables <- unique(variables)
  missing_cols <- setdiff(c(target, variables), names(data))
  if (length(missing_cols) > 0) stop("Not found in data: ", paste(missing_cols, collapse = ", "))

  # Drop incomplete rows once, so every candidate set uses the same rows and folds
  keep <- complete.cases(data[, c(target, variables)])
  if (any(!keep)) {
    message("Dropping ", sum(!keep), " row(s) with missing values.")
    data <- data[keep, , drop = FALSE]
  }
  groups     <- if ("id_patient" %in% names(data)) data$id_patient else seq_len(nrow(data))
  fold_index <- make_patient_folds(data[[target]], groups, repeats = cv_repeats)

  # RF performance of one variable set. tune = FALSE uses one fixed setting (much faster).
  score_set <- function(vars) {
    grid <- if (tune) NULL else data.frame(mtry = max(1, floor(sqrt(length(vars)))),
                                           splitrule = "gini", min.node.size = 10)
    res <- suppressMessages(train_binary_models(
      data, target, vars, positive_class,
      fold_index = fold_index, which_models = "rf", rf_grid = grid))
    res$comparison
  }
  sign <- if (selection_metric == "Brier") -1 else 1   # lower Brier is better

  selected    <- character(0)
  best_so_far <- -Inf
  steps <- candidates <- list()
  n_steps <- if (is.null(max_vars)) length(variables) else min(max_vars, length(variables))

  for (step in seq_len(n_steps)) {
    remaining <- setdiff(variables, selected)
    message(sprintf("Step %d: trying %d candidate(s)", step, length(remaining)))

    comp <- do.call(rbind, lapply(remaining, function(v) {
      data.frame(step = step, set_size = step, candidate_variable = v, score_set(c(selected, v)))
    }))
    candidates[[step]] <- comp

    best <- comp[which.max(sign * comp[[selection_metric]]), ]
    gain <- sign * best[[selection_metric]] - best_so_far
    if (step > 1 && gain < min_improvement) {
      message(sprintf("  Best gain %.4f < min_improvement - stopping.", gain))
      break
    }
    best_so_far <- sign * best[[selection_metric]]
    selected    <- c(selected, best$candidate_variable)
    message(sprintf("  + %s (%s = %.4f)", best$candidate_variable,
                    selection_metric, best[[selection_metric]]))

    steps[[step]] <- data.frame(
      best[c("step", "set_size")],
      variable_added = best$candidate_variable,
      selected_set   = paste(selected, collapse = " + "),
      best[!names(best) %in% c("step", "set_size", "candidate_variable")])
  }

  list(combined_steps        = data.frame(target = target, do.call(rbind, steps), row.names = NULL),
       combined_candidates   = data.frame(target = target, do.call(rbind, candidates), row.names = NULL),
       selected_vars         = selected,
       selection_metric_used = selection_metric,
       fold_index            = fold_index)
}