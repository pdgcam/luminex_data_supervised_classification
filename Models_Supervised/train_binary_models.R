# ---- Patient-grouped, stratified, repeated CV folds ----
# Patients (not rows) are assigned to folds, so no patient is ever in train and test
# at once. Patients are stratified by which classes their rows contain (only negative,
# only positive, or both), so this also works when a patient's outcome changes over time.
# Pass groups = row numbers for ordinary stratified CV.
make_patient_folds <- function(y, groups, k = 5, repeats = 1, seed = 42) {
  pats <- unique(groups)
  strat <- sapply(
    split(as.character(y), factor(groups, levels = pats)),
    function(v) paste(sort(unique(v)), collapse = "+")
  )
  set.seed(seed)
  unlist(lapply(seq_len(repeats), function(r) {
    pat_fold <- caret::createFolds(factor(strat), k = k, list = FALSE)
    row_fold <- pat_fold[match(groups, pats)]
    setNames(
      lapply(seq_len(k), function(f) which(row_fold != f)),
      sprintf("Fold%d.Rep%02d", seq_len(k), r)
    )
  }), recursive = FALSE)
}


train_binary_models <- function(data, target,
                                positive_class = NULL, k = 5,
                                cv_repeats = 10, ratio_or_cs = NULL,
                                which_models = c("glmnet", "rf", "naive_bayes"),
                                rf_grid = NULL, n_boot = 2000) {
  which_models <- match.arg(which_models, several.ok = TRUE)

  # ---- Inputs ----
  variables <- setdiff(names(data), c(target, "id_patient"))
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

  # pos class to be provided for DENV vs CHICK 
  # for rest, neg = no infection, pos = infection
  # positive class = level 1
  positive_class <- if (is.null(positive_class)) levels(y)[1] else make.names(positive_class)
  if (!positive_class %in% levels(y)) {
    stop("positive_class '", positive_class, "' not in: ", paste(levels(y), collapse = ", "))
  }
  model_data[[target]] <- relevel(y, ref = positive_class)
  target_levels <- levels(model_data[[target]])
  message("Positive class: ", target_levels[1], " | Negative: ", target_levels[2])


  # ---- Folds (shared across models) ----
  ratio_or_cs <- match.arg(ratio_or_cs, c("ratio", "cs"))

  groups <- if (ratio_or_cs == "ratio" || is.null(id_patient)) {
  seq_len(nrow(model_data))   # row-level folds: every row is its own group
  } else {
  id_patient     # patient-level folds: all rows of a patient stay together
  }

  fold_index <- make_patient_folds(model_data[[target]], k = k, groups, repeats = cv_repeats)

  # Every held-out fold needs both classes, or AUC/AUPRC can't be computed
  one_class <- sapply(fold_index, function(tr) length(unique(model_data[[target]][-tr])) < 2)
  if (any(one_class)) {
    yy <- model_data[[target]]
    n_pat <- sapply(levels(yy), function(l) length(unique(groups[yy == l])))
    n_row <- table(yy)
    stop(
      sum(one_class), " of ", length(fold_index), " test folds contain only one class. ",
      "Patients with >= 1 row of each class: ",
      paste(names(n_pat), n_pat, sep = " = ", collapse = ", "),
      " (rows: ", paste(names(n_row), n_row, sep = " = ", collapse = ", "), "). ",
      "5-fold CV needs at least 5 patients in each class."
    )
  }

  # ---- Per-fold metrics: ROC (= AUC) / Sens / Spec + AUPRC + Brier ----
  binary_summary <- function(data, lev, model = NULL) {
    p <- data[[lev[1]]] # prob of positive class
    y <- as.integer(data$obs == lev[1]) # 1 = positive
    c(caret::twoClassSummary(data, lev),
      AUPRC = MLmetrics::PRAUC(p, y),
      Brier = mean((p - y)^2)
    )
  }

  ctrl <- caret::trainControl(
    method          = "repeatedcv",
    number          = k,
    repeats         = cv_repeats,
    index           = fold_index,
    summaryFunction = binary_summary,
    classProbs      = TRUE,
    savePredictions = "final"
  )

  # ---- Train models ----
  form <- as.formula(paste(target, "~ ."))
  n_vars <- length(variables)

  fit <- function(method, ...) {
    message("Training ", method)
    set.seed(42)
    caret::train(form,
      data = model_data, method = method,
      metric = "ROC", trControl = ctrl, ...
    )
  }

  builders <- list(
    glmnet = function() {
      if (n_vars < 2) {
        fit("bayesglm", preProcess = c("center", "scale")) # glmnet needs >= 2 predictors
      } else {
        fit("glmnet",
          preProcess = c("center", "scale"),
          tuneGrid = expand.grid(
            alpha = seq(0, 1, length.out = 5),
            lambda = 10^seq(-3, 1, length.out = 20)
          )
        )
      }
    },
    rf = function() {
      if (is.null(rf_grid)) {
        rf_grid <- expand.grid(
          mtry = unique(pmax(1, floor(c(sqrt(n_vars), n_vars / 3, n_vars / 2, n_vars)))),
          splitrule = c("gini", "extratrees"),
          min.node.size = c(1, 5, 10),
          stringsAsFactors = FALSE
        )
      }
      fit("ranger", importance = "permutation", tuneGrid = rf_grid)
    },
    naive_bayes = function() {
      fit("naive_bayes",
        preProcess = c("center", "scale"),
        tuneGrid = expand.grid(
          usekernel = c(TRUE, FALSE),
          laplace = c(0, 1),
          adjust = c(0.5, 1, 1.5)
        )
      )
    }
  )

  models <- lapply(setNames(nm = which_models), function(m) builders[[m]]())
  display_names <- c(glmnet = "GLMnet", rf = "Random Forest", naive_bayes = "Naive Bayes")


  # ---- Metrics: pooled AUC + per-fold mean AUC, AUPRC, Brier ----
  #   pooled: computed once on all out-of-fold predictions
  #               (with repeats, each row's predictions are averaged first)
  #   mean_fold: mean of the per-fold values caret computed
  #               (model$resample = one row per held-out fold, best tune only)

  # ---- Out-of-fold predictions ----
  # savePredictions = "final" already keeps only the best tune, so no filtering needed
  prob_cols <- paste0("prob_", target_levels)

  predictions <- do.call(rbind, lapply(names(models), function(m) {
    p <- models[[m]]$pred
    data.frame(
      Model = display_names[[m]],
      rowIndex = p$rowIndex,
      obs_class = p$obs,
      pred_class = p$pred,
      setNames(p[target_levels], prob_cols),
      row.names = NULL
    )
  }))
  predictions$id_patient <- id_patient[predictions$rowIndex]

  # Mean of the per-fold metrics (model$resample = one row per held-out fold, best tune only).
  # CI for the fold mean uses the Nadeau-Bengio correction: folds share training data,
  # so the naive SD/sqrt(n) is far too narrow. Var = SD^2 * (1/n_folds + n_test/n_train).
  test_train_ratio <- mean(sapply(fold_index, function(tr) (nrow(model_data) - length(tr)) / length(tr)))
  fold_metrics <- function(model) {
    r <- model$resample
    roc <- r$ROC[!is.na(r$ROC)]
    n <- length(roc)
    m <- mean(roc)
    se <- sd(roc) * sqrt(1 / n + test_train_ratio)
    tq <- qt(0.975, df = n - 1)
    c(
      ROC_mean_fold = m,
      ROC_sd_fold = sd(roc),
      ROC_mean_fold_low = max(0, m - tq * se),
      ROC_mean_fold_high = min(1, m + tq * se),
      AUPRC_mean_fold = mean(r$AUPRC, na.rm = TRUE),
      Brier_mean_fold = mean(r$Brier, na.rm = TRUE),
      n_folds_used = n
    )
  }


  # -- Pooled summary --
  # AUC from Mann-Whitney ranks: identical to pROC's AUC (ties count 1/2),
  # and fast enough to recompute thousands of times in the bootstrap.
  auc_rank <- function(p, y) {
    n1 <- sum(y == 1)
    n0 <- sum(y == 0)
    if (n1 == 0 || n0 == 0) return(NA_real_)
    (sum(rank(p)[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
  }

  # Patient-level (cluster) bootstrap: resample whole patients with replacement and
  # keep all of each sampled patient's rows, so within-patient correlation is respected.
  # Same seed for every model, so all models are scored on the same resamples.
  cluster_boot_ci <- function(p, y, cluster, n_boot, seed = 42) {
    rows_by_cluster <- split(seq_along(p), cluster)
    set.seed(seed)
    boots <- replicate(n_boot, {
      pick <- sample.int(length(rows_by_cluster), replace = TRUE)
      i <- unlist(rows_by_cluster[pick], use.names = FALSE)
      auc_rank(p[i], y[i])
    })
    if (mean(is.na(boots)) > 0.05) {
      warning(
        round(100 * mean(is.na(boots))),
        "% of bootstrap resamples contained only one class and were dropped."
      )
    }
    quantile(boots, c(0.025, 0.975), na.rm = TRUE, names = FALSE)
  }

  # cs:    rows within a patient are correlated -> patient-level bootstrap CI
  # ratio: rows are independent -> DeLong CI from pROC
  pooled_metrics <- function(label) {
    d <- predictions[predictions$Model == label, ]
    d <- aggregate(d[prob_cols[1]], by = d[c("rowIndex", "obs_class")], FUN = mean)
    p <- d[[prob_cols[1]]]
    y <- as.integer(d$obs_class == target_levels[1]) # 1 = positive

    if (ratio_or_cs == "cs") {
      auc <- auc_rank(p, y)
      ci <- cluster_boot_ci(p, y, cluster = groups[d$rowIndex], n_boot = n_boot)
    } else {
      r <- pROC::roc(d$obs_class, p,
        levels = rev(target_levels), direction = "<", quiet = TRUE
      )
      ci3 <- as.numeric(pROC::ci.auc(r))
      auc <- ci3[2]
      ci <- ci3[c(1, 3)]
    }

    c(
      ROC_pooled = auc,
      ROC_pooled_low = ci[1],
      ROC_pooled_high = ci[2],
      AUPRC_pooled = MLmetrics::PRAUC(p, y),
      Brier_pooled = mean((p - y)^2)
    )
  }

  comparison <- data.frame(
    model_key = names(models),
    Model = unname(display_names[names(models)]),
    t(sapply(names(models), function(m) {
      c(pooled_metrics(display_names[[m]]), fold_metrics(models[[m]]))
    })),
    ROC_pooled_ci_method = if (ratio_or_cs == "cs") "patient bootstrap" else "DeLong",
    row.names = NULL
  )

  # ---- Return ----
  list(
    models = models,
    predictions = predictions,
    comparison = comparison,
    variables_used = variables,
    target_used = target,
    target_levels = target_levels,
    n_resamples = length(fold_index),
    glmnet_method = models$glmnet$method
  )
}

