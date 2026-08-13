
# extract variable importance for each antigen
library(ggplot2)
library(dplyr)
library(caret)
library(purrr)
library(glmnet)
library(dplyr)
library(purrr)
library(ggplot2)
library(ggrepel)
library(DALEX)


load_binary_result <- function(comparison, isotype, dataset_type) {
  path <- file.path("Results", "Binary_Classification", dataset_type,
                     paste0(comparison, "_", isotype, ".rds"))
  if (!file.exists(path)) stop(paste("File not found:", path))
  readRDS(path)
}

# plot 1 - feature importance 
extract_varimp <- function(fit_result, model_name, n_permutations = 50) {
  model_obj <- fit_result$models[[model_name]]
  if (is.null(model_obj)) return(NULL)

  positive_class <- fit_result$target_levels[1]

  # caret::train() stores the training data it fit on (post any resampling
  # setup, pre-preprocessing — preProcess is reapplied automatically at predict time)
  train_data <- model_obj$trainingData
  y_fac  <- train_data$.outcome
  X      <- train_data[, setdiff(names(train_data), ".outcome"), drop = FALSE]
  y_num  <- as.numeric(y_fac == positive_class)

  # same predict_function works for all 4 model types because they're all
  # caret 'train' objects — this is what makes the importance comparable
  pred_fun <- function(object, newdata) {
    predict(object, newdata = newdata, type = "prob")[[positive_class]]
  }

  explainer <- DALEX::explain(
    model            = model_obj,
    data             = X,
    y                = y_num,
    predict_function = pred_fun,
    label            = model_name,
    verbose          = FALSE
  )

  vi <- DALEX::model_parts(
    explainer,
    type          = "difference",        # dropout_loss - loss of full model
    loss_function = DALEX::loss_one_minus_auc,
    B             = n_permutations,       # number of permutation repeats
    N             = NULL                  # use all rows, not a subsample
  )

  vi_df <- as.data.frame(vi)
  vi_df <- vi_df[!vi_df$variable %in% c("_baseline_", "_full_model_"), ]

  vi_df %>%
    group_by(variable) %>%
    summarise(Importance = mean(dropout_loss, na.rm = TRUE), .groups = "drop") %>%
    rename(Antigen = variable) %>%
    mutate(Model = model_name)
}



dataset_types <- c("ratio", "cross_sectional") 
isotypes <- c("IgG", "IgA", "IgM", "avidity")

feature_importance_list <- list()

for (dataset_type in dataset_types) {
  for (isotype in isotypes) {
    fit_result <- load_binary_result("dengue_vs_not", isotype, dataset_type)
    model_imps <- lapply(names(fit_result$models), function(m) extract_varimp(fit_result, m))
    model_imps <- bind_rows(model_imps)

    agg_imp <- model_imps %>%
      group_by(Model) %>%
      mutate(
        # normalized: each model's importances rescaled 0-1 by its own max,
        # so a model with small absolute drops doesn't get automatically
        # outweighted by a model with larger absolute drops
        Importance_norm = Importance / max(Importance, na.rm = TRUE),
        # rank: 1 = most important antigen *within that model*; ties averaged
        Importance_rank = rank(-Importance, ties.method = "average")
      ) %>%
      ungroup() %>%
      group_by(Antigen) %>%
      summarise(
        Importance_mean_raw  = mean(Importance, na.rm = TRUE),      # what you had before
        Importance_mean_norm = mean(Importance_norm, na.rm = TRUE), # scale-corrected
        Importance_mean_rank = mean(Importance_rank, na.rm = TRUE), # lower = more important
        .groups = "drop"
      ) %>%
      mutate(isotype = isotype, dataset_type = dataset_type)

    key <- paste(dataset_type, isotype, sep = "_")
    feature_importance_list[[key]] <- agg_imp
  }
}

all_importance <- bind_rows(feature_importance_list)
importance_plots <- list()

# pick which aggregation to plot: "Importance_mean_norm" (recommended default),
# "Importance_mean_rank" (lower = more important - flip sign if using this),
# or "Importance_mean_raw" (your original, kept for comparison)
importance_metric <- "Importance_mean_norm"

for (dtype in dataset_types) {
  for (iso in isotypes) {

    df_sub <- all_importance %>%
      filter(dataset_type == dtype, isotype == iso)

    p <- ggplot(df_sub, aes(x = reorder(Antigen, -.data[[importance_metric]]),
                             y = .data[[importance_metric]])) +
      geom_col(fill = "#1d486b") +
      labs(
        title = paste("Feature Importance -", iso, paste0("(", dtype, ", Dengue vs Not)")),
        x = "Antigen",
        y = "Importance (normalized, averaged across models)"
      ) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 20),
            axis.text.y  = element_text(size = 20),
            axis.title.y = element_text(size = 20),
            aspect.ratio = 0.7)

    key <- paste(dtype, iso, sep = "_")
    importance_plots[[key]] <- p
  }
}

# example access
importance_plots[["ratio_IgG"]]
importance_plots[["cross_sectional_IgG"]]




# plot 2 - backward selection
backward_selection_avg_model <- function(dataset_type, isotype, antigen_set = flavi_antigens) {

  binom_flavi    <- data_with_binomial_targets_flavi_list[[dataset_type]][[isotype]]
  dengue_data    <- binom_flavi$data$dengue
  positive_class <- binom_flavi$positive_class_map$dengue

  get_avg_auc <- function(vars) {
    res <- train_binary_models(
      data           = dengue_data,
      target         = "dengue_target",
      variables      = vars,
      positive_class = positive_class,
      metrics        = c("ROC", "AUPRC", "Brier")
    )
    list(
      avg_auc   = mean(res$comparison$ROC, na.rm = TRUE),
      sd_auc    = sd(res$comparison$ROC, na.rm = TRUE),
      per_model = setNames(res$comparison$ROC, res$comparison$Model)
    )
  }

  # ---- baseline: pull from saved result rather than retraining ----
  baseline_result <- load_binary_result("dengue_vs_not", isotype, dataset_type)
  baseline_avg_auc <- mean(baseline_result$comparison$ROC, na.rm = TRUE)
  baseline_sd_auc  <- sd(baseline_result$comparison$ROC, na.rm = TRUE)

  cat(sprintf("\n=== Baseline (%s, %s): all %d antigens ===\n",
              isotype, dataset_type, length(antigen_set)))
  cat(sprintf("Baseline avg AUC across models: %.3f (sd %.3f)\n",
              baseline_avg_auc, baseline_sd_auc))

  drop_results <- lapply(antigen_set, function(ag) {
    cat(sprintf("Dropping: %s\n", ag))
    vars_dropped <- setdiff(antigen_set, ag)
    res <- get_avg_auc(vars_dropped)

    data.frame(
      antigen_dropped = ag,
      avg_AUC         = res$avg_auc,
      sd_AUC          = res$sd_auc,
      delta_AUC       = res$avg_auc - baseline_avg_auc   # negative = antigen was useful
    )
  })

  drop_df <- bind_rows(drop_results) %>%
    mutate(baseline_avg_AUC = baseline_avg_auc, isotype = isotype, dataset_type = dataset_type)

  list(baseline_avg_auc = baseline_avg_auc, baseline_sd_auc = baseline_sd_auc, drop_df = drop_df)
}

save_backward_selection <- function(result, dataset_type, isotype,
                                     out_dir = file.path("Results", "Backward_Selection")) {

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  cat(sprintf("Saved: %s\n", out_path))

  invisible(result)
}

# Run and save 
backward_avg_results <- list()

for (dataset_type in dataset_types) {
  for (isotype in isotypes) {

    cat(sprintf("\n>>> Running backward selection: %s / %s\n", dataset_type, isotype))

    result <- backward_selection_avg_model(dataset_type = dataset_type, isotype = isotype)
    save_backward_selection(result, dataset_type = dataset_type, isotype = isotype)

    key <- paste(dataset_type, isotype, sep = "_")
    backward_avg_results[[key]] <- result
  }
}

# load results
load_backward_selection <- function(dataset_type, isotype,
                                     out_dir = file.path("Results", "Backward_Selection")) {
  path <- file.path(out_dir, paste0("backward_avg_", dataset_type, "_", isotype, ".rds"))
  if (!file.exists(path)) stop(paste("File not found:", path))
  readRDS(path)
}

res <- load_backward_selection("ratio", "IgG")
# plot
plot_backward_selection <- function(dataset_type, isotype) {

  res <- load_backward_selection(dataset_type, isotype)

  ggplot(res$drop_df, aes(x = reorder(antigen_dropped, avg_AUC), y = avg_AUC)) +
    geom_point(aes(color = delta_AUC), size = 3) +
    geom_errorbar(aes(ymin = avg_AUC - sd_AUC, ymax = avg_AUC + sd_AUC), width = 0.2, alpha = 0.4) +
    geom_hline(yintercept = res$baseline_avg_auc, linetype = "dashed", color = "red") +
    scale_color_gradient2(low = "firebrick", mid = "grey80", high = "steelblue", midpoint = 0,
                           name = "Δ AUC\n(vs baseline)") +
    coord_flip() +
    labs(
      title    = paste0("Backward Selection (avg AUC across models): ", isotype,
                         " (", dataset_type, "), Dengue vs Not"),
      x = "Antigen Dropped", y = "Mean AUC across models (model without this antigen)"
    ) +
    theme_bw() +
    theme(axis.text.x = element_text(size = 20),
          axis.title.x = element_text(size = 20),
          axis.text.y      = element_text(size = 20),
          axis.title.y     = element_text(size = 20), aspect.ratio = 1.4)
}

# instant, no retraining
plot_backward_selection("ratio", "IgG")
plot_backward_selection("cross_sectional", "IgG")




# Run NEW Backward selection just for one isotype and dataset type 
# Run and save - single combination only (ratio, IgG) for timing/testing
# plot 2 - true backward (stepwise) selection
NEW_backward_selection_avg_model <- function(dataset_type, isotype, antigen_set = flavi_antigens,
                                          min_vars = 2, stop_on_auc_drop = NULL) {

  binom_flavi    <- data_with_binomial_targets_flavi_list[[dataset_type]][[isotype]]
  dengue_data    <- binom_flavi$data$dengue
  positive_class <- binom_flavi$positive_class_map$dengue

  get_avg_auc <- function(vars) {
    res <- train_binary_models(
      data           = dengue_data,
      target         = "dengue_target",
      variables      = vars,
      positive_class = positive_class,
      metrics        = c("ROC", "AUPRC", "Brier")
    )
    list(
      avg_auc   = mean(res$comparison$ROC, na.rm = TRUE),
      sd_auc    = sd(res$comparison$ROC, na.rm = TRUE),
      per_model = setNames(res$comparison$ROC, res$comparison$Model)
    )
  }

  # ---- baseline: pull from saved result rather than retraining ----
  baseline_result  <- load_binary_result("dengue_vs_not", isotype, dataset_type)
  baseline_avg_auc <- mean(baseline_result$comparison$ROC, na.rm = TRUE)
  baseline_sd_auc  <- sd(baseline_result$comparison$ROC, na.rm = TRUE)

  cat(sprintf("\n=== Baseline (%s, %s): all %d antigens ===\n",
              isotype, dataset_type, length(antigen_set)))
  cat(sprintf("Baseline avg AUC across models: %.3f (sd %.3f)\n",
              baseline_avg_auc, baseline_sd_auc))

  current_vars    <- antigen_set
  current_auc     <- baseline_avg_auc
  elimination_log <- list()
  round_num       <- 0

  # keep dropping until we hit min_vars, or (optionally) until dropping
  # anything further would cost too much AUC
  while (length(current_vars) > min_vars) {
    round_num <- round_num + 1
    cat(sprintf("\n--- Round %d: %d antigens remaining ---\n", round_num, length(current_vars)))

    # try dropping each remaining antigen, one at a time, from the CURRENT set
    round_results <- lapply(current_vars, function(ag) {
      cat(sprintf("  Trying drop: %s\n", ag))
      vars_dropped <- setdiff(current_vars, ag)
      res <- get_avg_auc(vars_dropped)
      data.frame(
        antigen_dropped = ag,
        avg_AUC         = res$avg_auc,
        sd_AUC          = res$sd_auc,
        delta_AUC       = res$avg_auc - current_auc   # relative to THIS round's starting point
      )
    })

    round_df <- bind_rows(round_results)

    # pick the antigen whose removal hurts least (or helps most) -
    # i.e. the row with the highest resulting avg_AUC
    best_idx    <- which.max(round_df$avg_AUC)
    best_drop   <- round_df$antigen_dropped[best_idx]
    new_auc     <- round_df$avg_AUC[best_idx]
    auc_change  <- new_auc - current_auc

    cat(sprintf("  -> Dropping '%s' (AUC %.3f -> %.3f, delta %+.3f)\n",
                best_drop, current_auc, new_auc, auc_change))

    # optional stopping rule: stop BEFORE dropping if it would cost too much AUC
    if (!is.null(stop_on_auc_drop) && auc_change < -stop_on_auc_drop) {
      cat(sprintf("  Stopping: further removal would drop AUC by more than %.3f\n",
                  stop_on_auc_drop))
      break
    }

    round_df$round          <- round_num
    round_df$n_vars_before  <- length(current_vars)
    round_df$vars_remaining <- paste(current_vars, collapse = ";")
    round_df$antigen_removed_this_round <- best_drop
    round_df$auc_after_removal          <- new_auc

    elimination_log[[round_num]] <- round_df

    # commit: remove the chosen antigen, update running AUC, move to next round
    current_vars <- setdiff(current_vars, best_drop)
    current_auc  <- new_auc
  }

  elimination_df <- bind_rows(elimination_log) %>%
    mutate(isotype = isotype, dataset_type = dataset_type,
           baseline_avg_AUC = baseline_avg_auc)

  # the actual elimination PATH - one row per round, in the order antigens were dropped
  elimination_path <- elimination_df %>%
    filter(antigen_dropped == antigen_removed_this_round) %>%
    dplyr::select(round, antigen_removed_this_round, n_vars_before,
                  auc_after_removal, isotype, dataset_type) %>%
    mutate(delta_from_baseline = auc_after_removal - baseline_avg_auc)

  list(
    baseline_avg_auc = baseline_avg_auc,
    baseline_sd_auc  = baseline_sd_auc,
    elimination_df   = elimination_df,     # full detail: every candidate tried, every round
    elimination_path = elimination_path,   # the ranked drop order you actually want
    final_vars       = current_vars        # whatever's left when the loop stops
  )
}

backward_avg_results <- list()
dataset_type <- "ratio"
isotype <- "IgG"

cat(sprintf("\n>>> Running backward selection: %s / %s\n", dataset_type, isotype))

start_time <- Sys.time()
result <- NEW_backward_selection_avg_model(dataset_type = dataset_type, isotype = isotype)
elapsed <- Sys.time() - start_time
cat(sprintf("\nElapsed time: %.1f %s\n", as.numeric(elapsed), units(elapsed)))

key <- paste(dataset_type, isotype, sep = "_")
backward_avg_results[[key]] <- result

dir.create(file.path("Results", "Backward_Selection_Sequential"),
           recursive = TRUE, showWarnings = FALSE)

saveRDS(result,
        file.path("Results", "Backward_Selection_Sequential",
                  "backward_sequential_ratio_IgG.rds"))
result <- readRDS(here("Results/Backward_Selection_Sequential/backward_sequential_ratio_IgG.rds"))

result$final_vars

ggplot(result$elimination_path, aes(x = reorder(antigen_removed_this_round, -round),
                 y = auc_after_removal)) +
    geom_point(aes(color = delta_from_baseline), size = 4) +
    geom_line(aes(group = 1), color = "#343434", alpha = 0.7) +
    geom_hline(yintercept = result$baseline_avg_auc, linetype = "dashed", color = "red") +
    scale_color_gradient2(low = "#f6a5a5", mid = "#f00202", high = "#1d3f5b", midpoint = 0,
                           name = "Δ AUC\n(vs baseline)") +
    labs(
      x = "Antigen Removed", y = "Mean AUC across models (after this antigen's removal)"
    ) +
    theme_bw() +
    theme(axis.text.x    = element_text(size = 20, angle = 90, hjust = 1),
          axis.title.x   = element_text(size = 20),
          axis.text.y    = element_text(size = 20),
          axis.title.y   = element_text(size = 20),
          aspect.ratio   = 0.7)

setdiff(baseline_result$variables_used, flavi_antigens)
setdiff(flavi_antigens, baseline_result$variables_used)

a <- get_avg_auc(flavi_antigens)$avg_auc
b <- get_avg_auc(flavi_antigens)$avg_auc
c <- get_avg_auc(flavi_antigens)$avg_auc
c(a, b, c)


  






# SHAP 
install.packages("fastshap")
library(fastshap)


extract_shap_tree <- function(fit_result, model_name) {
  if (!model_name %in% c("rf", "xgb")) {
    stop("This function only supports tree-based models: 'rf' or 'xgb'")
  }

  model_obj <- fit_result$models[[model_name]]
  if (is.null(model_obj)) return(NULL)

  positive_class <- fit_result$target_levels[1]
  train_data <- model_obj$trainingData
  X_raw <- train_data[, setdiff(names(train_data), ".outcome"), drop = FALSE]

  X_proc <- if (!is.null(model_obj$preProcess)) {
    predict(model_obj$preProcess, X_raw)
  } else {
    X_raw
  }

  if (model_name == "xgb") {
    booster <- model_obj$finalModel
    X_mat <- as.matrix(X_proc)
    contrib <- predict(booster, X_mat, predcontrib = TRUE)
    contrib <- contrib[, setdiff(colnames(contrib), "BIAS"), drop = FALSE]
    sv <- shapviz::shapviz(contrib, X = X_proc)

  } else if (model_name == "rf") {
    pred_fun_rf <- function(object, newdata) {
      predict(object, newdata = newdata, type = "prob")[[positive_class]]
    }

    shap_vals <- fastshap::explain(
      model_obj,
      X            = X_proc,
      pred_wrapper = pred_fun_rf,
      nsim         = 50,     # Monte Carlo repeats - more = less noisy, slower
      adjust       = TRUE
    )

    sv <- shapviz::shapviz(as.matrix(shap_vals), X = X_proc)
  }

  data.frame(
    Antigen    = colnames(sv$S),
    Importance = colMeans(abs(sv$S)),
    Model      = model_name,
    sv         = I(list(sv))
  )
}

shap_imps <- lapply(c("rf", "xgb"), function(m) extract_shap_tree(fit_result, m))
shap_imps <- bind_rows(shap_imps)





















  # correlation between antigens

library(dplyr)
library(ggplot2)
library(reshape2)   # for melt

plot_antigen_correlation <- function(dataset_type, isotype, antigen_set = flavi_antigens) {

  # pull the same source data used for modelling (pre-target-selection, logged values)
  df_logged <- get_source_df(dataset_type, isotype)

  antigen_data <- df_logged %>%
    dplyr::select(all_of(antigen_set)) %>%
    as.data.frame()

  cor_mat <- cor(antigen_data, use = "pairwise.complete.obs", method = "spearman")

  cor_df <- reshape2::melt(cor_mat, varnames = c("Antigen1", "Antigen2"), value.name = "correlation")

  # keep antigens in their original (input) order - no clustering
  cor_df <- cor_df %>%
    mutate(
      Antigen1 = factor(Antigen1, levels = antigen_set),
      Antigen2 = factor(Antigen2, levels = antigen_set)
    )

  print(cor_df)

  ggplot(cor_df, aes(x = Antigen1, y = Antigen2, fill = correlation)) +
    geom_tile() +
    geom_text(aes(label = sprintf("%.1f", correlation)), size = 3) +
    scale_fill_gradient2(low = "firebrick", mid = "white", high = "steelblue",
                          midpoint = 0, limits = c(-1, 1), name = "Spearman\ncorrelation") +
    labs(
      title = paste("Antigen Correlation Matrix -", isotype, paste0("(", dataset_type, ")")),
      x = NULL, y = NULL
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 15),
      axis.text.y = element_text(size = 15),
      panel.grid = element_blank(), aspect.ratio = 1
    )
}

# example
plot_antigen_correlation(dataset_type = "ratio", isotype = "IgG")


















# OLD


ggplot(igg_backward_avg$drop_df, aes(x = reorder(antigen_dropped, avg_AUC), y = avg_AUC)) +
  geom_point(aes(color = delta_AUC), size = 3) +
  geom_errorbar(aes(ymin = avg_AUC - sd_AUC, ymax = avg_AUC + sd_AUC), width = 0.2, alpha = 0.4) +
  geom_hline(yintercept = igg_backward_avg$baseline_avg_auc, linetype = "dashed", color = "red") +
  scale_color_gradient2(low = "firebrick", mid = "grey80", high = "steelblue", midpoint = 0,
                         name = "Δ AUC\n(vs baseline)") +
  coord_flip() +
  labs(
    title    = "Backward Selection (avg AUC across models)",
    subtitle = "Dashed red line = baseline avg AUC (from saved model); error bars = SD across the 4 models",
    x = "Antigen Dropped", y = "Mean AUC across models (model without this antigen)"
  ) +
  theme_bw()


igg_backward <- backward_selection_single_drop(dataset_type = "ratio", isotype = "IgG")

ggplot(igg_backward$drop_df, aes(x = reorder(antigen_dropped, AUC), y = AUC)) +
  geom_point(aes(color = delta_AUC), size = 3) +
  geom_hline(yintercept = igg_backward$baseline_auc, linetype = "dashed", color = "red") +
  scale_color_gradient2(low = "firebrick", mid = "grey80", high = "steelblue", midpoint = 0,
                         name = "Δ AUC\n(vs baseline)") +
  coord_flip() +
  labs(
    title    = "Backward Selection: Effect of Dropping Each Antigen (IgG, Dengue vs Not)",
    subtitle = "Dashed red line = baseline AUC with all antigens included",
    x = "Antigen Dropped", y = "AUC (model without this antigen)"
  ) +
  theme_bw()
