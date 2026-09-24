
# extract variable importance for each antigen
library(ggplot2)
library(dplyr)
library(caret)
library(purrr)
library(glmnet)
library(purrr)
library(ggplot2)
library(ggrepel)
library(DALEX)
library(reshape2)   # for melt
library(scales)
# --- Source functions
source(here('Models/train_binary_models.R'))
source(here('Models/train_multinomial_models.R'))
source(here('Functions.R'))
source(here('Models/forward_selection.R'))

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




# ----- FORWARD SELECTION

# read data 

# --- Source functions
source(here('Models/train_binary_models.R'))
source(here('Models/train_multinomial_models.R'))
source(here('Functions.R'))


isotypes <- c("IgG", "IgA", "IgM", "avidity")

flavi_antigens <- c("DENV1_DIII","DENV1_NS1","DENV1_VLP","SHERPADES_DENV1_DIII",
"DENV2_DIII","DENV2_NS1","DENV2_VLP","SHERPADES_DENV2_DIII",
"DENV3_DIII","DENV3_NS1","DENV3_VLP", "SHERPADES_DENV3_DIII",
"DENV4_DIII","DENV4_NS1","DENV4_VLP", "SHERPADES_DENV4_DIII",
"JEV_E", "JEV_NS1", "SHERPADES_JEV_DIII",
"YFV_E", "YFV_NS1", "SHERPADES_YFV_DIII",
"WNV_DIII","WNV_NS1","SHERPADES_WNV_DIII",
"ZIKV_NS1","ZIKV_VLP","ZIKVAS_DIII","ZIKVSU_NS1","SHERPADES_ZIKV_DIII")

length(flavi_antigens) # 30 antigens

#  Read: Flavi vs not (binary) and dengue serotype (multinomial) forward selection

data_with_binomial_targets_flavi_list <- readRDS(here("Data", "model_prepared", "data_with_binomial_targets_flavi_list.rds"))
data_with_multinomial_targets_list <- readRDS(here("Data", "model_prepared", "data_with_multinomial_targets_list.rds"))


ratio_IgG_flavi_vs_not <- data_with_binomial_targets_flavi_list$ratio$IgG$data$flavi
ratio_IgA_flavi_vs_not <- data_with_binomial_targets_flavi_list$ratio$IgA$data$flavi
ratio_IgM_flavi_vs_not <- data_with_binomial_targets_flavi_list$ratio$IgM$data$flavi
ratio_avidity_flavi_vs_not <- data_with_binomial_targets_flavi_list$ratio$avidity$data$flavi


ratio_IgG_dengue_serotype <- data_with_multinomial_targets_list$ratio$IgG$data$dengue_serotype
ratio_IgA_dengue_serotype <- data_with_multinomial_targets_list$ratio$IgA$data$dengue_serotype
ratio_IgM_dengue_serotype <- data_with_multinomial_targets_list$ratio$IgM$data$dengue_serotype
ratio_avidity_dengue_serotype <- data_with_multinomial_targets_list$ratio$avidity$data$dengue_serotype




forward_selection_ratio_IgG_flavi_vs_not<- forward_selection(
   data  = ratio_IgG_flavi_vs_not,
   variables  = flavi_antigens,
   target  = "flavi_target",
   positive_class = data_with_binomial_targets_flavi_list$ratio$IgG$positive_class_map$flavi,
   metrics   = c("ROC", "AUPRC", "Brier"),
   selection_metric = "auto",    
   selection_model  = NULL,      
   max_vars   = NULL,      
   min_improvement  = -Inf          
 )

names(foward_selection_ratio_IgG_flavi_vs_not)


forward_selection_ratio_IgG_flavi_vs_not$selected_vars
forward_selection_ratio_IgG_flavi_vs_not$combined_steps                 # step-by-step performance (winner only)
forward_selection_ratio_IgG_flavi_vs_not$combined_candidates            # every candidate tried at every step


names(forward_selection_ratio_IgG_flavi_vs_notforward_selection_ratio_IgG_flavi_vs_not)

forward_selection_ratio_IgG_flavi_vs_not$combined_steps


plot_forward_selection <- function(forward_results,
                                    metric = NULL,
                                    target = NULL,
                                    spread = c("sd", "range", "none"),
                                    show_models = FALSE) {
  
  spread <- match.arg(spread)

  df <- forward_results$combined_steps
  if (nrow(df) == 0) stop("combined_steps is empty -- nothing to plot.")

  if (is.null(metric)) metric <- forward_results$selection_metric_used
  if (identical(metric, "auto")) {
    stop("selection_metric was 'auto' (differs by target type) -- please pass `metric` explicitly, e.g. metric = 'ROC' or metric = 'AUC_Micro'.")
  }

  if (!is.null(target)) df <- df[df$target == target, ]

  if (!metric %in% names(df)) {
    stop(sprintf("Metric '%s' not present in results. Available metric columns: %s",
                 metric,
                 paste(setdiff(names(df),
                               c("target", "step", "set_size", "variable_added",
                                 "selected_set", "Model", "candidate_variable",
                                 "engine")),
                       collapse = ", ")))
  }

  # ---- Average across models within each step ----
  ens <- df %>%
    dplyr::group_by(target, step, set_size, variable_added) %>%
    dplyr::summarise(
      mean_metric = mean(.data[[metric]], na.rm = TRUE),
      sd_metric   = stats::sd(.data[[metric]], na.rm = TRUE),
      min_metric  = min(.data[[metric]], na.rm = TRUE),
      max_metric  = max(.data[[metric]], na.rm = TRUE),
      n_models    = sum(!is.na(.data[[metric]])),
      models_used = paste(sort(Model[!is.na(.data[[metric]])]), collapse = ", "),
      .groups = "drop"
    ) %>%
    dplyr::arrange(target, step)

  # sd is NA when only one model contributed -- treat as zero width
  ens$sd_metric[is.na(ens$sd_metric)] <- 0

  if (any(ens$n_models != max(ens$n_models))) {
    warning(paste0(
      "Number of models differs between steps (",
      paste(sort(unique(ens$n_models)), collapse = " vs "),
      "). The mean is then taken over different model sets, so steps are not ",
      "strictly comparable -- check `models_used` in the returned data."
    ))
  }

  ens$step_label <- factor(
    paste0(ens$step, ". +", ens$variable_added),
    levels = unique(paste0(ens$step, ". +", ens$variable_added))
  )

  # ---- Build plot ----
  p <- ggplot2::ggplot(ens, ggplot2::aes(x = step_label, y = mean_metric, group = target))

  # faint per-model curves underneath, if asked for
  if (show_models) {
    df_models <- df
    df_models$step_label <- factor(
      paste0(df_models$step, ". +", df_models$variable_added),
      levels = levels(ens$step_label)
    )
    p <- p + ggplot2::geom_line(
      data = df_models,
      ggplot2::aes(x = step_label, y = .data[[metric]], group = Model),
      inherit.aes = FALSE,
      colour = "grey70", linewidth = 0.4, alpha = 0.8
    )
  }

  if (spread == "sd") {
    p <- p + ggplot2::geom_ribbon(
      ggplot2::aes(ymin = mean_metric - sd_metric, ymax = mean_metric + sd_metric),
      fill = "steelblue", alpha = 0.20
    )
  } else if (spread == "range") {
    p <- p + ggplot2::geom_ribbon(
      ggplot2::aes(ymin = min_metric, ymax = max_metric),
      fill = "steelblue", alpha = 0.20
    )
  }

  ribbon_label <- switch(spread,
    sd    = "ribbon = +/- 1 SD across models",
    range = "ribbon = min-max across models",
    none  = NULL
  )

  p <- p +
    ggplot2::geom_line(linewidth = 0.9, colour = "steelblue4") +
    ggplot2::geom_point(size = 2, colour = "steelblue4") +
    ggplot2::labs(
      x        = "Antigen added (cumulative)",
      y        = paste0("Mean ", metric, " across models"),
      title    = paste("Forward selection -- ensemble mean", metric),
      subtitle = ribbon_label
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

  if (length(unique(ens$target)) > 1) {
    p <- p + ggplot2::facet_wrap(~ target, scales = "free_x")
  }

  # attach the aggregated table so you can inspect / export the numbers
  attr(p, "ensemble_data") <- ens
  p
}

quartz()
plot_forward_selection(forward_selection_ratio_IgG_flavi_vs_not, metric = "ROC")


set.seed(42)
fold_index <- caret::createMultiFolds(
  factor(ratio_IgG_flavi_vs_not$flavi_target),   # use your actual target column
  k = 5, times = 10
)

length(fold_index)   # expect 50




# does WVN really work that well??
chk <- train_binary_models(
  data = ratio_IgG_flavi_vs_not, target = "flavi_target",
  variables = c("WNV_NS1", "YFV_E", "SHERPADES_YFV_DIII"), positive_class = "positive",
  fold_index = fold_index
)

chk$comparison




dat  <- antigen_data_ratio_IgG_flavi_vs_not
keep <- flavi_antigens

pairs_df <- t(combn(keep, 2)) |> as.data.frame() |>
  setNames(c("a", "b"))

long <- lapply(seq_len(nrow(pairs_df)), function(i) {
  a <- pairs_df$a[i]; b <- pairs_df$b[i]
  data.frame(pair = paste(a, "vs", b), x = dat[[a]], y = dat[[b]])
}) |> bind_rows()

long$pair <- factor(long$pair, levels = unique(long$pair))  # keeps combn order


antigen_pairs <- ggplot(long, aes(x, y)) +
  geom_point(size = 0.4, alpha = 0.4) +
  facet_wrap(~ pair, scales = "free", ncol = 12) +
  theme_minimal(base_size = 6) +
  theme(strip.text = element_text(size = 4), axis.title = element_blank())


ggsave(here::here("Results/antigen_pairs.pdf"), antigen_pairs,
       width = 30, height = 60, limitsize = FALSE)




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
