
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
library(here)
# --- Source functions
source(here('Models/NEW_train_binary_models.R'))
source(here('Models/train_multinomial_models.R'))
source(here('Models/RF_forward_selection.R'))

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
isotypes <- c("IgG", "IgA", "IgM")

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
source(here('Models/RF_forward_selection.R'))

isotypes <- c("IgG", "IgA", "IgM")

flavi_antigens <- c("DENV1_DIII","DENV1_NS1","DENV1_VLP","SHERPADES_DENV1_DIII",
"DENV2_DIII","DENV2_NS1","DENV2_VLP","SHERPADES_DENV2_DIII",
"DENV3_DIII","DENV3_NS1","DENV3_VLP", "SHERPADES_DENV3_DIII",
"DENV4_DIII","DENV4_NS1","DENV4_VLP", "SHERPADES_DENV4_DIII",
"JEV_E", "JEV_NS1", "SHERPADES_JEV_DIII",
"YFV_E", "YFV_NS1", "SHERPADES_YFV_DIII",
"WNV_DIII","WNV_NS1","SHERPADES_WNV_DIII",
"ZIKV_NS1","ZIKV_VLP","ZIKV_DIII","ZIKVSU_NS1","SHERPADES_ZIKV_DIII")

alpha_antigens <- c("CHIKV_E2", "CHIKV_NSP123", "CHIKV_VLP", "SHERPADES_CHIKV_E2", 
                    "MAYV_E2" , "SHERPADES_MAYV_E2",
                    "ONNV_E2", "ONNV_VLP",
                    "RR" , "SHERPADES_RR")



#  Read: Flavi vs not (binary) and dengue serotype (multinomial) forward selection

data_with_binomial_targets_flavi_list <- readRDS(here("Data", "model_prepared", "data_with_binomial_targets_flavi_list.rds"))
data_with_binomial_targets_chik_list <- readRDS(here("Data", "model_prepared", "data_with_binomial_targets_chik_list.rds"))
data_with_multinomial_targets_list <- readRDS(here("Data", "model_prepared", "data_with_multinomial_targets_list.rds"))


# ratio binomial 
ratio_IgG_flavi_vs_not <- data_with_binomial_targets_flavi_list$ratio$IgG$data$flavi
ratio_IgG_dengue_vs_not <- data_with_binomial_targets_flavi_list$ratio$IgG$data$dengue
ratio_IgG_chik_vs_not <- data_with_binomial_targets_chik_list$ratio$IgG$data$chik
# cross-sectional binomial 
cs_IgG_flavi_vs_not <- data_with_binomial_targets_flavi_list$cross_sectional$IgG$data$flavi
cs_IgG_dengue_vs_not <- data_with_binomial_targets_flavi_list$cross_sectional$IgG$data$dengue
cs_IgG_chik_vs_not <- data_with_binomial_targets_chik_list$cross_sectional$IgG$data$chik


# ratio multinomial 
ratio_IgG_dengue_serotype <- data_with_multinomial_targets_list$ratio$IgG$data$dengue_serotype
# cross-sectional binomial 
cs_IgG_dengue_serotype <- data_with_multinomial_targets_list$cross_sectional$IgG$data$dengue_serotype

# forward selection - only Random Forest 
# ---- Positive classes -------------------------------------------------------
pc_ratio_flavi  <- data_with_binomial_targets_flavi_list$ratio$IgG$positive_class_map$flavi
pc_ratio_dengue <- data_with_binomial_targets_flavi_list$ratio$IgG$positive_class_map$dengue
pc_ratio_chik   <- data_with_binomial_targets_chik_list$ratio$IgG$positive_class_map$chik
 
pc_cs_flavi     <- data_with_binomial_targets_flavi_list$cross_sectional$IgG$positive_class_map$flavi
pc_cs_dengue    <- data_with_binomial_targets_flavi_list$cross_sectional$IgG$positive_class_map$dengue
pc_cs_chik      <- data_with_binomial_targets_chik_list$cross_sectional$IgG$positive_class_map$chik


# =============================================================================
# RATIO
# =============================================================================

forward_selection_ratio_IgG_flavi_vs_not <- forward_selection_rf(
  data             = ratio_IgG_flavi_vs_not,
  target           = "flavi_target",
  variables        = flavi_antigens,
  positive_class   = pc_ratio_flavi,
  metrics          = c("ROC", "AUPRC", "Brier"),
  selection_metric = "auto",
  max_vars         = NULL,
  min_improvement  = -Inf
)
 
forward_selection_ratio_IgG_dengue_vs_not <- forward_selection_rf(
  data             = ratio_IgG_dengue_vs_not,   # was ratio_IgG_dengue_serotype - see note
  target           = "dengue_target",
  variables        = flavi_antigens,
  positive_class   = pc_ratio_dengue,
  metrics          = c("ROC", "AUPRC", "Brier"),
  selection_metric = "auto",
  max_vars         = NULL,
  min_improvement  = -Inf
)
 
forward_selection_ratio_IgG_chik_vs_not <- forward_selection_rf(
  data             = ratio_IgG_chik_vs_not,
  target           = "chik_target",
  variables        = alpha_antigens,
  positive_class   = pc_ratio_chik,
  metrics          = c("ROC", "AUPRC", "Brier"),
  selection_metric = "auto",
  max_vars         = NULL,
  min_improvement  = -Inf
)
 
forward_selection_ratio_IgG_dengue_serotype <- forward_selection_rf(
  data             = ratio_IgG_dengue_serotype,
  target           = "dengue_serotype_target",
  variables        = flavi_antigens,
  metrics          = c("AUC_Macro", "AUC_Micro", "Brier"),
  selection_metric = "auto",
  max_vars         = NULL,
  min_improvement  = -Inf
)
 
 
# =============================================================================
# CROSS-SECTIONAL
# =============================================================================
 
forward_selection_cs_IgG_flavi_vs_not <- forward_selection_rf(
  data             = cs_IgG_flavi_vs_not,
  target           = "flavi_target",
  variables        = flavi_antigens,
  positive_class   = pc_cs_flavi,
  metrics          = c("ROC", "AUPRC", "Brier"),
  selection_metric = "auto",
  max_vars         = NULL,
  min_improvement  = -Inf
)

forward_selection_cs_IgG_dengue_vs_not <- forward_selection_rf(
  data             = cs_IgG_dengue_vs_not,
  target           = "dengue_target",
  variables        = flavi_antigens,
  positive_class   = pc_cs_dengue,
  metrics          = c("ROC", "AUPRC", "Brier"),
  selection_metric = "auto",
  max_vars         = NULL,
  min_improvement  = -Inf
)

cs_IgG_chik_vs_not

forward_selection_cs_IgG_chik_vs_not <- forward_selection_rf(
  data             = cs_IgG_chik_vs_not,
  target           = "chik_target",
  variables        = alpha_antigens,
  positive_class   = pc_cs_chik,
  metrics          = c("ROC", "AUPRC", "Brier"),
  selection_metric = "auto",
  max_vars         = NULL,
  min_improvement  = -Inf
)
 
forward_selection_cs_IgG_dengue_serotype <- forward_selection_rf(
  data             = cs_IgG_dengue_serotype,
  target           = "dengue_serotype_target",
  variables        = flavi_antigens,
  metrics          = c("AUC_Macro", "AUC_Micro", "Brier"),
  selection_metric = "auto",
  max_vars         = NULL,
  min_improvement  = -Inf
)





# ---- NEW RANDOM FOREST CS FORWARD SELECTION 
fs <- forward_selection_rf(binom_flavi$data$dengue, "dengue_target",
                           variables = setdiff(names(binom_flavi$data$dengue), c("dengue_target", "id_patient")),
                           positive_class = binom_flavi$positive_class_map$dengue)



fs_chik <- forward_selection_rf(cs_IgG_chik_vs_not, "chik_target",
                           variables = setdiff(names(cs_IgG_chik_vs_not$chik_target), c("chik_target", "id_patient")),
                           positive_class = pc_cs_chik)

dir.create(file.path("Results", "NEW_RF_Forward_Selection"), recursive = TRUE, showWarnings = FALSE)



saveRDS(
  fs,
  here("Results", "NEW_RF_Forward_Selection",
       "forward_selection_cs_IgG_dengue_vs_not.rds")
)


 
# read forward selection 
forward_selection_ratio_IgG_flavi_vs_not <- readRDS(here('Results/Forward_Selection/forward_selection_ratio_IgG_flavi_vs_not.rds'))
forward_selection_ratio_IgG_dengue_vs_not <- readRDS(here('Results/Forward_Selection/forward_selection_ratio_IgG_dengue_vs_not.rds'))
forward_selection_ratio_IgG_chik_vs_not <- readRDS(here('Results/Forward_Selection/forward_selection_ratio_IgG_chik_vs_not.rds'))


forward_selection_cs_IgG_dengue_vs_not <- readRDS(here('Results/RF_Forward_Selection/forward_selection_cs_IgG_dengue_vs_not.rds'))


forward_selection_cs_IgG_flavi_vs_not <- readRDS(
  here("Results", "RF_Forward_Selection",
       "forward_selection_cs_IgG_flavi_vs_not.rds"))
 
forward_selection_cs_IgG_dengue_vs_not <- readRDS(
  here("Results", "RF_Forward_Selection",
       "forward_selection_cs_IgG_dengue_vs_not.rds"))
 
forward_selection_cs_IgG_chik_vs_not <- readRDS(
  here("Results", "RF_Forward_Selection",
       "forward_selection_cs_IgG_chik_vs_not.rds"))
 
forward_selection_cs_IgG_dengue_serotype <- readRDS(
  here("Results", "RF_Forward_Selection",
       "forward_selection_cs_IgG_dengue_serotype.rds"))
 

plot_forward_selection_rf <- function(forward_results,
                                      metric = NULL,
                                      target = NULL,
                                      spread = c("sd", "none")) {

  spread <- match.arg(spread)

  df <- forward_results$combined_steps
  if (nrow(df) == 0) stop("combined_steps is empty -- nothing to plot.")

  if (is.null(metric)) metric <- forward_results$selection_metric_used
  if (!is.null(target)) df <- df[df$target == target, ]

  if (!metric %in% names(df)) {
    stop(sprintf("Metric '%s' not present in results. Available: %s",
                 metric,
                 paste(grep("_sd$", setdiff(names(df),
                       c("target", "step", "set_size", "variable_added",
                         "selected_set", "model_key", "Model",
                         "candidate_variable")),
                       value = TRUE, invert = TRUE), collapse = ", ")))
  }

  sd_col <- paste0(metric, "_sd")
  has_sd <- spread == "sd" && sd_col %in% names(df)
  if (spread == "sd" && !has_sd) {
    warning(sprintf("No '%s' column - drawing the line without a ribbon.", sd_col))
  }

  # one row per step (RF only), so just order and label
  df <- df %>%
    dplyr::arrange(target, step) %>%
    dplyr::mutate(
      value      = .data[[metric]],
      sd_value   = if (has_sd) .data[[sd_col]] else 0,
      step_label = factor(paste0(step, ". +", variable_added),
                          levels = unique(paste0(step, ". +", variable_added)))
    )

  p <- ggplot2::ggplot(df, ggplot2::aes(x = step_label, y = value, group = target))

  if (has_sd) {
    p <- p + ggplot2::geom_ribbon(
      ggplot2::aes(ymin = value - sd_value, ymax = value + sd_value),
      fill = "steelblue", alpha = 0.20)
  }

  p <- p +
    ggplot2::geom_line(linewidth = 0.9, colour = "steelblue4") +
    ggplot2::geom_point(size = 2, colour = "steelblue4") +
    ggplot2::labs(
      x        = "Antigen added",
      y        = metric,
      title    = paste("Forward selection (Random Forest) ", metric),
      subtitle = if (has_sd) "ribbon = +/- 1 SD across resamples" else NULL) +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

  if (length(unique(df$target)) > 1) {
    p <- p + ggplot2::facet_wrap(~ target, scales = "free_x")
  }

  attr(p, "plot_data") <- df
  p
}


# ---- Usage ------------------------------------------------------------------
# plot_forward_selection_rf(forward_selection_cs_IgG_flavi_vs_not)



quartz()
p <- plot_forward_selection_rf(fs, metric = "ROC")
ggplot2::ggsave(here("Results/NEW_RF_Forward_Selection/forward_selection_roc.png"),
                p, width = 12, height = 6, dpi = 300)

quartz()
plot_forward_selection(forward_selection_ratio_IgG_chik_vs_not, metric = "ROC")
quartz()
plot_forward_selection(forward_selection_ratio_IgG_dengue_serotype, metric = "AUC_Macro")


