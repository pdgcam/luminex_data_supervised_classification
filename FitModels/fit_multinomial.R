
# Multinomial Classification 
results_dengue_serotype_list   <- list(ratio = list(), cross_sectional = list())
results_dengue_serotype_neg_list <- list(ratio = list(), cross_sectional = list())


for (dataset_type in dataset_types) {
  for (isotype in isotypes) {

    multinom <- data_with_multinomial_targets_list[[dataset_type]][[isotype]]

    # Dengue serotype (infected-only: DENV1-4)
    results_dengue_serotype <- train_multinomial_models(
      data = multinom$data$dengue_serotype,
      target = "dengue_serotype_target",
      variables = NULL,
      metrics = c("AUC_Macro", "AUC_Micro", "AUPRC_Macro", "Brier", "StratBrier")
    )

    # Dengue serotype + negative class
    #results_dengue_serotype_neg <- train_multinomial_models(
     # data = multinom$data$dengue_serotype_neg,
     # variables = NULL,
      #target = "dengue_serotype_neg_target",
     # metrics = c("ROC", "AUPRC", "Brier", "StratBrier")
    #)

    results_dengue_serotype_list[[dataset_type]][[isotype]]     <- results_dengue_serotype
    #results_dengue_serotype_neg_list[[dataset_type]][[isotype]] <- results_dengue_serotype_neg
  }
}





# save multinomial
dir.create("Results/Multinomial_Classification", showWarnings = FALSE, recursive = TRUE)
base_dir_multinom <- "Results/Multinomial_Classification"

for (dataset_type in c("ratio", "cross_sectional")) {
  dir.create(file.path(base_dir_multinom, dataset_type), recursive = TRUE, showWarnings = FALSE)
}

# ---- Bundle the two result lists together, keyed by question name ----
all_results_multinom <- list(
  dengue_serotype     = results_dengue_serotype_list
 # dengue_serotype_neg = results_dengue_serotype_neg_list
)

# ---- Save full result objects (.rds) + collect comparison tables ----
comparison_rows_multinom <- list()

for (question in names(all_results_multinom)) {
  for (dataset_type in names(all_results_multinom[[question]])) {
    for (isotype in names(all_results_multinom[[question]][[dataset_type]])) {

      result_obj <- all_results_multinom[[question]][[dataset_type]][[isotype]]
      if (is.null(result_obj)) next

      label <- paste(question, dataset_type, isotype, sep = " / ")
      cat("\n==== ", label, " ====\n", sep = "")

      out_path <- file.path(base_dir_multinom, dataset_type,
                            paste0(question, "_", isotype, ".rds"))
      t_save <- Sys.time()
      saveRDS(result_obj, out_path)
      cat(sprintf("  saveRDS %.1fs (%.0f MB)\n",
                  as.numeric(difftime(Sys.time(), t_save, units = "secs")),
                  file.size(out_path) / 1e6))

      comp <- result_obj$comparison
      comp$question     <- question
      comp$dataset_type <- dataset_type
      comp$isotype      <- isotype
      comparison_rows_multinom[[paste(question, dataset_type, isotype, sep = "_")]] <- comp
    }
  }
}

# ---- Combine all comparison tables into one long df ----
combined_comparisons_multinom <- dplyr::bind_rows(comparison_rows_multinom) %>%
  dplyr::select(
    dataset_type, question, isotype, Model,
    dplyr::any_of(c("AUC_Macro_pooled", "AUC_Micro_pooled",
                    "AUC_Macro_mean", "AUC_Macro_sd",
                    "AUC_Micro_mean", "AUC_Micro_sd",
                    "AUPRC_Macro_mean", "Brier_mean", "StratBrier_mean",
                    "n_predictions", "coverage")),
    dplyr::everything()
  )

for (ty in c("ratio", "cross_sectional")) {
  write.csv(
    combined_comparisons_multinom %>% dplyr::filter(dataset_type == ty),
    file.path(base_dir_multinom, ty, "all_comparisons.csv"),
    row.names = FALSE
  )
}

# ---- Summary across the four models, per question/isotype ----
auc_summary_multinomial <- combined_comparisons_multinom %>%
  dplyr::group_by(dataset_type, question, isotype) %>%
  dplyr::summarise(
    mean_ROC_macro        = mean(AUC_Macro_pooled, na.rm = TRUE),
    mean_ROC_micro  = mean(AUC_Micro_pooled, na.rm = TRUE),
    min_ROC         = min(AUC_Micro_pooled,  na.rm = TRUE),
    max_ROC         = max(AUC_Micro_pooled,  na.rm = TRUE),
    range_ROC       = max(AUC_Micro_pooled, na.rm = TRUE) -
                      min(AUC_Micro_pooled, na.rm = TRUE),
    sd_ROC          = stats::sd(AUC_Micro_pooled, na.rm = TRUE),  # across models
    best_model      = Model[which.max(AUC_Micro_pooled)],
    worst_model     = Model[which.min(AUC_Micro_pooled)],
    # split-to-split variability of macro AUC, averaged over the four models
    mean_ROC_by_rep = mean(AUC_Macro_mean, na.rm = TRUE),
    mean_sd_by_rep  = mean(AUC_Macro_sd,   na.rm = TRUE),
    min_coverage    = min(coverage, na.rm = TRUE),
    n_models        = dplyr::n(),
    .groups = "drop"
  ) %>%
  dplyr::arrange(dataset_type, question, isotype)

write.csv(auc_summary_multinomial,
          file.path(base_dir_multinom, "multinomial_auc_summary_across_models.csv"),
          row.names = FALSE)
View(auc_summary_multinomial)




# read data and average AUCs across model 
all_compairsons_ratio <- read.csv(here("Results/Binary_Classification/ratio/all_comparisons.csv"))

auc_avg_ratio <- all_compairsons_ratio %>%
  group_by(dataset_type, question, isotype) %>%
  summarise(
    mean_ROC = mean(ROC, na.rm = TRUE),
    n_models = dplyr::n(),
    .groups  = "drop"
  ) %>%
  arrange(question, isotype)

auc_avg_ratio


all_compairsons_cross_sectional <- read.csv(here("Results/Binary_Classification/cross_sectional/all_comparisons.csv"))

auc_avg_cross_sectional <- all_compairsons_cross_sectional %>%
  group_by(dataset_type, question, isotype) %>%
  summarise(
    mean_ROC = mean(ROC, na.rm = TRUE),
    n_models = dplyr::n(),
    .groups  = "drop"
  ) %>%
  arrange(question, isotype)

auc_avg_cross_sectional
