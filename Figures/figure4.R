library(caret)
library(here)
library(corrplot)

# --- Source functions
source(here('Models/train_binary_models.R'))
source(here('Models/train_multinomial_models.R'))
source(here('Functions.R'))



# ---- Import prepossessed datasets ---- 
ratio_df_logged <- readRDS('Results/logged_ratio_df.rds')
table(ratio_df_logged$target)


flavi_antigens <- c("DENV1_DIII","DENV1_NS1","DENV1_VLP","SHERPADES_DENV1_DIII",
"DENV2_DIII","DENV2_NS1","DENV2_VLP","SHERPADES_DENV2_DIII",
"DENV3_DIII","DENV3_NS1","DENV3_VLP", "SHERPADES_DENV3_DIII",
"DENV4_DIII","DENV4_NS1","DENV4_VLP", "SHERPADES_DENV4_DIII",
"JEV_E", "JEV_NS1", "SHERPADES_JEV_DIII",
"YFV_E", "YFV_NS1", "SHERPADES_YFV_DIII",
"WNV_DIII","WNV_NS1","SHERPADES_WNV_DIII",
"ZIKV_NS1","ZIKV_VLP","ZIKVAS_DIII","ZIKVSU_NS1","SHERPADES_ZIKV_DIII")

alpha_antigens <- c("CHIKV_E2", "CHIKV_NSP123", "CHIKV_VLP", "SHERPADES_CHIKV_E2", 
                    "MAYV_E2" , "SHERPADES_MAYV_E2",
                    "ONNV_E2", "ONNV_VLP",
                    "RR" , "SHERPADES_RR")


# subset flavi and alpha antigens 
ratio_df_logged_flavi <- ratio_df_logged %>%
  dplyr::select(target, id_patient, all_of(flavi_antigens))

ratio_df_logged_alpha <- ratio_df_logged %>%
  dplyr::select(target, id_patient, all_of(alpha_antigens))

# Define targets / classification question
data_with_binomial_targets <- select_targets(
  preprocessed_data = ratio_df_logged_flavi,
  targets = c("flavi", "dengue"),
  drop_original_target = TRUE,
  min_samples = 2
)


data_with_binomial_targets_chik <- select_targets(
  preprocessed_data = ratio_df_logged_alpha,
  targets = c("chik", "dengue_chik"),
  drop_original_target = TRUE,
  min_samples = 2
)


data_with_multinomial_targets <- select_targets(
  preprocessed_data = ratio_df_logged_flavi,
  targets = c("dengue_serotype", "dengue_serotype_neg"),
  drop_original_target = TRUE, 
  min_samples = 2
)

# View distribution of targets
table(data_with_binomial_targets$data$flavi$flavi)
table(data_with_binomial_targets$data$flavi$flavi)
table(data_with_binomial_targets$data$dengue$dengue)
table((data_with_binomial_targets_chik$data$dengue_chik$dengue_chik))
table((data_with_binomial_targets_chik$data$chik$chik))

data_with_multinomial_targets$dengue_serotype <- na.omit(data_with_multinomial_targets$dengue_serotype)
data_with_multinomial_targets$dengue_serotype_neg <- na.omit(data_with_multinomial_targets$dengue_serotype_neg)
table(data_with_multinomial_targets$data$dengue_serotype$dengue_serotype)
table(data_with_multinomial_targets$data$dengue_serotype_neg$dengue_serotype_neg)



# ---Binomial Classification
# Flavi vs not 
results_flavi_vs_not <- train_binary_models(
  data = data_with_binomial_targets$data$flavi,
  target = "flavi",
  positive_class = data_with_binomial_targets$positive_class_map$flavi,
  variables = NULL,  
  metrics = c("ROC", "AUPRC", "Brier"))

results_flavi_vs_not$comparison


# Dengue vs not 
results_dengue_vs_not <- train_binary_models(
  data = data_with_binomial_targets$data$dengue,
  target = "dengue",
  positive_class = data_with_binomial_targets$positive_class_map$dengue,
  variables = NULL, 
  metrics = c("ROC", "AUPRC", "Brier"))

results_dengue_vs_not$comparison

# Chik vs not 
results_chik_vs_not <- train_binary_models(
  data = data_with_binomial_targets_chik$data$chik,
  target = "chik",
  positive_class = data_with_binomial_targets_chik$positive_class_map$chik,
  variables = NULL, 
  metrics = c("ROC", "AUPRC", "Brier"),
)

results_chik_vs_not$comparison

# Dengue vs chik 
results_dengue_vs_chik <- train_binary_models(
  data = data_with_binomial_targets_chik$data$dengue_chik,
  target = "dengue_chik",
  positive_class = data_with_binomial_targets_chik$positive_class_map$dengue_chik,
  variables = NULL, 
  metrics = c("ROC", "AUPRC", "Brier"),
)


# results
results_flavi_vs_not$comparison
results_dengue_vs_not$comparison
results_dengue_vs_chik$comparison
results_chik_vs_not$comparison



# --- Multinomial Results 

# Given infection, classify stereotype + serotype and neg 
results_dengue_serotype <- train_multinomial_models(
  data = data_with_multinomial_targets$data$dengue_serotype,
  target = "dengue_serotype",
  variables = NULL,
  metrics = c("ROC", "AUPRC", "Brier", "StratBrier"))

results_dengue_serotype_neg <-  train_multinomial_models(
  data = data_with_multinomial_targets$data$dengue_serotype_neg,
  target = "dengue_serotype_neg",
  variables = NULL,
  metrics = c("ROC", "AUPRC", "Brier", "StratBrier"))


results_dengue_serotype$comparison
results_dengue_serotype_neg$comparison

# Save results 
saveRDS(results_flavi_vs_not, 'Results/NewClassificationResults/results_flavi_vs_not.rds')
saveRDS(results_dengue_vs_not, 'Results/NewClassificationResults/results_dengue_vs_not.rds')
saveRDS(results_dengue_vs_chik, 'Results/NewClassificationResults/results_dengue_vs_chik.rds')
saveRDS(results_chik_vs_not, 'Results/NewClassificationResults/results_chik_vs_not.rds')
saveRDS(results_dengue_serotype, 'Results/NewClassificationResults/results_dengue_serotype.rds')
saveRDS(results_dengue_serotype_neg, 'Results/NewClassificationResults/results_dengue_serotype_neg.rds')

# Micro-average "takes imbalance into account" in the sense that the resulting performance is based on the proportion of every class
# i.e.the performance of a large class has more impact on the result than of a small class.
# Macro-average "doesn't take imbalance into account" in the sense that the resulting performance 
# is a simple average over the classes, so every class is given equal weight independently from their proportion.







# ---- Best from binomial results ----
best_binomial <- bind_rows(
  results_flavi_vs_not$combined_comparison,
  results_dengue_vs_not$combined_comparison,
  results_dengue_vs_chik$combined_comparison
) %>%
  group_by(target) %>%
  slice_max(AUROC, n = 1) %>%
  ungroup() %>%
  rename(AUC = AUROC) %>%
  mutate(type = "binary")


# ---- Best from multinomial results ----
best_multinomial <- results_dengue_serotype$combined_comparison %>%
  group_by(target) %>%
  slice_max(AUC_Micro, n = 1) %>%
  ungroup() %>%
  rename(AUC = AUC_Micro) %>%
  mutate(type = "multinomial")


# ---- Combine into summary table ----
best_models <- bind_rows(
  best_binomial %>% dplyr::select(type, target, Model, AUC, AUPRC, Brier),
  best_multinomial %>% dplyr::select(type, target, Model, AUC = AUC, 
                               AUPRC = AUPRC_Micro, Brier)
) %>%
  arrange(type, target)


saveRDS(best_models, 'Results/multivariate_best_models.rds')


