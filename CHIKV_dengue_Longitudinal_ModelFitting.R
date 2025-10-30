# ---- Script for binomial analysis: CHIKV vs denuge

#--- source functions
source(here('Functions.R'))

#  ---- import preprocessed datasets ---- 
preprocessed_data <- readRDS("luminex_processed_data_with_chik.rds")
# extract ratio data (longitudinal analysis)
ratio_data <- preprocessed_data$ratio
# Drop rows with any NA
ratio_data <- na.omit(ratio_data)
table(ratio_data$Target)


# ---- Run Model on multiple target  (binomial) ----
data_with_binomial_targets_chik <- select_targets(
  preprocessed_data = ratio_data,
  targets = c("dengue_chik"),
  drop_original_target = TRUE
)


# Drop nas
data_with_binomial_targets_chik$dengue_chik <- na.omit(data_with_binomial_targets_chik$dengue_chik)

# View distribution of targets
table(data_with_binomial_targets_chik$dengue_chik$dengue_chik)

binomial_modeling_results <- train_multiple_targets(
  data_list = data_with_binomial_targets_chik,
  variables = NULL,  # Uses all columns except target
  k_fold = 5,
  metrics = c("AUROC", "AUPRC", "Brier")
)


# --- save results --- 
saveRDS(binomial_modeling_results,  here("Desktop/supervised_learning_flavi/Model_Results/chik_longitudinal_binomial_modeling_results.rds"))

