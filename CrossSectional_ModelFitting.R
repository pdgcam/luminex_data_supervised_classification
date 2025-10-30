#--- source functions
source(here('Functions.R'))

# extract last infection (cross-sectional) data
last_infection_data <- preprocessed_data$last_infection
# Drop rows with any NA
last_infection_data <- na.omit(last_infection_data)
table(last_infection_data$Target)

# ---- Select binomial targets / classification questions  ----
cross_sectional_data_with_binomial_targets <- select_targets(
  preprocessed_data = last_infection_data,
  targets = c("flavi", "dengue"),
  drop_original_target = TRUE
)

# Drop Nas
cross_sectional_data_with_binomial_targets$flavi <- na.omit(cross_sectional_data_with_binomial_targets$flavi)
cross_sectional_data_with_binomial_targets$dengue <- na.omit(cross_sectional_data_with_binomial_targets$dengue)

# View distribution of targets
table(cross_sectional_data_with_binomial_targets$flavi$flavi)
table(cross_sectional_data_with_binomial_targets$dengue$dengue)

# Fit binomial models
cross_sectional_binomial_modeling_results <- train_multiple_targets(
  data_list = cross_sectional_data_with_binomial_targets,
  variables = NULL,  
  k_fold = 5,
  metrics = c("AUROC", "AUPRC", "Brier", "StratBrier")
)


# ---- Select multinomial targets  ----
cross_sectional_data_with_multinomial_targets <- select_targets(
  preprocessed_data = last_infection_data,
  targets = c("dengue_serotype", "dengue_serotype_neg"),
  drop_original_target = TRUE
)

# Drop Nas
cross_sectional_data_with_multinomial_targets$dengue_serotype <- 
        na.omit(cross_sectional_data_with_multinomial_targets$dengue_serotype)
cross_sectional_data_with_multinomial_targets$dengue_serotype_neg <- 
  na.omit(cross_sectional_data_with_multinomial_targets$dengue_serotype_neg)


cross_sectional_multinomial_modeling_results <- train_multinomial_models(
  data =cross_sectional_data_with_multinomial_targets$dengue_serotype,
  target =  "dengue_serotype",
  variables = NULL, 
  k_fold = "LOOCV", 
  metrics = c("AUROC", "AUPRC", "Brier", "StratBrier"))


cross_sectional_multinomial_modeling_dengue_serotype_neg <- train_multinomial_models(
  data =cross_sectional_data_with_multinomial_targets$dengue_serotype_neg,
  target =  "dengue_serotype_neg",
  variables = NULL, 
  k_fold = "LOOCV", 
  metrics = c("AUROC", "AUPRC", "Brier", "StratBrier"))


# ---- Run univariate analysis ---- 
# Get dengue and zika variables (or variable of choice) 
serotype_vars <- grep("DENV|ZIKV", names(last_infection_data), value = TRUE)

cross_sectional_univariate_results <- train_multiple_targets_univariate(
  data_list = cross_sectional_data_with_binomial_targets,
  variables = serotype_vars,  
  k_fold = 5,
  metrics = c("AUROC", "AUPRC", "Brier", "StratBrier"),
  univariate = TRUE
)

# ---- Save all results ----
saveRDS(cross_sectional_binomial_modeling_results, here("Desktop/supervised_learning_flavi/Model_Results/cross_sectional_binomial_modeling_results.rds"))
saveRDS(cross_sectional_multinomial_modeling_dengue_serotype, here("Desktop/supervised_learning_flavi/Model_Results/cross_sectional_multinomial_modeling_dengue_serotype.rds"))
saveRDS(cross_sectional_multinomial_modeling_dengue_serotype_neg, here("Desktop/supervised_learning_flavi/Model_Results/cross_sectional_multinomial_modeling_dengue_serotype_neg.rds"))
write_csv(cross_sectional_univariate_results$combined_comparison, 
          here("Desktop/supervised_learning_flavi/Model_Results/cross_sectional_univariate_results.csv"))
