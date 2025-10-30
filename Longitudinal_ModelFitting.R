# ----- Script to fit to Binary and Multinational models to longitudinal (post /pre ratio) lumindex data ----- #

# ---- Import librarys 
library(tidyverse)
library(nnet) 
library(randomForest)
library(pROC)
library(MLmetrics)
library(DescTools)
library(glmnet)
library(corrplot)
library(MASS) 
library(dplyr)
library(purrr)
library(tidyr)
library(ggplot2)
library(scales)
library(lme4)
library(ggrepel)
library(ranger)
library(pROC)
library(C50)
library(naivebayes)
library(neuralnet)
library(caret)
library(e1071)
library(PRROC)
library(here)

# --- source functions ----
source(here('/Users/ap2488/Desktop/supervised_learning_flavi/new_code_scripts/Functions.R'))

# ---- Import prepossessed datasets ---- 
preprocessed_data <- readRDS("/Users/ap2488/Desktop/supervised_learning_flavi/preprocessed_data/luminex_processed_data.rds")

# If using already preprocessed data, skip this 
# Extract ratio data (longitudinal analysis)
ratio_data <- preprocessed_data$ratio
# Drop rows with any NA
ratio_data <- na.omit(ratio_data)
table(ratio_data$Target)


# ---- Run on multiple (binomial) targets  ----
# Define targets / classification question
data_with_binomial_targets <- select_targets(
  preprocessed_data = ratio_data,
  targets = c("flavi", "dengue"),
  drop_original_target = TRUE
)

# Drop Nas
data_with_binomial_targets$flavi <- na.omit(data_with_binomial_targets$flavi)
data_with_binomial_targets$dengue <- na.omit(data_with_binomial_targets$dengue)
data_with_binomial_targets$zika <- na.omit(data_with_binomial_targets$zika)

# View distribution of targets
table(data_with_binomial_targets$flavi$flavi)
table(data_with_binomial_targets$dengue$dengue)


# Fit binomial models
binomial_modeling_results <- train_multiple_targets(
  data_list = data_with_binomial_targets,
  variables = NULL,  # Uses all columns except target
  k_fold = 5,
  metrics = c("AUROC", "AUPRC", "Brier")
)



# ---- Run on multinomial targets  ----
#Define targets 
data_with_multinomial_targets <- select_targets(
  preprocessed_data = ratio_data,
  targets = c("dengue_serotype", "dengue_serotype_neg"),
  drop_original_target = TRUE
)

# Drop Nas
data_with_multinomial_targets$dengue_serotype <- na.omit(data_with_multinomial_targets$dengue_serotype)
data_with_multinomial_targets$dengue_serotype_neg <- na.omit(data_with_multinomial_targets$dengue_serotype_neg)

# View distribution of targets
table(data_with_multinomial_targets$dengue_serotype$dengue_serotype)
table(data_with_multinomial_targets$dengue_serotype_neg$dengue_serotype_neg)


# Given infection, classify stereotype 
dengue_serotype_results <- train_multinomial_models(
  data = data_with_multinomial_targets$dengue_serotype,
  target = "dengue_serotype",
  variables = NULL,
  k_fold = "LOOCV", 
  metrics = c("AUROC", "AUPRC", "Brier", "StratBrier"))


# Classify DENV1 vs DENV2 vs DENV3 vs DEN4 vs Neg 
dengue_serotype_neg_results <- train_multinomial_models(
  data = data_with_multinomial_targets$dengue_serotype_neg,
  target =  "dengue_serotype_neg",
  variables = NULL, 
  k_fold = "LOOCV", 
  metrics = c("AUROC", "AUPRC", "Brier", "StratBrier"))


# ----  Run univariate analysis (one variable at a time)
# Get dengue and zika serotype variables (or single variable of interest)
serotype_vars <- grep("DENV|ZIKV", names(ratio_data), value = TRUE)

univariate_results <- train_multiple_targets_univariate(
  data_list = data_with_binomial_targets,
  variables = serotype_vars,  
  k_fold = 5,
  metrics = c("AUROC", "AUPRC", "Brier"),
  univariate  = TRUE
)


# --- Save results --- 
saveRDS(binomial_modeling_results,  here("longitudinal_binomial_modeling_results.rds"))
saveRDS(dengue_serotype_results, here("dengue_serotype_results.rds"))
saveRDS(dengue_serotype_neg_results, here("dengue_serotype_neg_results.rds"))
write_csv(univariate_results$combined_comparison, here("flavi_univariate_results.csv"))
