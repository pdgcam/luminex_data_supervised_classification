library(caret)
library(here)
library(corrplot)
packageVersion("xgboost")
packageVersion("caret")
library(dplyr)



# install remotes if you don't have it
install.packages("remotes")

# downgrade xgboost to a known-compatible version
remotes::install_version("xgboost", version = "1.7.8.1", repos = "https://cran.r-project.org")

# --- Source functions
source(here('Models/train_binary_models.R'))
source(here('Models/train_multinomial_models.R'))
source(here('Functions.R'))

packageVersion("xgboost")

isotypes <- c("IgG", "IgA", "IgM", "avidity")

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


all_antigens <- c(flavi_antigens, alpha_antigens)

# ---- Import prepossessed datasets ---- 
# -- Ratio (Longitudnal)
logged_ratio_df_IgG <- readRDS("Data/by_isotype/logged_ratio_df_IgG.rds")
logged_ratio_df_IgA <- readRDS("Data/by_isotype/logged_ratio_df_IgA.rds")
logged_ratio_df_IgM <- readRDS("Data/by_isotype/logged_ratio_df_IgM.rds")
logged_ratio_df_avidity <- readRDS("Data/by_isotype/logged_ratio_df_avidity.rds")


# -- Cross-sectional 
cross_sectional_df_IgG <- readRDS("Data/by_isotype/cross_sectional_df_IgG.rds")
cross_sectional_df_IgA <- readRDS("Data/by_isotype/cross_sectional_df_IgA.rds")
cross_sectional_df_IgM <- readRDS("Data/by_isotype/cross_sectional_df_IgM.rds")
cross_sectional_df_avidity <- readRDS("Data/by_isotype/cross_sectional_df_avidity.rds")

# cross-section have 3 more samples (these didnt have all timepoints - so had to be dropped from ratio)

#logged cross-sectional
logged_cross_sectional_df_IgG <- cross_sectional_df_IgG %>%
  mutate(across(all_of(all_antigens), log10))

logged_cross_sectional_df_IgA <-cross_sectional_df_IgA %>%
  mutate(across(all_of(all_antigens), log10))

logged_cross_sectional_df_IgM <- cross_sectional_df_IgM %>%
  mutate(across(all_of(all_antigens), log10))

logged_cross_sectional_df_avidity <- cross_sectional_df_avidity %>%
  mutate(across(all_of(all_antigens), log10))

# dataset types we'll loop over, and a helper to fetch the right df
dataset_types <- c("ratio", "cross_sectional")

get_source_df <- function(dataset_type, isotype) {
  prefix <- if (dataset_type == "ratio") "logged_ratio_df_" else "logged_cross_sectional_df_"
  get(paste0(prefix, isotype))
}


# --- Define output target 
data_with_binomial_targets_flavi_list <- list(ratio = list(), cross_sectional = list())
data_with_binomial_targets_chik_list  <- list(ratio = list(), cross_sectional = list())
data_with_multinomial_targets_list    <- list(ratio = list(), cross_sectional = list())

for (dataset_type in dataset_types) {
  for (isotype in isotypes) {

    df_logged <- get_source_df(dataset_type, isotype)

    # subset flavi and alpha antigens
    df_logged_flavi <- df_logged %>%
      dplyr::select(target, id_patient, all_of(flavi_antigens))

    df_logged_alpha <- df_logged %>%
      dplyr::select(target, id_patient, all_of(alpha_antigens))

    # binomial - flavi / dengue
    binom_flavi <- select_targets(
      preprocessed_data = df_logged_flavi,
      targets = c("flavi", "dengue"),
      drop_original_target = TRUE,
      min_samples = 2
    )

    # binomial - chik / dengue_chik
    binom_chik <- select_targets(
      preprocessed_data = df_logged_alpha,
      targets = c("chik", "dengue_chik"),
      drop_original_target = TRUE,
      min_samples = 2
    )

    # multinomial - dengue serotype
    multinom <- select_targets(
      preprocessed_data = df_logged_flavi,
      targets = c("dengue_serotype", "dengue_serotype_neg"),
      drop_original_target = TRUE,
      min_samples = 2
    )

    data_with_binomial_targets_flavi_list[[dataset_type]][[isotype]] <- binom_flavi
    data_with_binomial_targets_chik_list[[dataset_type]][[isotype]]  <- binom_chik
    data_with_multinomial_targets_list[[dataset_type]][[isotype]]    <- multinom
  }
}

# View distribution of targets - binomial (only looking at IgG)
table(data_with_binomial_targets_flavi_list$ratio$IgG$data$dengue$dengue_target)
table(data_with_binomial_targets_flavi_list$ratio$IgG$data$flavi$flavi_target)
table(data_with_binomial_targets_chik_list$ratio$IgG$data$dengue_chik$dengue_chik_target)
table(data_with_binomial_targets_chik_list$ratio$IgG$data$chik$chik_target)

# View distribution of targets - binomial (only looking at IgG)
table(data_with_binomial_targets_flavi_list$cross_sectional$IgG$data$dengue$dengue_target)
table(data_with_binomial_targets_flavi_list$cross_sectional$IgG$data$flavi$flavi_target)
table(data_with_binomial_targets_chik_list$cross_sectional$IgG$data$dengue_chik$dengue_chik_target)
table(data_with_binomial_targets_chik_list$cross_sectional$IgG$data$chik$chik_target)

# Multinomial - also only looking at IgG
table(data_with_multinomial_targets_list$ratio$IgG$data$dengue_serotype$dengue_serotype_target)
table(data_with_multinomial_targets_list$ratio$IgG$data$dengue_serotype_neg$dengue_serotype_neg_target)


results_flavi_vs_not_list   <- list(ratio = list(), cross_sectional = list())
results_dengue_vs_not_list  <- list(ratio = list(), cross_sectional = list())
results_chik_vs_not_list    <- list(ratio = list(), cross_sectional = list())
results_dengue_vs_chik_list <- list(ratio = list(), cross_sectional = list())


# Binomial Classification 
for (dataset_type in dataset_types) {
  for (isotype in isotypes) {

    binom_flavi <- data_with_binomial_targets_flavi_list[[dataset_type]][[isotype]]
    binom_chik  <- data_with_binomial_targets_chik_list[[dataset_type]][[isotype]]

    # Flavi vs not
    results_flavi_vs_not <- train_binary_models(
      data = binom_flavi$data$flavi,
      target = "flavi_target",
      positive_class = binom_flavi$positive_class_map$flavi,
      variables = NULL,
      metrics = c("ROC", "AUPRC", "Brier")
    )

    # Dengue vs not
    results_dengue_vs_not <- train_binary_models(
      data = binom_flavi$data$dengue,
      target = "dengue_target",
      positive_class = binom_flavi$positive_class_map$dengue,
      variables = NULL,
      metrics = c("ROC", "AUPRC", "Brier")
    )

    # Chik vs not
    results_chik_vs_not <- train_binary_models(
      data = binom_chik$data$chik,
      target = "chik_target",
      positive_class = binom_chik$positive_class_map$chik,
      variables = NULL,
      metrics = c("ROC", "AUPRC", "Brier")
    )

    # Dengue vs chik
    results_dengue_vs_chik <- train_binary_models(
      data = binom_chik$data$dengue_chik,
      target = "dengue_chik_target",
      positive_class = binom_chik$positive_class_map$dengue_chik,
      variables = NULL,
      metrics = c("ROC", "AUPRC", "Brier")
    )

    results_flavi_vs_not_list[[dataset_type]][[isotype]]   <- results_flavi_vs_not
    results_dengue_vs_not_list[[dataset_type]][[isotype]]  <- results_dengue_vs_not
    results_chik_vs_not_list[[dataset_type]][[isotype]]    <- results_chik_vs_not
    results_dengue_vs_chik_list[[dataset_type]][[isotype]] <- results_dengue_vs_chik
  }
}



# Multinomial Classification 
results_dengue_serotype_list     <- list(ratio = list(), cross_sectional = list())
results_dengue_serotype_neg_list <- list(ratio = list(), cross_sectional = list())

for (dataset_type in dataset_types) {
  for (isotype in isotypes) {

    multinom <- data_with_multinomial_targets_list[[dataset_type]][[isotype]]

    # Dengue serotype (infected-only: DENV1-4)
    results_dengue_serotype <- train_multinomial_models(
      data = multinom$data$dengue_serotype,
      target = "dengue_serotype_target",
      variables = NULL,
      metrics = c("ROC", "AUPRC", "Brier", "StratBrier")
    )

    # Dengue serotype + negative class
    results_dengue_serotype_neg <- train_multinomial_models(
      data = multinom$data$dengue_serotype_neg,
      target = "dengue_serotype_neg_target",
      variables = NULL,
      metrics = c("ROC", "AUPRC", "Brier", "StratBrier")
    )

    results_dengue_serotype_list[[dataset_type]][[isotype]]     <- results_dengue_serotype
    results_dengue_serotype_neg_list[[dataset_type]][[isotype]] <- results_dengue_serotype_neg
  }
}


# ----Flavi vs Not 

# Get best performing tuned model for each isotype - ratios
results_flavi_vs_not_list$ratio$IgG$comparison
results_flavi_vs_not_list$ratio$IgA$comparison
results_flavi_vs_not_list$ratio$IgM$comparison
results_flavi_vs_not_list$ratio$avidity$comparison

# Cross-sectional
results_flavi_vs_not_list$cross_sectional$IgG$comparison
results_flavi_vs_not_list$cross_sectional$IgA$comparison
results_flavi_vs_not_list$cross_sectional$IgM$comparison
results_flavi_vs_not_list$cross_sectional$avidity$comparison


# Save all results (binomial)
dir.create("Results/Binary_Classification", showWarnings = FALSE, recursive = TRUE)
base_dir <- "Results/Binary_Classification"

for (dataset_type in c("ratio", "cross_sectional")) {
  dir.create(file.path(base_dir, dataset_type), recursive = TRUE, showWarnings = FALSE)
}

# ---- Bundle the four result lists together, keyed by question name ----
all_results <- list(
  flavi_vs_not   = results_flavi_vs_not_list,
  dengue_vs_not  = results_dengue_vs_not_list,
  chik_vs_not    = results_chik_vs_not_list,
  dengue_vs_chik = results_dengue_vs_chik_list
)

# ---- Save full result objects (.rds) + collect comparison tables ----
comparison_rows <- list()

for (question in names(all_results)) {
  for (dataset_type in names(all_results[[question]])) {
    for (isotype in names(all_results[[question]][[dataset_type]])) {

      result_obj <- all_results[[question]][[dataset_type]][[isotype]]

      # save the full result object (models, predictions, comparison, aucs, etc.)
      out_path <- file.path(base_dir, dataset_type, paste0(question, "_", isotype, ".rds"))
      saveRDS(result_obj, out_path)

      # collect comparison table for the combined CSV
      comp <- result_obj$comparison
      comp$question     <- question
      comp$dataset_type <- dataset_type
      comp$isotype      <- isotype

      comparison_rows[[paste(question, dataset_type, isotype, sep = "_")]] <- comp
    }
  }
}

# ---- Combine all comparison tables into one long df, split by dataset_type ----
combined_comparisons <- dplyr::bind_rows(comparison_rows) %>%
  dplyr::select(dataset_type, question, isotype, Model, dplyr::everything())

for (dataset_type in c("ratio", "cross_sectional")) {
  write.csv(
    combined_comparisons %>% dplyr::filter(dataset_type == !!dataset_type),
    file.path(base_dir, dataset_type, "all_comparisons.csv"),
    row.names = FALSE
  )
}



# save multinomial
dir.create("Results/Multinomial_Classification", showWarnings = FALSE, recursive = TRUE)
base_dir_multinom <- "Results/Multinomial_Classification"

for (dataset_type in c("ratio", "cross_sectional")) {
  dir.create(file.path(base_dir_multinom, dataset_type), recursive = TRUE, showWarnings = FALSE)
}

# ---- Bundle the two result lists together, keyed by question name ----
all_results_multinom <- list(
  dengue_serotype     = results_dengue_serotype_list,
  dengue_serotype_neg = results_dengue_serotype_neg_list
)

# ---- Save full result objects (.rds) + collect comparison tables ----
comparison_rows_multinom <- list()

for (question in names(all_results_multinom)) {
  for (dataset_type in names(all_results_multinom[[question]])) {
    for (isotype in names(all_results_multinom[[question]][[dataset_type]])) {

      result_obj <- all_results_multinom[[question]][[dataset_type]][[isotype]]

      # save the full result object (models, predictions, comparison/oof_metrics, etc.)
      out_path <- file.path(base_dir_multinom, dataset_type, paste0(question, "_", isotype, ".rds"))
      saveRDS(result_obj, out_path)

      # collect comparison table for the combined CSV
      comp <- result_obj$comparison
      comp$question     <- question
      comp$dataset_type <- dataset_type
      comp$isotype      <- isotype

      comparison_rows_multinom[[paste(question, dataset_type, isotype, sep = "_")]] <- comp
    }
  }
}

# ---- Combine all comparison tables into one long df, split by dataset_type ----
combined_comparisons_multinom <- dplyr::bind_rows(comparison_rows_multinom) %>%
  dplyr::select(dataset_type, question, isotype, Model, dplyr::everything())

for (dataset_type in c("ratio", "cross_sectional")) {
  write.csv(
    combined_comparisons_multinom %>% dplyr::filter(dataset_type == !!dataset_type),
    file.path(base_dir_multinom, dataset_type, "all_comparisons.csv"),
    row.names = FALSE
  )
}



# read data and average AUCs across model 

all_compairsons_ratio <- read.csv(here("Results/Binary_Classification/ratio/all_comparisons.csv"))
head(all_compairsons_ratio)
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


# multinomial 

all_compairsons_ratio_multinomial <- read.csv(here("Results/Multinomial_Classification/ratio/all_comparisons.csv"))

auc_avg_ratio_multinomial <- all_compairsons_ratio_multinomial %>%
  group_by(dataset_type, question, isotype) %>%
  summarise(
    mean_ROC = mean(AUC_Micro, na.rm = TRUE),
    n_models = dplyr::n(),
    .groups  = "drop"
  ) %>%
  arrange(question, isotype)

auc_avg_ratio_multinomial

all_compairsons_cross_sectional_multinomial <- read.csv(here("Results/Multinomial_Classification/cross_sectional/all_comparisons.csv"))

auc_avg_cross_sectional_multinomial <- all_compairsons_cross_sectional_multinomial %>%
  group_by(dataset_type, question, isotype) %>%
  summarise(
    mean_ROC = mean(AUC_Micro, na.rm = TRUE),
    n_models = dplyr::n(),
    .groups  = "drop"
  ) %>%
  arrange(question, isotype)

auc_avg_cross_sectional_multinomial



# tag each summary with its classification type (dataset_type already
# distinguishes ratio vs cross_sectional, so we just need this)
auc_avg_ratio$comparison_type                  <- "Binary"
auc_avg_cross_sectional$comparison_type        <- "Binary"
auc_avg_ratio_multinomial$comparison_type      <- "Multinomial"
auc_avg_cross_sectional_multinomial$comparison_type <- "Multinomial"

all_auc <- bind_rows(
  auc_avg_ratio,
  auc_avg_cross_sectional,
  auc_avg_ratio_multinomial,
  auc_avg_cross_sectional_multinomial
)

plot_auc_bar <- function(df, chance_line = 0.5) {
  ggplot(df, aes(x = question, y = mean_ROC, fill = isotype)) +
    geom_col(
      position  = position_dodge(width = 0.8),
      width     = 0.8,
      linewidth = 0.3
    ) +
    geom_text(
      aes(label = sprintf("%.2f", mean_ROC)),
      position = position_dodge(width = 0.8),
      vjust    = -0.5,
      size     = 5,
      colour   = "black"
    ) +
    scale_y_continuous(limits = c(0, 1), expand = expansion(mult = c(0, 0.1))) +
    scale_x_discrete(
      labels = function(x) tools::toTitleCase(gsub("_", " ", x))
    ) +
    scale_fill_manual(
      values = c(
        "IgG"     = "#8F94B7",
        "IgA"     = "#227DAA",
        "IgM"     = "#2D3047",
        "avidity" = "#71bbd8"
      )
    ) +
    facet_grid(dataset_type ~ comparison_type,  scales = "fixed") +
    labs(
      x    = NULL,
      y    = "Mean AUC (averaged across models)",
      fill = "Isotype"
    ) +
    theme_minimal() +
    theme(
      axis.text.x      = element_text(hjust = 0.5, size = 20, colour = "black"),
      axis.text.y      = element_text(size = 20),
      axis.title.y     = element_text(size = 20), 
      legend.text      = element_text(size = 20),
      panel.grid.minor = element_blank(),
      strip.text       = element_text(size = 20),
      strip.background = element_rect(fill = "#ffffff", colour = NA),
      legend.position  = "bottom", aspect.ratio = 0.5
    )
}




plot_auc_bar(auc_avg_ratio)
plot_auc_bar(auc_avg_cross_sectional)
plot_auc_bar(auc_avg_ratio_multinomial)
plot_auc_bar(auc_avg_cross_sectional_multinomial)
