library(caret)
library(here)
packageVersion("xgboost")
packageVersion("caret")
library(dplyr)



# --- Source functions
#source(here('Models/train_binary_models.R'))
source(here('Models/train_multinomial_models.R'))
source(here('Models/bootstrap_AUC.R'))
source(here('Models/NEW_train_binary_models.R'))



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


all_antigens <- c(flavi_antigens, alpha_antigens)

# ---- Import prepossessed datasets ---- 
# -- Ratio (Longitudnal)
logged_ratio_df_IgG <- readRDS("Data/by_isotype_gm/logged_ratio_df_IgG.rds")
logged_ratio_df_IgA <- readRDS("Data/by_isotype_gm/logged_ratio_df_IgA.rds")
logged_ratio_df_IgM <- readRDS("Data/by_isotype_gm/logged_ratio_df_IgM.rds")


# -- Cross-sectional (IgG and IgM)
flavi_cross_sectional_df_IgG <- readRDS("Data/by_isotype/flavi_cross_sectional_logged_IgG.rds")
flavi_cross_sectional_df_IgM <- readRDS("Data/by_isotype/flavi_cross_sectional_logged_IgM.rds")

nrow(flavi_cross_sectional_df_IgG)
table(flavi_cross_sectional_df_IgG$target)

alpha_cross_sectional_df_IgG <- readRDS("Data/by_isotype/alpha_cross_sectional_logged_IgA.rds")
alpha_cross_sectional_df_IgM <- readRDS("Data/by_isotype/alpha_cross_sectional_logged_IgM.rds")

nrow(alpha_cross_sectional_df_IgG)
table(alpha_cross_sectional_df_IgG$target)


# dataset types we'll loop over, and a helper to fetch the right df
dataset_types <- c("ratio", "cross_sectional")

get_source_df <- function(dataset_type, isotype, panel) {

  antigens <- if (panel == "flavi") flavi_antigens else alpha_antigens

  path <- if (dataset_type == "cross_sectional") {
    sprintf("Data/by_isotype/%s_cross_sectional_logged_%s.rds", panel, isotype)
  } else {
    sprintf("Data/by_isotype/logged_ratio_df_%s.rds", isotype)
  }

  if (!file.exists(path)) stop("File not found: ", path)

  readRDS(path) %>%
    dplyr::select(target, id_patient, all_of(antigens))
}

select_targets <- function(preprocessed_data,
                           targets = c("flavi", "dengue", 
                                       "zika", "chik",
                                       "dengue_zika",  "dengue_chik",
                                       "dengue_serotype", "dengue_serotype_neg"), 
                                       drop_original_target = TRUE,  min_samples = 2) {
  
  if (!"target" %in% names(preprocessed_data)) {
    stop('Column "Target" not found in preprocessed_data.')
  }
  
  # Validate all targets
  valid_targets <- c("flavi", "dengue", 
                     "zika", "chik",
                     "dengue_zika", 
                     "dengue_serotype", "dengue_serotype_neg",
                     "dengue_chik")
  
  invalid <- setdiff(targets, valid_targets)
  if (length(invalid) > 0) {
    stop("Invalid target(s): ", paste(invalid, collapse = ", "))
  }
  
  # mappings
  mapping_list <- list(
    flavi               = c(DENV1 = 1, DENV2 = 1, DENV3 = 1, DENV4 = 1, ZIKV = 1),
    dengue              = c(DENV1 = 1, DENV2 = 1, DENV3 = 1, DENV4 = 1),
    zika                = c(ZIKV = 1),
    chik                = c(CHIKV = 1),
    dengue_zika         = c(DENV1 = 1, DENV2 = 1, DENV3 = 1, DENV4 = 1, ZIKV = 2),
    dengue_serotype     = c(DENV1 = 1, DENV2 = 2, DENV3 = 3, DENV4 = 4),
    dengue_serotype_neg = c(DENV1 = 1, DENV2 = 2, DENV3 = 3, DENV4 = 4),
    dengue_chik        = c(DENV1 = 1, DENV2 = 1, DENV3 = 1, DENV4 = 1, CHIKV = 2))
  
  # factor levels/labels
  label_specs <- list(
    flavi               = list(levels = c(0, 1),           labels = c("negative", "positive")),
    dengue              = list(levels = c(0, 1),           labels = c("negative", "positive")),
    zika                = list(levels = c(0, 1),           labels = c("negative", "positive")),
    chik                = list(levels = c(0, 1),           labels = c("negative", "positive")),
    dengue_zika         = list(levels = c(0, 1, 2),        labels = c("negative", "dengue", "zika")),
    dengue_serotype     = list(levels = c(1, 2, 3, 4),     labels = c("DENV1", "DENV2", "DENV3", "DENV4")),
    dengue_serotype_neg = list(levels = c(0, 1, 2, 3, 4),  labels = c("negative", "DENV1", "DENV2", "DENV3", "DENV4")),
    dengue_chik         = list(levels = c(1, 2),           labels = c("dengue", "CHIKV"))
  )

  positive_class_map <- list(
  flavi       = "positive",
  dengue      = "positive",
  zika        = "positive",
  chik        = "positive",
  dengue_zika = "dengue",
  dengue_chik = "dengue"
)
  
  
  # List to store results
  data_with_target <- list()
  
  # Process each target
  for (target in targets) {
    df_copy <- preprocessed_data

    target_col <- paste0(target, "_target")
    
    # Build the chosen target column
    mp <- mapping_list[[target]]
    tgt_chr <- as.character(df_copy$target)
    
   # Map positives, everything else is negative (0)
    df_copy[[target_col]] <- ifelse(
      tgt_chr %in% names(mp),
      as.numeric(dplyr::recode(tgt_chr, !!!mp, .default = NA_real_)),
      0
    )
    
    # Remove negatives for classifications that only include infected samples (eg: given infection, classify dengue vs chik)
    if (target %in% c("dengue_chik", "dengue_serotype")) {
      df_copy <- df_copy[df_copy[[target_col]] %in% c(1, 2, 3, 4), ]  # keep only mapped positives
    }
    # convert to factor with labels
    spec <- label_specs[[target]]
    df_copy[[target_col]] <-
      factor(df_copy[[target_col]], 
             levels = spec$levels, 
             labels = spec$labels)

    # Drop classes with fewer than min_samples
    class_counts <- table(df_copy[[target_col]])
    valid_classes <- names(class_counts[class_counts >= min_samples])
    
    if (length(valid_classes) < length(class_counts)) {
      dropped <- names(class_counts[class_counts < min_samples])
      cat(sprintf("Target '%s': dropping class(es) with < %d samples: %s\n",
                  target, min_samples, paste(dropped, collapse = ", ")))
      df_copy <- df_copy[df_copy[[target_col]] %in% valid_classes, ]
      df_copy[[target_col]] <- droplevels(df_copy[[target_col]])  #
    }
    
    
    # drop original Target column if requested
    if (drop_original_target) {
      if ("target" %in% names(df_copy)) {
        df_copy[["target"]] <- NULL
      }
    }
    
    data_with_target[[target]] <- df_copy
  }
  return(list(
  data = data_with_target,
  positive_class_map = positive_class_map
))
}


# for cross-sectional datasets, we dont want zika to be coded as 'neg' in the deng vs not classification 
drop_zikv_negatives <- function(res, source_df, targets = c("flavi", "dengue")) {
  for (tg in targets) {
    if (is.null(res$data[[tg]])) next

    df  <- res$data[[tg]]
    col <- paste0(tg, "_target")

    # patients whose original label was ZIKV but who got encoded as negative
    zikv_ids <- source_df$id_patient[source_df$target == "ZIKV"]
    bad <- df$id_patient %in% zikv_ids & df[[col]] == "negative"

    if (any(bad)) {
      cat(sprintf("Target '%s': dropping %d ZIKV sample(s) encoded as negative\n",
                  tg, sum(bad)))
      df <- df[!bad, ]
      df[[col]] <- droplevels(df[[col]])
      res$data[[tg]] <- df
    }
  }
  res
}

# View distribution of targets - binomial (only looking at IgG)
table(data_with_binomial_targets_flavi_list$ratio$IgG$data$dengue$dengue_target)
table(data_with_binomial_targets_flavi_list$ratio$IgG$data$flavi$flavi_target)
table(data_with_binomial_targets_chik_list$ratio$IgG$data$dengue_chik$dengue_chik_target)
table(data_with_binomial_targets_chik_list$ratio$IgG$data$chik$chik_target)


# View distribution of targets - binomial (only looking at IgG)
# NOTE - CHIKV VS DENV is not correctly run in the select target function - NEED TO FIX!!
table(data_with_binomial_targets_flavi_list$cross_sectional$IgG$data$dengue$dengue_target)
table(data_with_binomial_targets_flavi_list$cross_sectional$IgG$data$flavi$flavi_target)
table(data_with_binomial_targets_chik_list$cross_sectional$IgG$data$chik$chik_target)

# Multinomial - also only looking at IgG
table(data_with_multinomial_targets_list$ratio$IgG$data$dengue_serotype$dengue_serotype_target)
table(data_with_multinomial_targets_list$ratio$IgG$data$dengue_serotype_neg$dengue_serotype_neg_target)

table(data_with_multinomial_targets_list$cross_sectional$IgG$data$dengue_serotype$dengue_serotype_target)
table(data_with_multinomial_targets_list$cross_sectional$IgG$data$dengue_serotype_neg$dengue_serotype_neg_target)



dir.create("Data/model_prepared_gm", showWarnings = FALSE, recursive = TRUE)
 
saveRDS(data_with_binomial_targets_flavi_list,
        file.path("Data/model_prepared_gm", "data_with_binomial_targets_flavi_list.rds"))
saveRDS(data_with_binomial_targets_chik_list,
        file.path("Data/model_prepared_gm", "data_with_binomial_targets_chik_list.rds"))
saveRDS(data_with_multinomial_targets_list,
        file.path("Data/model_prepared_gm", "data_with_multinomial_targets_list.rds"))
 


results_flavi_vs_not_list   <- list(ratio = list(), cross_sectional = list())
results_dengue_vs_not_list  <- list(ratio = list(), cross_sectional = list())
results_chik_vs_not_list    <- list(ratio = list(), cross_sectional = list())
results_dengue_vs_chik_list <- list(ratio = list(), cross_sectional = list())



# --- NEW BINARY CLASSIFICATION ----
for (dataset_type in dataset_types) {
  for (isotype in isotypes) {
    binom_flavi <- data_with_binomial_targets_flavi_list[[dataset_type]][[isotype]]
    binom_chik  <- data_with_binomial_targets_chik_list[[dataset_type]][[isotype]]

    results_flavi_vs_not_list[[dataset_type]][[isotype]] <- train_binary_models(
      binom_flavi$data$flavi, "flavi_target",
      positive_class = binom_flavi$positive_class_map$flavi, cv_repeats = 10)

    results_dengue_vs_not_list[[dataset_type]][[isotype]] <- train_binary_models(
      binom_flavi$data$dengue, "dengue_target",
      positive_class = binom_flavi$positive_class_map$dengue, cv_repeats = 10)

    results_chik_vs_not_list[[dataset_type]][[isotype]] <- train_binary_models(
      binom_chik$data$chik, "chik_target",
      positive_class = binom_chik$positive_class_map$chik, cv_repeats = 10)
  }
}



# -- Save binomial

# ----Flavi vs Not 
# Get best performing tuned model for each isotype - ratios
results_flavi_vs_not_list$ratio$IgG$comparison
results_dengue_vs_not_list$ratio$IgG$comparison
results_chik_vs_not_list$ratio$IgG$comparison

results_flavi_vs_not_list$ratio$IgM$comparison
results_dengue_vs_not_list$ratio$IgM$comparison


# Cross-sectional
results_flavi_vs_not_list$cross_sectional$IgG$comparison
results_dengue_vs_not_list$cross_sectional$IgG$comparison
#results_flavi_vs_not_list$cross_sectional$IgA$comparison
results_flavi_vs_not_list$cross_sectional$IgM$comparison
results_dengue_vs_not_list$cross_sectional$IgM$comparison

results_chik_vs_not_list$cross_sectional$IgG$comparison

# Save all results (binomial)
dir.create("Results/NEW_Binary_Classification", showWarnings = FALSE, recursive = TRUE)
base_dir <- "Results/NEW_Binary_Classification"

for (dataset_type in c("ratio", "cross_sectional")) {
  dir.create(file.path(base_dir, dataset_type), recursive = TRUE, showWarnings = FALSE)
}



# --- read results files 
base_dir  <- here("Results/NEW_Binary_Classification")
types     <- c("ratio", "cross_sectional")
questions <- c("flavi_vs_not", "dengue_vs_not", "chik_vs_not", "dengue_vs_chik")
isotypes  <- c("IgG", "IgA", "IgM", "avidity")



# ---- Bundle the four result lists together, keyed by question name ----
all_results <- list(
  flavi_vs_not   = results_flavi_vs_not_list,
  dengue_vs_not  = results_dengue_vs_not_list,
  chik_vs_not    = results_chik_vs_not_list,
  dengue_vs_chik = results_dengue_vs_chik_list
)


# ---- Save each result object and collect comparison tables ----
comparison_rows <- list()

for (question in names(all_results)) {
  for (dataset_type in names(all_results[[question]])) {
    for (isotype in names(all_results[[question]][[dataset_type]])) {

      result_obj <- all_results[[question]][[dataset_type]][[isotype]]
      if (is.null(result_obj)) next

      label <- paste(question, dataset_type, isotype, sep = " / ")
      cat("\n==== ", label, " ====\n", sep = "")

      out_path <- file.path(base_dir, dataset_type,
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
      comp$n_resamples  <- result_obj$n_resamples
      comparison_rows[[paste(question, dataset_type, isotype, sep = "_")]] <- comp
    }
  }
}

# ---- Combine into one long df ----
combined_comparisons <- dplyr::bind_rows(comparison_rows) %>%
  dplyr::select(dataset_type, question, isotype, Model,
                ROC, ROC_low, ROC_high,
                ROC_mean_fold, ROC_sd_fold, n_resamples,
                dplyr::everything())

for (ty in c("ratio", "cross_sectional")) {
  write.csv(
    combined_comparisons %>% dplyr::filter(dataset_type == ty),
    file.path(base_dir, ty, "all_comparisons.csv"),
    row.names = FALSE
  )
}

# ---- Summary across the four models, per question/isotype ----
# ---- Summary across the four models, per question/isotype ----
auc_summary_binomial <- combined_comparisons %>%
  dplyr::group_by(dataset_type, question, isotype) %>%
  dplyr::summarise(
    # Pooled out-of-fold AUC (DeLong CI)
    mean_ROC        = mean(ROC, na.rm = TRUE),
    min_ROC         = min(ROC, na.rm = TRUE),
    max_ROC         = max(ROC, na.rm = TRUE),
    range_ROC       = max(ROC, na.rm = TRUE) - min(ROC, na.rm = TRUE),
    sd_ROC          = stats::sd(ROC, na.rm = TRUE),
    best_model      = Model[which.max(ROC)],
    worst_model     = Model[which.min(ROC)],
    # envelope of the per-model CIs - NOT a CI for the mean
    env_low         = min(ROC_low,  na.rm = TRUE),
    env_high        = max(ROC_high, na.rm = TRUE),

    # Mean per-fold AUC (Nadeau-Bengio corrected CI)
    mean_ROC_fold   = mean(ROC_mean_fold, na.rm = TRUE),
    min_ROC_fold    = min(ROC_mean_fold,  na.rm = TRUE),
    max_ROC_fold    = max(ROC_mean_fold,  na.rm = TRUE),
    best_model_fold = Model[which.max(ROC_mean_fold)],
    best_ROC_fold   = max(ROC_mean_fold, na.rm = TRUE),
    best_fold_low   = ROC_fold_low[which.max(ROC_mean_fold)],   # CI of the best model
    best_fold_high  = ROC_fold_high[which.max(ROC_mean_fold)],
    # envelope of the per-model corrected CIs - NOT a CI for the mean
    env_fold_low    = min(ROC_fold_low,  na.rm = TRUE),
    env_fold_high   = max(ROC_fold_high, na.rm = TRUE),

    n_models        = dplyr::n(),
    n_resamples     = dplyr::first(n_resamples),
    .groups = "drop"
  ) %>%
  dplyr::arrange(dataset_type, question, isotype)

write.csv(auc_summary_binomial,
          file.path(base_dir, "binomial_auc_summary_across_models.csv"),
          row.names = FALSE)

View(auc_summary_binomial)





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
      axis.title.y     = element_text(size = 18), 
      legend.text      = element_text(size = 20),
      panel.grid.minor = element_blank(),
      strip.text       = element_text(size = 20),
      strip.background = element_rect(fill = "#ffffff", colour = NA),
      legend.position  = "bottom", aspect.ratio = 0.5
    )
}



quartz()
plot_auc_bar(auc_avg_ratio)
plot_auc_bar(auc_avg_cross_sectional)
plot_auc_bar(auc_avg_ratio_multinomial)
plot_auc_bar(auc_avg_cross_sectional_multinomial)





# extra code 

results_flavi_vs_not_list$ratio$IgG$models$glmnet$bestTune
results_flavi_vs_not_list$ratio$IgG$target_levels
m <- results_flavi_vs_not_list$ratio$IgG$models$glmnet
coef(m$finalModel, s = m$bestTune$lambda)

res <- results_flavi_vs_not_list$ratio$IgG
aggregate(prob_positive ~ obs_class, data = res$predictions, FUN = mean)

res$comparison

library(dplyr)
res$predictions %>%
  filter(Model == "GLMnet") %>%
  group_by(rowIndex, obs_class) %>%
  summarise(p = mean(prob_positive), .groups = "drop") %>%
  ggplot(aes(obs_class, p)) + geom_boxplot() + geom_jitter(width = 0.2, alpha = 0.4)

lapply(res$models[c("rf","pls")], function(m) m$bestTune)
res$models$rf$results[match_best_row(res$models$rf$results, res$models$rf$bestTune, "rf"), ]


res$predictions %>%
  filter(Model == "GLMnet") %>%
  group_by(rowIndex, obs_class) %>%
  summarise(p = mean(prob_positive), .groups = "drop") %>%
  with(roc(obs_class, p, levels = c("negative","positive"), direction = "<")) %>%
  ci.auc()

res$predictions %>%
  filter(Model == "GLMnet") %>%
  group_by(rowIndex, obs_class) %>%
  summarise(p = mean(prob_positive), .groups = "drop") %>%
  filter((obs_class == "positive" & p < 0.58) | (obs_class == "negative" & p > 0.55))

plot(m$finalModel, xvar = "lambda", label = TRUE)   # coefficient paths vs log λ
abline(v = log(m$bestTune$lambda), lty = 2)




# Binomial Classification 
for (dataset_type in dataset_types) {
  for (isotype in isotypes) {

    binom_flavi <- data_with_binomial_targets_flavi_list[[dataset_type]][[isotype]]
    binom_chik  <- data_with_binomial_targets_chik_list[[dataset_type]][[isotype]]

    # Flavi vs not
    d <- binom_flavi$data$flavi
    set.seed(42)
    folds <- caret::groupKFold(d$id_patient, k = 5)
    results_flavi_vs_not <- train_binary_models(
      data = binom_flavi$data$flavi,
      target = "flavi_target",
      positive_class = binom_flavi$positive_class_map$flavi,
      variables = NULL,cv_repeats = 10,  fold_index = folds,
      metrics = c("ROC", "AUPRC", "Brier")
    )

    # Dengue vs not
    d <- binom_flavi$data$dengue
    set.seed(42)
    folds <- caret::groupKFold(d$id_patient, k = 5)
    results_dengue_vs_not <- train_binary_models(
      data = binom_flavi$data$dengue,
      target = "dengue_target",
      positive_class = binom_flavi$positive_class_map$dengue,
      variables = NULL, cv_repeats = 10, fold_index = folds,
      metrics = c("ROC", "AUPRC", "Brier")
    )

    # Chik vs not
    d <- binom_chik$data$chik
    set.seed(42)
    folds <- caret::groupKFold(d$id_patient, k = 5)
    results_chik_vs_not <- train_binary_models(
      data = binom_chik$data$chik,
      target = "chik_target",
      positive_class = binom_chik$positive_class_map$chik,
      variables = NULL, cv_repeats = 10,    fold_index = folds,
      metrics = c("ROC", "AUPRC", "Brier")
    )

    # Dengue vs chik
    #results_dengue_vs_chik <- train_binary_models(
     # data = binom_chik$data$dengue_chik,
     # target = "dengue_chik_target",
     # positive_class = binom_chik$positive_class_map$dengue_chik,
     # variables = NULL,
     # metrics = c("ROC", "AUPRC", "Brier")
    #)

    results_flavi_vs_not_list[[dataset_type]][[isotype]]   <- results_flavi_vs_not
    results_dengue_vs_not_list[[dataset_type]][[isotype]]  <- results_dengue_vs_not
    results_chik_vs_not_list[[dataset_type]][[isotype]]    <- results_chik_vs_not
    #results_dengue_vs_chik_list[[dataset_type]][[isotype]] <- results_dengue_vs_chik
  }
}
