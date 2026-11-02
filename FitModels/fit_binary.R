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
