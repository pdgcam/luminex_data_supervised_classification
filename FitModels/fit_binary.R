library(caret)
library(here)
packageVersion("xgboost")
packageVersion("caret")
library(dplyr)



# --- Source functions
source(here('Models_Supervised/train_binary_models.R'))



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


# -- Ratio (Longitudnal) DENV
logged_ratio_df_IgG_denv <- readRDS("Data/by_isotype_gm/denv/logged_ratio_df_IgG.rds")
logged_ratio_df_IgA_denv <- readRDS("Data/by_isotype_gm/denv/logged_ratio_df_IgA.rds")
logged_ratio_df_IgM_denv <- readRDS("Data/by_isotype_gm/denv/logged_ratio_df_IgM.rds")

logged_ratio_df_IgG_chik <- readRDS("Data/by_isotype_gm/chik/logged_ratio_df_IgG.rds")
logged_ratio_df_IgA_chik <- readRDS("Data/by_isotype_gm/chik/logged_ratio_df_IgA.rds")
logged_ratio_df_IgM_chik <- readRDS("Data/by_isotype_gm/chik/logged_ratio_df_IgM.rds")

logged_ratio_df_IgG_denv_chik <- readRDS("Data/by_isotype_gm/denv_vs_chik/logged_ratio_df_IgG.rds")
logged_ratio_df_IgA_denv_chik <- readRDS("Data/by_isotype_gm/denv_vs_chik/logged_ratio_df_IgA.rds")
logged_ratio_df_IgM_denv_chik <- readRDS("Data/by_isotype_gm/denv_vs_chik/logged_ratio_df_IgM.rds")




# --- Cross-sectional (IgG and IgM)
flavi_cross_sectional_df_IgG <- readRDS("Data/by_isotype_cross_sectional/flavi_cross_sectional_logged_IgG.rds")
flavi_cross_sectional_df_IgM <- readRDS("Data/by_isotype_cross_sectional/flavi_cross_sectional_logged_IgM.rds")
flavi_cross_sectional_df_IgA <- readRDS("Data/by_isotype_cross_sectional/flavi_cross_sectional_logged_IgA.rds")


alpha_cross_sectional_df_IgG <- readRDS("Data/by_isotype_cross_sectional/alpha_cross_sectional_logged_IgG.rds")
alpha_cross_sectional_df_IgM <- readRDS("Data/by_isotype_cross_sectional/alpha_cross_sectional_logged_IgM.rds")
alpha_cross_sectional_df_IgA <- readRDS("Data/by_isotype_cross_sectional/alpha_cross_sectional_logged_IgA.rds")



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


# ---- Lookup of the source dfs already loaded above ----
source_dfs <- list(
  ratio = list(
    denv = list(IgG = logged_ratio_df_IgG_denv,
                IgA = logged_ratio_df_IgA_denv,
                IgM = logged_ratio_df_IgM_denv),
    chik = list(IgG = logged_ratio_df_IgG_chik,
                IgA = logged_ratio_df_IgA_chik,
                IgM = logged_ratio_df_IgM_chik),
    denv_chik = list(IgG = logged_ratio_df_IgG_denv_chik,
                     IgA = logged_ratio_df_IgA_denv_chik,
                     IgM = logged_ratio_df_IgM_denv_chik)
  ),
  cross_sectional = list(
    flavi = list(IgG = flavi_cross_sectional_df_IgG,
                 IgA = flavi_cross_sectional_df_IgA,
                 IgM = flavi_cross_sectional_df_IgM),
    alpha = list(IgG = alpha_cross_sectional_df_IgG,
                 IgA = alpha_cross_sectional_df_IgA,
                 IgM = alpha_cross_sectional_df_IgM)
  )
)

binomial_flavi_targets <- c("flavi", "dengue")
data_with_binomial_targets_flavi_list <- list(ratio = list(), cross_sectional = list())
data_with_binomial_targets_chik_list  <- list(ratio = list(), cross_sectional = list())

for (iso in isotypes) {
  cat(sprintf("\n===== %s =====\n", iso))

  # ---- Binomial: flavi / dengue ----
  # Ratio
  data_with_binomial_targets_flavi_list$ratio[[iso]] <-
    select_targets(source_dfs$ratio$denv[[iso]], targets = binomial_flavi_targets)

  # Cross-sectional (drop ZIKV encoded as negative)
  src_cs_flavi <- source_dfs$cross_sectional$flavi[[iso]]
  res_cs_flavi <- select_targets(src_cs_flavi, targets = binomial_flavi_targets)
  data_with_binomial_targets_flavi_list$cross_sectional[[iso]] <-
    drop_zikv_negatives(res_cs_flavi, src_cs_flavi, targets = binomial_flavi_targets)

  # ---- Binomial: chik / dengue vs chik ----
  # Ratio: chik vs not from the chik df, dengue vs chik from the denv_vs_chik df
  res_chik    <- select_targets(source_dfs$ratio$chik[[iso]],      targets = "chik")
  res_denv_ch <- select_targets(source_dfs$ratio$denv_chik[[iso]], targets = "dengue_chik")
  res_chik$data$dengue_chik <- res_denv_ch$data$dengue_chik
  data_with_binomial_targets_chik_list$ratio[[iso]] <- res_chik

  # Cross-sectional: chik vs not only
  # TODO: dengue_chik for cross-sectional still needs fixing (see note below)
  data_with_binomial_targets_chik_list$cross_sectional[[iso]] <-
    select_targets(source_dfs$cross_sectional$alpha[[iso]], targets = "chik")

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



dir.create("Data/binary_model_prepared", showWarnings = FALSE, recursive = TRUE)
 
saveRDS(data_with_binomial_targets_flavi_list,
        file.path("Data/binary_model_prepared", "data_with_binomial_targets_flavi_list.rds"))
saveRDS(data_with_binomial_targets_chik_list,
        file.path("Data/binary_model_prepared", "data_with_binomial_targets_chik_list.rds"))
saveRDS(data_with_multinomial_targets_list,
        file.path("Data/binary_model_prepared", "data_with_multinomial_targets_list.rds"))
 


results_flavi_vs_not_list   <- list(ratio = list(), cross_sectional = list())
results_dengue_vs_not_list  <- list(ratio = list(), cross_sectional = list())
results_chik_vs_not_list    <- list(ratio = list(), cross_sectional = list())
results_dengue_vs_chik_list <- list(ratio = list(), cross_sectional = list())



# ---  BINARY CLASSIFICATION ----
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


# Save all results (binomial)
base_dir <- "Results/Binary_Classification"

binary_results <- list(
  flavi_vs_not   = results_flavi_vs_not_list,
  dengue_vs_not  = results_dengue_vs_not_list,
  chik_vs_not    = results_chik_vs_not_list,
  dengue_vs_chik = results_dengue_vs_chik_list
)

for (dataset_type in dataset_types) {
  out_dir <- file.path(base_dir, dataset_type)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  for (comparison_name in names(binary_results)) {
    for (isotype in isotypes) {
      res <- binary_results[[comparison_name]][[dataset_type]][[isotype]]
      if (is.null(res)) next  # skip models that weren't trained

      # Full results object
      saveRDS(res, file.path(out_dir, sprintf("%s_%s.rds", comparison_name, isotype)))

      # Model comparison table as CSV for quick viewing
      if (is.data.frame(res$comparison)) {
        write.csv(res$comparison,
                  file.path(out_dir, sprintf("%s_%s_comparison.csv", comparison_name, isotype)),
                  row.names = FALSE)
      }
    }
  }
}

# Combined file per comparison (all dataset types and isotypes together)
for (comparison_name in names(binary_results)) {
  saveRDS(binary_results[[comparison_name]],
          file.path(base_dir, sprintf("results_%s_list.rds", comparison_name)))
}



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
results_flavi_vs_not_list$cross_sectional$IgA$comparison
results_flavi_vs_not_list$cross_sectional$IgM$comparison
results_dengue_vs_not_list$cross_sectional$IgM$comparison

results_chik_vs_not_list$cross_sectional$IgG$comparison

