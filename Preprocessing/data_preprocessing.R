# --- Script for preprocessing data and visualising tire dynamics 
# Import libraries
library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(matrixStats)
library(stringr)
library(data.table)
library(dplyr)
library(scales)
library(purrr)
library(tidyr)
library(stringr)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(dplyr)
library(here)
library(terra)
library(exactextractr)
library(raster)
library(readxl)
library(knitr)
library(patchwork)
library(tidyverse)

#--- source functions
source(here('/Users/ap2488/Documents/GitHub/luminex_data_supervised_classification/Functions.R'))


# Functions
convert_hai_to_numeric <- function(x) {
  case_when(
    is.na(x) ~ NA_real_,
    x == "NA" ~ NA_real_,
    x == "<10" ~ 5,  # Treat <10 as 5 (undetectable HI)
    TRUE ~ as.numeric(as.character(x))
  )
}

#--- log transform of HI values
log_transform <- function (titre) {
  1 + (log2(titre / 10) / log2(2))
}


# --- post/pre dataset + cross-sectional dataset
prepare_luminex_datasets <- function(raw_data, patient_mapping, antigen_cols, pre_threshold = -1) {
  
  # Validate that antigen columns exist in the data
  missing_cols <- setdiff(antigen_cols, names(raw_data))
  if (length(missing_cols) > 0) {
    stop("Missing columns in data: ", paste(missing_cols, collapse = ", "))
  }

   # ---- Dataset 1: Post/Pre Ratio ----
  ratio_df <- raw_data %>%
    # Label pre vs post
    mutate(timepoint = ifelse(days_since_infection <= pre_threshold, 'pre', 'post')) %>%
    # Average across all pre and post samples for each anitgen and each patient 
    group_by(id_patient, timepoint) %>%
    summarise(across(all_of(antigen_cols), ~mean(.x, na.rm = TRUE)), .groups = 'drop') %>%
    # Pivot to get pre and post columns
    pivot_wider(
      names_from = timepoint,
      values_from = all_of(antigen_cols),
      names_sep = "_"
    ) %>%
    # Calculate ratios
    mutate(across(
      .cols = ends_with("_post"),
      .fns = ~. / get(sub("_post$", "_pre", cur_column())),
      .names = "{sub('_post$', '', .col)}"
    )) %>%
    # Keep only ratio columns + id
    dplyr::select(id_patient, all_of(antigen_cols)) %>%
    # Add target
    left_join(patient_mapping, by = "id_patient") %>%
    # Remove rows with any NA
    na.omit() %>%
    # Convert to dataframe with rownames
    as.data.frame()
  
  # Dataset 2: Last draw (ie draw at max days since infection)
  cross_sectional_data <- raw_data %>%
    group_by(id_patient) %>%
    filter(days_since_infection == max(days_since_infection, na.rm = TRUE)) %>%
    ungroup() %>%
    # Keep only selected antigen columns
    dplyr::select(id_patient, all_of(antigen_cols)) %>%
    # Add target
    left_join(patient_mapping, by = "id_patient") %>%
    na.omit() %>%
    as.data.frame()
  
 
  # Return both datasets
  return(list(
    ratio = ratio_df,
    cross_sectional_data = cross_sectional_data
  ))
}



# validation and random datasets (cebu)
# validation subset = PCR confirmed cases
validation_subset <- read.csv("/Users/ap2488/Desktop/supervised_learning_flavi/MIA_DataBaseOut_ValidationSet.csv")
random_subset <- read.csv("/Users/ap2488/Desktop/supervised_learning_flavi/MIA_DataBaseOut_RandomSubset.csv")
cebu_mutiple_antigens <- read_excel("/Users/ap2488/Desktop/supervised_learning_flavi/db_philippines_IgG_IgA_IgM_avidity.xlsx")

# align patient IDs / PCR cols across datasets
cebu_mutiple_antigens$id_patient <- gsub("_", "-", cebu_mutiple_antigens$id_patient)
length(intersect(validation_subset$ids, cebu_mutiple_antigens$id_patient)) #39 samples intersect


# pivot to get RAU as main data
cebu_pivot <- cebu_mutiple_antigens %>%
  pivot_wider(
    names_from = antigen,
    values_from = RAU
  )

View(cebu_pivot)


# --- Add col: days_since_infection 
cebu_pivot_days_since_inf  <- cebu_pivot %>%
  group_by(id_patient) %>%
  arrange(date_sample) %>%
  mutate(
    # Count number of samples/timepoints per patient
    n_samples = n(),
    # First PCR test date (regardless of result - includes "negative")
    infection_date = first(date_sample[!is.na(PCR)]),
    # Days since first PCR test
    days_since_infection = as.numeric(difftime(date_sample, infection_date, units = "days"))
  ) %>%
  ungroup()


# --- Use HI threshold (1.6) to remove subclinical infection 
# use HI to informed if pcr neg == no infection 
HI_ratio_df <- cebu_pivot_days_since_inf %>%
  mutate(across(starts_with("HAI_"), convert_hai_to_numeric)) %>%
  arrange(id_patient, days_since_infection) %>%
  group_by(id_patient) %>%
  mutate(
    # --- log transform HI
    log_HI_DENV1 = if_else(HAI_DENV1 > 0 ,log_transform(HAI_DENV1), 
                                 NA_real_),
    log_HI_DENV2 = if_else(HAI_DENV2 > 0, log_transform(HAI_DENV2), 
                                 NA_real_),
    log_HI_DENV3 = if_else(HAI_DENV3 > 0, log_transform(HAI_DENV3), 
                                 NA_real_),
    log_HI_DENV4 = if_else(HAI_DENV4 > 0 ,log_transform(HAI_DENV4), 
                                 NA_real_),
    log_HI_CHIKV = if_else(HAI_CHIKV > 0 ,log_transform(HAI_CHIKV), 
    NA_real_),
    
     # --- mean across serotypes 
    mean_log_HI = rowMeans(
      cbind(log_HI_DENV1, log_HI_DENV2, 
            log_HI_DENV3, log_HI_DENV4),
      na.rm = TRUE
    ),
    
    # --- ratio T2 / T1
    prev_days_since_infection = lag(days_since_infection),
    time_diff = days_since_infection - lag(days_since_infection),
    log_mean_HI_ratio_deng = mean_log_HI - lag(mean_log_HI), 
    log_mean_HI_ratio_chik = log_HI_CHIKV - lag(log_HI_CHIKV)
  ) %>%
  ungroup()  %>% dplyr::select(-prev_days_since_infection, -time_diff)



# use ratio threshold == 1.6
HI_ratio_df <- HI_ratio_df %>%
  group_by(id_patient) %>%
  mutate(
    # Check if patient was EVER PCR positive
    ever_pcr_positive = any(!is.na(PCR) & PCR != "negative"),
    
    infection_status = case_when(
      # check subclinical infection (regardless of PCR)
      log_mean_HI_ratio_deng >= log2(1.6) &  log_mean_HI_ratio_chik >= log2(1.6) ~ 'subclinical_deng_chik',
      log_mean_HI_ratio_deng >= log2(1.6) ~ 'subclinical_dengue',
      log_mean_HI_ratio_chik >= log2(1.6) ~ 'subclinical_chik',
      # Then check if ever PCR positive (without subclinical HI rise)
      ever_pcr_positive ~ 'true_positive',
      # Remaining PCR-negative patients with no HI rise
      log_mean_HI_ratio_deng < log2(1.6) & log_mean_HI_ratio_chik < log2(1.6) ~ 'true_negative',
      # Handle NAs (first timepoint, missing HI data)
      TRUE ~ NA_character_
    )
  ) %>% ungroup() %>% dplyr::select(-ever_pcr_positive)


table(HI_ratio_df$infection_status)

# drop rows (timepoints) that are suspected to be subclinical using the 1.6 threshold 
clean_HI_ratio_df <- HI_ratio_df %>%
  filter(!grepl("subclinical", infection_status) | is.na(infection_status))

nrow(clean_HI_ratio_df)
nrow(HI_ratio_df)

# Identify patient IDs that were excluded 
excluded_rows <- HI_ratio_df %>%
  filter(grepl("subclinical", infection_status))
View(excluded_rows)



#  ---- Antigens and isotypes interest 
antigen_cols <- c(
  "CHIKV_E2","CHIKV_NSP123","CHIKV_VLP","SHERPADES_CHIKV_E2",
  "DENV1_DIII","DENV1_NS1","DENV1_VLP","SHERPADES_DENV1_DIII",
  "DENV2_DIII","DENV2_NS1","DENV2_VLP", "SHERPADES_DENV2_DIII",
  "DENV3_DIII","DENV3_NS1","DENV3_VLP", "SHERPADES_DENV3_DIII",
  "DENV4_DIII","DENV4_NS1","DENV4_VLP",  "SHERPADES_DENV4_DIII",
  "JEV_E","JEV_NS1","SHERPADES_JEV_DIII",
  "MAYV_E2","SHERPADES_MAYV_E2",
  "ONNV_E2","ONNV_VLP",
  "RR", "SHERPADES_RR",
  "WNV_DIII","WNV_NS1","SHERPADES_WNV_DIII",
  "YFV_E","YFV_NS1", "SHERPADES_YFV_DIII",
  "ZIKV_NS1","ZIKV_VLP","ZIKVAS_DIII", "ZIKVSU_NS1","SHERPADES_ZIKV_DIII"
)

isotypes <- c("IgG", "IgA", "IgM", "avidity") # not using avidity currently 



# rename + IgG, IgA and IgM isotypes
final_preprocessed_data <- clean_HI_ratio_df
final_preprocessed_data <- final_preprocessed_data %>% filter(isotype %in% isotypes)

View(final_preprocessed_data)


# keep unlogged data 
final_preprocessed_data_raw <- final_preprocessed_data

# Create the logged version as a separate dataframe
final_preprocessed_data_log <- final_preprocessed_data
final_preprocessed_data_log[antigen_cols] <- log2(final_preprocessed_data_log[antigen_cols])



# Create patient-level mapping
patient_pcr_mapping <- final_preprocessed_data_raw %>%
  group_by(id_patient) %>%
  summarise(
    # Check what PCR results this patient has
    has_zikv = any(PCR == 'ZIKV', na.rm = TRUE),
    has_chikv = any(PCR == 'CHIKV', na.rm = TRUE),
    has_denv = any(PCR %in% c('DENV1', 'DENV2', 'DENV3', 'DENV4'), na.rm = TRUE),
    
    # Get first PCR result for reference
    patient_pcr = first(na.omit(PCR)),
    
    # Get infection status
    patient_infection_status = first(na.omit(infection_status))
  ) %>%
  mutate(
    # Map to target label with ZIKV priority
    target = case_when(
      # Prioritize ZIKV if present
      has_zikv ~ 'ZIKV',
      
      # Then DENV
      has_denv ~ patient_pcr,
      
      # Then CHIKV
      has_chikv ~ 'CHIKV',
      
      # True negatives
      patient_infection_status == 'true_negative' ~ 'negative',
      
      TRUE ~ NA_character_
    )
  ) %>% dplyr::select(id_patient, target)  

# negative == NO PCR confirmed infection + subclinical (during study) removed 
table(patient_pcr_mapping$target)
target_counts <- as.data.frame(table(patient_pcr_mapping$target))
colnames(target_counts) <- c("Target", "Count")

# save counts per target 
write.csv(target_counts, "Results/target_counts.csv", row.names = FALSE)

# save preprocessed data for downstream analysis 
write.csv(final_preprocessed_data_raw,
"Results/raw_preprocessed_cebu_data.csv")

write.csv(final_preprocessed_data_log,
"Results/logged_preprocessed_cebu_data.csv")



# --- calculate  post / pre ratios and extract crosssectional data (use unloged data)
processed_dfs <- prepare_luminex_datasets(final_preprocessed_data_raw, patient_pcr_mapping, antigen_cols, pre_threshold = -1)


# Save each dataset separately (easier to load individually later)
saveRDS(processed_dfs$ratio,"Results/ratio_df.rds")
saveRDS(processed_dfs$cross_sectional_data, "Results/cross_sectional_df.rds")
