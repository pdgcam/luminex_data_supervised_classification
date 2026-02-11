# --- Script for preprocessing data to run CHIKV vs Dengue
install.packages("gtsummary")
install.packages("rlang")
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


#--- source functions
source(here('/Users/ap2488/Documents/GitHub/updates_uminex_data_supervised_classification/Functions.R'))



# ---  function to convert HAI values to numeric
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
  1 + (log(titre / 10) / log(2))
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


# --- Add col: days_since_infection 
cebu_pivot_days_since_inf  <- cebu_pivot %>%
  group_by(id_patient) %>%
  arrange(date_sample) %>%
  mutate(
    # Count number of samples/timepoints per patient
    n_samples = n(),
    # Count actual infections (positive PCRs only)
    n_infections = sum(!is.na(PCR) & PCR != "negative"),
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
    
     # --- mean across serotypes 
    mean_log_HI = rowMeans(
      cbind(log_HI_DENV1, log_HI_DENV2, 
            log_HI_DENV3, log_HI_DENV4),
      na.rm = TRUE
    ),
    
    # --- ratio T2 / T1
    prev_days_since_infection = lag(days_since_infection),
    time_diff = days_since_infection - lag(days_since_infection),
    log_mean_HI_ratio = mean_log_HI - lag(mean_log_HI) 
  ) %>%
  ungroup()


# use ratio threshold == 1.6
HI_ratio_df <- HI_ratio_df %>%
  group_by(id_patient) %>%
  mutate(
    # Check if patient was EVER PCR positive
    ever_pcr_positive = any(!is.na(PCR) & PCR != "negative"),
    
    # Only flag subclinical infections in patients who were always PCR negative
    infection_status = case_when(
      # If ever PCR positive, mark as PCR positive (regardless of HI)
      ever_pcr_positive ~ 'PCR_positive',
      
      # For always-negative patients, check HI ratio for subclinical infection
      !ever_pcr_positive & log_mean_HI_ratio >= log2(1.6) ~ 'subclinical_infection',
      !ever_pcr_positive & log_mean_HI_ratio < log2(1.6) ~ 'true_negative',
      
      # Handle NAs (first timepoint, missing HI data)
      TRUE ~ NA_character_
    )
  ) %>%
  ungroup()


# Identify patients to exclude
patients_to_exclude <- HI_ratio_df %>%
  group_by(id_patient) %>%
  summarise(
    has_subclinical = any(infection_status == 'subclinical_infection', na.rm = TRUE),
    all_na = all(is.na(infection_status))
  ) %>%
  filter(has_subclinical | all_na) %>%
  pull(id_patient)

# Clean dataframe
clean_HI_ratio_df <- HI_ratio_df %>%
  filter(!id_patient %in% patients_to_exclude)
length(unique(clean_HI_ratio_df$id_patient))

# excluded IDs: 
# "CPC-C-0277-00" #all HI log mean = NA unknown if subclinical
# "CPC-C-0329-00" #all HI log mean = NA unknown if subclinical
# "CPC-C-0068-00" # subclinical using log(1.6) threshold


table(clean_HI_ratio_df$PCR)

# Create patient-level mapping
patient_pcr_mapping <- clean_HI_ratio_df %>%
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

table(patient_pcr_mapping$target)

# table of samples per pathogen (Fig1)
patient_pcr_mapping %>%
  count(target) %>%
  mutate(percent = round(100 * n / sum(n), 1),
         `n (%)` = paste0(n, " (", percent, "%)")) %>%
  select(`PCR Target` = target, `n (%)`) %>%
  kable(caption = "Distribution of PCR Targets")



all_preprocessed_dfs <- process_luminex_data(clean_HI_ratio_df, patient_pcr_mapping, pre_threshold = -1)

# ---- save pre-processed dfs ----
saveRDS(all_processed_dfs_with_chik, here("luminex_processed_data_with_chik.rds"))

