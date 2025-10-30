# --- Script for preprocessing data to run CHIKV vs Dengue

#--- source functions
source(here('Functions.R'))

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

# import Luminex / MSD data
# assumes data has: MFI values, target infection, patient ID, days since infection
# acute and convalescent samples have been removed 
raw_data_with_chik <- read_csv('luminex_data_with_hai_chik.csv')


# use HI to informed if pcr neg == no infection 
df_with_HI_ratio_chik <- raw_data_with_chik %>%
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
df_with_HI_ratio_chik <- df_with_HI_ratio_chik %>% 
  mutate(new_target = case_when(
    log_mean_HI_ratio < 1.6  ~ 'no_infection', 
    TRUE ~ NA_character_  
  ))


# --- Keep CHIKV as PCR target for dengue vs CHIKV 
df_with_HI_ratio_chik <- df_with_HI_ratio_chik %>%
  mutate(mapped_target = case_when(
    Target %in% c('DENV1', 'DENV2', 'DENV3', 'DENV4', 'CHIKV') ~ Target,
    Target %in% c('ZIKV','negative') ~ new_target,
    TRUE ~ NA_character_
  ))
# Create patient-level mapping based on all rows per patient
patient_mapping_chik <- df_with_HI_ratio_chik %>%
  group_by(id_patient) %>%
  summarise(
    # For new_target: take the first non-NA value if any exists
    patient_new_target = first(na.omit(new_target)),
    # For Target: take the first non-NA value if any exists
    patient_target = first(na.omit(Target))
  ) %>%
  mutate(
    # Apply your mapping logic at patient level
    mapped_target = case_when(
      patient_target %in% c('DENV1', 'DENV2', 'DENV3', 'DENV4', 'CHIKV') ~ patient_target,
      patient_target %in% c('ZIKV', 'negative') ~ patient_new_target,
      TRUE ~ NA_character_
    )
  )

all_processed_dfs_with_chik <- process_luminex_data(raw_data_with_chik, patient_mapping_chik, pre_threshold = -1)

# ---- save pre-processed dfs ----
saveRDS(all_processed_dfs_with_chik, here("luminex_processed_data_with_chik.rds"))

