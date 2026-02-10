# --- Script for preprocessing data to run CHIKV vs Dengue
install.packages("readxl")
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


# cebu data with multiple antigens 
head(cebu_mutiple_antigens)
colnames(cebu_mutiple_antigens)

# Simpler version - just antigen as columns
cebu_pivot <- cebu_mutiple_antigens %>%
  pivot_wider(
    names_from = antigen,
    values_from = RAU
  )
View(cebu_pivot)
class(cebu_pivot$date_sample)
table(cebu_pivot$PCR)


cebu_pivot_days_since_inf  <- cebu_pivot %>%
  group_by(id_patient) %>%
  arrange(date_sample) %>%
  mutate(
    # Count actual infections (positive PCRs only)
    n_infections = sum(!is.na(PCR) & PCR != "negative"),
    
    # First PCR test date (regardless of result - includes "negative")
    infection_date = first(date_sample[!is.na(PCR)]),
    
    # Days since first PCR test
    days_since_infection = as.numeric(difftime(date_sample, infection_date, units = "days"))
  ) %>%
  ungroup()
View(cebu_pivot_days_since_inf)
head(cebu_pivot_days_since_inf$id_patient)

validate <- cebu_pivot_days_since_inf %>%
  filter(id_patient %in% c("CPC-C-0010-00", "CPC-C-0047-00")) %>%
  dplyr::select(id_patient, date_sample, PCR, infection_date, days_since_infection, n_infections) %>%
  arrange(id_patient, date_sample)
View(validate)

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

