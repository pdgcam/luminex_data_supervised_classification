
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
library(dplyr)
library(here)
library(exactextractr)
library(raster)
library(readxl)
library(knitr)
library(patchwork)
library(tidyverse)
library(standardize)



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

# --- post/pre dataset + cross-sectional dataset (remove any samples between day 0 and day 30)
prepare_luminex_datasets <- function(raw_data, patient_mapping, antigen_cols, pre_threshold = -1, post_threshold = 30) {
  
  # Validate that antigen columns exist in the data
  missing_cols <- setdiff(antigen_cols, names(raw_data))
  if (length(missing_cols) > 0) {
    stop("Missing columns in data: ", paste(missing_cols, collapse = ", "))
  }


   # ---- Dataset 1: Post/Pre Ratio ----
  ratio_df <- raw_data %>%
  filter(days_since_infection != 0) %>%
  mutate(timepoint = case_when(
          days_since_infection <= pre_threshold  ~ "pre",
          days_since_infection >  post_threshold ~ "post",
          TRUE ~ NA_character_)) %>%
  group_by(id_patient, timepoint) %>%
    summarise(across(all_of(antigen_cols), ~exp(mean(log(.x), na.rm = TRUE))), .groups = 'drop') %>%
  pivot_wider(
    names_from = timepoint,
    values_from = all_of(antigen_cols),
    names_sep = "_"
  )

  # Compute ratios 
  post_cols <- paste0(antigen_cols, "_post")
  pre_cols  <- paste0(antigen_cols, "_pre")

  ratio_matrix <- ratio_df[post_cols] / ratio_df[pre_cols]
  names(ratio_matrix) <- antigen_cols

  ratio_df <- bind_cols(
    ratio_df %>% dplyr::select(id_patient),
    ratio_matrix
  ) %>%
    left_join(patient_mapping, by = "id_patient") %>%
    na.omit() %>%
    as.data.frame()
  

  # Return both datasets
  return(list(
    ratio = ratio_df
  ))
}

# thresholds
DENV_THRESHOLD <- log2(1.6)   # >=1.6 rise in mean log2 HAI between consecutive timepoints
CHIK_THRESHOLD <- log2(8)     # >=8-fold rise in PRNT vs enrolment baseline
DENV_SEROTYPES <- c("DENV1", "DENV2", "DENV3", "DENV4")

isotypes     <- c("IgG", "IgA", "IgM", "avidity")
antigen_cols <- c(
  "CHIKV_E2","CHIKV_NSP123","CHIKV_VLP","SHERPADES_CHIKV_E2",
  "DENV1_DIII","DENV1_NS1","DENV1_VLP","SHERPADES_DENV1_DIII",
  "DENV2_DIII","DENV2_NS1","DENV2_VLP","SHERPADES_DENV2_DIII",
  "DENV3_DIII","DENV3_NS1","DENV3_VLP","SHERPADES_DENV3_DIII",
  "DENV4_DIII","DENV4_NS1","DENV4_VLP","SHERPADES_DENV4_DIII",
  "JEV_E","JEV_NS1","SHERPADES_JEV_DIII",
  "MAYV_E2","SHERPADES_MAYV_E2",
  "ONNV_E2","ONNV_VLP",
  "RR","SHERPADES_RR",
  "WNV_DIII","WNV_NS1","SHERPADES_WNV_DIII",
  "YFV_E","YFV_NS1","SHERPADES_YFV_DIII",
  "ZIKV_NS1","ZIKV_VLP","ZIKV_DIII", "ZIKVSU_NS1","SHERPADES_ZIKV_DIII"
)


# --- Read data 
# validation and random datasets (cebu) - validation subset = PCR confirmed cases
validation_subset <- read.csv(here("Data/MIA_DataBaseOut_ValidationSet.csv"))
random_subset <- read.csv(here("Data/MIA_DataBaseOut_RandomSubset.csv"))
cebu_mutiple_antigens <- read_excel(here("Data/db_philippines_IgG_IgA_IgM_avidity.xlsx"))
cpc_gps <- read.csv(here("Data/CPC_GPS.csv"))


# align patient IDs / PCR cols across datasets
cebu_mutiple_antigens$id_patient <- gsub("_", "-", cebu_mutiple_antigens$id_patient)
length(intersect(validation_subset$ids, cebu_mutiple_antigens$id_patient)) #39 samples intersect

unique(cebu_mutiple_antigens$antigen)

# using ZIKVAS_DIII as ZIKV_DIII
# for ZIKA NS1 using ZIKV_NS1 (not using - ZIKVSU_NS1)
cebu_mutiple_antigens <- cebu_mutiple_antigens %>%
  mutate(antigen = dplyr::recode(antigen, "ZIKVAS_DIII" = "ZIKV_DIII"))


# pivot to get RAU as main data
cebu_pivot <- cebu_mutiple_antigens %>%
  pivot_wider(
    names_from = antigen,
    values_from = RAU
  )

View(cebu_pivot)
table(cebu_pivot$PCR)



# --- Add col: days_since_infection 
cebu_pivot_days_since_inf <- cebu_pivot %>%
  group_by(id_patient) %>%
  arrange(date_sample) %>%
  mutate(
    n_samples = n(),
    first_positive_date = dplyr::first(date_sample[PCR %in% c("DENV1", "DENV2", "DENV3", "DENV4", "ZIKV", "CHIKV")],
                                 default = as.Date(NA)),                            
    first_negative_date = dplyr::first(date_sample[PCR %in% "negative"],
                                 default = as.Date(NA)),
    infection_date = if_else(!is.na(first_positive_date), first_positive_date, first_negative_date),
    days_since_infection = as.numeric(difftime(date_sample, infection_date, units = "days"))
  ) %>%
  ungroup() %>%
  dplyr::select(-first_positive_date, -first_negative_date)



# patients with no PCR sample at all -> days_since_infection is NA throughout
cebu_pivot_days_since_inf %>%
  dplyr::group_by(id_patient) %>%
  dplyr::summarise(has_pcr = any(!is.na(PCR)), .groups = "drop") %>%
  dplyr::count(has_pcr)


pcr_history <- cebu_pivot_days_since_inf %>%
  group_by(id_patient) %>%
  summarise(
    # %in% so the many NA PCR rows give FALSE, not NA
    has_symptomatic_denv = any(PCR %in% DENV_SEROTYPES),
    has_symptomatic_chik = any(PCR %in% "CHIKV"),
    .groups = "drop"
  )
 

dplyr::count(pcr_history, has_symptomatic_denv, has_symptomatic_chik)



# isotype-long subsets
no_denv_pcr  <- cebu_pivot_days_since_inf %>%
  semi_join(filter(pcr_history, !has_symptomatic_denv), by = "id_patient")

no_chikv_pcr <- cebu_pivot_days_since_inf %>%
  semi_join(filter(pcr_history, !has_symptomatic_chik), by = "id_patient")
 
n_distinct(no_denv_pcr$id_patient)    # patients with no DENV +ve PCR during study 
n_distinct(no_chikv_pcr$id_patient)   # patients with no CHIKV +ve PCR during study 



# get one sample per patient for HAI and PRNT 
# rather than duplicated, which is the case in original data because we have 4 isotypes per sample)
sample_per_patient <- function(df) {
  out <- df %>%
    dplyr::mutate(
      dplyr::across(starts_with("HAI_"), convert_hai_to_numeric),
      PRNT_CHIKV = convert_hai_to_numeric(PRNT_CHIKV)
    ) %>%
    dplyr::distinct(id_patient, id_sample, date_sample, days_since_infection, PCR,
             HAI_DENV1, HAI_DENV2, HAI_DENV3, HAI_DENV4, PRNT_CHIKV)
 
  # must be one row per sample
  dups <- out %>% dplyr::count(id_patient, id_sample) %>% filter(n > 1)
  if (nrow(dups) > 0) stop("Titres differ between isotype rows for ",
                           nrow(dups), " patients.")
 
  # one sample per patient per day (otherwise can't compute fold-rise)
  same_day <- out %>% dplyr::count(id_patient, date_sample) %>% filter(n > 1)
  if (nrow(same_day) > 0) warning(nrow(same_day), " patient-dates have >1 sample.")
 
  out
}



denv_HAIs  <- sample_per_patient(no_denv_pcr)
chik_PRNTs <- sample_per_patient(no_chikv_pcr)

View(denv_HAIs)
exists(".vsc.attach")

# DENV Subclinicals 
denv_subclinicals_df <- denv_HAIs %>%
  arrange(id_patient, date_sample) %>%     # date_sample is always present
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
    days_gap  = as.numeric(difftime(date_sample, lag(date_sample), units = "days")),
    log_mean_HI_ratio_deng = mean_log_HI - lag(mean_log_HI),
 
    denv_subclinical = coalesce(log_mean_HI_ratio_deng >= DENV_THRESHOLD, FALSE)
  ) %>%
  ungroup() %>%
  dplyr::select(id_patient, id_sample, mean_log_HI, date_sample,days_gap, 
                log_mean_HI_ratio_deng, denv_subclinical)

table(denv_subclinicals_df$denv_subclinical)


# CHIKV subclincicals
chik_baseline <- chik_PRNTs %>%
  filter(!is.na(PRNT_CHIKV)) %>%
  group_by(id_patient) %>%
  slice_min(date_sample, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  transmute(
    id_patient,
    baseline_PRNT_CHIKV = PRNT_CHIKV,
    baseline_date_chik  = date_sample,
    baseline_days_chik  = days_since_infection
  )
 
chik_subclinicals_df <- chik_PRNTs %>%
  left_join(chik_baseline, by = "id_patient") %>%
  mutate(
    log_PRNT_CHIKV             = if_else(PRNT_CHIKV > 0, log_transform(PRNT_CHIKV), NA_real_),
    baseline_log_PRNT_CHIKV    = log_transform(baseline_PRNT_CHIKV),
    log_fold_rise_chik         = log_PRNT_CHIKV - baseline_log_PRNT_CHIKV,
    chik_baseline_seropositive = baseline_PRNT_CHIKV >= 10,   # >=10 = past CHIKV infection
 
    # baseline must be <10 (no past CHIKV) AND a >=8-fold rise from it
    chik_subclinical = coalesce(!chik_baseline_seropositive, FALSE) &
                   coalesce(log_fold_rise_chik >= CHIK_THRESHOLD, FALSE)
  ) %>%
  dplyr::select(id_patient, id_sample, log_PRNT_CHIKV, baseline_PRNT_CHIKV,
                baseline_days_chik, chik_baseline_seropositive,
                log_fold_rise_chik, chik_subclinical)


table(chik_subclinicals_df$chik_subclinical)



subclinical_samples <- full_join(
  denv_subclinicals_df %>%
    filter(denv_subclinical) %>%
    dplyr::select(id_patient, id_sample, days_gap,
                  mean_log_HI, log_mean_HI_ratio_deng),
  chik_subclinicals_df %>%
    filter(chik_subclinical) %>%
    dplyr::select(id_patient, id_sample, log_PRNT_CHIKV,
                  baseline_PRNT_CHIKV, log_fold_rise_chik),
  by = c("id_patient", "id_sample")
) %>%
  mutate(
    reason = case_when(
      !is.na(log_mean_HI_ratio_deng) & !is.na(log_fold_rise_chik) ~ "subclinical_deng_chik",
      !is.na(log_mean_HI_ratio_deng)   ~ "subclinical_dengue",
      TRUE  ~ "subclinical_chik"
    ),
    denv_fold_rise = 2^log_mean_HI_ratio_deng,   # human-readable
    chik_fold_rise = 2^log_fold_rise_chik
  ) %>%
  left_join(
    cebu_pivot_days_since_inf %>%
      distinct(id_patient, id_sample, date_sample, days_since_infection),
    by = c("id_patient", "id_sample")
  ) %>%
  arrange(id_patient, date_sample)
 
dplyr::count(subclinical_samples, reason)

# remove subclinical samples 
clean_data <- cebu_pivot_days_since_inf %>%
  anti_join(subclinical_samples, by = c("id_patient", "id_sample"))

 
# rows removed should be n_excluded_samples * rows_per_sample
nrow(cebu_pivot_days_since_inf) - nrow(clean_data)
nrow(subclinical_samples)


patient_pcr_mapping <- clean_data %>%
  group_by(id_patient) %>%
  summarise(
    has_zikv   = any(PCR %in% "ZIKV"),
    has_chikv  = any(PCR %in% "CHIKV"),
    has_denv   = any(PCR %in% DENV_SEROTYPES),
    # first DENV serotype specifically (not just the first PCR of any virus)
    denv_pcr   = first(c(na.omit(PCR[PCR %in% DENV_SEROTYPES]), NA_character_)),
    .groups = "drop"
  ) %>%
  mutate(
    # subclinical samples are already gone, so "no PCR+ of any virus" = negative
    target = case_when(
      has_zikv  ~ "ZIKV",
      has_denv  ~ denv_pcr,
      has_chikv ~ "CHIKV",
      TRUE      ~ "negative"
    )
  ) %>%
  dplyr::select(id_patient, target)

table(patient_pcr_mapping$target)
target_counts <- as.data.frame(table(patient_pcr_mapping$target))
colnames(target_counts) <- c("Target", "Count")

# save counts per target 
write.csv(target_counts, "Results/target_counts.csv", row.names = FALSE)


out_dir <- "Data/by_isotype_gm"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)


# geometric mean ratios 
preprocessed_by_isotype <- map(set_names(isotypes), function(iso) {
 
  raw_iso <- clean_data %>% filter(isotype == iso)
 
  log_iso <- raw_iso
  log_iso[antigen_cols] <- log10(log_iso[antigen_cols])
 
  # ratios + cross-sectional built from UNLOGGED data
  dfs <- prepare_luminex_datasets(raw_iso, patient_pcr_mapping,
                                  antigen_cols, pre_threshold = -1, post_threshold = 30)

 
  logged_ratio <- dfs$ratio
  logged_ratio[antigen_cols] <- log10(logged_ratio[antigen_cols])


  write.csv(raw_iso, sprintf("%s/raw_preprocessed_cebu_%s.csv", out_dir, iso),    row.names = FALSE)
  write.csv(log_iso, sprintf("%s/logged_preprocessed_cebu_%s.csv", out_dir, iso), row.names = FALSE)

  saveRDS(dfs$ratio,  sprintf("%s/ratio_df_%s.rds", out_dir, iso))
  saveRDS(logged_ratio,  sprintf("%s/logged_ratio_df_%s.rds", out_dir, iso))


  list(raw = raw_iso, log = log_iso,
       ratio = dfs$ratio, logged_ratio = logged_ratio)
})

saveRDS(preprocessed_by_isotype, file.path(out_dir, "preprocessed_by_isotype.rds"))
 

igg_ratio <- readRDS(file.path(out_dir, "ratio_df_IgG.rds"))
iga_ratio <- readRDS(file.path(out_dir, "ratio_df_IgA.rds"))
igm_ratio <- readRDS(file.path(out_dir, "ratio_df_IgM.rds"))













