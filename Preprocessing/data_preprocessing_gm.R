
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
library(readxl)


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
prepare_ratio_datasets <- function(raw_data, patient_mapping, antigen_cols, pre_threshold = -1, post_threshold = 30) {
  
  # Validate that antigen columns exist in the data
  missing_cols <- setdiff(antigen_cols, names(raw_data))
  if (length(missing_cols) > 0) {
    stop("Missing columns in data: ", paste(missing_cols, collapse = ", "))
  }
    # removea any samples between day 0 and 30 
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
#DENV_THRESHOLD <- log2(1.6)   # >=1.6 rise in mean log2 HAI between first and last sample
DENV_THRESHOLD <- 1.6 # is it 1.6 or log2(1.6) - need to check??!!
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

denv_HAIs_no_denv_PCR  <- sample_per_patient(no_denv_pcr)
chik_PRNTs_no_chikv_PCR <- sample_per_patient(no_chikv_pcr)


# DENV Subclinicals  - look at no_denv_pcr subset (patients with no symptomatic DENV infection during study)
denv_subclinicals_df <- denv_HAIs_no_denv_PCR %>%
  arrange(id_patient, date_sample) %>%
  group_by(id_patient) %>%
  mutate(
    log_HI_DENV1 = if_else(HAI_DENV1 > 0, log_transform(HAI_DENV1), NA_real_),
    log_HI_DENV2 = if_else(HAI_DENV2 > 0, log_transform(HAI_DENV2), NA_real_),
    log_HI_DENV3 = if_else(HAI_DENV3 > 0, log_transform(HAI_DENV3), NA_real_),
    log_HI_DENV4 = if_else(HAI_DENV4 > 0, log_transform(HAI_DENV4), NA_real_),

    mean_log_HI = rowMeans(
      cbind(log_HI_DENV1, log_HI_DENV2, log_HI_DENV3, log_HI_DENV4),
      na.rm = TRUE
    ),
    # all-NA rows come back as NaN, not NA - make that explicit
    mean_log_HI = if_else(is.nan(mean_log_HI), NA_real_, mean_log_HI),
    has_HI = !is.na(mean_log_HI),

    first_log_HI = dplyr::first(mean_log_HI[has_HI], default = NA_real_),
    last_log_HI  = dplyr::last(mean_log_HI[has_HI],  default = NA_real_),

    enrollment_seroneg    = first_log_HI < 1,
    enrollment_serostatus = if_else(enrollment_seroneg, "seronegative", "seropositive"),

    first_HI_date = dplyr::first(date_sample[has_HI], default = as.Date(NA)),
    last_HI_date  = dplyr::last(date_sample[has_HI],  default = as.Date(NA)),
    days_gap      = as.numeric(difftime(last_HI_date, first_HI_date, units = "days")),

    log_mean_HI_ratio_deng = last_log_HI - first_log_HI,

    denv_status = case_when(
      is.na(log_mean_HI_ratio_deng) | is.na(days_gap) | days_gap == 0 ~ "insufficient_data",
      enrollment_seroneg & log_mean_HI_ratio_deng <  DENV_THRESHOLD   ~ "seronegative_no_rise",
      enrollment_seroneg & log_mean_HI_ratio_deng >= DENV_THRESHOLD   ~ "seronegative_start_subclinical",
      log_mean_HI_ratio_deng >= DENV_THRESHOLD                        ~ "seropositive_start_subclinical",
      .default = "seropositive_no_rise"
    )
  ) %>%
  ungroup() %>%
  dplyr::select(id_patient, id_sample, mean_log_HI, date_sample, days_gap,
                days_since_infection, PCR, enrollment_serostatus,
                log_mean_HI_ratio_deng, denv_status)


# CHIKV subclincicals
chik_baseline <- chik_PRNTs_no_chikv_PCR %>%
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
 
chik_subclinicals_df <- chik_PRNTs_no_chikv_PCR %>%
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




DENV_SUBCLINICAL <- c("seronegative_start_subclinical", "seropositive_start_subclinical")
FLAVI_TARGETS <- c(DENV_SEROTYPES, "ZIKV")

# --- patient-level subclinical flags (one row per patient) ---
denv_subclinical_ids <- denv_subclinicals_df %>%
  filter(denv_status %in% DENV_SUBCLINICAL) %>%
  distinct(id_patient) %>%
  pull(id_patient)

chik_subclinical_ids <- chik_subclinicals_df %>%
  filter(chik_subclinical) %>%
  distinct(id_patient) %>%
  pull(id_patient)

length(denv_subclinical_ids)
length(chik_subclinical_ids)
length(intersect(denv_subclinical_ids, chik_subclinical_ids))


patient_pcr_mapping <- cebu_pivot_days_since_inf %>%
  group_by(id_patient) %>%
  summarise(
    has_zikv  = any(PCR %in% "ZIKV"),
    has_chikv = any(PCR %in% "CHIKV"),
    has_denv  = any(PCR %in% DENV_SEROTYPES),
    denv_pcr  = first(c(na.omit(PCR[PCR %in% DENV_SEROTYPES]), NA_character_)),
    .groups = "drop"
  ) %>%
  mutate(
    target = case_when(
      has_zikv  ~ "ZIKV",
      has_denv  ~ denv_pcr,
      has_chikv ~ "CHIKV",
      TRUE      ~ "negative"
    )
  ) %>%
  dplyr::select(id_patient, target)

print(table(patient_pcr_mapping$target))

target_counts <- patient_pcr_mapping %>%
  dplyr::count(target) %>%
  dplyr::mutate(
    target = factor(target, levels = c("DENV1", "DENV2", "DENV3", "DENV4", "ZIKV", "CHIKV", "negative"))
  ) %>%
  arrange(target)
print(target_counts)

# save counts per target 
write.csv(target_counts, "Results/target_counts.csv", row.names = FALSE)



denv_data <- cebu_pivot_days_since_inf %>%
  left_join(patient_pcr_mapping, by = "id_patient") %>%
  filter(target %in% DENV_SEROTYPES | !(id_patient %in% denv_subclinical_ids))

chik_data <- cebu_pivot_days_since_inf %>%
  left_join(patient_pcr_mapping, by = "id_patient") %>%
  filter(target == "CHIKV" | !(id_patient %in% chik_subclinical_ids))

denv_vs_chik_data <- cebu_pivot_days_since_inf %>%
  left_join(patient_pcr_mapping, by = "id_patient") %>%
  filter(
    (target %in% DENV_SEROTYPES | !(id_patient %in% denv_subclinical_ids)) &
    (target == "CHIKV"          | !(id_patient %in% chik_subclinical_ids))
  )



purrr::imap(
  list(denv = denv_data, chik = chik_data, denv_vs_chik = denv_vs_chik_data),
  ~ tibble(set = .y, n_rows = nrow(.x), n_patients = n_distinct(.x$id_patient))
) %>% bind_rows()

denv_data %>% dplyr::distinct(id_patient, target) %>% dplyr::count(target)
chik_data %>% dplyr::distinct(id_patient, target) %>% dplyr::count(target)
denv_vs_chik_data %>% dplyr::distinct(id_patient, target) %>% dplyr::count(target)



base_out_dir <- "Data/by_isotype_gm"

analysis_sets <- list(
  denv         = denv_data,
  chik         = chik_data,
  denv_vs_chik = denv_vs_chik_data
)

preprocessed_by_set <- purrr::imap(analysis_sets, function(set_data, set_name) {

  out_dir <- file.path(base_out_dir, set_name)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  by_isotype <- purrr::map(rlang::set_names(isotypes), function(iso) {

    raw_iso <- set_data %>% dplyr::filter(isotype == iso)

    log_iso <- raw_iso
    log_iso[antigen_cols] <- log10(log_iso[antigen_cols])

    dfs <- prepare_ratio_datasets(
      raw_iso, patient_pcr_mapping, antigen_cols,
      pre_threshold = -1, post_threshold = 30
    )

    logged_ratio <- dfs$ratio
    logged_ratio[antigen_cols] <- log10(logged_ratio[antigen_cols])

    write.csv(raw_iso, sprintf("%s/raw_preprocessed_cebu_%s.csv", out_dir, iso),
              row.names = FALSE)
    write.csv(log_iso, sprintf("%s/logged_preprocessed_cebu_%s.csv", out_dir, iso),
              row.names = FALSE)

    saveRDS(dfs$ratio,    sprintf("%s/ratio_df_%s.rds", out_dir, iso))
    saveRDS(logged_ratio, sprintf("%s/logged_ratio_df_%s.rds", out_dir, iso))

    list(raw = raw_iso, log = log_iso,
         ratio = dfs$ratio, logged_ratio = logged_ratio)
  })

  saveRDS(by_isotype, file.path(out_dir, "preprocessed_by_isotype.rds"))
  by_isotype
})

saveRDS(preprocessed_by_set, file.path(base_out_dir, "preprocessed_by_set.rds"))
