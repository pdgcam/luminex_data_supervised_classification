# --- Script for preprocessing data and visualising tire dynamics 
# --- Preprocessing for cross-sectional samples 
# mostly the same as data_preprocessing.R
# but here we need to make sure to remove potetinal subclincals at the start of the study 

# --- also here we are using geometric mean ratio across all four serotypes of DENV 
#  Unless otherwise specified, we classified an interval as a no-infection event if this ratio was 1 or smaller. 



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
library(here)
library(exactextractr)
library(raster)
library(readxl)
library(knitr)
library(patchwork)
library(tidyverse)
library(standardize)



DENV_SEROTYPES <- c("DENV1", "DENV2", "DENV3", "DENV4")
isotypes     <- c("IgG", "IgA", "IgM")
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




# --- DENV Seronegative at baseline 

hai_cols <- c("HAI_DENV1", "HAI_DENV2", "HAI_DENV3", "HAI_DENV4")

# one row per sample, with HAI status attached
sample_hai <- cebu_pivot_days_since_inf %>%
  distinct(id_patient, .keep_all = TRUE) %>%
  mutate(across(all_of(hai_cols), ~ trimws(as.character(.x))),
         hai_neg = if_all(all_of(hai_cols), ~ .x == "<10"))

# patients whose EARLIEST sample is negative across all four serotypes
flavi_baseline_neg_ids <- sample_hai %>%
  group_by(id_patient) %>%
  arrange(date_sample, .by_group = TRUE) %>%
  slice(1) %>%
  ungroup() %>%
  filter(hai_neg) %>%
  pull(id_patient)

length(flavi_baseline_neg_ids)

# ----- CHIKV Seronegative at baseline 
# one row per sample, with CHIKV PRNT status attached
sample_prnt <- cebu_pivot_days_since_inf %>%
  distinct(id_sample, .keep_all = TRUE) %>%
  mutate(PRNT_CHIKV = trimws(as.character(PRNT_CHIKV)),
         prnt_neg = PRNT_CHIKV == "<10")

# patients whose EARLIEST sample is CHIKV PRNT negative
alpha_baseline_neg_ids <- sample_prnt %>%
  group_by(id_patient) %>%
  arrange(days_since_infection, date_sample, .by_group = TRUE) %>%
  slice(1) %>%
  ungroup() %>%
  filter(prnt_neg) %>%
  pull(id_patient)


# --- patient mapping (PCRs)
patient_pcr_mapping <- cebu_pivot_days_since_inf %>%
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


cebu_pivot_days_since_inf %>%
  filter(PCR == "DENV1") %>%
  distinct(id_patient, id_sample) %>%
  dplyr::count(id_patient, name = "n_samples") %>%
  dplyr::count(n_samples)


# two IDs only have a single infection - so DENV1 goes down from 12 to 10 PCR confirmed 
singletons <- cebu_pivot_days_since_inf %>%
  filter(PCR == "DENV1") %>%
  distinct(id_patient, id_sample) %>%
  dplyr::count(id_patient) %>%
  filter(n == 1) %>%
  pull(id_patient) 



# --- build cross-sectional datasets 
# --- positives == PCR confirmed infections 
# --- negatives == baseline neg (flavi or alpha) 

build_cross_sectional_datasets <- function(dat, antigen_cols, pos_targets, neg_ids, iso,
                               pcr_mapping = patient_pcr_mapping) {

  d <- dat %>% filter(isotype == iso)

  pos_map <- pcr_mapping %>% filter(target %in% pos_targets)

  # NEGATIVE rows: baseline draw of the baseline-negative patients.
  # Keep only draws that are genuinely pre-infection (or from never-infected patients).
  neg <- d %>%
    filter(id_patient %in% neg_ids,
           is.na(days_since_infection) | days_since_infection < 0) %>%
    group_by(id_patient) %>%
    slice_min(date_sample, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    mutate(target = "negative")


  # POSITIVE rows: last post-infection draw of the PCR-confirmed patients
  pos <- d %>%
    filter(id_patient %in% pos_map$id_patient,
           !is.na(days_since_infection), days_since_infection > 0) %>%
    group_by(id_patient) %>%
    slice_max(days_since_infection, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    left_join(pos_map, by = "id_patient")

  bind_rows(neg, pos) %>%
  dplyr::select(id_patient, id_sample, days_since_infection, target, all_of(antigen_cols)) %>%
  tidyr::drop_na(all_of(antigen_cols)) %>%     # instead of na.omit()
  as.data.frame()
}

flavi_IgG <- build_cross_sectional_datasets(cebu_pivot_days_since_inf, flavi_antigens,
                                c(DENV_SEROTYPES, "ZIKV"), flavi_baseline_neg_ids, "IgG")
flavi_IgA <- build_cross_sectional_datasets(cebu_pivot_days_since_inf, flavi_antigens,
                                c(DENV_SEROTYPES, "ZIKV"), flavi_baseline_neg_ids, "IgA")
flavi_IgM <- build_cross_sectional_datasets(cebu_pivot_days_since_inf, flavi_antigens,
                                c(DENV_SEROTYPES, "ZIKV"), flavi_baseline_neg_ids, "IgM")

alpha_IgG <- build_cross_sectional_datasets(cebu_pivot_days_since_inf, alpha_antigens,
                                "CHIKV", alpha_baseline_neg_ids, "IgG")
alpha_IgA <- build_cross_sectional_datasets(cebu_pivot_days_since_inf, alpha_antigens,
                                "CHIKV", alpha_baseline_neg_ids, "IgA")
alpha_IgM <- build_cross_sectional_datasets(cebu_pivot_days_since_inf, alpha_antigens,
                                "CHIKV", alpha_baseline_neg_ids, "IgM")


nrow(flavi_IgG)
View(flavi_IgG)

#log10 
log_antigens <- function(df, antigen_cols) {
  df[antigen_cols] <- log10(df[antigen_cols])
  df
}

flavi_IgG_log <- log_antigens(flavi_IgG, flavi_antigens)
flavi_IgA_log <- log_antigens(flavi_IgA, flavi_antigens)
flavi_IgM_log <- log_antigens(flavi_IgM, flavi_antigens)

alpha_IgG_log <- log_antigens(alpha_IgG, alpha_antigens)
alpha_IgA_log <- log_antigens(alpha_IgA, alpha_antigens)
alpha_IgM_log <- log_antigens(alpha_IgM, alpha_antigens)



saveRDS(flavi_IgG_log, "Data/by_isotype/flavi_cross_sectional_logged_IgG.rds")
saveRDS(flavi_IgA_log, "Data/by_isotype/flavi_cross_sectional_logged_IgA.rds")
saveRDS(flavi_IgM_log, "Data/by_isotype/flavi_cross_sectional_logged_IgM.rds")
saveRDS(alpha_IgG_log, "Data/by_isotype/alpha_cross_sectional_logged_IgG.rds")
saveRDS(alpha_IgA_log, "Data/by_isotype/alpha_cross_sectional_logged_IgA.rds")
saveRDS(alpha_IgM_log, "Data/by_isotype/alpha_cross_sectional_logged_IgM.rds")






flavi_baseline_neg <- cebu_pivot_days_since_inf %>%
  filter(id_patient %in% baseline_neg_ids,
         isotype %in% c("IgG", "IgA", "IgM")) %>%
  left_join(dplyr::select(sample_hai, id_sample, hai_neg), by = "id_sample") %>%
  pivot_longer(all_of(flavi_antigens),
               names_to = "antigen", values_to = "RAU") %>%
  filter(!is.na(RAU)) %>%
  mutate(
    isotype = factor(isotype, levels = c("IgG", "IgA", "IgM")),
    target = case_when(
      grepl("^SHERPADES_", antigen) ~ "DIII (SHERPADES)",
      grepl("_DIII$", antigen)      ~ "DIII",
      grepl("_NS1$", antigen)       ~ "NS1",
      grepl("_VLP$", antigen)       ~ "VLP",
      grepl("_E$|_E2$", antigen)    ~ "E",
      TRUE                          ~ "other"
    )
  ) %>%
  arrange(id_patient, antigen, days_since_infection)

View(flavi_baseline_neg)

plot_flavi_traj <- function(dat, iso, ids = NULL) {

  stopifnot(all(c("days_since_infection", "RAU", "target",
                  "antigen", "isotype", "id_patient") %in% names(dat)))

  d <- dplyr::filter(dat, isotype == iso)
  if (!is.null(ids)) d <- dplyr::filter(d, id_patient %in% ids)
  if (nrow(d) == 0) stop("No rows for isotype ", iso)

  n_pat <- dplyr::n_distinct(d$id_patient)

  ggplot(d, aes(.data$days_since_infection, .data$RAU,
                colour = .data$target, group = .data$antigen)) +
    geom_line(linewidth = 0.3, alpha = 0.7) +
    geom_point(size = 0.7) +
    facet_wrap(~ id_patient, ncol = min(8, n_pat)) +
    scale_y_log10() +
    labs(x = "Days since infection", y = paste0("Flavivirus ", iso, " (RAU)"),
         title = paste0(iso, " — HAI <10 across DENV1-4 at baseline (n = ", n_pat, ")")) +
    theme_minimal(base_size = 8) +
    theme(panel.grid.minor = element_blank(),
          strip.text = element_text(face = "bold", size = 7),
          legend.position = "bottom")
}

quartz()
plot_flavi_traj(flavi_baseline_neg, "IgG", ids = baseline_neg_ids)