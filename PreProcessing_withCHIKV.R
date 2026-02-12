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
  ungroup()  %>% dplyr::select(-prev_days_since_infection, -time_diff)


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
  ) %>% ungroup() %>% dplyr::select(-ever_pcr_positive)


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


antigen_cols <- c(
  "CCHF_NP", "CHIKV_E2","CHIKV_NSP123",
  "CHIKV_VLP","DENV1_DIII","DENV1_NS1",
  "DENV1_VLP","DENV2_DIII","DENV2_NS1",
  "DENV2_VLP","DENV3_DIII","DENV3_NS1",
  "DENV3_VLP", "DENV4_DIII","DENV4_NS1",
  "DENV4_VLP","JEV_E","JEV_NS1", "MAYV_E2",
  "ONNV_E2","ONNV_VLP","OROV_Gc","OROV_Gn",
  "OROV_NP_Clinisciences","OROV_NP_NA",
  "RR","RVFV","USUV_NS1","WNV_DIII",
  "WNV_NS1","YFV_E","YFV_NS1","ZIKV_NS1",
  "ZIKV_VLP","ZIKVAS_DIII", "ZIKVSU_NS1",
  "SHERPADES_CHIKV_E2", "SHERPADES_DENV1_DIII",
  "SHERPADES_DENV2_DIII","SHERPADES_DENV3_DIII",
  "SHERPADES_DENV4_DIII","SHERPADES_JEV_DIII",
  "SHERPADES_MAYV_E2","SHERPADES_RR",
  "SHERPADES_RVFV","SHERPADES_USUV_DIII",
  "SHERPADES_WNV_DIII","SHERPADES_YFV_DIII",
  "SHERPADES_ZIKV_DIII"
)


# log antigen cols for downstream analysis 
clean_HI_ratio_df[antigen_cols] <- 
  log10(clean_HI_ratio_df[antigen_cols])


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




# Define antigen groups
flavivirus_antigens <- c(
  "DENV1_DIII", "DENV1_VLP", "DENV1_NS1",
  "DENV2_DIII", "DENV2_VLP", "DENV2_NS1",
  "DENV3_DIII", "DENV3_VLP", "DENV3_NS1",
  "DENV4_DIII", "DENV4_VLP", "DENV4_NS1",
  "ZIKVAS_DIII", "ZIKV_VLP", "ZIKV_NS1", "ZIKVSU_NS1",
  "SHERPADES_DENV1_DIII", "SHERPADES_DENV2_DIII", 
  "SHERPADES_DENV3_DIII", "SHERPADES_DENV4_DIII",
  "SHERPADES_ZIKV_DIII"
)

alphavirus_antigens <- c(
  "CHIKV_E2", "CHIKV_NSP123", "CHIKV_VLP",
  "ONNV_E2", "ONNV_VLP",
  "MAYV_E2",
  "RR"
)

# Filter for IgG and prepare data # Filter for IgG and prepare data
plot_data <- clean_HI_ratio_df %>%
  filter(isotype == "IgG") %>%
  select(id_patient, id_sample, days_since_infection, PCR, 
         all_of(c(flavivirus_antigens, alphavirus_antigens))) %>%
  pivot_longer(cols = c(all_of(flavivirus_antigens), all_of(alphavirus_antigens)),
               names_to = "antigen") 



# Separate pathogen and antigen type
plot_data <- plot_data %>%
  mutate(
    pathogen = case_when(
      str_detect(antigen, "^SHERPADES_") ~ str_extract(antigen, "(?<=SHERPADES_)[^_]+"),
      str_detect(antigen, "^ZIKVSU") ~ "ZIKV",
      str_detect(antigen, "^ZIKVAS") ~ "ZIKV",
      TRUE ~ str_extract(antigen, "^[^_]+")
    ),
    antigen_type = case_when(
      str_detect(antigen, "^SHERPADES_") ~ "SHERPADES_DIII",
      TRUE ~ str_extract(antigen, "[^_]+$")
    ),
    virus_family = case_when(
      pathogen %in% c("DENV1", "DENV2", "DENV3", "DENV4", "ZIKV") ~ "Flavivirus",
      pathogen %in% c("CHIKV", "ONNV", "MAYV", "RRV") ~ "Alphavirus"
    )
  )

# Create time categories for x-axis (like screenshot)
plot_data <- plot_data %>%
  mutate(
    time_category = case_when(
      days_since_infection < -30 ~ "pre-infection",
      days_since_infection >= -30 & days_since_infection < 0 ~ as.character(days_since_infection),
      days_since_infection >= 0 & days_since_infection <= 30 ~ as.character(days_since_infection),
      days_since_infection > 30 ~ ">30"
    ),
    # Create numeric position for plotting
    x_position = case_when(
      days_since_infection < -30 ~ -35,
      days_since_infection >= -30 & days_since_infection < 0 ~ days_since_infection,
      days_since_infection >= 0 & days_since_infection <= 30 ~ days_since_infection,
      days_since_infection > 30 ~ 35
    )
  )

# Create homologous/heterologous indicator
plot_data <- plot_data %>%
  mutate(
    response_type = case_when(
      pathogen == PCR ~ "Homologous",
      TRUE ~ "Heterologous"
    )
  )

plot_data <- plot_data %>%
  group_by(id_patient) %>%
  fill(PCR, .direction = "downup") %>%
  ungroup()

plot_data <- plot_data %>%
  mutate(
    color_group = if_else(pathogen == PCR, PCR, "Other")
  )

# Color palette
pcr_colors <- c(
  "DENV1" = "#012b48",
  "DENV2" = "#0396f8", 
  "DENV3" = "#693bf1",
  "DENV4" = "#de5a7b",
  "ZIKV" = "#21737c",
  "Other" = "grey90"
)

plot_antigen_dynamics <- function(data, pathogen_subset = NULL) {
  
  if (!is.null(pathogen_subset)) {
    data <- data %>% filter(pathogen %in% pathogen_subset)
  }
  
  ggplot(data, aes(x = x_position, y = value, group = id_patient)) +
    
    geom_line(aes(color = color_group),
              alpha = 0.7, linewidth = 0.5) +
    
    geom_point(aes(color = color_group),
               alpha = 0.8, size = 1.3) +
    
    facet_grid(antigen_type ~ pathogen, scales = "free_y") +
    
    scale_color_manual(values = pcr_colors, name = "PCR Confirmed") +
    
    scale_x_continuous(
      breaks = c(-35, seq(-20, 30, 10), 35),
      labels = c("pre-infection", seq(-20, 30, 10), ">30")
    ) +
    
    geom_vline(xintercept = 0, linetype = "dashed", color = "black", alpha = 0.3) +
    
    labs(x = "Days since PCR+ infection",
         y = "Antibody Titre (log10)",
         title = "IgG Antibody Dynamics Around Infection") +
    
    theme_bw() +
    theme(
       # Facet strips
        strip.text = element_text(size = 12),
        
        # Axis text
        axis.text = element_text(size = 11),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_text(size = 13),
        
        # Legend
        legend.position = "bottom",
        legend.text = element_text(size = 11),
        legend.title = element_text(size = 11),
        legend.direction = "horizontal",
        legend.margin = margin(20, 0, 0, 0),
        
        # Panel spacing (important for faceted plots)
        panel.spacing = unit(1.2, "lines"),
        
        # Grid styling
        panel.grid.minor = element_blank(),
        
        # Outer box
        plot.background = element_rect(
          fill = "white", 
          color = "black", 
          linewidth = 0.5
        ),
        
        plot.margin = margin(10, 5, 15, 5)
    )
}


# Create Flavivirus plot
flavi_plot <- plot_antigen_dynamics(
  data = plot_data %>% filter(pathogen %in% c("DENV1", "DENV2", "DENV3", "DENV4", "ZIKV"))
)
print(flavi_plot)
# Create Alphavirus plot  
alpha_plot <- plot_antigen_dynamics(
  data = plot_data %>% filter(pathogen %in% c("CHIKV", "ONNV", "MAYV", "RRV"))
)
print(alpha_plot)
# Save plots
ggsave("/Users/ap2488/Desktop/supervised_learning_flavi/FinalLuminexClassification/flavivirus_IgG_dynamics.png", 
flavi_plot, width = 16, height = 10)

# save file 
write.csv(clean_HI_ratio_df,
"/Users/ap2488/Desktop/supervised_learning_flavi/FinalLuminexClassification/clean_HI_ratio_df.csv")
