library(grid)
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(patchwork)
library(ggh4x)

# import ratio df 
logged_preprocessed_cebu_data <- read.csv('Results/logged_preprocessed_cebu_data.csv')
colnames(logged_preprocessed_cebu_data)
head(logged_preprocessed_cebu_data)

# Define antigen groups
dengue_zika_antigens <- c(
  "DENV1_DIII", "DENV1_VLP", "DENV1_NS1","SHERPADES_DENV1_DIII",
  "DENV2_DIII", "DENV2_VLP", "DENV2_NS1", "SHERPADES_DENV2_DIII", 
  "DENV3_DIII", "DENV3_VLP", "DENV3_NS1", "SHERPADES_DENV3_DIII",
  "DENV4_DIII", "DENV4_VLP", "DENV4_NS1", "SHERPADES_DENV4_DIII",
  "ZIKVAS_DIII", "ZIKV_VLP", "ZIKV_NS1", "ZIKVSU_NS1", "SHERPADES_ZIKV_DIII"
)

chik_onnv_mayv_antigens <- c(
  "CHIKV_E2", "CHIKV_NSP123", "CHIKV_VLP","SHERPADES_CHIKV_E2",
  "ONNV_E2", "ONNV_VLP",
  "MAYV_E2" , "SHERPADES_MAYV_E2"
)

pcr_colours <- c(
  "DENV1" = "#012b48",
  "DENV2" = "#0396f8",
  "DENV3" = "#4525a4",
  "DENV4" = "#de5a7b",
  "ZIKV"  = "#21737c",
  "CHIKV" = "#7d2102",
  "Other" = "#cbcbcb"
)



# --- Prepare data for a given isotype 
prepare_antibody_data <- function(data, isotype, antigens = c(dengue_zika_antigens, chik_onnv_mayv_antigens)) {
  
  data %>%
    filter(isotype == !!isotype) %>%
    dplyr::select(id_patient, id_sample, days_since_infection, PCR, all_of(antigens)) %>%
    pivot_longer(
      cols = all_of(antigens),
      names_to = "antigen",
      values_to = "value"
    ) %>%
    mutate(
      pathogen = case_when(
        str_detect(antigen, "^SHERPADES_") ~ str_extract(antigen, "(?<=SHERPADES_)[^_]+"),
        str_detect(antigen, "^ZIKVSU") ~ "ZIKV",
        str_detect(antigen, "^ZIKVAS") ~ "ZIKV",
        TRUE ~ str_extract(antigen, "^[^_]+")
      ),
      antigen_type = case_when(
        str_detect(antigen, "^SHERPADES_") ~ "DIII-SHERPADES",
        TRUE ~ str_extract(antigen, "[^_]+$")
      ),
      virus_family = case_when(
        pathogen %in% c("DENV1", "DENV2", "DENV3", "DENV4", "ZIKV") ~ "Flavivirus",
        pathogen %in% c("CHIKV", "ONNV", "MAYV", "RRV") ~ "Alphavirus",
        TRUE ~ NA_character_
      ),
      days_since_infection = as.numeric(days_since_infection)
    ) %>%
    group_by(id_patient) %>%
    fill(PCR, .direction = "downup") %>%
    ungroup() %>%
    mutate(
      time_bin = case_when(
        days_since_infection < 0 ~ "pre",
        days_since_infection > 30 ~ "post",
        TRUE ~ as.character(days_since_infection)
      ),
      x_position = case_when(
        days_since_infection < 0 ~ -30,
        days_since_infection > 30 ~ 30,
        TRUE ~ days_since_infection
      ),
      color_group = if_else(pathogen == PCR, PCR, "Other")
    ) %>%
    group_by(id_patient, antigen, time_bin, x_position) %>%
    summarise(
      value = mean(value, na.rm = TRUE),
      PCR = first(PCR),
      color_group = first(color_group),
      antigen_type = first(antigen_type),
      pathogen = first(pathogen),
      virus_family = first(virus_family),
      .groups = "drop"
    ) %>%
    arrange(color_group == "Other", color_group)
}

# --- single plot for one isotype + virus family 
plot_antibody_dynamics <- function(data,
                                   isotype) {
  
  
  all_pathogens    <- c("DENV1", "DENV2", "DENV3", "DENV4", "ZIKV", "CHIKV", "ONNV", "MAYV")
  confirmed        <- c("DENV1", "DENV2", "DENV3", "DENV4", "ZIKV", "CHIKV")
  alpha_pathogens <- c("CHIKV", "ONNV", "MAYV", "RRV")
  flavi_pathogens <- c("DENV1", "DENV2", "DENV3", "DENV4", "ZIKV")

  y_label <- if (isotype == "avidity") "Avidity" else "Antibody Titre (log2)"

  foreground_data <- data %>%
  filter(
    PCR %in% confirmed,
    pathogen == PCR
  ) %>%
  mutate(
    facet_pathogen = PCR,
    color_group = PCR
  )

  background_data <- data %>%
    filter(
      PCR %in% confirmed,
      pathogen %in% flavi_pathogens,
      pathogen != PCR
    ) %>%
    mutate(
      facet_pathogen = PCR,
      color_group = "Other"
    )

 

  # Foreground data: only confirmed pathogens
  facet_data <- data %>%
    filter(pathogen %in% confirmed) %>%
    mutate(
      facet_pathogen = pathogen
    )

  foreground_data <- facet_data %>%
    filter(color_group != "Other")

  ggplot(facet_data, aes(x = x_position, y = value, group = interaction(id_patient, pathogen, antigen))) +
    geom_line(
      data = background_data,
      color = "grey90",
      alpha = 0.7,
      linewidth = 0.5
    ) +
    geom_point(
      data = background_data,
      color = "grey90",
      alpha = 0.7,
      size = 2
    ) +
    geom_line(
      data = foreground_data,
      aes(color = color_group),
      alpha = 0.8,
      linewidth = 0.7
    ) +
    geom_point(
      data = foreground_data,
      aes(color = color_group),
      alpha = 0.8,
      size = 2
    ) +
    facet_grid(antigen_type ~ facet_pathogen, scales = "free_y", drop = TRUE) +
    scale_color_manual(values = pcr_colours, name = "PCR Confirmed") +
    scale_x_continuous(
      breaks = c(-30, 0, 15, 30),
      labels = c("Pre-infection", "0", "15", "> 30 days \n since infection")
    ) +
    labs(
      x = "Days since PCR+ve infection",
      y = y_label
    ) +
    theme_bw() +
    theme(
      strip.text = element_text(size = 20),
      strip.background = element_rect(fill = "#ffffff"),
      axis.text = element_text(size = 20),
      axis.text.x = element_text(size = 20),
      axis.text.y = element_text(size = 20),
      axis.title = element_text(size = 20),
      legend.position = "bottom",
      legend.text = element_text(size = 20),
      legend.title = element_text(size = 20),
      legend.direction = "horizontal",
      legend.margin = margin(0, 0, 0, 0),
      panel.spacing.x = unit(1, "lines"),
      panel.spacing.y = unit(1, "lines"),
      panel.grid.minor = element_blank(),
      plot.background = element_rect(fill = "white", color = "black", linewidth = 0.5)
    )
}



igg_data <- prepare_antibody_data(logged_preprocessed_cebu_data, "IgG", antigens = c(dengue_zika_antigens, chik_onnv_mayv_antigens))
igm_data <- prepare_antibody_data(logged_preprocessed_cebu_data, "IgM", antigens = c(dengue_zika_antigens, chik_onnv_mayv_antigens))
iga_data <- prepare_antibody_data(logged_preprocessed_cebu_data, "IgA", antigens = c(dengue_zika_antigens, chik_onnv_mayv_antigens))
avidity_data <- prepare_antibody_data(logged_preprocessed_cebu_data, "avidity", antigens = c(dengue_zika_antigens, chik_onnv_mayv_antigens))


# Plot separately
p_igg <- plot_antibody_dynamics(igg_data, "IgG") +
  force_panelsizes(rows = unit(3, "cm"), cols = unit(3, "cm"))

quartz()
print(p_igg)

# save IgG 
ggsave("Results/dynamics_IgG_.png", 
p_igg, width = 25, height = 12)

# save IgG 
ggsave("Results/dynamics_alphavirus_IgG_.png", 
p_alpha_igg, width = 8, height = 10)

#save combined 
ggsave("Results/Figure2.png", 
p_flavi_igg + p_alpha_igg, width = 8, height = 10)

# other isotype plots 
igm_flavi <- plot_antibody_dynamics(igm_data, "IgM", "Flavivirus")
iga_flavi <- plot_antibody_dynamics(iga_data, "IgA", "Flavivirus")
avidity_flavi <- plot_antibody_dynamics(avidity_data, "Avidity", "Flavivirus")


print(igm_flavi)
print(iga_flavi)
print(avidity_flavi)

