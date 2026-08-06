library(grid)
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(patchwork)
library(ggh4x)

# import ratio df -  log10(post / pre)
igg_data <- read.csv("Data/by_isotype/logged_preprocessed_cebu_IgG.csv")
iga_data <- read.csv("Data/by_isotype/logged_preprocessed_cebu_IgA.csv")
igm_data <- read.csv("Data/by_isotype/logged_preprocessed_cebu_IgM.csv")
avidity_data <- read.csv("Data/by_isotype/logged_preprocessed_cebu_avidity.csv")


sentivity_igg_data <- read.csv()


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
prepare_antibody_data <- function(data, 
                                  antigens = c(dengue_zika_antigens, chik_onnv_mayv_antigens), 
                                  day0_as_post = FALSE) {
  data %>%
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
        day0_as_post & days_since_infection == 0 ~ "post", # for sensitivity analysis - acute (day 0) sample in post infection
        TRUE ~ as.character(days_since_infection)
      ),
      x_position = case_when(
        days_since_infection < 0 ~ -30,
        days_since_infection > 30 ~ 30,
        day0_as_post & days_since_infection == 0 ~ 30,
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
                                   isotype,
                                   subtitle = NULL,
                                   # explicit pairing of structurally-equivalent antigen types into shared rows
                                   antigen_class_map = c(
                                     "DIII"           = "DIII / E2",     # flavi Domain III <-> alpha E2 (structural)
                                     "E2"             = "DIII / E2",
                                     "NS1"            = "NS1 / NSP123",  # flavi NS1 <-> alpha NSP123 (non-structural)
                                     "NSP123"         = "NS1 / NSP123",
                                     "DIII-SHERPADES" = "SHERPADES-DIII",   # shortened so the strip label doesn't clip
                                     "SHERPA-DIII"    = "SHERPADES-DIII",
                                     "VLP"            = "VLP"
                                   ),
                                   row_order = c("DIII / E2", "SHERPADES-DIII", "NS1 / NSP123", "VLP"),
                                   col_order = c("CHIKV", "DENV1", "DENV2", "DENV3", "DENV4", "ZIKV")) {

  confirmed       <- c("DENV1", "DENV2", "DENV3", "DENV4", "ZIKV", "CHIKV")
  flavi_pathogens <- c("DENV1", "DENV2", "DENV3", "DENV4", "ZIKV")
  alpha_pathogens <- c("CHIKV", "ONNV", "MAYV", "RRV")

  y_label <- if (isotype == "avidity") "Avidity" else paste0(isotype, " Antibody Titre (log2)")

  facet_data <- data %>%
    filter(PCR %in% confirmed) %>%
    mutate(
      # facet is the patient's confirmed infection, NOT the antigen's own target
      facet_pathogen = PCR,
      facet_family = case_when(
        facet_pathogen %in% flavi_pathogens ~ "Flavivirus",
        facet_pathogen %in% alpha_pathogens ~ "Alphavirus",
        TRUE ~ NA_character_
      )
    ) %>%
    # keep only antigens from the same family as the facet
    filter(virus_family == facet_family) %>%
    mutate(
      antigen_class = unname(ifelse(
        antigen_type %in% names(antigen_class_map),
        antigen_class_map[antigen_type],
        antigen_type                       # unmapped types fall through unchanged
      )),
      antigen_class = factor(
        antigen_class,
        levels = intersect(c(row_order, sort(unique(antigen_class))), unique(antigen_class))
      ),
      facet_pathogen = factor(
        facet_pathogen,
        levels = intersect(col_order, unique(facet_pathogen))
      )
    )

  # coloured: antigen matches this facet's confirmed pathogen
  foreground_data <- facet_data %>% filter(pathogen == facet_pathogen)
  # grey: all other same-family antigens
  background_data <- facet_data %>% filter(pathogen != facet_pathogen)

  ggplot(facet_data, aes(x = x_position, y = value,
                         group = interaction(id_patient, pathogen, antigen))) +
    geom_line(data = background_data,  color = "grey80", alpha = 0.5, linewidth = 0.5) +
    geom_point(data = background_data, color = "grey80", alpha = 0.5, size = 1.5) +
    geom_line(data = foreground_data,  aes(color = color_group), alpha = 0.8, linewidth = 0.7) +
    geom_point(data = foreground_data, aes(color = color_group), alpha = 0.8, size = 1.8) +
    facet_grid(
      antigen_class ~ facet_pathogen,
      scales = "free_y", drop = TRUE,
      labeller = labeller(antigen_class = label_wrap_gen(12))
    ) +
    scale_color_manual(values = pcr_colours, name = "PCR Confirmed") +
    scale_x_continuous(
      breaks = c(-30, 0, 30),
      labels = c("Pre", "0", ">30")     # short so repeated axes don't overlap
    ) +
    labs(x = "Days since PCR+ve infection", y = y_label, subtitle = subtitle) +
    theme_bw() +
    theme(
      strip.text.x     = element_text(size = 16),
      strip.text.y     = element_text(size = 12, angle = 0, hjust = 0),  # horizontal = no rotation clipping
      strip.background = element_rect(fill = "#ffffff"),
      axis.text        = element_text(size = 12),
      axis.title       = element_text(size = 16),
      legend.position  = "bottom",
      legend.text      = element_text(size = 14),
      legend.title     = element_text(size = 14),
      legend.direction = "horizontal",
      legend.margin    = margin(0, 0, 0, 0),
      panel.spacing.x  = unit(1.5, "lines"),
      panel.spacing.y  = unit(0.8, "lines"),
      panel.grid.minor = element_blank(),
      plot.subtitle = element_text(size = 14, face = "italic"),
      plot.margin      = margin(t = 10, r = 25, b = 10, l = 10),
      plot.background  = element_rect(fill = "white", color = "black", linewidth = 0.5)
    )
}


prepared_igg_data <- prepare_antibody_data(igg_data, antigens = c(dengue_zika_antigens, chik_onnv_mayv_antigens))
prepared_igm_data <- prepare_antibody_data(igm_data, antigens = c(dengue_zika_antigens, chik_onnv_mayv_antigens))
prepared_iga_data <- prepare_antibody_data(iga_data, antigens = c(dengue_zika_antigens, chik_onnv_mayv_antigens))
prepared_avidity_data <- prepare_antibody_data(avidity_data, antigens = c(dengue_zika_antigens, chik_onnv_mayv_antigens))



# Plot separately
p_igg <- plot_antibody_dynamics(prepared_igg_data, "IgG") +
  force_panelsizes(rows = unit(3, "cm"), cols = unit(3, "cm"))

p_igm <- plot_antibody_dynamics(prepared_igm_data, "IgM") +
  force_panelsizes(rows = unit(3, "cm"), cols = unit(3, "cm"))

p_iga <- plot_antibody_dynamics(prepared_iga_data, "IgA") +
  force_panelsizes(rows = unit(3, "cm"), cols = unit(3, "cm"))

p_avidity <- plot_antibody_dynamics(prepared_avidity_data, "avidity") +
  force_panelsizes(rows = unit(3, "cm"), cols = unit(3, "cm"))

print(p_igg)
print(p_igm)
print(p_iga)
print(p_avidity)

# --- create output folder
dir.create("Results/Fig2", recursive = TRUE, showWarnings = FALSE)

# --- save all four
ggsave("Results/Fig2/dynamics_IgG.png", p_igg,  width = 12, height = 10)
ggsave("Results/Fig2/dynamics_IgM.png", p_igm,  width = 12, height = 10)
ggsave("Results/Fig2/dynamics_IgA.png", p_iga,  width = 12, height = 10)
ggsave("Results/Fig2/dynamics_avidity.png", p_avidity, width = 12, height = 10)


# sensitivity analysis

# --- Prepare sensitivity data (day 0 averaged into post-infection)
prepared_igg_data_sens  <- prepare_antibody_data(igg_data,day0_as_post = TRUE)
prepared_igm_data_sens  <- prepare_antibody_data(igm_data,  day0_as_post = TRUE)
prepared_iga_data_sens  <- prepare_antibody_data(iga_data, day0_as_post = TRUE)
prepared_avidity_data_sens <- prepare_antibody_data(avidity_data, day0_as_post = TRUE)


# Plot separately
p_igg_sens <- plot_antibody_dynamics(prepared_igg_data_sens, "IgG") +
  force_panelsizes(rows = unit(3, "cm"), cols = unit(3, "cm"))

p_igm_sens <- plot_antibody_dynamics(prepared_igm_data_sens, "IgM") +
  force_panelsizes(rows = unit(3, "cm"), cols = unit(3, "cm"))

p_iga_sens <- plot_antibody_dynamics(prepared_iga_data_sens, "IgA") +
  force_panelsizes(rows = unit(3, "cm"), cols = unit(3, "cm"))

p_avidity_sens <- plot_antibody_dynamics(prepared_avidity_data_sens, "avidity") +
  force_panelsizes(rows = unit(3, "cm"), cols = unit(3, "cm"))

print(p_igg_sens)
print(p_igm_sens)
print(p_iga_sens)
print(p_avidity_sens)

# --- Save
dir.create("Results/Fig2_sensitivity", recursive = TRUE, showWarnings = FALSE)

ggsave("Results/Fig2_sensitivity/dynamics_IgG_day0_as_post.png",     p_igg_sens, width = 12, height = 10)
ggsave("Results/Fig2_sensitivity/dynamics_IgM_day0_as_post.png",     p_igm_sens, width = 12, height = 10)
ggsave("Results/Fig2_sensitivity/dynamics_IgA_day0_as_post.png",     p_iga_sens, width = 12, height = 10)
ggsave("Results/Fig2_sensitivity/dynamics_avidity_day0_as_post.png", p_avidity_sens, width = 12, height = 10)
