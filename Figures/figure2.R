
library(grid)

# import ratio df 
logged_preprocessed_cebu_data <- read.csv('Results/logged_preprocessed_cebu_data.csv')

colnames(logged_preprocessed_cebu_data)



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
  "MAYV_E2" 
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
prepare_antibody_data <- function(data, isotype, antigens = c(flavivirus_antigens, alphavirus_antigens)) {
  
  data %>%
    filter(isotype == !!isotype) %>%
    dplyr::select(id_patient, id_sample, days_since_infection, PCR,
                  all_of(antigens)) %>%
    pivot_longer(
      cols     = all_of(antigens),
      names_to = "antigen"
    ) %>%
    mutate( # extract pathogen
      pathogen = case_when(
        str_detect(antigen, "^SHERPADES_") ~ str_extract(antigen, "(?<=SHERPADES_)[^_]+"),
        str_detect(antigen, "^ZIKVSU")     ~ "ZIKV",
        str_detect(antigen, "^ZIKVAS")     ~ "ZIKV",
        TRUE                               ~ str_extract(antigen, "^[^_]+")
      ), # extract antigen type
      antigen_type = case_when(
        str_detect(antigen, "^SHERPADES_") ~ "DIII-SHERPADES",
        TRUE                               ~ str_extract(antigen, "[^_]+$")
      ),
      virus_family = case_when(
        pathogen %in% c("DENV1", "DENV2", "DENV3", "DENV4", "ZIKV") ~ "Flavivirus",
        pathogen %in% c("CHIKV", "ONNV", "MAYV", "RRV")  ~ "Alphavirus"
      ),
      days_since_infection = as.numeric(days_since_infection),
      x_position = case_when(
        days_since_infection < 0 ~ -30, # all pre-infection samples will be plotted at -30
        days_since_infection > 30  ~ 30, # ALL post >30 days → 30
        TRUE  ~ days_since_infection) # 0 to 30 kept as-is
    ) %>%                                        
    group_by(id_patient) %>%
    fill(PCR, .direction = "downup") %>%
    ungroup() %>%
    mutate(                                     
      color_group = if_else(pathogen == PCR, PCR, "Other")
    ) %>%

    # If >1 sample maps to the same x_position (e.g. >30 days, or <-30 days),
    # take the average value across those samples per patient/antigen
    group_by(id_patient, antigen, x_position) %>%
    summarise(value = mean(value, na.rm = TRUE),
           PCR = first(PCR),
           color_group = first(color_group),
           antigen_type = first(antigen_type),
           pathogen = first(pathogen),
           virus_family = first(virus_family),
          .groups = "drop") %>%
    arrange(color_group == "Other", color_group) # ← grey lines drawn first
}

# --- single plot for one isotype + virus family 
plot_antibody_dynamics <- function(data,
                                   isotype,
                                   flavi_or_alpha = c("Flavivirus", "Alphavirus")) {
  
  flavi_or_alpha <- match.arg(flavi_or_alpha)
  
  # All pathogens per family (used for background grey lines)
  all_pathogens <- list(
    Flavivirus = c("DENV1", "DENV2", "DENV3", "DENV4", "ZIKV"),
    Alphavirus = c("CHIKV", "ONNV", "MAYV")
  )

  # Pathogens with confirmed infections (used for facet columns)
  confirmed_pathogens <- list(
  Flavivirus = c("DENV1", "DENV2", "DENV3", "DENV4", "ZIKV"),
  Alphavirus = c("CHIKV")
  )

  y_label <- if (isotype == "avidity") "Avidity" else "Antibody Titre (log2)"

  confirmed <- confirmed_pathogens[[flavi_or_alpha]]
  all_family <- all_pathogens[[flavi_or_alpha]]
  

  # Background: all pathogens in family, in grey
  background_data <- data %>%
    filter(pathogen %in% all_family) %>%
    mutate(color_group = "Other") %>%
    tidyr::crossing(facet_pathogen = confirmed)

  # Foreground data: only confirmed pathogens
  facet_data <- data %>%
    filter(pathogen %in% confirmed) %>%
    mutate(
      facet_pathogen = pathogen
    )

  
  ggplot(data = facet_data,
         aes(x = x_position, y = value, group = id_patient)) +
    
    # Layer 1: all family pathogens in grey (background)
    geom_line(data  = background_data,
              color = "grey90", alpha = 0.7, linewidth = 0.5) +
    geom_point(data = background_data,
               color = "grey90", alpha = 0.7, size = 2) +
    
    # Layer 2: PCR-confirmed coloured lines on top
    geom_line(data = ~filter(.x, color_group != "Other"),
              aes(color = color_group), alpha = 0.8, linewidth = 0.7) +
    geom_point(data = ~filter(.x, color_group != "Other"),
               aes(color = color_group), alpha = 0.8, size = 2) +
    
    facet_grid(antigen_type ~ facet_pathogen, scales = "free_y", drop = TRUE) +

    scale_color_manual(values = pcr_colours, name = "PCR Confirmed") +
    scale_x_continuous(
      breaks = c(-30, -15, 0, 15,  30),
      labels = c("<-30", -15, 0, 15, ">30")
    ) +
    
    labs(
      x     = "Days since PCR+ve infection",
      y     = y_label) +
    
    theme_bw() +
    theme(
      strip.text       = element_text(size = 20),
      strip.background = element_rect(fill = "#ffffff"),
      axis.text        = element_text(size = 20),
      axis.text.x      = element_text(size = 20),
      axis.text.y      = element_text(size = 20),
      axis.title       = element_text(size = 20),
      legend.position  = "bottom",
      legend.text      = element_text(size = 20),
      legend.title     = element_text(size = 20),
      legend.direction = "horizontal",
      legend.margin    = margin(0, 0, 0, 0),
      panel.spacing.x  = unit(1, "lines"), 
      panel.spacing.y  = unit(1, "lines"),
      panel.grid.minor = element_blank(),
      plot.background  = element_rect(fill = "white", color = "black", linewidth = 0.5)
    )
    
}



igg_data <- prepare_antibody_data(logged_preprocessed_cebu_data, "IgG", antigens = c(flavivirus_antigens, alphavirus_antigens))
igm_data <- prepare_antibody_data(logged_preprocessed_cebu_data, "IgM", antigens = c(flavivirus_antigens, alphavirus_antigens))
iga_data <- prepare_antibody_data(logged_preprocessed_cebu_data, "IgA", antigens = c(flavivirus_antigens, alphavirus_antigens))
avidity_data <- prepare_antibody_data(logged_preprocessed_cebu_data, "avidity", antigens = c(flavivirus_antigens, alphavirus_antigens))

# Plot separately
igg_flavi <- plot_antibody_dynamics(igg_data, "IgG", "Flavivirus")
igg_alpha <- plot_antibody_dynamics(igg_data, "IgG", "Alphavirus")

print(igg_flavi)
print(igg_alpha)

# save IgG 
ggsave("Results/flavivirus_IgG_dynamics.png", 
igg_flavi, width = 25, height = 12)

# save IgG 
ggsave("Results/alphavirus_IgG_dynamics.png", 
igg_alpha, width = 8, height = 10)


# other isotype plots 
igm_flavi <- plot_antibody_dynamics(igm_data, "IgM", "Flavivirus")
iga_flavi <- plot_antibody_dynamics(iga_data, "IgA", "Flavivirus")
avidity_flavi <- plot_antibody_dynamics(avidity_data, "Avidity", "Flavivirus")

print(iga_flavi)
print(igm_flavi)
print(avidity_flavi)


logged_preprocessed_cebu_data$isotype %>% unique()
