library(grid)
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(patchwork)
library(ggh4x)

# import preprocessed data (log10 RAU)
igg_data <- read.csv("Data/by_isotype/logged_preprocessed_cebu_IgG.csv")
iga_data <- read.csv("Data/by_isotype/logged_preprocessed_cebu_IgA.csv")
igm_data <- read.csv("Data/by_isotype/logged_preprocessed_cebu_IgM.csv")
avidity_data <- read.csv("Data/by_isotype/logged_preprocessed_cebu_avidity.csv")



# Define antigen groups
dengue_zika_antigens <- c(
  "DENV1_DIII", "DENV1_VLP", "DENV1_NS1","SHERPADES_DENV1_DIII",
  "DENV2_DIII", "DENV2_VLP", "DENV2_NS1", "SHERPADES_DENV2_DIII", 
  "DENV3_DIII", "DENV3_VLP", "DENV3_NS1", "SHERPADES_DENV3_DIII",
  "DENV4_DIII", "DENV4_VLP", "DENV4_NS1", "SHERPADES_DENV4_DIII",
  "ZIKV_DIII", "ZIKV_VLP", "ZIKV_NS1", "SHERPADES_ZIKV_DIII"
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


family_defs <- list(
  Flavivirus = list(
    confirmed = c("DENV1", "DENV2", "DENV3", "DENV4", "ZIKV"),
    row_order = c("DIII", "SHERPADES-DIII", "NS1", "VLP"),
    col_order = c("DENV1", "DENV2", "DENV3", "DENV4", "ZIKV")
  ),
  Alphavirus = list(
    confirmed = c("CHIKV"),                    # PCR-confirmable alphaviruses
    row_order = c("E2", "NSP123", "VLP"),
    col_order = c("CHIKV")
  )
)

# --- Prepare data for a given isotype 
prepare_antibody_data <- function(data, 
                                  antigens = c(dengue_zika_antigens, chik_onnv_mayv_antigens), 
                                  isotype = NULL,
                                  day0_as_post = FALSE) {
  data %>%
    dplyr::select(id_patient, id_sample, isotype, days_since_infection, PCR, all_of(antigens)) %>%
    tidyr::pivot_longer(
      cols = all_of(antigens),
      names_to = "antigen",
      values_to = "value"
    ) %>%
    mutate(
      pathogen = case_when(
        str_detect(antigen, "^SHERPADES_") ~ str_extract(antigen, "(?<=SHERPADES_)[^_]+"),
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
    dplyr::group_by(id_patient) %>%
    tidyr::fill(PCR, .direction = "downup") %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
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
      color_group = dplyr::if_else(pathogen == PCR, PCR, "Other")
    ) %>%
    dplyr::group_by(id_patient, isotype, antigen, time_bin, x_position) %>%
    dplyr::summarise(
      value = mean(value, na.rm = TRUE),
      PCR =  dplyr::first(PCR),
      color_group = dplyr::first(color_group),
      antigen_type = dplyr::first(antigen_type),
      pathogen = dplyr::first(pathogen),
      virus_family = dplyr::first(virus_family),
      .groups = "drop"
    ) 
}

# --- single plot for one isotype + virus family 
# isotypes = NULL  -> every isotype in the data, as nested facet rows
# isotypes = "IgG" -> single isotype, plain antigen rows
plot_antibody_dynamics <- function(data,
                                   family,
                                   isotypes = NULL, 
                                   subtitle = NULL,
                                   defs = family_defs) {

  fam <- defs[[family]]
  isotypes <- isotypes %||% unique(data$isotype)
  multi    <- length(isotypes) > 1


   d <- data %>%
    dplyr::filter(isotype %in% isotypes,
                  PCR %in% fam$confirmed,     # facet columns: this family's infections
                  virus_family == family) %>% # antigens from the same family only
    dplyr::mutate(
      facet_pathogen = factor(PCR, levels = fam$col_order),
      isotype        = factor(isotype, levels = isotypes),
      # row_order first, then any unlisted types alphabetically after
      antigen_type   = factor(antigen_type,
                              levels = intersect(c(fam$row_order, sort(unique(antigen_type))),
                                                 unique(antigen_type)))
    )
 
  y_label <- if (identical(isotypes, "avidity")) "Avidity"
             else if (multi) "Antibody Titre (log2)"
             else paste0(isotypes, " Antibody Titre (log2)")

    f <- if (multi) isotype + antigen_type ~ facet_pathogen
       else    antigen_type ~ facet_pathogen
 
 
  fg <- d %>% dplyr::filter(pathogen == facet_pathogen)  # coloured: matches the facet
  bg <- d %>% dplyr::filter(pathogen != facet_pathogen)  # grey: other same-family antigens
 
  ggplot(d, aes(x = x_position, y = value,
                group = interaction(id_patient, pathogen, antigen, isotype))) +
    geom_line(data = bg,  color = "grey80", alpha = 0.5, linewidth = 0.5) +
    geom_point(data = bg, color = "grey80", alpha = 0.5, size = 1.5) +
    geom_line(data = fg,  aes(color = PCR), alpha = 0.8, linewidth = 0.7) +
    geom_point(data = fg, aes(color = PCR), alpha = 0.8, size = 1.8) +
    ggh4x::facet_nested(
      f, scales = "free_y", drop = TRUE,
      labeller = labeller(antigen_type = label_wrap_gen(12)),
      strip = ggh4x::strip_nested(bleed = FALSE),
      nest_line = element_line(colour = "grey50")
    ) +
    scale_color_manual(values = pcr_colours, name = "PCR Confirmed") +
    scale_x_continuous(breaks = c(-30, 0, 30), labels = c("Pre", "0", ">30")) +
    labs(x = "Days since PCR+ve infection", y = y_label,
         subtitle = subtitle %||% family) +
    theme_bw() +
    theme(
      strip.text.x     = element_text(size = 16),
      strip.text.y     = element_text(size = 12, angle = 0, hjust = 0),
      strip.background = element_rect(fill = "#ffffff"),
      axis.text        = element_text(size = 12),
      axis.title       = element_text(size = 16),
      panel.spacing.x  = unit(1.5, "lines"),
      legend.position = "none",
      panel.spacing.y  = unit(0.8, "lines"),
      panel.grid.minor = element_blank(),
      plot.subtitle    = element_text(size = 14, face = "italic"),
      plot.margin      = margin(t = 10, r = 25, b = 10, l = 10),
      plot.background  = element_rect(fill = "white", color = "black", linewidth = 0.5)
    )
}



all_antigens <- c(dengue_zika_antigens, chik_onnv_mayv_antigens)

prepared_data <- prepare_antibody_data(
  dplyr::bind_rows(igg_data, igm_data, iga_data),
  antigens = all_antigens
)


isos <- c("IgG", "IgA", "IgM")

fig_flavi <- purrr::imap(purrr::set_names(isos), function(iso, i) {
  p <- plot_antibody_dynamics(prepared_data, "Flavivirus", isotypes = iso) +
    labs(subtitle = NULL, title = iso) +
    theme(
      plot.title = ggtext::element_textbox(
        fill = "grey85", colour = "black", width = unit(1, "npc"),
        padding = margin(4, 8, 4, 8), margin = margin(b = 4),
        halign = 0.5, size = 15, face = "bold"
      ),
      plot.background = element_blank()
    )
  if (iso != tail(isos, 1)) p <- p + theme(axis.title.x = element_blank(),
                                           axis.text.x  = element_blank())
  if (iso != isos[1])       p <- p + theme(strip.text.x = element_blank())
  p
})


flavi_dynamics <- patchwork::wrap_plots(fig_flavi, ncol = 1)  



fig_alpha <- purrr::imap(purrr::set_names(isos), function(iso, i) {
  p <- plot_antibody_dynamics(prepared_data, "Alphavirus", isotypes = iso) +
    labs(subtitle = NULL, title = iso) +
    theme(
      plot.title = ggtext::element_textbox(
        fill = "grey85", colour = "black", width = unit(1, "npc"),
        padding = margin(4, 8, 4, 8), margin = margin(b = 4),
        halign = 0.5, size = 15, face = "bold"
      ),
      plot.background = element_blank()
    )
  if (iso != tail(isos, 1)) p <- p + theme(axis.title.x = element_blank(),
                                           axis.text.x  = element_blank())
  if (iso != isos[1])       p <- p + theme(strip.text.x = element_blank())
  p
})

alpha_dynamics <- patchwork::wrap_plots(fig_alpha, ncol = 1)  


# --- create output folder
dir.create("Results/Fig2", recursive = TRUE, showWarnings = FALSE)
ggsave("Results/Fig2/dynamics_flavi.png", flavi_dynamics,
       width = 12, height = 14, dpi = 300, limitsize = FALSE)

ggsave("Results/Fig2/dynamics_alpha.png", alpha_dynamics,
       width = 5, height = 14, dpi = 300, limitsize = FALSE)















# ---- OLD CODE ------

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
quartz()
print(p_igg)
quartz()
print(p_igm)
print(p_iga)
print(p_avidity)



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
