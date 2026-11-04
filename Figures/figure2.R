library(grid)
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(patchwork)
library(ggh4x)

# FLAVI import preprocessed data (log10 RAU)
igg_data_denv <- read.csv("Data/by_isotype_gm/denv/logged_preprocessed_cebu_IgG.csv")
iga_data_denv <- read.csv("Data/by_isotype_gm/denv/logged_preprocessed_cebu_IgA.csv")
igm_data_denv <- read.csv("Data/by_isotype_gm/denv/logged_preprocessed_cebu_IgM.csv")



# ALPHA import preprocessed data (log10 RAU)
igg_data_alpha <- read.csv("Data/by_isotype_gm/chik/logged_preprocessed_cebu_IgG.csv")
iga_data_alpha <- read.csv("Data/by_isotype_gm/chik/logged_preprocessed_cebu_IgA.csv")
igm_data_alpha <- read.csv("Data/by_isotype_gm/chik/logged_preprocessed_cebu_IgM.csv")


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
    confirmed = c("CHIKV"),
    row_order = c("E2", "NSP123", "VLP"),
    col_order = c("CHIKV")
  )
)

# --- Prepare data for a given isotype 
prepare_antibody_data <- function(data, 
                                  antigens = c(dengue_zika_antigens, chik_onnv_mayv_antigens)) {
  data %>%
    dplyr::select(id_patient, id_sample, days_since_infection, isotype, PCR = target, all_of(antigens)) %>%
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
      str_detect(antigen, "^SHERPADES_") ~ paste0(str_extract(antigen, "[^_]+$"), "-SHERPADES"),
      TRUE ~ str_extract(antigen, "[^_]+$")
      ),
      virus_family = case_when(
        pathogen %in% c("DENV1", "DENV2", "DENV3", "DENV4", "ZIKV") ~ "Flavivirus",
        pathogen %in% c("CHIKV", "ONNV", "MAYV", "RRV") ~ "Alphavirus",
        TRUE ~ NA_character_
      ),
      days_since_infection = as.numeric(days_since_infection),
    ) %>%
    dplyr::mutate(
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

  fam      <- defs[[family]]
  isotypes <- isotypes %||% unique(data$isotype)
  multi    <- length(isotypes) > 1

  d <- data %>%
    dplyr::filter(isotype %in% isotypes,
                  PCR %in% fam$confirmed,
                  virus_family == family) %>%
    dplyr::mutate(
      facet_pathogen = factor(PCR, levels = fam$col_order),
      isotype        = factor(isotype, levels = isotypes),
      antigen_type   = factor(antigen_type,
                              levels = intersect(c(fam$row_order, sort(unique(antigen_type))),
                                                 unique(antigen_type)))
    )

  y_label <- if (identical(isotypes, "avidity")) "Avidity"
             else if (multi) "Antibody Titre (log10)"
             else paste0(isotypes, " Antibody Titre (log10)")

  f <- if (multi) antigen_type ~ isotype + facet_pathogen
       else       antigen_type ~ facet_pathogen

  fg <- d %>% dplyr::filter(pathogen == facet_pathogen)
  bg <- d %>% dplyr::filter(pathogen != facet_pathogen)

  ggplot(d, aes(x = x_position, y = value,
                group = interaction(id_patient, pathogen, antigen, isotype))) +
    geom_line(data = bg,  color = "grey80", alpha = 0.5, linewidth = 0.5) +
    geom_point(data = bg, color = "grey80", alpha = 0.5, size = 1.5) +
    geom_line(data = fg,  aes(color = PCR), alpha = 0.8, linewidth = 0.7) +
    geom_point(data = fg, aes(color = PCR), alpha = 0.8, size = 1.8) +
    ggh4x::facet_nested(
      f, scales = "free_y", drop = TRUE,
            labeller = labeller(antigen_type = function(x) stringr::str_replace(x, "-", "\n")),
      strip = ggh4x::strip_nested(bleed = FALSE),
      nest_line = element_line(colour = "grey50")
    ) + 
    scale_color_manual(values = pcr_colours, name = "PCR Confirmed") +
    scale_x_continuous(breaks = c(-30, 0, 30), labels = c("Pre", "0", ">30"),
                       expand = expansion(mult = 0.15)) +
    labs(x = "Days since PCR+ve infection", y = y_label) +
    theme_bw() +
    theme(
      strip.text.x       = element_text(size = 16),
      strip.text.y       = element_text(size = 16, angle = 0, hjust = 0),
      strip.background.x = element_rect(fill = "#ffffff", colour = NA),
      strip.background.y = element_blank(),
      axis.text          = element_text(size = 14),
      axis.title         = element_text(size = 16),
      panel.spacing.x    = unit(1.2, "lines"),
      panel.spacing.y    = unit(0.8, "lines"),
      legend.position    = "none",
      panel.grid.minor   = element_blank(),
      plot.subtitle      = element_text(size = 14, face = "italic"),
      plot.margin        = margin(t = 10, r = 25, b = 10, l = 10),
      plot.background    = element_rect(fill = "white", color = "black", linewidth = 0.5)
    )
}

all_antigens <- c(dengue_zika_antigens, chik_onnv_mayv_antigens)

prepared_data <- prepare_antibody_data(
  dplyr::bind_rows(igg_data_denv, igm_data_denv, iga_data_denv,
                   igg_data_alpha, igm_data_alpha, iga_data_alpha),
  antigens = all_antigens
)

isos <- c("IgG", "IgA", "IgM")

flavi_dynamics <- plot_antibody_dynamics(prepared_data, "Flavivirus", isotypes = isos)
alpha_dynamics <- plot_antibody_dynamics(prepared_data, "Alphavirus", isotypes = isos)

dir.create("Results/Fig2", recursive = TRUE, showWarnings = FALSE)

ggsave("Results/Fig2/dynamics_flavi.png", flavi_dynamics,
       width = 25, height = 8, dpi = 300, bg = "white", limitsize = FALSE)

ggsave("Results/Fig2/dynamics_alpha.png", alpha_dynamics,
       width = 10, height = 7, dpi = 300, bg = "white", limitsize = FALSE)








prepared_data %>%
  dplyr::filter(PCR == "ZIKV", pathogen == "ZIKV") %>%
  dplyr::count(isotype, time_bin, antigen_type) %>%
  tidyr::pivot_wider(names_from = isotype, values_from = n)


prepared_data %>%
  dplyr::filter(PCR == "ZIKV", is.nan(value)) %>%
  dplyr::count(isotype, time_bin)


dplyr::bind_rows(
  igg_data_denv %>% dplyr::mutate(.src = "IgG"),
  iga_data_denv %>% dplyr::mutate(.src = "IgA"),
  igm_data_denv %>% dplyr::mutate(.src = "IgM")
) %>%
  dplyr::filter(PCR == "ZIKV") %>%
  dplyr::count(.src, id_patient, days_since_infection) %>%
  tidyr::pivot_wider(names_from = .src, values_from = n)


zikv_id <- prepared_data %>%
  dplyr::filter(PCR == "ZIKV") %>%
  dplyr::pull(id_patient) %>% unique()
zikv_id

dplyr::bind_rows(
  igg_data_denv %>% dplyr::mutate(.src = "IgG"),
  iga_data_denv %>% dplyr::mutate(.src = "IgA"),
  igm_data_denv %>% dplyr::mutate(.src = "IgM")
) %>%
  dplyr::filter(id_patient %in% zikv_id) %>%
  dplyr::count(.src, days_since_infection, PCR)


dplyr::bind_rows(
  igg_data_denv %>% dplyr::mutate(.src = "IgG"),
  iga_data_denv %>% dplyr::mutate(.src = "IgA"),
  igm_data_denv %>% dplyr::mutate(.src = "IgM")
) %>%
  dplyr::filter(stringr::str_detect(id_patient, stringr::fixed(zikv_id[1]))) %>%
  dplyr::distinct(.src, id_patient) %>%
  dplyr::mutate(raw = sprintf("[%s]", id_patient))   # exposes whitespace
