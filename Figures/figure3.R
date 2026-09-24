library(ggtext)
library(DescTools)
library(ggforce)
library(stringr)
library(dplyr)
library(patchwork)
library(corrplot)


ratio_df_IgG <- readRDS('Data/by_isotype/ratio_df_IgG.rds')

ratio_df_IgG_gmean <- readRDS('Data/by_isotype_gmean/ratio_df_IgG.rds')

pcr_target <- c("DENV1", "DENV2", "DENV3", "DENV4", "ZIKV")

flavi_antigens <- c("DENV1_DIII","DENV1_NS1","DENV1_VLP","SHERPADES_DENV1_DIII",
"DENV2_DIII","DENV2_NS1","DENV2_VLP","SHERPADES_DENV2_DIII",
"DENV3_DIII","DENV3_NS1","DENV3_VLP", "SHERPADES_DENV3_DIII",
"DENV4_DIII","DENV4_NS1","DENV4_VLP", "SHERPADES_DENV4_DIII",
"JEV_E", "JEV_NS1", "SHERPADES_JEV_DIII",
"YFV_E", "YFV_NS1", "SHERPADES_YFV_DIII",
"WNV_DIII","WNV_NS1","SHERPADES_WNV_DIII",
"ZIKV_DIII", "ZIKV_NS1","ZIKV_VLP","SHERPADES_ZIKV_DIII")

antigen_blocks <- list(
  "VLP / E" = c("VLP", "E"),
  "NS1"  = "NS1",
  "DIII (SHERPADES)"  = "SHERPADES"
)



# --- Functions ---
# antigen column name for a given pathogen + panel suffix
get_antigen_col <- function(pathogen, suffix) {
    if (suffix == "DIII" && pathogen == "ZIKV") return("ZIKVAS_DIII")
    if (suffix == "SHERPADES")   return(paste0("SHERPADES_", pathogen, "_DIII"))
    paste0(pathogen, "_", suffix)
  }

  # reverse map: which pathogen does this antigen column belong to?
get_antigen_pathogen <- function(col_name) {
    if (str_detect(col_name, "^SHERPADES_")) return(str_extract(col_name, "(?<=^SHERPADES_)[^_]+"))
    if (str_detect(col_name, "^ZIKVAS_"))    return("ZIKV")
    str_extract(col_name, "^[^_]+")
  }
 
  # pretty antigen label for axes / facets
fmt_antigen <- function(x) dplyr::case_when(
    str_detect(x, "^SHERPADES_") ~ str_replace(x, "^SHERPADES_([^_]+)_DIII$", "SHERPADES \\1 DIII"),
    str_detect(x, "^ZIKVAS_")    ~ "ZIKV DIII",
    TRUE                         ~ str_replace(x, "_", " ")
  )


# -- prepare data ---
prep_data <- function(data, 
                      antigen, 
                      flavi_targets = pcr_target,
                      antigens = flavi_antigens) {

if (identical(antigen, "SHERPADES")) {
    filtered_antigens <- antigens[str_detect(antigens, "^SHERPADES_.*_DIII$")]
  } else {
    pattern <- paste0("_(", paste(antigen, collapse = "|"), ")$")
    filtered_antigens <- antigens[
      str_detect(antigens, pattern) & !str_detect(antigens, "^SHERPADES_")
    ]
  }


  if (length(filtered_antigens) == 0) {
    stop("No antigens matched antigen = '", antigen, "'.\n  Available: ",
         paste(antigens, collapse = ", "))
  }
  
  missing_cols <- setdiff(filtered_antigens, names(data))
  
  if (length(missing_cols) > 0) {
    stop("These expected columns are not present in `data`: ",
         paste(missing_cols, collapse = ", "))
  }
  
  # antigen -> pathogen lookup
  antigen_pathogen <- vapply(filtered_antigens, get_antigen_pathogen, character(1))

  plot_data <- data %>%
    filter(target %in% flavi_targets) %>%
    dplyr::select(id_patient, target, all_of(filtered_antigens)) %>%
    pivot_longer(
      cols      = all_of(filtered_antigens),
      names_to  = "antigen",
      values_to = "ratio"
    ) %>%
    mutate(
      infecting_target = target,
      antigen_pathogen = unname(antigen_pathogen[antigen]),
      reactivity = if_else(infecting_target == antigen_pathogen,
                           "Homologous", "Heterologous"),
      antigen_label = factor(fmt_antigen(antigen),
                             levels = fmt_antigen(filtered_antigens)),
      target_label  = factor(infecting_target, levels = flavi_targets)
    )

  list(
    plot_data    = plot_data,
    antigen_cols = filtered_antigens
  )}



denv_serotypes <- c("DENV1", "DENV2", "DENV3", "DENV4")

antigen_blocks <- list(
  "VLP / E" = c("VLP", "E"),
  "NS1"                = "NS1",
  "DIII (SHERPADES)"   = "SHERPADES"
)

res_list <- lapply(antigen_blocks, function(a) prep_data(ratio_df_IgG, a))

all_ratios <- unlist(lapply(res_list, function(r) r$plot_data$ratio))
x_lim <- c(min(all_ratios, na.rm = TRUE) * 0.8,
           max(all_ratios, na.rm = TRUE) * 1.2)

make_panel <- function(res, title, suffixes, show_x = FALSE) {

  homol_cells <- res$plot_data %>%
    filter(reactivity == "Homologous") %>%
    distinct(target_label, antigen_label)

  p <- ggplot(res$plot_data, aes(x = ratio, fill = reactivity)) +
        geom_rect(data = homol_cells,
              xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf,
              fill = "#fffcfc", inherit.aes = FALSE) +
    geom_vline(xintercept = 1, linetype = "dashed", colour = "red") +
    geom_histogram(binwidth = 0.25, colour = NA) +
    facet_grid(target_label ~ antigen_label, switch = "y") +
    scale_x_log10(breaks = c(0.1, 1, 10, 100),
                  labels = c("0.1", "1", "10", "100")) +
    scale_y_continuous(breaks = scales::breaks_pretty(n = 3),
                       expand = expansion(mult = c(0, 0.05))) +
    scale_fill_manual(
      values = c(Homologous = "palevioletred", Heterologous = "steelblue"),
      breaks = c("Homologous", "Heterologous")
    ) +
    coord_cartesian(xlim = x_lim) +
    labs(title = title, x = NULL, y = "Count") +
    theme_bw(base_size = 8) +
    theme(panel.grid.minor   = element_blank(),
          strip.background   = element_rect(fill = "white"),
          panel.spacing.y    = unit(0.8, "lines"),
          strip.text.x     = element_text(size = 16),
          strip.text.y     = element_text(size = 16, angle = 0, hjust = 0),
          strip.text.y.left  = element_text(angle = 0),
          axis.text.x        = element_text(size = 16),
          axis.text.y        = element_text(size = 16),
           axis.title         = element_text(size = 16),
          plot.title         = element_text(face = "bold", size = 9),
          legend.position    = "none")

  if (!show_x) {
    p <- p + theme(axis.text.x = element_blank(),
                   axis.ticks.x = element_blank())
  }
  p
}

panels <- lapply(seq_along(antigen_blocks), function(i) {
  make_panel(res_list[[i]],
             title    = "",
             suffixes = antigen_blocks[[i]],
             show_x   = i == length(antigen_blocks))
})


combined_igg_ratio_panel <- wrap_plots(panels, ncol = 1) 


dir.create("Results/Fig3", recursive = TRUE, showWarnings = FALSE)
ggsave("Results/Fig3/ratio_histograms_IgG.png", combined_igg_ratio_panel,
       width = 23, height = 15, dpi = 300)



# ---- Summary Nos : how many id_patient have ratio < 1 - mean thath post titre is greater than pre titre
flavi_gt1_summary <- ratio_long %>%
  filter(
    target  %in% pcr_target,          # same target filter as prep_data()
    antigen %in% flavi_antigens,
    component %in% c("VLP", "E", "NS1", "DIII (SHERPADES)")
  ) %>%
  mutate(
    antigen = factor(
      if_else(component %in% c("VLP", "E"), "VLP / E", component),
      levels = names(antigen_blocks)
    ),
    reactivity = factor(
      if_else(target == virus, "Homologous", "Heterologous"),
      levels = c("Homologous", "Heterologous")
    )
  ) %>%
  group_by(isotype, antigen, reactivity) %>%
  summarise(
    n               = n(),
    n_gt_1          = sum(ratio > 1),
    prop_gt_1       = mean(ratio > 1),
    median_ratio    = median(ratio),
    n_samples      = n_distinct(id_patient),
    .groups = "drop"
  ) %>%
  arrange(isotype, antigen, reactivity)

flavi_gt1_summary %>% print(n = Inf)


flavi_gt1_wide <- flavi_gt1_summary %>%
  dplyr::mutate(summary = sprintf("%d/%d (%.0f%%)", n_gt_1, n, 100 * prop_gt_1)) %>%
  dplyr::select(isotype, antigen, reactivity, summary) %>%
  tidyr::pivot_wider(names_from = reactivity, values_from = summary)

flavi_gt1_wide





# --- correlation plot 
colnames(ratio_df_IgG_gmean)
names(ratio_df_IgG_gmean)[names(ratio_df_IgG_gmean) == "ZIKVAS_DIII"] <- "ZIKV_DIII"


ratio_df_IgG_flavi <- ratio_df_IgG %>%
    dplyr::select(all_of(flavi_antigens)) %>%
    as.data.frame()

cor_mat <- cor(ratio_df_IgG_flavi, use = "pairwise.complete.obs", method = "pearson")

cor_df <- reshape2::melt(cor_mat, varnames = c("Antigen1", "Antigen2"),
               value.name = "correlation", na.rm = TRUE)


grp_order   <- c("DIII", "SHERPADES", "NS1", "VLP", "E")
virus_order <- c("DENV1", "DENV2", "DENV3", "DENV4", "ZIKV", "JEV", "WNV", "YFV")

grp <- ifelse(grepl("^SHERPADES", flavi_antigens),
              "SHERPADES",
              sub(".*_", "", flavi_antigens))

vir <- sub("^SHERPADES_", "", flavi_antigens)
vir <- sub("_(DIII|NS1|VLP|E)$", "", vir)

ord <- order(match(grp, grp_order), match(vir, virus_order))
antigen_levels <- flavi_antigens[ord]

cor_df$Antigen1 <- factor(cor_df$Antigen1, levels = antigen_levels)
cor_df$Antigen2 <- factor(cor_df$Antigen2, levels = antigen_levels)

i <- as.integer(cor_df$Antigen1)
j <- as.integer(cor_df$Antigen2)

cor_half <- cor_df[j >= i, ]   # flip to j <= i for the other triangle
                               # use > / < to also drop the diagonal


quartz()
ggplot(cor_half, aes(x = Antigen1, y = Antigen2, fill = correlation)) +
  geom_tile() +
  scale_x_discrete(drop = FALSE) +
  scale_y_discrete(drop = FALSE) +
  scale_fill_gradient(low = "steelblue", high = "firebrick",
                      limits = c(0, 1), oob = scales::squish,
                      name = "Spearman\ncorrelation") +
  coord_fixed() + 
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 10),
    axis.text.y = element_text(size = 10),
    panel.grid = element_blank(),
    axis.title = element_blank()
  )















# ---- OLD Plots ------ 
# Including all ratio plots (diagonal + off-diagonal) for each isotype
# Forest plot 

plot_ratio_forest <- function(res,
                              pad = 0.5,
                              max_breaks = 8) {

  stopifnot(is.list(res), !is.null(res$forest_data))

  forest_data <- res$forest_data
  conf_level  <- if (is.null(res$conf_level)) 0.95 else res$conf_level

  facet_rows <- vars(target_label)     # <-- changed: facet by infecting pathogen
  y_var      <- "antigen_label"        # <-- changed: y-axis = antigen tested

  off_diag <- forest_data %>% dplyr::filter(panel_type == "off_diagonal")

  x_min  <- min(c(off_diag$gm_lower, off_diag$log_gm_rel, 0), na.rm = TRUE) - pad
  x_max  <- max(c(off_diag$gm_upper, off_diag$log_gm_rel, 0), na.rm = TRUE) + pad
  brks   <- scales::breaks_pretty(n = max_breaks)(c(x_min, x_max))

  ggplot(forest_data, aes(x = log_gm_rel, y = .data[[y_var]])) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
    geom_errorbar(
      data        = off_diag,
      aes(xmin    = gm_lower, xmax = gm_upper),
      orientation = "y",
      width       = 0.25,
      linewidth   = 0.5
    ) +
    geom_point(aes(shape = panel_type, colour = panel_type, fill = panel_type), size = 3) +
    scale_shape_manual(values = c(diagonal = 21, off_diagonal = 19), guide = "none") +
    scale_colour_manual(values = c(diagonal = "#e972a7", off_diagonal = "#024a9c"), guide = "none") +
    scale_fill_manual(values = c(diagonal = "#e972a7", off_diagonal = "#024a9c"), guide = "none") +
    scale_y_discrete(limits = rev) +
    facet_col(facet_rows, scales = "free_y", space = "free", strip.position = "top") +
    scale_x_continuous(breaks = brks, labels = brks) +
    coord_cartesian(xlim = c(x_min, x_max), clip = "off") +
    labs(
      x  = expression(log[2]~"(Geometric mean relative ratio)"),
      y  = ""
    ) +
    theme_bw(base_size = 12) +
    theme(
      strip.placement  = "outside",
      strip.text.x  = element_text(size = 20, angle = 0),
      panel.border   = element_rect(colour = "black", fill = NA, linewidth = 0.5),
      panel.spacing.y    = unit(0.4, "lines"),
      strip.background   = element_rect(fill = "#ffffff", colour = "black",
                                        linewidth = 0.5),
      axis.text   = element_text(size = 20),
      axis.title.x  = element_text(size = 20),
      axis.title.y   = element_text(size = 20),
      plot.title.position = "plot",
      plot.subtitle  = element_text(size = 20, hjust = 0),
      axis.line.x = element_line(colour = "black"),
      plot.margin  = margin(t = 10, r = 20, b = 10, l = 10),
      plot.background    = element_rect(fill = "white", colour = NA)
    )
}


plot_ratio_histogram <- function(res,
                                 label_fmt = "GM = %.2f",
                                 ratio_label_fmt = "Ratio = %.2f",
                                 bins = 10) {
 
  stopifnot(is.list(res), !is.null(res$plot_data), !is.null(res$mean_data))
 
  #diagonal -> geometric mean (GM)
  #off-diagonal -> relative ratio, GM[current antigen] / GM[infecting pathogen]
  gm_labels <- res$mean_data %>%
    dplyr::filter(!is.na(gm)) %>%
    mutate(
      label = dplyr::case_when(
        panel_type == "diagonal"     ~ sprintf(label_fmt,       log2(gm)),
        panel_type == "off_diagonal" ~ sprintf(ratio_label_fmt, log2(gm)),
        TRUE                         ~ NA_character_
      )
    )

  antigen_labeller <- as_labeller(fmt_antigen_multiline)
  # diagonal data == GM 
  # off diagonal data == (GM[current antigen] / GMp[infecting pathogen]) == relative ratio

  # values are plotted on a log2 scale - 
  # diagonal == log2(GM)
  # off-diagonal == log2(relative ratio)
  # if relative ratio > 1 --> cross reactive and log2(relative ratio) is pos else, log2(relative ratio) is neg
  ggplot(res$plot_data, aes(x = log2(value), fill = panel_type)) +
    geom_histogram(bins = bins, alpha = 0.8) +
    geom_vline(xintercept = 0, colour = "red", linetype = "dashed",
               linewidth = 0.5, alpha = 0.8) +
    geom_label(
      data          = gm_labels,
      aes(label     = label),
      x             = -Inf,
      y             = Inf,
      inherit.aes   = FALSE,
      hjust         = -0.1,
      vjust         = 2.3,
      size          = 5,
      fill          = "white",
      label.padding = unit(0.15, "lines"),
      label.size    = 0.3
    ) +
    facet_grid(
      rows     = vars(target_label),
      cols     = vars(antigen),
      labeller = labeller(antigen = antigen_labeller)
    ) +
    scale_fill_manual(values = c(diagonal = "#e972a7", off_diagonal = "#024a9c")) +
    theme_minimal() +
    theme(
      strip.text.x     = element_text(size = 20),
      strip.text.y     = element_text(size = 20),
      axis.line        = element_line(colour = "black", linewidth = 0.7),
      panel.grid       = element_blank(),
      panel.spacing    = unit(1, "lines"),
      legend.position  = "none",
      axis.text        = element_text(size = 20),
      axis.title.x     = element_blank(),
      axis.title.y     = element_blank(),
      panel.background = element_rect(fill = "#ffffff", colour = NA),
      plot.background  = element_rect(fill = "white", colour = NA),
      axis.ticks.x     = element_line(colour = "black", linewidth = 0.5),
      axis.ticks.y     = element_line(colour = "black", linewidth = 0.5),
      panel.border     = element_rect(colour = "black", fill = NA, linewidth = 0.3),
      aspect.ratio     = 1
    )
}



# testing gmean ratios 
ratio_df_IgG <- readRDS("Data/by_isotype/ratio_df_IgG.rds")

# DENV1 
denv1_pcr_pos <- ratio_df_IgG %>% filter(target == "DENV1")
nrow(denv1_pcr_pos)


denv1_pcr_pos_VLP <- Gmean(denv1_pcr_pos$DENV1_VLP, conf.level = 0.95)
log2(denv1_pcr_pos_VLP)
denv1_pcr_pos_NS1 <- Gmean(denv1_pcr_pos$DENV1_NS1, conf.level = 0.95)
log2(denv1_pcr_pos_NS1)
denv1_pcr_pos_D111 <- Gmean(denv1_pcr_pos$DENV1_DIII, conf.level = 0.95)
log2(denv1_pcr_pos_D111)
denv1_pcr_pos_SHERPADES <- Gmean(denv1_pcr_pos$SHERPADES_DENV1_DIII, conf.level = 0.95)
log2(denv1_pcr_pos_SHERPADES)


denv1_pcr_pos_denv2_titres <- Gmean(denv1_pcr_pos$SHERPADES_DENV2_DIII)
denv1_pcr_pos_denv2_titres

denv1_pcr_pos_denv3_titres <- Gmean(denv1_pcr_pos$SHERPADES_DENV3_DIII)
denv1_pcr_pos_denv3_titres

denv1_pcr_pos_denv4_titres <- Gmean(denv1_pcr_pos$SHERPADES_DENV4_DIII)
denv1_pcr_pos_denv4_titres


# ratio and relative ratios 
denv1_pcr_pos_denv1_titres  # homologus
denv1_pcr_pos_denv2_titres / denv1_pcr_pos_denv1_titres # heterlogous or off diagonal
denv1_pcr_pos_denv3_titres / denv1_pcr_pos_denv1_titres # heterlogous or off diagonal
denv1_pcr_pos_denv4_titres / denv1_pcr_pos_denv1_titres # heterlogous or off diagonal

#mean(ratio_pairs)


# logged ratio dfs
logged_ratio_df_IgG <- readRDS("Data/by_isotype/logged_ratio_df_IgG.rds")


logged_denv1_pcr_pos <- logged_ratio_df_IgG %>% filter(target == "DENV1")
logged_denv1_pcr_pos


logged_denv1_pcr_pos_denv1_titres <- mean(logged_denv1_pcr_pos$SHERPADES_DENV1_DIII)
logged_denv1_pcr_pos_denv1_titres

logged_denv1_pcr_pos_denv2_titres <- mean(logged_denv1_pcr_pos$SHERPADES_DENV2_DIII)
logged_denv1_pcr_pos_denv2_titres

logged_denv1_pcr_pos_denv3_titres <- mean(logged_denv1_pcr_pos$SHERPADES_DENV3_DIII)
logged_denv1_pcr_pos_denv3_titres

logged_denv1_pcr_pos_denv4_titres <- mean(logged_denv1_pcr_pos$SHERPADES_DENV4_DIII)
logged_denv1_pcr_pos_denv4_titres



# ratio and relative ratios 
logged_denv1_pcr_pos_denv1_titres
logged_denv1_pcr_pos_denv2_titres / logged_denv1_pcr_pos_denv1_titres
logged_denv1_pcr_pos_denv3_titres / logged_denv1_pcr_pos_denv1_titres
logged_denv1_pcr_pos_denv4_titres / logged_denv1_pcr_pos_denv1_titres


