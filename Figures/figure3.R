library(ggtext)
install.packages("DescTools")
library(DescTools)
library(ggforce)

ratio_df_IgG <- readRDS('Data/by_isotype/ratio_df_IgG.rds')
ratio_df_IgA <- readRDS('Data/by_isotype/ratio_df_IgA.rds')
ratio_df_IgM <- readRDS('Data/by_isotype/ratio_df_IgM.rds')
ratio_df_avidity <- readRDS('Data/by_isotype/ratio_df_avidity.rds')


pcr_target <- c("DENV1", "DENV2", "DENV3", "DENV4", "ZIKV")

dengue_zika_antigens <- c(
  "DENV1_DIII", "DENV1_VLP", "DENV1_NS1",
  "DENV2_DIII", "DENV2_VLP", "DENV2_NS1",
  "DENV3_DIII", "DENV3_VLP", "DENV3_NS1",
  "DENV4_DIII", "DENV4_VLP", "DENV4_NS1",
  "ZIKVAS_DIII", "ZIKV_VLP", "ZIKV_NS1", "ZIKVSU_NS1",
  "SHERPADES_DENV1_DIII", "SHERPADES_DENV2_DIII", 
  "SHERPADES_DENV3_DIII", "SHERPADES_DENV4_DIII",
  "SHERPADES_ZIKV_DIII"
)


#   log_gm_rel = mean( log2( b_i / a_i ) )  ==  log2( GM(B) / GM(A) )
#
#   where  a = homologous antigen (matching the infecting serotype)  [reference]
#          b = the antigen of the current panel                     [test]
#
#   These two forms are algebraically identical because the mean of differences
#   equals the difference of means:
#       mean(log2 b - log2 a) = mean(log2 b) - mean(log2 a)
#                             = log2(GM(B)) - log2(GM(A))
#                             = log2(GM(B) / GM(A))
#   This holds for paired (element-wise) ratios and for all n_a x n_b pairwise
#   combinations alike -- the point estimate is the same. We use the PAIRED form
#   because the resulting confidence interval is the correct within-subject one.
# ================================================


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

  fmt_antigen_multiline <- function(x) dplyr::case_when(
  str_detect(x, "^SHERPADES_") ~ str_replace(x, "^SHERPADES_([^_]+)_DIII$", "SHERPADES\n\\1\nDIII"),
  str_detect(x, "^ZIKVAS_")    ~ "ZIKV\nDIII",
  TRUE                         ~ str_replace(x, "_", "\n")
)
  
  gmean <- function(x, conf_level = 0.95) {
    x <- x[is.finite(x) & x > 0]                   # GM is undefined for <= 0
    if (length(x) == 0) {
      return(c(mean = NA_real_, lwr.ci = NA_real_, upr.ci = NA_real_))
    }
    out <- DescTools::Gmean(x, conf.level = conf_level)
    if (is.null(names(out))) {
      return(c(mean = unname(out[1]), lwr.ci = NA_real_, upr.ci = NA_real_))
    }
    c(mean   = unname(out[["mean"]]),
      lwr.ci = unname(out[["lwr.ci"]]),
      upr.ci = unname(out[["upr.ci"]]))
  }
  
  # log2 GMR + t-based CI, computed on the log2 scale. Returns exactly ONE row.
  log2_gmr <- function(ratios, conf.level =  0.95) {
    ratios <- ratios[is.finite(ratios) & ratios > 0]
    n <- length(ratios)
    if (n == 0) {
      return(tibble(n = 0L, log_gm_rel = NA_real_,
                    gm_lower = NA_real_, gm_upper = NA_real_))
    }
    lr <- log2(ratios)
    m  <- mean(lr)
    s  <- stats::sd(lr)
    # n < 2, or zero variance (e.g. the diagonal, where every ratio is 1)
    if (n < 2 || !is.finite(s) || s == 0) {
      return(tibble(n = n, log_gm_rel = m,
                    gm_lower = NA_real_, gm_upper = NA_real_))
    }
    se <- s / sqrt(n)
    tc <- stats::qt(1 - (1 - conf.level) / 2, df = n - 1)
    tibble(n          = n,
           log_gm_rel = m,
           gm_lower   = m - tc * se,
           gm_upper   = m + tc * se)
  }
 


# prepare data and plot 
prep_data <- function(data, 
                      antigen, 
                      flavi_targets = pcr_target,
                      antigen_pool = dengue_zika_antigens,
                      conf_level = 0.95) {
if (antigen == "SHERPADES") {
    filtered_antigens <- antigen_pool[str_detect(antigen_pool, "^SHERPADES_.*_DIII$")]
    target_to_col <- setNames(
      vapply(flavi_targets, get_antigen_col, character(1), suffix = "SHERPADES"),
      flavi_targets
    )
  } else {
    filtered_antigens <- antigen_pool[
      str_detect(antigen_pool, paste0(antigen, "$")) &
        !str_detect(antigen_pool, "^SHERPADES_")
    ]
    target_to_col <- setNames(
      vapply(flavi_targets, function(t) get_antigen_col(t, antigen), character(1)),
      flavi_targets
    )
  }
  if (length(filtered_antigens) == 0) {
    stop("No antigens matched antigen = '", antigen, "'.\n  Available: ",
         paste(antigen_pool, collapse = ", "))
  }
  missing_cols <- setdiff(unique(c(filtered_antigens, unname(target_to_col))),
                          names(data))
  if (length(missing_cols) > 0) {
    stop("These expected columns are not present in `data`: ",
         paste(missing_cols, collapse = ", "))
  }

plot_data_list <- list()
mean_data_list <- list()
forest_data_list <- list() 


# Map each pcr_target to its actual column name in filtered_antigens
for (infecting_pathogen in flavi_targets) {
  
  diagonal_data  <- data %>% filter(target == infecting_pathogen)
  infecting_col  <- target_to_col[[infecting_pathogen]]
  
  for (current_antigen in filtered_antigens) {
    key <- paste(infecting_pathogen, current_antigen, sep = "_")
    antigen_pathogen <- get_antigen_pathogen(current_antigen)
    is_diagonal <- identical(antigen_pathogen, infecting_pathogen)
    panel_type <- if (is_diagonal) "diagonal" else "off_diagonal"

    if (is_diagonal){
      vals <- diagonal_data[[infecting_col]] # homologous - absolute titres 
    } else {
      vals <- diagonal_data[[current_antigen]] / diagonal_data[[infecting_col]] # relative ratio
    }
    
      
    plot_data_list[[key]] <- tibble(
      value            = vals,
      infecting_target = infecting_pathogen,
      antigen          = current_antigen,
      panel_type       = panel_type
    )


    # Geometric mean
    gm_value <- gmean(vals, conf_level = conf_level)
    
    mean_data_list[[key]] <- tibble(
        infecting_target = infecting_pathogen,
        antigen          = current_antigen,
        panel_type       = panel_type,
        n                = length(vals),
        gm               = gm_value[["mean"]],
        gm_lwr           = gm_value[["lwr.ci"]],
        gm_upr           = gm_value[["upr.ci"]],
        label            = round(gm_value[["mean"]], 2)
      )

    ratios <- diagonal_data[[current_antigen]] / diagonal_data[[infecting_col]]

    forest_data_list[[key]] <- log2_gmr(ratios) %>%
      mutate(
        infecting_target = infecting_pathogen,
        antigen          = current_antigen,
        panel_type       = panel_type) %>%
      dplyr::select(infecting_target, antigen, panel_type,
                    n, log_gm_rel, gm_lower, gm_upper)}}

  antigen_levels <- fmt_antigen(filtered_antigens)
  target_levels  <- flavi_targets


  add_labels <- function(df) {
    df %>% mutate(
      antigen_label = factor(fmt_antigen(antigen), levels = antigen_levels),
      target_label  = factor(infecting_target,  levels = target_levels)
    )
  }
 
  plot_data <- bind_rows(plot_data_list) %>% add_labels()
  mean_data <- bind_rows(mean_data_list) %>% add_labels()

  
  forest_data <- bind_rows(forest_data_list) %>%
    dplyr::filter(!is.na(log_gm_rel)) %>%
    add_labels() %>%
    mutate(
      # homologous sits at exactly 0 on the log2 scale, by construction
      label_text = if_else(
        panel_type == "diagonal",
        "0",
        sprintf("%.2f (%.2f - %.2f)", log_gm_rel, gm_lower, gm_upper)
      )
    )
  list(
    plot_data    = plot_data,
    mean_data    = mean_data,
    forest_data  = forest_data,
    antigen_cols = filtered_antigens,
    conf_level   = conf_level
  )
}


# ---- Plots 
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

antigen_types <- c("VLP", "NS1", "DIII", "SHERPADES")



ratio_dfs <- list(
  IgG     = ratio_df_IgG,
  IgA     = ratio_df_IgA,
  IgM     = ratio_df_IgM,
  avidity = ratio_df_avidity
)

prep_data_list <- list()
histogram_list <- list()
forest_list <- list()

dir.create("Results/Fig3", showWarnings = FALSE)

for (iso in names(ratio_dfs)) {
  for (antigen in antigen_types) {

    key <- paste(iso, antigen, sep = "_")   # composite key: nothing overwrites
    message("Processing: ", key)

    res <- prep_data(ratio_dfs[[iso]], antigen = antigen)
    hist_p  <- plot_ratio_histogram(res)
    forest_p <- plot_ratio_forest(res)

    prep_data_list[[key]] <- res
    histogram_list[[key]] <- hist_p
    forest_list[[key]]    <- forest_p

    n_row <- dplyr::n_distinct(res$plot_data$target_label)
    n_col <- dplyr::n_distinct(res$plot_data$antigen)
    tag   <- gsub("[^A-Za-z0-9_-]", "_", key)

    ggsave(file.path("Results/Fig3", sprintf("hist_%s.png", tag)), hist_p,
           width = n_col * 2.2 + 1, height = n_row * 2.2 + 1,
           units = "in", limitsize = FALSE)

    ggsave(file.path("Results/Fig3", sprintf("forest_%s.png", tag)), forest_p,
           width = 13, height = 16, units = "in", limitsize = FALSE)
  }
}




ratio_df_IgG <- readRDS("Data/by_isotype/ratio_df_IgG.rds")

View(ratio_df_IgG)

denv1_pcr_pos <- ratio_df_IgG %>% filter(target == "DENV1")
denv1_pcr_pos

denv1_pcr_pos_denv1_titres <- Gmean(denv1_pcr_pos$DENV1_VLP, conf.level = 0.95)
denv1_pcr_pos_denv1_titres

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



heatmap_data <- mean_data %>%
    mutate(log_odds = log2(gm))

# heat map 
heatmap <- ggplot(heatmap_data, aes(x = antigen, y = infecting_target, fill = log_odds)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = round(log_odds, 1)), size = 3.5, fontface = "bold") +
    scale_fill_gradient2(
      low      = "#024a9c",
      mid      = "white",
      high     = "#e972a7",
      midpoint = 0,
      limits   = c(-5, 5),
      name     = "log2 odds\n(up/down)"
    ) +
    scale_x_discrete(
      position = "top",
      labels = function(x) case_when(
      str_detect(x, "^SHERPADES_") ~ str_replace(x, "^SHERPADES_([^_]+)_DIII$", "SHERPADES\n\\1 DIII"),
      str_detect(x, "^ZIKVAS_")    ~ "ZIKV\nDIII",
      TRUE                         ~ str_replace(x, "_", "\n")
    )) + 
    scale_y_discrete(limits = rev) +
    labs(x = "Antigen", y = "Infecting pathogen") +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x      = element_text(size = 16),
      axis.text.y      = element_text(size = 16),
      axis.title       = element_text(size = 16),
      legend.title     = element_text(size = 16),
      legend.text      = element_text(size = 16),
      panel.grid       = element_blank(),
      plot.background  = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA)
    )





# This function computes summary statistics for ZIKV cross-reactivity ratios
# across all DENV-infected samples (DENV1-4).
# For each DENV-infected group, we calculate the ratio:
#   ZIKV antigen signal / infecting DENV antigen signal
# This tells us how strongly DENV-infected patients' antibodies
# cross-react with ZIKV antigen, relative to their own infecting virus.
# All ratios from DENV1, DENV2, DENV3, and DENV4 infected groups are
# pooled together, since there is only one ZIKV-infected sample —
# not enough to compute meaningful within-group statistics.
# Output: mean, IQR, min/max, and proportion of ratios above/below 1
#   ratio < 1 → weaker reaction to ZIKV than to infecting DENV serotype
#   ratio = 1 → equal cross-reactivity
#   ratio > 1 → stronger reaction to ZIKV than to infecting DENV serotype



compute_zikv_column_stats <- function(data, antigen_suffix, flavi_targets = pcr_target) {
  
  all_zikv_ratios <- c()
       # Map antigen names to column prefixes for specific suffixes
    antigen_col_map <- list(
      DIII = c(DENV1 = "DENV1", DENV2 = "DENV2", DENV3 = "DENV3", DENV4 = "DENV4", ZIKV = "ZIKVAS")
    )

    get_col_name <- function(antigen, suffix) {
    if (suffix %in% names(antigen_col_map) && antigen %in% names(antigen_col_map[[suffix]])) {
      paste0(antigen_col_map[[suffix]][[antigen]], "_", suffix)
    } else {
      paste0(antigen, "_", suffix)
    }
  }
  
  # Loop through all non-ZIKV infecting targets (off-diagonal ZIKV column)
  for (infecting_target in flavi_targets[flavi_targets != "ZIKV"]) {
    
    target_data <- data %>% filter(target == infecting_target)
    
    infecting_col <- get_col_name(infecting_target, antigen_suffix)
    zikv_col      <- get_col_name("ZIKV", antigen_suffix)  # "ZIKVAS_DIII"
    
    valid_rows <- complete.cases(target_data[, c(infecting_col, zikv_col)])
    
    if (sum(valid_rows) > 0) {
      relative_ratios <- target_data[valid_rows, zikv_col] / target_data[valid_rows, infecting_col]
      all_zikv_ratios <- c(all_zikv_ratios, relative_ratios)
    }
  }
  
  # --- Summary statistics ---
  if (length(all_zikv_ratios) > 0) {
    
    overall_mean <- mean(all_zikv_ratios, na.rm = TRUE)
    overall_q25  <- quantile(all_zikv_ratios, 0.25, na.rm = TRUE)
    overall_q75  <- quantile(all_zikv_ratios, 0.75, na.rm = TRUE)
    overall_iqr  <- overall_q75 - overall_q25
    
    below_1 <- sum(all_zikv_ratios < 1, na.rm = TRUE)
    above_1 <- sum(all_zikv_ratios > 1, na.rm = TRUE)
    n       <- length(all_zikv_ratios)
    
    cat("\nOVERALL LAST COLUMN STATISTICS (All DENV-ZIKV comparisons):\n")
    cat(sprintf("Total samples: %d\n",        n))
    cat(sprintf("Mean: %.3f\n",               overall_mean))
    cat(sprintf("IQR: %.3f - %.3f (range: %.3f)\n", overall_q25, overall_q75, overall_iqr))
    cat(sprintf("Min: %.3f\n",                min(all_zikv_ratios)))
    cat(sprintf("Max: %.3f\n",                max(all_zikv_ratios)))
    cat(sprintf("Values < 1.0: %d/%d (%.1f%%)\n", below_1, n, 100 * below_1 / n))
    cat(sprintf("Values > 1.0: %d/%d (%.1f%%)\n", above_1, n, 100 * above_1 / n))
  }
}

compute_zikv_column_stats(ratio_df, "VLP")
compute_zikv_column_stats(ratio_df, "NS1")
compute_zikv_column_stats(ratio_df, "DIII")