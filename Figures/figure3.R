library(ggtext)

preprocessed_data_raw <- read.csv('Results/raw_preprocessed_cebu_data.csv')
View(preprocessed_data_raw)
ratio_df <- readRDS('Results/ratio_df.rds')
ratio_df_logged <- readRDS('Results/logged_ratio_df.rds')
View(ratio_df)
View(ratio_df_logged)

pcr_target <- c("DENV1", "DENV2", "DENV3", "DENV4", "ZIKV")

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

create_custom_bins <- function(data, min_val = 0.75, max_val = 1.25, n_middle_bins = 3) {
 

  data <- data[!is.na(data)]
  
  middle_edges <- seq(min_val, max_val, length.out = n_middle_bins + 1)
  
  counts <- c()
  bin_centers <- c()
  bin_names <- c()
  
  low_count <- sum(data <= min_val)
  counts <- c(counts, low_count)
  bin_centers <- c(bin_centers, log2(min_val))
  bin_names <- c(bin_names, "low")
  
  for (i in 1:(length(middle_edges) - 1)) {
    mask <- (data > middle_edges[i]) & (data <= middle_edges[i + 1])
    counts <- c(counts, sum(mask))
    center_val <- (middle_edges[i] + middle_edges[i + 1]) / 2
    bin_centers <- c(bin_centers, log2(center_val))
    bin_names <- c(bin_names, paste0("mid_", i))
  }
  
  high_count <- sum(data >= max_val)
  counts <- c(counts, high_count)
  bin_centers <- c(bin_centers, log2(max_val))
  bin_names <- c(bin_names, "high")
  
  data.frame(
    bin = bin_names,
    count = counts,
    bin_center = bin_centers
  )
}

ratio_plot_data <- function(data, antigen, flavi_targets = pcr_target) {
  
  plot_data_list <- list()
  mean_data_list <- list()

  get_antigen_col <- function(pathogen, suffix) {
    if (suffix == "DIII" && pathogen == "ZIKV") {
      return("ZIKVAS_DIII")
    }
    if (suffix == "SHERPADES") {
      return(paste0("SHERPADES_", pathogen, "_DIII"))
    }
    return(paste0(pathogen, "_", suffix))
  }

  get_antigen_pathogen <- function(col_name) {
    if (str_detect(col_name, "^SHERPADES_")) {
      return(str_extract(col_name, "(?<=^SHERPADES_)[^_]+"))
    }
    if (str_detect(col_name, "^ZIKVAS_")) {
      return("ZIKV")
    }
    return(str_extract(col_name, "^[^_]+"))
  }

  if (antigen == "SHERPADES") {
    filtered_antigens <- flavivirus_antigens[str_detect(flavivirus_antigens, "^SHERPADES_.*_DIII$")]
    target_to_col <- setNames(
      sapply(flavi_targets, get_antigen_col, suffix = "SHERPADES"),
      flavi_targets
    )
  } else {
    filtered_antigens <- flavivirus_antigens[
      str_detect(flavivirus_antigens, paste0(antigen, "$")) &
      !str_detect(flavivirus_antigens, "^SHERPADES_")
    ]
    target_to_col <- setNames(
      sapply(flavi_targets, function(t) get_antigen_col(t, antigen)),
      flavi_targets
    )
  }

   # Map each pcr_target to its actual column name in filtered_antigens

  for (infecting_pathogen in flavi_targets) {
    
    diagonal_data    <- data %>% filter(target == infecting_pathogen)
    off_diagonal_data <- data %>% filter(target != infecting_pathogen)
    infecting_col  <- target_to_col[[infecting_pathogen]]
    
    for (current_antigen in filtered_antigens) {
      antigen_pathogen <- get_antigen_pathogen(current_antigen)
      
      if (antigen_pathogen == infecting_pathogen) {
        # Diagonal
        ratios     <- diagonal_data[[infecting_col]]
        ratios     <- ratios[!is.na(ratios)]
        bins_df    <- create_custom_bins(ratios, 0.75, 1.25, 3)
        panel_type <- "diagonal"
        
      } else {
        # Off-diagonal
        non_infecting_col <- current_antigen
        relative_ratios   <- off_diagonal_data[[non_infecting_col]] / off_diagonal_data[[infecting_col]]
        relative_ratios   <- relative_ratios[!is.na(relative_ratios)]
        bins_df           <- create_custom_bins(relative_ratios, 0.75, 1.25, 3)
        panel_type        <- "off_diagonal"
      }
      
      # Store plot data
      plot_data_list[[paste(infecting_pathogen, current_antigen, sep = "_")]] <- bins_df %>%
        mutate(
          infecting_target = infecting_pathogen,
          antigen          = current_antigen,
          panel_type       = panel_type
        )
      
      # Store mean data
      mean_data_list[[paste(infecting_pathogen, current_antigen, sep = "_")]] <- tibble(
        infecting_target = infecting_pathogen,
        antigen          = current_antigen,
        panel_type       = panel_type,
        x                = mean(if (panel_type == "diagonal") ratios else relative_ratios, na.rm = TRUE),
        y                = max(bins_df$count),
        label            = round(x, 2)
      )
    }
  }
  
  list(
    plot_data = bind_rows(plot_data_list),
    mean_data = bind_rows(mean_data_list)
  )
}

plot_ratio_grid <- function(data, antigen_suffix, flavi_targets = pcr_target) {
  
  plot_inputs <- ratio_plot_data(data, antigen_suffix, flavi_targets)
  
  antigen_labeller <- as_labeller(function(x) {
  case_when(
    str_detect(x, "^SHERPADES_") ~ str_replace(x, "^SHERPADES_([^_]+)_DIII$", "SHERPADES\n\\1\nDIII"),
    str_detect(x, "^ZIKVAS_")    ~ "ZIKV\nDIII",
    TRUE                         ~ str_replace(x, "_", "\n")
  )
})

plot_data <- plot_inputs$plot_data
  
  ggplot(plot_data, aes(x = bin_center, y = count, fill = panel_type)) +
    geom_col(width = 0.4, alpha = 0.8) +
    # Add vertical line at x=0 (ratio = 1 on log2 scale) for reference
    geom_vline(xintercept = 0, color = "red", linetype = "dashed", linewidth = 0.5, alpha = 0.8) +
    geom_label(
        data = plot_inputs$mean_data %>% filter(!is.na(label)),
        aes(x = log2(x), y = y, label = label),
        x = -Inf,        # left edge of every panel
        y = Inf,
        inherit.aes = FALSE,
        hjust = 0,
        vjust = 1.3,
        size = 5,
        fill = "white",        
        label.padding = unit(0.15, "lines"),  
        label.size = 0.3       
        ) +
    facet_grid(
      rows = vars(infecting_target),
      cols = vars(antigen),
      labeller = labeller(antigen = antigen_labeller)) +
    scale_fill_manual(
      values = c(
        diagonal = "#e972a7",
        off_diagonal = "#024a9c"
      )
    ) + scale_x_continuous( # using log2 scale 
      limits = c(-2.2, 2.2),
      breaks = c(-2, -1, 0, 1, 2),
      labels = c("1/4", "1/2", "1", "2", "4")
    ) +
    theme_minimal(base_size = 10) +
    theme(
      strip.text.x = element_text(size = 20),
      strip.text.y = element_text(size = 20),
      axis.line = element_line(color = "black", linewidth = 0.7),
      panel.grid = element_blank(),
      panel.spacing = unit(0.1, "lines"),
      legend.position = "none",
      axis.text = element_text(size = 20),
      axis.text.x = element_text(size = 20, hjust = 0.5),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      panel.background = element_rect(fill = "#ffffff", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      axis.ticks.x = element_line(color = "black", size = 0.5),
      axis.ticks.y = element_line(color = "black", size = 0.5),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.3),
      aspect.ratio = 1

    )
}


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

vlp_ratio_plot <- plot_ratio_grid(ratio_df_logged, antigen_suffix = "VLP")
ns1_ratio_plot <- plot_ratio_grid(ratio_df_logged, antigen_suffix = "NS1")
DIII_ratio_plot <- plot_ratio_grid(ratio_df_logged, antigen_suffix = "DIII")
sherpades_DIII_ratio_plot <- plot_ratio_grid(ratio_df_logged, antigen_suffix = "SHERPADES")

vlp_ratio_plot
ns1_ratio_plot
DIII_ratio_plot
sherpades_DIII_ratio_plot




# save VLP 
ggsave("Results/vlp_ratio_plot.png", 
vlp_ratio_plot, width = 12, height = 8)

# save NS1
ggsave("Results/ns1_ratio_plot.png", 
ns1_ratio_plot, width = 12, height = 8)

# save DIII
ggsave("Results/DIII_ratio_plot.png", 
DIII_ratio_plot, width = 12, height = 8)


# save SHERPADES DIII
ggsave("Results/sherpades_DIII_ratio_plot.png", 
sherpades_DIII_ratio_plot, width = 12, height = 8)



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

compute_zikv_column_stats(ratio_df, "VLP")
compute_zikv_column_stats(ratio_df, "NS1")
compute_zikv_column_stats(ratio_df, "DIII")
