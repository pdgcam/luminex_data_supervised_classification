ratio_df <- readRDS('Results/ratio_df.rds')


pcr_target <- c("DENV1", "DENV2", "DENV3", "DENV4", "ZIKV")
antigen_labels_map <- c(
  DENV1 = "IgG DENV-1",
  DENV2 = "IgG DENV-2",
  DENV3 = "IgG DENV-3",
  DENV4 = "IgG DENV-4",
  ZIKV  = "IgG ZIKV"
)

create_custom_bins <- function(data, min_val = 0.75, max_val = 1.25, n_middle_bins = 3) {
  if (length(data) == 0) {
    return(data.frame(
      bin = character(),
      count = numeric(),
      bin_center = numeric()
    ))
  }
  
  data <- data[!is.na(data)]
  
  if (length(data) == 0) {
    return(data.frame(
      bin = character(),
      count = numeric(),
      bin_center = numeric()
    ))
  }
  
  middle_edges <- seq(min_val, max_val, length.out = n_middle_bins + 1)
  
  counts <- c()
  bin_centers <- c()
  bin_names <- c()
  
  low_count <- sum(data <= min_val)
  counts <- c(counts, low_count)
  bin_centers <- c(bin_centers, log2(min_val - 0.1))
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
  bin_centers <- c(bin_centers, log2(max_val + 0.1))
  bin_names <- c(bin_names, "high")
  
  data.frame(
    bin = bin_names,
    count = counts,
    bin_center = bin_centers
  )
}

prepare_ratio_plot_data <- function(data, antigen_suffix, flavi_targets = pcr_target) {
  plot_data_list <- list()
  mean_data_list <- list()
  
  for (infecting_target in flavi_targets) {
    for (antigen in flavi_targets) {
      
      target_data <- data %>%
        filter(target == infecting_target)
      
      if (infecting_target == antigen) {
        antigen_col <- paste0(antigen, "_", antigen_suffix)
        ratios <- target_data[[antigen_col]]
        ratios <- ratios[!is.na(ratios)]
        
        bins_df <- create_custom_bins(ratios, 0.75, 1.25, 3)
        
        if (nrow(bins_df) > 0) {
          bins_df <- bins_df %>%
            mutate(
              infecting_target = infecting_target,
              antigen = antigen,
              panel_type = "diagonal"
            )
          plot_data_list[[length(plot_data_list) + 1]] <- bins_df
        }
        
        mean_data_list[[length(mean_data_list) + 1]] <- data.frame(
          infecting_target = infecting_target,
          antigen = antigen,
          label = paste0("Mean: ", sprintf("%.2f", mean(ratios, na.rm = TRUE)))
        )
        
      } else {
        infecting_col <- paste0(infecting_target, "_", antigen_suffix)
        non_infecting_col <- paste0(antigen, "_", antigen_suffix)
        
        valid_rows <- complete.cases(target_data[, c(infecting_col, non_infecting_col)])
        
        if (sum(valid_rows) > 0) {
          relative_ratios <- target_data[valid_rows, non_infecting_col] / target_data[valid_rows, infecting_col]
          
          bins_df <- create_custom_bins(relative_ratios, 0.75, 1.25, 3)
          
          if (nrow(bins_df) > 0) {
            bins_df <- bins_df %>%
              mutate(
                infecting_target = infecting_target,
                antigen = antigen,
                panel_type = "off_diagonal"
              )
            plot_data_list[[length(plot_data_list) + 1]] <- bins_df
          }
          
          mean_data_list[[length(mean_data_list) + 1]] <- data.frame(
            infecting_target = infecting_target,
            antigen = antigen,
            label = paste0("Mean: ", sprintf("%.2f", mean(relative_ratios, na.rm = TRUE)))
          )
        } else {
          mean_data_list[[length(mean_data_list) + 1]] <- data.frame(
            infecting_target = infecting_target,
            antigen = antigen,
            label = NA_character_
          )
        }
      }
    }
  }
  
  plot_data <- bind_rows(plot_data_list)
  mean_data <- bind_rows(mean_data_list)
  
  plot_data$infecting_target <- factor(plot_data$infecting_target, levels = flavi_targets)
  plot_data$antigen <- factor(plot_data$antigen, levels = flavi_targets)
  
  mean_data$infecting_target <- factor(mean_data$infecting_target, levels = flavi_targets)
  mean_data$antigen <- factor(mean_data$antigen, levels = flavi_targets)
  
  mean_data <- mean_data %>%
    mutate(
      x = -2.1,
      y = Inf
    )
  
  list(plot_data = plot_data, mean_data = mean_data)
}

plot_ratio_grid <- function(data,
                            antigen_suffix,
                            flavi_targets = pcr_target,
                            antigen_labels = antigen_labels_map) {
  
  plot_inputs <- prepare_ratio_plot_data(data, antigen_suffix, flavi_targets)
  
  ggplot(plot_inputs$plot_data, aes(x = bin_center, y = count, fill = panel_type)) +
    geom_col(width = 0.2, alpha = 0.7, color = "black", linewidth = 0.3) +
    geom_vline(xintercept = 0, color = "red", linetype = "dashed", linewidth = 0.5, alpha = 0.8) +
    geom_text(
      data = plot_inputs$mean_data %>% filter(!is.na(label)),
      aes(x = x, y = y, label = label),
      inherit.aes = FALSE,
      hjust = 0,
      vjust = 1.1,
      size = 4
    ) +
    facet_grid(
      rows = vars(infecting_target),
      cols = vars(antigen),
      labeller = labeller(antigen = antigen_labels)
    ) +
    scale_fill_manual(
      values = c(
        diagonal = "#FF6B6B",
        off_diagonal = "#4A90E2"
      )
    ) +
    scale_x_continuous(
      limits = c(-2.2, 2.2),
      breaks = c(-2, -1, 0, 1, 2),
      labels = c("1/4", "1/2", "1", "2", "4")
    ) +
    labs(
      x = "Ratio relative to infecting serotype",
      y = NULL
    ) +
    theme_minimal(base_size = 10) +
    theme(
      strip.text.x = element_text(size = 20),
      strip.text.y = element_text(size = 20),
      axis.line = element_line(color = "black", linewidth = 0.7),
      panel.grid = element_blank(),
      panel.spacing = unit(0.15, "lines"),
      legend.position = "none",
      axis.text = element_text(size = 20),
      axis.title.x = element_text(size = 20),
      panel.background = element_rect(fill = "#ffffff", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      axis.ticks.x = element_line(color = "black", size = 0.5),
      axis.ticks.y = element_line(color = "black", size = 0.5),
      panel.border = element_blank()
    )
}


vlp_ratio_plot <- plot_ratio_grid(ratio_df, antigen_suffix = "VLP")
# save VLP 
ggsave("Results/vlp_ratio_plot.png", 
vlp_ratio_plot, width = 10, height = 10)

plot_ratio_grid(ratio_df, antigen_suffix = "NS1")
plot_ratio_grid(ratio_df, antigen_suffix = "EDIII")

