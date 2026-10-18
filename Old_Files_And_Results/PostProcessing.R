#--- Script for plotting results ----

# -- longitudinal
# multivariate binomial 
longitudinal_binomial_modeling_results <- readRDS("longitudinal_binomial_modeling_results.rds")
chik_longitudinal_binomial_modeling_results <- readRDS("chik_longitudinal_binomial_modeling_results.rds")
# multivariate multinomial 
dengue_serotype_results <- readRDS("dengue_serotype_results.rds")
dengue_serotype_neg_results <- readRDS("dengue_serotype_neg_results.rds")
# univariate binomial 
flavi_univariate_results <- read_csv("flavi_univariate_results.csv")


#-- cross-sectional 
# multivariate binomial 
cross_sectional_binomial_modeling_results <- readRDS("cross_sectional_binomial_modeling_results.rds")
# multivariate multinomial 
cross_sectional_multinomial_modeling_dengue_serotype <- readRDS("cross_sectional_multinomial_modeling_dengue_serotype.rds")
cross_sectional_multinomial_modeling_dengue_serotype_neg <- readRDS("cross_sectional_multinomial_modeling_dengue_serotype_neg.rds")
# univariate binomial 
cross_sectional_univariate_results <- read_csv("cross_sectional_univariate_results.csv")


model_palette <- c(
  "GLMnet"                 = "#003554",
  "Random Forest"          = "#006494",
  "SVM"                    = "#00a6fb", 
  "NaiveBayes"             = "#90e0ef"
)

# --- Multivariate binomial plots 
plot_model_comparison_binomial <- function(results_df, target_name = NULL, metrics = c("AUROC" ,
                                                                             "AUPRC",
                                                                             "Brier")) {
  results_long <- results_df %>%
    pivot_longer(cols = all_of(metrics), names_to = "Metric", values_to = "Value") %>%
    filter(!is.na(Value))
  
  if (!is.null(target_name) && "Target" %in% colnames(results_long)) {
    results_long <- results_long %>% filter(Target == target_name)
    plot_title <- paste(target_name)
  } else {
    plot_title <- "Model Comparison"
  }
  
  ggplot(results_long, aes(x = Model, y = Value, fill = Model)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
    geom_text(aes(label = sprintf("%.3f", Value)), 
              position = position_dodge(width = 0.8), 
              vjust = -0.5, size = 4, fontface = "bold") +
    facet_wrap(~ Metric, scales = "free_y") +
    coord_cartesian(ylim = c(0, 1)) +
    scale_fill_manual(values = model_palette) +
    theme_minimal(base_size = 13) +
    theme(
      axis.text.x = element_text(angle = 30, hjust = 1, size = 11),
      axis.title = element_text(size = 13, face = "bold"),
      strip.text = element_text(size = 13, face = "bold"),
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      legend.position = "none"
    ) +
    labs(title = plot_title, x = "Model", y = "Value")
}

# --- Multivaraite multinomial plots 
plot_model_comparison_multinomial <- function(results_df, target_name = NULL, metrics = c("AUC_Micro" ,
                                                                             "AUPRC_Micro",
                                                                             "Brier", "StratBrier")) {
  results_long <- results_df %>%
    pivot_longer(cols = all_of(metrics), names_to = "Metric", values_to = "Value") %>%
    filter(!is.na(Value))
  
  if (!is.null(target_name) && "Target" %in% colnames(results_long)) {
    results_long <- results_long %>% filter(Target == target_name)
    plot_title <- paste(target_name)
  } else {
    plot_title <- "Model Comparison"
  }
  
  ggplot(results_long, aes(x = Model, y = Value, fill = Model)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
    geom_text(aes(label = sprintf("%.3f", Value)), 
              position = position_dodge(width = 0.8), 
              vjust = -0.5, size = 4, fontface = "bold") +
    facet_wrap(~ Metric, scales = "free_y") +
    coord_cartesian(ylim = c(0, 1)) +
    scale_fill_manual(values = model_palette) +
    theme_minimal(base_size = 13) +
    theme(
      axis.text.x = element_text(angle = 30, hjust = 1, size = 11),
      axis.title = element_text(size = 13, face = "bold"),
      strip.text = element_text(size = 13, face = "bold"),
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      legend.position = "none"
    ) +
    labs(title = plot_title, x = "Model", y = "Value")
}

# --- Longitudinal
# Binomial
plot_model_comparison_binomial(longitudinal_binomial_modeling_results$combined_comparison, target_name = "flavi")
plot_model_comparison_binomial(longitudinal_binomial_modeling_results$combined_comparison, target_name = "dengue")
plot_model_comparison_binomial(chik_longitudinal_binomial_modeling_results$combined_comparison, target_name = "dengue_chik")
# Multinomial
plot_model_comparison_multinomial(dengue_serotype_results$comparison)
plot_model_comparison_multinomial(dengue_serotype_neg_results$comparison)


# --- cross-sectional 
# Binomial
plot_model_comparison_binomial(cross_sectional_binomial_modeling_results$combined_comparison, target_name = "flavi")
plot_model_comparison_binomial(cross_sectional_binomial_modeling_results$combined_comparison, target_name = "dengue")
# Multinomial
plot_model_comparison_multinomial(cross_sectional_multinomial_modeling_dengue_serotype$comparison)
plot_model_comparison_multinomial(cross_sectional_multinomial_modeling_dengue_serotype_neg$comparison)


# --- Longitudinal : univariate
univariate_plot_data_longitudinal <- flavi_univariate_results %>%
  dplyr::select(Target, Variable, Model, AUC) %>%
  pivot_longer(cols = c(AUC), names_to = "Metric", values_to = "Value")
# Order Variable by mean AUC (descending)
univariate_plot_data_longitudinal <- univariate_plot_data_longitudinal %>%
  group_by(Variable) %>%
  mutate(mean_auc = mean(Value, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(Variable = reorder(Variable, mean_auc))
ggplot(univariate_plot_data_longitudinal, aes(x = Value, y = Variable, color = Model, shape = Model)) +
  geom_point(size = 3, alpha = 0.8) +
  facet_grid(Metric ~ Target, scales = "free_y") +
  scale_x_continuous(limits = c(0, 1)) +
  theme_minimal(base_size = 13) +
  theme(
    axis.text.y = element_text(size = 9),
    strip.text = element_text(size = 13, face = "bold"),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5)
  ) +
  labs(
    title = "Model Performance by Variable",
    x = "AUC",
    y = "Variable"
  )


# --- Cross-sectional : univariate
univariate_plot_data_cross_sectional <- cross_sectional_univariate_results %>%
  dplyr::select(Target, Variable, Model, AUC) %>%
  pivot_longer(cols = c(AUC), names_to = "Metric", values_to = "Value")

# Order Variable by mean AUC (descending)
univariate_plot_data_cross_sectional <- univariate_plot_data_cross_sectional %>%
  group_by(Variable) %>%
  mutate(mean_auc = mean(Value, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(Variable = reorder(Variable, mean_auc))

ggplot(univariate_plot_data_cross_sectional, aes(x = Value, y = Variable, color = Model, shape = Model)) +
  geom_point(size = 3, alpha = 0.8) +
  facet_grid(Metric ~ Target, scales = "free_y") +
  scale_x_continuous(limits = c(0, 1)) +
  theme_minimal(base_size = 13) +
  theme(
    axis.text.y = element_text(size = 9),
    strip.text = element_text(size = 13, face = "bold"),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5)
  ) +
  labs(
    title = "Model Performance by Variable (Cross-Sectional)",
    x = "AUC",
    y = "Variable"
  )

