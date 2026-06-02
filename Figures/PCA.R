library(ggrepel)


# marco PCA results
arb_pca_model <- readRDS('/Users/ap2488/Desktop/supervised_learning_flavi/MarcoCode/my_luminex_arb_pca_model.rds')
pca_model <- readRDS('/Users/ap2488/Desktop/supervised_learning_flavi/MarcoCode/my_luminex_pca_model.rds')

names(arb_pca_model)
names(pca_model)
names(pca_model$pca)
head(pca_model$pca$sdev)
head(pca_model$pca$rotation)
head(pca_model$pca$center)
head(pca_model$pca$scale)
head(pca_model$pca$x)
arb_pca_model$pca$x

pca_model$center
pca_model$scale

arb_pca_model$pca
arb_pca_model$center
arb_pca_model$scale


# combine scores back with target/id_patient for later overlay
pca_scores_marco <- as.data.frame(arb_pca_model$pca$x) 
quartz()
plot(arb_pca_model$pca$x[,1], arb_pca_model$pca$x[,2], col = "blue", pch = 16, xlab = "PC1", ylab = "PC2", main = "PCA Scores")

plot_cols <- c(
  DENV = "#8a0101",
  ZIKV = "#0d3392",
  CHIKV = "#008080",
  negative = "#747272"
)

names(arb_pca_model)

# --- Loadings (variables) ---
loadings <- as.data.frame(arb_pca_model$rotation[, 1:2])
loadings$variable <- rownames(loadings)

# Scale loadings to scores range for overlay
scale_factor <- max(abs(pca_scores_marco[, c("PC1","PC2")])) / max(abs(loadings[, 1:2])) * 0.7
loadings[, 1:2] <- loadings[, 1:2] * scale_factor

loadings$magnitude <- sqrt(loadings$PC1^2 + loadings$PC2^2)


# --- Variance explained ---
var_exp <- round(summary(arb_pca_model)$importance[2, 1:2] * 100, 1)
quartz()
# --- Plot ---
ggplot(pca_scores_marco, aes(x = PC1, y = PC2, colour = target)) +
  geom_point(alpha = 0.7, size = 2) +
  # Draw all arrows, label only top ones
  geom_segment(data = loadings,
               aes(x = 0, y = 0, xend = PC1, yend = PC2),
               arrow = arrow(length = unit(0.15, "cm")),
               colour = "grey50", linewidth = 0.4,
               inherit.aes = FALSE) +
  geom_text_repel(data = loadings,
                  aes(x = PC1, y = PC2, label = variable),
                  colour = "grey20", size = 3,
                  max.overlaps = Inf,
                  box.padding = 0.4,
                  segment.colour = "grey60",
                  segment.size = 0.3,
                  inherit.aes = FALSE) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey70") +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey70") +
  theme_minimal(base_size = 13)
quartz()
print(p)













ratio_df <- readRDS('Results/ratio_df.rds')

flavi_antigens <- c("DENV1_DIII","DENV1_NS1","DENV1_VLP","SHERPADES_DENV1_DIII",
"DENV2_DIII","DENV2_NS1","DENV2_VLP","SHERPADES_DENV2_DIII",
"DENV3_DIII","DENV3_NS1","DENV3_VLP", "SHERPADES_DENV3_DIII",
"DENV4_DIII","DENV4_NS1","DENV4_VLP", "SHERPADES_DENV4_DIII",
"JEV_E", "JEV_NS1", "SHERPADES_JEV_DIII",
"YFV_E", "YFV_NS1", "SHERPADES_YFV_DIII",
"WNV_DIII","WNV_NS1","SHERPADES_WNV_DIII",
"ZIKV_NS1","ZIKV_VLP","ZIKVAS_DIII","ZIKVSU_NS1","SHERPADES_ZIKV_DIII")

alpha_antigens <- c("CHIKV_E2", "CHIKV_NSP123", "CHIKV_VLP", "SHERPADES_CHIKV_E2", 
                    "MAYV_E2" , "SHERPADES_MAYV_E2",
                    "ONNV_E2", "ONNV_VLP",
                    "RR" , "SHERPADES_RR")
     


all_antigens <- c(flavi_antigens, alpha_antigens)


# subset flavi antigens 
ratio_df_arbo <- ratio_df %>%
  dplyr::select(target, id_patient, all_of(c(flavi_antigens, alpha_antigens)))


# scale only the antigen columns and keep the metadata columns unchanged
feature_matrix <- ratio_df_arbo %>%
   dplyr::select(all_of(c(flavi_antigens, alpha_antigens)))
scaled_features <- scale(feature_matrix)

# run PCA on the scaled antigen data only
pca <- prcomp(scaled_features)
summary(pca)
quartz()
biplot(pca)


# Get the scaling parameters from your original scale() call
scale_center <- attr(scaled_features, "scaled:center")
scale_scale <- attr(scaled_features, "scaled:scale")

# combine scores back with target/id_patient for later overlay
pca_scores <- as.data.frame(pca$x) %>%
  dplyr::bind_cols(ratio_df_arbo %>% dplyr::select(target, id_patient))



# Save everything to share with Anchita
saveRDS(list(
pca = pca, center = scale_center, scale = scale_scale
), file = here("Results/cebu_pca.rds"))




# Plot the first two principal component scores, grouping by DENV, ZIKV, CHIKV, negative
pca_scores <- pca_scores %>%
  dplyr::mutate(
    target_group = dplyr::case_when(
      grepl("^DENV[1-4]", target, ignore.case = TRUE) ~ "DENV",
      grepl("^ZIKV", target, ignore.case = TRUE) ~ "ZIKV",
      grepl("^CHIKV", target, ignore.case = TRUE) ~ "CHIKV",
      target %in% c("negative", "NEGATIVE", "Negative", "neg") ~ "negative",
      TRUE ~ "other"
    ),
    target_group = factor(target_group, levels = c("DENV", "ZIKV", "CHIKV", "negative", "other"))
  )

plot_cols <- c(
  DENV = "#8a0101",
  ZIKV = "#0d3392",
  CHIKV = "#008080",
  negative = "#747272"
)

# --- Loadings (variables) ---
loadings <- as.data.frame(pca$rotation[, 1:2])
loadings$variable <- rownames(loadings)

# Scale loadings to scores range for overlay
scale_factor <- max(abs(pca_scores[, c("PC1","PC2")])) / max(abs(loadings[, 1:2])) * 0.7
loadings[, 1:2] <- loadings[, 1:2] * scale_factor

loadings$magnitude <- sqrt(loadings$PC1^2 + loadings$PC2^2)


# --- Variance explained ---
var_exp <- round(summary(pca)$importance[2, 1:2] * 100, 1)
quartz()
# --- Plot ---
ggplot(pca_scores, aes(x = PC1, y = PC2, colour = target)) +
  geom_point(alpha = 0.7, size = 2) +
  # Draw all arrows, label only top ones
  geom_segment(data = loadings,
               aes(x = 0, y = 0, xend = PC1, yend = PC2),
               arrow = arrow(length = unit(0.15, "cm")),
               colour = "grey50", linewidth = 0.4,
               inherit.aes = FALSE) +
  geom_text_repel(data = loadings,
                  aes(x = PC1, y = PC2, label = variable),
                  colour = "grey20", size = 3,
                  max.overlaps = Inf,
                  box.padding = 0.4,
                  segment.colour = "grey60",
                  segment.size = 0.3,
                  inherit.aes = FALSE) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey70") +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey70") +
  theme_minimal(base_size = 13)
quartz()
print(p)

ggsave(
  filename = here("Results/cebu_pca_biplot.png"),
  plot = p,
  width = 15,
  height = 6)




scale_factor <- max(abs(pca_scores[, c("PC1","PC2")])) /
                max(abs(loadings[, c("PC1","PC2")])) * 0.85

loadings_scaled <- loadings %>%
  mutate(PC1_s = PC1 * scale_factor,
         PC2_s = PC2 * scale_factor,
         magnitude = sqrt(PC1^2 + PC2^2))   # for filtering

# ── 2. Keep only top N most influential variables (reduces clutter) ────────
top_loadings <- loadings_scaled %>%
  slice_max(magnitude, n = 40)   # adjust n as needed

# ── 3. Variance explained labels (add if you have it) ─────────────────────
# var_pct <- summary(pca)$importance[2, ] * 100


quartz()

ggplot(pca_scores, aes(x = PC1, y = PC2, colour = target)) +
  geom_point(alpha = 0.7, size = 2) +

  # All loadings as faint grey arrows (context)
  geom_segment(data = loadings_scaled,
               aes(x = 0, y = 0, xend = PC1_s, yend = PC2_s),
               arrow = arrow(length = unit(0.12, "cm"), type = "closed"),
               colour = "grey80", linewidth = 0.3,
               inherit.aes = FALSE) +

  # Top loadings with bold arrows + labels
  geom_segment(data = top_loadings,
               aes(x = 0, y = 0, xend = PC1_s, yend = PC2_s),
               arrow = arrow(length = unit(0.2, "cm"), type = "closed"),
               colour = "firebrick", linewidth = 0.6,
               inherit.aes = FALSE) +
  geom_text_repel(data = top_loadings,
                  aes(x = PC1_s, y = PC2_s, label = variable),
                  colour = "firebrick", size = 3, fontface = "bold",
                  max.overlaps = Inf,
                  box.padding   = 0.5,
                  point.padding = 0.2,
                  segment.colour = "firebrick", segment.size = 0.3,
                  inherit.aes = FALSE) +

  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey70") +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey70") +

  # Dual axes: scores (bottom/left) + loadings (top/right)
  scale_x_continuous(
    sec.axis = sec_axis(~ . / scale_factor, name = "PC1 loading")
  ) +
  scale_y_continuous(
    sec.axis = sec_axis(~ . / scale_factor, name = "PC2 loading")
  ) +

  labs(
    x = "PC1",   # replace with paste0("PC1 (", round(var_pct[1],1), "%)")
    y = "PC2",
    title = "PCA Biplot",
    colour = "Target"
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "right")
