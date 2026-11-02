
library(posterior)
library(tidyverse)
library(cmdstanr)
library(ggraph)
library(lhs)
library(mvtnorm)
library(matrixcalc)
library(emdbook)
library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(matrixStats)
library(stringr)
library(bayesplot)
library(data.table)
library(loo)
library(dplyr)
library(tidyr)
library(here)
library(patchwork)
library(mixR)
library(writexl)
library(tidybayes)
library(GGally)
library(ggplot2)


# Clear conda's settings FIRST, before loading cmdstanr
Sys.unsetenv(c("CMDSTAN", "CXX", "CC", "CXXFLAGS", "CPPFLAGS", "LDFLAGS", "CFLAGS"))

library(cmdstanr)

source(here('FitModels/MultiSeroFunctions.R'))

# Recompile model
model_path = here('Models_Unsupervised/MultiSero.stan')
mod = cmdstan_model(model_path, pedantic = FALSE, force_recompile = TRUE)


# read cross-sectional IgG data 
flavi_IgG_log <- readRDS("Data/by_isotype/flavi_cross_sectional_logged_IgG.rds")
nrow(flavi_IgG_log)


# rank samples within each patient: smallest days_since_infection = first
sample_num <- ave(flavi_IgG_log$days_since_infection, flavi_IgG_log$id_patient,
                 FUN = function(x) rank(x, ties.method = "first", na.last = "keep"))

flavi_IgG_log$visit <- c("first", "second", "third", "fourth")[sample_num]


# --- all flavis --- 
flavi_antigens <- c("DENV1_DIII","DENV1_NS1","DENV1_VLP","SHERPADES_DENV1_DIII",
"DENV2_DIII","DENV2_NS1","DENV2_VLP","SHERPADES_DENV2_DIII",
"DENV3_DIII","DENV3_NS1","DENV3_VLP", "SHERPADES_DENV3_DIII",
"DENV4_DIII","DENV4_NS1","DENV4_VLP", "SHERPADES_DENV4_DIII",
"JEV_E", "JEV_NS1", "SHERPADES_JEV_DIII",
"YFV_E", "YFV_NS1", "SHERPADES_YFV_DIII",
"WNV_DIII","WNV_NS1","SHERPADES_WNV_DIII",
"ZIKV_NS1","ZIKV_VLP","ZIKV_DIII","ZIKVSU_NS1","SHERPADES_ZIKV_DIII")

alpha_antigens <- c("CHIKV_E2", "CHIKV_NSP123", "CHIKV_VLP", "SHERPADES_CHIKV_E2", 
                    "MAYV_E2" , "SHERPADES_MAYV_E2",
                    "ONNV_E2", "ONNV_VLP",
                    "RR" , "SHERPADES_RR")


flavi_NS1_antigens <- c("DENV1_NS1","DENV2_NS1","DENV3_NS1","DENV4_NS1","ZIKV_NS1",
                        "JEV_NS1","YFV_NS1","WNV_NS1")

denv_NS1_antigens <- c("DENV1_NS1","DENV2_NS1","DENV3_NS1","DENV4_NS1")   
zikv_NS1_antigens <- c("ZIKV_NS1")


# --- 
flavi_IgG_log_NS1_subset <- flavi_IgG_log[, colnames(flavi_IgG_log) %in% flavi_NS1_antigens]
rownames(flavi_IgG_log_NS1_subset) <- paste(flavi_IgG_log$id_patient, flavi_IgG_log$visit, sep = "_")


flavi_vs_not_present_antigens <- c("DENV1_NS1", "DENV2_NS1", "DENV3_NS1", "DENV4_NS1", "ZIKV_NS1")
denv_vs_not_present_antigens <- c("DENV1_NS1", "DENV2_NS1", "DENV3_NS1", "DENV4_NS1")
denv_serotype_present_antigens <- c("DENV1_NS1", "DENV2_NS1", "DENV3_NS1", "DENV4_NS1")

denv_NS1_antigens <- c("DENV1_NS1","DENV2_NS1","DENV3_NS1","DENV4_NS1")   
zikv_NS1_antigens <- c("ZIKV_NS1")

denv_IgG_log_NS1_subset <- as.data.frame(flavi_IgG_log_NS1_subset[, denv_antigens])


# -- Plot data ---
diag_dens <- function(data, mapping, ...) {
  ggplot(data, mapping) +
    geom_histogram(aes(y = after_stat(density)), bins = 10,
                   fill = "grey75", colour = "white") +
    scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05)))
}

ggpairs(
  denv_IgG_log_NS1_subset[, denv_antigens],
  lower = list(continuous = wrap("points", alpha = 0.8, size = 1.8)),
  diag  = list(continuous = diag_dens),
  upper = list(continuous = wrap("cor", method = "spearman"))
)


# Fit multi-sero
preprocessed_data_flavi_IgG_log <- prepare_multiplex_sero_data(
  data = flavi_IgG_log_NS1_subset, 
  pathogens = flavi_NS1_antigens,
  present_pathogens = denv_NS1_antigens
)


saveRDS(preprocessed_data_flavi_IgG_log, here('Data/UnsupervisedModel/preprocessed_data_flavi_IgG_log.rds'))

#--- Fit full model 
ini <- init(preprocessed_data_flavi_IgG_log$data, nChains = 3)

 
fit <- mod$sample(
  data = preprocessed_data_flavi_IgG_log$data, 
  chains = 3, 
  iter_sampling = 3000, 
  refresh = 100, 
  iter_warmup = 1000, 
  parallel_chains = 3,
  init = ini,
  save_cmdstan_config=TRUE
)

#saveRDS(fit, '/Users/ap2488/Desktop/Cameroon_Analysis_2025/16thDEC_CHIK+ONNV_MultiSeroFit.rds')
fit$save_object(here('Results/MultiSero_fits.rds'))

# Read RDS files
fit <- readRDS(here('Results/MultiSero_fits.rds'))
fit

# extract chains
chains <- fit$draws(format='df')
chains_df <- as.data.frame(chains)



quartz()
# Plot trace plots with all chains clearly visible
color_scheme_set("mix-blue-red")
p1 <- mcmc_trace(fit$draws(c("seroAll", "lp__")))
p2 <- mcmc_trace(fit$draws(c("mu0", "mu1")))
p3 <- mcmc_trace(fit$draws(c('sd0','sd1')))
p4 <- mcmc_trace(fit$draws(c('phi','rho00')))
print(p1 + p2 + p3 + p4)
print(p1)


n_keep <- 500
keep_rows <- round(seq(1, nrow(chains_df), length.out = n_keep))
chains_thin <- chains_df[keep_rows, ]

chains_2_3 <- chains_df %>%
  filter(.chain %in% c(2, 3))

distfits <- plot_fits(chains_df, preprocessed_data_flavi_IgG_log$data, pathogens=preprocessed_data_flavi_IgG_log$pathogens, 
                    show_crossreactive_for = seq_along(preprocessed_data_flavi_IgG_log$pathogens))

quartz()
distfits$fit
distfits$fitPN


fit$diagnostic_summary()                      # divergences, treedepth
fit$summary(c("mu", "sigma", "theta"))        # Rhat, ESS (use your parameter names)
mcmc_trace(fit$draws(c("mu", "sigma")))       # chains swapping places = label switching

# Save posterior estimates from full model
# extract posterior estimates
phi <- extract_phi(chains_df, preprocessed_data_flavi_IgG_log$data, pathogens=preprocessed_data_flavi_IgG_log$pathogens)
mu <- extract_mu(chains_df, preprocessed_data_flavi_IgG_log$data, pathogens=preprocessed_data_flavi_IgG_log$pathogens)
sero <- extract_sero(chains_df, preprocessed_data_flavi_IgG_log$data, pathogens=preprocessed_data_flavi_IgG_log$pathogens)
sds <- extract_sd(chains_df, preprocessed_data_flavi_IgG_log$data, pathogens=preprocessed_data_flavi_IgG_log$pathogens)

print(sero)
print(mu)


# --- Cluster assignment based on max probability - For INLA analysis 
N  <- preprocessed_data_flavi_IgG_log$data$N
nC <- preprocessed_data_flavi_IgG_log$data$nC
draws_post <- as_draws_df(fit$draws("post_prob"))
prob_matrix <- matrix(NA_real_, nrow = N, ncol = nC)
for (n in 1:N) {
  for (c in 1:nC) {
    prob_matrix[n, c] <- mean(draws_post[[sprintf("post_prob[%d,%d]", n, c)]])
  }
}

cluster_label <- max.col(prob_matrix, ties.method = "first") #returns leftmost if two cols have equal (and max) prob
cluster_prob <- apply(prob_matrix, 1, max)
table(cluster_label)

length(cluster_label)

cluster_df <- flavi_IgG_log|>
  dplyr::mutate(
    cluster = cluster_label,
    cluster_prob = cluster_prob
  )
View(cluster_df)


