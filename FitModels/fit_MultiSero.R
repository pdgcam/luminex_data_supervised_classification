
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

denv_IgG_log_NS1_subset <- as.data.frame(flavi_IgG_log_NS1_subset[, denv_NS1_antigens])


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

# get marginal prob  of being seropositive for each pathogen
neg <- which(rowSums(preprocessed_data_flavi_IgG_log$data$infM) == 0)
N <- preprocessed_data_flavi_IgG_log$data$N
p_denv <- 1 - unclass(fit$draws(paste0("post_prob[", seq_len(N), ",", neg, "]"),
                                format = "draws_matrix"))

# use PCR confirmed infections as truth to calculate AUC for DENV vs not DENV                 
truth_denv_vs_not <- ifelse(flavi_IgG_log$target %in% c("DENV1","DENV2","DENV3","DENV4"), 1L,
         ifelse(flavi_IgG_log$target %in% c("negative","ZIKV"), 0L, NA_integer_))

auc <- function(p, y) {
  r <- rank(p)                                  # average ranks handle ties
  n1 <- sum(y == 1); n0 <- sum(y == 0)
  (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}
auc_draws <- apply(p_denv, 1, auc, y = truth_denv_vs_not)
quantile(auc_draws, c(0.025, 0.5, 0.975))


denv_labels <- c("DENV1", "DENV2", "DENV3", "DENV4")

# plot posterior predictive distribution for each pathogen
d <- flavi_IgG_log[, denv_NS1_antigens]
d$p_denv <- colMeans(p_denv)
d$group  <- factor(ifelse(grepl("^DENV", flavi_IgG_log$target), "DENV", "No DENV"),
                   levels = c("No DENV", "DENV"))

pts <- function(data, mapping, ...) {
  ggplot(data, mapping) +
    geom_point(aes(colour = p_denv, shape = group), size = 2.8, stroke = 1) +
    scale_colour_gradientn(colours = c("#2166AC", "#B2182B"),
                           limits = c(0, 1), name = "P(DENV+)") +
    scale_shape_manual(values = c("No DENV" = 1, "DENV" = 16), name = "True status")
}


hist_diag <- function(data, mapping, ...) {
  ggplot(data, mapping) +
    geom_histogram(aes(fill = group), binwidth = 0.4, alpha = 0.6,
                   position = "identity", colour = "white", linewidth = 0.2) +
    scale_fill_manual(values = c("No DENV" = "#2166AC", "DENV" = "#B2182B")) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05)))
}

ggpairs(d, columns = 1:4,
        columnLabels = sub("_NS1", "", denv_labels),
        lower = list(continuous = pts),
        diag  = list(continuous = hist_diag),
        upper = list(continuous = wrap("cor", method = "spearman")),
        legend = c(2, 1)) +
  theme_minimal(base_size = 16) +
  theme(legend.position  = "right",
        panel.grid.minor = element_blank(),
        panel.spacing    = unit(8, "pt"),
        axis.text        = element_text(size = 16),
        strip.text       = element_text(size = 16),
        legend.title     = element_text(size = 16),
        legend.text      = element_text(size = 16))


