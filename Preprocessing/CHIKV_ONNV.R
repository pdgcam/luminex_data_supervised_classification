


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
model_path_censored = here('Models_Unsupervised/MultiSero_Censored.stan')
mod = cmdstan_model(model_path, pedantic = FALSE, force_recompile = TRUE)
mod_censored = cmdstan_model(model_path_censored, pedantic = FALSE, force_recompile = TRUE)

alpha_IgG_log <- readRDS("Data/by_isotype_cross_sectional/alpha_cross_sectional_logged_IgG.rds")
head(alpha_IgG_log)

table(alpha_IgG_log$target)

# Antigens to plot
alpha_antigens <- c("CHIKV_E2", "CHIKV_VLP", "ONNV_E2", "ONNV_VLP")
alpha_labels   <- c("CHIKV E2", "CHIKV VLP", "ONNV E2", "ONNV VLP")

# Check your target labels first and adjust the patterns below if needed:
# unique(alpha_IgG_log$target)
d <- alpha_IgG_log[, alpha_antigens]
d$group <- with(alpha_IgG_log, ifelse(
  grepl("CHIKV", target, ignore.case = TRUE), "CHIKV",
  ifelse(grepl("ONNV", target, ignore.case = TRUE), "ONNV",
  ifelse(grepl("negative", target, ignore.case = TRUE), "Negative", "Other"))))
d$group <- droplevels(factor(d$group, levels = c("Negative", "CHIKV", "ONNV", "Other")))

grp_cols   <- c("Negative" = "grey55", "CHIKV" = "#B2182B",
                "ONNV" = "#2166AC", "Other" = "#E08214")
grp_shapes <- c("Negative" = 1, "CHIKV" = 16, "ONNV" = 17, "Other" = 15)

# Lower panels: scatter coloured and shaped by true status
pts <- function(data, mapping, ...) {
  ggplot(data, mapping) +
    geom_point(aes(colour = group, shape = group),
               size = 2.8, stroke = 1, alpha = 0.8) +
    scale_colour_manual(values = grp_cols, name = "PCR positive") +
    scale_shape_manual(values = grp_shapes, name = "PCR positive")
}

# Diagonal: overlapping histograms by true status
hist_diag <- function(data, mapping, ...) {
  ggplot(data, mapping) +
    geom_histogram(aes(fill = group), binwidth = 0.4, alpha = 0.6,
                   position = "identity", colour = "white", linewidth = 0.2) +
    scale_fill_manual(values = grp_cols, name = "PCR positive") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05)))
}
quartz()
ggpairs(d, columns = 1:4,
        columnLabels = alpha_labels,
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



# Antigens to plot
pairs_antigens <- c("CHIKV_E2", "ONNV_VLP")
pairs_labels   <- c("CHIKV E2", "ONNV VLP")

# Check your target labels first and adjust the patterns below if needed:
# unique(alpha_IgG_log$target)
d <- alpha_IgG_log[, pairs_antigens]
d$group <- with(alpha_IgG_log, ifelse(
  grepl("CHIKV", target, ignore.case = TRUE), "CHIKV",
  ifelse(grepl("ONNV", target, ignore.case = TRUE), "ONNV",
  ifelse(grepl("negative", target, ignore.case = TRUE), "Negative", "Other"))))
d$group <- droplevels(factor(d$group, levels = c("Negative", "CHIKV", "ONNV", "Other")))

grp_cols   <- c("Negative" = "grey55", "CHIKV" = "#B2182B",
                "ONNV" = "#2166AC", "Other" = "#E08214")
grp_shapes <- c("Negative" = 1, "CHIKV" = 16, "ONNV" = 17, "Other" = 15)

# Lower panels: scatter coloured and shaped by true status
pts <- function(data, mapping, ...) {
  ggplot(data, mapping) +
    geom_point(aes(colour = group, shape = group),
               size = 2.8, stroke = 1, alpha = 0.8) +
    scale_colour_manual(values = grp_cols, name = "PCR positive") +
    scale_shape_manual(values = grp_shapes, name = "PCR positive")
}

# Diagonal: overlapping histograms by true status
hist_diag <- function(data, mapping, ...) {
  ggplot(data, mapping) +
    geom_histogram(aes(fill = group), binwidth = 0.4, alpha = 0.6,
                   position = "identity", colour = "white", linewidth = 0.2) +
    scale_fill_manual(values = grp_cols, name = "PCR positive") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05)))
}

quartz()
ggpairs(d, columns = 1:2,
        columnLabels = pairs_labels,
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


antigens <- c("CHIKV_E2", "ONNV_VLP", "MAYV_E2")
alpha_IgG_log_subset <- alpha_IgG_log[, c(antigens)]

max(alpha_IgG_log_subset$CHIKV_E2)
max(alpha_IgG_log_subset$ONNV_VLP)
max(alpha_IgG_log_subset$MAYV_E2)

table(alpha_IgG_log$CHIKV_E2 >= 3.999)
table(censored = alpha_IgG_log$CHIKV_E2 >= 3.999, PCR = alpha_IgG_log$target)

# Fit multi-sero
preprocessed_data_alpha_IgG_log <- prepare_multiplex_sero_data(
  data = alpha_IgG_log_subset, 
  pathogens = antigens,
  present_pathogens = "CHIKV_E2"
)

# Fit multi-sero
preprocessed_data_alpha_IgG_log_ONNV <- prepare_multiplex_sero_data(
  data = alpha_IgG_log_subset, 
  pathogens = antigens,
  present_pathogens = "ONNV_VLP"
)

add_censoring <- function(dat, j, tol = 1e-6) {
  ymat <- as.matrix(dat$y)                  # individuals in rows, antigens in columns
  stopifnot(nrow(ymat) == dat$N)
  U <- max(ymat[, j])
  dat$j_cens <- j
  dat$U      <- U
  dat$cens   <- as.integer(ymat[, j] >= U - tol)
  dat
}

dat_chikv <- add_censoring(preprocessed_data_alpha_IgG_log$data,
                           j = which(preprocessed_data_alpha_IgG_log$pathogens == "CHIKV_E2"))
dat_onnv  <- add_censoring(preprocessed_data_alpha_IgG_log_ONNV$data,
                           j = which(preprocessed_data_alpha_IgG_log_ONNV$pathogens == "CHIKV_E2"))
sum(dat_chikv$cens); sum(dat_onnv$cens)   # should be 11 each
dat_chikv$U; dat_onnv$U                   # should be 4


saveRDS(preprocessed_data_alpha_IgG_log, here('Data/UnsupervisedModel/preprocessed_data_alpha_IgG_log.rds'))
saveRDS(preprocessed_data_alpha_IgG_log_ONNV, here('Data/UnsupervisedModel/preprocessed_data_alpha_IgG_log_ONNV.rds'))

#--- Fit full model 
ini <- init(preprocessed_data_alpha_IgG_log$data, nChains = 3)


 
fit_censored <- mod_censored$sample(
  data = dat_chikv, 
  chains = 3, 
  iter_sampling = 3000, 
  refresh = 100, 
  iter_warmup = 1000, 
  parallel_chains = 3,
  init = ini,
  save_cmdstan_config=TRUE
)


fit_ONNV_censored <- mod_censored$sample(
  data = dat_onnv, 
  chains = 3, 
  iter_sampling = 3000, 
  refresh = 100, 
  iter_warmup = 1000, 
  parallel_chains = 3,
  init = ini,
  save_cmdstan_config=TRUE
)


#saveRDS(fit, '/Users/ap2488/Desktop/Cameroon_Analysis_2025/16thDEC_CHIK+ONNV_MultiSeroFit.rds')
fit$save_object(here('Results/MultiSero_fits_CHIKV.rds'))
fit_ONNV$save_object(here('Results/MultiSero_fits_ONNV.rds'))

#saveRDS(fit, '/Users/ap2488/Desktop/Cameroon_Analysis_2025/16thDEC_CHIK+ONNV_MultiSeroFit.rds')
fit_censored$save_object(here('Results/Censored_MultiSero_fits_CHIKV.rds'))
fit_ONNV_censored$save_object(here('Results/Censored_MultiSero_fits_ONNV.rds'))

fit_censored <- readRDS(here('Results/Censored_MultiSero_fits_CHIKV.rds'))
fit_ONNV_censored <- readRDS(here('Results/Censored_MultiSero_fits_ONNV.rds'))

# extract chains
chains <- fit_censored$draws(format='df')
chains_df <- as.data.frame(chains)

# extract chains
chains_ONNV <- fit_ONNV_censored$draws(format='df')
chains_df_ONNV <- as.data.frame(chains_ONNV)


quartz()
# Plot trace plots with all chains clearly visible
color_scheme_set("mix-blue-red")
p1 <- mcmc_trace(fit_censored$draws(c("seroAll", "lp__")))
p2 <- mcmc_trace(fit_censored$draws(c("mu0", "mu1")))
p3 <- mcmc_trace(fit_censored$draws(c('sd0','sd1')))
p4 <- mcmc_trace(fit_censored$draws(c('phi','rho00')))
print(p1 + p2 + p3 + p4)


quartz()
# Plot trace plots with all chains clearly visible
color_scheme_set("mix-blue-red")
p1 <- mcmc_trace(fit_ONNV_censored$draws(c("seroAll", "lp__")))
p2 <- mcmc_trace(fit_ONNV_censored$draws(c("mu0", "mu1")))
p3 <- mcmc_trace(fit_ONNV_censored$draws(c('sd0','sd1')))
p4 <- mcmc_trace(fit_ONNV_censored$draws(c('phi','rho00')))
print(p1 + p2 + p3 + p4)


distfits <- plot_fits(chains_df, preprocessed_data_alpha_IgG_log$data, pathogens=preprocessed_data_alpha_IgG_log$pathogens, 
                    show_crossreactive_for = seq_along(preprocessed_data_alpha_IgG_log$pathogens))

quartz()
distfits$fit
distfits$fitPN



distfits_onnv <- plot_fits(chains_df_ONNV, preprocessed_data_alpha_IgG_log_ONNV$data, pathogens=preprocessed_data_alpha_IgG_log_ONNV$pathogens, 
                    show_crossreactive_for = seq_along(preprocessed_data_alpha_IgG_log_ONNV$pathogens))

quartz()
distfits_onnv$fit
distfits_onnv$fitPN



loo_chikv <- fit_censored$loo(variables = "log_lik", cores = 3)
loo_onnv  <- fit_ONNV_censored$loo(variables = "log_lik", cores = 3)
 
print(loo_chikv)   # check Pareto k diagnostics: k > 0.7 means LOO is unreliable
print(loo_onnv)


ll_df <- data.frame(
  loglik = c(as.vector(fit_censored$draws("sumloglik", format = "matrix")),
             as.vector(fit_ONNV_censored$draws("sumloglik", format = "matrix"))),
  model  = rep(c("CHIKV", "ONNV"),
               times = c(ndraws(fit_censored$draws("sumloglik")),
                         ndraws(fit_ONNV_censored$draws("sumloglik"))))
)

ll_sum <- ll_df %>%
  group_by(model) %>%
  summarise(median = median(loglik),
            lo     = quantile(loglik, 0.025),
            hi     = quantile(loglik, 0.975))
print(ll_sum)

quartz()
ggplot(ll_sum, aes(median, model, colour = model)) +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.15, linewidth = 0.8) +
  geom_point(size = 4) +
  scale_colour_manual(values = c(CHIKV = "#B2182B", ONNV = "#2166AC"), guide = "none") +
  labs(x = "Total log-likelihood (median, 95% CrI)", y = NULL) +
  theme_minimal(base_size = 16) +
  theme(panel.grid.minor = element_blank())


loo_chikv <- fit_censored$loo(variables = "log_lik", cores = 3)
loo_onnv  <- fit_ONNV_censored$loo(variables = "log_lik", cores = 3)

loo_compare(list(CHIKV = loo_chikv, ONNV = loo_onnv))   # elpd_diff, se_diff
loo_model_weights(list(CHIKV = loo_chikv, ONNV = loo_onnv), method = "stacking")

