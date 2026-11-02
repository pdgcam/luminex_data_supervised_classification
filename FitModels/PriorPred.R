
# read preprocessed data 
flavi_IgG_log <- readRDS("Data/by_isotype/flavi_cross_sectional_logged_IgG.rds")

colnames(flavi_IgG_log)

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


# rank samples within each patient: smallest days_since_infection = first
sample_num <- ave(flavi_IgG_log$days_since_infection, flavi_IgG_log$id_patient,
                 FUN = function(x) rank(x, ties.method = "first", na.last = "keep"))

flavi_IgG_log$visit <- c("first", "second", "third", "fourth")[sample_num]

flavi_IgG_log_NS1_subset <- flavi_IgG_log[, colnames(flavi_IgG_log) %in% flavi_NS1_antigens]
rownames(flavi_IgG_log_NS1_subset) <- paste(flavi_IgG_log$id_patient, flavi_IgG_log$visit, sep = "_")


present_antigens <- c("DENV1_NS1", "DENV2_NS1", "DENV3_NS1", "DENV4_NS1", "ZIKV_NS1")

# prior predictive:
nP  <- 8
nPp <- 5 
N   <- 1000
n_draws <- 100 

# store simulated titers
prior_pred <- list()

for(i in 1:n_draws) {
  
  # --- sample from priors ---
  sero  <- rbeta(nPp, 3, 2)
  mu0 <- rnorm(nP, 0.5, 0.3)
  mu1 <- rnorm(nPp, 3.5, 0.2)      
  sd0 <- rnorm(nP, 0.6, 0.1)     
  sd1   <- rnorm(1, 0.7, 0.2)     
  rho00 <- rbeta(1, 2, 2)
  
  # cross-reactivity matrix (nPp x nP)
  n_phi <- (nP * nPp) - nPp
  phi <- rexp(n_phi, rate = 8)
  CR    <- matrix(0, nrow = nPp, ncol = nP)
  idx   <- 1
  for(p in 1:nPp) for(p2 in 1:nP) {
    if(p != p2) { CR[p, p2] <- phi[idx]; idx <- idx + 1 }
  }
  
  sd0_vec <- sd0      
  sd1_vec <- rep(sd1, nP)  # expand scalar to vector
  
  tryCatch({
    sim <- sim_multisero(N = N, nP = nP, nPp = nPp,
                         sero = sero, mu0 = mu0, mu1 = mu1,
                         sd0 = sd0_vec, sd1 = sd1_vec,
                         CR = CR, corr00 = rho00)
    
    prior_pred[[i]] <- sim$titers[, 1:nP]
  }, error = function(e) {message("Error at draw ", i, ": ", e$message)})
  
  
}

# --- plot prior predictive distributions ---
prior_pred <- Filter(Negate(is.null), prior_pred)

antigens <- c(present_antigens, setdiff(names(flavi_IgG_log_NS1_subset), present_antigens))
flavi_IgG_log_NS1_subset <- flavi_IgG_log_NS1_subset[, antigens]

quartz()
# --- 2) all panels in one figure ---
op <- par(mfrow = n2mfrow(nP + 1),        # grid sized automatically; +1 panel for the legend
          mar = c(4, 4, 2.5, 1))

for (p in seq_len(nP)) {
  antigen <- antigens[p]

  all_titers <- unlist(lapply(prior_pred, function(x) x[, p]))
  obs <- flavi_IgG_log_NS1_subset[[antigen]]
  obs <- obs[!is.na(obs)]

  obs_dens <- density(obs)
  h <- hist(all_titers, breaks = 150, plot = FALSE)

  xlim <- range(c(all_titers, obs), na.rm = TRUE)
  ylim <- c(0, max(c(h$density, obs_dens$y)) * 1.05)

  hist(all_titers, breaks = 150, freq = FALSE,
       xlim = xlim, ylim = ylim,
       xlab = "log titer",
       main = if (p <= nPp) paste0(antigen, " (present)") else antigen,
       col = adjustcolor("steelblue", 0.4), border = "white")
  lines(obs_dens, col = "tomato", lwd = 2)
}

# legend in the empty last panel
plot.new()
legend("center",
       legend = c("prior predictive", "observed"),
       pch    = c(15, NA), pt.cex = 2,
       lty    = c(NA, 1),  lwd = 2,
       col    = c(adjustcolor("steelblue", 0.4), "tomato"),
       bty    = "n")
