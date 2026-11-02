
# ---- Generate infection status matrix
inf_matrix <- function(N_pathogen, pres=rep(1,N_pathogen)){
  
  # list of possible outcomes for each pathogen
  combos <- list()
  for(c in 1:N_pathogen) combos[[c]] <- c(0,1)
  
  # matrix of all possible infection status combinations
  m <- expand.grid(combos)
  colnames(m) <- letters[1:N_pathogen]
  
  # remove positives of absent pathogens
  if(sum(pres)<N_pathogen){
    for(abs in which(pres==0)) m <- m[m[,abs]==0, ]
  }
  
  return(m)
}

# ---- Build data for inputting into stan model 
prepare_multiplex_sero_data <- function(data, pathogens,present_pathogens) {
  
  # Validate inputs
  if (!all(present_pathogens %in% pathogens)) {
    stop("All present_pathogens must be in the pathogens list")
  }
  
  # Reorder pathogens (present first, then non-present)
  nonpres <- pathogens[!pathogens %in% present_pathogens]
  pathogens_ordered <- c(present_pathogens, nonpres)
  
  # Compile data for model fitting
  stan_data <- list()
  
  # Antibody titer data (log-transformed)
  stan_data$y <- cbind(data[, c(present_pathogens, nonpres)])
  stan_data$y <- stan_data$y
  
  # Pathogen presence indicators
  stan_data$pres <- c(rep(1, length(present_pathogens)), rep(0, length(nonpres)))
  
  # Basic dimensions
  stan_data$N <- nrow(stan_data$y)
  stan_data$nP <- ncol(stan_data$y)
  stan_data$nPp <- sum(stan_data$pres)
  
  # Infection matrix and related variables
  stan_data$infM <- inf_matrix(stan_data$nP, pres = stan_data$pres)
  stan_data$nC <- nrow(stan_data$infM)
  
  # Create indexes for model fitting
  npos <- rowSums(stan_data$infM)
  wpos <- wneg <- matrix(0, ncol = stan_data$nP, nrow = stan_data$nC)
  
  for (c in 1:nrow(stan_data$infM)) {
    for (p in 1:stan_data$nP) {
      if (npos[c] > 0) {
        wpos[c, 1:npos[c]] <- which(stan_data$infM[c, ] == 1)
      }
      if (npos[c] < stan_data$nP) {
        wneg[c, 1:(stan_data$nP - npos[c])] <- which(stan_data$infM[c, ] == 0)
      }
    }
  }
  
  stan_data$npos <- npos  # N pos pathogens per infection status
  stan_data$wpos <- wpos  # index pos pathogens per infection status
  stan_data$wneg <- wneg  # index neg pathogens per infection status
  
  # Return list with data and model
  return(list(
    data = stan_data,
    pathogens = pathogens_ordered
  ))
}



# ---- Generate infection status matrix (updated to have different questions explicitly)
inf_matrix_NEW <- function(question, antigens, denv = NULL, zikv = NULL, chikv = NULL) {

  row <- function(...) {
    x <- c(...)
    if (!all(x %in% antigens))
      stop("Unknown antigen(s): ", paste(setdiff(x, antigens), collapse = ", "))
    as.integer(antigens %in% x)
  }

  m <- switch(question,
    "DENV_vs_not"   = rbind(neg = row(), DENV = row(denv)),
    "Flavi_vs_not"  = rbind(neg = row(), Flavi = row(denv, zikv)),
    "CHIKV_vs_not"  = rbind(neg = row(), CHIKV = row(chikv)),
    "DENV_vs_CHIKV" = rbind(neg   = row(),
                            DENV  = row(denv),
                            CHIKV = row(chikv),
                            Both  = row(denv, chikv)),
    "DENV_serotype" = {
      s <- as.matrix(expand.grid(rep(list(c(0L, 1L)), length(denv))))
      m <- matrix(0L, nrow(s), length(antigens))
      m[, match(denv, antigens)] <- s
      rownames(m) <- c("neg", apply(s[-1, , drop = FALSE], 1, function(r)
                         paste(denv[r == 1], collapse = "+")))
      m
    },
    stop("Unknown question: ", question)
  )

  colnames(m) <- antigens
  m
}

# ---- Build data for inputting into stan model 
prepare_multiplex_sero_data_NEW<- function(data, infM) {

  pathogens <- colnames(infM)
  miss <- setdiff(pathogens, names(data))
  if (length(miss)) stop("Columns not in data: ", paste(miss, collapse = ", "))

  pres <- as.integer(colSums(infM) > 0)

  y  <- as.matrix(data[, pathogens, drop = FALSE])
  nP <- ncol(infM); nC <- nrow(infM)

  npos <- as.integer(rowSums(infM))
  wpos <- wneg <- matrix(0L, nrow = nC, ncol = nP)
  for (c in seq_len(nC)) {
    p <- which(infM[c, ] == 1); n <- which(infM[c, ] == 0)
    if (length(p)) wpos[c, seq_along(p)] <- p
    if (length(n)) wneg[c, seq_along(n)] <- n
  }

  stan_data <- list(
    y = y, N = nrow(y), nP = nP, nPp = sum(pres), pres = pres,
    infM = unname(infM), nC = nC,
    npos = npos, nneg = nP - npos, wpos = wpos, wneg = wneg,
    neg_row = which(npos == 0)
  )

  list(data = stan_data,
       pathogens = pathogens,
       class_labels = rownames(infM))
}


init <- function(data, nChains){
  n_phi <- data$nP * data$nPp - data$nPp
  
  lapply(seq_len(nChains), function(i) list(
    sero  = array(runif(data$nPp, 0.45, 0.75)),   # beta(2,3): IQR ~0.26–0.55
    mu0   = array(runif(data$nP, -0.20, 0.20)),   # normal(0, 0.3): IQR ~±0.20
    mu1   = array(runif(data$nPp, 2.80, 3.20)),   # normal(3, 0.3)
    sd0   = array(runif(data$nP,  0.50, 0.67)),   # normal(0.6, 0.1), kept > 0
    sd1   = runif(1, 0.55, 0.85),                 # normal(0.7, 0.2)
    phi   = array(runif(n_phi, 0.04, 0.17)),      # exponential(8): IQR ~0.036–0.173
    rho00 = runif(1, 0.35, 0.65)                  # beta(2,2)
  ))
}

#----- Extract prevalence estimates
extract_sd <- function(chains, data, pathogens){
  
  sd0s <- data.frame(pathogen=pathogens, par='sd0', med=NA, ciL=NA, ciU=NA)
  sd1s <- data.frame(pathogen='all', par='sd1', med=NA, ciL=NA, ciU=NA)
  
  for(p in 1:data$nP) sd0s[p, 3:5] <- quantile(chains[,paste('sd0[', paste(p,']', sep=''), sep='')], c(0.5,0.025,0.975))
  sd1s[1, 3:5] <- quantile(chains[,'sd1'], c(0.5,0.025,0.975))
  
  sds <- rbind(sd0s, sd1s)
  return(sds)
}


#----- Extract covariance matrices per iteration
extract_covM <- function(chains, data){
  
  iter <- length(chains$lp__)
  covM <- list()
  for(i in 1:iter){
    covM[[i]] <- list()
    for(c in 1:data$nC){
      x <- matrix(NA, ncol=data$nP, nrow=data$nP)
      for(p in 1:data$nP) for(p2 in 1:data$nP){
        
        y <- paste(paste(c,p,sep=','),p2,sep=',')
        name <- paste0("covM[", y, "]")
        x[p,p2] <- chains[[name]][i]        
      }
      covM[[i]][[c]] <- x
    }
  }
  return(covM)
}

#----- Extract gaussian sds
extract_sds <- function(chains, data, pathogens){
  
  sds <- data.frame(pathogen=rep(pathogens, 2), 
                    par=rep(c('sd0','sd1'), each=length(pathogens)), 
                    med=NA, ciL=NA, ciU=NA)
  
  for(p in 1:data$nP) sds[p, 3:5] <- quantile(chains[,paste('sd0[', paste(p,']', sep=''), sep='')], c(0.5,0.025,0.975))
  for(p in 1:data$nP) sds[p+data$nP, 3:5] <- quantile(chains[,paste('sd1[', paste(p,']', sep=''), sep='')], c(0.5,0.025,0.975))
  
  return(sds)
}


#----- Extract gaussian means
extract_mu <- function(chains, data, pathogens){
  
  
  # label combination positives
  pos <- rep('neg', data$nC)
  for(c in 1:data$nC){
    np <- sum(data$infM[c,])
    if(np==1) pos[c] <- pathogens[which(data$infM[c,]==1)]
    else if(np>1) pos[c] <- paste(pathogens[which(data$infM[c,]==1)], sep='&', collapse="&")
  }
  
  # all gaussian means
  mus0 <- data.frame(pg=rep(NA,length(which(data$infM==0))), pos=NA, med=NA, ciL=NA, ciU=NA)
  mus1 <- data.frame(pg=rep(NA,length(which(data$infM==1))), pos=NA, med=NA, ciL=NA, ciU=NA)
  ix0 <- 1
  ix1 <- 1
  for(c in 1:data$nC) for(p in 1:data$nP){
    if(data$infM[c,p]==0){
      mus0$pg[ix0] <- pathogens[p]
      mus0$pos[ix0] <- pos[c] 
      y <- paste(c,p, sep=',')
      mus0[ix0,3:5] <- quantile(chains[,paste(paste('mu[',y,sep=''),']',sep='')], c(0.5,0.025,0.975))
      ix0 <- ix0+1
    }else{
      mus1$pg[ix1] <- pathogens[p]
      mus1$pos[ix1] <- pos[c] 
      y <- paste(c,p, sep=',')
      mus1[ix1,3:5] <- quantile(chains[,paste(paste('mu[',y,sep=''),']',sep='')], c(0.5,0.025,0.975))
      ix1 <- ix1+1
    }
  }
  colnames(mus0)[1] <- 'antigen'
  colnames(mus1)[1] <- 'antigen'
  
  return(list(mus0=mus0, mus1=mus1))
}

extract_phi <- function(chains, data, pathogens){
  
  phi <- data.frame(pos=NA, neg=NA, med=NA, ciL=NA, ciU=NA)
  ind <- 1
  for(p in 1:data$nPp) for(p2 in 1:data$nP){
    if(!p==p2){
      phi[ind,1:2] <- c(pathogens[p], pathogens[p2])
      y <- str_replace_all(toString(c(p,p2))," ","")
      phi[ind,3:5] <- quantile(chains[,paste('CR[', paste(y, ']', sep=''), sep='')], c(0.5,0.025,0.975))
      ind <- ind+1
    }
  }
  
  rho <- data.frame(pars=c('rho00'), med=NA, ciL=NA, ciU=NA)
  rho[1,2:4] <- quantile(chains[,paste('rho00')], c(0.5,0.025,0.975))
  
  return(list(phi=phi,rho=rho))
  
}

extract_sero <- function(chains, data, pathogens){
  
  sero <- data.frame(pathogen=pathogens, med=NA, ciL=NA, ciU=NA)
  
  for(p in 1:data$nP) sero[p,2:4] <- quantile(chains[,paste('seroAll[', paste(p,']', sep=''), sep='')], c(0.5,0.025,0.975))
  sero <- sero[!sero$med==0, ]
  
  return(sero)
}


#----- Plot gaussian distribution fits
plot_fits <- function(chains, data, pathogens, show_crossreactive_for = NULL){
  
  iter <- length(chains$lp__)
  covM <- extract_covM(chains, data)
  
  # simulate multivariate gaussians per combination
  yy <- yp <- list()
  for(p in 1:data$nP) yp[[p]] <- matrix(NA, nrow=512, ncol=iter)
  ypN <- ypP <- ypNcross <- yp
  for(i in 1:iter){
    for(c in 1:data$nC){
      
      # simulate gaussion for combination c, iteration i
      g <- paste('mu[',c, sep='')
      nn <- ceiling(sum(data$N*chains[,paste('theta[',paste(c,']',sep=''), sep='')][i]))
      muu <- vector()
      for(p in 1:data$nP) muu[p] <- chains[i,paste(paste(g,p,sep=','),']',sep='')]
      yy[[c]] <- as.data.frame(rmvnorm(nn, mean=muu,sigma=covM[[i]][[c]]))
      yy[[c]]$C <- c
    }
    
    # density distributions per pathogen
    yc <- do.call('rbind',yy)
    for(p in 1:data$nP){
      
      yp[[p]][,i] <- density(yc[,p], bw=0.01, from=-2.5, to=12)$y
      
      # Check if we should show cross-reactive for this pathogen
      show_cr <- !is.null(show_crossreactive_for) && p %in% show_crossreactive_for
       
      ## Positive component
      if(data$pres[p] == 1){
        z <- which(data$infM[, p] == 1)
        pw <- numeric(length(z))
        for(s in seq_along(z)) {
          pw[s] <- chains[i, paste0("theta[", z[s], "]")]
        }
        propP <- sum(data$N * pw) / data$N
        
        ypP[[p]][, i] <- density(
          yc[yc$C %in% z, p],
          bw = 0.01, from = -2.5, to = 12
        )$y * propP
      } else {
        propP <- 0
        ypP[[p]][, i] <- rep(0, 512)
      }

      
      if(show_cr){
        # Split negatives into true negatives and cross-reactive negatives
        neg_idx <- which(data$infM[,p]==0)
        
        # Cross-reactive negatives: negative for pathogen p but positive for at least one other pathogen
        cross_reactive_idx <- neg_idx[rowSums(data$infM[neg_idx, , drop=FALSE]) > 0]
        # True negatives: negative for all pathogens
        true_neg_idx <- neg_idx[rowSums(data$infM[neg_idx, , drop=FALSE]) == 0]
        
        # Calculate proportions
        pw_cross <- vector()
        pw_true <- vector()
        if(length(cross_reactive_idx) > 0){
          for(s in 1:length(cross_reactive_idx)) pw_cross[s] <- chains[i,paste(paste('theta[',paste(cross_reactive_idx[s]),sep=''),']',sep='')]
          propNcross <- sum(data$N*pw_cross)/data$N
        } else {
          propNcross <- 0
        }
        if(length(true_neg_idx) > 0){
          for(s in 1:length(true_neg_idx)) pw_true[s] <- chains[i,paste(paste('theta[',paste(true_neg_idx[s]),sep=''),']',sep='')]
          propNtrue <- sum(data$N*pw_true)/data$N
        } else {
          propNtrue <- 0
        }
        
        # Density for cross-reactive negatives
        if(length(cross_reactive_idx) > 0){
          ypNcross[[p]][,i] <- density(yc[yc$C %in% cross_reactive_idx,p], bw=0.01, from=-2.5, to=12)$y * propNcross
        } else {
          ypNcross[[p]][,i] <- rep(0, 512)
        }
        
        # Density for true negatives
        if(length(true_neg_idx) > 0){
          ypN[[p]][,i] <- density(yc[yc$C %in% true_neg_idx,p], bw=0.01, from=-2.5, to=12)$y * propNtrue
        } else {
          ypN[[p]][,i] <- rep(0, 512)
        }
       } else {
        neg_idx <- which(data$infM[, p] == 0)
        ypN[[p]][, i] <- density(
          yc[yc$C %in% neg_idx, p],
          bw = 0.01, from = -2.5, to = 12
        )$y * (1 - propP)
        
        ypNcross[[p]][, i] <- rep(0, 512)
      }
    }
   }

  # quantiles of density distributions
  dpq <- dpqP <- dpqN <- dpqNcross <- list()
  titer <- density(yc[,1], bw=0.01, from=-2.5, to=12)$x
  for(p in 1:data$nP){
    dpq[[p]] <- as.data.frame(rowQuantiles(yp[[p]], probs=c(0.5,0.025,0.975)))
    dpqN[[p]] <- as.data.frame(rowQuantiles(ypN[[p]], probs=c(0.5,0.025,0.975)))
    dpqNcross[[p]] <- as.data.frame(rowQuantiles(ypNcross[[p]], probs=c(0.5,0.025,0.975)))
    dpqP[[p]] <- as.data.frame(rowQuantiles(ypP[[p]], probs=c(0.5,0.025,0.975)))
    dpq[[p]]$titer <- dpqN[[p]]$titer <- dpqNcross[[p]]$titer <- dpqP[[p]]$titer <- titer
    dpq[[p]]$pathogen <- dpqN[[p]]$pathogen <- dpqNcross[[p]]$pathogen <- dpqP[[p]]$pathogen <- pathogens[p]
  }
  dpq <- do.call('rbind',dpq)
  dpqN <- do.call('rbind',dpqN)
  dpqNcross <- do.call('rbind',dpqNcross)
  dpqP <- do.call('rbind',dpqP)
  colnames(dpq)[1:3] <- colnames(dpqN)[1:3] <- colnames(dpqNcross)[1:3] <- colnames(dpqP)[1:3] <- c('med','ciL','ciU')
  
  # compile data for plotting 
  dta <- as.data.frame(data$y)
  colnames(dta) <- pathogens
  dta <- tidyr::gather(dta, key='pathogen', value='t')
  
  # Set factor levels to match input order
  dta$pathogen <- factor(dta$pathogen, levels = pathogens)
  dpq$pathogen <- factor(dpq$pathogen, levels = pathogens)
  dpqN$pathogen <- factor(dpqN$pathogen, levels = pathogens)
  dpqNcross$pathogen <- factor(dpqNcross$pathogen, levels = pathogens)
  dpqP$pathogen <- factor(dpqP$pathogen, levels = pathogens)
  
  # overall fit
  fitD <- ggplot()+ geom_histogram(data=dta, aes(t,y=..density..), bins=10, fill='#dedede', col='#dedede')+
    theme_minimal() + 
    theme(text=element_text(size=20))+
    geom_line(data=dpq, aes(titer,med),col='#008b7f')+ facet_wrap(~pathogen, scales='free_y')+ xlab('titer')+
    geom_ribbon(data=dpq, aes(x=titer,y=med,ymin=ciL,ymax=ciU), fill='#02d0bf', alpha=0.4)
  
  # pos-neg fit with cross-reactive negatives
  fitDPN <- ggplot()+ geom_histogram(data=dta, aes(t,y=..density..), bins=10, fill='#dedede', col='#dedede')+
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      axis.line = element_line(color = "black", linewidth = 0.7),
      axis.ticks.x = element_line(color = "black", size = 0.5),
      axis.ticks.y = element_line(color = "black", size = 0.5),
      axis.text = element_text(size = 12),
      axis.text.x = element_text(size = 12),
      axis.title = element_text(size = 12),
      strip.text = element_text(size = 12),
      aspect.ratio = 1,
      legend.title = element_blank(),
      legend.text = element_text(size = 12)
    )+
      # Lines
      geom_line(data = dpqN,
                aes(titer, med, colour = "True negative"),
                linewidth = 1) +
      geom_line(data = dpqNcross,
                aes(titer, med, colour = "Cross-reactive negative"),
                linewidth = 1) +
      geom_line(data = dpqP,
                aes(titer, med, colour = "Positive"),
                linewidth = 1) +

      # Ribbons 
      geom_ribbon(data = dpqN,
                  aes(x = titer, y = med,
                      ymin = ciL, ymax = ciU,
                      fill = "True negative"),
                  alpha = 0.3) +
      geom_ribbon(data = dpqNcross,
                  aes(x = titer, y = med,
                      ymin = ciL, ymax = ciU,
                      fill = "Cross-reactive negative"),
                  alpha = 0.3) +
      geom_ribbon(data = dpqP,
                  aes(x = titer, y = med,
                      ymin = ciL, ymax = ciU,
                      fill = "Positive"),
                  alpha = 0.5) +
    facet_wrap(~pathogen,
        labeller = labeller(
          pathogen = c(
            "ONNV_VLP_log" = "ONNV VLP",
            "CHIKV_sE2_log" = "CHIKV E2",
            "MAYV_E2_log"  = "MAYV E2"
          )
        )
      ) +
    labs(x = "Log(titer)", y = "Density") +
    # Manual colours
    scale_colour_manual(
      values = c(
        "True negative" = "#823e29",
        "Cross-reactive negative" = "#134611",
        "Positive" = "#EF81A5"
      )
    ) +
    scale_fill_manual(
      values = c(
        "True negative" = "#823e29",
        "Cross-reactive negative" = "#134611",
        "Positive" = "#EF81A5"
      )
    )+ guides(fill = "none",  colour = "none")
  
  # return plots
  return(list(fit=fitD,fitPN=fitDPN)) 
}


plot_titres_coloured_by_clusters <- function(meta_data) {

  cluster_labels <- c("1" = "ONNV VLP -ve, CHIKV E2 -ve", "2" = "ONNV VLP +ve", "3" = "CHIKV E2 +ve")
  cluster_colours <- c("1" = "#C56244", "2" = "#648767", "3" = "#715F95")

    p <- meta_data %>%
        filter(!is.na(cluster)) %>%
        ggplot(aes(x = log(ONNV_VLP),
                  y = log(CHIKV_sE2),
                  color = factor(cluster))) +
    geom_point(alpha = 0.3, size = 2) +
    scale_color_manual(values = cluster_colours, labels = cluster_labels) +
    labs(
      x = "Log(ONNV VLP MFI)",
      y = "Log(CHIKV E2 MFI)") +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 20),
      panel.grid = element_blank(),
      axis.text.x = element_text(size = 24),
      axis.text.y = element_text(size = 24),
      axis.title = element_text(size = 20),
      legend.position = 'bottom',
      legend.title = element_blank(),
      legend.text = element_text(size = 20),
      panel.grid.major.x = element_blank(),
      aspect.ratio = 1,
    )
  return(p)
}



#----- Function to simulate multivariate Gaussian data
sim_multisero <- function(N, nP, nPp, sero, mu0, mu1, sd0, sd1, CR, corr00){
  
  #--- infection status matrix
  pres <- c(rep(1,nPp), rep(0,nP-nPp)) # index for which pathogens are present (1) or absent (0)
  infM <- inf_matrix(nP, pres=pres) # matrix of possible infection statuses
  
  #--- Gaussian weights
  w <- infM
  w[,which(pres==0)] <- 1
  for(p in 1:nPp) for(c in 1:nrow(infM)){
    if(infM[c,p]==0) w[c,p] <- 1 - sero[p]
    else w[c,p] <- sero[p]
  }
  ws <- rowProds(as.matrix(w)) # weights
  Nws <- round(N*ws) # N individuals per infection status
  
  #--- Gaussian means & standard deviations
  mus <- infM
  sigs <- infM
  for(c in 1:nrow(infM)){
    
    if(sum(infM[c,])==0){ # neg to all
      
      mus[c,] <- mu0 
      sigs[c,] <- sd0
      
    }else if(sum(infM[c,])==1){ # pos to just 1
      
      # which are pos/neg
      wp <- which(infM[c,]==1)
      wn <- which(infM[c,]==0)
      
      # pos
      mus[c,wp] <- mu0[wp] + mu1[wp]
      sigs[c,wp] <- sd1[wp]
      
      # neg with cr
      for(n in wn) mus[c,n] <- mu0[n]+ mu1[wp]*CR[wp,n]
      for(n in wn) sigs[c,n] <- sqrt(sd0[n]^2 + (sd1[wp]*CR[wp,n])^2)
      
    }else if(sum(infM[c,])>1){ # pos to > 1
      
      # which are pos/neg
      wp <- which(infM[c,]==1)
      npos <- sum(infM[c,])
      wn <- which(infM[c,]==0)
      
      # pos
      for(p in wp) mus[c,p] <- mu0[p] + mu1[p]
      for(p in wp) sigs[c,p] <- sd1[p]
      
      # neg with cr
      for(n in wn){
        mus[c,n] <- mu0[n]
        for(p in wp) mus[c,n] <- mus[c,n] + mu1[p]*CR[p,n]
        tmp <- rep(0,ncol(infM))
        tmp2 <- rep(0,((npos*(npos-1))/2))
        for(p in wp) tmp[p] <- (sd1[p]*CR[p,n])^2
        sigs[c,n] <- sqrt(sd0[n]^2 + sum(tmp))
      }
    }
  }
  
  #--- covariance matrices
  covs <- list()
  for(c in 1:nrow(infM)) covs[[c]] <- matrix(0, ncol=nP,nrow=nP)
  for(c in 1:nrow(infM)){
    
    for(p in 1:nP) for(p2 in p:nP){
      
      if(p==p2){ # Gaussian variance
        
        covs[[c]][p,p2] <- sigs[c,p]^2 
        
      }else if(sum(infM[c,])==0){ # negative to all pathogens
        
        covs[[c]][p,p2] <- covs[[c]][p2,p] <- corr00*sigs[c,p]*sigs[c,p2]
        
      }else if(infM[c,p]+infM[c,p2]==0){ # negative to both where others are pos
        
        npos <- sum(infM[c,])
        wp <- which(infM[c,]==1)
        
        if(npos==1){
          cv00 <- corr00*sd0[p]*sd0[p2] + CR[wp,p]*CR[wp,p2]*sd1[wp]^2
        }else{
          vrs <- rep(0,npos)
          for(x in 1:npos) vrs[x] = CR[wp[x],p]*CR[wp[x],p2]*sd1[wp[x]]^2
          cv00 <- corr00*sd0[p]*sd0[p2] + sum(vrs)
        }
        covs[[c]][p,p2] <- covs[[c]][p2,p] <- cv00
        
        
      }else if(infM[c,p]+infM[c,p2]==1){ # negative & positive
        
        pos <- c(p,p2)[which(infM[c,c(p,p2)]==1)]
        neg <- c(p,p2)[which(infM[c,c(p,p2)]==0)]
        
        covs[[c]][p,p2] <- covs[[c]][p2,p] <- CR[pos,neg]*sd1[pos]^2
        
      }
    }
  }
  
  # check the covariance matrixes are semi definite
  Qcovs <- list()
  for(c in 1:nrow(infM)) Qcovs[[c]] <- is.positive.semi.definite(covs[[c]])
  
  #--- simulate multivariate Gaussians
  gs <- list()
  for(c in 1:nrow(infM)){
    if(Nws[c]>0){
      gs[[c]] <- as.data.frame(rmvnorm(Nws[c], as.numeric(mus[c,]), covs[[c]]))
      gs[[c]]$status <- c
      gs[[c]][,LETTERS[seq(1,nP)]] <- NA
      for(p in 1:nP) gs[[c]][,nP+1+p] <- infM[c,p]
      gs[[c]]$stat <- paste('{',paste(infM[c,],collapse=','),'}',sep='')
      
    }
  }
  gs <- do.call('rbind',gs)
  
  
  #--- return data
  return(list(titers=gs, mus=mus, sigs=sigs, covs=covs, covs_check=Qcovs))
}