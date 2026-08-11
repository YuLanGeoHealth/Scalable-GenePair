#analysis results
library("dplyr")
library(ggdist)
library(tidyr)
library(parallelMCMCcombine)
library(coda)
library(metafor)
library(tibble)

sim_table <- list()
M = 10 #subsets
#pb <- txtProgressBar(min = 0, max = 99, style = 3)
for(sim in (1:100)){
  ####meta analysis
  d <- 3 # dimension of the parameter space 
  sampT <- 10000 # number of subset posterior samples
  theta <- array(NA,c(d,sampT,M)) # M total number of subsets
  #d=2 beta[2] age differences and beta[3] sp distance
  
  df_RR = data.frame(matrix( 
    vector(), 0, 3, dimnames=list(c(), c("Group","yi","vi"))), 
    stringsAsFactors=F)
  summary_list <- list()  # create empty list before loop
  name = c("AgeDifference","SpatialDistance","AgeSum")
  for(v_id in (1:3)){
    key_name = name[v_id]
    for(i in (1:M)){
      group = i-1
      data_source <- paste0("./data/sim",sim,"_subset",M,"_",group,".RData")
      data_res <- paste0("./results/result_m", sprintf("%03d", sim), "_n", group, ".rds")
      results = readRDS(data_res)
      load(data_source)
      
      #reverse SD for interpret
      theta[1,,i] <- (readRDS(data_res)$beta[2,seq(10001, 60000, by = 5)])
      theta[2,,i] <- (readRDS(data_res)$beta[3,seq(10001, 60000, by = 5)])
      theta[3,,i] <- (readRDS(data_res)$gamma[1,seq(10001, 60000, by = 5)])
      
      eff<- mean(theta[v_id,,i])
      SE <- sd(theta[v_id,,i])
      # Calculate variance of the log-relative risk
      variance <- SE^2
      df_RR[i,]<-c(paste0("group_",group),eff,variance)
      
    }
    
    res <- rma(as.numeric(yi), as.numeric(vi), data=df_RR, method="FE",
               slab=df_RR$Group)
    
    # Extract results
    summary_table <- tibble(
      key = key_name,  # or extract from model context
      value = exp(res$b),   # exponentiate to get OR
      .lower = exp(res$ci.lb),
      .upper = exp(res$ci.ub),
      .width = 0.95,
      .point = "estimate",
      .interval = "ci"
    )
    summary_list[[v_id]] <- summary_table
    
  }
  
  sum_meta <- bind_rows(summary_list)
  sum_meta$source = "Meta"
  ####meta analysis end
  
  #####################semiparamDPE#####################
  ## estimate (mean) standard deviations for each parameter across the subsets
  norm.var.est <- rep(0,d)
  for(i in 1:d){
    for(s in 1:M){
      norm.var.est[i] <- norm.var.est[i] + var(theta[i,,s])
    }
  }
  norm.sd.est <- sqrt(norm.var.est/M)
  
  ## Compute the diagonal of the optimal bandwidth
  ## matrix according to Silverman's rule
  h_opt1 = (4/(d+2))^(1/(4+d)) * (sampT^(-1/(4+d))) * norm.sd.est
  ## Combine samples. The bandwidth matrix is fixed:
  full.theta1 <- exp(semiparamDPE( subchain = theta, bandw = h_opt1 * 2, anneal = FALSE))
  
  
  ## combine samples:
  cov_fin <- t(rbind(full.theta1[1,], full.theta1[2,],full.theta1[3,]))
  colnames(cov_fin) <- c("AgeDifference","SpatialDistance","AgeSum")
  cov_fin <- data.frame(cov_fin) %>% gather()
  
  # check convergency
  #plot(full.theta1[1,], type = "l",  ylab = "ageD") 
  #plot(full.theta1[2,], type = "l",  ylab = "spD") 
  #plot(full.theta1[3,], type = "l",  ylab = "ageSum") 
  #analyze combined results
  sum4 <- cov_fin %>% group_by(key) %>%ggdist::median_hdci(.width = c(0.95))%>%
    mutate(source = "DPE1")
  
  ## Compute the diagonal of the optimal bandwidth
  ## matrix for the method that uses annealing
  h_opt2 = (4/(d+2))^(1/(4+d)) * norm.sd.est
  ## Combine samples. The bandwidth matrix will be annealed:
  full.theta2 <- exp(semiparamDPE(subchain = theta, bandw = h_opt2 * 2, anneal = TRUE))
  cov_fin <- t(rbind(full.theta2[1,], full.theta2[2,],full.theta2[3,]))
  colnames(cov_fin) <- c("AgeDifference","SpatialDistance","AgeSum")
  cov_fin <- data.frame(cov_fin) %>% gather()
  # check convergency
  #plot(full.theta2[1,], type = "l",  ylab = "ageD") 
  #plot(full.theta2[2,], type = "l",  ylab = "spD") 
  #plot(full.theta2[3,], type = "l",  ylab = "ageSum") 
  #analyze combined results
  sum5 <- cov_fin %>% group_by(key) %>%ggdist::median_hdci(.width = c(0.95))%>%
    mutate(source = "DPE2")
  #########################################################
  
  #####################sampleAvg#####################
  ## combine samples:
  full.theta <- exp(sampleAvg(subchain=theta, shuff=FALSE))
  cov_fin <- t(rbind(full.theta[1,], full.theta[2,],full.theta[3,]))
  colnames(cov_fin) <- c("AgeDifference","SpatialDistance","AgeSum")
  cov_fin <- data.frame(cov_fin) %>% gather()
  
  #analyze combined results
  sum1 <- cov_fin %>% group_by(key) %>%ggdist::median_hdci(.width = c(0.95))%>%
    mutate(source = "sampleAvg")
  #########################################################
  
  #####################consensusMCindep#####################
  ## combine samples:
  full.theta <- exp(consensusMCindep(subchain=theta, shuff=FALSE))
  cov_fin <- t(rbind(full.theta[1,], full.theta[2,],full.theta[3,]))
  colnames(cov_fin) <- c("AgeDifference","SpatialDistance","AgeSum")
  cov_fin <- data.frame(cov_fin) %>% gather()
  
  #analyze combined results
  sum2 <- cov_fin %>% group_by(key) %>%ggdist::median_hdci(.width = c(0.95))%>%
    mutate(source = "MCindep")
  #########################################################
  
  #####################consensusMCcov#####################
  ## combine samples:
  full.theta <- exp(consensusMCcov(subchain=theta, shuff=FALSE))
  cov_fin <- t(rbind(full.theta[1,], full.theta[2,],full.theta[3,]))
  colnames(cov_fin) <- c("AgeDifference","SpatialDistance","AgeSum")
  cov_fin <- data.frame(cov_fin) %>% gather()
  
  #analyze combined results
  sum3 <- cov_fin %>% group_by(key) %>%ggdist::median_hdci(.width = c(0.95))%>%
    mutate(source = "MCcov")
  #########################################################
  
  all_sum <- bind_rows(sum_meta,sum1, sum2, sum3, sum4, sum5)
  
  truth <- tibble::tibble(
    key = c("AgeDifference", "AgeSum", "SpatialDistance"),
    true_value = c(0.752, 0.754, 0.516) #the reference values from the benchmark dataset
  )
  
  comparison <- all_sum %>%
  left_join(truth, by = "key") %>%
  mutate(
    bias = value - true_value,
    coverage = true_value >= .lower & true_value <= .upper,
    ci_width = .upper - .lower
  ) %>%
  select(key, source, value, true_value, bias, coverage, ci_width)
  
  # Bias per key
  by_key <- comparison %>%
  group_by(source, key) %>%
  summarise(
    bias = mean(bias),
    coverage = mean(coverage),         # Proportion of TRUEs = coverage rate
    ci_width = mean(ci_width),
    .groups = "drop"
  )
  
  by_key <- by_key %>%
    mutate(source = factor(source, levels = c(setdiff(unique(source), "Meta"), "Meta")))
  
  # Overall mean bias per source
  mean_bias <- by_key %>%
  group_by(source) %>%
  summarise(
    bias = mean(bias),
    coverage = mean(coverage),
    ci_width = mean(ci_width),
    .groups = "drop") %>%
  mutate(key = "mean_bias")  # fake key
  
  # Combine
  plot_data <- bind_rows(by_key, mean_bias)
  # Reorder 'key' so that 'mean_bias' is last
  plot_data <- plot_data %>%
    mutate(key = factor(key, levels = c(setdiff(unique(key), "mean_bias"), "mean_bias")))
  
    
  plot_data$simulation = sim
  sim_table[[sim+1]] <- plot_data
 # setTxtProgressBar(pb, i)
}

final_bias <- bind_rows(sim_table)
write.csv2(final_bias,"final_bias_update.csv")
#close(pb)
