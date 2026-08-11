
# Simulate processing of divide-and-conquer approach
#simulation of scalable Genepair#

args <- commandArgs(trailingOnly = TRUE)
jobid <- as.integer(args[1])

set.seed(jobid+1014)  # set a unique seed for each run

start_time <- Sys.time()
cat(paste("Running job", jobid, "with seed", jobid, "\n"))

library(mnormt)  #Multivariate normal distribution
#library(ggplot2)
library(dplyr)
library(GenePair)
library(Matrix) #z matrix
M = 10 # subset of one simus
####################################################
#Setting the global data values:

TB_cluster <- read.csv(".../xx.csv")
TB_cluster <-TB_cluster[!is.na(TB_cluster$Response_X_M),]
full_data <- TB_cluster

max <- full_data %>% distinct(Response_X_M,Response_Y_M, .keep_all = TRUE) %>%
  select(Response_X_M, Response_Y_M)

max_spatial_dists<-as.matrix(dist(max,
                                  diag = TRUE,
                                  upper = TRUE))
max_dis <- max(max_spatial_dists) 

####################################################
 
results <- readRDS("xxx.rds") 

sample_size = 60000
beta_1 <- mean(results$beta[1,seq(10001, sample_size,by = 5)]) 
beta_2 <- mean(results$beta[2,seq(10001, sample_size,by = 5)])
beta_3 <- mean(results$beta[3,seq(10001, sample_size,by = 5)])
beta_true <- c(beta_1,beta_2,beta_3) 

gamma_true <- mean(results$gamma[1,seq(10001, sample_size,by = 5)])

tau2_true <- mean(results$tau2[seq(10001, sample_size,by = 5)])
sigma2_zeta_true <- mean(results$sigma2_zeta[seq(10001, sample_size,by = 5)])
phi_true <- mean(results$phi[seq(10001, sample_size,by = 5)])

n_truth = 700

simu_indices <- sample(1:nrow(full_data), n_truth)  # Randomly select 700 subjects
simu__data <- full_data[simu_indices,c("age","Response_X_M","Response_Y_M") ]
simu__data <- simu__data[, c("age","Response_X_M","Response_Y_M")]
simu__data <- simu__data %>% 
  mutate(ID = row_number())

df_random <- simu__data
df_random$group = 0
subeset_size = n_truth%/%M

for(i in 1:(M-1)){
  tmp<-df_random[df_random$group==0,]
  group <- sample(tmp$ID, subeset_size, replace = FALSE)
  for(j in 1:nrow(df_random)){
    if(df_random$ID[j] %in% group){df_random$group[j]=i}
  }
}

table(df_random$group)

TB_cluster <- df_random
table(TB_cluster$group)

#standardize covariates before division######################

#Setting the global data values:
n<-nrow(TB_cluster)  #Number of individuals
n_star<-0  #Number of paired responses
for(j in 1:(n - 1)){ 
  for(k in (j + 1):n){
    n_star<-n_star + 1
  }
}

#Unique locations (one for each individual)
locs <- TB_cluster %>% distinct(Response_X_M,Response_Y_M, .keep_all = TRUE) %>%
  select(Response_X_M, Response_Y_M)

spatial_dists<-as.matrix(dist(locs,
                              diag = TRUE,
                              upper = TRUE))

spatial_dists<-spatial_dists/max_dis

#Unique locations ID
locs <- locs %>% 
  mutate(locID = row_number())

#add unique locID to all cases
TB_cluster <- inner_join(TB_cluster, locs, by = c("Response_X_M", "Response_Y_M"))
# Select only distinct columns
TB_cluster <- TB_cluster %>%
  select(ID, age, Response_X_M, Response_Y_M, locID,group)

v<-matrix(0,
          nrow=n,
          ncol=nrow(locs))
for(i in 1:n){
  v[i,TB_cluster[i,"locID"]]=1
}


#x_pair
x_pair_full <- matrix(1, 
                      nrow = n_star, 
                      ncol = 7)  

counter<-1
for(j in 1:(n - 1)){
  for(k in (j + 1):n){
    x_pair_full[counter, 1] <- j  # or TB_cluster$ID[j] 
    x_pair_full[counter, 2] <- k  # or TB_cluster$ID[k]  
    # x_pair_full[counter, 3] intercept
    x_pair_full[counter, 4]<-abs(TB_cluster$age[j] - TB_cluster$age[k])
    x_pair_full[counter, 5]<-abs(sqrt((TB_cluster$Response_X_M[j]-TB_cluster$Response_X_M[k])^2+
                                        (TB_cluster$Response_Y_M[j] - TB_cluster$Response_Y_M[k])^2))
    x_pair_full[counter, 6]<-(TB_cluster$age[k] + TB_cluster$age[j]) # combined age x_ind
    counter<-counter + 1
    
  }
}

sd_x1 <- sd(x_pair_full[,4])
sd_x2 <- sd(x_pair_full[,5])
sd_ind <- sd(x_pair_full[,6])

x_pair_full[,4:6]<-scale(x_pair_full[,4:6])

z<-matrix(0,
          nrow = n_star,
          ncol = n)
counter<-1
for(j in 1:(n - 1)){
  for(k in (j + 1):n){
    
    z[counter, j]<-1
    z[counter, k]<-1
    counter<-counter +
      1
    
  }
}

Sigma_true<-tau2_true*exp(-phi_true*spatial_dists)

unique_locs = nrow(spatial_dists)
all_locs <- TB_cluster$locID
eta_simu_unique<-rmnorm(n = 1, 
                 mean = rep(0.00, times = unique_locs), 
                 varcov = Sigma_true)
eta_simu <- eta_simu_unique[all_locs] 

theta_simu_expanded <-  eta_simu +  rnorm(n = length(all_locs), mean = 0.00, 
                                   sd = sqrt(sigma2_zeta_true))

theta_simu_expanded <- theta_simu_expanded - mean(theta_simu_expanded)

mu_simu<-x_pair_full[,3:5]%*%beta_true +
  x_pair_full[,6, drop = FALSE]%*%gamma_true + 
  z%*%theta_simu_expanded


p_simu<-1.00/(1.00 + exp(-mu_simu))

clustered_indicators_full <-rbinom(n = n_star, 
                                   size = 1, 
                                   prob = p_simu)

x_pair_full[,7] = clustered_indicators_full #simulated outcome for the full dataset
x_pair_full <- data.frame(x_pair_full)
colnames(x_pair_full) <- c("id1", "id2","intercept", "age_diff", "distance","combined_age", "clustered_indicators")

#####################################

#
for(i in 1:M){
  start_time_subset <- Sys.time()
  group_index = i-1
  
  subL= TB_cluster[TB_cluster$group==group_index ,]
  subset_ID<- TB_cluster[TB_cluster$group==group_index ,]$ID
  subset_locs<- unique(TB_cluster[TB_cluster$group==group_index ,]$locID)
  subset_pairs <- subset(x_pair_full, id1 %in% subset_ID & id2 %in% subset_ID)
  
  #Setting the global data values:
  n<-nrow(subL)  #Number of individuals
  n_star<-0  #Number of paired responses
  for(j in 1:(n - 1)){ 
    for(k in (j + 1):n){
      n_star<-n_star + 
        1
    }
  }
  
  #unique locations
  locs <- subL %>% distinct(Response_X_M,Response_Y_M, .keep_all = TRUE) %>%
    select(Response_X_M, Response_Y_M)
  
  spatial_dists<-as.matrix(dist(locs,
                                diag = TRUE,
                                upper = TRUE))
  
  spatial_dists<-spatial_dists/max_dis
  
  #extract subset of V
  v_subset <- v[subset_ID, subset_locs, drop = FALSE]
  
  
  #x_pair
  x_pair <- subset_pairs[,3:5]
  
  
  #x_ind
  x_ind <- subset_pairs[,6, drop = FALSE]
  
  #z
  z<-matrix(0, 
            nrow = n_star, 
            ncol = n)  
  
  counter<-1
  
  for(j in 1:(n - 1)){
    for(k in (j + 1):n){
      
      z[counter, j]<-1
      z[counter, k]<-1
      counter<-counter +
        1
      
    }
  }
  
  clustered_indicators<- subset_pairs[,7]
  
  filepath = paste0("xxx/",jobid,"_subset",M,"_",group_index,".RData")
  
  save(list=c("clustered_indicators", "x_pair", "x_ind","z","spatial_dists","v_subset",
              "sd_x1","sd_x2","sd_ind"),
       file=filepath)
 
  end_time_subset <- Sys.time()
  elapsed_subset <- end_time_subset - start_time_subset
  cat(paste("Finished job",jobid, "_subset_",M,"_",group_index, "\n"))
  cat(paste("Runtime:", round(as.numeric(elapsed_subset, units = "mins"), 2), "mins\n"))
  
}

end_time <- Sys.time()
elapsed <- end_time - start_time

cat(paste("Finished job", jobid, "\n"))
cat(paste("Runtime:", round(as.numeric(elapsed, units = "mins"), 2), "mins\n"))

