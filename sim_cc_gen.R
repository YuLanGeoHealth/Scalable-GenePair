
# Simulate processing for Case-control approach
#simulation of scalable Genepair#
#require(devtools)
#devtools::install_github("warrenjl/GenePair")
args <- commandArgs(trailingOnly = TRUE)
jobid <- as.integer(args[1])

set.seed(jobid+1014)  # set a unique seed for each run

start_time <- Sys.time()
cat(paste("Running job", jobid, "with seed", jobid, "\n"))

library(mnormt)  #Multivariate normal distribution
library(dplyr)
library(Matrix) #z matrix

####################################################
#Setting the global data values:
TB_cluster <- read.csv("XXX.csv")
TB_cluster <-TB_cluster[!is.na(TB_cluster$Response_X_M),]
full_data <- TB_cluster


max <- full_data %>% distinct(Response_X_M,Response_Y_M, .keep_all = TRUE) %>%
  select(Response_X_M, Response_Y_M)

max_spatial_dists<-as.matrix(dist(max,
                                  diag = TRUE,
                                  upper = TRUE))
max_dis <- max(max_spatial_dists) 

####################################################
#Setting the values for the statistical model parameters: Fixed
results <- readRDS("XXX.rds") #GenePair results from the benchmark dataset

# Assign unique locations to the individuals

sample_size <- 60000 #MCMC sampling sizes
beta_1 <- mean(results$beta[1,seq(10001, sample_size,by = 5)]) #mean
beta_2 <- mean(results$beta[2,seq(10001, sample_size,by = 5)])
beta_3 <- mean(results$beta[3,seq(10001, sample_size,by = 5)])
beta_true <- c(beta_1,beta_2,beta_3) 

gamma_true <- mean(results$gamma[1,seq(10001, sample_size,by = 5)])

tau2_true <- mean(results$tau2[seq(10001, sample_size,by = 5)])
sigma2_zeta_true <- mean(results$sigma2_zeta[seq(10001, sample_size,by = 5)])
phi_true <- mean(results$phi[seq(10001, sample_size,by = 5)])

simu__data <- full_data[, c("age","Response_X_M","Response_Y_M")] %>%
  sample_n(700)
simu__data <- simu__data %>% 
  mutate(ID = row_number())

TB_cluster <- simu__data

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
  
spatial_dists_full<-as.matrix(dist(locs,
                              diag = TRUE,
                              upper = TRUE))
  
spatial_dists_full<-spatial_dists_full/max_dis
  
#Unique locations ID
locs <- locs %>% 
  mutate(locID = row_number())
  
#add unique locID to all cases
TB_cluster <- inner_join(TB_cluster, locs, by = c("Response_X_M", "Response_Y_M"))
# Select only distinct columns
TB_cluster <- TB_cluster %>%
  select(ID, age, Response_X_M, Response_Y_M, locID)
  
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

pairs <- combn(n, 2) # 2 x n_star; each column is one pair (j, k)

# one row per pair (i indexes rows), two entries per row
i <- rep(seq_len(ncol(pairs)), each = 2L)  # 1,1,2,2,3,3, ...
j <- c(pairs)                              # column-wise: j1,k1,j2,k2,...

z <- sparseMatrix(i = i, j = j, x = 1L,
                         dims = c(ncol(pairs), n))

Sigma_true<-tau2_true*exp(-phi_true*spatial_dists_full)

unique_locs = nrow(spatial_dists_full)

eta_simu_unique<-rmnorm(n = 1, 
                        mean = rep(0.00, times = unique_locs), 
                        varcov = Sigma_true)
all_locs <- TB_cluster$locID
eta_simu <- eta_simu_unique[all_locs] 
theta_simu_expanded <-  eta_simu +  rnorm(n = length(all_locs), mean = 0.00, 
                                          sd = sqrt(sigma2_zeta_true))

theta_simu_expanded <- theta_simu_expanded - mean(theta_simu_expanded)

mu_simu<-as.matrix(x_pair_full[,3:5]%*%beta_true +
  x_pair_full[,6, drop = FALSE]%*%gamma_true + 
  z%*%theta_simu_expanded)

p_simu<-1.00/(1.00 + exp(-mu_simu))
  
#simulated clusters
clustered_indicators_full <-rbinom(n = n_star, 
                             size = 1, 
                              prob = p_simu)

x_pair_full[,7] = clustered_indicators_full #simulated outcome for the full dataset
x_pair_full <- data.frame(x_pair_full)
colnames(x_pair_full) <- c("id1", "id2","intercept", "age_diff", "distance","combined_age", "clustered_indicators")

#####################################

####### CASE–CONTROL SAMPLING OF PAIRS ########################

# how many controls per case?
m_controls_per_case <- 10   

y_full <- x_pair_full$clustered_indicators

# indices of cases (Y_ij = 1) and controls (Y_ij = 0)
case_idx    <- which(y_full == 1)
control_idx <- which(y_full == 0)

# keep ALL cases
n_cases <- length(case_idx)

# total number of controls to sample
n_controls  <- n_cases * m_controls_per_case

control_idx_sub <- sample(control_idx, n_controls, replace = FALSE)

# combine and optionally shuffle
keep_idx <- c(case_idx, control_idx_sub)

# subset to case–control dataset
x_pair_cc <- x_pair_full[keep_idx, ]
clustered_indicators_cc <- x_pair_cc$clustered_indicators

# effective number of pairs in case–control data
n_star_cc <- nrow(x_pair_cc)
cat("Case–control dataset has", n_star_cc, "pairs\n")
cat("  Cases (Y=1):  ", sum(x_pair_cc$clustered_indicators == 1), "\n")
cat("  Controls (Y=0):", sum(x_pair_cc$clustered_indicators == 0), "\n")

############################################################
## Build objects for a single GenePair dataset (no subsets)
############################################################

# extract the individuals that appear in ANY kept pair
keep_ids <- unique(c(x_pair_cc$id1, x_pair_cc$id2))

subL <- TB_cluster[TB_cluster$ID %in% keep_ids, ]
subset_ID <- subL$ID
subset_locs <- unique(subL$locID)

n <- nrow(subL)   # number of individuals in this dataset

# subset spatial distance matrix and v consistently using original locID indexing
spatial_dists <- spatial_dists_full[subset_locs, subset_locs, drop = FALSE]
v_subset      <- v[subset_ID, subset_locs, drop = FALSE]


############################################################
## Construct x_pair, x_ind, and z for the retained pairs
############################################################

subset_pairs <- x_pair_cc[
  x_pair_cc$id1 %in% subset_ID &
    x_pair_cc$id2 %in% subset_ID, ]

# x_pair = the columns GenePair uses
x_pair <- as.matrix(subset_pairs[, c("intercept", "age_diff", "distance")])

# x_ind (combined age categorical)
x_ind <- as.matrix(subset_pairs[, "combined_age", drop = FALSE])

# Build z matrix
n_star <- nrow(subset_pairs)
z <- matrix(0, nrow = n_star, ncol = n)

for (row in 1:n_star) {
  i <- subset_pairs$id1[row]
  j <- subset_pairs$id2[row]
  # map IDs to row numbers of subL
  ii <- match(i, subset_ID)
  jj <- match(j, subset_ID)
  z[row, ii] <- 1
  z[row, jj] <- 1
}

clustered_indicators <- subset_pairs$clustered_indicators

filepath <- paste0(
  "xx/",
  "xx",jobid,".RData")


save(list = c("clustered_indicators", "x_pair", "x_ind", "z",
              "spatial_dists", "v_subset", "sd_x1", "sd_x2","sd_ind"),
     file = filepath)


end_time <- Sys.time()
elapsed <- end_time - start_time

cat(paste("Finished job", jobid, "\n"))
cat(paste("Runtime:", round(as.numeric(elapsed, units = "mins"), 2), "mins\n"))


