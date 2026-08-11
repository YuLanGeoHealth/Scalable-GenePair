
set.seed(931014)
library(dplyr)
library(GenePair)
library(Matrix) #z matrix
####################################################

# how many controls per case?
m_controls_per_case <- 10   # for example, set m here


#Setting the global data values:
TB_cluster <- read.csv("xxx.csv")
TB_cluster <-TB_cluster[!is.na(TB_cluster$Response_X_M),]
TB_cluster <- TB_cluster %>% mutate(ID = row_number())

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
  
spatial_dists_ori<-as.matrix(dist(locs,
                              diag = TRUE,
                              upper = TRUE))


max_dis <- max(spatial_dists_ori) 
  
spatial_dists_full <- spatial_dists_ori/max_dis
  
#Unique locations ID
locs <- locs %>% 
  mutate(locID = row_number())
  
#add unique locID to all cases
TB_cluster <- inner_join(TB_cluster, locs, by = c("Response_X_M", "Response_Y_M"))
  
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
    x_pair_full[counter, 1] <- j 
    x_pair_full[counter, 2] <- k 
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

unique_locs = nrow(spatial_dists_full)

  
#clustered_indicators
clustered_indicators<-matrix(0, 
                             nrow = n_star, 
                             ncol = 1)  

counter<-1
for(j in 1:(n - 1)){
  for(k in (j + 1):n){
    clustered_indicators[counter, 1]<-ifelse(TB_cluster$clusterID[j] == TB_cluster$clusterID[k] & 
                                               !is.na(TB_cluster$clusterID[k])& 
                                               !is.na(TB_cluster$clusterID[j]), 1, 0)
    counter<-counter +1
    
  }
}


x_pair_full[,7] = clustered_indicators #simulated outcome for the full dataset
x_pair_full <- data.frame(x_pair_full)
colnames(x_pair_full) <- c("id1", "id2","intercept", "age_diff", "distance","combined_age", "clustered_indicators")

#####################################

####### CASE–CONTROL SAMPLING OF PAIRS ########################
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
  "../",
  "cc_input_m",m_controls_per_case,".RData")

save(list = c("clustered_indicators", "x_pair", "x_ind", "z",
              "spatial_dists", "v_subset", "sd_x1", "sd_x2","sd_ind"),
     file = filepath)




