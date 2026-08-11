

set.seed(654) 
library("dplyr")
#library(GenePair)
library(ggdist)
library(tidyr)
# library(ggplot2)
TB_cluster <- read.csv("xxx.csv")
TB_cluster <-TB_cluster[!is.na(TB_cluster$Response_X_M),]
TB_cluster <- TB_cluster %>% mutate(ID = row_number())

#Unique locations (one for each individual)
locs <- TB_cluster %>% distinct(Response_X_M,Response_Y_M, .keep_all = TRUE) %>%
  select(Response_X_M, Response_Y_M)

spatial_dists_ori<-as.matrix(dist(locs,
                                  diag = TRUE,
                                  upper = TRUE))


max_dis <- max(spatial_dists_ori) 

df_random <- TB_cluster

set.seed(654) # set 2
df_random$group = 0
for (i in 1:5) {
  available <- which(df_random$group == 0)
  
  selected <- sample(available, size = 700, replace = FALSE)
  df_random$group[selected] <- i
}

table(df_random$group)

for(i in 0:5){
  
  group_index = i
  subL <- df_random[df_random$group==group_index ,]
  n<-nrow(subL)  #Number of individuals
  n_star<-0  #Number of paired responses
  for(j in 1:(n - 1)){ 
    for(k in (j + 1):n){
      n_star<-n_star + 1
    }
  }
  

#Unique locations (one for each individual)
locs <- subL %>% distinct(Response_X_M,Response_Y_M, .keep_all = TRUE) %>%
  select(Response_X_M, Response_Y_M)

spatial_dists<-as.matrix(dist(locs,
                              diag = TRUE,
                              upper = TRUE))
spatial_dists<-spatial_dists/max_dis

#Unique locations ID
locs <- locs %>% 
  mutate(locID = row_number())

#add unique locID to all cases
subL <- inner_join(subL, locs, by = c("Response_X_M", "Response_Y_M"))
# Select only distinct columns
subL <- subL %>%
  select(ID, age, Response_X_M, Response_Y_M, locID,group,clusterID)

v_subset<-matrix(0,
                 nrow=n,
                 ncol=nrow(locs))
for(i in 1:n){
  v_subset[i,subL[i,"locID"]]=1
}


#x_pair
x_pair<-matrix(1, 
               nrow = n_star, 
               ncol = 3)  

counter<-1
for(j in 1:(n - 1)){
  for(k in (j + 1):n){
    
    x_pair[counter, 2]<-abs(subL$age[j] - subL$age[k])
    #ncol =2 betas for age differences
    
    x_pair[counter, 3]<-abs(sqrt((subL$Response_X_M[j]-subL$Response_X_M[k])^2+
                                   (subL$Response_Y_M[j] - subL$Response_Y_M[k])^2))
    #ncol =3 betas for spatial differences
    
    counter<-counter + 1
    
  }
}
sd_x1 <- sd(x_pair[,2])
sd_x2 <- sd(x_pair[,3])


x_pair[,2:3]<-scale(x_pair[,2:3])

#x_ind, only for continues variables
x1<-subL$age
x_ind<-matrix(0.00, 
              nrow = n_star, 
              ncol = 1)  
counter<-1
for(j in 1:(n - 1)){
  for(k in (j + 1):n){
    
    x_ind[counter, 1]<- x1[k] + x1[j]
    counter<-counter + 1
    
  }
}
sd_ind <- sd(x_ind)
x_ind<-scale(x_ind)

#z
z<-matrix(0, 
          nrow = n_star, 
          ncol = n) 
counter<-1
for(j in 1:(n - 1)){
  for(k in (j + 1):n){
    
    z[counter, j]<-1
    z[counter, k]<-1
    counter<-counter + 1
    
  }
}

#clustered_indicators
clustered_indicators<-matrix(0, 
                             nrow = n_star, 
                             ncol = 1)  

counter<-1
for(j in 1:(n - 1)){
  for(k in (j + 1):n){
    clustered_indicators[counter, 1]<-ifelse(subL$clusterID[j] == subL$clusterID[k] & 
                                               !is.na(subL$clusterID[k])& 
                                               !is.na(subL$clusterID[j]), 1, 0)
    counter<-counter +1
    
  }
}

filepath = paste0("xxx",group_index,".RData")

save(list=c("clustered_indicators", "x_pair", "x_ind","z","spatial_dists","v_subset", "sd_x1", "sd_x2","sd_ind"),
     file=filepath)
}

