###################################
########## Analysis file ##########
###################################

####################
##### Packages #####
####################

## "ggh4x" package.
# Check whether "ggh4x" package is installed and install if necessary.
if(!require(ggh4x)) install.packages("ggh4x")
# Push "ggh4x" package to library.
library(ggh4x)

## "ggarrange" package.
# Check whether "ggarrange" package is installed and install if necessary.
if(!require(ggarrange)) install.packages("ggarrange")
# Push "ggarrange" package to library.
library(ggarrange)

## "ergm" package.
# Check whether "ergm" package is installed and install if necessary.
if(!require(ergm)) install.packages("ergm")
# Push "ergm" package to library.
library(ergm)

## "data.table" package.
# Check whether "data.table" package is installed and install if necessary.
if(!require(data.table)) install.packages("data.table")
# Push "data.table" package to library.
library(data.table)

## "plyr" package.
# Check whether "plyr" package is installed and install if necessary.
if(!require(plyr)) install.packages("plyr")
# Push "plyr" package to library.
library(plyr)

## "ggplot2" package.
# Check whether "ggplot2" package is installed and install if necessary.
if(!require(ggplot2)) install.packages("ggplot2")
# Push "ggplot2" package to library.
library(ggplot2)

## "sna" package.
# Check whether "sna" package is installed and install if necessary.
if(!require(sna)) install.packages("sna")
# Push "sna" package to library.
library(sna)

#####################
##### Load data #####
#####################

### Set working directory.
setwd("D:\\Academia\\Master MSBBSS\\Thesis\\Research archive JWG Simons\\Ad A. Data storage\\2. Data files\\Empirical data")
### Import network files. 
ties_het0 <- as.matrix(read.table("not8.net"))
ties_het1 <- as.matrix(read.table("not16.net"))
### Import covariate files.
cov_het0 <- read.table("cov8.dat", col.names = c("girl", "identcGMC", "identcCMC", "identcCM", "percgirl"))
cov_het1 <- read.table("cov16.dat", col.names = c("girl", "identcGMC", "identcCMC", "identcCM", "percgirl"))

### Set working directory.
setwd("D:\\Academia\\Master MSBBSS\\Thesis\\Research archive JWG Simons\\Ad A. Data storage\\3. Computer code\\Empirical computer code")
### Import "myEnvironment.RData".
load("myEnvironment.RData")

### Set working directory.
setwd("D:\\Academia\\Master MSBBSS\\Thesis\\Research archive JWG Simons\\Ad A. Data storage\\2. Data files\\Simulation study data")
### Import "simdata.RData".
load("simdata.RData")

setwd("D:\\Academia\\Master MSBBSS\\Thesis\\Simulation folder\\simdata files\\New folder 2")
load("simdata_1_1.RData")

#############################################################################################
##### Five point summary of the 8th and 16th Vermeij classroom relationship nominations #####
#############################################################################################

net8 <- network(ties_het0)
net16 <- network(ties_het1)

network.size(net8) # Network 8th size.
network.size(net16) # Network 16th size.
gden(net8) # Network 8th density.
gden(net16) # Network 16th density.
components(net8) # Network 8th components.
components(net16) # Network 16th components.
net8_geodist <- geodist(component.largest(net8, result = "graph"))
max(net8_geodist$gdist) # Diameter 8th network.
net16_geodist <- geodist(component.largest(net16, result = "graph"))
max(net16_geodist$gdist) # Diameter 16th network.
gtrans(net8) # Clustering coefficient 8th network.
gtrans(net16) # Clustering coefficient 16th network.
length(isolates(net8)) # Number of isolates in 8th network.
length(isolates(net16)) # Number of isolates in 16th network.

#######################################################
##### Target statistics of population level ERGMs #####
#######################################################

sample(true_ergm_results_het0$target.stats)
sample(true_ergm_results_het1$target.stats)

##########################################################################################################################
##### Mean (%) convergence failure for the structural, true and dense triadic models by samplesize and heterogeneity #####
##########################################################################################################################

fail_1 <- vector(length = 124)
fail_2 <- vector(length = 124)
fail_3 <- vector(length = 124)
fail_4 <- vector(length = 124)
fail_5 <- vector(length = 124)
fail_6 <- vector(length = 124)
fail_7 <- vector(length = 124)
fail_8 <- vector(length = 124)
fail_9 <- vector(length = 124)
fail_10 <- vector(length = 124)
fail_11 <- vector(length = 124)
fail_12 <- vector(length = 124)

for (i in 1:124){
  
  fail_1[i] <- sum_cell_25[[i]]$gwodegree_results_het0$next_counter
  fail_2[i] <- sum_cell_25[[i]]$cov_results_het0$next_counter
  fail_3[i] <- sum_cell_25[[i]]$dtriad_results_het0$next_counter
  fail_4[i] <- sum_cell_25[[i]]$gwodegree_results_het1$next_counter
  fail_5[i] <- sum_cell_25[[i]]$cov_results_het1$next_counter
  fail_6[i] <- sum_cell_25[[i]]$dtriad_results_het1$next_counter
  fail_7[i] <- sum_cell_75[[i]]$gwodegree_results_het0$next_counter
  fail_8[i] <- sum_cell_75[[i]]$cov_results_het0$next_counter
  fail_9[i] <- sum_cell_75[[i]]$dtriad_results_het0$next_counter
  fail_10[i] <- sum_cell_75[[i]]$gwodegree_results_het1$next_counter
  fail_11[i] <- sum_cell_75[[i]]$cov_results_het1$next_counter
  fail_12[i] <- sum_cell_75[[i]]$dtriad_results_het1$next_counter
  
}

## Mean.
mean(fail_1)
mean(fail_2)
mean(fail_3)
mean(fail_4)
mean(fail_5)
mean(fail_6)
mean(fail_7)
mean(fail_8)
mean(fail_9)
mean(fail_10)
mean(fail_11)
mean(fail_12)

## Percentage.
mean(fail_1) / 25
mean(fail_2) / 25
mean(fail_3) / 25
mean(fail_4) / 25
mean(fail_5) / 25
mean(fail_6) / 25
mean(fail_7) / 75
mean(fail_8) / 75
mean(fail_9) / 75 
mean(fail_10) / 75
mean(fail_11) / 75
mean(fail_12) / 75

#########################
##### Bias and RMSE #####
#########################

### Bias and RMSE meta ERGM.

## Bias.

bias_true_25_hom <- matrix(nrow = 124, ncol = 6)
for (i in 1 : nrow(bias_true_25_hom)){
  
  bias_true_25_hom[i, ] <- (c(sum_cell_25[[1]]$cov_results_het0$ma_model[[1]]$b, sum_cell_25[[1]]$cov_results_het0$ma_model[[2]]$b,
                              sum_cell_25[[1]]$cov_results_het0$ma_model[[3]]$b, sum_cell_25[[1]]$cov_results_het0$ma_model[[4]]$b,
                              sum_cell_25[[1]]$cov_results_het0$ma_model[[5]]$b, sum_cell_25[[1]]$cov_results_het0$ma_model[[6]]$b) - true_ergm_results_het0$coef)
  
}

bias_struct_25_hom <- matrix(nrow = 124, ncol = 4)
for (i in 1 : nrow(bias_struct_25_hom)){
  
  bias_struct_25_hom[i, ] <- (c(sum_cell_25[[i]]$gwodegree_results_het0$ma_model[[1]]$b, sum_cell_25[[i]]$gwodegree_results_het0$ma_model[[2]]$b,
                                sum_cell_25[[i]]$gwodegree_results_het0$ma_model[[3]]$b, sum_cell_25[[i]]$gwodegree_results_het0$ma_model[[4]]$b) - true_ergm_results_het0$coef[1:4]) 
  
}

bias_dtriad_25_hom <- matrix(nrow = 124, ncol = 6)
for (i in 1 : nrow(bias_dtriad_25_hom)){
  
  bias_dtriad_25_hom[i, ] <- (c(sum_cell_25[[i]]$dtriad_results_het0$ma_model[[1]]$b, sum_cell_25[[i]]$dtriad_results_het0$ma_model[[2]]$b,
                                sum_cell_25[[i]]$dtriad_results_het0$ma_model[[3]]$b, sum_cell_25[[i]]$dtriad_results_het0$ma_model[[4]]$b,
                                sum_cell_25[[i]]$dtriad_results_het0$ma_model[[5]]$b, sum_cell_25[[i]]$dtriad_results_het0$ma_model[[6]]$b) - true_ergm_results_het0$coef) 
  
}

bias_true_75_hom <- matrix(nrow = 124, ncol = 6)
for (i in 1 : nrow(bias_true_75_hom)){
  
  bias_true_75_hom[i, ] <- (c(sum_cell_75[[i]]$cov_results_het0$ma_model[[1]]$b, sum_cell_75[[i]]$cov_results_het0$ma_model[[2]]$b,
                              sum_cell_75[[i]]$cov_results_het0$ma_model[[3]]$b, sum_cell_75[[i]]$cov_results_het0$ma_model[[4]]$b,
                              sum_cell_75[[i]]$cov_results_het0$ma_model[[5]]$b, sum_cell_75[[i]]$cov_results_het0$ma_model[[6]]$b) - true_ergm_results_het0$coef) 
  
}

bias_struct_75_hom <- matrix(nrow = 124, ncol = 4)
for (i in 1 : nrow(bias_struct_75_hom)){
  
  bias_struct_75_hom[i, ] <- (c(sum_cell_75[[i]]$gwodegree_results_het0$ma_model[[1]]$b, sum_cell_75[[i]]$gwodegree_results_het0$ma_model[[2]]$b,
                                sum_cell_75[[i]]$gwodegree_results_het0$ma_model[[3]]$b, sum_cell_75[[i]]$gwodegree_results_het0$ma_model[[4]]$b) - true_ergm_results_het0$coef[1:4]) 
  
}

bias_dtriad_75_hom <- matrix(nrow = 124, ncol = 6)
for (i in 1 : nrow(bias_dtriad_75_hom)){
  
  bias_dtriad_75_hom[i, ] <- (c(sum_cell_75[[i]]$dtriad_results_het0$ma_model[[1]]$b, sum_cell_75[[i]]$dtriad_results_het0$ma_model[[2]]$b,
                                sum_cell_75[[i]]$dtriad_results_het0$ma_model[[3]]$b, sum_cell_75[[i]]$dtriad_results_het0$ma_model[[4]]$b,
                                sum_cell_75[[i]]$dtriad_results_het0$ma_model[[5]]$b, sum_cell_75[[i]]$dtriad_results_het0$ma_model[[6]]$b) - true_ergm_results_het0$coef)
  
}

bias_true_25_het0 <- matrix(nrow = 124, ncol = 6)
for (i in 1 : nrow(bias_true_25_het0)){
  
  bias_true_25_het0[i, ] <- (c(sum_cell_25[[i]]$cov_results_het1$ma_model[[1]]$b, sum_cell_25[[i]]$cov_results_het1$ma_model[[2]]$b,
                               sum_cell_25[[i]]$cov_results_het1$ma_model[[3]]$b, sum_cell_25[[i]]$cov_results_het1$ma_model[[4]]$b,
                               sum_cell_25[[i]]$cov_results_het1$ma_model[[5]]$b, sum_cell_25[[i]]$cov_results_het1$ma_model[[6]]$b) - true_ergm_results_het0$coef) 
  
}

bias_true_25_het1 <- matrix(nrow = 124, ncol = 6)
for (i in 1 : nrow(bias_true_25_het1)){
  
  bias_true_25_het1[i, ] <- (c(sum_cell_25[[i]]$cov_results_het1$ma_model[[1]]$b, sum_cell_25[[i]]$cov_results_het1$ma_model[[2]]$b,
                               sum_cell_25[[i]]$cov_results_het1$ma_model[[3]]$b, sum_cell_25[[i]]$cov_results_het1$ma_model[[4]]$b,
                               sum_cell_25[[i]]$cov_results_het1$ma_model[[5]]$b, sum_cell_25[[i]]$cov_results_het1$ma_model[[6]]$b) - true_ergm_results_het1$coef)
  
}

bias_struct_25_het0 <- matrix(nrow = 124, ncol = 4)
for (i in 1 : nrow(bias_struct_25_het0)){
  
  bias_struct_25_het0[i, ] <- (c(sum_cell_25[[i]]$gwodegree_results_het1$ma_model[[1]]$b, sum_cell_25[[i]]$gwodegree_results_het1$ma_model[[2]]$b,
                                 sum_cell_25[[i]]$gwodegree_results_het1$ma_model[[3]]$b, sum_cell_25[[i]]$gwodegree_results_het1$ma_model[[4]]$b) - true_ergm_results_het0$coef[1:4])
  
}

bias_struct_25_het1 <- matrix(nrow = 124, ncol = 4)
for (i in 1 : nrow(bias_struct_25_het1)){
  
  bias_struct_25_het1[i, ] <- (c(sum_cell_25[[i]]$gwodegree_results_het1$ma_model[[1]]$b, sum_cell_25[[i]]$gwodegree_results_het1$ma_model[[2]]$b,
                                 sum_cell_25[[i]]$gwodegree_results_het1$ma_model[[3]]$b, sum_cell_25[[i]]$gwodegree_results_het1$ma_model[[4]]$b) - true_ergm_results_het1$coef[1:4]) 
  
}

bias_dtriad_25_het0 <- matrix(nrow = 124, ncol = 6)
for (i in 1 : nrow(bias_dtriad_25_het0)){
  
  bias_dtriad_25_het0[i, ] <- (c(sum_cell_25[[i]]$dtriad_results_het1$ma_model[[1]]$b, sum_cell_25[[i]]$dtriad_results_het1$ma_model[[2]]$b,
                                 sum_cell_25[[i]]$dtriad_results_het1$ma_model[[3]]$b, sum_cell_25[[i]]$dtriad_results_het1$ma_model[[4]]$b,
                                 sum_cell_25[[i]]$dtriad_results_het1$ma_model[[5]]$b, sum_cell_25[[i]]$dtriad_results_het1$ma_model[[6]]$b) - true_ergm_results_het0$coef)
  
}

bias_dtriad_25_het1 <- matrix(nrow = 124, ncol = 6)
for (i in 1 : nrow(bias_dtriad_25_het1)){
  
  bias_dtriad_25_het1[i, ] <- (c(sum_cell_25[[i]]$dtriad_results_het1$ma_model[[1]]$b, sum_cell_25[[i]]$dtriad_results_het1$ma_model[[2]]$b,
                                 sum_cell_25[[i]]$dtriad_results_het1$ma_model[[3]]$b, sum_cell_25[[i]]$dtriad_results_het1$ma_model[[4]]$b,
                                 sum_cell_25[[i]]$dtriad_results_het1$ma_model[[5]]$b, sum_cell_25[[i]]$dtriad_results_het1$ma_model[[6]]$b) - true_ergm_results_het1$coef)
  
}

bias_true_75_het0 <- matrix(nrow = 124, ncol = 6)
for (i in 1 : nrow(bias_true_75_het0)){
  
  bias_true_75_het0[i, ] <- (c(sum_cell_75[[i]]$cov_results_het1$ma_model[[1]]$b, sum_cell_75[[i]]$cov_results_het1$ma_model[[2]]$b,
                               sum_cell_75[[i]]$cov_results_het1$ma_model[[3]]$b, sum_cell_75[[i]]$cov_results_het1$ma_model[[4]]$b,
                               sum_cell_75[[i]]$cov_results_het1$ma_model[[5]]$b, sum_cell_75[[i]]$cov_results_het1$ma_model[[6]]$b) - true_ergm_results_het0$coef) 
  
}

bias_true_75_het1 <- matrix(nrow = 124, ncol = 6)
for (i in 1 : nrow(bias_true_75_het1)){
  
  bias_true_75_het1[i, ] <- (c(sum_cell_75[[i]]$cov_results_het1$ma_model[[1]]$b, sum_cell_75[[i]]$cov_results_het1$ma_model[[2]]$b,
                               sum_cell_75[[i]]$cov_results_het1$ma_model[[3]]$b, sum_cell_75[[i]]$cov_results_het1$ma_model[[4]]$b,
                               sum_cell_75[[i]]$cov_results_het1$ma_model[[5]]$b, sum_cell_75[[i]]$cov_results_het1$ma_model[[6]]$b) - true_ergm_results_het1$coef)
  
}

bias_struct_75_het0 <- matrix(nrow = 124, ncol = 4)
for (i in 1 : nrow(bias_struct_75_het0)){
  
  bias_struct_75_het0[i, ] <- (c(sum_cell_75[[i]]$gwodegree_results_het1$ma_model[[1]]$b, sum_cell_75[[i]]$gwodegree_results_het1$ma_model[[2]]$b,
                                 sum_cell_75[[i]]$gwodegree_results_het1$ma_model[[3]]$b, sum_cell_75[[i]]$gwodegree_results_het1$ma_model[[4]]$b) - true_ergm_results_het0$coef[1:4])
  
}

bias_struct_75_het1 <- matrix(nrow = 124, ncol = 4)
for (i in 1 : nrow(bias_struct_75_het1)){
  
  bias_struct_75_het1[i, ] <- (c(sum_cell_75[[i]]$gwodegree_results_het1$ma_model[[1]]$b, sum_cell_75[[i]]$gwodegree_results_het1$ma_model[[2]]$b,
                                 sum_cell_75[[i]]$gwodegree_results_het1$ma_model[[3]]$b, sum_cell_75[[i]]$gwodegree_results_het1$ma_model[[4]]$b) - true_ergm_results_het1$coef[1:4])
  
}

bias_dtriad_75_het0 <- matrix(nrow = 124, ncol = 6)
for (i in 1 : nrow(bias_dtriad_75_het0)){
  
  bias_dtriad_75_het0[i, ] <- (c(sum_cell_75[[i]]$dtriad_results_het1$ma_model[[1]]$b, sum_cell_75[[i]]$dtriad_results_het1$ma_model[[2]]$b,
                                 sum_cell_75[[i]]$dtriad_results_het1$ma_model[[3]]$b, sum_cell_75[[i]]$dtriad_results_het1$ma_model[[4]]$b,
                                 sum_cell_75[[i]]$dtriad_results_het1$ma_model[[5]]$b, sum_cell_75[[i]]$dtriad_results_het1$ma_model[[6]]$b) - true_ergm_results_het0$coef) 
  
}

bias_dtriad_75_het1 <- matrix(nrow = 124, ncol = 6)
for (i in 1 : nrow(bias_dtriad_75_het1)){
  
  bias_dtriad_75_het1[i, ] <- (c(sum_cell_75[[i]]$dtriad_results_het1$ma_model[[1]]$b, sum_cell_75[[i]]$dtriad_results_het1$ma_model[[2]]$b,
                                 sum_cell_75[[i]]$dtriad_results_het1$ma_model[[3]]$b, sum_cell_75[[i]]$dtriad_results_het1$ma_model[[4]]$b,
                                 sum_cell_75[[i]]$dtriad_results_het1$ma_model[[5]]$b, sum_cell_75[[i]]$dtriad_results_het1$ma_model[[6]]$b) - true_ergm_results_het1$coef)
  
}

colnames(bias_dtriad_25_hom) <- c("edges", "mutual", "gwesp", "gwodegree", "nodeocov", "nodematch")
colnames(bias_dtriad_25_het0) <- c("edges", "mutual", "gwesp", "gwodegree", "nodeocov", "nodematch")
colnames(bias_dtriad_25_het1) <- c("edges", "mutual", "gwesp", "gwodegree", "nodeocov", "nodematch")
colnames(bias_dtriad_75_hom) <- c("edges", "mutual", "gwesp", "gwodegree", "nodeocov", "nodematch")
colnames(bias_dtriad_75_het0) <- c("edges", "mutual", "gwesp", "gwodegree", "nodeocov", "nodematch")
colnames(bias_dtriad_75_het1) <- c("edges", "mutual", "gwesp", "gwodegree", "nodeocov", "nodematch")

bias_dtriad_25_hom <- melt(bias_dtriad_25_hom)
bias_dtriad_25_het0 <- melt(bias_dtriad_25_het0)
bias_dtriad_25_het1 <- melt(bias_dtriad_25_het1)
bias_dtriad_75_hom <- melt(bias_dtriad_75_hom)
bias_dtriad_75_het0 <- melt(bias_dtriad_75_het0)
bias_dtriad_75_het1 <- melt(bias_dtriad_75_het1)

bias_dtriad_25_hom <- cbind(sampsize = "25", het = "Homogeneous sample", comp = "Population model 1", bias_dtriad_25_hom)
bias_dtriad_25_het0 <- cbind(sampsize = "25", het = "Heterogeneous sample", comp = "Population model 1", bias_dtriad_25_het0)
bias_dtriad_25_het1 <- cbind(sampsize = "25", het = "Heterogeneous sample", comp = "Population model 2", bias_dtriad_25_het1)
bias_dtriad_75_hom <- cbind(sampsize = "75", het = "Homogeneous sample", comp = "Population model 1", bias_dtriad_75_hom)
bias_dtriad_75_het0 <- cbind(sampsize = "75", het = "Heterogeneous sample", comp = "Population model 1", bias_dtriad_75_het0)
bias_dtriad_75_het1 <- cbind(sampsize = "75", het = "Heterogeneous sample", comp = "Population model 2", bias_dtriad_75_het1)

colnames(bias_struct_25_hom) <- c("edges", "mutual", "gwesp", "gwodegree")
colnames(bias_struct_25_het0) <- c("edges", "mutual", "gwesp", "gwodegree")
colnames(bias_struct_25_het1) <- c("edges", "mutual", "gwesp", "gwodegree")
colnames(bias_struct_75_hom) <- c("edges", "mutual", "gwesp", "gwodegree")
colnames(bias_struct_75_het0) <- c("edges", "mutual", "gwesp", "gwodegree")
colnames(bias_struct_75_het1) <- c("edges", "mutual", "gwesp", "gwodegree")

bias_struct_25_hom <- melt(bias_struct_25_hom)
bias_struct_25_het0 <- melt(bias_struct_25_het0)
bias_struct_25_het1 <- melt(bias_struct_25_het1)
bias_struct_75_hom <- melt(bias_struct_75_hom)
bias_struct_75_het0 <- melt(bias_struct_75_het0)
bias_struct_75_het1 <- melt(bias_struct_75_het1)

bias_struct_25_hom <- cbind(sampsize = "25", het = "Homogeneous sample", comp = "Population model 1", bias_struct_25_hom)
bias_struct_25_het0 <- cbind(sampsize = "25", het = "Heterogeneous sample", comp = "Population model 1", bias_struct_25_het0)
bias_struct_25_het1 <- cbind(sampsize = "25", het = "Heterogeneous sample", comp = "Population model 2", bias_struct_25_het1)
bias_struct_75_hom <- cbind(sampsize = "75", het = "Homogeneous sample", comp = "Population model 1", bias_struct_75_hom)
bias_struct_75_het0 <- cbind(sampsize = "75", het = "Heterogeneous sample", comp = "Population model 1", bias_struct_75_het0)
bias_struct_75_het1 <- cbind(sampsize = "75", het = "Heterogeneous sample", comp = "Population model 2", bias_struct_75_het1)

colnames(bias_true_25_hom) <- c("edges", "mutual", "gwesp", "gwodegree", "nodeocov", "nodematch")
colnames(bias_true_25_het0) <- c("edges", "mutual", "gwesp", "gwodegree", "nodeocov", "nodematch")
colnames(bias_true_25_het1) <- c("edges", "mutual", "gwesp", "gwodegree", "nodeocov", "nodematch")
colnames(bias_true_75_hom) <- c("edges", "mutual", "gwesp", "gwodegree", "nodeocov", "nodematch")
colnames(bias_true_75_het0) <- c("edges", "mutual", "gwesp", "gwodegree", "nodeocov", "nodematch")
colnames(bias_true_75_het1) <- c("edges", "mutual", "gwesp", "gwodegree", "nodeocov", "nodematch")

bias_true_25_hom <- melt(bias_true_25_hom)
bias_true_25_het0 <- melt(bias_true_25_het0)
bias_true_25_het1 <- melt(bias_true_25_het1)
bias_true_75_hom <- melt(bias_true_75_hom)
bias_true_75_het0 <- melt(bias_true_75_het0)
bias_true_75_het1 <- melt(bias_true_75_het1)

bias_true_25_hom <- cbind(sampsize = '25', het = "Homogeneous sample", comp = "Population model 1", bias_true_25_hom)
bias_true_25_het0 <- cbind(sampsize = "25", het = "Heterogeneous sample", comp = "Population model 1", bias_true_25_het0)
bias_true_25_het1 <- cbind(sampsize = "25", het = "Heterogeneous sample", comp = "Population model 2", bias_true_25_het1)
bias_true_75_hom <- cbind(sampsize = "75", het = "Homogeneous sample", comp = "Population model 1", bias_true_75_hom)
bias_true_75_het0 <- cbind(sampsize = "75", het = "Heterogeneous sample", comp = "Population model 1", bias_true_75_het0)
bias_true_75_het1 <- cbind(sampsize = "75", het = "Heterogeneous sample", comp = "Population model 2", bias_true_75_het1)

dat_bias_dtriad <- rbind.fill(data.frame(bias_dtriad_25_hom), data.frame(bias_dtriad_25_het0),  data.frame(bias_dtriad_25_het1),
                              data.frame(bias_dtriad_75_hom), data.frame(bias_dtriad_75_het0), data.frame(bias_dtriad_75_het1))

dat_bias_true <- rbind.fill(data.frame(bias_true_25_hom), data.frame(bias_true_25_het0), data.frame(bias_true_25_het1),
                            data.frame(bias_true_75_hom), data.frame(bias_true_75_het0), data.frame(bias_true_75_het1))

dat_bias_structural <- rbind.fill(data.frame(bias_struct_25_hom), data.frame(bias_struct_25_het0), data.frame(bias_struct_25_het1),
                                  data.frame(bias_struct_75_hom), data.frame(bias_struct_75_het0), data.frame(bias_struct_75_het1))

dat_bias_dtriad <- dat_bias_dtriad[dat_bias_dtriad[, "value"] > - 5, ] # None removed.
dat_bias_dtriad <- dat_bias_dtriad[dat_bias_dtriad[, "value"] < 5, ] # 2944, 2951 removed.
dat_bias_true <- dat_bias_true[dat_bias_true[, "value"] > - 5, ] # 4011 removed.
dat_bias_true <- dat_bias_true[dat_bias_true[, "value"] < 5, ] # 714, 3267 removed.
dat_bias_structural <- dat_bias_structural[dat_bias_structural[, "value"] > - 5, ] # 1313 removed.
dat_bias_structural <- dat_bias_structural[dat_bias_structural[, "value"] < 5, ] # 817, 990, 1486, 2957 removed.

dat_bias_dtriad$het <- factor(dat_bias_dtriad$het, levels = c('Homogeneous sample','Heterogeneous sample'), ordered = TRUE)
dat_bias_true$het <- factor(dat_bias_true$het, levels = c('Homogeneous sample','Heterogeneous sample'), ordered = TRUE)
dat_bias_structural$het <- factor(dat_bias_structural$het, levels = c('Homogeneous sample','Heterogeneous sample'), ordered = TRUE)
dat_bias_dtriad$X2 <- factor(dat_bias_dtriad$X2, levels = c("edges", "mutual", "gwesp", "gwodegree", "nodeocov", "nodematch")) 
dat_bias_true$X2 <- factor(dat_bias_true$X2, levels = c("edges", "mutual", "gwesp", "gwodegree", "nodeocov", "nodematch")) 
dat_bias_structural$X2 <- factor(dat_bias_structural$X2, levels = c("edges", "mutual", "gwesp", "gwodegree")) 

ggplot(dat_bias_dtriad, aes(x = X2, y = value, fill = sampsize)) +
  geom_boxplot() +
  facet_nested(~ het + comp) +
  labs(x = "Parameter", y = "Bias") +
  scale_fill_discrete(name = "Sample size") +
  theme(axis.text.x = element_text(angle = 40, vjust = 1, hjust = 1))

ggplot(dat_bias_true, aes(x = X2, y = value, fill = sampsize)) +
  geom_boxplot() +
  facet_nested(~ het + comp) +
  labs(x = "Parameter", y = "Bias") +
  scale_fill_discrete(name = "Sample size") +
  theme(axis.text.x = element_text(angle = 40, vjust = 1, hjust = 1)) 

ggplot(dat_bias_structural, aes(x = X2, y = value, fill = sampsize)) +
  geom_boxplot() +
  facet_nested(~ het + comp) +
  labs(x = "Parameter", y = "Bias") +
  scale_fill_discrete(name = "Sample size") +
  theme(axis.text.x = element_text(angle = 40, vjust = 1, hjust = 1)) 

## RMSE.
rmse_true_25_hom <- matrix(nrow = 124, ncol = 6)
for (i in 1 : nrow(rmse_true_25_hom)){
  
  rmse_true_25_hom[i, ] <- sqrt(((c(sum_cell_25[[i]]$cov_results_het0$ma_model[[1]]$b, sum_cell_25[[i]]$cov_results_het0$ma_model[[2]]$b,
                                    sum_cell_25[[i]]$cov_results_het0$ma_model[[3]]$b, sum_cell_25[[i]]$cov_results_het0$ma_model[[4]]$b,
                                    sum_cell_25[[i]]$cov_results_het0$ma_model[[5]]$b, sum_cell_25[[i]]$cov_results_het0$ma_model[[6]]$b) - true_ergm_results_het0$coef) / true_ergm_results_het0$coef) ** 2)
  
}

rmse_struct_25_hom <- matrix(nrow = 124, ncol = 4)
for (i in 1 : nrow(rmse_struct_25_hom)){
  
  rmse_struct_25_hom[i, ] <- sqrt(((c(sum_cell_25[[i]]$gwodegree_results_het0$ma_model[[1]]$b, sum_cell_25[[i]]$gwodegree_results_het0$ma_model[[2]]$b,
                                      sum_cell_25[[i]]$gwodegree_results_het0$ma_model[[3]]$b, sum_cell_25[[i]]$gwodegree_results_het0$ma_model[[4]]$b) - true_ergm_results_het0$coef[1:4]) / true_ergm_results_het0$coef[1:4]) ** 2)
  
}

rmse_dtriad_25_hom <- matrix(nrow = 124, ncol = 6)
for (i in 1 : nrow(rmse_dtriad_25_hom)){
  
  rmse_dtriad_25_hom[i, ] <- sqrt(((c(sum_cell_25[[i]]$dtriad_results_het0$ma_model[[1]]$b, sum_cell_25[[i]]$dtriad_results_het0$ma_model[[2]]$b,
                                      sum_cell_25[[i]]$dtriad_results_het0$ma_model[[3]]$b, sum_cell_25[[i]]$dtriad_results_het0$ma_model[[4]]$b,
                                      sum_cell_25[[i]]$dtriad_results_het0$ma_model[[5]]$b, sum_cell_25[[i]]$dtriad_results_het0$ma_model[[6]]$b) - true_ergm_results_het0$coef) / true_ergm_results_het0$coef) ** 2)
  
}

rmse_true_75_hom <- matrix(nrow = 124, ncol = 6)
for (i in 1 : nrow(rmse_true_75_hom)){
  
  rmse_true_75_hom[i, ] <- sqrt(((c(sum_cell_75[[i]]$cov_results_het0$ma_model[[1]]$b, sum_cell_75[[i]]$cov_results_het0$ma_model[[2]]$b,
                                    sum_cell_75[[i]]$cov_results_het0$ma_model[[3]]$b, sum_cell_75[[i]]$cov_results_het0$ma_model[[4]]$b,
                                    sum_cell_75[[i]]$cov_results_het0$ma_model[[5]]$b, sum_cell_75[[i]]$cov_results_het0$ma_model[[6]]$b) - true_ergm_results_het0$coef) / true_ergm_results_het0$coef)  ** 2)
  
}

rmse_struct_75_hom <- matrix(nrow = 124, ncol = 4)
for (i in 1 : nrow(rmse_struct_75_hom)){
  
  rmse_struct_75_hom[i, ] <- sqrt(((c(sum_cell_75[[i]]$gwodegree_results_het0$ma_model[[1]]$b, sum_cell_75[[i]]$gwodegree_results_het0$ma_model[[2]]$b,
                                      sum_cell_75[[i]]$gwodegree_results_het0$ma_model[[3]]$b, sum_cell_75[[i]]$gwodegree_results_het0$ma_model[[4]]$b) - true_ergm_results_het0$coef[1:4]) / true_ergm_results_het0$coef[1:4]) ** 2)
  
}

rmse_dtriad_75_hom <- matrix(nrow = 124, ncol = 6)
for (i in 1 : nrow(rmse_dtriad_75_hom)){
  
  rmse_dtriad_75_hom[i, ] <- sqrt(((c(sum_cell_75[[i]]$dtriad_results_het0$ma_model[[1]]$b, sum_cell_75[[i]]$dtriad_results_het0$ma_model[[2]]$b,
                                      sum_cell_75[[i]]$dtriad_results_het0$ma_model[[3]]$b, sum_cell_75[[i]]$dtriad_results_het0$ma_model[[4]]$b,
                                      sum_cell_75[[i]]$dtriad_results_het0$ma_model[[5]]$b, sum_cell_75[[i]]$dtriad_results_het0$ma_model[[6]]$b) - true_ergm_results_het0$coef) / true_ergm_results_het0$coef) ** 2)
  
}

rmse_true_25_het0 <- matrix(nrow = 124, ncol = 6)
for (i in 1 : nrow(rmse_true_25_het0)){
  
  rmse_true_25_het0[i, ] <- sqrt(((c(sum_cell_25[[i]]$cov_results_het1$ma_model[[1]]$b, sum_cell_25[[i]]$cov_results_het1$ma_model[[2]]$b,
                                     sum_cell_25[[i]]$cov_results_het1$ma_model[[3]]$b, sum_cell_25[[i]]$cov_results_het1$ma_model[[4]]$b,
                                     sum_cell_25[[i]]$cov_results_het1$ma_model[[5]]$b, sum_cell_25[[i]]$cov_results_het1$ma_model[[6]]$b) - true_ergm_results_het0$coef) / true_ergm_results_het0$coef) ** 2)
  
}

rmse_true_25_het1 <- matrix(nrow = 124, ncol = 6)
for (i in 1 : nrow(rmse_true_25_het1)){
  
  rmse_true_25_het1[i, ] <- sqrt(((c(sum_cell_25[[i]]$cov_results_het1$ma_model[[1]]$b, sum_cell_25[[i]]$cov_results_het1$ma_model[[2]]$b,
                                     sum_cell_25[[i]]$cov_results_het1$ma_model[[3]]$b, sum_cell_25[[i]]$cov_results_het1$ma_model[[4]]$b,
                                     sum_cell_25[[i]]$cov_results_het1$ma_model[[5]]$b, sum_cell_25[[i]]$cov_results_het1$ma_model[[6]]$b) - true_ergm_results_het1$coef) / true_ergm_results_het1$coef) ** 2)
  
}

rmse_struct_25_het0 <- matrix(nrow = 124, ncol = 4)
for (i in 1 : nrow(rmse_struct_25_het0)){
  
  rmse_struct_25_het0[i, ] <- sqrt(((c(sum_cell_25[[i]]$gwodegree_results_het1$ma_model[[1]]$b, sum_cell_25[[i]]$gwodegree_results_het1$ma_model[[2]]$b,
                                       sum_cell_25[[i]]$gwodegree_results_het1$ma_model[[3]]$b, sum_cell_25[[i]]$gwodegree_results_het1$ma_model[[4]]$b) - true_ergm_results_het0$coef[1:4]) / true_ergm_results_het0$coef[1:4]) ** 2)
  
}

rmse_struct_25_het1 <- matrix(nrow = 124, ncol = 4)
for (i in 1 : nrow(rmse_struct_25_het1)){
  
  rmse_struct_25_het1[i, ] <- sqrt(((c(sum_cell_25[[i]]$gwodegree_results_het1$ma_model[[1]]$b, sum_cell_25[[i]]$gwodegree_results_het1$ma_model[[2]]$b,
                                       sum_cell_25[[i]]$gwodegree_results_het1$ma_model[[3]]$b, sum_cell_25[[i]]$gwodegree_results_het1$ma_model[[4]]$b) - true_ergm_results_het1$coef[1:4]) / true_ergm_results_het1$coef[1:4]) ** 2)
  
}

rmse_dtriad_25_het0 <- matrix(nrow = 124, ncol = 6)
for (i in 1 : nrow(rmse_dtriad_25_het0)){
  
  rmse_dtriad_25_het0[i, ] <- sqrt(((c(sum_cell_25[[i]]$dtriad_results_het1$ma_model[[1]]$b, sum_cell_25[[i]]$dtriad_results_het1$ma_model[[2]]$b,
                                       sum_cell_25[[i]]$dtriad_results_het1$ma_model[[3]]$b, sum_cell_25[[i]]$dtriad_results_het1$ma_model[[4]]$b,
                                       sum_cell_25[[i]]$dtriad_results_het1$ma_model[[5]]$b, sum_cell_25[[i]]$dtriad_results_het1$ma_model[[6]]$b) - true_ergm_results_het0$coef) / true_ergm_results_het0$coef) ** 2)
  
}

rmse_dtriad_25_het1 <- matrix(nrow = 124, ncol = 6)
for (i in 1 : nrow(rmse_dtriad_25_het1)){
  
  rmse_dtriad_25_het1[i, ] <- sqrt(((c(sum_cell_25[[i]]$dtriad_results_het1$ma_model[[1]]$b, sum_cell_25[[i]]$dtriad_results_het1$ma_model[[2]]$b,
                                       sum_cell_25[[i]]$dtriad_results_het1$ma_model[[3]]$b, sum_cell_25[[i]]$dtriad_results_het1$ma_model[[4]]$b,
                                       sum_cell_25[[i]]$dtriad_results_het1$ma_model[[5]]$b, sum_cell_25[[i]]$dtriad_results_het1$ma_model[[6]]$b) - true_ergm_results_het1$coef) / true_ergm_results_het1$coef) ** 2)
  
}

rmse_true_75_het0 <- matrix(nrow = 124, ncol = 6)
for (i in 1 : nrow(rmse_true_75_het0)){
  
  rmse_true_75_het0[i, ] <- sqrt(((c(sum_cell_75[[i]]$cov_results_het1$ma_model[[1]]$b, sum_cell_75[[i]]$cov_results_het1$ma_model[[2]]$b,
                                     sum_cell_75[[i]]$cov_results_het1$ma_model[[3]]$b, sum_cell_75[[i]]$cov_results_het1$ma_model[[4]]$b,
                                     sum_cell_75[[i]]$cov_results_het1$ma_model[[5]]$b, sum_cell_75[[i]]$cov_results_het1$ma_model[[6]]$b) - true_ergm_results_het0$coef) / true_ergm_results_het0$coef)** 2)
  
}

rmse_true_75_het1 <- matrix(nrow = 124, ncol = 6)
for (i in 1 : nrow(rmse_true_75_het1)){
  
  rmse_true_75_het1[i, ] <- sqrt(((c(sum_cell_75[[i]]$cov_results_het1$ma_model[[1]]$b, sum_cell_75[[i]]$cov_results_het1$ma_model[[2]]$b,
                                     sum_cell_75[[i]]$cov_results_het1$ma_model[[3]]$b, sum_cell_75[[i]]$cov_results_het1$ma_model[[4]]$b,
                                     sum_cell_75[[i]]$cov_results_het1$ma_model[[5]]$b, sum_cell_75[[i]]$cov_results_het1$ma_model[[6]]$b) - true_ergm_results_het1$coef) / true_ergm_results_het1$coef) ** 2)
  
}

rmse_struct_75_het0 <- matrix(nrow = 124, ncol = 4)
for (i in 1 : nrow(rmse_struct_75_het0)){
  
  rmse_struct_75_het0[i, ] <- sqrt(((c(sum_cell_75[[i]]$gwodegree_results_het1$ma_model[[1]]$b, sum_cell_75[[i]]$gwodegree_results_het1$ma_model[[2]]$b,
                                       sum_cell_75[[i]]$gwodegree_results_het1$ma_model[[3]]$b, sum_cell_75[[i]]$gwodegree_results_het1$ma_model[[4]]$b) - true_ergm_results_het0$coef[1:4]) / true_ergm_results_het0$coef[1:4]) ** 2)
  
}

rmse_struct_75_het1 <- matrix(nrow = 124, ncol = 4)
for (i in 1 : nrow(rmse_struct_75_het1)){
  
  rmse_struct_75_het1[i, ] <- sqrt(((c(sum_cell_75[[i]]$gwodegree_results_het1$ma_model[[1]]$b, sum_cell_75[[i]]$gwodegree_results_het1$ma_model[[2]]$b,
                                       sum_cell_75[[i]]$gwodegree_results_het1$ma_model[[3]]$b, sum_cell_75[[i]]$gwodegree_results_het1$ma_model[[4]]$b) - true_ergm_results_het1$coef[1:4]) / true_ergm_results_het1$coef[1:4]) ** 2)
  
}

rmse_dtriad_75_het0 <- matrix(nrow = 124, ncol = 6)
for (i in 1 : nrow(rmse_dtriad_75_het0)){
  
  rmse_dtriad_75_het0[i, ] <- sqrt(((c(sum_cell_75[[i]]$dtriad_results_het1$ma_model[[1]]$b, sum_cell_75[[i]]$dtriad_results_het1$ma_model[[2]]$b,
                                       sum_cell_75[[i]]$dtriad_results_het1$ma_model[[3]]$b, sum_cell_75[[i]]$dtriad_results_het1$ma_model[[4]]$b,
                                       sum_cell_75[[i]]$dtriad_results_het1$ma_model[[5]]$b, sum_cell_75[[i]]$dtriad_results_het1$ma_model[[6]]$b) - true_ergm_results_het0$coef) / true_ergm_results_het0$coef) ** 2)
  
}

rmse_dtriad_75_het1 <- matrix(nrow = 124, ncol = 6)
for (i in 1 : nrow(rmse_dtriad_75_het1)){
  
  rmse_dtriad_75_het1[i, ] <- sqrt(((c(sum_cell_75[[i]]$dtriad_results_het1$ma_model[[1]]$b, sum_cell_75[[i]]$dtriad_results_het1$ma_model[[2]]$b,
                                       sum_cell_75[[i]]$dtriad_results_het1$ma_model[[3]]$b, sum_cell_75[[i]]$dtriad_results_het1$ma_model[[4]]$b,
                                       sum_cell_75[[i]]$dtriad_results_het1$ma_model[[5]]$b, sum_cell_75[[i]]$dtriad_results_het1$ma_model[[6]]$b) - true_ergm_results_het1$coef) / true_ergm_results_het1$coef) ** 2)
  
}

colnames(rmse_dtriad_25_hom) <- c("edges", "mutual", "gwesp", "gwodegree", "nodeocov", "nodematch")
colnames(rmse_dtriad_25_het0) <- c("edges", "mutual", "gwesp", "gwodegree", "nodeocov", "nodematch")
colnames(rmse_dtriad_25_het1) <- c("edges", "mutual", "gwesp", "gwodegree", "nodeocov", "nodematch")
colnames(rmse_dtriad_75_hom) <- c("edges", "mutual", "gwesp", "gwodegree", "nodeocov", "nodematch")
colnames(rmse_dtriad_75_het0) <- c("edges", "mutual", "gwesp", "gwodegree", "nodeocov", "nodematch")
colnames(rmse_dtriad_75_het1) <- c("edges", "mutual", "gwesp", "gwodegree", "nodeocov", "nodematch")

rmse_dtriad_25_hom <- melt(as.data.table(rmse_dtriad_25_hom))
rmse_dtriad_25_het0 <- melt(as.data.table(rmse_dtriad_25_het0))
rmse_dtriad_25_het1 <- melt(as.data.table(rmse_dtriad_25_het1))
rmse_dtriad_75_hom <- melt(as.data.table(rmse_dtriad_75_hom))
rmse_dtriad_75_het0 <- melt(as.data.table(rmse_dtriad_75_het0))
rmse_dtriad_75_het1 <- melt(as.data.table(rmse_dtriad_75_het1))

rmse_dtriad_25_hom <- cbind(sampsize = "25", het = "Homogeneous sample", comp = "Population model 1", rmse_dtriad_25_hom)
rmse_dtriad_25_het0 <- cbind(sampsize = "25", het = "Heterogeneous sample", comp = "Population model 1", rmse_dtriad_25_het0)
rmse_dtriad_25_het1 <- cbind(sampsize = "25", het = "Heterogeneous sample", comp = "Population model 2", rmse_dtriad_25_het1)
rmse_dtriad_75_hom <- cbind(sampsize = "75", het = "Homogeneous sample", comp = "Population model 1", rmse_dtriad_75_hom)
rmse_dtriad_75_het0 <- cbind(sampsize = "75", het = "Heterogeneous sample", comp = "Population model 1", rmse_dtriad_75_het0)
rmse_dtriad_75_het1 <- cbind(sampsize = "75", het = "Heterogeneous sample", comp = "Population model 2", rmse_dtriad_75_het1)

colnames(rmse_struct_25_hom) <- c("edges", "mutual", "gwesp", "gwodegree")
colnames(rmse_struct_25_het0) <- c("edges", "mutual", "gwesp", "gwodegree")
colnames(rmse_struct_25_het1) <- c("edges", "mutual", "gwesp", "gwodegree")
colnames(rmse_struct_75_hom) <- c("edges", "mutual", "gwesp", "gwodegree")
colnames(rmse_struct_75_het0) <- c("edges", "mutual", "gwesp", "gwodegree")
colnames(rmse_struct_75_het1) <- c("edges", "mutual", "gwesp", "gwodegree")

rmse_struct_25_hom <- melt(as.data.table(rmse_struct_25_hom))
rmse_struct_25_het0 <- melt(as.data.table(rmse_struct_25_het0))
rmse_struct_25_het1 <- melt(as.data.table(rmse_struct_25_het1))
rmse_struct_75_hom <- melt(as.data.table(rmse_struct_75_hom))
rmse_struct_75_het0 <- melt(as.data.table(rmse_struct_75_het0))
rmse_struct_75_het1 <- melt(as.data.table(rmse_struct_75_het1))

rmse_struct_25_hom <- cbind(sampsize = "25", het = "Homogeneous sample", comp = "Population model 1", rmse_struct_25_hom)
rmse_struct_25_het0 <- cbind(sampsize = "25", het = "Heterogeneous sample", comp = "Population model 1", rmse_struct_25_het0)
rmse_struct_25_het1 <- cbind(sampsize = "25", het = "Heterogeneous sample", comp = "Population model 2", rmse_struct_25_het1)
rmse_struct_75_hom <- cbind(sampsize = "75", het = "Homogeneous sample", comp = "Population model 1", rmse_struct_75_hom)
rmse_struct_75_het0 <- cbind(sampsize = "75", het = "Heterogeneous sample", comp = "Population model 1", rmse_struct_75_het0)
rmse_struct_75_het1 <- cbind(sampsize = "75", het = "Heterogeneous sample", comp = "Population model 2", rmse_struct_75_het1)

colnames(rmse_true_25_hom) <- c("edges", "mutual", "gwesp", "gwodegree", "nodeocov", "nodematch")
colnames(rmse_true_25_het0) <- c("edges", "mutual", "gwesp", "gwodegree", "nodeocov", "nodematch")
colnames(rmse_true_25_het1) <- c("edges", "mutual", "gwesp", "gwodegree", "nodeocov", "nodematch")
colnames(rmse_true_75_hom) <- c("edges", "mutual", "gwesp", "gwodegree", "nodeocov", "nodematch")
colnames(rmse_true_75_het0) <- c("edges", "mutual", "gwesp", "gwodegree", "nodeocov", "nodematch")
colnames(rmse_true_75_het1) <- c("edges", "mutual", "gwesp", "gwodegree", "nodeocov", "nodematch")

rmse_true_25_hom <- melt(as.data.table(rmse_true_25_hom))
rmse_true_25_het0 <- melt(as.data.table(rmse_true_25_het0))
rmse_true_25_het1 <- melt(as.data.table(rmse_true_25_het1))
rmse_true_75_hom <- melt(as.data.table(rmse_true_75_hom))
rmse_true_75_het0 <- melt(as.data.table(rmse_true_75_het0))
rmse_true_75_het1 <- melt(as.data.table(rmse_true_75_het1))

rmse_true_25_hom <- cbind(sampsize = '25', het = "Homogeneous sample", comp = "Population model 1", rmse_true_25_hom)
rmse_true_25_het0 <- cbind(sampsize = "25", het = "Heterogeneous sample", comp = "Population model 1", rmse_true_25_het0)
rmse_true_25_het1 <- cbind(sampsize = "25", het = "Heterogeneous sample", comp = "Population model 2", rmse_true_25_het1)
rmse_true_75_hom <- cbind(sampsize = "75", het = "Homogeneous sample", comp = "Population model 1", rmse_true_75_hom)
rmse_true_75_het0 <- cbind(sampsize = "75", het = "Heterogeneous sample", comp = "Population model 1", rmse_true_75_het0)
rmse_true_75_het1 <- cbind(sampsize = "75", het = "Heterogeneous sample", comp = "Population model 2", rmse_true_75_het1)

dat_rmse_dtriad <- rbind.fill(data.frame(rmse_dtriad_25_hom), data.frame(rmse_dtriad_25_het0), data.frame(rmse_dtriad_25_het1),
                              data.frame(rmse_dtriad_75_hom), data.frame(rmse_dtriad_75_het0), data.frame(rmse_dtriad_75_het1))

dat_rmse_structural <- rbind.fill(data.frame(rmse_struct_25_hom), data.frame(rmse_struct_25_het0), data.frame(rmse_struct_25_het1),
                                  data.frame(rmse_struct_75_hom), data.frame(rmse_struct_75_het0), data.frame(rmse_struct_75_het1))

dat_rmse_true <- rbind.fill(data.frame(rmse_true_25_hom), data.frame(rmse_true_25_het0), data.frame(rmse_true_25_het1),
                            data.frame(rmse_true_75_hom), data.frame(rmse_true_75_het0), data.frame(rmse_true_75_het1))
dat_rmse_dtriad <- dat_rmse_dtriad[dat_rmse_dtriad[, "value"] < 10, ] # None removed.
dat_rmse_structural <- dat_rmse_structural[dat_rmse_structural[, "value"] < 10, ] # 990, 1486, 2461 and 2957 removed.
dat_rmse_true <- dat_rmse_true[dat_rmse_true[, "value"] < 10, ] # None removed.

dat_rmse_dtriad$het <- factor(dat_rmse_dtriad$het, levels = c('Homogeneous sample','Heterogeneous sample'), ordered = TRUE)
dat_rmse_structural$het <- factor(dat_rmse_structural$het, levels = c('Homogeneous sample','Heterogeneous sample'), ordered = TRUE)
dat_rmse_true$het <- factor(dat_rmse_true$het, levels = c('Homogeneous sample','Heterogeneous sample'), ordered = TRUE)

ggplot(dat_rmse_dtriad, aes(x = variable, y = value * 100, fill = sampsize)) +
  geom_boxplot() +
  facet_nested(~ het + comp) +
  labs(x = "Parameter", y = "RMSE (in %)") +
  scale_fill_discrete(name = "Sample size") +
  scale_x_discrete(labels = c('Edges','Mutual','Gwesp', "Gwodegree", "Nodeocov", "Nodematch")) +
  theme(axis.text.x = element_text(angle = 40, vjust = 1, hjust = 1)) +
  ylim(0, 300)

ggplot(dat_rmse_structural, aes(x = variable, y = value * 100, fill = sampsize)) +
  geom_boxplot() +
  facet_nested(~ het + comp) +
  labs(x = "Parameter", y = "RMSE (in %)") +
  scale_fill_discrete(name = "Sample size") +
  scale_x_discrete(labels = c('Edges','Mutual','Gwesp', "Gwodegree")) +
  theme(axis.text.x = element_text(angle = 40, vjust = 1, hjust = 1)) +
  ylim(0, 300)

ggplot(dat_rmse_true, aes(x = variable, y = value * 100, fill = sampsize)) +
  geom_boxplot() +
  facet_nested(~ het + comp) +
  labs(x = "Parameter", y = "RMSE (in %)") +
  scale_fill_discrete(name = "Sample size") +
  scale_x_discrete(labels = c('Edges','Mutual','Gwesp', "Gwodegree", "Nodeocov", "Nodematch")) +
  theme(axis.text.x = element_text(angle = 40, vjust = 1, hjust = 1)) +
  ylim(0, 300)

### AIC and SIC meta ERGM.
diff_true_struct_25_het0 <- vector()
diff_true_dtriad_25_het0 <- vector()
for (i in 1:124){
  diff_true_struct_25_het0[i] <- sum_cell_25[[i]]$cov_results_het0$aic$mean_aic - sum_cell_25[[i]]$gwodegree_results_het0$aic$mean_aic
  diff_true_dtriad_25_het0[i] <- sum_cell_25[[i]]$cov_results_het0$aic$mean_aic - sum_cell_25[[i]]$dtriad_results_het0$aic$mean_aic
}

diff_true_struct_25_het1 <- vector()
diff_true_dtriad_25_het1 <- vector()
for (i in 1:124){
  diff_true_struct_25_het1[i] <- sum_cell_25[[i]]$cov_results_het1$aic$mean_aic - sum_cell_25[[i]]$gwodegree_results_het1$aic$mean_aic
  diff_true_dtriad_25_het1[i] <- sum_cell_25[[i]]$cov_results_het1$aic$mean_aic - sum_cell_25[[i]]$dtriad_results_het1$aic$mean_aic
}

diff_true_struct_75_het0 <- vector()
diff_true_dtriad_75_het0 <- vector()
for (i in 1:124){
  diff_true_struct_75_het0[i] <- sum_cell_75[[i]]$cov_results_het0$aic$mean_aic - sum_cell_75[[i]]$gwodegree_results_het0$aic$mean_aic
  diff_true_dtriad_75_het0[i] <- sum_cell_75[[i]]$cov_results_het0$aic$mean_aic - sum_cell_75[[i]]$dtriad_results_het0$aic$mean_aic
  
}

diff_true_struct_75_het1 <- vector()
diff_true_dtriad_75_het1 <- vector()
for (i in 1:124){
  diff_true_struct_75_het1[i] <- sum_cell_75[[i]]$cov_results_het1$aic$mean_aic - sum_cell_75[[i]]$gwodegree_results_het1$aic$mean_aic
  diff_true_dtriad_75_het1[i] <- sum_cell_75[[i]]$cov_results_het1$aic$mean_aic - sum_cell_75[[i]]$dtriad_results_het1$aic$mean_aic
}

diff_true_struct_25_het0 <- cbind(sampsize = "25", het = "Homogeneous sample", comp = "true_struct", diff_true_struct_25_het0)
diff_true_dtriad_25_het0 <- cbind(sampsize = "25", het = "Homogeneous sample", comp = "true_dtriad", diff_true_dtriad_25_het0)
diff_true_struct_25_het1 <- cbind(sampsize = "25", het = "Heterogeneous sample", comp = "true_struct", diff_true_struct_25_het1)
diff_true_dtriad_25_het1 <- cbind(sampsize = "25", het = "Heterogeneous sample", comp = "true_dtriad", diff_true_dtriad_25_het1)
diff_true_struct_75_het0 <- cbind(sampsize = "75", het = "Homogeneous sample", comp = "true_struct", diff_true_struct_75_het0)
diff_true_dtriad_75_het0 <- cbind(sampsize = "75", het = "Homogeneous sample", comp = "true_dtriad", diff_true_dtriad_75_het0)
diff_true_struct_75_het1 <- cbind(sampsize = "75", het = "Heterogeneous sample", comp = "true_struct", diff_true_struct_75_het1)
diff_true_dtriad_75_het1 <- cbind(sampsize = "75", het = "Heterogeneous sample", comp = "true_dtriad", diff_true_dtriad_75_het1)

colnames(diff_true_struct_25_het0)[4] <- "value"
colnames(diff_true_dtriad_25_het0)[4] <- "value"
colnames(diff_true_struct_25_het1)[4] <- "value"
colnames(diff_true_dtriad_25_het1)[4] <- "value"
colnames(diff_true_struct_75_het0)[4] <- "value"
colnames(diff_true_dtriad_75_het0)[4] <- "value"
colnames(diff_true_struct_75_het1)[4] <- "value"
colnames(diff_true_dtriad_75_het1)[4] <- "value"

dat_aic <- rbind(data.frame(diff_true_struct_25_het0), data.frame(diff_true_dtriad_25_het0), data.frame(diff_true_struct_25_het1),
                 data.frame(diff_true_dtriad_25_het1), data.frame(diff_true_struct_75_het0), data.frame(diff_true_dtriad_75_het0),
                 data.frame(diff_true_struct_75_het1), data.frame(diff_true_dtriad_75_het1))
dat_aic[ , 4] <- round(as.numeric(dat_aic[, 4]), 4)
dat_aic <- dat_aic[dat_aic[, 4] > - 100, ] # 120 252 260 296 321 322 370 443 485 635 712 719 737 788 849 857 865 867 899 removed.
dat_aic <- dat_aic[dat_aic[, 4] < 100, ] # 94 218 365 481 489 787 911 removed.
dat_aic <- dat_aic[complete.cases(dat_aic), ]

dat_aic$het <- factor(dat_aic$het, levels = c('Homogeneous sample','Heterogeneous sample'),ordered = TRUE)
dat_aic$comp <- factor(dat_aic$comp, levels = c('true_struct','true_dtriad'),ordered = TRUE)

ggplot(dat_aic, aes(x = comp, y = value, fill = sampsize)) +
  geom_boxplot() +
  facet_wrap(~ het) +
  labs(x = "Model relative to true model", y = "AIC Difference") +
  scale_fill_discrete(name = "Sample size") +
  scale_x_discrete(labels = c("Structural", "Dense triadic"))

diff_true_struct_25_het0 <- vector()
diff_true_dtriad_25_het0 <- vector()
for (i in 1:124){
  diff_true_struct_25_het0[i] <- sum_cell_25[[i]]$cov_results_het0$bic$mean_bic - sum_cell_25[[i]]$gwodegree_results_het0$bic$mean_bic
  diff_true_dtriad_25_het0[i] <- sum_cell_25[[i]]$cov_results_het0$bic$mean_bic - sum_cell_25[[i]]$dtriad_results_het0$bic$mean_bic
}

diff_true_struct_25_het1 <- vector()
diff_true_dtriad_25_het1 <- vector()
for (i in 1:124){
  diff_true_struct_25_het1[i] <- sum_cell_25[[i]]$cov_results_het1$bic$mean_bic - sum_cell_25[[i]]$gwodegree_results_het1$bic$mean_bic
  diff_true_dtriad_25_het1[i] <- sum_cell_25[[i]]$cov_results_het1$bic$mean_bic - sum_cell_25[[i]]$dtriad_results_het1$bic$mean_bic
}

diff_true_struct_75_het0 <- vector()
diff_true_dtriad_75_het0 <- vector()
for (i in 1:124){
  diff_true_struct_75_het0[i] <- sum_cell_75[[i]]$cov_results_het0$bic$mean_bic - sum_cell_75[[i]]$gwodegree_results_het0$bic$mean_bic
  diff_true_dtriad_75_het0[i] <- sum_cell_75[[i]]$cov_results_het0$bic$mean_bic - sum_cell_75[[i]]$dtriad_results_het0$bic$mean_bic
  
}

diff_true_struct_75_het1 <- vector()
diff_true_dtriad_75_het1 <- vector()
for (i in 1:124){
  diff_true_struct_75_het1[i] <- sum_cell_75[[i]]$cov_results_het1$bic$mean_bic - sum_cell_75[[i]]$gwodegree_results_het1$bic$mean_bic
  diff_true_dtriad_75_het1[i] <- sum_cell_75[[i]]$cov_results_het1$bic$mean_bic - sum_cell_75[[i]]$dtriad_results_het1$bic$mean_bic
}

diff_true_struct_25_het0 <- cbind(sampsize = "25", het = "Homogeneous sample", comp = "true_struct", diff_true_struct_25_het0)
diff_true_dtriad_25_het0 <- cbind(sampsize = "25", het = "Homogeneous sample", comp = "true_dtriad", diff_true_dtriad_25_het0)
diff_true_struct_25_het1 <- cbind(sampsize = "25", het = "Heterogeneous sample", comp = "true_struct", diff_true_struct_25_het1)
diff_true_dtriad_25_het1 <- cbind(sampsize = "25", het = "Heterogeneous sample", comp = "true_dtriad", diff_true_dtriad_25_het1)
diff_true_struct_75_het0 <- cbind(sampsize = "75", het = "Homogeneous sample", comp = "true_struct", diff_true_struct_75_het0)
diff_true_dtriad_75_het0 <- cbind(sampsize = "75", het = "Homogeneous sample", comp = "true_dtriad", diff_true_dtriad_75_het0)
diff_true_struct_75_het1 <- cbind(sampsize = "75", het = "Heterogeneous sample", comp = "true_struct", diff_true_struct_75_het1)
diff_true_dtriad_75_het1 <- cbind(sampsize = "75", het = "Heterogeneous sample", comp = "true_dtriad", diff_true_dtriad_75_het1)

colnames(diff_true_struct_25_het0)[4] <- "value"
colnames(diff_true_dtriad_25_het0)[4] <- "value"
colnames(diff_true_struct_25_het1)[4] <- "value"
colnames(diff_true_dtriad_25_het1)[4] <- "value"
colnames(diff_true_struct_75_het0)[4] <- "value"
colnames(diff_true_dtriad_75_het0)[4] <- "value"
colnames(diff_true_struct_75_het1)[4] <- "value"
colnames(diff_true_dtriad_75_het1)[4] <- "value"

dat_bic <- rbind(data.frame(diff_true_struct_25_het0), data.frame(diff_true_dtriad_25_het0), data.frame(diff_true_struct_25_het1),
                 data.frame(diff_true_dtriad_25_het1), data.frame(diff_true_struct_75_het0), data.frame(diff_true_dtriad_75_het0),
                 data.frame(diff_true_struct_75_het1), data.frame(diff_true_dtriad_75_het1))
dat_bic[ , 4] <- round(as.numeric(dat_bic[, 4]), 4)
dat_bic <- dat_bic[dat_bic[, 4] > - 100, ] # 120 252 260 296 321 322 370 443 485 635 712 719 737 788 849 857 865 867 899 removed.
dat_bic <- dat_bic[dat_bic[, 4] < 100, ] # 94 218 365 481 489 787 911 removed.
dat_bic <- dat_bic[complete.cases(dat_bic), ]

dat_bic$het <- factor(dat_bic$het, levels = c('Homogeneous sample','Heterogeneous sample'),ordered = TRUE)
dat_bic$comp <- factor(dat_bic$comp, levels = c('true_struct','true_dtriad'),ordered = TRUE)

ggplot(dat_bic, aes(x = comp, y = value, fill = sampsize)) +
  geom_boxplot() +
  facet_wrap(~ het) +
  labs(x = "Model relative to true model", y = "SIC Difference") +
  scale_fill_discrete(name = "Sample size") +
  scale_x_discrete(labels = c("Structural", "Dense triadic"))

###################################
##### Auxilary statistics GOF #####
###################################

##########################################
### Condition: 75 - structural - het0. ###
##########################################

coefs <- c(sum_cell_75[[1]]$gwodegree_results_het0$ma_model[[1]]$b, sum_cell_75[[1]]$gwodegree_results_het0$ma_model[[2]]$b,
           sum_cell_75[[1]]$gwodegree_results_het0$ma_model[[3]]$b, sum_cell_75[[1]]$gwodegree_results_het0$ma_model[[4]]$b)
simnet <- simulate(network(25) ~ edges + mutual + gwesp(0.25, fixed = T) + gwodegree(0.10, fixed = T), coef = coefs)
meta_ergm_fits <- gof(object = simnet ~  edges + mutual + gwesp(0.25, fixed = T) + gwodegree(0.10, fixed = T),
                      GOF =~ model + idegree + odegree + distance + espartners + dspartners + triadcensus,
                      coef = coefs, control.gof.formula(seed = 3791, nsim = 1000))

### In-degree.
idegree <- matrix(nrow = 63, ncol = 25)
for (i in 1 : 63){
  idegree[i, ] <- cell_75[[1]]$gwodegree_results_het0$ergm_fits[[i]]$obs.ideg
}

dat1 <- melt(data.table(meta_ergm_fits$sim.ideg))
dat2 <- melt(data.table(idegree))
dat2[ , 1] <- as.numeric(unlist(dat2[ , 1]))
dat1 <- as.data.frame(dat1)
dat2 <- as.data.frame(dat2)

struct_het0_indegree <- ggplot(dat1, aes(x = variable, y = value)) +
  geom_boxplot() +
  geom_jitter(data = dat2, aes(x = variable, y = value), position = position_jitter(0.1), color = "blue", shape = 18, alpha = 0.3, inherit.aes = FALSE) +
  labs(x = "In-degree", y = "Frequency") +
  theme(axis.text.x = element_text(angle = 40, vjust = 1, hjust=1))

### Out-degree.
odegree <- matrix(nrow = 63, ncol = 25)
for (i in 1 : 63){
  odegree[i, ] <- cell_75[[1]]$gwodegree_results_het0$ergm_fits[[i]]$obs.odeg
}
dat1 <- melt(data.table(meta_ergm_fits$sim.odeg))
dat2 <- melt(data.table(odegree))
dat2[ , 1] <- as.numeric(unlist(dat2[ , 1]))
dat1 <- as.data.frame(dat1)
dat2 <- as.data.frame(dat2)

struct_het0_outdegree <- ggplot(dat1, aes(x = variable, y = value)) +
  geom_boxplot() +
  geom_jitter(data = dat2, aes(x = variable, y = value), position = position_jitter(0.1), color = "blue", shape = 18, alpha = 0.3, inherit.aes = FALSE) +
  labs(x = "Out-degree", y = "Frequency") +
  theme(axis.text.x = element_text(angle = 40, vjust = 1, hjust=1))+
  ylim(-1, 16)

### Geodesic distance.
dist <- matrix(nrow = 63, ncol = 25)
for (i in 1 : 63){
  dist[i, ] <- cell_75[[1]]$gwodegree_results_het0$ergm_fits[[i]]$obs.dist
}
dat1 <- melt(data.table(meta_ergm_fits$sim.dist))
dat2 <- melt(data.table(dist))
dat2[ , 1] <- as.numeric(unlist(dat2[ , 1]))
dat1 <- as.data.frame(dat1)
dat2 <- as.data.frame(dat2)

struct_het0_dist <- ggplot(dat1, aes(x = variable, y = value)) +
  geom_boxplot() +
  geom_jitter(data = dat2, aes(x = variable, y = value), position = position_jitter(0.1), color = "blue", shape = 18, alpha = 0.3, inherit.aes = FALSE) +
  labs(x = "Geodesic distance", y = "Frequency") +
  theme(axis.text.x = element_text(angle = 40, vjust = 1, hjust=1)) +
  ylim(-5, 450)

### Triad census.
triad_census <- matrix(nrow = 63, ncol = 16)
for (i in 1 : 63){
  triad_census[i, ] <- cell_75[[1]]$gwodegree_results_het0$ergm_fits[[i]]$obs.triadcensus
}
dat1 <- melt(data.table(meta_ergm_fits$sim.triadcensus))
dat2 <- melt(data.table(triad_census))
dat2[ , 1] <- as.numeric(unlist(dat2[ , 1]))
dat1 <- as.data.frame(dat1)
dat2 <- as.data.frame(dat2)

struct_het0_triadc <- ggplot(dat1, aes(x = variable, y = value)) +
  geom_boxplot() +
  geom_jitter(data = dat2, aes(x = variable, y = value), position = position_jitter(0.1), color = "blue", shape = 18, alpha = 0.3, inherit.aes = FALSE) +
  labs(x = "Triad-census", y = "Frequency") +
  theme(axis.text.x = element_text(angle = 40, vjust = 1, hjust=1)) +
  ylim(-5, 1300)

##########################################
### Condition: 75 - structural - het1. ###
##########################################

coefs <- c(sum_cell_75[[1]]$gwodegree_results_het1$ma_model[[1]]$b, sum_cell_75[[1]]$gwodegree_results_het1$ma_model[[2]]$b,
           sum_cell_75[[1]]$gwodegree_results_het1$ma_model[[3]]$b, sum_cell_75[[1]]$gwodegree_results_het1$ma_model[[4]]$b)
simnet <- simulate(network(25) ~ edges + mutual + gwesp(0.25, fixed = T) + gwodegree(0.10, fixed = T), coef = coefs)
meta_ergm_fits <- gof(object = simnet ~  edges + mutual + gwesp(0.25, fixed = T) + gwodegree(0.10, fixed = T),
                      GOF =~ model + idegree + odegree + distance + espartners + dspartners + triadcensus,
                      coef = coefs, control.gof.formula(seed = 3791, nsim = 1000))

### In-degree.
idegree <- matrix(nrow = 57, ncol = 25)
for (i in 1 : 57){
  idegree[i, ] <- cell_75[[1]]$gwodegree_results_het1$ergm_fits[[i]]$obs.ideg
}
dat1 <- melt(data.table(meta_ergm_fits$sim.ideg))
dat2 <- melt(data.table(idegree))
dat2[ , 1] <- as.numeric(unlist(dat2[ , 1]))
dat1 <- as.data.frame(dat1)
dat2 <- as.data.frame(dat2)

struct_het1_indegree <- ggplot(dat1, aes(x = variable, y = value)) +
  geom_boxplot() +
  geom_jitter(data = dat2, aes(x = variable, y = value), position = position_jitter(0.1), color = "blue", shape = 18, alpha = 0.3, inherit.aes = FALSE) +
  labs(x = "In-degree", y = "Frequency") +
  theme(axis.text.x = element_text(angle = 40, vjust = 1, hjust=1))+
  ylim(-1, 13)

### Out-degree.
odegree <- matrix(nrow = 57, ncol = 25)
for (i in 1 : 57){
  odegree[i, ] <- cell_75[[1]]$gwodegree_results_het1$ergm_fits[[i]]$obs.odeg
}
dat1 <- melt(data.table(meta_ergm_fits$sim.odeg))
dat2 <- melt(data.table(odegree))
dat2[ , 1] <- as.numeric(unlist(dat2[ , 1]))
dat1 <- as.data.frame(dat1)
dat2 <- as.data.frame(dat2)

struct_het1_outdegree <- ggplot(dat1, aes(x = variable, y = value)) +
  geom_boxplot() +
  geom_jitter(data = dat2, aes(x = variable, y = value), position = position_jitter(0.1), color = "blue", shape = 18, alpha = 0.3, inherit.aes = FALSE) +
  labs(x = "Out-degree", y = "Frequency") +
  theme(axis.text.x = element_text(angle = 40, vjust = 1, hjust=1))+
  ylim(-1, 16)

### Geodesic distance.
dist <- matrix(nrow = 57, ncol = 25)
for (i in 1 : 57){
  dist[i, ] <- cell_75[[1]]$gwodegree_results_het1$ergm_fits[[i]]$obs.dist
}
dat1 <- melt(data.table(meta_ergm_fits$sim.dist))
dat2 <- melt(data.table(dist))
dat2[ , 1] <- as.numeric(unlist(dat2[ , 1]))
dat1 <- as.data.frame(dat1)
dat2 <- as.data.frame(dat2)

struct_het1_dist <- ggplot(dat1, aes(x = variable, y = value)) +
  geom_boxplot() +
  geom_jitter(data = dat2, aes(x = variable, y = value), position = position_jitter(0.1), color = "blue", shape = 18, alpha = 0.3, inherit.aes = FALSE) +
  labs(x = "Geodesic distance", y = "Frequency") +
  theme(axis.text.x = element_text(angle = 40, vjust = 1, hjust=1))+
  ylim(-5, 500)

### Triad census.
triad_census <- matrix(nrow = 57, ncol = 16)
for (i in 1 : 57){
  triad_census[i, ] <- cell_75[[1]]$gwodegree_results_het1$ergm_fits[[i]]$obs.triadcensus
}
dat1 <- melt(data.table(meta_ergm_fits$sim.triadcensus))
dat2 <- melt(data.table(triad_census))
dat2[ , 1] <- as.numeric(unlist(dat2[ , 1]))
dat1 <- as.data.frame(dat1)
dat2 <- as.data.frame(dat2)

struct_het1_triadc <- ggplot(dat1, aes(x = variable, y = value)) +
  geom_boxplot() +
  geom_jitter(data = dat2, aes(x = variable, y = value), position = position_jitter(0.1), color = "blue", shape = 18, alpha = 0.3, inherit.aes = FALSE) +
  labs(x = "Triad-census", y = "Frequency") +
  theme(axis.text.x = element_text(angle = 40, vjust = 1, hjust=1)) +
  ylim(-5, 1500)

##########################################
### Condition: 75 - true - het0. ###
##########################################

net <- network(25)
net %v% "sex" <- cov_het0$girl # Add sex actor variable.
net %v% "identclass" <- cov_het0$identcGMC # Add class membership actor variable.

coefs <- c(sum_cell_75[[1]]$cov_results_het0$ma_model[[1]]$b, sum_cell_75[[1]]$cov_results_het0$ma_model[[2]]$b,
           sum_cell_75[[1]]$cov_results_het0$ma_model[[3]]$b, sum_cell_75[[1]]$cov_results_het0$ma_model[[4]]$b,
           sum_cell_75[[1]]$cov_results_het0$ma_model[[5]]$b, sum_cell_75[[1]]$cov_results_het0$ma_model[[6]]$b)
simnet <- simulate(net ~ edges + mutual + gwesp(0.25, fixed = T) + gwodegree(0.10, fixed = T) + nodeocov("identclass") + nodematch("sex"), coef = coefs)
meta_ergm_fits <- gof(object = simnet ~  edges + mutual + gwesp(0.25, fixed = T) + gwodegree(0.10, fixed = T) + nodeocov("identclass") + nodematch("sex"),
                      GOF =~ model + idegree + odegree + distance + espartners + dspartners + triadcensus,
                      coef = coefs, control.gof.formula(seed = 3791, nsim = 1000))

### In-degree.
idegree <- matrix(nrow = 75, ncol = 25)
for (i in 1 : 75){
  idegree[i, ] <- cell_75[[1]]$cov_results_het0$ergm_fits[[i]]$obs.ideg
}
dat1 <- melt(data.table(meta_ergm_fits$sim.ideg))
dat2 <- melt(data.table(idegree))
dat2[ , 1] <- as.numeric(unlist(dat2[ , 1]))
dat1 <- as.data.frame(dat1)
dat2 <- as.data.frame(dat2)

true_het0_indegree <- ggplot(dat1, aes(x = variable, y = value)) +
  geom_boxplot() +
  geom_jitter(data = dat2, aes(x = variable, y = value), position = position_jitter(0.1), color = "blue", shape = 18, alpha = 0.3, inherit.aes = FALSE) +
  labs(x = "In-degree", y = "Frequency") +
  theme(axis.text.x = element_text(angle = 40, vjust = 1, hjust = 1))+
  ylim(-1, 17)

### Out-degree.
odegree <- matrix(nrow = 75, ncol = 25)
for (i in 1 : 75){
  odegree[i, ] <- cell_75[[1]]$cov_results_het0$ergm_fits[[i]]$obs.odeg
}
dat1 <- melt(data.table(meta_ergm_fits$sim.odeg))
dat2 <- melt(data.table(odegree))
dat2[ , 1] <- as.numeric(unlist(dat2[ , 1]))
dat1 <- as.data.frame(dat1)
dat2 <- as.data.frame(dat2)

true_het0_outdegree <- ggplot(dat1, aes(x = variable, y = value)) +
  geom_boxplot() +
  geom_jitter(data = dat2, aes(x = variable, y = value), position = position_jitter(0.1), color = "blue", shape = 18, alpha = 0.3, inherit.aes = FALSE) +
  labs(x = "Out-degree", y = "Frequency") +
  theme(axis.text.x = element_text(angle = 40, vjust = 1, hjust = 1))+
  ylim(-1, 16)

### Geodesic distance.
dist <- matrix(nrow = 75, ncol = 25)
for (i in 1 : 75){
  dist[i, ] <- cell_75[[1]]$cov_results_het0$ergm_fits[[i]]$obs.dist
}
dat1 <- melt(data.table(meta_ergm_fits$sim.dist))
dat2 <- melt(data.table(dist))
dat2[ , 1] <- as.numeric(unlist(dat2[ , 1]))
dat1 <- as.data.frame(dat1)
dat2 <- as.data.frame(dat2)

true_het0_dist <- ggplot(dat1, aes(x = variable, y = value)) +
  geom_boxplot() +
  geom_jitter(data = dat2, aes(x = variable, y = value), position = position_jitter(0.1), color = "blue", shape = 18, alpha = 0.3, inherit.aes = FALSE) +
  labs(x = "Geodesic distance", y = "Frequency") +
  theme(axis.text.x = element_text(angle = 40, vjust = 1, hjust = 1))+
  ylim(-5, 450)

### Triad census.
triad_census <- matrix(nrow = 75, ncol = 16)
for (i in 1 : 75){
  triad_census[i, ] <- cell_75[[1]]$cov_results_het0$ergm_fits[[i]]$obs.triadcensus
}
dat1 <- melt(data.table(meta_ergm_fits$sim.triadcensus))
dat2 <- melt(data.table(triad_census))
dat2[ , 1] <- as.numeric(unlist(dat2[ , 1]))
dat1 <- as.data.frame(dat1)
dat2 <- as.data.frame(dat2)

true_het0_triadc <- ggplot(dat1, aes(x = variable, y = value)) +
  geom_boxplot() +
  geom_jitter(data = dat2, aes(x = variable, y = value), position = position_jitter(0.1), color = "blue", shape = 18, alpha = 0.3, inherit.aes = FALSE) +
  labs(x = "Triad census", y = "Frequency") +
  theme(axis.text.x = element_text(angle = 40, vjust = 1, hjust = 1)) +
  ylim(-5, 1300)

##########################################
### Condition: 75 - true - het1. ###
##########################################

net <- network(25)
net %v% "sex" <- cov_het1$girl # Add sex actor variable.
net %v% "identclass" <- cov_het1$identcGMC # Add class membership actor variable.

coefs <- c(sum_cell_75[[1]]$cov_results_het1$ma_model[[1]]$b, sum_cell_75[[1]]$cov_results_het1$ma_model[[2]]$b,
           sum_cell_75[[1]]$cov_results_het1$ma_model[[3]]$b, sum_cell_75[[1]]$cov_results_het1$ma_model[[4]]$b,
           sum_cell_75[[1]]$cov_results_het1$ma_model[[5]]$b, sum_cell_75[[1]]$cov_results_het1$ma_model[[6]]$b)
simnet <- simulate(net ~ edges + mutual + gwesp(0.25, fixed = T) + gwodegree(0.10, fixed = T) + nodeocov("identclass") + nodematch("sex"), coef = coefs)
meta_ergm_fits <- gof(object = simnet ~  edges + mutual + gwesp(0.25, fixed = T) + gwodegree(0.10, fixed = T) + nodeocov("identclass") + nodematch("sex"),
                      GOF =~ model + idegree + odegree + distance + espartners + dspartners + triadcensus,
                      coef = coefs, control.gof.formula(seed = 3791, nsim = 1000))

### In-degree.
idegree <- matrix(nrow = 75, ncol = 25)
for (i in 1 : 75){
  idegree[i, ] <- cell_75[[1]]$cov_results_het1$ergm_fits[[i]]$obs.ideg
}
dat1 <- melt(data.table(meta_ergm_fits$sim.ideg))
dat2 <- melt(data.table(idegree))
dat2[ , 1] <- as.numeric(unlist(dat2[ , 1]))
dat1 <- as.data.frame(dat1)
dat2 <- as.data.frame(dat2)

true_het1_indegree <- ggplot(dat1, aes(x = variable, y = value)) +
  geom_boxplot() +
  geom_jitter(data = dat2, aes(x = variable, y = value), position = position_jitter(0.1), color = "blue", shape = 18, alpha = 0.3, inherit.aes = FALSE) +
  labs(x = "In-degree", y = "Frequency") +
  theme(axis.text.x = element_text(angle = 40, vjust = 1, hjust = 1))+
  ylim(-1, 13)

### Out-degree.
odegree <- matrix(nrow = 75, ncol = 25)
for (i in 1 : 75){
  odegree[i, ] <- cell_75[[1]]$cov_results_het1$ergm_fits[[i]]$obs.odeg
}
dat1 <- melt(data.table(meta_ergm_fits$sim.odeg))
dat2 <- melt(data.table(odegree))
dat2[ , 1] <- as.numeric(unlist(dat2[ , 1]))
dat1 <- as.data.frame(dat1)
dat2 <- as.data.frame(dat2)

true_het1_outdegree <- ggplot(dat1, aes(x = variable, y = value)) +
  geom_boxplot() +
  geom_jitter(data = dat2, aes(x = variable, y = value), position = position_jitter(0.1), color = "blue", shape = 18, alpha = 0.3, inherit.aes = FALSE) +
  labs(x = "Out-degree", y = "Frequency") +
  theme(axis.text.x = element_text(angle = 40, vjust = 1, hjust = 1))+
  ylim(-1, 16)

### Geodesic distance.
dist <- matrix(nrow = 75, ncol = 25)
for (i in 1 : 75){
  dist[i, ] <- cell_75[[1]]$cov_results_het1$ergm_fits[[i]]$obs.dist
}
dat1 <- melt(data.table(meta_ergm_fits$sim.dist))
dat2 <- melt(data.table(dist))
dat2[ , 1] <- as.numeric(unlist(dat2[ , 1]))
dat1 <- as.data.frame(dat1)
dat2 <- as.data.frame(dat2)

true_het1_dist <- ggplot(dat1, aes(x = variable, y = value)) +
  geom_boxplot() +
  geom_jitter(data = dat2, aes(x = variable, y = value), position = position_jitter(0.1), color = "blue", shape = 18, alpha = 0.3, inherit.aes = FALSE) +
  labs(x = "Geodesic distance", y = "Frequency") +
  theme(axis.text.x = element_text(angle = 40, vjust = 1, hjust = 1))+
  ylim(-5, 500)

### Triad census.
triad_census <- matrix(nrow = 75, ncol = 16)
for (i in 1 : 75){
  triad_census[i, ] <- cell_75[[1]]$cov_results_het1$ergm_fits[[i]]$obs.triadcensus
}
dat1 <- melt(data.table(meta_ergm_fits$sim.triadcensus))
dat2 <- melt(data.table(triad_census))
dat2[ , 1] <- as.numeric(unlist(dat2[ , 1]))
dat1 <- as.data.frame(dat1)
dat2 <- as.data.frame(dat2)

true_het1_triadc <- ggplot(dat1, aes(x = variable, y = value)) +
  geom_boxplot() +
  geom_jitter(data = dat2, aes(x = variable, y = value), position = position_jitter(0.1), color = "blue", shape = 18, alpha = 0.3, inherit.aes = FALSE) +
  labs(x = "Triad-census", y = "Frequency") +
  theme(axis.text.x = element_text(angle = 40, vjust = 1, hjust = 1)) +
  ylim(-5, 1500)

##########################################
##### Condition: 75 - dtriad - het0. #####
##########################################

coefs <- c(sum_cell_75[[1]]$dtriad_results_het0$ma_model[[1]]$b, sum_cell_75[[1]]$dtriad_results_het0$ma_model[[2]]$b,
           sum_cell_75[[1]]$dtriad_results_het0$ma_model[[3]]$b, sum_cell_75[[1]]$dtriad_results_het0$ma_model[[4]]$b,
           sum_cell_75[[1]]$dtriad_results_het0$ma_model[[5]]$b, sum_cell_75[[1]]$dtriad_results_het0$ma_model[[6]]$b)

simnet <- simulate(net ~ edges + mutual + gwesp(0.25, fixed = T) + gwodegree(0.10, fixed = T) + nodeocov("identclass") + nodematch("sex"), coef = coefs)
meta_ergm_fits <- gof(object = simnet ~  edges + mutual + gwesp(0.25, fixed = T) + gwodegree(0.10, fixed = T) + nodeocov("identclass") + nodematch("sex"),
                      GOF =~ model + idegree + odegree + distance + espartners + dspartners + triadcensus,
                      coef = coefs, control.gof.formula(seed = 3791, nsim = 1000))

### In-degree.
idegree <- matrix(nrow = 38, ncol = 25)
for (i in 1 : 38){
  idegree[i, ] <- cell_75[[1]]$dtriad_results_het0$ergm_fits[[i]]$obs.ideg
}
dat1 <- melt(data.table(meta_ergm_fits$sim.ideg))
dat2 <- melt(data.table(idegree))
dat2[ , 1] <- as.numeric(unlist(dat2[ , 1]))
dat1 <- as.data.frame(dat1)
dat2 <- as.data.frame(dat2)

dtriad_het0_indegree <- ggplot(dat1, aes(x = variable, y = value)) +
  geom_boxplot() +
  geom_jitter(data = dat2, aes(x = variable, y = value), position = position_jitter(0.1), color = "blue", shape = 18, alpha = 0.3, inherit.aes = FALSE) +
  labs(x = "In-Degree", y = "Frequency") +
  theme(axis.text.x = element_text(angle = 40, vjust = 1, hjust = 1))+
  ylim(-1, 17)

### Out-degree.
odegree <- matrix(nrow = 38, ncol = 25)
for (i in 1 : 38){
  odegree[i, ] <- cell_75[[1]]$dtriad_results_het0$ergm_fits[[i]]$obs.odeg
}
dat1 <- melt(data.table(meta_ergm_fits$sim.odeg))
dat2 <- melt(data.table(odegree))
dat2[ , 1] <- as.numeric(unlist(dat2[ , 1]))
dat1 <- as.data.frame(dat1)
dat2 <- as.data.frame(dat2)

dtriad_het0_outdegree <- ggplot(dat1, aes(x = variable, y = value)) +
  geom_boxplot() +
  geom_jitter(data = dat2, aes(x = variable, y = value), position = position_jitter(0.1), color = "blue", shape = 18, alpha = 0.3, inherit.aes = FALSE) +
  labs(x = "Out-degree", y = "Frequency") +
  theme(axis.text.x = element_text(angle = 40, vjust = 1, hjust = 1)) +
  ylim(-1, 16)

### Geodesic distance.
dist <- matrix(nrow = 38, ncol = 25)
for (i in 1 : 38){
  dist[i, ] <- cell_75[[1]]$dtriad_results_het0$ergm_fits[[i]]$obs.dist
}
dat1 <- melt(data.table(meta_ergm_fits$sim.dist))
dat2 <- melt(data.table(dist))
dat2[ , 1] <- as.numeric(unlist(dat2[ , 1]))
dat1 <- as.data.frame(dat1)
dat2 <- as.data.frame(dat2)

dtriad_het0_dist <- ggplot(dat1, aes(x = variable, y = value)) +
  geom_boxplot() +
  geom_jitter(data = dat2, aes(x = variable, y = value), position = position_jitter(0.1), color = "blue", shape = 18, alpha = 0.3, inherit.aes = FALSE) +
  labs(x = "Geodesic distance", y = "Frequency") +
  theme(axis.text.x = element_text(angle = 40, vjust = 1, hjust = 1)) +
  ylim(-5, 450)

### Triad census.
triad_census <- matrix(nrow = 38, ncol = 16)
for (i in 1 : 38){
  triad_census[i, ] <- cell_75[[1]]$dtriad_results_het0$ergm_fits[[i]]$obs.triadcensus
}
dat1 <- melt(data.table(meta_ergm_fits$sim.triadcensus))
dat2 <- melt(data.table(triad_census))
dat2[ , 1] <- as.numeric(unlist(dat2[ , 1]))
dat1 <- as.data.frame(dat1)
dat2 <- as.data.frame(dat2)

dtriad_het0_triadc <- ggplot(dat1, aes(x = variable, y = value)) +
  geom_boxplot() +
  geom_jitter(data = dat2, aes(x = variable, y = value), position = position_jitter(0.1), color = "blue", shape = 18, alpha = 0.3, inherit.aes = FALSE) +
  labs(x = "Triad-census", y = "Frequency") +
  theme(axis.text.x = element_text(angle = 40, vjust = 1, hjust = 1)) +
  ylim(-5, 1300)

#############################################
###### Condition: 75 - dtriad - het1. #######
#############################################

coefs <- c(sum_cell_75[[1]]$dtriad_results_het1$ma_model[[1]]$b, sum_cell_75[[1]]$dtriad_results_het1$ma_model[[2]]$b,
           sum_cell_75[[1]]$dtriad_results_het1$ma_model[[3]]$b, sum_cell_75[[1]]$dtriad_results_het1$ma_model[[4]]$b,
           sum_cell_75[[1]]$dtriad_results_het1$ma_model[[5]]$b, sum_cell_75[[1]]$dtriad_results_het1$ma_model[[6]]$b)

simnet <- simulate(net ~ edges + mutual + gwesp(0.25, fixed = T) + gwodegree(0.10, fixed = T) + nodeocov("identclass") + nodematch("sex"), coef = coefs)
meta_ergm_fits <- gof(object = simnet ~  edges + mutual + gwesp(0.25, fixed = T) + gwodegree(0.10, fixed = T) + nodeocov("identclass") + nodematch("sex"),
                      GOF =~ model + idegree + odegree + distance + espartners + dspartners + triadcensus,
                      coef = coefs, control.gof.formula(seed = 3791, nsim = 1000))

### In-degree.
idegree <- matrix(nrow = 47, ncol = 25)
for (i in 1 : 47){
  idegree[i, ] <- cell_75[[1]]$dtriad_results_het1$ergm_fits[[i]]$obs.ideg
}
dat1 <- melt(data.table(meta_ergm_fits$sim.ideg))
dat2 <- melt(data.table(idegree))
dat2[ , 1] <- as.numeric(unlist(dat2[ , 1]))
dat1 <- as.data.frame(dat1)
dat2 <- as.data.frame(dat2)

dtriad_het1_indegree <- ggplot(dat1, aes(x = variable, y = value)) +
  geom_boxplot() +
  geom_jitter(data = dat2, aes(x = variable, y = value), position = position_jitter(0.1), color = "blue", shape = 18, alpha = 0.3, inherit.aes = FALSE) +
  labs(x = "In-degree", y = "Frequency") +
  theme(axis.text.x = element_text(angle = 40, vjust = 1, hjust = 1))+
  ylim(-1, 13)

### Out-degree.
odegree <- matrix(nrow = 47, ncol = 25)
for (i in 1 : 47){
  odegree[i, ] <- cell_75[[1]]$dtriad_results_het1$ergm_fits[[i]]$obs.odeg
}
dat1 <- melt(data.table(meta_ergm_fits$sim.odeg))
dat2 <- melt(data.table(odegree))
dat2[ , 1] <- as.numeric(unlist(dat2[ , 1]))
dat1 <- as.data.frame(dat1)
dat2 <- as.data.frame(dat2)

dtriad_het1_outdegree <- ggplot(dat1, aes(x = variable, y = value)) +
  geom_boxplot() +
  geom_jitter(data = dat2, aes(x = variable, y = value), position = position_jitter(0.1), color = "blue", shape = 18, alpha = 0.3, inherit.aes = FALSE) +
  labs(x = "Out-degree", y = "Frequency") +
  theme(axis.text.x = element_text(angle = 40, vjust = 1, hjust = 1)) +
  ylim(-1, 16)

### Geodesic distance.
dist <- matrix(nrow = 47, ncol = 25)
for (i in 1 : 47){
  dist[i, ] <- cell_75[[1]]$dtriad_results_het1$ergm_fits[[i]]$obs.dist
}
dat1 <- melt(data.table(meta_ergm_fits$sim.dist))
dat2 <- melt(data.table(dist))
dat2[ , 1] <- as.numeric(unlist(dat2[ , 1]))
dat1 <- as.data.frame(dat1)
dat2 <- as.data.frame(dat2)

dtriad_het1_dist <- ggplot(dat1, aes(x = variable, y = value)) +
  geom_boxplot() +
  geom_jitter(data = dat2, aes(x = variable, y = value), position = position_jitter(0.1), color = "blue", shape = 18, alpha = 0.3, inherit.aes = FALSE) +
  labs(x = "Geodesic distance", y = "Frequency") +
  theme(axis.text.x = element_text(angle = 40, vjust = 1, hjust = 1))+
  ylim(-5, 500)

### Triad census.
triad_census <- matrix(nrow = 47, ncol = 16)
for (i in 1 : 47){
  triad_census[i, ] <- cell_75[[1]]$dtriad_results_het1$ergm_fits[[i]]$obs.triadcensus
}

dat1 <- melt(data.table(meta_ergm_fits$sim.triadcensus))
dat2 <- melt(data.table(triad_census))
dat2[ , 1] <- as.numeric(unlist(dat2[ , 1]))
dat1 <- as.data.frame(dat1)
dat2 <- as.data.frame(dat2)

dtriad_het1_triadc <- ggplot(dat1, aes(x = variable, y = value)) +
  geom_boxplot() +
  geom_jitter(data = dat2, aes(x = variable, y = value), position = position_jitter(0.1), color = "blue", shape = 18, alpha = 0.3, inherit.aes = FALSE) +
  labs(x = "Triad-census", y = "Frequency") +
  theme(axis.text.x = element_text(angle = 40, vjust = 1, hjust = 1)) +
  ylim(-5, 1500)

### Drawing final plots ###. 
text1 <- "Structural model"
text2 <- "True model"
text3 <- "Dense triadic model"
# Create a text grob
tgrob1 <- text_grob(text1, size = 11)
tgrob2 <- text_grob(text2, size = 11)
tgrob3 <- text_grob(text3, size = 11)
# Draw the text
plot_1 <- as_ggplot(tgrob1) + theme(plot.margin = margin(2.5, -0.5, 0 ,0, "cm"))
plot_2 <- as_ggplot(tgrob2) + theme(plot.margin = margin(2.5, -0.5, 0 ,0, "cm"))
plot_3 <- as_ggplot(tgrob3) + theme(plot.margin = margin(2.5, -0.5, 0 ,0, "cm"))

### Arranged homogeneous sample.
ggpubr::ggarrange(struct_het0_indegree, true_het0_indegree, dtriad_het0_indegree,
                  struct_het0_outdegree, true_het0_outdegree, dtriad_het0_outdegree,
                  struct_het0_dist, true_het0_dist, dtriad_het0_dist,
                  struct_het0_triadc, true_het0_triadc, dtriad_het0_triadc,
                  ncol = 3, nrow = 4)
### Arranged heterogeneous sample.
ggpubr::ggarrange(struct_het1_indegree, true_het1_indegree, dtriad_het1_indegree,
                  struct_het1_outdegree, true_het1_outdegree, dtriad_het1_outdegree,
                  struct_het1_dist, true_het1_dist, dtriad_het1_dist,
                  struct_het1_triadc, true_het1_triadc, dtriad_het1_triadc,
                  ncol = 3, nrow = 4)
