#################################################
########## Code Master Thesis Project  ##########
#################################################

####################
##### Packages #####
####################

## "network" package.
# Check whether "network" package is installed and install if necessary.
if(!require(network)) install.packages("network")
# Push "network" package to library.
library(network)

## "ergm" package.
# Check whether "ergm" package is installed and install if necessary.
if(!require(ergm)) install.packages("ergm")
# Push "ergm" package to library.
library(ergm)

############################################
##### Importing Empricial Network Data #####
############################################

### Set working directory.
setwd("D:\\Academia\\Master MSBBSS\\Thesis\\Simulation folder")

### Import network files. 
ties_het0 <- as.matrix(read.table("not8.net"))
ties_het1 <- as.matrix(read.table("not16.net"))
### Import covariate files.
cov_het0 <- read.table("cov8.dat", col.names = c("girl", "identcGMC", "identcCMC", "identcCM", "percgirl")) 
cov_het1 <- read.table("cov16.dat", col.names = c("girl", "identcGMC", "identcCMC", "identcCM", "percgirl")) 

### Transform the matrices into network objects and attach sex variables to the network objects:
net_het0 <- network(ties_het0, directed = T) # Specify as directed network.
net_het0 %v% "sex" <- cov_het0$girl # Add sex actor variable. 
net_het0 %v% "identclass" <- cov_het0$identcGMC # Add class membership actor variable.

net_het1 <- network(ties_het1, directed = T) # Specify as directed network.
net_het1 %v% "sex" <- cov_het1$girl # Add sex actor variable. 
net_het1 %v% "identclass" <- cov_het1$identcGMC # Add class membership actor variable.

############################################################
##### Fitting an Exponential Random Graph Model (ERGM) #####
############################################################

### Set two random seeds.
seed1 <- 3791 # One ring to rule them all.
seed2 <- 2 * seed1 # Or two?

#### First true model #####
### Specifying true ERGM configuration.   
true_ergm_het0 <- net_het0 ~ edges + mutual + gwesp(0.25, fixed = T) + gwodegree(0.10, fixed = T) + nodeocov("identclass") + nodematch("sex")
### Fitting true ERG model.
## ERG model results.
# Fitting the ERG model.
true_ergm_results_het0 <- ergm(formula = true_ergm_het0, control = control.ergm(MCMC.interval = 10 * 1024, MCMC.burnin = (10 * 1024) * 16, MCMC.samplesize = 10 * 1024, 
                                                                                seed = seed1))

# Print ERG model estimation results.
summary(object = true_ergm_results_het0)

### Inspecting goodness-of-fit (GOF). In- and outdegree, triad-census, minimum geodesic distance.
## Calculate model and auxillary GOF.
true_ergm_gof_het0 <- gof(object = true_ergm_results_het0, GOF =~ idegree + odegree + triadcensus + distance + model, 
                          control.gof.formula(seed = seed2, nsim = 1000))
## Plot result.
plot(x = true_ergm_gof_het0) 

#### Second true model ####
### Specifying true ERGM configuration.   
true_ergm_het1 <- net_het1 ~ edges() + mutual() + gwesp(0.25, fixed = T) + gwodegree(0.10, fixed = T) + nodeocov("identclass") + nodematch("sex")
#search.ergmTerms()
### Fitting true ERG model.
## ERG model results.
# Fitting the ERG model.
true_ergm_results_het1 <- ergm(formula = true_ergm_het1, control = control.ergm(MCMC.interval = 10 * 1024, MCMC.burnin = (10 * 1024) * 16, MCMC.samplesize = 10 * 1024, 
                                                                                seed = seed1))

# Print ERG model estimation results.
summary(object = true_ergm_results_het1)

### Inspecting goodness-of-fit (GOF). In- and outdegree, triad-census, minimum geodesic distance.
## Calculate model and auxillary GOF.
true_ergm_gof_het1 <- gof(object = true_ergm_results_het1, GOF =~ idegree + odegree + triadcensus + distance + model, 
                          control.gof.formula(seed = seed2, nsim = 1000))
## Plot result.
plot(x = true_ergm_gof_het1) 